/*
 * This file is part of MPlayer.
 *
 * MPlayer is free software; you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation; either version 2 of the License, or
 * (at your option) any later version.
 *
 * MPlayer is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License along
 * with MPlayer; if not, write to the Free Software Foundation, Inc.,
 * 51 Franklin Street, Fifth Floor, Boston, MA 02110-1301 USA.
 */

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include "mp_image.h"
#include "mp_msg.h"
#include "vf.h"
#include "img_format.h"

#include "libvo/fastmemcpy.h"
#include "libavutil/common.h"

//////////////////////////////////////////////////

struct vf_priv_s {
    int w, h;
    int* inv_cache;
    unsigned int fmt;
};


const int NumSegments = 11;

// From OVR_Stereo.cpp, evaluates Catmull-Rom spline based on provided K values
static float EvalCatmullRom10Spline(float const *K, float scaledVal)
{
    float t, omt, res;
    float p0, p1;
    float m0, m1;
    int k;

    float scaledValFloor = floorf(scaledVal);
    scaledValFloor = fmaxf(0.0f, fminf((float)(NumSegments - 1), scaledValFloor));
    t = scaledVal - scaledValFloor;
    k = (int)scaledValFloor;

    switch (k)
    {
        case 0:
            // Curve starts at 1.0 with gradient K[1]-K[0]
            p0 = 1.0f;
            m0 = (K[1] - K[0]);    // general case would have been (K[1]-K[-1])/2
            p1 = K[1];
            m1 = 0.5f * (K[2] - K[0]);
            break;
        default:
            // General case
            p0 = K[k];
            m0 = 0.5f * (K[k + 1] - K[k - 1]);
            p1 = K[k + 1];
            m1 = 0.5f * (K[k + 2] - K[k]);
            break;
        case NumSegments - 2:
            // Last tangent is just the slope of the last two points.
            p0 = K[NumSegments - 2];
            m0 = 0.5f * (K[NumSegments - 1] - K[NumSegments - 2]);
            p1 = K[NumSegments - 1];
            m1 = K[NumSegments - 1] - K[NumSegments - 2];
            break;
        case NumSegments - 1:
            // Beyond the last segment it's just a straight line
            p0 = K[NumSegments - 1];
            m0 = K[NumSegments - 1] - K[NumSegments - 2];
            p1 = p0 + m0;
            m1 = m0;
            break;
    }

    omt = 1.0f - t;
    res = (p0 * (1.0f + 2.0f *   t) + m0 *   t) * omt * omt
            + (p1 * (1.0f + 2.0f * omt) - m1 * omt) *   t *   t;

    return res;
}

// From OVR_DeviceConstants.h
enum DistortionEqnType
{
    Distortion_Poly4 = 0,
    Distortion_RecipPoly4 = 1,
    Distortion_CatmullRom10 = 2,
};

// Also based on function from OVR_Stereo.cpp
static float DistortionFnScaleRadiusSquared(enum DistortionEqnType Eqn, float const *K, float MaxR, float const CA0, float const CA1, float rsq)
{
    float scale = 1.0f;
    switch (Eqn)
    {
        case Distortion_Poly4:
            // This version is deprecated! Prefer one of the other two.
            scale = (K[0] + rsq * (K[1] + rsq * (K[2] + rsq * K[3])));
            break;
        case Distortion_RecipPoly4:
            scale = 1.0f / (K[0] + rsq * (K[1] + rsq * (K[2] + rsq * K[3])));
            break;
        case Distortion_CatmullRom10:{
            // A Catmull-Rom spline through the values 1.0, K[1], K[2] ... K[10]
            // evenly spaced in R^2 from 0.0 to MaxR^2
            // K[0] controls the slope at radius=0.0, rather than the actual value.
            float scaledRsq = (float)(NumSegments - 1) * rsq / (MaxR * MaxR);
            scale = EvalCatmullRom10Spline(K, scaledRsq);
        }break;
    }
    scale *= 1.0f + CA0 + CA1 * rsq;
    return scale;
}

// Computes inverse of DistortionFnScaleRadiusSquared function using binary search
// Function is monotonic increasing so this ought to work, although might be slow
static float DistortionFnScaleRadiusSquaredInv(enum DistortionEqnType Eqn, float const *K, float MaxR, float const CA0, float const CA1, float rsq)
{
    float low_guess = 0.0f, high_guess = 10.0f;
    // The "high_guess > 0.00001" is needed for the singular case where zero is the solution
    // With the relative error at 0.001 I observed a dot in the center on some frames, so lowered to 0.0001
    while ((high_guess - low_guess) / low_guess > 0.0001 && high_guess > 0.00001) {
        float mid_guess = (low_guess + high_guess) / 2.0f;
        float scale = DistortionFnScaleRadiusSquared(Eqn, K, MaxR, CA0, CA1, mid_guess);
        float mid_guess_value = scale * scale * mid_guess;
        if (rsq < mid_guess_value) {
            high_guess = mid_guess;
        }
        else {
            low_guess = mid_guess;
        }
    }
    return (low_guess + high_guess) / 2.0f;
}

#define NUM_EYES 2
#define NUM_CHANNELS 4


static unsigned int getfmt(unsigned int outfmt){
    switch(outfmt){
        case IMGFMT_BGRA:
            return outfmt;
        default:
            return 0;
    }
}


static int config_props(vf_instance_t *vf)
{
    int64_t w, h;
    int ret;


    w = vf->priv->w;//unwarpvr->w;
    h = vf->priv->h;//unwarpvr->h;


    // Initialize inv_cache
    {
        int i, j, eye_count;
        float MetersPerTanAngleAtCenter;
        float screenWidthMeters;
        float screenHeightMeters;
        float LensCenterXOffset; // For left eye, determined by physical parameters
        float TanEyeAngleScaleX, TanEyeAngleScaleY, DevicePPDInCenterX, DevicePPDInCenterY;
        enum DistortionEqnType Eqn = Distortion_CatmullRom10;
        float K[11];
        float MaxR = 1.0f;
        float ChromaticAberration[4];
        int DeviceResX, DeviceResY;
        int channel;
        int in_linesize;

        // Create temporary input frame just so we can get its linesize
        mp_image_t *dmpi;
        dmpi=vf_get_image(vf->next, vf->priv->fmt,
                MP_IMGTYPE_TEMP, MP_IMGFLAG_ACCEPT_STRIDE | MP_IMGFLAG_PREFER_ALIGNED_STRIDE,
                w, h);
        in_linesize = dmpi->stride[0]; //w * NUM_CHANNELS; //in->linesize[0];
        free_mp_image(dmpi);

        // Distortion varies by SDK version but never by cup type or eye relief (for DK2 in 0.4.2)
        const float K_DK2[] = { 1.003f, 1.02f, 1.042f, 1.066f, 1.094f, 1.126f, 1.162f, 1.203f, 1.25f, 1.31f, 1.38f };
        // ChromaticAbberation varies by eye relief and lerps between the following two arrays
        const float ChromaticAberrationMin[] = { -0.0112f, -0.015f, 0.0187f, 0.015f };
        const float ChromaticAberrationMax[] = { -0.015f, -0.02f, 0.025f, 0.02f };
        memmove(K, K_DK2, sizeof(K));
        int eye_relief_dial=3;
        for (i = 0; i < sizeof(ChromaticAberration) / sizeof(*ChromaticAberration); i++) {
            ChromaticAberration[i] = ChromaticAberrationMin[i] + eye_relief_dial / 10.0f * (ChromaticAberrationMax[i] - ChromaticAberrationMin[i]);
        }

        MetersPerTanAngleAtCenter = 0.036f;
        screenWidthMeters = 0.12576f;
        screenHeightMeters = 0.07074f;
        LensCenterXOffset = -0.00986003876f;
        DeviceResX = 1920; DeviceResY = 1080;

        DevicePPDInCenterX = MetersPerTanAngleAtCenter / screenWidthMeters * DeviceResX;
        DevicePPDInCenterY = MetersPerTanAngleAtCenter / screenHeightMeters * DeviceResY;

        float scale_in_width = 1, scale_in_height = 1;
        float ppd = 0;
        if (ppd != 0.0f) {
            scale_in_width *= (ppd * 53.1301f) / DevicePPDInCenterX; // 53.1301 deg = tan(0.5) - (tan-0.5)
            scale_in_height *= (ppd * 53.1301f) / DevicePPDInCenterY;
        }

        // As computed in CalculateDistortionRenderDesc() distortion.TanEyeAngleScale in OVR_Stereo.cpp
        TanEyeAngleScaleX = 0.25f * screenWidthMeters / MetersPerTanAngleAtCenter;
        TanEyeAngleScaleY = 0.5f * screenHeightMeters / MetersPerTanAngleAtCenter;

        vf->priv->inv_cache = malloc(w * h * NUM_CHANNELS * sizeof(int));
        if (!vf->priv->inv_cache)
        {
            //av_log(ctx, AV_LOG_ERROR, "unwarpvr: Out of memory allocating cache\n");
            return -1;//AVERROR(EINVAL);
        }
        for (i = 0; i < h * w * NUM_CHANNELS; i++) {
            vf->priv->inv_cache[i] = -1;
        }
        for (eye_count = 0; eye_count < NUM_EYES; eye_count++) {
            int mono_input=1;
            int one_eye_multiplier = 1;
            float lensCenterXOffsetEye;
            int in_eye = eye_count;
            int out_eye = eye_count;
            int in_width_per_eye = w;
            if (mono_input)
                in_eye = 0;
            int forward_warp=1;
            lensCenterXOffsetEye = ((!forward_warp && in_eye) || (forward_warp && out_eye)) ? -LensCenterXOffset : LensCenterXOffset;


            for (i = 0; i < h; i++) {
                for (j = 0; j < w / 2 * one_eye_multiplier; j++) {
                    float ndcx, ndcy, tanx_distorted, tany_distorted, rsq;
                    float scale[NUM_CHANNELS];
                    float scale_width=1, scale_height=1;

                    ndcx = ((-1.0f + 2.0f * j / (w / 2 * one_eye_multiplier)) * one_eye_multiplier) / scale_width - lensCenterXOffsetEye;
                    ndcy = (-1.0f + 2.0f * i / h) / scale_height;
                    tanx_distorted = ndcx * TanEyeAngleScaleX;
                    tany_distorted = ndcy * TanEyeAngleScaleY;
                    rsq = tanx_distorted*tanx_distorted + tany_distorted*tany_distorted;
                    scale[0] = DistortionFnScaleRadiusSquared(Eqn, K, MaxR, ChromaticAberration[0], ChromaticAberration[1], rsq);
                    scale[1] = DistortionFnScaleRadiusSquared(Eqn, K, MaxR, 0, 0, rsq);
                    scale[2] = DistortionFnScaleRadiusSquared(Eqn, K, MaxR, ChromaticAberration[2], ChromaticAberration[3], rsq);
                    scale[3] = 1;
                    for (channel = 0; channel < NUM_CHANNELS; channel++) {
                        float x, y;
                        int srcj, srci;
                        float tanx, tany, rt_ndcx, rt_ndcy;
                        int output_idx = (i*w + eye_count*w / 2 + j)*NUM_CHANNELS + channel;

                        tanx = tanx_distorted * scale[channel];
                        tany = tany_distorted * scale[channel];

                        rt_ndcx = tanx / TanEyeAngleScaleX;
                        rt_ndcy = tany / TanEyeAngleScaleY;

                        x = (rt_ndcx * scale_in_width / 2.0f * DeviceResX / 2) + (in_width_per_eye / 2.0f);
                        y = (rt_ndcy * scale_in_height / 2.0f * DeviceResY) + (h / 2.0f);

                        srcj = (int)x;
                        srci = (int)y;

                        if (srci >= 0 && srcj >= 0 && srci < h && srcj < in_width_per_eye)
                            vf->priv->inv_cache[output_idx] = (srci*in_linesize) + (in_eye * in_width_per_eye + srcj)*NUM_CHANNELS + channel;
                    }
                }
            }
        }
    }

    return 0;

}


/////////////////////////////////////////////////




static int config(struct vf_instance *vf,
       int width, int height, int d_width, int d_height,
       unsigned int flags, unsigned int outfmt)
{
    if (vf->priv->w < 0 || width < vf->priv->w)
        vf->priv->w = width;
    if (vf->priv->h < 0 || height < vf->priv->h)
        vf->priv->h = height;
//        mp_msg(MSGT_VFILTER,MSGL_WARN,"rectangle: bad position/width/height - rectangle area is out of the original!\n");
    vf->priv->fmt=getfmt(outfmt);
    config_props(vf);

    return vf_next_config(vf, width, height, d_width, d_height, flags, outfmt);
}

static int control(struct vf_instance *vf, int request, void *data)
{
    return vf_next_control(vf, request, data);
}

static int put_image(struct vf_instance *vf, mp_image_t* mpi, double pts){
    mp_image_t* dmpi;
    dmpi = vf_get_image(vf->next, mpi->imgfmt, MP_IMGTYPE_TEMP,
                        MP_IMGFLAG_ACCEPT_STRIDE | MP_IMGFLAG_PREFER_ALIGNED_STRIDE,
                        mpi->w, mpi->h);

    uint8_t *dst = dmpi->planes[0];
    int *inv_cache_p = vf->priv->inv_cache;
    for (int y = 0; y < dmpi->height; y++) {
        for (int x = 0; x < dmpi->width*NUM_CHANNELS; x++, inv_cache_p++) {
            dst[x] = (*inv_cache_p == -1) ? 0 : mpi->planes[0][*inv_cache_p];
        }
        dst+=dmpi->stride[0];
    }


    return vf_next_put_image(vf, dmpi, pts);
}


static int query_format(struct vf_instance *vf, unsigned int outfmt){
    unsigned int fmt=getfmt(outfmt);
    if(!fmt) return 0;
    return vf_next_query_format(vf,fmt) & (~VFCAP_CSP_SUPPORTED_BY_HW);
}


static int vf_open(vf_instance_t *vf, char *args) {
    vf->query_format=query_format;
    vf->config = config;
    vf->control = control;
    vf->put_image = put_image;
    vf->priv = malloc(sizeof(struct vf_priv_s));
    vf->priv->w = -1;
    vf->priv->h = -1;
    return 1;
}

const vf_info_t vf_info_unwarpvr = {
    "warp image for oculus rift",
    "unwarpvr",
    "eVRydayVR",
    "",
    vf_open,
    NULL
};
