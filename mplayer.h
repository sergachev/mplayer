
#ifndef __MPLAYER_MAIN
#define __MPLAYER_MAIN

#include "libvo/sub.h"
#include "subreader.h"

extern int use_gui;
extern char* current_module;
extern int fullscreen;
extern int vcd_track;

extern char ** audio_fm_list;
extern char ** video_fm_list;
extern char * video_driver;
extern char * audio_driver;
extern float  audio_delay;

extern int osd_level;
extern int osd_visible;

extern char * font_name;
extern float  font_factor;
extern float movie_aspect;

extern char * sub_name;
extern float  sub_delay;
extern float  sub_fps;
extern int    sub_auto;
extern int    sub_pos;
extern int    sub_unicode;
extern subtitle* subtitles;
extern subtitle* vo_sub;

extern char * filename;

extern int stream_cache_size;
extern int flip;
extern int force_ni;
extern int index_mode;
extern int frame_dropping;

extern int auto_quality;

extern void exit_player(char* how);

#endif
