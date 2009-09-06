#ifndef MI_RANDOM_TESTING_H
#define MI_RANDOM_TESTING_H

#include <stdio.h>

void MISetRandomPlaybackMode(FILE* logfile, bool playback);
long GetRand(const char* msg, long limit);

#endif
