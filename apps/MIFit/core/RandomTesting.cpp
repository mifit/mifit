#include "nonguilib.h"
#include <stdlib.h>
#include "RandomTesting.h"

#ifdef _WIN32
#define random rand
#endif

static FILE* LOGFILE = 0;
bool PlayBackMode = false;
void MISetRandomPlaybackMode(FILE* logfile, bool playback) {
  LOGFILE = logfile;
  PlayBackMode = playback;
}

long GetRand(const char* msg, long limit) {
  if (PlayBackMode) {
    long r;
    char buf[1024];
    do {
      buf[0] = ' ';
      fgets(buf, 1024, LOGFILE);
      if (buf[0] != 'M' && sscanf(buf, "%ld", &r) == 1) {
        return r;
      }
    } while (buf[0] == 'M');

    Logger::log("Out of data in playback!\n");
    return 0;
  }

  long r = 0;
  if (limit > 0) {
    r = random()%limit;
  }
  if (LOGFILE) {
    fprintf(LOGFILE, "%8ld of %8ld; %s\n", r, limit, msg ? msg : "");
    fflush(LOGFILE);
  }
  return r;
}

