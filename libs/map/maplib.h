#ifndef MI_MAPLIB_H
#define MI_MAPLIB_H

#include <string>

//private #include "sfcalc.h"

int MIMapFactor(int ntest, int prime, int even, int inc);
unsigned int MIMapInitializeScatteringFactorTables(
    const std::string &crystal_data_dir,
    const std::string &molimage_home_dir);
unsigned int MIMapFreeScatteringFactorTables();

#include "CMapHeaderBase.h"
#include "EMapBase.h"
#include "CMapHeaderBase.h"
#include "MAP_POINT.h"
#include "MapSettingsBase.h"
#include "InterpBox.h"
#include "maptypes.h"
#include "CREFL.h"

// private #include "fssubs.h"
// private #include "fft.h"
// private #include "MINATOM.h"
// private #include "PEAK.h"
// private #include "rescalc.h"
// private #include "sfcalc_data.h"
// private #include "sfcalc.h"



#endif // ifndef MI_MAPLIB_H
