#include "EMapBase.h"
#include "MapSettingsBase.h"

void MapSettingsBase::load()
{
    static float defaultMapLevel[5];
    static int defaultMapLevelColor[5];
    defaultMapLevel[0] = 50.0;
    defaultMapLevel[1] = 100.0;
    defaultMapLevel[2] = 150.0;
    defaultMapLevel[3] = 200.0;
    defaultMapLevel[4] = 250.0;

    // note colors won't be used in a non-gui application, so we could
    // probably just skip this.
    defaultMapLevelColor[0] = 21; //Colors::MAP1;
    defaultMapLevelColor[1] = 22; //Colors::MAP2;
    defaultMapLevelColor[2] = 23; //Colors::MAP3;
    defaultMapLevelColor[3] = 24; //Colors::MAP4;
    defaultMapLevelColor[4] = 25; //Colors::MAP5;

    ContourMethod = MAP_CUBE;
    BlobRadius = 1.7f;
    Radius = 8.0f;

    for (int i = 0; i < 5; ++i)
    {
        MapLevel[i] = defaultMapLevel[i];
        MapLevelOn[i] = true;
        MapLevelColor[i] = defaultMapLevelColor[i];
    }
    maplinewidth = 1.0f;
    m_radiusmax = 50;
}

