#include "EMapBase.h"
#include "MapSettingsBase.h"

MapSettingsBase::MapSettingsBase()
{
    static const float defaultMapLevel[5] = { 50.0f, 100.0f, 150.0f, 200.0f, 250.0f };
    static const int defaultMapLevelColor[5] =
    {
        // note colors won't be used in a non-gui application, so we could
        // probably just skip this.
        21, //Colors::MAP1;
        22, //Colors::MAP2;
        23, //Colors::MAP3;
        24, //Colors::MAP4;
        25  //Colors::MAP5;
    };

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

