#include "core/corelib.h"
#include "MapSettings.h"

#include "core/MIConfig.h"
#include "core/Colors.h"
#include <util/utillib.h>


static bool GetHardcodedDefaults(MapSettingsBase &settings, int styleno)
{
    //these are styleno-independent
    settings.m_radiusmax = 50;
    settings.Radius = 8.0;

    settings.ContourMethod = MAP_CUBE;
    settings.BlobRadius = 1.7f;
    settings.maplinewidth = 1.0;

    //these are styleno-dependent
    switch (styleno)
    {
    case 0:  // Blue Map 1,2,3,4,5 sigma
        settings.MapLevelColor[0] = Colors::MAP1;
        settings.MapLevelColor[1] = Colors::MAP2;
        settings.MapLevelColor[2] = Colors::MAP3;
        settings.MapLevelColor[3] = Colors::MAP4;
        settings.MapLevelColor[4] = Colors::MAP5;
        settings.MapLevelOn[0] = true;
        settings.MapLevelOn[1] = true;
        settings.MapLevelOn[2] = true;
        settings.MapLevelOn[3] = true;
        settings.MapLevelOn[4] = true;
        settings.MapLevel[0] = 50.0f;
        settings.MapLevel[1] = 100.0f;
        settings.MapLevel[2] = 150.0f;
        settings.MapLevel[3] = 200.0f;
        settings.MapLevel[4] = 250.0f;
        break;
    case 1:  // Single Level Blue 1 sigma
        settings.MapLevelColor[0] = Colors::MAP1;
        settings.MapLevelColor[1] = Colors::MAP2;
        settings.MapLevelColor[2] = Colors::MAP3;
        settings.MapLevelColor[3] = Colors::MAP4;
        settings.MapLevelColor[4] = Colors::MAP5;
        settings.MapLevelOn[0] = true;
        settings.MapLevelOn[1] = false;
        settings.MapLevelOn[2] = false;
        settings.MapLevelOn[3] = false;
        settings.MapLevelOn[4] = false;
        settings.MapLevel[0] = 50.0f;
        settings.MapLevel[1] = 100.0f;
        settings.MapLevel[2] = 150.0f;
        settings.MapLevel[3] = 200.0f;
        settings.MapLevel[4] = 250.0f;
        break;
    case 2: // Green Map 1,2,3,4,5 Level
        settings.MapLevelColor[0] = Colors::MAP6;
        settings.MapLevelColor[1] = Colors::MAP7;
        settings.MapLevelColor[2] = Colors::MAP8;
        settings.MapLevelColor[3] = Colors::MAP9;
        settings.MapLevelColor[4] = Colors::MAP10;
        settings.MapLevelOn[0] = true;
        settings.MapLevelOn[1] = true;
        settings.MapLevelOn[2] = true;
        settings.MapLevelOn[3] = true;
        settings.MapLevelOn[4] = true;
        settings.MapLevel[0] = 50.0f;
        settings.MapLevel[1] = 100.0f;
        settings.MapLevel[2] = 150.0f;
        settings.MapLevel[3] = 200.0f;
        settings.MapLevel[4] = 250.0f;
        break;
    case 3: // Single Level Green 1 sigma
        settings.MapLevelColor[0] = Colors::MAP6;
        settings.MapLevelColor[1] = Colors::MAP7;
        settings.MapLevelColor[2] = Colors::MAP8;
        settings.MapLevelColor[3] = Colors::MAP9;
        settings.MapLevelColor[4] = Colors::MAP10;
        settings.MapLevelOn[0] = true;
        settings.MapLevelOn[1] = false;
        settings.MapLevelOn[2] = false;
        settings.MapLevelOn[3] = false;
        settings.MapLevelOn[4] = false;
        settings.MapLevel[0] = 50.0f;
        settings.MapLevel[1] = 100.0f;
        settings.MapLevel[2] = 150.0f;
        settings.MapLevel[3] = 200.0f;
        settings.MapLevel[4] = 250.0f;
        break;
    case 4: // Difference Map -4,-3,3,4,5 sigma
        settings.MapLevelColor[0] = Colors::MAP4;
        settings.MapLevelColor[1] = Colors::MAP5;
        settings.MapLevelColor[2] = Colors::MAP1;
        settings.MapLevelColor[3] = Colors::MAP2;
        settings.MapLevelColor[4] = Colors::MAP3;
        settings.MapLevelOn[0] = true;
        settings.MapLevelOn[1] = true;
        settings.MapLevelOn[2] = true;
        settings.MapLevelOn[3] = true;
        settings.MapLevelOn[4] = true;
        settings.MapLevel[0] = -200.0f;
        settings.MapLevel[1] = -150.0f;
        settings.MapLevel[2] = 150.0f;
        settings.MapLevel[3] = 200.0f;
        settings.MapLevel[4] = 250.0f;
        break;
    case 5: // Difference Map -3 red,+3 blue sigma
        settings.MapLevelColor[0] = Colors::MAP4;
        settings.MapLevelColor[1] = Colors::MAP5;
        settings.MapLevelColor[2] = Colors::MAP1;
        settings.MapLevelColor[3] = Colors::MAP2;
        settings.MapLevelColor[4] = Colors::MAP3;
        settings.MapLevelOn[0] = false;
        settings.MapLevelOn[1] = true;
        settings.MapLevelOn[2] = true;
        settings.MapLevelOn[3] = false;
        settings.MapLevelOn[4] = false;
        settings.MapLevel[0] = -200.0f;
        settings.MapLevel[1] = -150.0f;
        settings.MapLevel[2] = 150.0f;
        settings.MapLevel[3] = 200.0f;
        settings.MapLevel[4] = 250.0f;
        break;
    default:
        return false;
        break;
    }
    return true;
}



void MapSettings::saveStyle(int i, const MapSettingsBase& settings)
{
    std::string key;
    std::string stylestr("");
    if (i!=-1)
        stylestr = format("%d", i);
    const char *sty = stylestr.c_str();

    key = format("MapSettings%s/m_radiusmax", sty);
    MIConfig::Instance()->Write(key, (long)settings.m_radiusmax);
    key = format("MapSettings%s/Radius", sty);
    MIConfig::Instance()->Write(key, settings.Radius);
    key = format("MapSettings%s/ContourMethod", sty);
    MIConfig::Instance()->Write(key, (long)settings.ContourMethod);
    key = format("MapSettings%s/BlobRadius", sty);
    MIConfig::Instance()->Write(key, settings.BlobRadius);
    key = format("MapSettings%s/maplinewidth", sty);
    MIConfig::Instance()->Write(key, settings.maplinewidth);

    for (int i = 0; i < 5; ++i)
    {
        key = format("MapSettings%s/MapLevelColor%d", sty, i);
        MIConfig::Instance()->Write(key, (long)settings.MapLevelColor[i]);
    }

    for (int i = 0; i < 5; ++i)
    {
        key = format("MapSettings%s/MapLevelOn%d", sty, i);
        MIConfig::Instance()->Write(key, settings.MapLevelOn[i]);
    }

    for (int i = 0; i < 5; ++i)
    {
        key = format("MapSettings%s/MapLevel%d", sty, i);
        MIConfig::Instance()->Write(key, (double)settings.MapLevel[i]);
    }
}

void MapSettings::loadStyle(int i, MapSettingsBase& settings)
{

    // first get the hard-coded defaults
    GetHardcodedDefaults(settings, i);

    // now try in the user prefs file, using the hardcoded values as defaults
    // in case they're not in the settings file
    std::string key;
    std::string stylestr("");
    if (i!=-1)
        stylestr = format("%d", i);
    const char *sty = stylestr.c_str();

    long longValue;
    double d;

    key = format("MapSettings%s/m_radiusmax", sty);
    MIConfig::Instance()->Read(key, &longValue, (long)settings.m_radiusmax);
    settings.m_radiusmax = longValue;
    key = format("MapSettings%s/Radius", sty);
    MIConfig::Instance()->Read(key, &d, settings.Radius);
    settings.Radius = (float) d;
    key = format("MapSettings%s/ContourMethod", sty);
    MIConfig::Instance()->Read(key, &longValue, (long)settings.ContourMethod);
    settings.ContourMethod = longValue;
    key = format("MapSettings%s/BlobRadius", sty);
    MIConfig::Instance()->Read(key, &d, settings.BlobRadius);
    settings.BlobRadius = (float) d;
    key = format("MapSettings%s/maplinewidth", sty);
    MIConfig::Instance()->Read(key, &d, settings.maplinewidth);
    settings.maplinewidth = d;

    for (int i = 0; i < 5; ++i)
    {
        key = format("MapSettings%s/MapLevelColor%d", sty, i);
        MIConfig::Instance()->Read(key, &longValue, (long)settings.MapLevelColor[i]);
        settings.MapLevelColor[i] = longValue;
    }
    for (int i = 0; i < 5; ++i)
    {
        key = format("MapSettings%s/MapLevelOn%d", sty, i);
        MIConfig::Instance()->Read(key, &settings.MapLevelOn[i], settings.MapLevelOn[i]);
    }
    for (int i = 0; i < 5; ++i)
    {
        key = format("MapSettings%s/MapLevel%d", sty, i);
        MIConfig::Instance()->Read(key, &d, settings.MapLevel[i]);
        settings.MapLevel[i] = d;
    }
}

void MapSettings::save(const MapSettingsBase& settings)
{
    saveStyle(-1, settings);
}

void MapSettings::load(MapSettingsBase& settings)
{
    loadStyle(-1, settings);
}

bool GetMapSettingsForStyle(int styleno, MapSettingsBase &settings, bool hardcodedOnly)
{
    if (hardcodedOnly)
    {
        return GetHardcodedDefaults(settings, styleno);
    }

    MapSettings::loadStyle(styleno, settings);
    return true;
}
