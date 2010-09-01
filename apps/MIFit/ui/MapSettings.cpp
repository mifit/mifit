#include "core/corelib.h"
#include "MapSettings.h"

#include "core/Colors.h"
#include <util/utillib.h>
#include <QSettings>

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
    QString key;
    QString sty("");
    if (i!=-1)
        sty = QString::number(i);

    QSettings qsettings;
    key = QString("MapSettings%1/m_radiusmax").arg(sty);
    qsettings.setValue(key, settings.m_radiusmax);
    key = QString("MapSettings%1/Radius").arg(sty);
    qsettings.setValue(key, settings.Radius);
    key = QString("MapSettings%1/ContourMethod").arg(sty);
    qsettings.setValue(key, settings.ContourMethod);
    key = QString("MapSettings%1/BlobRadius").arg(sty);
    qsettings.setValue(key, settings.BlobRadius);
    key = QString("MapSettings%1/maplinewidth").arg(sty);
    qsettings.setValue(key, settings.maplinewidth);

    for (int i = 0; i < 5; ++i)
    {
        key = QString("MapSettings%1/MapLevelColor%2").arg(sty).arg(i);
        qsettings.setValue(key, settings.MapLevelColor[i]);
    }

    for (int i = 0; i < 5; ++i)
    {
        key = QString("MapSettings%1/MapLevelOn%2").arg(sty).arg(i);
        qsettings.setValue(key, settings.MapLevelOn[i]);
    }

    for (int i = 0; i < 5; ++i)
    {
        key = QString("MapSettings%1/MapLevel%2").arg(sty).arg(i);
        qsettings.setValue(key, settings.MapLevel[i]);
    }
}

void MapSettings::loadStyle(int i, MapSettingsBase& settings)
{

    // first get the hard-coded defaults
    GetHardcodedDefaults(settings, i);

    // now try in the user prefs file, using the hardcoded values as defaults
    // in case they're not in the settings file
    QString key;
    QString sty("");
    if (i!=-1)
        sty = QString::number(i);

    QSettings qsettings;
    key = QString("MapSettings%1/m_radiusmax").arg(sty);
    settings.m_radiusmax = qsettings.value(key, settings.m_radiusmax).toInt();
    key = QString("MapSettings%1/Radius").arg(sty);
    settings.Radius = (float) qsettings.value(key, settings.Radius).toDouble();
    key = QString("MapSettings%1/ContourMethod").arg(sty);
    settings.ContourMethod = qsettings.value(key, settings.ContourMethod).toInt();
    key = QString("MapSettings%1/BlobRadius").arg(sty);
    settings.BlobRadius = (float) qsettings.value(key, settings.BlobRadius).toDouble();
    key = QString("MapSettings%1/maplinewidth").arg(sty);
    settings.maplinewidth = qsettings.value(key, settings.maplinewidth).toDouble();

    for (int i = 0; i < 5; ++i)
    {
        key = QString("MapSettings%1/MapLevelColor%2").arg(sty).arg(i);
        settings.MapLevelColor[i] = qsettings.value(key, settings.MapLevelColor[i]).toInt();
    }
    for (int i = 0; i < 5; ++i)
    {
        key = QString("MapSettings%1/MapLevelOn%2").arg(sty).arg(i);
        settings.MapLevelOn[i] = qsettings.value(key, settings.MapLevelOn[i]).toInt();
    }
    for (int i = 0; i < 5; ++i)
    {
        key = QString("MapSettings%1/MapLevel%2").arg(sty).arg(i);
        settings.MapLevel[i] = qsettings.value(key, settings.MapLevel[i]).toDouble();
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
