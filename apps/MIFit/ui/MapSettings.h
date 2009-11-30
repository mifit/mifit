#ifndef mifit_ui_MapSettings_h
#define mifit_ui_MapSettings_h

#include <map/maplib.h>

class MapSettings : public MapSettingsBase
{
public:
    virtual void save();
    virtual void load();

    void saveStyle(int i);
    void loadStyle(int i);

};

const unsigned int REGULAR_MAP_SETTINGS = 0;
const unsigned int DIFFERENCE_MAP_SETTINGS = 4;

// get settings for a particular style
//   will first get hardcoded defaults, and then overwrite with any user preferences
//   unless hardcodedOnly is true, in which case it will stop with the hardcoded values
bool GetMapSettingsForStyle(int styleno, MapSettingsBase &settings, bool hardcodedOnly = false);

// MIData / MapSettings conversion
void MapSettingsToData(const MapSettingsBase &settings, MIData &data, float mapmin = -250.0f, float mapmax = 250.0f);
void DataToMapSettings(MIData &data, MapSettingsBase &settings, float &mapmin, float &mapmax);


#endif // ifndef mifit_ui_MapSettings_h
