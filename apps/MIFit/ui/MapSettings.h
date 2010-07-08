#ifndef mifit_ui_MapSettings_h
#define mifit_ui_MapSettings_h

#include <map/maplib.h>

class MapSettings
{
    MapSettings(const MapSettings&);
    MapSettings& operator=(const MapSettings&);
public:
    static void save(const MapSettingsBase& settings);
    static void load(MapSettingsBase& settings);

    static void saveStyle(int i, const MapSettingsBase& settings);
    static void loadStyle(int i, MapSettingsBase& settings);

};

const unsigned int REGULAR_MAP_SETTINGS = 0;
const unsigned int DIFFERENCE_MAP_SETTINGS = 4;

// get settings for a particular style
//   will first get hardcoded defaults, and then overwrite with any user preferences
//   unless hardcodedOnly is true, in which case it will stop with the hardcoded values
bool GetMapSettingsForStyle(int styleno, MapSettingsBase &settings, bool hardcodedOnly = false);


#endif // ifndef mifit_ui_MapSettings_h
