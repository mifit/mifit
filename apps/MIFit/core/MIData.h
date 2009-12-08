#ifndef mifit_core_MIData_h
#define mifit_core_MIData_h

#include <map>
#include <string>
#include <vector>

class MIDatum
{
public:
    MIDatum();

    bool b;
    int i;
    unsigned int u;
    short s;
    float f;
    double d;
    std::string str;
    std::vector<std::string> strList;

    // set this type instead of .u or .i for choices which should be
    // limited to a certain range by the randomizer
    unsigned int radio;
    unsigned int radio_count;
    std::vector<std::string> radio_labels;

    static const std::string INVALID_STRING;

};

typedef std::map<std::string, MIDatum> MIData;

#endif // ifndef mifit_core_MIData_h
