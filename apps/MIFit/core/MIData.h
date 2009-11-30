#ifndef mifit_core_MIData_h
#define mifit_core_MIData_h

#include <string>
#include <vector>
#include <map>

#include <float.h>
#include "stdio.h"

class MIDatum
{
public:
    MIDatum();

    // returns number of characters read, or 0 if error
    unsigned int ReadDatum(const std::string &dat);
    bool WriteDatum(std::string &str) const;

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

    bool isColorIndex;
    bool isColor;
    unsigned char color[3];

    static const std::string INVALID_STRING;

};

typedef std::map<std::string, MIDatum> MIData;
typedef std::map<std::string, MIDatum>::iterator MIDataIter;
typedef std::map<std::string, MIDatum>::const_iterator MIDataConstIter;

bool StringToMIData(const std::string &str, MIData &data, std::string *errorString = NULL);
std::string MIDataToString(MIData &data);


#endif // ifndef mifit_core_MIData_h
