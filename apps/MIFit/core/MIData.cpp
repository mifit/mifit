#include "MIData.h"

#include <cfloat>
#include <climits>

const std::string MIDatum::INVALID_STRING("INVALID_STRING");

MIDatum::MIDatum()
{
    b = 0;
    i = INT_MIN;
    u = UINT_MAX;
    s = SHRT_MIN;
    f = FLT_MIN;
    d = DBL_MIN;
    radio = UINT_MAX;
    radio_count = UINT_MAX;
    str = INVALID_STRING;
    strList.clear();
}
