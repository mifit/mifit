#ifndef MI_CMapHeader_H
#define MI_CMapHeader_H

#include <map/maplib.h>

class CMapHeader : public CMapHeaderBase
{
public:
    CMapHeader();
    virtual bool LoadCrystal(const char*);
    virtual bool SaveCrystal(const std::string&);
};

#endif
