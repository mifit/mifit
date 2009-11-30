#ifndef mifit_mimolio_h
#define mifit_mimolio_h

#include <chemlib/chemlib.h>

class MIMolIO : public chemlib::MIMolIOBase
{
public:
    MIMolIO();
    bool Write(chemlib::MIMolInfo &mol, const std::string &filename, int writerIndex = -1) const;
    bool Read(chemlib::MIMolInfo &mol, const std::string &filename, int readerIndex = -1) const;
};

#endif
