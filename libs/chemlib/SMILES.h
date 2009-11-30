#ifndef mifit_fios_SMILES_h
#define mifit_fios_SMILES_h

#include <cstdio>
#include <string>
#include <map>
#include <vector>
#include "MIMolIOBase.h"

namespace chemlib
{

    class SMILES : public Reader
    {
    public:
        SMILES();
        virtual ~SMILES();
        virtual std::string getDescription()
        {
            return "SMILES (*.smi)";
        }

        virtual std::string getExtension()
        {
            return "*.smi";
        }

        virtual bool Read(FILE *fp, MIMolInfo &mol);
        bool Read(const std::string &smiles, MIMolInfo &mol);
    };

}

#endif // ifndef mifit_fios_SMILES_h
