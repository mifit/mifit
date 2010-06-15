#ifndef MIMolIOBase_h
#define MIMolIOBase_h

#include <cstdio>
#include <vector>
#include <string>

#include "model.h"
#include "MIAtom_fwd.h"
#include "Residue_fwd.h"
#include "Bond.h"
#include "ANGLE.h"
#include "PLANE.h"
#include "CHIRALDICT.h"
#include "TORSION.h"

namespace chemlib
{

    class MIMolInfo
    {
    public:
        // readers and writers use these
        MIMolInfo() : res(0)
        {
        }

        // readers allocate and return a new residue for this,
        // writers expect this to be set to the residue to write
        Residue *res;                   // could make a vector, prob not worth the trouble

        std::vector<Bond> bonds;
        std::vector<ANGLE> angles;       // mmCIF only

        // some readers populate these: mmCIF,
        std::vector<TORSDICT> tordict;   // mmCIF only
        std::vector<PLANEDICT> planedict; // mmCIF only
        std::vector<CHIRALDICT> chiralsdict; // mmCIF only

        // some writers use these
        std::vector<TORSION> torsions; // mmCIF only
        std::vector<PLANE> planes;   // mmCIF only
        std::vector<CHIRAL> chirals; // mmCIF only
        std::vector<std::string> header; // PDB only
        std::vector<std::string> tail; // PDB only
    };

    class Reader
    {
    public:
        virtual ~Reader()
        {
        }

        virtual std::string getDescription() = 0;
        virtual std::string getExtension() = 0;
        virtual bool Read(FILE *fp, MIMolInfo &mi) = 0;
    };

    class Writer
    {
    public:
        virtual ~Writer()
        {
        }

        virtual std::string getDescription() = 0;
        virtual std::string getExtension() = 0;
        virtual bool Write(FILE *fp, MIMolInfo &mi) = 0;
    };

}

#endif // ifndef MIMolIOBase_h
