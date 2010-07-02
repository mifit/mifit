#ifndef mifit_fios_PDB_h
#define mifit_fios_PDB_h
#include "MIMolIOBase.h"


namespace chemlib
{
    class MIAtom;
    class Bond;
    class Residue;

    class PDB : public Reader, public Writer
    {

        void CreateBond(MIAtom *atom1, MIAtom *atom2, std::vector<Bond> &bonds);

    public:
        PDB();
        virtual ~PDB();
        virtual std::string getDescription()
        {
            return "PDB (*.pdb)";
        }

        virtual std::string getExtension()
        {
            return "*.pdb";
        }

        virtual bool Write(FILE *fp, MIMolInfo &mol);
        virtual bool Read(FILE *fp, MIMolInfo &mol);
    };

    void MISetIgnoreDummyAtomsOnLoad(bool ignore);
    Residue *LoadPDB(FILE *f, std::vector<Bond> *connects);
    bool SavePDB(FILE *fp, ResidueListIterator beginRes, ResidueListIterator endRes, Bond *Connects, int nConnects, bool mark_end = true, std::vector<std::string> *head = 0, std::vector<std::string> *tail = 0);
}

#endif // ifndef mifit_fios_PDB_h
