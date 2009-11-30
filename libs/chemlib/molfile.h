#ifndef mifit_fios_molfile_h
#define mifit_fios_molfile_h

#include <cstdio>
#include <string>
#include <map>
#include <vector>

#include "MIMolIOBase.h"

namespace chemlib
{
    class MIAtom;


    class molfile : public Reader, public Writer
    {

        bool SquirtAtoms(FILE *fp,  const RESIDUE &res, const std::vector<Bond> &bonds,
                         const std::vector<CHIRAL> &chirals);
        bool SquirtBonds(FILE *fp, const RESIDUE &res, const std::vector<Bond> &bonds);
        bool SlurpAtoms(FILE *fp, RESIDUE *res, int natoms);
        bool SlurpBonds(FILE *fp, std::vector<Bond> &bonds, const RESIDUE *res, int nbonds);
        bool SlurpProperties(FILE *fp, RESIDUE *res);
        bool SlurpChirals(RESIDUE *res, std::vector<CHIRALDICT> &chirals);
        const CHIRAL *SearchChirals(const MIAtom *atom, const std::vector<CHIRAL> &chirals);

        int ChargeCode(int charge);
        int ChiralCode(const MIAtom *patom, const RESIDUE &res, const std::vector<Bond> &bonds);
        int BondCode(unsigned char order);
        int ProcessCharge(int charge_code);
        void ProcessChiral(MIAtom *patom, int chiral_code);
        unsigned char ProcessOrder(int order_code);
        char ProcessStereo(int stereo_code, int bond_order);

    public:
        molfile();
        virtual ~molfile();
        virtual std::string getDescription()
        {
            return "MDL molfile (*.mol)";
        }

        virtual std::string getExtension()
        {
            return "*.mol";
        }

        virtual bool Write(FILE *fp, MIMolInfo &mol);
        virtual bool Read(FILE *fp, MIMolInfo &mol);
    };

    class MolPropertyLine
    {
    private:
        std::string _line;
        size_t _cur_pos;
        size_t _length;
    public:
        MolPropertyLine(const char *input);
        bool GetEntry(int &xAtom, int &value);
        int  NumEntries();
    };

}

#endif // ifndef mifit_fios_molfile_h
