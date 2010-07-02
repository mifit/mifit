#include "SMILES.h"
#include "MIMolIOBase.h"

#include <chemlib/Monomer.h>
#include "FirstToken.h"
#include "SmilesReader.h"

using namespace std;

namespace chemlib
{

SMILES::SMILES()
{
}

SMILES::~SMILES()
{
}

bool SMILES::Read(FILE *fp, MIMolInfo &mol)
{
    std::string error, smiles = MIFirstToken(fp);
    return Read(smiles, mol);
}

bool SMILES::Read(const std::string &smiles, MIMolInfo &mol)
{
    std::string error;

    // clear current mol
    MIMolInfo foo;
    mol = foo;
    Monomer *res = SmilesToMol(smiles, mol.bonds, error);
    if (res)
    {
        res->setSecstr('X');
        res->set_chain_id(' ');
        mol.beginRes = new Residue(*res);
    }
    return res != 0;
}

}
