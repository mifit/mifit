#include <vector>
#include <sstream>
#include <fstream>

#include <math/mathlib.h>
#include <chemlib/chemlib.h>
#include <conflib/conflib.h>
#include <util/utillib.h>

#include "HBond.h"


using namespace chemlib;
using namespace moldraw;

namespace HBondLookup
{

#define N_KNOWN_RESIDUES 24

bool IsKnown(const Residue &res)
{

    static char res_glossary[N_KNOWN_RESIDUES][MAXNAME] =
    {
        "ALA", "CYS", "ASP", "GLU", "PHE", "GLY", "HIS", "ILE", "LYS", "LEU",
        "MET", "ASN", "ASX", "PRO", "GLN", "GLX", "ARG", "SER", "THR", "VAL",
        "TRP", "TYR", "WAT", "HOH"
    };
    for (int i = 0; i < N_KNOWN_RESIDUES; ++i)
    {
        if (res.type()[0] == res_glossary[i][0]                   //Compare the first three chars
            && res.type()[1] == res_glossary[i][1]
            && res.type()[2] == res_glossary[i][2])
        {
            return true;
        }
    }

    return false;
}

bool IsKnownDonor(const MIAtom &atom, const Residue &res)
{
    static char res_glossary[N_KNOWN_RESIDUES][MAXNAME] =
    {
        "ALA", "CYS", "ASP", "GLU", "PHE", "GLY", "HIS", "ILE", "LYS", "LEU",
        "MET", "ASN", "ASX", "PRO", "GLN", "GLX", "ARG", "SER", "THR", "VAL",
        "TRP", "TYR", "WAT", "HOH"
    };
    int i;
    for (i = 0; i < N_KNOWN_RESIDUES; ++i)
    {
        if (res.type()[0] == res_glossary[i][0]                   //Compare the first three chars
            && res.type()[1] == res_glossary[i][1]
            && res.type()[2] == res_glossary[i][2])
        {
            break;
        }
    }

    if (i >= N_KNOWN_RESIDUES)
    {
        return false;
    }

    if (strcmp(atom.name(), "N") == 0
        && i < 22)
    {
        return true;
    }
    if (strcmp(atom.name(), "ND1") == 0          //histidine
        && i == 6)
    {
        return true;
    }
    if (strcmp(atom.name(), "NE2") == 0
        && i == 6)
    {
        return true;
    }
    if (strcmp(atom.name(), "NZ") == 0
        && i == 8)
    {
        return true;
    }
    if (strcmp(atom.name(), "ND2") == 0
        && (i == 11 || i == 12))
    {
        return true;
    }
    if (strcmp(atom.name(), "NE2") == 0
        && (i == 14 || i == 15))
    {
        return true;
    }
    if (strcmp(atom.name(), "NH1") == 0
        && i == 16)
    {
        return true;
    }
    if (strcmp(atom.name(), "NH2") == 0
        && i == 16)
    {
        return true;
    }
    if (strcmp(atom.name(), "OG") == 0
        && i == 17)
    {
        return true;
    }
    if (strcmp(atom.name(), "OG1") == 0
        && i == 18)
    {
        return true;
    }
    if (strcmp(atom.name(), "NE1") == 0
        && i == 20)
    {
        return true;
    }
    if (strcmp(atom.name(), "OH") == 0
        && i == 21)
    {
        return true;
    }
    if (strcmp(atom.name(), "O") == 0
        && (i == 22 || i == 23))
    {
        return true;
    }

    return false;
}

bool IsKnownAcceptor(const MIAtom &atom, const Residue &res)
{
    static char res_glossary[N_KNOWN_RESIDUES][MAXNAME] =
    {
        "ALA", "CYS", "ASP", "GLU", "PHE", "GLY", "HIS", "ILE", "LYS", "LEU",
        "MET", "ASN", "ASX", "PRO", "GLN", "GLX", "ARG", "SER", "THR", "VAL",
        "TRP", "TYR", "WAT", "HOH"
    };
    int i;
    for (i = 0; i < N_KNOWN_RESIDUES; ++i)
    {
        if (res.type()[0] == res_glossary[i][0]                   //Compare the first three chars
            && res.type()[1] == res_glossary[i][1]
            && res.type()[2] == res_glossary[i][2])
        {
            break;
        }
    }

    if (i >= N_KNOWN_RESIDUES)
    {
        return false;
    }

    if (strcmp(atom.name(), "O") == 0                //Peptide backbone
        && i < 22)
    {
        return true;
    }
    if (strcmp(atom.name(), "OT1") == 0              //Peptide backbone
        && i < 22)
    {
        return true;
    }
    if (strcmp(atom.name(), "OT2") == 0              //Peptide backbone
        && i < 22)
    {
        return true;
    }
    if (strcmp(atom.name(), "OTX") == 0              //Peptide backbone
        && i < 22)
    {
        return true;
    }
    if (strcmp(atom.name(), "OD1") == 0
        && (i == 2 || i == 11 || i == 12))
    {
        return true;
    }
    if (strcmp(atom.name(), "OD2") == 0
        && (i == 2 || i == 11 || i == 12))
    {
        return true;
    }
    if (strcmp(atom.name(), "OE1") == 0
        && (i == 3 || i == 14 || i == 15))
    {
        return true;
    }
    if (strcmp(atom.name(), "OE2") == 0
        && (i == 3 || i == 14 || i == 15))
    {
        return true;
    }
    if (strcmp(atom.name(), "ND1") == 0
        && i == 6)
    {
        return true;
    }
    if (strcmp(atom.name(), "NE2") == 0
        && i == 6)
    {
        return true;
    }
    if (strcmp(atom.name(), "OG") == 0
        && i == 17)
    {
        return true;
    }
    if (strcmp(atom.name(), "OG1") == 0
        && i == 18)
    {
        return true;
    }
    if (strcmp(atom.name(), "OH") == 0
        && i == 21)
    {
        return true;
    }
    if (strcmp(atom.name(), "O") == 0
        && (i == 22 || i == 23))
    {
        return true;
    }

    return false;
}

} //end namespace HBondLookup


/////////////////////////////////////////////////////////////////////////////
// Function:    Print
// Purpose:		Reports summary information for the hydrogen bond
// Input:       None
// Output:      Summary returned as a std::string
// Requires:
/////////////////////////////////////////////////////////////////////////////

std::string HBond::Print() const
{

    std::string s;

    s += donor->name();
    s += " ";
    //	s += donor_res->type();
    //	s += donor_res->name();
    //	s += donor_res->chain_id();
    s += acceptor->name();
    s += " ";
    //	s += acceptor_res->type();
    //	s += acceptor_res->name();
    //	s += acceptor_res->chain_id();
    s += "  ";

    std::string dist;
    dist = format("%5.2f", distance);
    s += dist;
    s += "\n";

    return s;
}

namespace moldraw
{

/////////////////////////////////////////////////////////////////////////////
// Function:    GetDonors
// Purpose:		Searches a sequence of residues for atoms capable of donating
//				a hydrogen bond
// Input:       A vector of residues to search
// Output:      A vector (passed by reference) of atoms that meet the criteria of a donor
// Requires:	The function IsDonor()
/////////////////////////////////////////////////////////////////////////////
void GetDonors(const std::vector <chemlib::Residue*> &residues,
               std::vector <chemlib::MIAtom*> &donor_atoms)
{

    std::vector<Residue*>::const_iterator res, resEnd = residues.end();
    MIAtom_const_iter atom, endAtom;
    for (res = residues.begin(); res != resEnd; ++res)
    {
        Residue &r = **res;
        const MIAtomList &atoms = r.atoms();
        endAtom = atoms.end();
        for (atom = atoms.begin(); atom != endAtom; ++atom)
        {
            if (IsDonor(**atom, r))
            {
                donor_atoms.push_back(*atom);
            }
        }
    }
}

/////////////////////////////////////////////////////////////////////////////
// Function:    GetAcceptors
// Purpose:		Searches a sequence of residues for atoms capable of accepting
//				a hydrogen bond
// Input:       A vector of residues to search
// Output:      A vector (passed by reference) of atoms that meet the criteria of an acceptor
// Requires:	The function IsAcceptor()
/////////////////////////////////////////////////////////////////////////////
void GetAcceptors(const std::vector <Residue*> &residues,
                  std::vector <MIAtom*> &acceptor_atoms)
{

    std::vector<Residue*>::const_iterator res, resEnd = residues.end();
    MIAtom_const_iter atom, endAtom;
    for (res = residues.begin(); res != resEnd; ++res)
    {
        Residue &r = **res;
        const MIAtomList &atoms = r.atoms();
        endAtom = atoms.end();
        for (atom = atoms.begin(); atom != endAtom; ++atom)
        {
            if (IsAcceptor(**atom, r))
            {
                acceptor_atoms.push_back(*atom);
            }
        }
    }
}

/////////////////////////////////////////////////////////////////////////////
// Function:	CheckHBond
// Purpose:		Checks whether a pair of atoms meet the geometric criteria
//				for forming a hydrogen bond
// Input:       An ordered pair of atom pointers, and an HBOND struct for writing
// Output:      Records the donor, acceptor, and distance of the hydrogen bond in a struct
//				Returns true if the atoms meet the geometric criteria
// Requires:	MAX_HBOND_LENGTH to be defined
/////////////////////////////////////////////////////////////////////////////
bool CheckHBond(MIAtom *donor, MIAtom *acceptor, HBond &hbond)
{

    float dist = AtomDist(*donor, *acceptor);                 //Measure the length
    if (dist > MAX_HBOND_LENGTH || dist < MIN_HBOND_LENGTH)
    {
        return false;
    }

    float h_angle;
    float h_a_dist;
    float bondlength = chemlib::IdealBondLength(donor->atomicnumber(), 1, 1);

    MIAtomList hpositions;
    if (HeavyDegree(donor) == 0)                //water or ammonium
    {
        h_angle = 180.0F;
        h_a_dist = dist - bondlength;
    }
    else if (HeavyDegree(donor) == 1)
    {
        const MIAtom *donor_nabor = MIAtom::GetHeavyNeighbor(donor);
        h_a_dist = CalcBestHDistance(donor, acceptor, donor_nabor, bondlength);
        h_angle = RAD2DEG * TriangleAngle(dist, h_a_dist, bondlength);
    }
    else
    {
        PositionHydrogens(donor, hpositions);
        if (hpositions.size() == 0)
        {
            return false;
        }
        float tmp_dist, min_dist = 1000.0F;
        for (unsigned int i = 0; i < hpositions.size(); ++i)
        {
            if ((tmp_dist = AtomDist(*hpositions[i], *acceptor)) < min_dist)
            {
                h_angle = CalcAtomAngle(*donor, *hpositions[i], *acceptor);
                min_dist = tmp_dist;
            }
        }
        h_a_dist = min_dist;
    }

    if (h_a_dist > 2.7F)
    {
        return false;
    }
    if (h_angle < 90.0F)
    {
        return false;
    }

    float a_angle;
    if (HeavyDegree(acceptor) == 1)
    {
        const MIAtom *aa = MIAtom::GetHeavyNeighbor(acceptor);
        a_angle = CalcAtomAngle(*donor, *acceptor, *aa);

        if (a_angle < 90.0F)
        {
            return false;
        }
    }

    hbond.distance = dist;
    hbond.donor = donor;                                    //store the atom ptrs and return true
    hbond.acceptor = acceptor;

    for (unsigned int i = 0; i < hpositions.size(); ++i)
    {
        delete hpositions[i];
    }
    return true;
}

} //namespace moldraw
