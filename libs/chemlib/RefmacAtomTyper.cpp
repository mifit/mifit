#include "RefmacType.h"
#include "RefmacAtomTyper.h"
#include "atom_util.h"
#include <chemlib/RESIDUE_.h>
#include "mol_util.h"

using namespace std;

namespace chemlib
{

RefmacAtomTyper::RefmacAtomTyper(const Residue &res, const vector<Bond> &bonds)
    : AtomTyper(res, bonds)
{
}

char*RefmacAtomTyper::AtomType(const MIAtom *a) const
{
    std::string name = a->name();
    const MIAtom *atom = m_mol.residues[0]->atomByName(name);

    switch (atom->atomicnumber())
    {
    case 1: return TypeHydrogen(atom);
    case 6: return TypeCarbon(atom);
    case 7: return TypeNitrogen(atom);
    case 8: return TypeOxygen(atom);
    case 14: return TypeSilicon(atom);
    case 15: return TypePhosphorus(atom);
    case 16: return TypeSulfur(atom);
    case 32: return TypeGermanium(atom);
    case 33: return TypeArsenic(atom);
    default: return TypeOther(atom);
    }
}

char*RefmacAtomTyper::TypeHydrogen(const MIAtom *atom) const
{
    if (atom->nabors().size())
    {
        return Name(RefmacType::H);
    }

    if (atom->nabors()[0]->atomicnumber() == 16)                //Only one type for
    {
        return Name(RefmacType::HSH1);                          //H's bound to sulfur
    }

    char nabor_type[RefmacType::MAXLENGTH];
    strncpy(nabor_type, AtomType(atom->nabors()[0]), RefmacType::MAXLENGTH);

    switch (Index(nabor_type))
    {
    case RefmacType::CSP1:
        return Name(RefmacType::H);
    case RefmacType::C1:
        return Name(RefmacType::HC1);
    case RefmacType::C2:
        return Name(RefmacType::HC2);
    case RefmacType::CR1H:
        return Name(RefmacType::HCR1);
    case RefmacType::CR15:
        return Name(RefmacType::HCR5);
    case RefmacType::CR16:
        return Name(RefmacType::HCR6);
    case RefmacType::CH1:
        return Name(RefmacType::HCH1);
    case RefmacType::CH2:
        return Name(RefmacType::HCH2);
    case RefmacType::CH3:
        return Name(RefmacType::HCH3);
    case RefmacType::NC1:
        return Name(RefmacType::HNC1);
    case RefmacType::NC2:
        return Name(RefmacType::HNC2);
    case RefmacType::NH1:
        return Name(RefmacType::HNH1);
    case RefmacType::NH2:
        return Name(RefmacType::HNH2);
    case RefmacType::NR15:
        return Name(RefmacType::HNR5);
    case RefmacType::NR16:
        return Name(RefmacType::HNR6);
    case RefmacType::NT1:
        return Name(RefmacType::HNT1);
    case RefmacType::NT2:
        return Name(RefmacType::HNT2);
    case RefmacType::NT3:
        return Name(RefmacType::HNT3);
    case RefmacType::OH1:
        return Name(RefmacType::HOH1);
    case RefmacType::OH2:
        return Name(RefmacType::HOH2);
    case RefmacType::OHA:
        return Name(RefmacType::HOHA);
    case RefmacType::OHB:
        return Name(RefmacType::HOHB);
    case RefmacType::OHC:
        return Name(RefmacType::HOHC);
    default:
        return Name(RefmacType::H);
    }
}

char*RefmacAtomTyper::TypeCarbon(const MIAtom *atom) const
{
    int n_hydrogen = GetNumHydrogens(*atom);
    int naBonds = CountAromaticBonds(*atom, m_mol.bonds);

    switch (atom->hybrid())
    {
    case 1:
        switch (n_hydrogen)
        {
        case 0:
            return Name(RefmacType::CSP);
        case 1:
            return Name(RefmacType::CSP1);
        }

    case 2:

        if (naBonds == 2 && atom->smallest_aromatic_ring() == 6)
        {
            switch (n_hydrogen)
            {
            case 0:
                return Name(RefmacType::CR6);
            case 1:
                return Name(RefmacType::CR16);
            }
        }

        if (naBonds == 2 && atom->smallest_aromatic_ring() == 5)
        {
            switch (n_hydrogen)
            {
            case 0:
                return Name(RefmacType::CR5);
            case 1:
                return Name(RefmacType::CR15);
            }
        }

        if (naBonds == 3 && AtSixSixFusion(*atom, m_mol.bonds))
        {
            return Name(RefmacType::CR66);
        }

        if (naBonds == 3 && AtFiveSixFusion(*atom, m_mol.bonds))
        {
            return Name(RefmacType::CR55);
        }

        if (naBonds == 3 && AtFiveFiveFusion(*atom, m_mol.bonds))
        {
            return Name(RefmacType::CR56);
        }

        switch (n_hydrogen)
        {
        case 0:
            return Name(RefmacType::C);
        case 1:
            return Name(RefmacType::C1);
        case 2:
            return Name(RefmacType::C2);
        }

    case 3:
        switch (n_hydrogen)
        {
        case 0:
            return Name(RefmacType::CT);
        case 1:
            return Name(RefmacType::CH1);
        case 2:
            return Name(RefmacType::CH2);
        case 3:
            return Name(RefmacType::CH3);
        default:
            return Name(RefmacType::CH3);           //Carbon in methane ends up here
        }
    }     //switch(atom->hybrid)

    return Name(RefmacType::C);
}     //TypeCarbon() function

char*RefmacAtomTyper::TypeNitrogen(const MIAtom *atom) const
{
    int n_hydrogen = GetNumHydrogens(*atom);
    int degree = atom->nabors().size() + atom->hcount();

    switch (atom->hybrid())
    {
    case 1:
        return Name(RefmacType::NS);

    case 2:

        if (atom->smallest_aromatic_ring() == 6 && n_hydrogen == 1)
        {
            return Name(RefmacType::NR16);
        }
        if (atom->smallest_aromatic_ring() == 6 && degree == 3)
        {
            return Name(RefmacType::NR6);
        }
        if (atom->smallest_aromatic_ring() == 6 && degree == 2)
        {
            return Name(RefmacType::NRD6);
        }


        if (atom->smallest_aromatic_ring() == 5 && n_hydrogen == 1)
        {
            return Name(RefmacType::NR15);
        }
        if (atom->smallest_aromatic_ring() == 5 && degree == 3)
        {
            return Name(RefmacType::NR5);
        }
        if (atom->smallest_aromatic_ring() == 5 && degree == 2)
        {
            return Name(RefmacType::NRD5);
        }

        switch (n_hydrogen)
        {
        case 0:
            return Name(RefmacType::N);
        case 1:
            return Name(RefmacType::NH1);
        case 2:
            return Name(RefmacType::NH2);
        }

    case 3:
        switch (n_hydrogen)
        {
        case 0:
            return Name(RefmacType::NT);
        case 1:
            return Name(RefmacType::NT1);
        case 2:
            return Name(RefmacType::NT2);
        case 3:
            return Name(RefmacType::NT3);
        default:
            return Name(RefmacType::NT3);           //ammonium ends up here
        }
    }     //switch(atom->hybrid)


    return Name(RefmacType::N);
}     //TypeNitrogen() function

char*RefmacAtomTyper::TypeOxygen(const MIAtom *atom) const
{
    int n_hydrogen = GetNumHydrogens(*atom);
    int hvy_degree = atom->nabors().size() + atom->hcount() - n_hydrogen;       //# of nonhydrogen nabors

    if (hvy_degree == 1 && atom->nabors()[0]->atomicnumber() == 15)
    {
        return Name(RefmacType::OP);
    }
    else if (hvy_degree == 1 && atom->nabors()[0]->atomicnumber() == 16)
    {
        return Name(RefmacType::OS);
    }
    else if (hvy_degree == 1 && atom->nabors()[0]->atomicnumber() == 5)
    {
        return Name(RefmacType::OB);
    }


    switch (atom->hybrid())
    {
    case 2:
        return Name(RefmacType::O);

    case 3:
        switch (n_hydrogen)
        {
        case 0:
            return Name(RefmacType::O2);
        case 1:
            return Name(RefmacType::OH1);
        case 2:
            return Name(RefmacType::OH2);
        }
    }     //switch(atom->hybrid)

    return Name(RefmacType::O);

}     //TypeOxygen() function

char*RefmacAtomTyper::TypeSilicon(const MIAtom *atom) const
{
    switch (atom->hybrid())
    {
    case 3:
        return Name(RefmacType::SI);
    default: return Name(RefmacType::SI1);
    }     //switch(atom->hybrid)
}     //TypeSilicon() function

char*RefmacAtomTyper::TypePhosphorus(const MIAtom *atom) const
{
    switch (atom->hybrid())
    {
    case 3:
        return Name(RefmacType::P);
    default: return Name(RefmacType::P1);
    }     //switch(atom->hybrid)
}     //TypePhosphorus() function

char*RefmacAtomTyper::TypeSulfur(const MIAtom *atom) const
{
    int n_hydrogen = GetNumHydrogens(*atom);

    if (n_hydrogen > 0)                                     // Mercapto/sulfhydryls/thiols
    {
        return Name(RefmacType::SH1);
    }

    int degree = atom->nabors().size();

    if (degree == 1 && atom->hybrid() == 2)               //Thiocarbonyls, isothiocyanates
    {
        return Name(RefmacType::S1);
    }

    if (degree == 2 && atom->hybrid() == 3)               //Sulfides, disulfides
    {
        return Name(RefmacType::S2);
    }

    if (degree == 3)                                        //Sulfoxides, sulfur trioxide
    {
        return Name(RefmacType::S3);
    }

    if (degree == 4 && atom->hybrid() == 3)               //Sulfites, sulfonamides, etc.
    {
        return Name(RefmacType::ST);
    }

    return Name(RefmacType::S);                         //Catch-all for sulfur with no
}                                                       //hydrogen, including sulfates

//I haven't found documentation for this.  Rather it's inferred by analogy
//with silicon.
char*RefmacAtomTyper::TypeGermanium(const MIAtom *atom) const
{
    switch (atom->hybrid())
    {
    case 3:
        return Name(RefmacType::GE);
    default: return Name(RefmacType::GE1);
    }
}

//I haven't found documentation for this.  Rather it's inferred by analogy
//with silicon.
char*RefmacAtomTyper::TypeArsenic(const MIAtom *atom) const
{
    switch (atom->hybrid())
    {
    case 3:
        return Name(RefmacType::AS);
    default: return Name(RefmacType::AS1);
    }
}

char*RefmacAtomTyper::TypeOther(const MIAtom *atom) const
{
    static char per_tab[NELEMENTS+1][RefmacType::MAXLENGTH+1] =
    {
        "DUM",
        "H", "HE",
        "LI", "BE", "B", "C", "N", "O", "F", "NE",
        "NA", "MG", "AL", "SI", "P", "S", "CL", "AR",
        "K", "CA", "SC", "TI", "V", "CR", "MN", "FE", "CO",
        "NI", "CU", "ZN", "GA", "GE", "AS", "SE", "BR", "KR",
        "RB", "SR", "Y", "ZR", "NB", "MO", "TC", "RU", "RH",
        "PD", "AG", "CD", "IN", "SN", "SB", "TE", "I", "XE",
        "CS", "BA", "LA",
        "CE", "PR", "ND", "DUM", "SM", "EU", "GD", "TB", "DY",
        "HO", "ER", "TM", "YB", "LU",
        "HF", "TA", "W", "RE", "OSE", "IR",
        "PT", "AU", "HG", "TL", "PB", "BI", "PO", "AT", "RN",
        "DUM", "DUM", "AC",
        "TH", "PA", "U", "NP",
    };

    return per_tab[atom->atomicnumber()];
}

char*RefmacAtomTyper::Name(unsigned int index) const
{
    static char types[RefmacType::MAXTYPE][RefmacType::MAXLENGTH+1] =
    {
        "DUM",
        "CSP", "CSP1", "C", "C1", "C2", "CR1", "CR2", "CR1H", //1-19 Carbon
        "CR15", "CR5", "CR56", "CR55", "CR16", "CR6", "CR66",
        "CH1", "CH2", "CH3", "CT",
        "NS", "N", "NC1", "NH1", "NC2", "NH2", "NC3", "NT", //20-41 Nitrogen
        "NT1", "NT2", "NT3", "NPA", "NPB", "NR5", "NR15",
        "NRD5", "NR56", "NR55", "NR6", "NR66", "NR16", "NRD6",
        "OS", "O", "O2", "OH1", "OH2", "OHA", "OHB", "OHC", //42-53 Oxygen
        "OC2", "OC", "OP", "OB",
        "P", "P1", "PS",                                     //54-56 Phosphorus
        "S", "S3", "S2", "S1", "ST", "SH1",               //57-62 Sulfur
        "H", "HCH", "HCH1", "HCH2", "HCH3", "HCR1", "HC1", //63-87 Hydrogen
        "HC2", "HCR5", "HCR6", "HNC1", "HNC2", "HNH1", "HNH2",
        "HNR5", "HNR6", "HNT1", "HNT2", "HNT3", "HOH1", "HOH2",
        "HOHA", "HOHB", "HOHC", "HSH1",

        "SI", "SI1", "GE", "GE1", "SN", "PB",         //88-93 Group 14
        "LI", "NA", "K", "RB", "CS",                   //94-98 Group 1
        "BE", "MG", "CA", "SR", "BA",                      //99-103 Group 2
        "SC", "Y", "LA", "CE", "PR", "ND",                //104-109 Group 3
        "SM", "EU", "GD", "TB", "DY", "HO", "ER", "TM", "YB", "LU", //110-119 Lanthanides
        "AC", "TH", "PA", "U", "NP",                   //120-124 Actinides
        "TI", "ZR", "HF",                                    //125-127 Group 4
        "V", "NB", "TA",                                     //128-130 Group 5
        "CR", "MO", "W",                                     //131-133 Group 6
        "MN", "TC", "RE",                                    //134-136 Group 7
        "FE", "RU", "OSE",                                   //137-139 Group 8
        "CO", "RH", "IR",                                    //140-142 Group 9
        "NI", "PD", "PT",                                    //143-145 Group 10
        "CU", "AG", "AU",                                    //146-148 Group 11
        "ZN", "CD", "HG",                                    //149-151 Group 12
        "B", "AL", "GA", "IN", "TL",                   //152-156 Group 13
        "AS", "AS1", "SB", "BI",                            //157-160 Group 15
        "SE", "TE", "PO",                                    //161-163 Group 16
        "F", "CL", "BR", "I", "AT",                        //164-168 Halogens
        "HE", "NE", "AR", "KR", "XE", "RN"                //169-174 Noble Gases
    };

    return types[index];
}

unsigned int RefmacAtomTyper::Index(const char *name) const
{
    static char types[RefmacType::MAXTYPE][RefmacType::MAXLENGTH+1] =
    {
        "DUM",
        "CSP", "CSP1", "C", "C1", "C2", "CR1", "CR2", "CR1H", //1-19 Carbon
        "CR15", "CR5", "CR56", "CR55", "CR16", "CR6", "CR66",
        "CH1", "CH2", "CH3", "CT",
        "NS", "N", "NC1", "NH1", "NC2", "NH2", "NC3", "NT", //20-41 Nitrogen
        "NT1", "NT2", "NT3", "NPA", "NPB", "NR5", "NR15",
        "NRD5", "NR56", "NR55", "NR6", "NR66", "NR16", "NRD6",
        "OS", "O", "O2", "OH1", "OH2", "OHA", "OHB", "OHC", //42-53 Oxygen
        "OC2", "OC", "OP", "OB",
        "P", "P1", "PS",                                     //54-56 Phosphorus
        "S", "S3", "S2", "S1", "ST", "SH1",               //57-62 Sulfur
        "H", "HCH", "HCH1", "HCH2", "HCH3", "HCR1", "HC1", //63-87 Hydrogen
        "HC2", "HCR5", "HCR6", "HNC1", "HNC2", "HNH1", "HNH2",
        "HNR5", "HNR6", "HNT1", "HNT2", "HNT3", "HOH1", "HOH2",
        "HOHA", "HOHB", "HOHC", "HSH1",

        "SI", "SI1", "GE", "GE1", "SN", "PB",         //88-93 Group 14
        "LI", "NA", "K", "RB", "CS"                        //94-98 Group 1
        "BE", "MG", "CA", "SR", "BA",                      //99-103 Group 2
        "SC", "IN", "LA", "CE", "PR", "ND"                //104-109 Group 3
        "SM", "EU", "GD", "TB", "DY", "HO", "ER", "TM", "YB", "LU" //110-119 Lanthanides
        "AC", "TH", "PA", "U", "NP",                   //120-124 Actinides
        "TI", "ZR", "HF",                                    //125-127 Group 4
        "V", "NB", "TA",                                     //128-130 Group 5
        "CR", "MO", "W",                                     //131-133 Group 6
        "MN", "TC", "RE",                                    //134-136 Group 7
        "FE", "RU", "OSE",                                   //137-139 Group 8
        "CO", "RH", "IR",                                    //140-142 Group 9
        "NI", "PD", "PT",                                    //143-145 Group 10
        "CU", "AG", "AU",                                    //146-148 Group 11
        "ZN", "CD", "HG",                                    //149-151 Group 12
        "B", "AL", "GA", "IN", "TL",                   //152-156 Group 13
        "AS", "AS1", "SB", "BI",                            //157-160 Group 15
        "SE", "TE", "PO",                                    //161-163 Group 16
        "F", "CL", "BR", "I", "AT",                        //164-168 Halogens
        "HE", "NE", "AR", "KR", "XE", "RN"                //169-174 Noble Gases
    };

    for (unsigned int i = 0; i < RefmacType::MAXTYPE; ++i)
    {
        if (!strncmp(name, types[i], RefmacType::MAXLENGTH))
        {
            return i;
        }
    }
    return 0;
}

}
