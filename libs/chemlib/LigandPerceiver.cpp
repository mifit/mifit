#include <vector>
#include <sstream>
#include <fstream>

#include "LigandPerceiver.h"
#include "mol_util.h"
#include "valence.h"
#include "CovalentGeom.h"

namespace chemlib
{

void LigandPerceiver::AssignImpHydrogens(Ligand *lig)
{
    std::vector<Monomer*>::iterator res;
    for (res = lig->residues.begin(); res != lig->residues.end(); ++res)
    {
        AssignImpHydrogens(**res, lig->bonds);
    }
}

void LigandPerceiver::AssignImpHydrogens(Monomer &res, const std::vector<Bond> bonds)
{
    unsigned int i;
    for (i = 0; i < res.atoms().size(); ++i)
    {
        AssignImpHydrogens(*res.atom(i), bonds);
    }
}

void LigandPerceiver::AssignImpHydrogens(MIAtom &atom, const std::vector<Bond> bonds)
{
    if (Is_Onium(atom, bonds))
    {
        return;
    }

    if (atom.atomicnumber() == 16 && atom.isaromatic())         //Don't add to aromatic sulfurs
    {
        atom.setHcount(0);
        return;
    }

    int v = CurrentValence(atom, bonds);

    //Guess for terminal atoms, since we sometimes call this function without
    //bond orders set
    unsigned char order;
    if (atom.nabors().size() == 1)
    {
        order = PredictBondOrder(atom.atomicnumber(),
                                 atom.nabors().front()->atomicnumber(),
                                 (float)AtomDist(atom, *atom.nabors().front()));
        if (order == DOUBLEBOND)
        {
            v = 2;
        }

        if (order == TRIPLEBOND)
        {
            v = 3;
        }
    }

    std::vector<int> pssbl_valences = GetValenceStates(atom.atomicnumber());

    std::vector<int>::iterator i, e;
    i = pssbl_valences.begin();
    e = pssbl_valences.end();

    while (i != e)
    {
        if (v <= *i)
        {
            atom.setHcount(*i - v);
            return;
        }
        ++i;
    }

    //If we reach here, either we haven't defined valence states for the atom (see valence.h)
    //or the atom already has a higher valence than any of the defined states.  In either
    //case, don't give it any implicit hydrogens.

    atom.setHcount(0);
}

void LigandPerceiver::AssignChirality(Ligand *lig)
{
    std::vector<Monomer*>::iterator res;
    for (res = lig->residues.begin(); res != lig->residues.end(); ++res)
    {
        AssignChirality(**res, lig);
    }
}

void LigandPerceiver::AssignChirality(Monomer &res, Ligand*)
{
    unsigned int i;
    for (i = 0; i < res.atoms().size(); ++i)
    {
        if (res.atom(i)->chiral_class() == CH_DEFAULT)
        {
            res.atom(i)->chiral_class(DefaultChiralClass(*res.atom(i)));
        }
    }
}

int LigandPerceiver::DefaultChiralClass(MIAtom &atom)
{
    int degree;
    degree = atom.bondnumbers().size() + atom.hcount();

    switch (degree)
    {
    case 2:
        return CH_ALLENE_LIKE;
    case 4:
        return CH_TETRAHEDRAL;
    case 5:
        return CH_TRIGONAL_BIPYRAMIDAL;
    case 6:
        return CH_OCTAHEDRAL;
    }
    //	cout << "WARNING: input of chirality information for atom with degree " << degree << endl
    //		 << "No known chiral geometry of degree " << degree << endl
    //		 << "Chiral specification ignored" << endl;

    return CH_NONE;
}

/*
   void LigandPerceiver::AssignAromaticBonds(Ligand *lig) {
    std::vector <Bond>::iterator bnd;
    for(bnd=lig->bonds.begin(); bnd!=lig->bonds.end(); ++bnd) {
        bnd->isaromatic = false;								//Set all the bonds to be
    }															//non-aromatic

    lig->FindCycles();

    std::vector <Cycle>::iterator rng;
    std::vector <MIAtom *>::iterator patm;

    bool is_aromatic;
    for(rng=lig->rings.begin(); rng!=lig->rings.end(); ++rng) {
        is_aromatic = true;
        for(patm=rng->atoms.begin(); patm!=rng->atoms.end(); ++patm) {
            if(!(*patm)->isaromatic){
                is_aromatic = false;
                break;
            }
        }
        rng->SetAromatic(is_aromatic);
    }
   }
 */
void LigandPerceiver::AssignHybridization(Ligand *lig)
{

    std::vector<Monomer*>::iterator res;
    for (res = lig->residues.begin(); res != lig->residues.end(); ++res)
    {
        AssignHybridization(**res, lig);
    }
}

void LigandPerceiver::AssignHybridization(Monomer &res, Ligand *lig)
{

    MIAtom *atm;

    unsigned int i;
    for (i = 0; i < res.atoms().size(); ++i)          //First "local" pass, using
    {
        AssignHybridization(res.atom(i), lig); //orders of the atom's own bonds
    }
    for (i = 0; i < res.atoms().size(); ++i)          //Second pass, adjusting the
    {
        atm = res.atom(i);                 //hybridizations for heteroatoms (N,O)
        //next to aromatics
        int degree = atm->nabors().size() + atm->hcount();

        if (atm->atomicnumber() == 6 && atm->formal_charge() > 0 && degree < 4)
        {
            atm->setHybrid(2);
        }
        else if (atm->atomicnumber() == 7 || atm->atomicnumber() == 8)
        {
            AdjustHybridization(atm, lig);
        }
        //Boron as a special case
        else if (atm->atomicnumber() == 5
                 && atm->formal_charge() == 0
                 && atm->nabors().size() < 4
                 && atm->hybrid() > 1)
        {
            atm->setHybrid(atm->hybrid()-1);                   //Adjust boron hybridization
        }
        //Beryllium as a special case
        else if (atm->atomicnumber() == 4 && atm->nabors().size() == 2)
        {
            atm->setHybrid(1);
        }
        else if (atm->atomicnumber() == 4 && atm->nabors().size() == 3)
        {
            atm->setHybrid(2);
        }
        else if (atm->atomicnumber() == 4 && atm->nabors().size() == 4)
        {
            atm->setHybrid(3);
        }
        //Phosphorus as a special case
        else if (atm->atomicnumber() == 15 && atm->nabors().size() == 5)
        {
            atm->setHybrid(4);
        }
        else if (atm->atomicnumber() == 15 && atm->nabors().size() == 4)
        {
            atm->setHybrid(3);
        }
        //Sulfur as a special case
        else if (atm->atomicnumber() == 16 && atm->nabors().size() == 6)
        {
            atm->setHybrid(5);
        }
        else if (atm->atomicnumber() == 16 && atm->nabors().size() == 4 && atm->hybrid() == 3)
        {
            atm->setHybrid(4);
        }
        else if (atm->atomicnumber() == 16 && atm->nabors().size() == 4 && atm->hybrid() < 3)
        {
            atm->setHybrid(3);
        }
        //Chlorine as a special case
        else if (atm->atomicnumber() == 17)
        {
            atm->setHybrid(3);
        }


    }
}

void LigandPerceiver::AssignHybridization(MIAtom *atm, Ligand *lig)
{
    int ndouble = 0;
    std::vector<int>::const_iterator xbond;
    Bond *bnd;

    if (atm->isaromatic())                //Count all aromatics as Sp2 hybridized
    {
        atm->setHybrid(2);
        return;
    }


    for (xbond = atm->bondnumbers().begin(); xbond != atm->bondnumbers().end(); ++xbond)
    {
        bnd = lig->GetBond(*xbond);
        switch (bnd->getOrder())
        {
        case TRIPLEBOND:                //Triple bonds imply Sp hybridization
            atm->setHybrid(1);
            return;
        case PARTIALDOUBLEBOND:
            atm->setHybrid(2);            //Count all aromatics as Sp2 hybridization
            return;                         //(some sulfurs get adjusted later)
        case DOUBLEBOND:
            ndouble++;                  //Count double bonds for rules below
        }
    }

    //If we've gotten this far, the atom is aliphatic and has only single and
    //double bonds.  The atom can have either 2, 1, or 0 double bonds.

    int degree = atm->nabors().size() + atm->hcount();

    switch (ndouble)
    {
    case 3:
        atm->setHybrid(2);                    //Trioxides (perchlorates get fixed)
        return;
    case 2:
        if (atm->atomicnumber() == 16          //Dioxides of sulfur, selenium, & tellurium
            || atm->atomicnumber() == 34      //Other cmpds of these (sulfates, sulfones...)
            || atm->atomicnumber() == 52)     //get fixed in the second pass
        {
            atm->setHybrid(2);
        }
        else if (degree < 3)
        {
            atm->setHybrid(1);                //Allenes, azides, and similar
        }
        else if (degree == 3)
        {
            atm->setHybrid(2);
        }
        else
        {
            atm->setHybrid(3);
        }
        return;
    case 1:
        atm->setHybrid(2);
        return;
    case 0:
        atm->setHybrid(3);
        return;
    default:
        atm->setHybrid(3);                    //4+ double bonds is at least sp3
    }
}

void LigandPerceiver::AdjustHybridization(MIAtom *atm, Ligand *lig)
{

    MIAtom *nabor;
    std::vector<int>::const_iterator xbond;
    Bond *bnd;

    for (xbond = atm->bondnumbers().begin(); xbond != atm->bondnumbers().end(); ++xbond)
    {
        bnd = lig->GetBond(*xbond);
        if (bnd->getAtom1() == atm)
        {
            nabor = bnd->getAtom2();
        }
        else if (bnd->getAtom2() == atm)
        {
            nabor = bnd->getAtom1();
        }

        //		if(nabor->atomicnumber == 6 && nabor->hybrid == 2) {	//Adjust if neighbor atom is an
        //			atm->hybrid = 2;									//sp2-hybridized carbon
        //			return;
        //		}

        if (nabor->hybrid() == 2)
        {
            atm->setHybrid(2);
        }
    }
}

void LigandPerceiver::AssignAtomGeom(Ligand *lig)
{
    std::vector<Monomer*>::iterator res;
    for (res = lig->residues.begin(); res != lig->residues.end(); ++res)
    {
        AssignAtomGeom(**res);
    }
}

void LigandPerceiver::AssignAtomGeom(Monomer &res)
{

    unsigned int i;
    for (i = 0; i < res.atoms().size(); ++i)                  //Just loop through atoms
    {
        AssignAtomGeom(*res.atom(i));
    }
}

void LigandPerceiver::AssignAtomGeom(MIAtom &atm)
{

    if (atm.hybrid() == 1)
    {
        atm.setGeom(LINEAR);
    }
    if (atm.hybrid() == 2)
    {
        atm.setGeom(TRIGONAL_PLANAR);
    }
    if (atm.hybrid() == 3)
    {
        atm.setGeom(TETRAHEDRAL);
    }
    if (atm.hybrid() == 4)                            //Sp3d hybridization
    {
        atm.setGeom(TRIGONAL_BIPYRAMIDAL);
    }
    if (atm.hybrid() == 5)                            //Sp3d2 hybridization
    {
        atm.setGeom(OCTAHEDRAL);
    }
}

bool Is_Onium(MIAtom &atom, std::vector<Bond> bonds)
{
    int d = atom.bondnumbers().size();
    int cv = CurrentValence(atom, bonds);

    if (d != cv)
    {
        return false;
    }

    if (atom.atomicnumber() == 7 && d == 4)               //ammonium
    {
        return true;
    }
    if (atom.atomicnumber() == 8 && d == 3)               //oxonium
    {
        return true;
    }
    if (atom.atomicnumber() == 9 && d == 2)               //fluronium
    {
        return true;
    }
    if (atom.atomicnumber() == 15 && d == 4)              //phosphonium
    {
        return true;
    }
    if (atom.atomicnumber() == 16 && d == 3)              //sulfonium
    {
        return true;
    }
    if (atom.atomicnumber() == 17 && d == 2)              //chloronium
    {
        return true;
    }
    if (atom.atomicnumber() == 33 && d == 4)              //arsonium
    {
        return true;
    }
    if (atom.atomicnumber() == 34 && d == 3)              //selenonium
    {
        return true;
    }
    if (atom.atomicnumber() == 35 && d == 2)              //bromonium
    {
        return true;
    }
    if (atom.atomicnumber() == 51 && d == 4)              //stibonium
    {
        return true;
    }
    if (atom.atomicnumber() == 52 && d == 3)              //telluronium
    {
        return true;
    }
    if (atom.atomicnumber() == 53 && d == 2)              //iodonium
    {
        return true;
    }
    if (atom.atomicnumber() == 83 && d == 4)              //bismuthonium
    {
        return true;
    }

    return false;
}

/*

   void PrepPolarAtom(MIAtom &atom) {

    int d = atom.bondnumbers.size();

    if ((atom.atomicnumber == 7 && d == 4) ||
        (atom.atomicnumber == 8 && d == 3)) {
        atom.hcount = 0;
        return;
    }

    int v = PredictValenceFrom3D(atom);

    std::vector<int> pssbl_valences = GetValenceStates(atom.atomicnumber);

    std::vector<int>::iterator i, e;
    i = pssbl_valences.begin();
    e = pssbl_valences.end();

    while(i != e) {
        if(v <= *i) {
            atom.hcount = *i - v;
            return;
        }
 ++i;
    }

    //If we reach here, either we haven't defined valence states for the atom (see valence.h)
    //or the atom already has a higher valence than any of the defined states.  In either
    //case, don't give it any implicit hydrogens.

    atom.hcount = 0;
    return;

   }

   void PrepPolarAtoms(Residue &res) {
    for_each(res.atoms().begin(), res.atoms().end(), PrepAtom);
   }
 */
} //namespace chemlib
