#include <vector>
#include <sstream>
#include <fstream>

#include "RingSystem.h"
#include "Ligand.h"
#include "Dictionary.h"
#include "sequence_util.h"
#include "mol_util.h"

#include <list>

namespace chemlib
{

/////////////////////////////////////////////////////////////////////////////
// Function:    Clear
// Purpose:		Remove all data from the object
// Input:       None
// Output:      None
// Requires:
/////////////////////////////////////////////////////////////////////////////
void RingSystem::Clear()
{
    _atoms.clear();
    _bonds.clear();
    _aromatics.clear();

    int i, j;

    for (i = 0; i < _conn_table.num_rows(); ++i)
    {
        for (j = 0; j < _conn_table.num_cols(); ++j)
        {
            _conn_table[i][j] = false;
        }
    }

}

/////////////////////////////////////////////////////////////////////////////
// Function:    SetIndex
// Purpose:		Gives this ring system a number, its place in the vector of ring systems
// Input:       Integer value of the index of the ring system
// Output:      None
// Requires:
/////////////////////////////////////////////////////////////////////////////
void RingSystem::SetIndex(int index)
{
    _rsnumber = index;
}

/////////////////////////////////////////////////////////////////////////////
// Function:    SetMolecule
// Purpose:		Stores a pointer to the parent molecule of this ring system
// Input:       ptr to the parent mol
// Output:      None
// Requires:
/////////////////////////////////////////////////////////////////////////////
void RingSystem::SetMolecule(Ligand *mol)
{
    _lig = mol;
}

/////////////////////////////////////////////////////////////////////////////
// Function:    AddBond
// Purpose:		Include an existing bond in this ring system
// Input:       Pointer to the bond
// Output:      None
// Requires:	The bond is already present in the parent molecule
/////////////////////////////////////////////////////////////////////////////
void RingSystem::AddBond(Bond *bond)
{
    if (std::find(_bonds.begin(), _bonds.end(), bond) == _bonds.end())    //Don't duplicate
    {
        _bonds.push_back(bond);
    }
    bond->ring_system = _rsnumber;
}

/////////////////////////////////////////////////////////////////////////////
// Function:    AddAtom
// Purpose:		Include an existing atom in this ring system
// Input:       Pointer to the atom
// Output:      None
// Requires:	The atom is already present in the parent molecule
/////////////////////////////////////////////////////////////////////////////
void RingSystem::AddAtom(MIAtom *atom)
{
    if (std::find(_atoms.begin(), _atoms.end(), atom) == _atoms.end())    //Don't duplicate
    {
        _atoms.push_back(atom);
    }
    atom->set_ring_system(_rsnumber);
}

/////////////////////////////////////////////////////////////////////////////
// Function:    Contains
// Purpose:		Test whether an atom is in this ring system
// Input:       Pointer to the atom
// Output:      True if the atom is in this ring system, false otherwise
// Requires:
/////////////////////////////////////////////////////////////////////////////
bool RingSystem::Contains(const MIAtom *query) const
{
    return std::find(_atoms.begin(), _atoms.end(), query) != _atoms.end();
}

/////////////////////////////////////////////////////////////////////////////
// Function:    Contains
// Purpose:		Test whether a given bond is in this aromatic system
// Input:       Pointer to the bond
// Output:      True if the atom is in this ring system, false otherwise
// Requires:
/////////////////////////////////////////////////////////////////////////////
bool RingSystem::Contains(const Bond *query) const
{
    return std::find(_bonds.begin(), _bonds.end(), query) != _bonds.end();
}

//////////////////////////////////////////////////////////////////////////////
// Function:    GetAtomIndex
// Purpose:		Get the index of an atom used in the internal storage of the
//				ring system object (i.e. _atoms, _conn_table, and the local variables
//				in the ring-finding & search functions.)
// Input:       Pointer to an atom contained in the aromatic
// Output:      Position of the atom in the _atoms vector
// Requires:
/////////////////////////////////////////////////////////////////////////////
int RingSystem::GetAtomIndex(const MIAtom *query) const
{

    MIAtom_const_iter atom_location;

    atom_location = std::find(_atoms.begin(), _atoms.end(), query);

    if (atom_location == _atoms.end())
    {
        return -1;
    }
    else
    {
        return atom_location - _atoms.begin();
    }
}

//////////////////////////////////////////////////////////////////////////////
// Function:    IsAromatic
// Purpose:		Report whether the ring system contains any aromatic systems
// Input:       None
// Output:      False only if the ring system has no aromatics
// Requires:	The method DetectAromatics() has been run for the ring system
/////////////////////////////////////////////////////////////////////////////
bool RingSystem::IsAromatic() const
{
    if (_aromatics.empty())
    {
        return true;
    }
    else
    {
        return false;
    }
}

//////////////////////////////////////////////////////////////////////////////
// Function:    IsAllAromatic
// Purpose:		Report whether the ring system is entirely aromatic
// Input:       None
// Output:      True or false
// Requires:	The method DetectAromatics() has been run for the ring system
/////////////////////////////////////////////////////////////////////////////
bool RingSystem::IsAllAromatic() const
{
    if (_aromatics.size() == 1 && _aromatics.front().NumAtoms() == (int)_atoms.size())
    {
        return true;
    }
    else
    {
        return false;
    }
}

//////////////////////////////////////////////////////////////////////////////
// Function:    SetAllAromatic
// Purpose:		Set the ring system to be entirely aromatic, including setting
//				the bond orders and atom hybrids accordingly
// Input:       None
// Output:      None
// Requires:
/////////////////////////////////////////////////////////////////////////////
void RingSystem::SetAllAromatic()
{
    _aromatics.clear();
    Aromatic arom;

    std::vector<Bond*>::iterator bnd;
    for (bnd = _bonds.begin(); bnd != _bonds.end(); ++bnd)
    {
        (*bnd)->isaromatic = 1;
        (*bnd)->setOrder(PARTIALDOUBLEBOND);                //May need to revisit for pyrroles, furans, etc.
        arom.AddBond(*bnd);
    }

    MIAtom_iter atm;
    for (atm = _atoms.begin(); atm != _atoms.end(); ++atm)
    {
        (*atm)->setIsaromatic(1);
        (*atm)->setHybrid(2);
        arom.AddAtom(*atm);
    }

    arom.GenerateConnTable();

    _aromatics.push_back(arom);
}

/////////////////////////////////////////////////////////////////////////////
// Function:    Extend
// Purpose:		Recursive search for all atoms connected by cyclic bonds
// Input:       Ptr to an atom to use as a starting point
// Output:      Vector of atoms stored in this object
//				Vector of bonds stored in this object
//				Index of this ring system stored in the atom structs
//				Index of this ring system stored in the bond structs
// Requires:
/////////////////////////////////////////////////////////////////////////////
void RingSystem::Extend(MIAtom *atom)
{
    Bond *bond;
    MIAtom *nabor;
    // Loop over bonds to this atom
    for (unsigned int i = 0; i < atom->bondnumbers().size(); ++i)
    {

        bond = _lig->GetBond(atom->bondnumbers()[i]);

        if (!bond->iscyclic)                // Check if bond is cyclic
        {
            continue;
        }

        AddBond(bond);                      // Add bond to obiect

        nabor = _lig->GetNabor(atom, i);

        if (!this->Contains(nabor))
        {
            AddAtom(nabor);                 // Add atom to obiect
            Extend(nabor);                  // Recurse
        }
    }
}

/////////////////////////////////////////////////////////////////////////////
// Function:    AssignAromaticBonds
// Purpose:		Performs depth-first searches for to assign the aromaticity
//				of all bonds in the ring system
// Input:       Looks up connectivity from parent molecule
// Output:      Stores data within the bond objects object ("_isaromatic" flag)
// Requires:	Atoms have been labeled using the "isaromatic" flag
/////////////////////////////////////////////////////////////////////////////
void RingSystem::AssignAromaticBonds()
{
    std::vector<Bond*>::iterator bnd;
    MIAtom_const_iter atm;

    for (bnd = _bonds.begin(); bnd != _bonds.end(); ++bnd)
    {
        (*bnd)->isaromatic = 0;
    }

    std::list<MIAtom*> path1;
    std::list<MIAtom*> path2;
    MIAtom *root;

    for (bnd = _bonds.begin(); bnd != _bonds.end(); ++bnd)
    {
        for (atm = _atoms.begin(); atm != _atoms.end(); ++atm)
        {
            (*atm)->set_search_flag(0);
        }
        path1.clear();
        path2.clear();

        if ((*bnd)->getAtom1()->isaromatic() == 0         //Skip bonds that don't join
            || (*bnd)->getAtom2()->isaromatic() == 0)     //aromatic atoms
        {
            continue;
        }

        path1.push_back((*bnd)->getAtom1());        //Start with just the 2 atoms
        path2.push_back((*bnd)->getAtom2());        //in this bond
        (*bnd)->getAtom1()->set_search_flag(1);
        (*bnd)->getAtom2()->set_search_flag(2);

        //Now perform 2 breadth-first searches simultaneously from
        //both atoms in the bond.  If the two should meet, as detected
        //by the "search_flag", declare this bond cyclic

        while (path1.size() > 0 && path2.size() > 0)
        {
            root = path1.front();
            for (atm = root->nabors().begin(); atm != root->nabors().end(); ++atm)
            {
                if ((*atm)->isaromatic() == 0)
                {
                    continue;
                }
                if (this->Contains(*atm) == false)
                {
                    continue;
                }
                if (*atm == (*bnd)->getAtom2())         //Prevents backtracking across the
                {
                    continue;                   //original bond
                }
                if ((*atm)->search_flag() == 2)
                {
                    (*bnd)->isaromatic = 1;
                    (*bnd)->setOrder(PARTIALDOUBLEBOND);
                }
                else if ((*atm)->search_flag() == 0)
                {
                    (*atm)->set_search_flag(1);
                    path1.push_back(*atm);
                }
            }
            path1.pop_front();

            root = path2.front();
            for (atm = root->nabors().begin(); atm != root->nabors().end(); ++atm)
            {
                if ((*atm)->isaromatic() == 0)
                {
                    continue;
                }
                if (this->Contains(*atm) == false)
                {
                    continue;
                }
                if (*atm == (*bnd)->getAtom1())         //Prevents backtracking across the
                {
                    continue;                   //original bond
                }
                if ((*atm)->search_flag() == 1)
                {
                    (*bnd)->isaromatic = 1;
                    (*bnd)->setOrder(PARTIALDOUBLEBOND);
                }
                else if ((*atm)->search_flag() == 0)
                {
                    (*atm)->set_search_flag(2);
                    path2.push_back(*atm);
                }
            }
            path2.pop_front();

            if ((*bnd)->isaromatic == 1)
            {
                break;
            }
        } //Breadth-first search
    } //Loop over bonds

    for (atm = _atoms.begin(); atm != _atoms.end(); ++atm)
    {
        (*atm)->set_search_flag(0);
    }
}

void RingSystem::AssignAromaticAtoms()
{
    bool allSp2 = true;
    MIAtom_iter atm;
    for (atm = _atoms.begin(); atm != _atoms.end(); ++atm)
    {
        if (!(*atm)->isaromatic()
            && !((*atm)->hybrid() == 2 && (*atm)->formal_charge() == 0))
        {
            allSp2 = false;
            break;
        }
    }
    if (allSp2)
    {
        for (atm = _atoms.begin(); atm != _atoms.end(); ++atm)
        {
            (*atm)->setIsaromatic(1);
        }
    }
}

/////////////////////////////////////////////////////////////////////////////
// Function:    DetectAromatics
// Purpose:		Initiates a search for all aromatic systems within the ring
//				system.
// Input:       Looks up connectivity from parent molecule
// Output:      Stores aromatic system info within the object ("_aromatics")
// Requires:	Atoms have been labeled using the "isaromatic" flag
/////////////////////////////////////////////////////////////////////////////
void RingSystem::DetectAromatics()
{

    AssignAromaticAtoms();
    AssignAromaticBonds();

    Aromatic arom;
    bool *used = new bool[NumAtoms()];
    int xatom;


    for (xatom = 0; xatom < NumAtoms(); ++xatom)             //Initialize the array that tracks
    {
        used[xatom] = false;                            //which candidate atoms are already used
    }

    MIAtom_iter atm;
    for (atm = _atoms.begin(); atm != _atoms.end(); ++atm)
    {

        xatom = GetAtomIndex(*atm);

        if (used[xatom])
        {
            continue;
        }
        if (!(*atm)->isaromatic())
        {
            used[xatom] = true;
            continue;
        }

        arom.Clear();
        //		arom.SetMolecule(_lig);
        arom.AddAtom(*atm);
        used[xatom] = true;
        ExtendAromatic(arom, *atm, used);

        arom.GenerateConnTable();

        _aromatics.push_back(arom);
    }

    delete[] used;
}

/////////////////////////////////////////////////////////////////////////////
// Function:    ExtendAromatics
// Purpose:		Recursive search for all atoms connected by aromatic bonds
// Input:       Ptr to an atom to use as a starting point
//				Array to track which atoms have already been checked
// Output:      Populates an Aromatic object with atoms and bonds
// Requires:	Bonds and atoms have been labeled using the "isaromatic" flag
/////////////////////////////////////////////////////////////////////////////
void RingSystem::ExtendAromatic(Aromatic &arom, MIAtom *atom, bool *used)
{
    // Loop over bonds to this atom
    Bond *bond;
    MIAtom *nabor;
    int xnabor;
    for (unsigned int i = 0; i < atom->bondnumbers().size(); ++i)
    {

        bond = _lig->GetBond(atom->bondnumbers()[i]);     //Get pointer to the bond

        if (!bond->isaromatic
            || !bond->iscyclic)             //Confine search to aromatic bonds
        {
            continue;
        }

        arom.AddBond(bond);                 //Add bond to obiect

        nabor = _lig->GetNabor(atom, i);    //Get pointer to the adjoining atom
        xnabor = GetAtomIndex(nabor);       //Get index of the atom within the ringsys

        if (used[xnabor])
        {
            continue;
        }

        arom.AddAtom(nabor);                //Add atom to obiect
        used[xnabor] = true;
        ExtendAromatic(arom, nabor, used);  //Recurse
    }
}

/////////////////////////////////////////////////////////////////////////////
// Function:    PrintAromatics
// Purpose:		Loops through the aromatics in this ring system, printing
//				a summary of the atoms and bonds in each.
// Input:       None
// Output:      Appends the summaries to a std::string.
// Requires:	The method DetectAromatics() has been run for the ring system
/////////////////////////////////////////////////////////////////////////////
std::string RingSystem::PrintAromatics() const
{
    std::string s;
    std::vector<Aromatic>::const_iterator arom;
    for (arom = _aromatics.begin(); arom != _aromatics.end(); ++arom)
    {
        arom->Print(s);
    }
    return s;
}

/////////////////////////////////////////////////////////////////////////////
// Function:    Print
// Purpose:		Reports summary information for the aromatic system
// Input:       Looked up from members of the Aromatic class
// Output:      Summary appended to a std::string
// Requires:	The method DetectAromatics() has been run for the parent ring system
/////////////////////////////////////////////////////////////////////////////

std::string RingSystem::Print() const
{
    MIAtom_const_iter atm;
    std::vector<Bond*>::const_iterator bond;

    std::string s;
    s += "Printing Ring System...";
    s += "\n";

    for (atm = _atoms.begin(); atm != _atoms.end(); ++atm)
    {
        s += (*atm)->name();
        s += " ";
    }
    s += "\n";

    for (bond = _bonds.begin(); bond != _bonds.end(); ++bond)
    {
        s +=  (*bond)->getAtom1()->name();
        s += "-";
        s += (*bond)->getAtom2()->name();
        s += " ";
    }
    s += "\n";

    return s;
}

/////////////////////////////////////////////////////////////////////////////
// Function:    GeneratePlanes
// Purpose:		Loops through the aromatics in this ring system, generating Planes
//				for each.
// Input:       None
// Output:      Populates the LigDictionary object with Planes.
// Requires:	The method DetectAromatics() has been run for the ring system
/////////////////////////////////////////////////////////////////////////////
void RingSystem::GeneratePlanes(LigDictionary &dict)
{
    std::vector<Aromatic>::const_iterator arom;
    for (arom = _aromatics.begin(); arom != _aromatics.end(); ++arom)
    {
        arom->GeneratePlane(dict);
    }
}

/////////////////////////////////////////////////////////////////////////////
// Function:    GenerateImpropers
// Purpose:		Loops through the aromatics in this ring system, generating
//				all the Improper objects for each aromatic.
// Input:       None
// Output:      Populates the LigDictionary object with Impropers.
// Requires:	The method DetectAromatics() has been run for the ring system
/////////////////////////////////////////////////////////////////////////////
void RingSystem::GenerateImpropers(LigDictionary &dict)
{
    std::vector<Aromatic>::const_iterator arom;
    for (arom = _aromatics.begin(); arom != _aromatics.end(); ++arom)
    {
        arom->GenerateImpropers(dict);
    }
}

/////////////////////////////////////////////////////////////////////////////
// Function:    GenerateConnTable
// Purpose:		Stores a matrix of booleans that maps which pairs of atoms
//				in the ring system are connected.
// Input:       None
// Output:      Stores the matrix in the _conn_table member of this object
// Requires:
/////////////////////////////////////////////////////////////////////////////
void RingSystem::GenerateConnTable()
{
    _conn_table.newsize(NumAtoms(), NumAtoms());
    int i, j;

    for (i = 0; i < _conn_table.num_rows(); ++i)
    {
        for (j = 0; j < _conn_table.num_cols(); ++j)
        {
            _conn_table[i][j] = false;
        }
    }

    int xatom1, xatom2;
    std::vector<Bond*>::iterator bnd;
    for (bnd = _bonds.begin(); bnd != _bonds.end(); ++bnd)
    {
        xatom1 = GetAtomIndex((*bnd)->getAtom1());
        xatom2 = GetAtomIndex((*bnd)->getAtom2());

        if (xatom1 != -1 && xatom2 != -1)
        {
            _conn_table[xatom1][xatom2] = true;
            _conn_table[xatom2][xatom1] = true;
        }
    }
}

/////////////////////////////////////////////////////////////////////////////
// Function:    GenerateSRData
// Purpose:		Initiates calculation of the size of the smallest ring in which
//				each atom and bond in the ring system is contained
// Input:       None
// Output:      Stores the size of the smallest ring in each MIAtom (atom) struct
//				Stores the size of the smallest ring in each Bond (bond) struct
// Requires:	The method GenerateConnTable() has been run for the ring system
/////////////////////////////////////////////////////////////////////////////
void RingSystem::GenerateSRData()
{
    MIAtom_const_iter atm;
    for (atm = _atoms.begin(); atm != _atoms.end(); ++atm)
    {
        SmallestRing(*atm);
    }

    std::vector<Bond*>::const_iterator bnd;
    for (bnd = _bonds.begin(); bnd != _bonds.end(); ++bnd)
    {
        (*bnd)->smallest_ring_size = SmallestRing(*bnd);
    }

    //Set the values of the aromatic ring sizes for each aromatic system in
    //this ring system
    std::for_each(_aromatics.begin(), _aromatics.end(), void_mem_fun_ref(&Aromatic::GenerateSRData));
}

/////////////////////////////////////////////////////////////////////////////
// Function:    DetectFusedRing
// Purpose:		Determines whether the ring system consists of a single ring or
//				multiple fused rings.
// Input:       Looks up connectivity information
// Output:      Stores the size of the smallest ring in each MIAtom (atom) struct
//				Stores the size of the smallest ring in each Bond (bond) struct
// Requires:	The method SetRingFlags() has been run for the parent molecule
/////////////////////////////////////////////////////////////////////////////
void RingSystem::DetectFusedRing()
{
    _fused = false;

    MIAtom_const_iter atm;
    for (atm = _atoms.begin(); atm != _atoms.end(); ++atm)
    {
        if ((*atm)->CountCyclicBonds() > 2)
        {
            _fused = true;
            return;
        }
    }
}

/////////////////////////////////////////////////////////////////////////////
// Function:    GetNormal
// Purpose:     Calclulates a normal vector to the best-fit ring plane
// Input:       Internal lookup of coordinates from parent mol
// Output:      Writes x, y, and z values to the input vector
// Requires:
/////////////////////////////////////////////////////////////////////////////
void RingSystem::GetNormal(std::vector<double> &normal) const
{
    double p_displace;          //not needed, in this case

    normal.clear();             //Initialize the vector to be sure that we
    normal.push_back(0.0);      //have enough space to store the x, y, and z
    normal.push_back(0.0);      //values--and to be sure that &normal[0]
    normal.push_back(0.0);      //produces a valid pointer

    LSqrPlane(&normal[0], &p_displace);
}

/////////////////////////////////////////////////////////////////////////////
// Function:    Flatten
// Purpose:     Projects the atoms of the ring system into 2 dimensions, adjusting
//				the other atoms in the molecule to preserve bond lengths
// Input:       Internal lookup of coordinates from parent mol
// Output:      Writes new coordinates to the atom structs in the parent mol
// Requires:
/////////////////////////////////////////////////////////////////////////////
void RingSystem::Flatten()
{
    double p_normal[3];         //vector normal to best-fit plane
    double p_displace;          //displacement from origin of best-fit plane

    LSqrPlane(p_normal, &p_displace);

    double dev;                 //deviation from the plane of an atom

    std::vector<std::pair<MIAtom*, MIAtom*> > exocycs;
    std::vector<double> original_length;

    GetExocyclics(&exocycs);

    for (unsigned int i = 0; i < exocycs.size(); ++i)
    {
        original_length.push_back(AtomDist(*exocycs[i].first,
                                           *exocycs[i].second) );
    }

    MIAtom_iter atm;
    for (atm = _atoms.begin(); atm != _atoms.end(); ++atm)
    {
        dev = -DotVect((*atm)->x(), (*atm)->y(), (*atm)->z(),
                       p_normal[0], p_normal[1], p_normal[2]);
        dev += p_displace;

        AtomStep(*atm, p_normal, dev);
    }

    for (unsigned int i = 0; i < exocycs.size(); ++i)
    {
        if (exocycs[i].second->nabors().size() == 1)          //Check if this is a terminal
        {
            MoveIntoPlane(exocycs[i].second,                //bond.  If so, move it into
                          p_normal,                         //the same plane with the ring
                          p_displace);                      //atoms
        }
    }

    double bv[3];
    double new_length;

    for (unsigned int i = 0; i < exocycs.size(); ++i)
    {
        BondVector(exocycs[i].second,
                   exocycs[i].first,
                   bv);

        new_length = AtomDist(*exocycs[i].second,
                              *exocycs[i].first);

        ScaleVect(bv, -(original_length[i] - new_length) / new_length);

        _lig->TranslateFragment(exocycs[i].second, exocycs[i].first, bv);
    }
}

/////////////////////////////////////////////////////////////////////////////
// Function:    FixExocyclics
// Purpose:		Adjusts covalent bonds extending from a ring so that the look
//				normal...they extend directly out from the ring.
// Input:       None
// Output:      Updates coordinates
// Requires:
/////////////////////////////////////////////////////////////////////////////
void RingSystem::FixExocyclics()
{

    std::vector<std::pair<MIAtom*, MIAtom*> > exocycs;
    GetCovExocyclics(&exocycs);
    double target[3], blength;

    for (unsigned int i = 0; i < exocycs.size(); ++i)
    {
        if (exocycs[i].second->nabors().size() > 1)
        {
            continue;
        }
        blength = AtomDist(*exocycs[i].first, *exocycs[i].second);
        ExoDirection(exocycs[i].first, target);

        exocycs[i].second->setPosition((float)(exocycs[i].first->x() + target[0] * blength),
                                       (float)(exocycs[i].first->y() + target[1] * blength),
                                       (float)(exocycs[i].first->z() + target[2] * blength));
    }

}

/////////////////////////////////////////////////////////////////////////////
// Function:    GetExocyclics
// Purpose:		Populates a vector with atom pairs representing bond vectors of
//				all the ring-substituent bonds emanating from this ring system
// Input:       Internal lookup of atoms from molecule
// Output:      A vector of atom pairs, each with the substituent atom second
// Requires:
/////////////////////////////////////////////////////////////////////////////
void RingSystem::GetExocyclics(std::vector<std::pair <MIAtom*, MIAtom*> > *exocycs) const
{
    MIAtom *atom;
    std::pair <MIAtom*, MIAtom*> atom_pair;
    MIAtom_const_iter nabor;

    for (unsigned int i = 0; i < _atoms.size(); ++i)
    {
        atom = _atoms[i];

        for (nabor = atom->nabors().begin();
             nabor != atom->nabors().end();
             nabor++)
        {

            if (!this->Contains(*nabor))
            {
                atom_pair.first = atom;
                atom_pair.second = *nabor;
                exocycs->push_back(atom_pair);
            }
        }
    }
}

/////////////////////////////////////////////////////////////////////////////
// Function:    GetCovExocyclics
// Purpose:		Populates a vector with atom pairs representing bond vectors of
//				all the ring-substituent covalent bonds emanating from this ring system
// Input:       Internal lookup of atoms from molecule
// Output:      A vector of atom pairs, each with the substituent atom second
// Requires:
/////////////////////////////////////////////////////////////////////////////
void RingSystem::GetCovExocyclics(std::vector<std::pair <MIAtom*, MIAtom*> > *exocycs)
{
    MIAtom *atom;
    Bond *bond;
    std::pair <MIAtom*, MIAtom*> atom_pair;
    MIAtom_const_iter nabor;

    for (unsigned int i = 0; i < _atoms.size(); ++i)
    {
        atom = _atoms[i];

        for (nabor = atom->nabors().begin();
             nabor != atom->nabors().end();
             nabor++)
        {
            bond = _lig->GetBond(atom, *nabor);
            if (bond->getOrder() > 4)
            {
                continue;
            }

            if (!this->Contains(*nabor))
            {
                atom_pair.first = atom;
                atom_pair.second = *nabor;
                exocycs->push_back(atom_pair);
            }
        }
    }
}

/////////////////////////////////////////////////////////////////////////////
// Function:    ExoDirection
// Purpose:		Computes the vector that will extend directly out from a ring
//				atom.  Specifically, it points opposite the average vector
//				of all the cyclic bonds from the given atom.
// Input:       A ptr to the atom of interest
// Output:      A unit-length vector (in 3-D), written to "v"
// Requires:
/////////////////////////////////////////////////////////////////////////////
void RingSystem::ExoDirection(MIAtom *source, double *v) const
{
    v[0] = 0;
    v[1] = 0;
    v[2] = 0;
    float bond_vect[3];
    int nnabors = 0;

    MIAtom_const_iter nbr;
    for (nbr = source->nabors().begin(); nbr != source->nabors().end(); ++nbr)
    {

        if (!this->Contains(*nbr))
        {
            continue;
        }

        bond_vect[0] = (*nbr)->x() - source->x();
        bond_vect[1] = (*nbr)->y() - source->y();
        bond_vect[2] = (*nbr)->z() - source->z();

        NormVect(bond_vect);

        v[0] -= bond_vect[0];
        v[1] -= bond_vect[1];
        v[2] -= bond_vect[2];

        nnabors++;
    }
    NormVect(v);
}

/////////////////////////////////////////////////////////////////////////////
// Function:    LSqrPlane
// Purpose:		Calculates the best-fit plane to all the atoms in the ring system
// Input:       Internal lookup of coordinates from parent mol
// Output:      A unit vector normal to the best-fit plane
//				A scalar that multiplies the normal vector to extend it from the origin
//				to the best-fit plane (i.e. the signed distance from origin to plane)
// Requires:
/////////////////////////////////////////////////////////////////////////////
void RingSystem::LSqrPlane(double vm[3], double *d) const
{

    double xs[3], xxs[3][3], b[3][3];
    double a[3][3], vmi[3], bv[3];
    double zip = 1.0E-5;
    double orm, vm0, ratio0, ratio1, ratio2, rat01, rat02;
    int i, j, kk, nnn, k, n;

    n = _atoms.size();
    for (i = 0; i < 3; i++)
    {
        xs[i] = 0.0;
    }
    for (k = 0; k < n; k++)
    {
        /* zip is added to prevent numerical instability if
         * atoms in a plane = 0
         */
        xs[0] += _atoms[k]->x()+zip;
        xs[1] += _atoms[k]->y()+zip;
        xs[2] += _atoms[k]->z()+zip;
    }
    for (i = 0; i < 3; i++)
    {
        for (j = 0; j < 3; j++)
        {
            xxs[i][j] = 0.0;
        }
    }

    for (k = 0; k < n; k++)
    {
        xxs[0][0] += _atoms[k]->x() * _atoms[k]->x();
        xxs[0][1] += _atoms[k]->x() * _atoms[k]->y();
        xxs[0][2] += _atoms[k]->x() * _atoms[k]->z();
        xxs[1][0] += _atoms[k]->y() * _atoms[k]->x();
        xxs[1][1] += _atoms[k]->y() * _atoms[k]->y();
        xxs[1][2] += _atoms[k]->y() * _atoms[k]->z();
        xxs[2][0] += _atoms[k]->z() * _atoms[k]->x();
        xxs[2][1] += _atoms[k]->z() * _atoms[k]->y();
        xxs[2][2] += _atoms[k]->z() * _atoms[k]->z();
    }
    for (i = 0; i < 3; i++)
    {
        for (j = 0; j < 3; j++)
        {
            a[i][j] = xxs[i][j] - xs[i]*xs[j]/(float)n;
        }
    }

    /* evaluate matrix */
    b[0][0] = a[1][1] * a[2][2] - a[1][2] * a[2][1];
    b[1][0] = a[2][0] * a[1][2] - a[1][0] * a[2][2];
    b[2][0] = a[1][0] * a[2][1] - a[2][0] * a[1][1];
    b[0][1] = a[2][1] * a[0][2] - a[0][1] * a[2][2];
    b[1][1] = a[0][0] * a[2][2] - a[2][0] * a[0][2];
    b[2][1] = a[2][0] * a[0][1] - a[0][0] * a[2][1];
    b[0][2] = a[0][1] * a[1][2] - a[0][2] * a[1][1];
    b[1][2] = a[1][0] * a[0][2] - a[0][0] * a[1][2];
    b[2][2] = a[0][0] * a[1][1] - a[1][0] * a[0][1];

    /* choose the largest column vector of b as initial solution */
    bv[0] = b[0][0]*b[0][0] + b[1][0]*b[1][0] + b[2][0]*b[2][0];
    bv[1] = b[0][1]*b[0][1] + b[1][1]*b[1][1] + b[2][1]*b[2][1];
    bv[2] = b[0][0]*b[0][2] + b[1][2]*b[1][2] + b[2][2]*b[2][2];
    kk = 0;
    if (bv[1] > bv[0])
    {
        kk = 1;
    }
    if (bv[2] > bv[kk])
    {
        kk = 2;
    }
    vm0 = b[0][kk];
    for (i = 0; i < 3; i++)
    {
        vmi[i] = b[i][kk]/vm0;
    }
    /* solve to convergence by iteration of M(I)=B*M(I-1) */
    for (nnn = 0; nnn < 10; nnn++)
    {
        vm[0] = b[0][0]*vmi[0]+b[0][1]*vmi[1]+b[0][2]*vmi[2];
        vm[1] = b[1][0]*vmi[0]+b[1][1]*vmi[1]+b[1][2]*vmi[2];
        vm[2] = b[2][0]*vmi[0]+b[2][1]*vmi[1]+b[2][2]*vmi[2];
        ratio0 = vm[0]/vmi[0];
        ratio1 = vm[1]/vmi[1];
        ratio2 = vm[2]/vmi[2];
        rat01 = fabs(ratio1/ratio0-1.0);
        rat02 = fabs(ratio2/ratio0-1.0);
        if (rat01 < zip && rat02 < zip)
        {
            break;
        }
        else
        {
            for (i = 0; i < 3; i++)
            {
                vmi[i] = vm[i]/vm[0];
            }
        }
    } /* 100*/
    orm = 0.0;
    /* normalize the solution vectors */
    for (i = 0; i < 3; i++)
    {
        orm += vm[i]*vm[i];
    }
    orm = sqrt(orm);
    for (i = 0; i < 3; i++)
    {
        vm[i] /= orm;
    }
    *d = (vm[0]*xs[0]+vm[1]*xs[1]+vm[2]*xs[2])/(float)n;
}

/////////////////////////////////////////////////////////////////////////////
// Function:    Smallest Ring
// Purpose:		Initiates a calculation of the smallest ring in which the given atom
//				is contained.
// Input:       Ptr to the atom to use as the starting point
// Output:      Returns the size of the smallest ring
// Requires:	The method GenerateConnTable() has been run for the ring system
/////////////////////////////////////////////////////////////////////////////
int RingSystem::SmallestRing(const MIAtom *atom)
{
    int xatom;
    if ((xatom = GetAtomIndex(atom)) == -1)         //Check that this atom is in the ring
    {
        return -1;
    }

    if (!_fused)                       //If the ring system consists of only one ring
    {
        return _atoms.size();           //return the size of that ring
    }


    int *path = new int[_atoms.size()];     //Get space for the arrays used in the
    bool *used = new bool[_atoms.size()];   //search
    InitializeArray(used, _atoms.size(), false);
    int min = _atoms.size() + 1;

    path[0] = xatom;                        //Initialize the search arrays, starting
    used[xatom] = true;                     //with the input atom

    ExtendPath(path, used, 1, min);         //Perform the recursive, depth-first search

    delete[] path;
    delete[] used;

    return min;
}

/////////////////////////////////////////////////////////////////////////////
// Function:    Smallest Ring
// Purpose:		Initiates a calculation of the smallest ring in which the given bond
//				is contained.
// Input:       Ptr to the bond to use as the starting point
// Output:      Returns the size of the smallest ring
// Requires:	The method GenerateConnTable() has been run for the ring system
/////////////////////////////////////////////////////////////////////////////
int RingSystem::SmallestRing(const Bond *bond)
{
    int xatom1, xatom2;

    if ((xatom1 = GetAtomIndex(bond->getAtom1())) == -1)    //Check that the 1st atom is in the ringsys
    {
        return -1;
    }
    if ((xatom2 = GetAtomIndex(bond->getAtom2())) == -1)    //Check that the 2nd atom is in the ringsys
    {
        return -1;
    }

    if (!_fused)                       //If the ring system consists of only one ring
    {
        return _atoms.size();           //return the size of that ring
    }


    int *path = new int[_atoms.size()];     //Get space for the arrays used in the
    bool *used = new bool[_atoms.size()];   //search
    InitializeArray(used, _atoms.size(), false);
    int min = _atoms.size();

    path[0] = xatom1;
    path[1] = xatom2;                       //Initialize the search arrays, starting
    used[xatom1] = true;                    //with the input atom
    used[xatom2] = true;

    ExtendPath(path, used, 2, min);         //Perform the recursive, depth-first search

    delete[] path;
    delete[] used;

    return min;
}

/////////////////////////////////////////////////////////////////////////////
// Function:    Smallest Ring
// Purpose:		Initiates a calculation of the smallest ring in which the given
//				triplet of atoms is contained.
// Input:       Ptr to a triplet of atoms to use as the starting point
// Output:      Returns the size of the smallest ring
// Requires:	The method GenerateConnTable() has been run for the ring system
//				That the atoms be in order (i.e. atom2 is bonded to atom1 and atom3)
/////////////////////////////////////////////////////////////////////////////
int RingSystem::SmallestRing(const MIAtom *atom1,
                             const MIAtom *atom2,
                             const MIAtom *atom3)
{
    int xatom1, xatom2, xatom3;

    if ((xatom1 = GetAtomIndex(atom1)) == -1)       //Check that the 1st atom is in the ringsys
    {
        return -1;
    }
    if ((xatom2 = GetAtomIndex(atom2)) == -1)       //Check that the 2nd atom is in the ringsys
    {
        return -1;
    }
    if ((xatom3 = GetAtomIndex(atom3)) == -1)       //Check that the 3rd atom is in the ringsys
    {
        return -1;
    }

    if (!_fused)                       //If the ring system consists of only one ring
    {
        return _atoms.size();           //return the size of that ring
    }


    int *path = new int[_atoms.size()];     //Get space for the arrays used in the
    bool *used = new bool[_atoms.size()];   //search
    InitializeArray(used, _atoms.size(), false);
    int min = _atoms.size();

    path[0] = xatom1;                       //Assumes the three atoms are in sequence
    path[1] = xatom2;
    path[2] = xatom3;                       //Initialize the search arrays, starting
    used[xatom1] = true;                    //with the input atoms
    used[xatom2] = true;
    used[xatom3] = true;

    ExtendPath(path, used, 3, min);         //Perform the recursive, depth-first search

    delete[] path;
    delete[] used;

    return min;
}

/////////////////////////////////////////////////////////////////////////////
// Function:    ExtendPath
// Purpose:		Recursively determines the smallest ring that can be formed from
//				a given set of atoms.
// Input:       Array "path" is a ordered sequence of atom indices
//				Array "used" is a set of booleans which tracks which atoms have been used
//				"depth" is the number of atoms in the path
//				"min" is the smallest ring identified so far
// Output:      Replaces input value of "min" with the size of the smallest ring
// Requires:	The method GenerateConnTable() has been run for the ring system
/////////////////////////////////////////////////////////////////////////////
void RingSystem::ExtendPath(int *path,          //Ordered array of atom indices in the path
                            bool *used,         //Array that tracks which atoms have been used
                            int depth,          //Number of atoms in current path
                            int &min) const      //Size of smallest ring found so far
{
    int i;
    for (i = 0; i < (int) _atoms.size(); ++i)           //Loop thru atoms in ringsys
    {
        if (_conn_table[path[depth-1]][i] == false)
        {
            continue;                           //Skip if this atom not connected to last
        }

        if (i == path[0] && depth > 2)          //Check if this completes a cycle
        {
            min = depth;
        }
        else if (used[i] == false && depth + 1 < min)
        {
            path[depth] = i;
            used[i] = true;
            ExtendPath(path, used, depth + 1, min);
            used[i] = false;                    //Reset for next bond
        }
    }                                           //End loop over atoms
}

/////////////////////////////////////////////////////////////////////////////
// Function:    CyclohexImpropers
// Purpose:
// Input:       A reference to the dictionary for the parent molecule
// Output:      Stores the Improper objects in the dictionary
// Requires:	The method KeySmallestRings() has been run for this RingSystem object
//				The method GenerateConnTable() has been run for this RingSystem object
/////////////////////////////////////////////////////////////////////////////

void RingSystem::CyclohexImpropers(LigDictionary &dict) const
{
    std::vector<Bond*>::const_iterator bnd;
    std::vector< MIAtomList > torsions;
    std::vector< MIAtomList >::iterator atms;

    std::vector<double> angles;                 //Set the ideal angles to those seen in a
    angles.push_back(58.0);                     //"chair"-shaped conformation
    angles.push_back(-58.0);

    Improper imp;
    int min_size;


    for (bnd = _bonds.begin(); bnd != _bonds.end(); ++bnd)
    {
        if ((*bnd)->getAtom1()->hybrid() != 3
            || (*bnd)->getAtom2()->hybrid() != 3)
        {
            continue;
        }
        if ((*bnd)->smallest_ring_size != 6)
        {
            continue;
        }

        if (((*bnd)->getAtom1())->CountCyclicBonds() > 2            //Skip bonds at ring fusions
            && ((*bnd)->getAtom2())->CountCyclicBonds() > 2)
        {
            continue;
        }

        min_size = _atoms.size() + 1;
        torsions.clear();
        EnumerateTorsions(*bnd, torsions);
        for (atms = torsions.begin(); atms != torsions.end(); ++atms)
        {

            if (Contains((*atms)[0])                        //Check if the end atoms are in
                && Contains((*atms)[3])                     //the ring system,
                && SmallestAliphaticRing((*atms)[0],
                                         (*atms)[1],        //Calc whether this torsions is in
                                         (*atms)[2],        //cyclohexane-type ring
                                         (*atms)[3]) == 6)
            {

                imp.ReInit(*atms, angles);
                dict.AddImproper(imp, false);
            }
        }                                                   //End loop over torsions around bond
    }                                                       //End loop over bonds in ring system
}

/////////////////////////////////////////////////////////////////////////////
// Function:    SmallestAliphaticRing
// Input:       Pointers to four consecutive atoms (in order)
// Output:      Returns the size (in atoms) of the smallest aliphatic ring
//				containing the atoms
// Requires:	The input args be ordered correctly.  That is, atom1 must be bonded
//				to atom2 and atom2 must be bonded to atom3, etc.
//				The method DetectAromatics() has been run for the parent ring system
//				The method GenerateConnTable() has been run for this object
/////////////////////////////////////////////////////////////////////////////

int RingSystem::SmallestAliphaticRing(const MIAtom *atom1,
                                      const MIAtom *atom2,
                                      const MIAtom *atom3,
                                      const MIAtom *atom4) const
{

    int xatom1, xatom2, xatom3, xatom4;

    if ((xatom1 = GetAtomIndex(atom1)) == -1)       //Check that the 1st atom is in the ringsys
    {
        return -1;
    }
    if ((xatom2 = GetAtomIndex(atom2)) == -1)       //Check that the 2nd atom is in the ringsys
    {
        return -1;
    }
    if ((xatom3 = GetAtomIndex(atom3)) == -1)       //Check that the 3rd atom is in the ringsys
    {
        return -1;
    }
    if ((xatom4 = GetAtomIndex(atom4)) == -1)       //Check that the 4th atom is in the ringsys
    {
        return -1;
    }

    if (atom1->hybrid() != 3)                 //Check that the 1st atom is aliphatic
    {
        return -1;
    }
    if (atom2->hybrid() != 3)                 //Check that the 2nd atom is aliphatic
    {
        return -1;
    }
    if (atom3->hybrid() != 3)                 //Check that the 3rd atom is aliphatic
    {
        return -1;
    }
    if (atom4->hybrid() != 3)                 //Check that the 4th atom is aliphatic
    {
        return -1;
    }

    int *path = new int[_atoms.size()];     //Get space for the arrays used in the
    bool *used = new bool[_atoms.size()];   //search
    InitializeArray(used, _atoms.size(), false);
    int min = _atoms.size();

    path[0] = xatom1;                       //Assumes the four atoms are in sequence
    path[1] = xatom2;
    path[2] = xatom3;                       //Initialize the search arrays, starting
    path[3] = xatom4;                       //with the input atoms
    used[xatom1] = true;
    used[xatom2] = true;
    used[xatom3] = true;
    used[xatom4] = true;

    ExtendAliphaticPath(path, used, 4, min);    //Perform the recursive, depth-first search

    delete[] path;
    delete[] used;

    return min;
}

/////////////////////////////////////////////////////////////////////////////
// Function:    ExtendAliphaticPath
// Purpose:		Recursive, depth-first search for the smallest cycle in a ring system
//				that consists only of sp3-hybridized atoms
// Input:       Sequence of indices of atoms representing a path thru the mol
//				Array of booleans tracking which atoms have been used in the path
//				Count of the number of atoms used in the path
// Output:      Stores the size of the smallest ring found so far in the integer
//				referred to in the final argument.
// Requires:	The method GenerateConnTable() has been run for this object
/////////////////////////////////////////////////////////////////////////////

void RingSystem::ExtendAliphaticPath(int *path,     //Ordered array of atom indices in the path
                                     bool *used,    //Array that tracks which atoms have been used
                                     int depth,         //Number of atoms in current path
                                     int &min) const     //Size of smallest ring found so far
{
    int i;
    for (i = 0; i < (int)_atoms.size(); ++i)           //Loop thru atoms in ringsys
    {
        if (_conn_table[path[depth-1]][i] == false)
        {
            continue;                           //Skip if this atom not connected to last
        }

        if (_atoms[i]->hybrid() != 3)
        {
            continue;
        }

        if (i == path[0] && depth > 2)          //Check if this completes a cycle
        {
            min = depth;
        }
        else if (used[i] == false && depth + 1 < min)
        {
            path[depth] = i;
            used[i] = true;
            ExtendAliphaticPath(path, used, depth + 1, min);
            used[i] = false;                    //Reset for next bond
        }
    }                                           //End loop over atoms
}

/////////////////////////////////////////////////////////////////////////////
// Function:    GetAromaticSys
// Purpose:		Identify which aromatic system, if any, a given atom belongs to.
// Input:       Ptr to the atom of interest
// Output:      Returns a ptr to the aromatic system, or zero if the atom is not
//				aromatic
// Requires:	The method DetectAromatics() has been run for the ring system
/////////////////////////////////////////////////////////////////////////////
Aromatic*RingSystem::GetAromaticSys(const MIAtom *atom)
{
    std::vector<Aromatic>::iterator arom;

    for (arom = _aromatics.begin(); arom != _aromatics.end(); ++arom)
    {
        if (arom->Contains(atom))
        {
            return &(*arom);
        }
    }
    return 0;
}

/////////////////////////////////////////////////////////////////////////////
// Function:    GetAromaticSys
// Purpose:		Identify which aromatic system, if any, a given bond belongs to.
// Input:       Ptr to the bond of interest
// Output:      Returns a ptr to the aromatic system, or zero if the bond is not
//				aromatic
// Requires:	The method DetectAromatics() has been run for the ring system
/////////////////////////////////////////////////////////////////////////////
Aromatic*RingSystem::GetAromaticSys(const Bond *bond)
{
    std::vector<Aromatic>::iterator arom;

    for (arom = _aromatics.begin(); arom != _aromatics.end(); ++arom)
    {
        if (arom->Contains(bond))
        {
            return &(*arom);
        }
    }
    return 0;
}

/////////////////////////////////////////////////////////////////////////////
// Function:    GetAromaticSys
// Purpose:		Identify which aromatic system, if any, a given triplet of
//				atoms belongs to.
// Input:       Ptrs to the 3 atoms of interest
// Output:      Returns a ptr to the aromatic system, or zero if the 3 atoms
//				are not in the same aromatic system.
// Requires:	The method DetectAromatics() has been run for the ring system
/////////////////////////////////////////////////////////////////////////////
Aromatic*RingSystem::GetAromaticSys(const MIAtom *atom1,
                                    const MIAtom *atom2,
                                    const MIAtom *atom3)
{
    std::vector<Aromatic>::iterator arom;

    for (arom = _aromatics.begin(); arom != _aromatics.end(); ++arom)
    {
        if (arom->Contains(atom1)
            && arom->Contains(atom2)
            && arom->Contains(atom3))
        {
            return &(*arom);
        }
    }
    return 0;
}

/////////////////////////////////////////////////////////////////////////////
// Function:    IsPlanar
// Purpose:		Assess whether the atoms of a ring system all lie within a plane,
//				generally used to guess aromaticity from 3D coordinates
// Input:       A cutoff value for the average deviation of the atoms
// Output:      True if the ring is planar, false otherwise
// Requires:
/////////////////////////////////////////////////////////////////////////////
bool RingSystem::IsPlanar(float tolerance) const
{
    double p_displace;

    std::vector<double> normal;         //Initialize the vector to be sure that we
    normal.push_back(0.0);              //have enough space to store the x, y, and z
    normal.push_back(0.0);              //values--and to be sure that &normal[0]
    normal.push_back(0.0);              //produces a valid pointer

    LSqrPlane(&normal[0], &p_displace);

    double d;
    double sum_d = 0.0;
    MIAtom_const_iter a, ae;
    ae = _atoms.end();
    a = _atoms.begin();
    while (a != ae)
    {
        d =  (*a)->x() * normal[0];
        d += (*a)->y() * normal[1];
        d += (*a)->z() * normal[2];
        d -= p_displace;

        sum_d += fabs(d);
        ++a;
    }

    return (sum_d / _atoms.size()) < tolerance;
}

}
