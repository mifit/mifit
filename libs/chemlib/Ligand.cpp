#include <string.h>
#include <numeric>
#include <algorithm>
#include <list>
#include <functional>
#include <iterator>

#include "Ligand.h"
#include "mol_util.h"
#include "sequence_util.h"
#include "transform_util.h"
#include "substituent_util.h"
#include "LigandPerceiver.h"
#include "atom_util.h"
#include <chemlib/RESIDUE_.h>

using namespace std;

namespace chemlib
{

Ligand::Ligand(int type)
{
    ModelType = type;
}

Ligand::Ligand(const RESIDUE &res, const std::vector<Bond> &orig_bonds)
{
    Residue *new_res = new Residue(res);
    residues.push_back(new_res);

    const MIAtomList &orig_atoms = res.atoms();      //Create a vector of ptrs to the atoms in the original residue
    const MIAtomList &new_atoms = residues.back()->atoms();

    std::vector<Bond> trimmed_bonds;
    TrimBonds(trimmed_bonds, orig_bonds, orig_atoms);       //Copy only bonds for which we have atoms

    Bond *bond;
    std::vector<Bond>::const_iterator bnd;
    for (bnd = trimmed_bonds.begin(); bnd != trimmed_bonds.end(); ++bnd)
    {
        bond = new Bond(*bnd, res.atoms(), new_atoms);                //Create a new Bond object
        AddBond(*bond);
        delete bond;
    }
}

Ligand::Ligand(const RESIDUE *res, const std::vector<Bond> &orig_bonds)
{
    Residue *new_res = new Residue(*res);
    residues.push_back(new_res);

    const MIAtomList &orig_atoms = res->atoms();     //Create a vector of ptrs to the atoms in the original residue
    const MIAtomList &new_atoms = residues.back()->atoms();

    std::vector<Bond> trimmed_bonds;
    TrimBonds(trimmed_bonds, orig_bonds, orig_atoms);       //Copy only bonds for which we have atoms

    Bond *bond;
    std::vector<Bond>::const_iterator bnd;
    for (bnd = trimmed_bonds.begin(); bnd != trimmed_bonds.end(); ++bnd)
    {
        bond = new Bond(*bnd, orig_atoms, new_atoms);               //Create a new Bond object
        AddBond(*bond);
        delete bond;
    }
}

Ligand::Ligand(const std::vector<Residue*> &orig_res, const std::vector<Bond> &orig_bonds)
{
    if (orig_res.empty())
    {
        return;
    }

    std::vector<Residue*>::const_iterator res;
    for (res = orig_res.begin(); res != orig_res.end(); ++res)
    {
        residues.push_back(new Residue(**res));
    }

    MIAtomList new_atoms;
    MIAtomList orig_atoms;
    GatherAtmPtrs(new_atoms, residues);
    GatherAtmPtrs(orig_atoms, orig_res);

    std::vector<Bond> trimmed_bonds;
    TrimBonds(trimmed_bonds, orig_bonds, orig_atoms);       //Copy only bonds for which we have atoms

    std::vector<Bond>::iterator bnd;
    for (bnd = trimmed_bonds.begin(); bnd != trimmed_bonds.end(); ++bnd)
    {
        bnd->UpdateAtoms(orig_atoms, new_atoms);
        AddBond(*bnd);
    }

}

Ligand::~Ligand()
{
    // clean up residues;
    std::vector<Residue*>::iterator ri;
    for (ri = residues.begin(); ri != residues.end(); ++ri)
    {
        Residue *res = *ri;
        delete res;
    }
}

void Ligand::Export(std::vector<RESIDUE*> &rdues, std::vector<Bond> &bnds)
{
    MIAtomList mifit_atms;
    MIAtomList orig_atms;

    std::vector<Residue*>::iterator ri;
    for (ri = residues.begin(); ri != residues.end(); ++ri)
    {
        Residue *res = *ri;
        for (unsigned int i = 0; i < res->atoms().size(); ++i)
        {
            orig_atms.push_back(res->atom(i));
        }
        RESIDUE *r = new RESIDUE(*res);
        copy(r->atoms().begin(), r->atoms().end(), back_inserter(mifit_atms));
        rdues.push_back(r);
    }

    std::vector<Bond>::iterator bnd;
    for (bnd = bonds.begin(); bnd != bonds.end(); ++bnd)
    {
        bnds.push_back(Bond(*bnd, orig_atms, mifit_atms));
    }
}

void Ligand::GetConstraints(const MIAtomList &mifit_atms,
                            const std::vector<RESIDUE*> &mifit_residues,
                            std::vector<Bond> &bond_lengths,
                            std::vector<ANGLE> &angles,
                            std::vector<TORSION> &torsions,
                            std::vector<TORSION> &impropers,
                            std::vector<PLANE> &planes,
                            std::vector<CHIRAL> &chirals)
{

    //Create a vector of ptrs to the MIAtom objects
    std::vector<const MIAtom*> orig_atms;
    std::vector<Residue*>::iterator ri;
    for (ri = residues.begin(); ri != residues.end(); ++ri)
    {
        Residue *res = *ri;
        for (unsigned int i = 0; i < res->atoms().size(); ++i)
        {
            orig_atms.push_back(res->atom(i));
        }
    }

    //Create a vector of ptrs to the residues objects
    std::vector<const Residue*> orig_residues;
    for (size_t i = 0; i< residues.size(); ++i) // was CpyResPtrs(orig_residues, residues);
    {
        orig_residues.push_back(residues[i]);
    }

    //In turn, create MIFit structs for each type of constraint
    std::vector<BondLength>::iterator bnd;
    for (bnd = geometry.Bonds.begin(); bnd != geometry.Bonds.end(); ++bnd)
    {
        bond_lengths.push_back(bnd->ToBond(mifit_atms, orig_atms));
    }

    std::vector<Angle>::iterator ang;
    for (ang = geometry.Angles.begin(); ang != geometry.Angles.end(); ++ang)
    {
        angles.push_back(ang->ToANGLE(mifit_atms, orig_atms));
    }

    std::vector<Improper>::iterator imp;
    for (imp = geometry.Impropers.begin(); imp != geometry.Impropers.end(); ++imp)
    {
        if (imp->DictTest())
        {
            impropers.push_back(imp->ToTORSION(mifit_atms, orig_atms,
                                               mifit_residues, orig_residues));
        }
    }

    std::vector<Torsion>::iterator tors;
    for (tors = geometry.Torsions.begin(); tors != geometry.Torsions.end(); ++tors)
    {
        torsions.push_back(tors->ToTORSION(mifit_atms, orig_atms,
                                           mifit_residues, orig_residues));
    }

    std::vector<Plane>::iterator pln;
    for (pln = geometry.Planes.begin(); pln != geometry.Planes.end(); ++pln)
    {
        planes.push_back(pln->ToPLANE(mifit_atms, orig_atms,
                                      mifit_residues, orig_residues));
    }

    std::vector<const MIAtom*> done_centers;
    std::vector<Chiral>::iterator chrl;
    for (chrl = geometry.Chirals.begin(); chrl != geometry.Chirals.end(); ++chrl)
    {
        if (find(done_centers.begin(), done_centers.end(), chrl->GetCenter()) == done_centers.end())
        {
            chirals.push_back(chrl->ToCHIRAL(mifit_atms, orig_atms));
            done_centers.push_back(chrl->GetCenter());
        }
    }
}

Residue*Ligand::AddRes(const std::string &type, const std::string &name, unsigned short linkage_type,
                       unsigned short chain_id, char secstr)
{

    Residue target;
    target.setType(type);
    target.setName(name);
    target.set_linkage_type(linkage_type);
    target.set_chain_id(chain_id);
    target.setSecstr(secstr);

    residues.push_back(new Residue(target));
    return residues.back();
}

Residue*Ligand::AddRes(const RESIDUE &res, const std::vector<Bond> &orig_bonds)
{
    Residue *new_res = new Residue(res);
    residues.push_back(new_res);

    const MIAtomList &orig_atoms = res.atoms();      //Create a vector of ptrs to the atoms in the original residue (The order of these is important.)
    const MIAtomList &new_atoms = residues.back()->atoms();

    //Add bonds between atoms in the new residue (intra-residue bonds)
    //Don't add any links between this residue and other residues
    Bond *bond;
    std::vector<Bond>::const_iterator bnd;
    for (bnd = orig_bonds.begin(); bnd != orig_bonds.end(); ++bnd)
    {

        if (bnd->getAtom1() == bnd->getAtom2())
        {
            continue;
        }

        //Check if this bond is between two atoms in this residue
        if (std::find(orig_atoms.begin(), orig_atoms.end(), bnd->getAtom1()) == orig_atoms.end())
        {
            continue;
        }
        if (std::find(orig_atoms.begin(), orig_atoms.end(), bnd->getAtom2()) == orig_atoms.end())
        {
            continue;
        }

        bond = new Bond(*bnd, orig_atoms, new_atoms);
        AddBond(*bond);
        //		bonds.push_back(*bond);
        delete bond;
    }

    return residues.back();
}

void Ligand::AddBond(Bond &bond)
{
    if (bond.getAtom1() == bond.getAtom2())
    {
        return;
    }

    bond.getAtom1()->addBondnumber(bonds.size());
    bond.getAtom2()->addBondnumber(bonds.size());
    bond.getAtom1()->addNabor(bond.getAtom2());
    bond.getAtom2()->addNabor(bond.getAtom1());
    bonds.push_back(bond);
}

void Ligand::AddBond(MIAtom *a1, MIAtom *a2)
{
    a1->addBondnumber(bonds.size());
    a2->addBondnumber(bonds.size());
    a1->addNabor(a2);
    a2->addNabor(a1);

    Bond bond;
    bond.setAtom1(a1);
    bond.setAtom2(a2);
    bond.type = B_CONNECT;
    bonds.push_back(bond);
}

void Ligand::AddBond(MIAtom *a1, MIAtom *a2, unsigned char bond_order, char stereo)
{
    a1->addBondnumber(bonds.size());
    a2->addBondnumber(bonds.size());
    a1->addNabor(a2);
    a2->addNabor(a1);

    Bond bond;
    bond.setAtom1(a1);
    bond.setAtom2(a2);
    bond.type = B_CONNECT;
    bond.setOrder(bond_order);
    bond.stereo = stereo;
    bonds.push_back(bond);
}

/////////////////////////////////////////////////////////////////////////////
// Function:    ClearBondOrder
// Purpose:		Sets all bonds to single, and all atoms to sp3.  Also removes
//				any implicit hydrogens.
// Input:       None
// Output:      None
// Requires:
/////////////////////////////////////////////////////////////////////////////
void Ligand::ClearBondOrders()
{
    for (size_t i = 0; i< bonds.size(); ++i)
    {
        bonds[i].ClearOrder();
    }

    for (size_t i = 0; i < residues.size(); ++i)
    {
        Residue *r = residues[i];
        for (int i = 0; i < r->atomCount(); ++i)
        {
            r->atom(i)->ClearHybrid();
        }
    }
}

static bool ReverseSortBy2nd(const std::pair< MIAtom_const_iter, int> &a1,
                             const std::pair< MIAtom_const_iter, int> &a2)
{
    return a1.second > a2.second;
}


void Ligand::GuessBondOrders()
{


    //*
    //Step 1: Iterate over all the atoms, guessing its hybridization from the bond
    //        angles or (for terminal atoms) the bond length.
    //*
    std::vector<Residue*>::iterator r, re;
    MIAtom_const_iter a, ae;

    int natoms = GetNumAtoms();

    re = residues.end();
    r = residues.begin();
    while (r != re)
    {

        ae = (*r)->atoms().end();
        a = (*r)->atoms().begin();
        while (a != ae)
        {
            switch ((*a)->nabors().size())
            {
            case 0: break;                                          //unconnected atoms
            case 1: if (natoms == 2)
                {
                    HybridizeTerminalAtom(**a, bonds);
                }
                break;
            default: HybridFromGeom(**a);
            }
            ++a;
        }
        ++r;
    }

    //*
    //Step 2: Iterate over all the ring systems, setting those that are planar (within a marg-of-err)
    //        to be aromatic, with sp2/aromatic atoms and partialdouble/aromatic bonds
    //*
    std::vector<RingSystem>::iterator rs, rse;
    rse = ringsystems.end();

    rs = ringsystems.begin();
    while (rs != rse)
    {
        if ((rs->NumAtoms() > 4) && rs->IsPlanar(MolFrom3D::PLANE_TOLERANCE))
        {
            rs->SetAllAromatic();
        }
        ++rs;
    }

    //*
    //Step 3: Loop over atoms, and iteratively increase bond orders until the bonds
    //        are consistent with the hybridization state assigned in step 2
    //*
    std::vector< std::pair<MIAtom_const_iter, int> > sorted_atoms;

    re = residues.end();
    r = residues.begin();
    while (r != re)
    {

        ae = (*r)->atoms().end();
        a = (*r)->atoms().begin();
        while (a != ae)
        {
            std::pair<MIAtom_const_iter, int>
            atom_to_sort(a, Electroneg((*a)->atomicnumber()));
            sorted_atoms.push_back(atom_to_sort);
            ++a;
        }
        ++r;
    }

    sort(sorted_atoms.begin(), sorted_atoms.end(), ReverseSortBy2nd);

    Bond *bond;
    MIAtom *partner;

    std::vector< std::pair<MIAtom_const_iter, int> >::iterator ap, ape;
    ape = sorted_atoms.end();
    ap = sorted_atoms.begin();

    while (ap != ape)
    {
        int degree = (*(ap->first))->nabors().size();

        //If a nitrogen appears to have all single bonds, skip it so we don't add
        //double bonds to nitrogens in amides, sulfonamides, anilines, etc.
        if ((*ap->first)->IsNitrogen() && PredictValenceFrom3D(*(*ap->first)) == degree)
        {
            ++ap;
            continue;
        }

        if ((*ap->first)->hybrid() == 1
            && UnusedValences(*(*ap->first), bonds) >= 2
            && (partner = GetTriplePartner(*ap->first, bonds)) != 0)
        {
            bond = GetBond(*ap->first, partner);
            bond->SetTriple();
        }
        else if ((*ap->first)->hybrid() == 2
                 && UnusedValences(*(*ap->first), bonds) >= 1
                 && CurrentValence(*(*ap->first), bonds) == degree
                 && (partner = GetDoublePartner(*ap->first, bonds)) != 0)
        {
            bond = GetBond(*ap->first, partner);
            bond->SetDouble();
        }
        ++ap;
    }
    LigandPerceiver lp;
    lp.AssignHybridization(this);
    lp.AssignImpHydrogens(this);
}

void Ligand::Connect(const Bond &connect)
{
    Bond bond;
    MIAtom *a1 = connect.getAtom1();
    MIAtom *a2 = connect.getAtom2();
    if (!AlreadyBonded(a1, a2))
    {
        a1->addNabor(a2);
        a2->addNabor(a1);
        a1->addBondnumber(bonds.size());
        a2->addBondnumber(bonds.size());

        bond.setAtom1(a1);
        bond.setAtom2(a2);
        bond.getAtom1()->setType(AtomType::BONDED);
        bond.getAtom2()->setType(AtomType::BONDED);
        bond.type = B_CONNECT;
        bonds.push_back(bond);
    }
}

void Ligand::FixAtomicNumbers()
{
    // fixes the atomicnumber field in a residue
    // sometimes this is a bit of a guess becuase it is not
    // clear what the atom type is by name alone!
    std::vector<Residue*>::iterator ri;
    unsigned int i;
    bool ispep, isnuc;
    //	char *name;
    //	name = (char *)malloc(sizeof(char) * MAXNAME);
    std::string name;
    for (ri = residues.begin(); ri != residues.end(); ++ri)
    {
        Residue *res = *ri;
        ispep = IsPeptide(*res) != 0;
        isnuc = IsNucleic(&(*res)) != 0;
        for (i = 0; i < res->atoms().size(); i++)
        {
            if (res->atom(i)->atomicnumber() <= 0 || res->atom(i)->atomicnumber() > NELEMENTS)
            {
                if (isdigit(res->atom(i)->name()[0]))
                {
                    name = " ";
                    name += res->atom(i)->name()[1];
                    res->atom(i)->setAtomicnumber(Atomic_Number(name.c_str()));
                }
                else if (ispep || isnuc)
                {
                    name = " ";
                    name += res->atom(i)->name()[0];
                    res->atom(i)->setAtomicnumber(Atomic_Number(name.c_str()));
                    if (strcmp(res->name().c_str(), "MET") == 0 && strncmp(res->atom(i)->name(), "SE", 2) == 0)
                    {
                        res->atom(i)->setAtomicnumber(Atomic_Number("SE"));
                    }
                }
                else
                {
                    name = res->atom(i)->name()[0];
                    name += res->atom(i)->name()[1];
                    res->atom(i)->setAtomicnumber(Atomic_Number(name.c_str()));
                    if (res->atom(i)->atomicnumber() == -1)
                    {
                        name = " ";
                        name += res->atom(i)->name()[0];
                        res->atom(i)->setAtomicnumber(Atomic_Number(name.c_str()));
                    }
                }
            }
        }
    }
}

int Ligand::FindRingSystems()
{
    ringsystems.clear();
    int nringsystems = 0;           //Number of ring systems

    //Initialize the variables in the atom and bond structures
    InitRingData();

    //Determine which bonds and atoms are cyclic and acyclic
    SetRingFlags();

    //Now cluster cyclic atoms into ring systems
    //by looping over all cyclic atoms
    unsigned int i;
    std::vector<Residue*>::iterator ri;
    RingSystem rs;
    for (ri = residues.begin(); ri != residues.end(); ++ri)
    {
        Residue *res = *ri;
        for (i = 0; i < res->atoms().size(); ++i)
        {
            if (!res->atom(i)->iscyclic())
            {
                continue;                           //Skip if this atom is not cyclic,
            }                                       //or if it has already been assigned
            if (res->atom(i)->ring_system() != -1)
            {
                continue;
            }

            //Init the Ring System object
            rs.Clear();
            rs.SetMolecule(this);
            rs.SetIndex(ringsystems.size());
            rs.AddAtom(res->atom(i));

            //Recursively add bonds and atoms to the ring system
            rs.Extend(res->atom(i));

            //Generate some of the internal data that the ring system object needs
            rs.GenerateConnTable();
            rs.DetectFusedRing();
            rs.DetectAromatics();
            rs.GenerateSRData();

            //Add this ringsystem to the molecule
            ringsystems.push_back(rs);
            nringsystems++;
        }                                           //close loop atoms
    }                                               //close loop residues
    return nringsystems;
}

void Ligand::InitRingData()
{
    unsigned int i;
    std::vector<Residue*>::iterator ri;
    for (ri = residues.begin(); ri != residues.end(); ++ri)
    {
        Residue *res = *ri;
        for (i = 0; i < res->atoms().size(); ++i)
        {
            res->atom(i)->setIscyclic(false);
            res->atom(i)->set_ring_system(-1);
            res->atom(i)->set_smallest_aromatic_ring(-1);
        }
    }

    std::vector<Bond>::iterator bnd;
    for (bnd = bonds.begin(); bnd != bonds.end(); ++bnd)
    {
        bnd->iscyclic = false;
        bnd->ring_system = -1;
        bnd->smallest_ring_size = -1;
        bnd->smallest_aromatic_ring = -1;
    }

}

void Ligand::SetRingFlags()
{
    std::list<MIAtom*> near_atom1;
    std::list<MIAtom*> near_atom2;

    MIAtom_const_iter atm;
    MIAtom *root;

    std::vector<Bond>::iterator bnd;
    for (bnd = bonds.begin(); bnd != bonds.end(); ++bnd)
    {
        ResetSearchFlags();                     //Clear the workspaces, which are the
        near_atom1.clear();                     //two lists of atoms ptrs and the
        near_atom2.clear();                     //search_flags in the atom structs

        near_atom1.push_back(bnd->getAtom1());      //Start with just the 2 atoms
        near_atom2.push_back(bnd->getAtom2());      //in this bond
        bnd->getAtom1()->set_search_flag(1);
        bnd->getAtom2()->set_search_flag(2);

        //Now perform 2 breadth-first searches simultaneously from
        //both atoms in the bond.  If the two should meet, as detected
        //by the "search_flag", declare this bond cyclic

        while (near_atom1.size() > 0 && near_atom2.size() > 0)
        {
            root = near_atom1.front();
            for (atm = root->nabors().begin(); atm != root->nabors().end(); ++atm)
            {
                if (*atm == bnd->getAtom2())        //Prevents backtracking across the
                {
                    continue;                   //original bond
                }
                if ((*atm)->search_flag() == 2)
                {
                    bnd->iscyclic = true;
                }
                else if ((*atm)->search_flag() == 0)
                {
                    (*atm)->set_search_flag(1);
                    near_atom1.push_back(*atm);
                }
            }
            near_atom1.pop_front();

            root = near_atom2.front();
            for (atm = root->nabors().begin(); atm != root->nabors().end(); ++atm)
            {
                if (*atm == bnd->getAtom1())        //Prevents backtracking across the
                {
                    continue;                   //original bond
                }
                if ((*atm)->search_flag() == 1)
                {
                    bnd->iscyclic = true;
                }
                else if ((*atm)->search_flag() == 0)
                {
                    (*atm)->set_search_flag(2);
                    near_atom2.push_back(*atm);
                }
            }
            near_atom2.pop_front();

            if (bnd->iscyclic == 1)
            {
                break;
            }
        } //Breadth-first search
    } //Loop over bonds

    ResetSearchFlags();


    //Set the "iscyclic" flags in the atom structure, using the fact that
    //any atom with a cyclic bond is cyclic
    unsigned int i, j;
    Bond *bond;
    std::vector<Residue*>::iterator ri;
    for (ri = residues.begin(); ri != residues.end(); ++ri)
    {
        Residue *res = *ri;
        for (i = 0; i < res->atoms().size(); ++i)
        {
            for (j = 0; j < res->atom(i)->bondnumbers().size(); ++j)
            {
                bond = GetBond(res->atom(i)->bondnumbers()[j]);
                if (bond != 0 && bond->iscyclic)
                {
                    res->atom(i)->setIscyclic(true);
                    break;
                }
            }
        }
    }
}

void Ligand::ResetSearchFlags()
{
    unsigned int i;
    std::vector<Residue*>::iterator ri;
    for (ri = residues.begin(); ri != residues.end(); ++ri)
    {
        Residue *res = *ri;
        for (i = 0; i < res->atoms().size(); ++i)
        {
            res->atom(i)->set_search_flag(0);
        }
    }
}

MIAtom*Ligand::MoreCentralAtm(MIAtom *atom1, MIAtom *atom2)
{
    std::vector <MIAtom*> near_atom1;
    std::vector <MIAtom*> near_atom2;

    int natoms = GetNumAtoms();

    near_atom1.reserve(natoms);
    near_atom2.reserve(natoms);

    DepthFirstSearch(atom1, atom2, near_atom1);
    DepthFirstSearch(atom2, atom1, near_atom2);

    if (near_atom2.size() > near_atom1.size())
    {
        return atom2;
    }
    else
    {
        return atom1;
    }
}

void Ligand::DepthFirstSearch(MIAtom *root, MIAtom *forbid, std::vector <MIAtom*> &aggregate)
{

    if (std::find(aggregate.begin(), aggregate.end(), root) == aggregate.end())
    {
        aggregate.push_back(root);                  //Add this atom to the aggregation
    }
    else
    {
        return;
    }

    std::vector <MIAtom*>::const_iterator nabor;
    for (nabor = root->nabors().begin(); nabor != root->nabors().end(); ++nabor)
    {
        if (*nabor != forbid)
        {
            DepthFirstSearch(*nabor, root, aggregate);
        }
    }
}

//**GetNabor															**//
//**Takes an atom pointer and the index (starting with zero) of the		**//
//**neighbor and uses the "bondnumbers" std::vector in the atom struct to	**//
//**return a pointer to the neighboring atom along that bond			**//
MIAtom*Ligand::GetNabor(const MIAtom *source, int nabor_index)
{
    Bond *bond = GetBond(source->bondnumbers()[nabor_index]);
    if (bond == 0)
    {
        return 0;
    }
    else if (bond->getAtom1() == source)
    {
        return bond->getAtom2();
    }
    else if (bond->getAtom2() == source)
    {
        return bond->getAtom1();
    }
    else
    {
        return 0;
    }
}

//Used to get a neighboring atom of the first argument.  MIAtom is arbitrary,
//but it cannot be the one specified in the second argument.
MIAtom*Ligand::GetNewNabor(MIAtom *root, MIAtom *forbid)
{
    std::vector <MIAtom*>::const_iterator nabor;
    for (nabor = root->nabors().begin(); nabor != root->nabors().end(); ++nabor)
    {
        if (*nabor != forbid)
        {
            return *nabor;
        }
    }
    return 0;
}

Bond*Ligand::GetBond(int index)
{
    if (index >= 0 && index < (int)bonds.size())
    {
        return &bonds[index];
    }
    else
    {
        return 0;
    }
}

const Bond*Ligand::GetBond(int index) const
{
    if (index >= 0 && index < (int)bonds.size())
    {
        return &bonds[index];
        //		return &(bonds.at(index));
    }
    else
    {
        return 0;
    }
}

Bond*Ligand::GetBond(const MIAtom *a1, const MIAtom *a2)
{
    for (unsigned int i = 0; i < bonds.size(); i++)
    {
        const Bond &b = bonds[i];
        if ((b.getAtom1() == a1 && b.getAtom2() == a2)
            || (b.getAtom1() == a2 && b.getAtom2() == a1) )
        {
            return &bonds[i];
        }
    }
    return 0;
}

bool Ligand::AlreadyBonded(const MIAtom *a1, const MIAtom *a2) const
{
    for (unsigned int i = 0; i < bonds.size(); i++)
    {
        const Bond &b = bonds[i];
        if ((b.getAtom1() == a1 && b.getAtom2() == a2)
            || (b.getAtom1() == a2 && b.getAtom2() == a1) )
        {
            return true;
        }
    }
    return false;
}

void Ligand::FreeBonds()
{
    //	nlinks=0;
    bonds.clear();
}

//Function to get all the bonds in a Ligand that contain at least one atom from
//the specified residue.
void Ligand::GetResidueBonds(const Residue *res, std::vector<Bond*> &resbonds)
{
    std::vector<int>::const_iterator xbond; //Iterator for looping over bonds to an atom
    Bond *bnd;              //Points to current bond
    MIAtom *atm;                //Points to current atom
    unsigned int i;                 //MIAtom counter

    for (i = 0; i < res->atoms().size(); ++i)
    {
        atm = res->atom(i);
        for (xbond = atm->bondnumbers().begin(); xbond != atm->bondnumbers().end(); ++xbond)
        {
            bnd = GetBond(*xbond);
            if (std::find(resbonds.begin(), resbonds.end(), bnd) == resbonds.end())
            {
                resbonds.push_back(bnd);                    //Check if we already have this
            }                                               //one, otherwise add it
        }
    }
}

//Function to get all the bonds in a molecule that are contained within the
//the specified residue.
void Ligand::GetResidueIntBonds(const Residue *res, std::vector<Bond> &resbonds)
{
    std::vector<int>::const_iterator xbond; //Iterator for looping over bonds to an atom

    MIAtom *atm;                //Points to current atom
    unsigned int i;                 //MIAtom counter

    std::vector<int> halfbonds;

    for (i = 0; i < res->atoms().size(); ++i)      //Loop over atoms in this residue,
    {
        atm = res->atom(i);                //then over bonds to that atom.
        for (xbond = atm->bondnumbers().begin(); xbond != atm->bondnumbers().end(); ++xbond)
        {
            halfbonds.push_back(*xbond);
        }
    }

    std::sort(halfbonds.begin(), halfbonds.end());

    xbond = halfbonds.begin();
    for (;;)
    {
        std::vector<int>::const_iterator hbond_end = halfbonds.end();
        if ((xbond = std::adjacent_find(xbond, hbond_end))
            != halfbonds.end())
        {
            resbonds.push_back(bonds[*xbond]);
            ++xbond;
        }
        else
        {
            break;
        }
    }
}

int Ligand::GetNumAtoms() const
{
    return std::accumulate(residues.begin(), residues.end(), 0, NumAtomsSum);
}

int Ligand::GetIndex(MIAtom *query, MIAtom **domain, int size)
{
    int i;
    for (i = 0; i < size; i++)
    {
        if (query == domain[i])
        {
            return i;
        }
    }
    return -1;
}

void Ligand::Translate(double xstep, double ystep, double zstep)
{
    std::vector<Residue*>::iterator ri;
    for (ri = residues.begin(); ri != residues.end(); ++ri)
    {
        Residue *res = *ri;
        for (unsigned int i = 0; i < res->atoms().size(); ++i)
        {
            res->atom(i)->translate((float)xstep,
                                    (float)ystep,
                                    (float)zstep);
        }
    }
}

void Ligand::Scale(double scale_factor)
{
    std::vector<Residue*>::iterator ri;
    for (ri = residues.begin(); ri != residues.end(); ++ri)
    {
        Residue *res = *ri;
        for (unsigned int i = 0; i < res->atoms().size(); ++i)
        {
            res->atom(i)->scale(scale_factor);
        }
    }
}

void Ligand::TranslateFragment(MIAtom *source, MIAtom *target, double *v)
{
    MIAtomList frag_atoms;
    DepthFirstSearch(source, target, frag_atoms);

    MIAtom_iter atm;
    for (atm = frag_atoms.begin();
         atm != frag_atoms.end();
         ++atm)
    {
        (*atm)->translate((float)v[0], (float)v[1], (float)v[2]);
    }
}

void Ligand::FlipBond(Bond &bond)
{
    MIAtomList near_atom1;
    MIAtomList near_atom2;

    int natoms = GetNumAtoms();
    near_atom1.reserve(natoms);
    near_atom2.reserve(natoms);

    DepthFirstSearch(bond.getAtom1(), bond.getAtom2(), near_atom1);
    DepthFirstSearch(bond.getAtom2(), bond.getAtom1(), near_atom2);

    //Since we could just as well flip either side, flip the
    //one with fewer atoms
    if (near_atom2.size() > near_atom1.size())
    {
        FlipAtoms(bond, near_atom1);
    }
    else
    {
        FlipAtoms(bond, near_atom2);
    }
}

void Ligand::Flatten()
{
    FlattenRings();
    Unroll();
    FitToXYPlane(this);
    CleanUp2D();
}

void Ligand::CleanUp2D()
{
    std::for_each(ringsystems.begin(),
                  ringsystems.end(),
                  void_mem_fun_ref(&RingSystem::FixExocyclics));

    //Could equalize lengths of acyclic bonds here

    //Could equalize covalent angles at acyclic junctions here
}

void Ligand::FlattenRings()
{
    std::for_each(ringsystems.begin(),
                  ringsystems.end(),
                  void_mem_fun_ref(&RingSystem::Flatten));
}

void Ligand::Unroll()
{
    std::vector<Bond>::const_iterator bnd;
    for (bnd = bonds.begin(); bnd != bonds.end(); ++bnd)
    {
        if (bnd->IsRotatable())
        {
            UnrollAtom(bnd->getAtom1(), bnd->getAtom2());
            UnrollAtom(bnd->getAtom2(), bnd->getAtom1());
        }
    }
}

void Ligand::UnrollAtom(MIAtom *origin, MIAtom *ref)
{
    //Translate whole molecule to align atom "origin" to (0,0,0)
    TranslateToOrigin(origin, this);

    //Rotate whole molecule to align atom "ref" to (X,0,0)
    RotateToXAxis(ref, this);

    std::vector< Substituent > subs;

    //Step 1: Create list of substituents, grouping those that
    //will be manipulated together
    GroupSubstituents(origin, ref, subs);

    //Step 2: Rotate each substituent into the XY-plane
    std::for_each(subs.begin(),
                  subs.end(),
                  void_mem_fun_ref(&Substituent::Squash));

    //Step 3: Sort the list of substituents in counter-clockwise order from
    //		  the positive x-axis
    SortSubsCounterClockwise(subs);

    //Step 4: Equalize spacing between substituents
    if (origin->nabors().size() > 2)
    {
        EqualizeSpacing(subs);
    }

}

void Ligand::GroupSubstituents(MIAtom *origin,
                               MIAtom *omit,
                               std::vector< Substituent > &subs) const
{
    Substituent sub;
    MIAtomList ring_atoms;
    MIAtom_const_iter nbr;

    //If degree(origin) == 3, handle both atoms together, so that the bond
    //angle between them remains the same
    if (origin->nabors().size() == 3)
    {
        sub.Clear();
        sub.SetOrigin(origin);
        sub.SetMolecule(this);
        for (nbr = origin->nabors().begin(); nbr != origin->nabors().end(); ++nbr)
        {
            if (*nbr == omit)
            {
                continue;
            }
            else
            {
                sub.AddBranch(*nbr);
            }
        }
        subs.push_back(sub);
        return;
    }

    //Loop over bonds, adding each atom to the std::vector of substituents,
    //so long as it is not connected by a cyclic bond
    ring_atoms.clear();
    for (nbr = origin->nabors().begin(); nbr != origin->nabors().end(); ++nbr)
    {
        if (*nbr == omit)
        {
            continue;
        }
        sub.Clear();
        sub.SetOrigin(origin);
        sub.SetMolecule(this);

        if (origin->iscyclic() && (*nbr)->iscyclic()
            && origin->ring_system() == (*nbr)->ring_system())
        {
            ring_atoms.push_back(*nbr);
        }
        else
        {
            sub.AddBranch(*nbr);
            subs.push_back(sub);
        }
    }

    //Add the all the atoms connected by cyclic bonds
    //as a single substituent
    if (!ring_atoms.empty())
    {
        sub.Clear();
        sub.SetOrigin(origin);
        sub.SetMolecule(this);
        for (nbr = ring_atoms.begin(); nbr != ring_atoms.end(); ++nbr)
        {
            sub.AddBranch(*nbr);
        }
        subs.push_back(sub);
    }
}

/////////////////////////////////////////////////////////////////////////////
// Function:    LSqrPlane
// Purpose:		Calculates the best-fit plane to all the atoms in the molecule
// Input:       Internal lookup of coordinates from parent mol
// Output:      A unit std::vector normal to the best-fit plane
//				A scalar that multiplies the normal std::vector to extend it from the origin
//				to the best-fit plane (i.e. the signed distance from origin to plane)
// Requires:
/////////////////////////////////////////////////////////////////////////////
void Ligand::LSqrPlane(double vm[3], double *d) const
{

    double xs[3], xxs[3][3], b[3][3];
    double a[3][3], vmi[3], bv[3];
    double zip = 1.0E-5;
    double orm, vm0, ratio0, ratio1, ratio2, rat01, rat02;
    int i, j, kk, nnn, n;
    std::vector<Residue*>::const_iterator ri;

    n = 0;
    for (i = 0; i < 3; i++)
    {
        xs[i] = 0.0;
    }
    for (ri = residues.begin(); ri != residues.end(); ++ri)
    {
        Residue *res = *ri;
        for (unsigned int k = 0; k < res->atoms().size(); ++k)
        {

            /* zip is added to prevent numerical instability if
             * atoms in a plane = 0
             */
            xs[0] += res->atom(k)->x()+zip;
            xs[1] += res->atom(k)->y()+zip;
            xs[2] += res->atom(k)->z()+zip;
            n++;
        }
    }

    for (i = 0; i < 3; i++)
    {
        for (j = 0; j < 3; j++)
        {
            xxs[i][j] = 0.0;
        }
    }

    for (ri = residues.begin(); ri != residues.end(); ++ri)
    {
        Residue *res = *ri;
        for (unsigned int k = 0; k < res->atoms().size(); ++k)
        {
            xxs[0][0] += res->atom(k)->x() * res->atom(k)->x();
            xxs[0][1] += res->atom(k)->x() * res->atom(k)->y();
            xxs[0][2] += res->atom(k)->x() * res->atom(k)->z();
            xxs[1][0] += res->atom(k)->y() * res->atom(k)->x();
            xxs[1][1] += res->atom(k)->y() * res->atom(k)->y();
            xxs[1][2] += res->atom(k)->y() * res->atom(k)->z();
            xxs[2][0] += res->atom(k)->z() * res->atom(k)->x();
            xxs[2][1] += res->atom(k)->z() * res->atom(k)->y();
            xxs[2][2] += res->atom(k)->z() * res->atom(k)->z();
        }
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

    /* choose the largest column std::vector of b as initial solution */
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

void Ligand::RandomizeCoords(std::pair<double, double> xbounds,
                             std::pair<double, double> ybounds,
                             std::pair<double, double> zbounds)
{

    if (xbounds.first > xbounds.second
        || ybounds.first > ybounds.second
        || zbounds.first > zbounds.second)
    {
        throw "Invalid bounds in Ligand::RandomizedCoords";
    }
    double xWidth = xbounds.second - xbounds.first;
    double yWidth = ybounds.second - ybounds.first;
    double zWidth = zbounds.second - zbounds.first;

    std::vector<Residue*>::iterator ri, re = residues.end();
    for (ri = residues.begin(); ri != re; ++ri)
    {
        Residue *res = *ri;
        for (unsigned int i = 0; i < res->atoms().size(); ++i)
        {
            res->atom(i)->setPosition((float)(xbounds.first + xWidth * (double) rand() / (double) RAND_MAX),
                                      (float)(ybounds.first + yWidth * (double) rand() / (double) RAND_MAX),
                                      (float)(zbounds.first + zWidth * (double) rand() / (double) RAND_MAX));
        }
    }
}

void Ligand::MapAtomPtrs(std::map<MIAtom*, int> &atom_map)
{

    std::vector<Residue*>::iterator ri, re = residues.end();
    int xAtom = 0;

    for (ri = residues.begin(); ri != re; ++ri)
    {
        Residue *res = *ri;
        int n = res->atoms().size();
        MIAtom *atm = res->atom(0);
        for (int i = 0; i < n; ++i, ++atm)
        {
            atom_map[atm] = xAtom++;
        }
    }
}

void Ligand::MapAtomPtrs(std::map<const MIAtom*, int> &atom_map) const
{

    std::vector<Residue*>::const_iterator ri, re = residues.end();
    int xAtom = 0;

    for (ri = residues.begin(); ri != re; ++ri)
    {
        Residue *res = *ri;
        int n = res->atoms().size();
        const MIAtom *atm = res->atom(0);
        for (int i = 0; i < n; ++i, ++atm)
        {
            atom_map[atm] = xAtom++;
        }
    }
}

}
