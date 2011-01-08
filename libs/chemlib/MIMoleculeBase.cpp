#include <algorithm>

#include <math/mathlib.h>

#include <chemlib/Residue.h>
#include "MIMoleculeBase.h"

#include "PDB.h"
#include "mol_util.h"

#include "MIMolIOBase.h"
#include "MIMolDictionary.h" // FIXME: remove this dependency?
#include "mol_util_private.h"


namespace chemlib
{

//#define MEMORY_CORRUPTION_DEBUG
#if defined(MEMORY_CORRUPTION_DEBUG)
static unsigned int DELETION_COUNT = 0;
static unsigned int magic_delete = UINT_MAX;
static std::map<MIMoleculeBase*, unsigned int> DELETED_MOLECULES;
#endif

MIMoleculeBase::MoleculeRefCountMap MIMoleculeBase::refCounts;

bool MIMoleculeBase::isValid(MIMoleculeBase *mol)
{
    bool result = false;
    if (mol != NULL && refCounts.find(mol) != refCounts.end())
    {
        result = refCounts[mol] > 0;
    }
#ifdef MEMORY_CORRUPTION_DEBUG
    if (!result && mol)
    {
        printf("Invalid molecule %p queried, was deletion #%d\n", mol, DELETED_MOLECULES[mol]);
    }
#endif
    return result;
}

MIMoleculeBase::MIMoleculeBase()
    : SymmResidues(NULL)
{
    ++refCounts[this];
    compound = "";

    Tatom1 = Tatom2 = NULL;
    residues = NULL;

    coords_changed = true;
    modified = true;
    strcpy(link_here, "C");  // protein link thru peptide bond
    strcpy(link_next, "N");  // in derived class you can change         to other bond type.

    nlinks = nresidues = 0;
}

MIMoleculeBase::MIMoleculeBase(Residue *reslist, const std::string &cmpd, Bond *conns, int nconns)
    : SymmResidues(NULL)
{
    ++refCounts[this];
    compound = cmpd;

    Tatom1 = Tatom2 = NULL;
    residues = reslist;

    strcpy(link_here, "C");  // protein link thru peptide bond
    strcpy(link_next, "N");  // in derived class you can change         to other bond type.

    nlinks = nresidues = 0;

    for (int i = 0; i < nconns; i++)
    {
        connects.push_back(conns[i]);
    }

    InitSeqPos();
    Build();

    coords_changed = true;
    modified = true;
}

MIMoleculeBase::~MIMoleculeBase()
{
    // follow residue list freeing as we go;
    if (residues)
    {
        FreeResidueList(residues);
        residues = 0;
    }
    --refCounts[this];
    if (refCounts[this] == 0)
    {
        refCounts.erase(this);
    }
    moleculeDeleted(this);
#if defined(MEMORY_CORRUPTION_DEBUG)
    if (DELETION_COUNT==magic_delete)
        printf("Set magic_delete in debugger and set breakpoint here.\n");
    DELETED_MOLECULES[this] = DELETION_COUNT++;
#endif
}

ResidueListIterator MIMoleculeBase::residuesBegin()
{
    return Residue::getIterator(residues);
}

ResidueListIterator MIMoleculeBase::residuesEnd()
{
    return Residue::getIterator(NULL);
}

ResidueListIterator MIMoleculeBase::symmResiduesBegin()
{
    return Residue::getIterator(SymmResidues);
}

ResidueListIterator MIMoleculeBase::symmResiduesEnd()
{
    return Residue::getIterator(NULL);
}


MIAtom*MIMoleculeBase::GetAtom(int natom)
{
    MIAtom_const_iter atom, endAtom;
    ResidueListIterator res = residuesBegin();
    for (; res != residuesEnd(); ++res)
    {
        const MIAtomList &atoms = res->atoms();
        endAtom = atoms.end();
        for (atom = atoms.begin(); atom != endAtom; ++atom)
        {
            MIAtom *a = *atom;
            if (a->atomnumber() == natom)
            {
                return a;
            }
        }
    }

    for (res = symmResiduesBegin(); res != symmResiduesEnd(); ++res)
    {
        const MIAtomList &atoms = res->atoms();
        endAtom = atoms.end();
        for (atom = atoms.begin(); atom != endAtom; ++atom)
        {
            MIAtom *a = *atom;
            if (a->atomnumber() == natom)
            {
                return a;
            }
        }
    }
    return NULL;
}

// this keept showing up as a performance bottleneck, so it's been
// somewhat optimized now
bool MIMoleculeBase::alreadybonded(MIAtom *a1, MIAtom *a2)
{
    size_t s = bonds.size();
    for (unsigned int i = 0; i < s; i++)
    {
        Bond &bond = bonds[i];
        MIAtom *b1 = bond.getAtom1();
        MIAtom *b2 = bond.getAtom2();
        if ( (b1==a1 && b2==a2) || (b1==a2 && b2==a1))
        {
            return true;
        }
    }
    return false;
}

//FIXME: this imposes a requirement that we have a dictionary loaded!
//perhaps we could fall back to using BondLimit() instead of
//MIMolDictionary->AreBonded if no dictionary is defined
int MIMoleculeBase::Build(bool symmetryAtomsOnly)
{
    int i, j /*, nlines=0*/;
    int na = 0;
    MIAtom *atom;
    MIAtom *atom2;
    float x, y, z;
    double dx, dy, dz, dist = 1.90;
    Residue *res = residues;
    Residue *nter = residues;
    short nextid;
    nresidues = 0;
    Bond bond;

    long nth = 0;

    nlinks = 0;
    if (!symmetryAtomsOnly)
    {
        bonds.clear();
    }
    InitSeqPos();
    FixChains();
    FixAtomicNumbers();
    MIMolDictionary *dict = MIGetDictionary();

    for (int ir = ((symmetryAtomsOnly) ? 1 : 0); ir < 2; ir++)
    {
        // first set all atoms bonded flag off
        if (ir == 1)
        {
            res = SymmResidues;
            nter = SymmResidues;
        }
        else
        {
            res = residues;
            nter = residues;
        }

        while (Monomer::isValid(res))
        {
            nresidues++;
            for (i = 0; i < res->atomCount(); i++)
            {
                res->atom(i)->removeType(AtomType::BONDED);
            }
            res = res->next();
        }

        // reset res and look for matches within each residue.
        if (ir == 1)
        {
            res = SymmResidues;
        }
        else
        {
            res = residues;
        }
        //xmax = ymax = zmax = -9999;
        //xmin = ymin = zmin = 9999;
        while (res != NULL)
        {

            /*   gather atoms in residue and "pointtovu" them */
            for (i = 1; i < res->atomCount(); i++)
            {
                for (j = 0; j < i; j++)
                {
                    if (dict->AreBonded(res->type().c_str(), res->atom(i), res->atom(j), bond) != 0)
                    {
                        bond.setAtom1(res->atom(i));
                        bond.setAtom2(res->atom(j));
                        bond.getAtom1()->addType(AtomType::BONDED);
                        bond.getAtom2()->addType(AtomType::BONDED);
                        if (ir == 0)
                        {
                            bond.type = B_NORMAL;
                            bonds.push_back(bond);
                        }
                        else
                        {
                            bond.type = B_SYMM;
                            symmetryBonds.push_back(bond);
                        }
                    }
                }
            }
            // find the N atom of next res and look for bond-
            int is_dna = IsDna(*res);
            if (is_dna == 1)
            {
                strcpy(link_here, "O3*"); // dna link thru phosphate bond
                strcpy(link_next, "P");
            }
            else if (is_dna == 2)
            {
                strcpy(link_here, "O3\'"); // dna link thru phosphate bond
                strcpy(link_next, "P");
            }
            else
            {
                strcpy(link_here, "C");    // protein link thru peptide bond
                strcpy(link_next, "N");    //
            }
            if (res->next() != NULL && (atom = atom_from_name(link_next, *res->next())) != NULL
                && (atom2 = atom_from_name(link_here, *res)) != NULL)
            {
                x = atom->x();
                y = atom->y();
                z = atom->z();
                dx = atom2->x() -x;
                dy = atom2->y() -y;
                dz = atom2->z() -z;
                dist = BondLimit(atom->name())+BondLimit(atom2->name());
                if (atom2->altloc() != ' ')
                {
                    /* check to see if OK to bond */
                    if (atom->altloc() != ' ')
                    {
                        if (atom->altloc() != atom2->altloc())
                        {
                            goto skip;
                        }
                    }
                }
                if ( (dx < dist || dx > -dist) && (dy < dist || dy > -dist)
                     && hypot(dz, hypot(dx, dy)) < dist)
                {
                    if (!alreadybonded(atom, atom2))
                    {
                        bond.setAtom1(atom);
                        bond.setAtom2(atom2);
                        bond.getAtom1()->setType(AtomType::BONDED);
                        bond.getAtom2()->setType(AtomType::BONDED);
                        if (ir == 0)
                        {
                            bond.type = B_NORMAL;
                            bonds.push_back(bond);
                        }
                        else
                        {
                            bond.type = B_SYMM;
                            symmetryBonds.push_back(bond);
                        }
                        //if(!checkbonds()) return( bonds.size() );
                    }
                }
            }
            // if the residue is last in chain look for cyclic peptide
skip:   if (res->next() == NULL)
            {
                nextid = -1;
            }
            else
            {
                nextid = res->next()->chain_id();
            }
            if ((res->chain_id()) != nextid)
            {
                if ((atom = atom_from_name(link_here, *res)) != NULL
                    && (atom2 = atom_from_name(link_next, *nter)) != NULL)
                {
                    x = atom->x();
                    y = atom->y();
                    z = atom->z();
                    dx = atom2->x() -x;
                    dy = atom2->y() -y;
                    dz = atom2->z() -z;
                    if (atom2->altloc() != ' ')
                    {
                        /* check to see if OK to bond */
                        if (atom->altloc() != ' ')
                        {
                            if (atom->altloc() != atom2->altloc())
                            {
                                goto skip2;
                            }
                        }
                    }
                    dist = BondLimit(atom->name())+BondLimit(atom2->name());
                    if ( (dx < dist || dx > -dist) && (dy < dist || dy > -dist)
                         && hypot(dz, hypot(dx, dy)) < dist)
                    {
                        if (!alreadybonded(atom, atom2))
                        {
                            bond.setAtom1(atom);
                            bond.setAtom2(atom2);
                            bond.getAtom1()->setType(AtomType::BONDED);
                            bond.getAtom2()->setType(AtomType::BONDED);
                            if (ir == 0)
                            {
                                bond.type = B_NORMAL;
                                bonds.push_back(bond);
                            }
                            else
                            {
                                bond.type = B_SYMM;
                                symmetryBonds.push_back(bond);
                            }
                            //if(!checkbonds()) return( bonds.size());
                        }
                    }
                }
skip2:      if (res->next())
                {
                    nter = res->next();
                }
            }
            // finally if only 1 atom and its a CA link it to next CA
            if (res->next() != NULL)
            {
                MIAtom *atom1;
                if (res->atomCount() <= 2 && (atom1 = atom_from_name("CA", *res)) != NULL
                    && (atom2 = atom_from_name("CA", *res->next())) != NULL
                    && (res->chain_id()) == (res->next()->chain_id()))
                {
                    if (!alreadybonded(res->atom(0), atom2))
                    {
                        bond.setAtom1(res->atom(0));
                        bond.setAtom2(atom2);
                        bond.getAtom1()->setType(AtomType::BONDED);
                        bond.getAtom2()->setType(AtomType::BONDED);
                        if (ir == 0)
                        {
                            bond.type = B_NORMAL;
                            bonds.push_back(bond);
                        }
                        else
                        {
                            bond.type = B_SYMM;
                            symmetryBonds.push_back(bond);
                        }
                    }
                }
            }
            nth++;

            res = res->next();
        }
        // Add the CONECT records to the list
        if (ir == 0)
        {
            for (i = 0; (unsigned int) i < connects.size(); i++)
            {
                Connect(connects[i]);
            }
        }
        // now find unbonded residues and issue a point
        if (ir == 1)
        {
            res = SymmResidues;
        }
        else
        {
            res = residues;
        }
        while (res != NULL)
        {
            for (i = 0; i < res->atomCount(); i++)
            {
                na++;
                if (!(res->atom(i)->type() & AtomType::BONDED))
                {
                    // add a point to list;
                    bond.setAtom1(res->atom(i));
                    bond.setAtom2(res->atom(i));
                    if (ir == 0)
                    {
                        bond.type = B_POINT;
                        bonds.push_back(bond);
                    }
                    else
                    {
                        bond.type = B_SYMM_POINT;
                        symmetryBonds.push_back(bond);
                    }
                }
            }
            res = res->next();
        }
    }
    if (!symmetryAtomsOnly)
    {
        moleculeChanged(this);
    }
    return (bonds.size());
}

bool MIMoleculeBase::BreakBond(MIAtom *a1, MIAtom *a2)
{
    unsigned int i = 0;
    int ret = 0;
    while (i < bonds.size())
    {
        //search for bond(s) with these atoms
        if ((bonds[i].getAtom1() == a1 && bonds[i].getAtom2() == a2)
            || (bonds[i].getAtom1() == a2 && bonds[i].getAtom2() == a1))
        {
            bonds.erase(bonds.begin()+i);
            ret++;
        }
        else
        {
            i++;
        }
    }
    i = 0;
    while (i < hbonds.size())
    {
        //search for bond(s) with these atoms
        if ((hbonds[i].getAtom1() == a1 && hbonds[i].getAtom2() == a2)
            || (hbonds[i].getAtom1() == a2 && hbonds[i].getAtom2() == a1))
        {
            hbonds.erase(bonds.begin()+i);
            ret++;
        }
        else
        {
            i++;
        }
    }
    i = 0;
    while (i < connects.size())
    {
        //search for bond(s) with these atoms
        if ((connects[i].getAtom1() == a1 && connects[i].getAtom2() == a2)
            || (connects[i].getAtom1() == a2 && connects[i].getAtom2() == a1))
        {
            connects.erase(connects.begin()+i);
            SetCoordsChanged(true);
            ret++;
        }
        else
        {
            i++;
        }
    }
    return ret != 0;
}

bool MIMoleculeBase::AddBond(MIAtom *a1, MIAtom *a2)
{
    Bond bond;
    if ((a1 == NULL) || (a2 == NULL))
    {
        return false;
    }
    bond.setAtom1(a1);
    bond.setAtom2(a2);
    bond.type = B_CONNECT;
    bonds.push_back(bond);
    connects.push_back(bond);
    SetCoordsChanged(true);
    return (true);
}

void MIMoleculeBase::FreeBonds()
{
    nlinks = 0;
    bonds.clear();
}

void MIMoleculeBase::Translate(float x, float y, float z, MIAtomList *atoms)
{
    MIAtom *atom;
    unsigned int i;
    for (i = 0; i < atoms->size(); i++)
    {
        atom = (MIAtom*) (*atoms)[i];
        atom->translate(x, y, z);
    }
    SetCoordsChanged(true);
}

void MIMoleculeBase::SymmLink()
{
    MIAtom *a1, *a2, *a3;
    Residue *firstres = SymmResidues, *nextres;
    Residue *start = firstres;
    Residue *end = NULL;
    char linkname[10];
    double maxdist;
    Residue *res = firstres;
    bool dna;
    Bond bond;
    int inrange = 0;

    symmetryBonds.clear();

    while (Monomer::isValid(res) && (nextres = res->next()) != NULL)
    {
        if (!inrange)
        {
            if (start == res)
            {
                if (res->next() != end)
                {
                    inrange = 1;
                    firstres = res;
                }
            }
        }
        if (res->chain_id() != nextres->chain_id())
        {
            nextres = firstres;
            firstres = res->next();
        }

        dna = IsDna(*res) != 0;
        if (dna)
        {
            strcpy(linkname, "P");
            maxdist = 8.0;
        }
        else
        {
            strcpy(linkname, "CA");
            maxdist = 4.5;
        }
        if (inrange && (a1 = atom_from_name(linkname, *res)) != NULL
            && (a2 = atom_from_name(linkname, *nextres)) != NULL)
        {
            if (!linked(a1, a2) && AtomDist(*a1, *a2) < maxdist)
            {
                // if already linked don't add it again
                bond.setAtom1(a1);
                bond.setAtom2(a2);
                bond.type = B_LINK;
                symmetryBonds.push_back(bond);
                nlinks++;
                //if(!checkbonds()) break;
                if (dna)
                {
                    if (!strcmp(res->type().c_str(), "A") || !strcmp(res->type().c_str(), "G"))
                    {
                        strcpy(linkname, "N1");
                    }
                    else
                    {
                        strcpy(linkname, "N3");
                    }
                    a3 = atom_from_name(linkname, *res);
                    if (a1 && a3 && !linked(a1, a3))
                    {
                        bond.setAtom1(a1);
                        bond.setAtom2(a3);
                        bond.type = B_LINK;
                        symmetryBonds.push_back(bond);
                        //a3->color = getcolor(res->next(),a3,docolor,color,colormethod,atoms);
                        //setradiustype(a3,radiustype,atoms);
                        nlinks++;
                        //if(!checkbonds())break;
                    }
                }
            }
            //a1->color = getcolor(res,a1,docolor,color,colormethod,atoms);
            //a2->color = getcolor(res->next(),a2,docolor,color,colormethod,atoms);
            //setradiustype(a1,radiustype,atoms);
            //setradiustype(a2,radiustype,atoms);
        }
        if (inrange)
        {
            if (end == res || end == res->next())
            {
                inrange = 0;
            }
        }
        res = res->next();
    }
}

bool MIMoleculeBase::linked(MIAtom *a1, MIAtom *a2)
{
    if (nlinks <= 0)
    {
        return false;
    }
    for (unsigned int i = 0; i < bonds.size(); i++)
    {
        if ( (bonds[i].type == B_LINK)
             && ((bonds[i].getAtom1() == a1 && bonds[i].getAtom2() == a2)
                 || (bonds[i].getAtom1() == a2 && bonds[i].getAtom2() == a1)))
        {
            return true;
        }
    }
    return false;
}

bool MIMoleculeBase::linked(MIAtom *a1)
{
    if (nlinks <= 0)
    {
        return false;
    }
    for (unsigned int i = 0; i < bonds.size(); i++)
    {
        if ( (bonds[i].type == B_LINK)
             && ((bonds[i].getAtom1() == a1)
                 || (bonds[i].getAtom2() == a1)))
        {
            return true;
        }
    }
    return false;
}

void MIMoleculeBase::BuildLinks()
{
    if (nlinks > 0 && (unsigned int) nlinks <= bonds.size())
    {
        // since the bonds get sorted the only way to
        // get rid of the old links is to rebuild
        ClearLinks();
    }
    nlinks = BuildCALinks(bonds, residues);
}

void MIMoleculeBase::FixChains()
{

    Residue *reslist = residues;

    if (reslist == NULL)
    {
        return;
    }
    // Reset linkage_types
    AminoOrNucleic(reslist, true);

    Residue *res = reslist;
    if (!Monomer::isValid(res))
    {
        return;
    }

    bool breakByDiscontinuity = true;
    bool breakByNonpeptide = false;
    updateFixChainOptions(&breakByDiscontinuity, &breakByNonpeptide);

    int chainno = 0;
    int state = 0;
    Residue *prevRes = NULL;
    int prevState = -1;
    do
    {
        if (!Monomer::isValid(res))
        {
            state = -1;
        }
        else
        {
            if (res->linkage_type() & PEPTIDE)
            {
                state = PEPTIDE;
            }
            else if (res->linkage_type() & NUCLEIC)
            {
                state = NUCLEIC;
            }
            else
            {
                state = 0;
            }
        }
        if (state != prevState)
        {
            if (breakByNonpeptide)
            {
                ++chainno;
            }
            if (prevState == PEPTIDE)
            {
                prevRes->add_linkage_type(CTERMINUS);
            }
            else if (prevState == NUCLEIC)
            {
                prevRes->add_linkage_type(LAST);
            }
            if (state == PEPTIDE)
            {
                res->add_linkage_type(NTERMINUS);
            }
            else if (state == NUCLEIC)
            {
                res->add_linkage_type(FIRST);
            }
        }
        else
        {
            if (state == PEPTIDE)
            {
                MIAtom *CA1 = atom_from_name("CA", *prevRes);
                MIAtom *CA2 = atom_from_name("CA", *res);
                if (CA1 != NULL && CA2 != NULL)
                {
                    if (AtomDist(*CA1, *CA2) >= 5.0)
                    {
                        prevRes->add_linkage_type(CTERMINUS);
                        res->add_linkage_type(NTERMINUS);
                        if (breakByDiscontinuity)
                        {
                            ++chainno;
                        }
                    }
                    else
                    {
                        res->add_linkage_type(MIDDLE);
                    }
                }
            }
            else if (state == NUCLEIC)
            {
                res->add_linkage_type(MIDDLE);
            }
        }
        if (state != -1)
        {
            res->set_chain_id((res->chain_id()&255) + 256*chainno);

            prevRes = res;
            prevState = state;
            res = res->next();
        }
    } while (state != -1);

}

void MIMoleculeBase::ClearSymmList()
{
    /*
       unsigned int i=0;
       while(i< bonds.size()){
        if(bonds[i].type == B_SYMM || bonds[i].type == B_SYMM_POINT){
            bonds.erase(bonds.begin()+i);
        } else {
            i++;
        }
       }
       if(SymmResidues)FreeResidueList(SymmResidues);
     */
    symmetryBonds.clear();
    symmetryToBeCleared(this);
    PurgeSymmetryResidues(SymmResidues);
    SymmResidues = NULL;
}

void MIMoleculeBase::Connect(Bond &connect)
{
    Bond bond;
    if (!alreadybonded(connect.getAtom1(), connect.getAtom2()))
    {
        bond.setAtom1(connect.getAtom1());
        bond.setAtom2(connect.getAtom2());
        bond.getAtom1()->setType(AtomType::BONDED);
        bond.getAtom2()->setType(AtomType::BONDED);
        bond.type = B_CONNECT;
        bonds.push_back(bond);
    }
    SetCoordsChanged(true);
}

bool MIMoleculeBase::ReplaceMainChain(Residue *where, Residue *with, int nres)
{
    // replaces the main chain atoms of where with atoms in with for nres residues
    int i, j;
    Residue *res_to = where, *res_from = with;
    // just checking...
    for (i = 0; i < nres; i++)
    {
        if (!res_to || !res_from)
        {
            return false;
        }
        res_to = res_to->next();
        res_from = res_from->next();
    }
    res_to = where;
    res_from = with;
    MIAtom *CA_from, *CA_to;
    MIAtom *N_from, *O_from, *HN_from, *C_from;
    float dx, dy, dz;
    const char *name;

    // go through each residue, replace mainchain atoms - translate the others by CA-CA vector
    for (j = 0; j < nres; j++)
    {
        CA_from = atom_from_name("CA", *res_from);
        CA_to = atom_from_name("CA", *res_to);
        if (CA_from && CA_to)
        {
            dx = CA_to->x() - CA_from->x();
            dy = CA_to->y() - CA_from->y();
            dz = CA_to->z() - CA_from->z();
            N_from = atom_from_name("N", *res_from);
            C_from = atom_from_name("C", *res_from);
            // shelx name
            HN_from = atom_from_name("H0", *res_from);
            O_from = atom_from_name("O", *res_from);
            for (i = 0; i < res_to->atomCount(); i++)
            {
                name = res_to->atom(i)->name();
                if (strcmp(name, "N") == 0 && N_from)
                {
                    res_to->atom(i)->copyPosition(*N_from);
                }
                else
                if (strcmp(name, "C") == 0 && C_from)
                {
                    res_to->atom(i)->copyPosition(*C_from);
                }
                else
                if (strcmp(name, "HN") == 0 && HN_from)
                {
                    res_to->atom(i)->copyPosition(*HN_from);
                }
                else
                if (strcmp(name, "CA") == 0 && CA_from)
                {
                    res_to->atom(i)->copyPosition(*CA_from);
                }
                else
                if (strcmp(name, "O") == 0 && O_from)
                {
                    res_to->atom(i)->copyPosition(*O_from);
                }
                else
                {
                    res_to->atom(i)->translate(dx, dy, dz);
                }
            }
        }
        res_to = res_to->next();
        res_from = res_from->next();
    }
    SetCoordsChanged(true);
    return true;
}

void MIMoleculeBase::SetSecStr(Residue *res1, Residue *res2, char sec_str)
{
    Residue *start = NULL, *end = NULL, *res;

    // fiind the first residue
    res = residuesBegin();
    while (res != NULL)
    {
        if (res == res1)
        {
            start = res1;
            end = res2;
            break;
        }
        if (res == res2)
        {
            start = res2;
            end = res1;
            break;
        }
        res = res->next();
    }
    if (!start)
    {
        return;
    }
    res = start;
    while (Monomer::isValid(res))
    {
        res->setSecstr(sec_str);
        if (res == end)
        {
            break;
        }
        res = res->next();
    }
}

void MIMoleculeBase::FixAtomicNumbers()
{
    // fixes the atomicnumber field in a residue
    // sometimes this is a bit of a guess becuase it is not
    // clear what the atom type is by name alone!
    Residue *res = residuesBegin();
    int i;
    bool ispep, isnuc;
    std::string name;
    while (Monomer::isValid(res))
    {
        ispep = IsPeptide(*res) != 0;
        isnuc = IsNucleic(res) != 0;
        for (i = 0; i < res->atomCount(); i++)
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
                    res->atom(i)->setAtomicnumber(Atomic_Number_Nformat(name));
                }
            }
        }
        res = res->next();
    }
}

bool MIMoleculeBase::Revert(const char *pathname)
{
    FILE *fp = fopen(pathname, "r");
    if (fp == NULL)
    {
        std::string s("Molecule::Revert:Can't open file: ");
        s += pathname;
        Logger::message(s);
        return false;
    }
    std::vector<Bond> connct;
    Residue *newlist = LoadPDB(fp, &connct);
    fclose(fp);
    if (newlist == NULL)
    {
        return false;
    }

    DeleteAllResidues();
    ClearSymmList();

    residues = newlist;
    connects.clear();
    connects = connct;
    Build();
    return true;
}

bool MIMoleculeBase::BuildCB(Residue *test)
{
    if (!Monomer::isValid(test))
    {
        return false;
    }
    if (strcmp(test->type().c_str(), "MRK") != 0)
    {
        return false;
    }
    Residue *res = residuesBegin(), *prev = NULL, *next = test->next();
    while (res != NULL)
    {
        if (res->next() == test)
        {
            prev = res;
            break;
        }
        res = res->next();
    }
    if (!next || !test || !prev)
    {
        return false;
    }
    MIAtom *CA = atom_from_name("CA", *test);
    MIAtom *CAnext = atom_from_name("CA", *next);
    MIAtom *CAprev = atom_from_name("CA", *prev);
    if (!CA || !CAnext || !CAprev)
    {
        return false;
    }
    float mx, my, mz;
    mx = (CAnext->x() + CAprev->x())/2.0F;
    my = (CAnext->y() + CAprev->y())/2.0F;
    mz = (CAnext->z() + CAprev->z())/2.0F;
    float vx =  (CA->x() - mx);
    float vy =  (CA->y() - my);
    float vz =  (CA->z() - mz);
    float d = (float)sqrt(vx*vx + vy*vy + vz*vz);
    if (d == 0.0)
    {
        return false;
    }
    vx *= 1.54F/d;
    vy *= 1.54F/d;
    vz *= 1.54F/d;
    MIAtom *CB = new MIAtom;
    CB->copyShallow(*CA);
    CB->setPosition(CA->x() + vx,
                    CA->y() + vy,
                    CA->z() + vz);
    CB->setName("CB");

    // since this function clears out all atoms which are not CA (if any),
    // we need to properly clean those atoms to avoid a memory and/or reference leak
    MIAtomList deaders;
    for (size_t i = 0; i < test->atoms().size(); ++i)
    {
        if (test->atom(i) != CA)
            deaders.push_back(test->atom(i));
    }
    if (deaders.size())
    {
        atomsToBeDeleted(this, deaders);
        for (size_t i = 0; i < deaders.size(); ++i)
            PurgeAtom(deaders[i]);
        atomsDeleted(this);
    }
    test->clearAtoms();

    test->addAtom(CA);
    test->addAtom(CB);
    test->setType("VEC");
    SetCoordsChanged(true);
    return true;
}

static int greaterthan(Residue *res1, Residue *res2)
{
    int chainid1 = res1->chain_id()&255;
    int chainid2 = res2->chain_id()&255;
    if (chainid1 == ' ')
    {
        chainid1 = 254;
    }
    if (chainid2 == ' ')
    {
        chainid2 = 254;
    }
    if ((chainid1 - chainid2) != 0)
    {
        return (chainid1 > chainid2);
    }
    return (atoi(res1->name().c_str()) > atoi(res2->name().c_str()));
}

Residue *MIMoleculeBase::AddWater(float x, float y, float z, bool rebuild)
{
    // adds a water to end of model
    Residue *res = residuesBegin();
    while (res->next() != NULL)
    {
        res = res->next();
    }
    Residue *water = new Residue();
    res->insertResidue(water);
    water->setName(std::string(ftoa(1+(int)(atof(res->name().c_str())))));
    water->setType(std::string("HOH"));
    water->set_chain_id('X');
    water->setSecstr('U');
    MIAtom *atom = new MIAtom;
    water->addAtom(atom);
    atom->setPosition(x, y, z);
    atom->setBValue(15.0);
    atom->setOcc(1.0);
    atom->setName("O");
    if (MIGetColorSetter())
    {
        (*MIGetColorSetter())(atom);
    }
    atom->setAtomicnumber(8);
    atom->setType(MIAtom::MIGetAtomTypeFromName(atom->name()));
    atom->setAltloc(' ');
    if (rebuild)
    {
        Build();
    }
    SetCoordsChanged(true);
    modified = true;
    return water;
}

static Residue *findChainEnd(Residue *res)
{
    if (res == NULL)
    {
        return res;
    }
    int chainId = res->chain_id()&255;
    Residue *next = res->next();
    while (next != NULL)
    {
        if ((next->chain_id()&255) != chainId)
        {
            return res;
        }
        res = next;
        next = res->next();
    }
    return res;
}

void MIMoleculeBase::SortChains()
{

    FixChains();
    int iterations = 0;
    int done = true;
    do
    {
        Residue *chainStart = residues;
        Residue *chainEnd = findChainEnd(chainStart);
        Residue *prevChainEnd = NULL;
        done = true;
        while (Monomer::isValid(chainStart))
        {
            Residue *nextChainStart = chainEnd->next();
            if (nextChainStart == NULL)
            {
                break;
            }
            Residue *nextChainEnd = findChainEnd(nextChainStart);
            if (greaterthan(chainStart, nextChainStart))
            {
                // Swap chains
                nextChainStart->removeFromList(nextChainEnd);
                if (prevChainEnd != NULL)
                {
                    prevChainEnd->insertResidue(nextChainStart);
                }
                else
                {
                    residues = nextChainStart;
                    nextChainEnd->insertResidue(chainStart);
                }

                prevChainEnd = nextChainEnd;
                done = false;
            }
            else
            {
                prevChainEnd = chainEnd;
            }
            chainStart = nextChainStart;
            chainEnd = nextChainEnd;
        }
        ++iterations;
    } while (!done && iterations < 1000);
    FixChains();
    SetCoordsChanged(true);
}

bool MIMoleculeBase::Contains(Residue *test)
{
    if (!Monomer::isValid(test))
    {
        return false;
    }
    Residue *res = residuesBegin();
    while (Monomer::isValid(res))
    {
        if (test == res)
        {
            return true;
        }
        res = res->next();
    }
    return false;
}

void MIMoleculeBase::DeleteAtoms(MIAtomList atoms)
{
    atomsToBeDeleted(this, atoms);
    MIAtom_iter iter;
    for (iter = atoms.begin(); iter != atoms.end(); ++iter)
    {
        doDeleteAtom(*iter);
    }
    atomsDeleted(this);
}

void MIMoleculeBase::DeleteAtom(MIAtom *a)
{
    MIAtomList av;
    av.push_back(a);
    atomsToBeDeleted(this, av);
    doDeleteAtom(a);
    atomsDeleted(this);
}

void MIMoleculeBase::doDeleteAtom(MIAtom *a)
{
    Residue *res = residue_from_atom(residuesBegin(), a);
    if (!Monomer::isValid(res))
    {
        return;
    }
    if (res->atomCount() == 1)
    {
        DeleteRes(res);
        return;
    }

    if (res->removeAtom(a))
    {
        PurgeAtom(a);
        SetCoordsChanged(true);
    }
}

size_t MIMoleculeBase::LengthChain(unsigned short chain_id)
{
    size_t n = 0;
    Residue *res = residuesBegin();
    while (Monomer::isValid(res))
    {
        if (res->chain_id() == chain_id)
        {
            n++;
        }
        res = res->next();
    }
    return n;
}

bool MIMoleculeBase::ClearTorsion(void)
{
    int i;
    Residue *res = residues;
    while (res != NULL)
    {
        for (i = 0; i < res->atomCount(); i++)
        {
            res->atom(i)->removeType(AtomType::TORSIONATOM);
        }
        res = res->next();
    }
    Tatom1 = Tatom2 = NULL;
    return false;
}

bool MIMoleculeBase::GetTorsionValue(char *str, Residue *res)
{
    static bool found = false;
    static MIAtom *atom1 = NULL;
    static MIAtom *atom2 = NULL;
    static MIAtom *atom3 = NULL;
    static MIAtom *atom4 = NULL;
    static char type[11];
    float val;

    if (!MIGetDictionary())
    {
        return false;
    }

    if (atom2 == Tatom1 && atom3 == Tatom2 && found)   //this is stupid since we'll always search this way...
    {
        val = (float)CalcAtomTorsion(atom1, atom2, atom3, atom4);
        sprintf(str, "%s %0.3f", type, val);
        return true;
    }
    found = false;
    TORSION *tors = NULL;

#ifndef _WIN32
    const std::set<std::string> &tornames = MIGetDictionary()->TorsNames;
    std::set<std::string>::iterator i(tornames.begin()), e(tornames.end());
#else
    std::set<std::string> &tornames = MIGetDictionary()->TorsNames;
    std::set<std::string>::iterator i, e;
    i = tornames.begin();
    e = tornames.end();
#endif

    for (; i != e; i++)
    {
        tors = MIGetDictionary()->getTORSION(res, i->c_str(), res->next());
        if (tors == NULL)
        {
            continue;
        }
        if ( (tors->getAtom2() == Tatom1 || tors->atom3 == Tatom1 )
             && (tors->getAtom2() == Tatom2 || tors->atom3 == Tatom2 ))
        {
            atom1 = tors->getAtom1();
            atom2 = Tatom1;
            atom3 = Tatom2;
            atom4 = tors->atom4;
            val = (float)CalcAtomTorsion(atom1, atom2, atom3, atom4);
            strncpy(type, tors->type, 11);
            free(tors);
            sprintf(str, "%s %0.3f", type, val);
            found = true;
            return true;
        }
        free(tors);
        tors = NULL;
    }
    return false;
}

void MIMoleculeBase::RotateTorsion(float alpha)
{
    if (!Tatom1 || !Tatom2)
        return;

    float mat[4][3];
    initrotate(
        Tatom1->x(),
        Tatom1->y(),
        Tatom1->z(),
        (Tatom2->x()-Tatom1->x()),
        (Tatom2->y()-Tatom1->y()),
        (Tatom2->z()-Tatom1->z()),
        alpha, mat);
    for (Residue *res = residues; res; res = res->next())
    {
        for (int i = 0; i < res->atomCount(); ++i)
        {
            if (res->atom(i)->type() & AtomType::TORSIONATOM)
            {
                float x = res->atom(i)->x();
                float y = res->atom(i)->y();
                float z = res->atom(i)->z();
                rotate(&x, &y, &z, mat);
                res->atom(i)->setPosition(x, y, z);
            }
        }
    }
    SetCoordsChanged(true);
}

size_t MIMoleculeBase::SplitAtoms(MIAtomList &atoms, bool torsion_only)
{
    size_t nsplit = 0;
    Residue *res;
    char old_alt, new_alt;
    float new_occ, old_occ;
    for (size_t i = 0; i < atoms.size(); i++)
    {
        if (torsion_only && !(atoms[i]->type()&AtomType::TORSIONATOM))
        {
            continue;
        }
        MIAtom *a = atoms[i];
        res = residue_from_atom(residuesBegin(), a);
        if (!res)
        {
            continue;
        }
        if (a->altloc() == ' ')
        {
            old_alt = 'A';
        }
        else
        {
            old_alt = a->altloc();
        }
        new_alt = old_alt+1;
        old_occ = a->occ() *0.6F;
        new_occ = a->occ() *0.4F;
        // a bit of trickery here - the old atoms are marked 'B' and the new atoms
        // are marked 'A' so that the user can fit the "new" B part leaving the
        // "old" A part where it is.
        a->setOcc(new_occ);
        a->setAltloc(new_alt);
        nsplit++;
        MIAtom *b = new MIAtom;
        b->copyShallow(*a);
        res->addAtom(b);
        if (!torsion_only)
        {
            a->translate(0.03F, 0.03F, 0.03F); // slight offset to make new part visible
        }
        if (MIGetColorSetter())
        {
            (*MIGetColorSetter())(a, a->altloc());
        }
        b->setOcc(old_occ);
        b->setAltloc(old_alt);
        // clear flags in A part
        if (torsion_only)
        {
            b->removeType(~AtomType::TORSIONATOM);
        }
        b->removeType(~AtomType::FITATOM);
    }
    SetCoordsChanged(true);
    return nsplit;
}

int MIMoleculeBase::searchbonds(MIAtomList *atoms)
{
    unsigned int i;
    int n = 0;
    MIAtom *a1, *a2;
    for (i = 0; i < bonds.size(); i++)
    {
        a1 = bonds[i].getAtom1();
        a2 = bonds[i].getAtom2();
        if ((a1->type()&AtomType::TORSIONATOM || a2->type()&AtomType::TORSIONATOM)
            && find(atoms->begin(), atoms->end(), a1) != atoms->end()
            && find(atoms->begin(), atoms->end(), a2) != atoms->end())
        {
            if (a1 != Tatom1 && a2 != Tatom1)
            {
                if (!(a1->type()&AtomType::TORSIONATOM))
                {
                    a1->addType(AtomType::TORSIONATOM);
                    n++;
                }
                if (!(a2->type()&AtomType::TORSIONATOM))
                {
                    a2->addType(AtomType::TORSIONATOM);
                    n++;
                }
            }
        }
    }
    return n;
}

int MIMoleculeBase::SetupTorsion(MIAtom *a1, MIAtom *a2, MIAtomList *atoms)
{
    //find the atoms bound to a2 and collect them for torsion
    int i, n, nt = 0;
    Residue *res = residues;
    Tatom1 = NULL;
    Tatom2 = NULL;
    while (res != NULL)
    {
        for (i = 0; i < res->atomCount(); i++)
        {
            res->atom(i)->removeType(AtomType::TORSIONATOM);
            if (res->atom(i) == a1 && find(atoms->begin(), atoms->end(), a1) != atoms->end())
            {
                Tatom1 = a1;
            }
            if (res->atom(i) == a2 && find(atoms->begin(), atoms->end(), a1) != atoms->end())
            {
                Tatom2 = a2;
            }
        }
        res = res->next();
    }
    if (!Tatom1 || !Tatom2)
    {
        return 0;
    }
    Tatom2->addType(AtomType::TORSIONATOM);
    i = 0;
    while (i < 10000)
    {
        n = searchbonds(atoms);
        if (n == 0)
        {
            break;
        }
        nt += n;
        i++;
    }
    Tatom2->removeType(AtomType::TORSIONATOM);
    return nt;
}

int MIMoleculeBase::SetupTorsion(MIAtom *a1, MIAtom *a2)
{
    //find the atoms bound to a2 and collect them for torsion
    int i, n, nt = 0;
    Residue *res = residues;
    Tatom1 = NULL;
    Tatom2 = NULL;
    while (res != NULL)
    {
        for (i = 0; i < res->atomCount(); i++)
        {
            res->atom(i)->removeType(AtomType::TORSIONATOM);
            if (res->atom(i) == a1)
            {
                Tatom1 = a1;
            }
            if (res->atom(i) == a2)
            {
                Tatom2 = a2;
            }
        }
        res = res->next();
    }
    if (!Tatom1 || !Tatom2)
    {
        return 0;
    }
    Tatom2->addType(AtomType::TORSIONATOM);
    i = 0;
    while (i < 10000)
    {
        n = searchbonds();
        if (n == 0)
        {
            break;
        }
        nt += n;
        i++;
    }
    Tatom2->removeType(AtomType::TORSIONATOM);
    return nt;
}

int MIMoleculeBase::searchbonds()
{
    unsigned int i;
    int n = 0;
    for (i = 0; i < bonds.size(); i++)
    {
        if (bonds[i].getAtom1()->type()&AtomType::TORSIONATOM || bonds[i].getAtom2()->type()&AtomType::TORSIONATOM)
        {
            if (bonds[i].getAtom1() != Tatom1 && bonds[i].getAtom2() != Tatom1)
            {
                if (!(bonds[i].getAtom1()->type()&AtomType::TORSIONATOM))
                {
                    bonds[i].getAtom1()->addType(AtomType::TORSIONATOM);
                    n++;
                }
                if (!(bonds[i].getAtom2()->type()&AtomType::TORSIONATOM))
                {
                    bonds[i].getAtom2()->addType(AtomType::TORSIONATOM);
                    n++;
                }
            }
        }
    }
    return n;
}

void MIMoleculeBase::DeleteRes(Residue *oldres)
{
    std::vector<Residue*> residues;
    residues.push_back(oldres);
    DeleteResidues(residues);
}

void MIMoleculeBase::DeleteResidues(std::vector<Residue*> residues)
{
    residuesToBeDeleted(this, residues);
    std::vector<Residue*>::iterator iter;
    for (iter = residues.begin(); iter != residues.end(); ++iter)
    {
        doDeleteRes(*iter);
    }
    residuesDeleted(this);
}

void MIMoleculeBase::DeleteAllResidues()
{
    std::vector<Residue*> r;
    Residue *i = residues;
    while (i != NULL)
    {
        r.push_back(i);
        i = i->next();
    }
    DeleteResidues(r);
}

void MIMoleculeBase::doDeleteRes(Residue *oldres)
{
    PurgeResidue(oldres);
    SetCoordsChanged(true);
}

void MIMoleculeBase::PurgeAllAtoms()
{
    bonds.clear();
    hbonds.clear();
    connects.clear();
}

void MIMoleculeBase::PurgeAtom(MIAtom *atom)
{
    PurgeReferences(atom);
    delete atom;
    SetCoordsChanged(true);
}

void MIMoleculeBase::PurgeSymmetryAtom(MIAtom *atom)
{
    delete atom;
    SetCoordsChanged(true);
}

//NOTE: the PurgeReferences functions just remove internal references to
//the given object, they don't delete the objects.  This is used in
//ReplaceRes, and as a delegate to do the work in PurgeResidue (which does
//actually do a deletion)
void MIMoleculeBase::PurgeReferences(Residue *oldres)
{
    for (int i = 0; i < oldres->atomCount(); i++)
    {
        PurgeReferences(oldres->atom(i));
    }
}

void MIMoleculeBase::PurgeReferences(MIAtom *atom)
{
    size_t j = 0;
    while (j < bonds.size())
    {
        if (bonds[j].getAtom1() == atom || bonds[j].getAtom2() == atom)
        {
            bonds.erase(bonds.begin()+j);
        }
        else
        {
            j++;
        }
    }
    j = 0;
    while (j < hbonds.size())
    {
        if (hbonds[j].getAtom1() == atom || hbonds[j].getAtom2() == atom)
        {
            hbonds.erase(hbonds.begin()+j);
        }
        else
        {
            j++;
        }
    }
    j = 0;
    while (j < connects.size())
    {
        //search for bond(s) with this atom
        if (connects[j].getAtom1() == atom || connects[j].getAtom2() == atom)
        {
            connects.erase(connects.begin()+j);
        }
        else
        {
            j++;
        }
    }
}

void MIMoleculeBase::PurgeResidue(Residue *oldres)
{
    Residue *res = residues;
    SetCoordsChanged(true);
    while (res != NULL)
    {
        if (res == oldres)
        {
            if (res == residues)
            {
                residues = res->next();
            }
            res->removeFromList();
            PurgeReferences(res);
            delete res;
            nresidues  -= 1;
            return;
        }
        res = res->next();
    }
}

void MIMoleculeBase::PurgeSymmetryResidues(Residue *res)
{
    int i;
    Residue *oldres;
    while (res != NULL)
    {
        for (i = 0; i < res->atomCount(); i++)
        {
            PurgeSymmetryAtom(res->atom(i));
        }
        oldres = res;
        res = res->next();
        oldres->clearAtoms();
        delete oldres;
    }
    res = NULL;
}

Residue*MIMoleculeBase::InsertRes(Residue *atres, const Residue *dictres, int where, unsigned short chain_id)
{
    Residue *result = doInsertRes(atres, dictres, where, chain_id);
    moleculeChanged(this);
    return result;
}

void MIMoleculeBase::InsertResidues(Residue *atResidue, Residue *residues, int where, unsigned short chain_id)
{
    if (residues != NULL)
    {
        Residue *res = doInsertRes(atResidue, residues, where, chain_id);
        residues = residues->next();
        while (residues != NULL)
        {
            res = doInsertRes(res, residues, 0, chain_id);
            residues = residues->next();
        }
    }
    moleculeChanged(this);
}

Residue*MIMoleculeBase::doInsertRes(Residue *atres, const Residue *dictres, int where, unsigned short chain_id)
{
    // values for where:
    // 0 - After atres
    // 1 - before atres
    // 2 - at head of model
    // 3 - at end of model
    if (atres == NULL && where < 2)
    {
        return NULL;
    }
    if (!dictres)
    {
        return NULL;
    }
    if (where < 0 || where > 3)
    {
        return NULL;
    }
    if (dictres->atomCount() < 1)
    {
        return NULL;
    }
    //if(natoms == 0)return NULL;

    Residue *res = residues, *prev_res = NULL, *new_res;
    if (where == 2)
    {
        atres = residues;
        where = 1;
    }
    if (where == 3)
    {
        if (res)
        {
            while (Monomer::isValid(res->next()) && !IsWater(res->next()))
            {
                res = res->next();
            }
        }
        atres = res;
        where = 0;
    }
    // allocate more space to copy dictres into
    //long insertat = natoms-1;
    new_res = new Residue(*dictres);
    if (!new_res)
    {
        Logger::message("Disaster! Unable to allocate one more residue in MIMoleculeBase::InsertRes! Save and restart!");
        return NULL;
    }

    new_res->setSecstr('U');
    new_res->setName1(singleletter(new_res->type().c_str()));
    new_res->clearPrefBonds();
    new_res->clearPrefAngles();

    if (chain_id == 0 && atres != NULL)
    {
        chain_id = atres->chain_id();
    }
    if (chain_id != 0)
    {
        new_res->set_chain_id(chain_id);
    }
    for (int i = 0; i < new_res->atomCount(); i++)
    {
        if (MIGetColorSetter())
        {
            (*MIGetColorSetter())(new_res->atom(i));
        }
    }

    if (atres == NULL && residues == NULL)
    {
        residues = new_res;
    }
    else
    {

        // find the insertion point
        res = residues;
        while (Monomer::isValid(res))
        {
            if (res == atres)
            {
                // insert residue
                if (where == 0)
                {
                    // insert after
                    atres->insertResidue(new_res);
                    new_res->set_linkage_type(MIDDLE);
                    if (new_res->next() == NULL)
                    {
                        new_res->set_linkage_type(CTERMINUS);
                    }
                }
                else
                {
                    // insert before
                    // check to see if at start of model
                    if (prev_res == NULL)
                    {
                        residues->set_linkage_type(MIDDLE);
                        new_res->insertResidue(residues);
                        residues = new_res;
                        new_res->set_linkage_type(NTERMINUS);
                    }
                    else
                    {
                        prev_res->insertResidue(new_res);
                    }
                }
                break;
            }
            prev_res = res;
            res = res->next();
        }
    }
    nresidues  += 1;
    SetCoordsChanged(true);
    return new_res;
}

void MIMoleculeBase::ReplaceRes(Residue *oldres, Residue *dictres)
{
    MIAtom N, O, C, CA, H;
    int i;
    Residue *res;
    if (!dictres || !oldres)
    {
        return;
    }
    // make sure residue is in this molecule
    for (res = residues; res; res = res->next())
    {
        if (res == oldres)
        {
            break;
        }
    }
    if (!res)
    {
        return;
    }
    //long oldnatoms = natoms;
    int oldispep = IsPeptide(*oldres);
    int dictispep = IsPeptide(*dictres);
    N.setMass(0);
    H.setMass(0);
    C.setMass(0);
    O.setMass(0);
    CA.setMass(0);

    if (oldispep && dictispep)
    {
        //save the old mainchain coordinates by saving copies of relevant atoms
        for (i = 0; i < oldres->atomCount(); i++)
        {
            if (!strcmp("N", oldres->atom(i)->name()))
            {
                N.copyShallow(*oldres->atom(i));
                N.setMass(1);
            }
            if (!strcmp("H", oldres->atom(i)->name()))
            {
                H.copyShallow(*oldres->atom(i));
                H.setMass(1);
            }
            else if (!strcmp("HN", oldres->atom(i)->name()))
            {
                H.copyShallow(*oldres->atom(i));
                H.setMass(1);
            }
            if (!strcmp("C", oldres->atom(i)->name()))
            {
                C.copyShallow(*oldres->atom(i));
                C.setMass(1);
            }
            if (!strcmp("O", oldres->atom(i)->name()))
            {
                O.copyShallow(*oldres->atom(i));
                O.setMass(1);
            }
            if (!strcmp("CA", oldres->atom(i)->name()))
            {
                CA.copyShallow(*oldres->atom(i));
                CA.setMass(1);
            }
        }
    }

    atomsToBeDeleted(this, oldres->atoms());
    PurgeReferences(oldres);
    for (size_t i = 0; i < oldres->atoms().size(); ++i)
        delete oldres->atom(i);
    oldres->clearAtoms(); //it's ok, we just deleted the old atoms:

    oldres->setType(dictres->type());
    oldres->setConfomer(dictres->confomer());
    oldres->setName1(dictres->name1());

    //restore old mainchain coordinates from saving copies of relevant atoms
    for (i = 0; i < dictres->atomCount(); i++)
    {
        MIAtom *newAtom = new MIAtom;
        newAtom->copyShallow(*dictres->atom(i));
        oldres->addAtom(newAtom);
        if (!strcmp("N", oldres->atom(i)->name()) && N.mass())
        {
            oldres->atom(i)->copyPosition(N);
            N.setMass(0);
        }
        if (!strcmp("H", oldres->atom(i)->name()) && H.mass())
        {
            oldres->atom(i)->copyPosition(H);
            H.setMass(0);
        }
        if (!strcmp("HN", oldres->atom(i)->name()) && H.mass())
        {
            oldres->atom(i)->copyPosition(H);
            H.setMass(0);
        }
        if (!strcmp("C", oldres->atom(i)->name()) && C.mass())
        {
            oldres->atom(i)->copyPosition(C);
            C.setMass(0);
        }
        if (!strcmp("O", oldres->atom(i)->name()) && O.mass())
        {
            oldres->atom(i)->copyPosition(O);
            O.setMass(0);
        }
        if (!strcmp("CA", oldres->atom(i)->name()) && CA.mass())
        {
            oldres->atom(i)->copyPosition(CA);
            CA.setMass(0);
        }
    }
    for (i = 0; i < oldres->atomCount(); i++)
    {
        if (MIGetColorSetter())
        {
            (*MIGetColorSetter())(oldres->atom(i));
        }
    }
    SetCoordsChanged(true);
}

void MIMoleculeBase::BuildHBonds()
{

    Residue *res = residues;
    Residue *list;
    MIAtom *a1;
    MIAtom *a2;
    int i, j;

    hbonds.clear();

    while (res != NULL)
    {
        for (i = 0; i < res->atomCount(); i++)
        {
            a1 = res->atom(i);
            if (a1->isHidden())
            {
                continue;
            }
            if (a1->name()[0] != 'N' && a1->name()[0] != 'O')
            {
                continue;
            }
            list = res->next();
            while (list != NULL)
            {
                for (j = 0; j < list->atomCount(); j++)
                {
                    a2 = list->atom(j);
                    if (a2->isHidden())
                    {
                        continue;
                    }
                    if (hbondable(*a1, *a2, *res, *list))
                    {
                        if (!AddHBond(a1, a2))
                        {
                            return;
                        }
                    }
                }
                list = list->next();
            }
        }
        res = res->next();
    }
}

bool MIMoleculeBase::AddHBond(MIAtom *a1, MIAtom *a2)
{
    Bond hbond;
    hbond.setAtom1(a1);
    hbond.setAtom2(a2);
    hbonds.push_back(hbond);
    return true;
}

void MIMoleculeBase::ClearLinks()
{
    // go through bonds deleting any that are of the the type B_LINK
    size_t i = 0;
    while (i < bonds.size())
    {
        if (bonds[i].type == B_LINK)
        {
            bonds.erase(bonds.begin()+i);
        }
        else
        {
            i++;
        }
    }
    nlinks = 0;
}

bool MIMoleculeBase::SavePDBFile(const char *savepath)
{
    FILE *fp = NULL;

    std::vector<std::string> headers = FileHead;
    FixHeaders(headers);

    // open for writing
    fp = fopen(savepath, "w");
    if (!fp)
    {
        Logger::message("Save failed! Unable to open output file. Disk full?");
        return false;
    }
    Bond *edges = NULL;
    if (connects.size() > 0)
    {
        edges = &connects[0];
    }

    bool ret = SavePDB(fp, residues, edges, connects.size(), true, &headers, &FileTail);

    if (fp)
    {
        fclose(fp);
    }
    SetCoordsChanged(!ret);
    return ret;
}

void MIMoleculeBase::SecStrFromAngles()
{
    Residue *res = residues;
    Residue *prev = NULL;
    int i = 0;
    float angle1, angle2;
    MIAtom *po, *pc, *ro, *rc, *no, *nc;
    if (res == NULL)
    {
        return;
    }
    while (res != NULL)
    {
        if (prev == NULL || res->next() == NULL)
        {
            res->setSecstr('C');
        }
        else
        {
            pc = atom_from_name("C", *prev);
            rc = atom_from_name("C", *res);
            nc = atom_from_name("C", *res->next());
            po = atom_from_name("O", *prev);
            ro = atom_from_name("O", *res);
            no = atom_from_name("O", *res->next());
            if (po && pc && rc && ro)
            {
                angle1 = (float)CalcAtomTorsion(po, pc, rc, ro);
            }
            else
            {
                angle1 = 90.0;
            }
            if (no && nc && rc && ro)
            {
                angle2 = (float)CalcAtomTorsion(ro, rc, nc, no);
            }
            else
            {
                angle2 = 90.0F;
            }
            if (angle1 > 180.0F)
            {
                angle1 = angle1 - 360.0F;
            }
            if (angle2 > 180.0F)
            {
                angle2 = angle2 - 360.0F;
            }
            angle1 = (float)fabs(angle1);
            angle2 = (float)fabs(angle2);
            if (angle1 < 40.0F && angle2 < 60.0F)
            {
                res->setSecstr('H');
            }
            else if (angle1 > 120.0F && angle2 > 100.0F)
            {
                res->setSecstr('S');
            }
            else
            {
                res->setSecstr('C');
            }
        }
        prev = res;
        res = res->next();
        i++;
    }

    res = residues->next();
    prev = res;
    // fill in gaps and get rid of widows
    while (res != NULL && res->next() != NULL)
    {
        if ((prev->secstr() == res->next()->secstr()) && (prev->secstr() != res->secstr()))
        {
            res->setSecstr(prev->secstr());
        }
        prev = res;
        res = res->next();
    }

    // get rid of helices shorter than 4
    Residue *start = NULL;
    int length = 0;
    res = residues;
    while (res != NULL)
    {
        if (res->secstr() == 'H')
        {
            if (!start)
            {
                start = res;
            }
            length++;
        }
        else
        {
            if (start)
            {
                if (length > 0 && length < 4)
                {
                    // erase the helix
                    Residue *r = start;
                    for (int i = 0; i < length; i++)
                    {
                        r->setSecstr('C');
                        r = r->next();
                    }
                }
            }
            start = NULL;
            length = 0;
        }
        res = res->next();
    }
}

bool MIMoleculeBase::SavePDBFileVisible(const char *savepath)
{
    FILE *fp = NULL;
    // open for writing
    fp = fopen(savepath, "w");
    if (!fp)
    {
        Logger::message("Save failed! Unable to open output file. Disk full?");
        return false;
    }
    bool ret = SavePDB(fp, residues, NULL, -1);
    if (fp)
    {
        fclose(fp);
    }
    return ret;
}

//here

void MIMoleculeBase::DeleteChain(Residue *chain)
{
    std::vector<Residue*> residuesToDelete;
    Residue *res = chain;
    short chain_id = chain->chain_id();
    while (Monomer::isValid(res) && res->chain_id() == chain_id)
    {
        residuesToDelete.push_back(res);
        res = res->next();
    }
    DeleteResidues(residuesToDelete);
}

void MIMoleculeBase::setResidueNames(std::vector<Residue*> &residues, const std::string &name)
{
    int rnum = atoi(name.c_str());
    for (unsigned int i = 0; i<residues.size(); ++i)
        residues[i]->setName(std::string(ftoa(rnum++)));
    moleculeChanged(this);
}

void MIMoleculeBase::setChainId(Residue *chain, char c)
{
    int chain_id = chain->chain_id();
    int chainno = chain_id/256;
    Residue *res = chain;
    while (Monomer::isValid(res) && (chain_id == res->chain_id()))
    {
        res->set_chain_id((c&255) + chainno* 256);
        res = res->next();
    }
    moleculeChanged(this);
}

void MIMoleculeBase::renumberChain(Residue *chain, int n)
{
    int chain_id = chain->chain_id();
    Residue *res = chain;
    while (Monomer::isValid(res) && (chain_id == res->chain_id()))
    {
        res->setName(std::string(ftoa(n)));
        n++;
        res = res->next();
    }
    moleculeChanged(this);
}

//this was in header of old version
int MIMoleculeBase::getnresidues()
{
    nresidues = 0;
    for (Residue *r = residues; r != NULL; r = r->next())
    {
        nresidues++;
    }
    return nresidues;
}

// this was in Mlculseq, but is needed by base class
void MIMoleculeBase::InitSeqPos()
{
    int n = 0;
    for (Residue *res = residues; Monomer::isValid(res); res = res->next())
    {
        res->setSeqpos(n);
        res->setName1(singleletter(res->type().c_str()));
        n++;
    }
}

void MIMoleculeBase::Center(float &x, float &y, float &z)
{

    float xmin = std::numeric_limits<float>::max();
    float xmax = -std::numeric_limits<float>::max();
    float ymin = xmin;
    float ymax = xmax;
    float zmin = xmin;
    float zmax = xmax;

    for (ResidueListIterator res = residuesBegin(); res != residuesEnd(); ++res)
    {
        for (int i = 0; i < res->atomCount(); ++i)
        {
            float ix = res->atom(i)->x();
            float iy = res->atom(i)->y();
            float iz = res->atom(i)->z();
            xmax = std::max(ix, xmax);
            ymax = std::max(iy, ymax);
            zmax = std::max(iz, zmax);
            xmin = std::min(ix, xmin);
            ymin = std::min(iy, ymin);
            zmin = std::min(iz, zmin);
        }
    }

    x = (xmin+xmax)/2;
    y = (ymin+ymax)/2;
    z = (zmin+zmax)/2;
}


}
