#include <ui/Logger.h>
#include <stdio.h>
#include <string>

#include "mmCIF.h"
#include "CifParser.h"
#include "CifData.h"
#include "RefmacAtomTyper.h"
#include <util/system.h>
#include <chemlib/Monomer.h>
#include "mol_util.h"
#include "atom_util.h"

using namespace std;

namespace chemlib
{

MITorsionWritePrompt *TORSION_WRITE_PROMPT = NULL;

MITorsionWritePrompt *MIGetTorsionWritePrompt()
{
    return TORSION_WRITE_PROMPT;
}

void MISetTorsionWritePrompt(MITorsionWritePrompt *p)
{
    TORSION_WRITE_PROMPT = p;
}

mmCIF::mmCIF()
{
    _squirt_torsions = 2;
    atyper = NULL;
}

mmCIF::~mmCIF()
{
    if (atyper != NULL)
    {
        delete atyper;
    }
}

bool mmCIF::Write(FILE *fp, MIMolInfo &mol)
{
    int i;
    Residue *res = mol.res;
    std::string resname; // Used to convert residue type to proper format
    std::string line; // Used to setup the output for an individual line item
    std::string poly;
    map<MIAtom*, int> atom_map; // Used to go from atom address to atom index
    std::string atom1, atom2; // Used for writing bonds
    atyper = new chemlib::RefmacAtomTyper(*res, mol.bonds);
    //Setup cif header info
    if (fprintf(fp, "# File created by MIFit\n") <= 0)
    {
        return false;
    }


    //convert residue type
    resname = res->type();
    resname = MIToUpper(resname);
    MIStringReplace(resname, " ", "_");
    if (fprintf(fp, "data_comp_list\nloop_\n_chem_comp.id\n_chem_comp.three_letter_code\n_chem_comp.name\n_chem_comp.group\n_chem_comp.number_atoms_all\n_chem_comp.number_atoms_nh\n_chem_comp.desc_level\n") <= 0)
    {
        return false;
    }

    // Count the hydrogens
    int nonHydrogens = 0;
    for (i = 0; i < res->atomCount(); i++)
    {
        if (!MIAtom::MIIsHydrogen(res->atom(i)))
        {
            nonHydrogens++;
        }
    }
    if (IsDna(*res) || IsPeptide(*res))
    {
        poly = "polymer";
    }
    else
    {
        poly = "non-polymer";
    }

    // Print in order
    // Residue name (i.e. RDC)
    // Residue name (i.e. Radicicol)
    // Type (polymer, non-polymer)
    // # Atoms
    // # Hydrogens

    char buf[1024];
    sprintf(buf, "%s  %s  '%s'  %s  %d  %d  .\n",
            resname.c_str(),
            resname.c_str(),
            resname.c_str(),
            poly.c_str(),
            res->atomCount(),
            nonHydrogens);
    line = buf;
    if (fprintf(fp, "%s", line.c_str()) <= 0)
    {
        return false;
    }
    if (fprintf(fp, "data_comp_%s\n", resname.c_str()) <= 0)
    {
        return false;
    }

    if (!SquirtAtoms(fp, atom_map, res))
    {
        return false;
    }
    if (!SquirtTree(fp, resname, atom_map, res, mol.bonds))
    {
        return false;
    }
    if (!SquirtBonds(fp, resname, atom_map, mol.bonds))
    {
        return false;
    }
    if (!SquirtAngles(fp, resname, atom_map, mol.angles))
    {
        return false;
    }
    if (!SquirtTorsions(fp, resname, atom_map, mol.torsions))
    {
        return false;
    }
    if (!SquirtChirality(fp, resname, atom_map, res, mol.bonds))
    {
        return false;
    }
    if (!SquirtPlanes(fp, resname, atom_map, mol.planes))
    {
        return false;
    }
    return true;
}

bool mmCIF::SquirtAtoms(FILE *fp, map<MIAtom*, int> &atom_map,
                        Residue *res)
{
    int i;
    float charge;
    std::string resname; // Used to convert residue type to proper format
    std::string line; // Used to setup the output for an individual line item
    std::string atom;

    //convert residue type
    resname = res->type();
    resname = MIToUpper(resname);
    MIStringReplace(resname, " ", "_");

    // Output atom information
    if (fprintf(fp, "loop_\n_chem_comp_atom.comp_id\n_chem_comp_atom.atom_id\n_chem_comp_atom.type_symbol\n_chem_comp_atom.type_energy\n_chem_comp_atom.partial_charge\n_chem_comp_atom.x\n_chem_comp_atom.y\n_chem_comp_atom.z\n") <= 0)
    {
        return false;
    }

    resname = res->type();
    resname = MIToUpper(resname);
    for (i = 0; i < res->atomCount(); i++)
    {
        atom_map[res->atom(i)] = i; // Add to the map

        // Print in order:
        // Res Name
        // AtomTypeName Index (to ensure unique name of atom)
        // Symbol (believe this is just the Atomic Symbol,
        //		i.e. C for carbon, N for Nitrogen etc)
        // Energy
        // Charge
        // x, y, z coords
        AtomName(res->atom(i), atom_map, atom);
        //atom.Printf("%s%d", Atomic_Name(res->atom(i)->atomicnumber), i);

        //Take the charge from either the formal charge or the partial charge
        if (res->atom(i)->charge() == 0.0)
        {
            charge = (float) res->atom(i)->formal_charge();
        }
        else
        {
            charge = res->atom(i)->charge();
        }

        char buf[1024];
        sprintf(buf, " %s %5s %5s %5s % 0.3f  % 0.3f  % 0.3f  % 0.3f\n",
                resname.c_str(),
                atom.c_str(),
                Atomic_Name(res->atom(i)->atomicnumber()),
                atyper->AtomType(res->atom(i)),
                charge,
                res->atom(i)->x(), res->atom(i)->y(), res->atom(i)->z());
        if (fprintf(fp, "%s", buf) <= 0)
        {
            return false;
        }
    }

    return true;
}

bool orderbysize(const Bond *e1, const Bond *e2)
{
    if (e1->ideal_length < e2->ideal_length)
    {
        return true;
    }
    return false;
}

void mmCIF::GetBonds(MIAtom *atom, vector<Bond*> &bondlist,
                     vector<Bond> &bonds)
{
    vector<Bond>::iterator o, r;
    o = bonds.begin();
    r = bonds.end();
    while (o != r)
    {
        if (o->getAtom1() == atom || o->getAtom2() == atom)
        {
            bondlist.push_back(&*o);
        }
        o++;
    }
}

void mmCIF::SquirtTreeFindNext(MIAtom *prev, MIAtom *curr, MIAtom **next,
                               map<MIAtom*, bool> &alreadyprinted, vector<Bond*> &sorted_neighbors,
                               map<Bond*, bool> &bond_map, vector<Bond*>::iterator &i)
{
    //Finding the next atom Rules:
    // 1) The next atom != to the prev atom
    // 2) The next atom isn't being printed by a previous recursion
    vector<Bond*>::iterator e;

    i = sorted_neighbors.begin();
    e = sorted_neighbors.end();
    while (i != e)
    {
        //It is the printing of a PREV that determins that an edge has been spit
        if ((*i)->getAtom1() == prev || (*i)->getAtom2() == prev)
        {
            bond_map.insert(pair<Bond*, bool>(*i, true));
            i++;
            continue;
        }
        if ((*i)->getAtom1() == curr
            && alreadyprinted.find((*i)->getAtom2()) == alreadyprinted.end())
        {
            *next = (*i)->getAtom2();
            break;
        }
        else if ((*i)->getAtom2() == curr
                 && alreadyprinted.find((*i)->getAtom1()) == alreadyprinted.end())
        {
            *next = (*i)->getAtom1();
            break;
        }
        i++;
    }
}

bool mmCIF::SquirtTreeAtom(FILE *fp, std::string &resname, MIAtom *patom,
                           MIAtom *catom, MIAtom *natom, map<MIAtom*, int> &atom_map,
                           bool start = false)
{
    std::string atom_curr, atom_prev, atom_next;

    if (!AtomName(patom, atom_map, atom_prev))
    {
        return false;
    }
    if (!AtomName(catom, atom_map, atom_curr))
    {
        return false;
    }
    if (!AtomName(natom, atom_map, atom_next))
    {
        return false;
    }
    if ((!start && fprintf(fp, " %s %5s %5s %5s",
                           resname.c_str(),
                           atom_curr.c_str(),
                           atom_prev.c_str(),
                           atom_next.c_str()) <= 0)
        || (start && fprintf(fp, " %s %5s %5s %5s",
                             resname.c_str(),
                             atom_curr.c_str(),
                             "n/a",
                             atom_next.c_str()) <= 0))
    {
        return false;
    }
    return true;
}

bool mmCIF::SquirtTree(FILE *fp, std::string &resname,
                       map<MIAtom*, int> &atom_map, Residue *res, vector<Bond> &bonds)
{
    int j;
    vector<Bond>::iterator o, r;
    vector<Bond*>::iterator i, e;
    vector<Bond*> sorted_neighbors;
    map<MIAtom*, bool> alreadyprinted;
    map<Bond*, bool> bond_map;
    MIAtom *next = NULL, *curr = res->atom(0);

    //Spit out the header
    if (fprintf(fp, "loop_\n_chem_comp_tree.comp_id\n_chem_comp_tree.atom_id\n_chem_comp_tree.atom_back\n_chem_comp_tree.atom_forward\n_chem_comp_tree.connect_type\n") <= 0)
    {
        return false;
    }

    //Find a valid first atom
    for (j = 0; j < res->atomCount(); j++)
    {
        if (MIAtom::MIIsHydrogen(res->atom(j)))
        {
            continue;
        }
        // Not a hydrogen, hence valid as first atom
        break;
    }
    curr = res->atom(j);
    // Setup initial atom output
    GetBonds(curr, sorted_neighbors, bonds);
    std::sort(sorted_neighbors.begin(), sorted_neighbors.end(), orderbysize);
    alreadyprinted.insert(pair<MIAtom*, bool>(curr, true)); //I am being handled
    e = sorted_neighbors.end();
    SquirtTreeFindNext(NULL, curr, &next, alreadyprinted, sorted_neighbors,
                       bond_map, i);

    // Spit out initial atom
    SquirtTreeAtom(fp, resname, NULL, curr, next, atom_map, true);
    if (fprintf(fp, " START\n") <= 0)
    {
        return false;
    }

    // Get the recursion rolling
    if ((*i)->getAtom1() == curr
        && alreadyprinted.find((*i)->getAtom2()) == alreadyprinted.end())
    {
        if (!CreateMonomerTree(fp, resname, curr, (*i)->getAtom2(), atom_map,
                               alreadyprinted, bonds, bond_map, false))
        {
            return false;
        }
    }
    else if ((*i)->getAtom2() == curr
             && alreadyprinted.find((*i)->getAtom1()) == alreadyprinted.end())
    {
        if (!CreateMonomerTree(fp, resname, curr, (*i)->getAtom1(), atom_map,
                               alreadyprinted, bonds, bond_map, false))
        {
            return false;
        }
    }
    i++;
    // Go through the remaining
    while (i != e)
    {
        if ((*i)->getAtom1() == curr
            && alreadyprinted.find((*i)->getAtom2()) == alreadyprinted.end())
        {
            if (!CreateMonomerTree(fp, resname, curr, (*i)->getAtom2(), atom_map,
                                   alreadyprinted, bonds, bond_map, true))
            {
                return false;
            }
        }
        else if ((*i)->getAtom2() == curr
                 && alreadyprinted.find((*i)->getAtom1()) == alreadyprinted.end())
        {
            if (!CreateMonomerTree(fp, resname, curr, (*i)->getAtom1(), atom_map,
                                   alreadyprinted, bonds, bond_map, true))
            {
                return false;
            }
        }
        i++;
    }
    //Place end tag on last output from atoms
    if (fprintf(fp, " END\n") <= 0)
    {
        return false;
    }

    //Check for bonds not written, squirt them if needed
    o = bonds.begin();
    r = bonds.end();
    for (; o != r; o++)
    {
        if (bond_map.find(&*o) == bond_map.end()) //An unprinted bond
        {
            if (o->getAtom1() == curr)
            {
                SquirtTreeAtom(fp, resname, curr, o->getAtom2(), NULL, atom_map);
                if (fprintf(fp, " ADD\n") <= 0)
                {
                    return false;
                }
            }
            if (o->getAtom2() == curr)
            {
                SquirtTreeAtom(fp, resname, curr, o->getAtom1(), NULL, atom_map);
                if (fprintf(fp, " ADD\n") <= 0)
                {
                    return false;
                }
            }
        }
    }
    return true;
}

bool mmCIF::CreateMonomerTree(FILE *fp, std::string &resname, MIAtom *prev,
                              MIAtom *curr, map<MIAtom*, int> &atom_map,
                              map<MIAtom*, bool> &alreadyprinted, vector<Bond> &bonds,
                              map<Bond*, bool> &bond_map, bool squirtend)
{
    MIAtom *next;
    vector<Bond*> sorted_neighbors;
    vector<Bond*>::iterator i, e;

    next = NULL;
    //Do I squirt out a . at the end of line or START for the first one
    if (squirtend && fprintf(fp, " .\n") <= 0)
    {
        return false;
    }

    //Setup the neighbor list
    GetBonds(curr, sorted_neighbors, bonds);
    std::sort(sorted_neighbors.begin(), sorted_neighbors.end(), orderbysize);

    //Find my next pointer (I'll just crap out if there is no next)
    alreadyprinted.insert(pair<MIAtom*, bool>(curr, true)); //I am being handled
    e = sorted_neighbors.end();
    SquirtTreeFindNext(NULL, curr, &next, alreadyprinted, sorted_neighbors,
                       bond_map, i);

    //Print myself
    SquirtTreeAtom(fp, resname, prev, curr, next, atom_map);
    //Do recursion
    //Recursion Rules:
    // 1) Do not recurse on myself,
    // 2) Do not recurse on the previous atom
    // 3) Do not recurse on atoms being handled by a previous recursion
    while (i != e)
    {
        if ((*i)->getAtom1() == prev || (*i)->getAtom2() == prev)
        {
            bond_map.insert(pair<Bond*, bool>(*i, true));
            i++;
            continue;
        }
        if ((*i)->getAtom1() == curr
            && alreadyprinted.find((*i)->getAtom2()) == alreadyprinted.end())
        {
            if (!CreateMonomerTree(fp, resname, curr, (*i)->getAtom2(), atom_map, alreadyprinted, bonds, bond_map, true))
            {
                return false;
            }
        }
        else if ((*i)->getAtom2() == curr
                 && alreadyprinted.find((*i)->getAtom1()) == alreadyprinted.end())
        {
            if (!CreateMonomerTree(fp, resname, curr, (*i)->getAtom1(), atom_map, alreadyprinted, bonds, bond_map, true))
            {
                return false;
            }
        }
        i++;
    }
    return true;
}

bool mmCIF::SquirtBonds(FILE *fp, std::string &resname,
                        map<MIAtom*, int> &atom_map, vector<Bond> &bonds)
{
    unsigned int i;
    std::string line; // Used to setup the output for an individual line item
    std::string atom1, atom2; // Used for writing bonds

    if (bonds.size() == 0)
    {
        return true; // No need to write

    }
    if (fprintf(fp, "loop_\n_chem_comp_bond.comp_id\n_chem_comp_bond.atom_id_1\n_chem_comp_bond.atom_id_2\n_chem_comp_bond.type\n_chem_comp_bond.value_dist\n_chem_comp_bond.value_dist_esd\n") <= 0)
    {
        return false;
    }
    // Print in order:
    // Res Name
    // Atom1 -- Atom2 (i.e. atom1 bonded to atom2)
    // Bond Type (i.e. single, double, etc.)
    // Bond Distance
    // Bond Estimated Standard Dev. (?) (for now defaulted to 0.020
    for (i = 0; i < bonds.size(); i++)
    {
        if (!AtomName(bonds[i].getAtom1(), atom_map, atom1))
        {
            return false;
        }
        if (!AtomName(bonds[i].getAtom2(), atom_map, atom2))
        {
            return false;
        }
        char buf[1024];
        sprintf(buf, " %s %5s %5s %10s  % 0.3f  0.020\n",
                resname.c_str(),
                atom1.c_str(), atom2.c_str(),
                GetBondOrder(bonds[i].getOrder()).c_str(),
                bonds[i].ideal_length
                );
        if (fprintf(fp, "%s", buf) <= 0)
        {
            return false;
        }

    }
    return true;
}

bool mmCIF::SquirtAngles(FILE *fp, std::string &resname,
                         map<MIAtom*, int> &atom_map, vector<ANGLE> &angles)
{
    unsigned int i;
    std::string line; // Used to setup the output for an individual line item
    std::string atom1, atom2, atom3;

    if (angles.size() == 0)
    {
        return true;
    }

    if (fprintf(fp, "loop_\n_chem_comp_angle.comp_id\n_chem_comp_angle.atom_id_1\n_chem_comp_angle.atom_id_2\n_chem_comp_angle.atom_id_3\n_chem_comp_angle.value_angle\n_chem_comp_angle.value_angle_esd\n") <= 0)
    {
        return false;
    }
    // Print in order:
    // Res Name
    // Atom1 -- Atom2 -- Atom3
    // Angle
    // Angle Estimated Standard Dev. (?) (for now defaulted to 3.000)
    for (i = 0; i < angles.size(); i++)
    {
        if (!AtomName(angles[i].getAtom1(), atom_map, atom1))
        {
            return false;
        }
        if (!AtomName(angles[i].getAtom2(), atom_map, atom2))
        {
            return false;
        }
        if (!AtomName(angles[i].atom3, atom_map, atom3))
        {
            return false;
        }
        char buf[1024];
        sprintf(buf, " %s %5s %5s %5s % 0.3f 3.000\n",
                resname.c_str(),
                atom1.c_str(), atom2.c_str(), atom3.c_str(),
                CalcAtomAngle(*angles[i].getAtom1(), *angles[i].getAtom2(), *angles[i].atom3));
        if (fprintf(fp, "%s", buf) <= 0)
        {
            return false;
        }
    }

    return true;
}

bool mmCIF::SquirtTorsions(FILE *fp, std::string &resname, map<MIAtom*, int> &atom_map, vector<TORSION> &torsions)
{
    unsigned int i;
    std::string line; // Used to setup the output for an individual line item
    std::string atom;
    std::string atom1, atom2, atom3, atom4;

    if (torsions.size() == 0)
    {
        return true;
    }
    if (_squirt_torsions == 0)
    {
        return true;
    }
    if (_squirt_torsions == 2)
    {
        // Ask user
        if (MIGetTorsionWritePrompt() == NULL)
        {
            Logger::message("A MIFit development error occured (mmCIF::SquirtTorsions). Please report this.");
            return false;
        }
        else
        {
            if (!(*MIGetTorsionWritePrompt())())
            {
                return true;
            }
        }
    }
    if (fprintf(fp, "loop_\n_chem_comp_tor.comp_id\n_chem_comp_tor.id\n_chem_comp_tor.atom_id_1\n_chem_comp_tor.atom_id_2\n_chem_comp_tor.atom_id_3\n_chem_comp_tor.atom_id_4\n_chem_comp_tor.value_angle\n_chem_comp_tor.value_angle_esd\n") <= 0)
    {
        return false;
    }

    //Print in order:
    // Res Name
    // Torsion Name
    // Atom1 -- Atom2 -- Atom3 -- Atom4
    // Angle (?) //spit out the first torsion angle
    // Angle Estimated Standard Dev. (?)
    for (i = 0; i < torsions.size(); i++)
    {
        if (!AtomName(torsions[i].getAtom1(), atom_map, atom1))
        {
            return false;
        }
        if (!AtomName(torsions[i].getAtom2(), atom_map, atom2))
        {
            return false;
        }
        if (!AtomName(torsions[i].atom3, atom_map, atom3))
        {
            return false;
        }
        if (!AtomName(torsions[i].atom4, atom_map, atom4))
        {
            return false;
        }
        char buf[1024];
        sprintf(buf, " %s CONST_%03d %5s %5s %5s %5s % 0.3f % 0.3f\n",
                resname.c_str(),
                i,
                atom1.c_str(), atom2.c_str(), atom3.c_str(), atom4.c_str(),
                torsions[i].ideal[0],
                0.0);
        if (fprintf(fp, "%s", buf) <= 0)
        {
            return false;
        }
    }
    return true;
}

bool mmCIF::SquirtChirality(FILE *fp, std::string &resname,
                            map<MIAtom*, int> &atom_map, Residue *res, vector<Bond> &bonds)
{
    //TODO: See e-mail from Keith on Friday 2:05:09pm on 5-20-05
    int i, id;
    std::string line; // Used to setup the output for an individual line item
    std::string list;
    std::string atom;
    std::string nabor0, nabor1, nabor2;
    vector<MIAtom*> nabors;


    chemlib::GuessBondOrders(res, bonds);
    list = chemlib::FindChiralCenters(res, bonds, false);

    if (list.length() == 0)
    {
        return true;
    }

    // Output header for chirals section
    std::string header;
    header = "loop_\n";
    header += "_chem_comp_chir.comp_id\n";
    header += "_chem_comp_chir.id\n";
    header += "_chem_comp_chir.atom_id_centre\n";
    header += "_chem_comp_chir.atom_id_1\n";
    header += "_chem_comp_chir.atom_id_2\n";
    header += "_chem_comp_chir.atom_id_3\n";
    header += "_chem_comp_chir.volume_sign\n";
    if (fprintf(fp, header.c_str()) <= 0)
    {
        return false;
    }

    for (i = 0, id = 1; i < res->atomCount(); i++)
    {
        if (res->atom(i)->chiral_class() == CH_NONE)
        {
            continue;
        }
        // Print in order:
        // Res Name
        // AtomName of chiral center(to ensure unique name of atom)
        // AtomName of first neighbor
        // AtomName of second neighbor
        // AtomName of third neighbor
        // "both", so that the chirality will be inferred by REFMAC (per meeting 6/2005)
        AtomName(res->atom(i), atom_map, atom);
        nabors.clear();
        chemlib::GetNabors(res->atom(i), bonds, nabors);
        if (nabors.size() < 3)
        {
            continue;
        }
        AtomName(nabors[0], atom_map, nabor0);
        AtomName(nabors[1], atom_map, nabor1);
        AtomName(nabors[2], atom_map, nabor2);

        //atom.Printf("%s%d", Atomic_Name(res->atom(i)->atomicnumber), i);
        char buf[1024];
        sprintf(buf, " %s      chir_%02d %5s %5s %5s %5s      both\n",
                resname.c_str(),
                id++,
                atom.c_str(),
                nabor0.c_str(),
                nabor1.c_str(),
                nabor2.c_str());
        if (fprintf(fp, "%s", buf) <= 0)
        {
            return false;
        }
    }

    return true;
}

bool mmCIF::SquirtPlanes(FILE *fp, std::string &resname,
                         map<MIAtom*, int> &atom_map, vector<PLANE> &planes)
{
    unsigned int i;
    int j;
    std::string line; // Used to setup the output for an individual line item
    std::string atom;

    if (planes.size() == 0)
    {
        return true;
    }

    if (fprintf(fp, "loop_\n_chem_comp_plane_atom.comp_id\n_chem_comp_plane_atom.plane_id\n_chem_comp_plane_atom.atom_id\n_chem_comp_plane_atom.dist_esd\n") <= 0)
    {
        return false;
    }

    //Print in order:
    // Res Name
    // Plane Name
    // Atom
    // Distance Estimated Standard Dev. (?) (defaulted to 0.020 for now)
    for (i = 0; i < planes.size(); i++)
    {
        for (j = 0; j < planes[i].natoms; j++)
        {
            if (!AtomName(planes[i].atoms[j], atom_map, atom))
            {
                return false;
            }
            char buf[1024];
            sprintf(buf, " %s plan-%02d %5s 0.020\n",
                    resname.c_str(),
                    i+1,
                    atom.c_str());
            if (fprintf(fp, "%s", buf) <= 0)
            {
                return false;
            }
        }
    }
    return true;
}

std::string mmCIF::GetBondOrder(unsigned char c)
{
    switch (c)
    {
    case NORMALBOND:
        return "single";
    case SINGLEBOND:
        return "single";
    case DOUBLEBOND:
        return "double";
    case TRIPLEBOND:
        return "triple";
    case PARTIALDOUBLEBOND:
        return "aromatic";
    case HYDROGENBOND:
    case IONICBOND:
    case METALLIGANDBOND:
        return "UNKNOWN";
    default:
        return "single";
    }
}

bool mmCIF::AtomName(MIAtom *atom, map<MIAtom*, int>&, std::string &str)
{
    if (atom == NULL)
    {
        str = ".";
        return true;
    }
    /*map<MIAtom *, int>::iterator ami;
       ami = atom_map.find(atom);
       if(ami == atom_map.end()){
       Logger::message("Invalid atom entry in bond list (not found in map)");
       return false;
       }
       string.Printf("%s%d", Atomic_Name(atom->atomicnumber), ami->second);*/
    str = atom->name();
    return true;
}

void mmCIF::WriteTorsions(bool val)
{
    if (val)
    {
        _squirt_torsions = 1;
    }
    if (!val)
    {
        _squirt_torsions = 0;
    }
}

bool mmCIF::Read(FILE *fp, MIMolInfo &mol)
{
    CifParser parser(fp);
    CifDataBlock block;
    CifLoop loop;

    // clear current mol
    MIMolInfo foo;
    mol = foo;
    mol.res = new Residue();

    map<std::string, TORSDICT> torsion_map;
    map<std::string, PLANEDICT> plane_map;
    map<std::string, CHIRALDICT> chiral_map;
    //Loop over data blocks in file
    while (parser.GetNextBlock(block) )
    {


        if (block.FindLoop("chem_comp", loop))
        {
            SlurpHeader(loop, mol.res);
        }

        if (!block.FindLoop("chem_comp_atom", loop))
        {
            continue;
        }
        //Get atoms
        if (!SlurpAtoms(loop, mol.res))
        {
            continue;
        }

        //Create atom map
        map<std::string, MIAtom*> atom_map;
        for (int i = 0; i < mol.res->atomCount(); ++i)
        {
            atom_map[mol.res->atom(i)->name()] = mol.res->atom(i); // Add to the map
        }

        //Get bonds
        if (!block.FindLoop("chem_comp_bond", loop))
        {
            continue;
        }
        if (!SlurpBonds(loop, atom_map, mol.bonds))
        {
            continue;
        }
        //Get angles
        if (block.FindLoop("chemp_comp_angle", loop))
        {
            SlurpAngles(loop, atom_map, mol.angles, mol.res);
        }

        //Get torsions
        if (block.FindLoop("chem_comp_tor", loop))
        {
            SlurpTorsions(loop, torsion_map);
        }
        if (block.FindLoop("chem_comp_tor_value", loop))
        {
            SlurpTorsionValues(loop, torsion_map);
        }

        //Get planes
        if (block.FindLoop("chem_comp_plane", loop))
        {
            SlurpPlanes(loop, plane_map);
        }
        if (block.FindLoop("chem_comp_plane_atom", loop))
        {
            SlurpPlaneAtoms(loop, plane_map);
        }

        //Get chirals
        if (block.FindLoop("chem_comp_chir", loop))
        {
            SlurpChirals(loop, chiral_map);
        }
        if (block.FindLoop("chem_comp_chir_atom", loop))
        {
            SlurpChiralAtoms(loop, chiral_map);
        }
        //Transfer data from the (local) maps to the vector arguments
        if (!torsion_map.empty())
        {
            AppendFromPairs(mol.tordict, torsion_map.begin(), torsion_map.end());
        }

        if (!plane_map.empty())
        {
            AppendFromPairs(mol.planedict, plane_map.begin(), plane_map.end());
        }

        if (!chiral_map.empty())
        {
            AppendFromPairs(mol.chiralsdict, chiral_map.begin(), chiral_map.end());
        }

    } //End loop over data blocks

    return true;
}

bool mmCIF::SlurpAtoms(CifLoop &loop, Residue *res)
{

    CifTokenizer toker(loop._values);
    int nCols = loop._names.size();

    AtomKeyIndices index(loop._names);

    int col = 0;
    MIAtom *atom;
    MIAtomList atoms;
    std::string token, restype;
    bool hasX = false;
    bool hasY = false;
    bool hasZ = false;
    while (toker.GetToken(token) )
    {
        if (col == 0)                       //Allocate a new atom each
        {
            atom = new MIAtom;              //time we start thru the columns
            atom->setAtomnumber(atoms.size() + 1);
            hasX = false;
            hasY = false;
            hasZ = false;
        }

        if (token == "?" || token == ".") //Skip values that are unknown
        {
            col++;                          //or undefined
            continue;
        }

        if (col == index.xName)
        {
            atom->setName(token.c_str());
        }
        else if (col == index.xSymbol)
        {
            atom->setAtomicnumber(Atomic_Number_Nformat(token));
        }
        else if (col == index.xCharge)
        {
            atom->setCharge((float)atof(token.c_str()));
        }
        else if (col == index.xX)
        {
            atom->setX((float)atof(token.c_str()));
            hasX = true;
        }
        else if (col == index.xY)
        {
            atom->setY((float)atof(token.c_str()));
            hasY = true;
        }
        else if (col == index.xZ)
        {
            atom->setZ((float)atof(token.c_str()));
            hasZ = true;
        }
        else if (col == index.xRes)
        {
            restype = token;
        }

        col++;
        if (col == nCols)
        {
            if (hasX && hasY && hasZ)
            {
                atoms.push_back(atom);
            }
            else
            {
                delete atom;
            }
            col = 0;
        }
    }

    if (atoms.empty())
    {
        return false;
    }
    res->setType(restype);
    for (size_t i = 0; i < atoms.size(); ++i)
    {
        res->addAtom(atoms[i]);
        //Give the atom a color
        if (MIGetColorSetter())
        {
            (*MIGetColorSetter())(res->atom(i));
        }
    }

    return true;
}

bool mmCIF::SlurpBonds(CifLoop &loop, map<std::string, MIAtom*> &atom_map, vector<Bond> &bonds)
{
    CifTokenizer toker(loop._values);
    int nCols = loop._names.size();

    BondKeyIndices index(loop._names);

    std::string token, res, name1, name2;
    int col = 0;
    Bond bond;
    bool hasAtom1 = false;
    bool hasAtom2 = false;
    while (toker.GetToken(token))
    {
        if (col == 0)                       //Initialize the bond each
        {
            bond.Clear();
            hasAtom1 = false;
            hasAtom2 = false;
        }

        if (token == "?" || token == ".") //Skip values that are unknown
        {
            col++;                          //or undefined
            continue;
        }

        if (col == index.xRes)
        {
            res = token;
        }
        else if (col == index.xAtom1)
        {
            MIAtom *atom = atom_map[token];
            if (atom != NULL)
            {
                bond.setAtom1(atom);
                hasAtom1 = true;
            }
        }
        else if (col == index.xAtom2)
        {
            MIAtom *atom = atom_map[token];
            if (atom != NULL)
            {
                bond.setAtom2(atom);
                hasAtom2 = true;
            }
        }
        else if (col == index.xOrder)
        {
            bond.setOrder(DecodeCifBondOrder(token));
        }
        else if (col == index.xLength)
        {
            bond.ideal_length = (float)atof(token.c_str());
        }
        else if (col == index.xTolerance)
        {
            bond.tolerance = (float)atof(token.c_str());
        }

        col++;
        if (col == nCols)
        {
            if (hasAtom1 && hasAtom2)
            {
                bonds.push_back(bond);
            }
            col = 0;
        }
    }
    return true;
}

bool mmCIF::SlurpAngles(CifLoop &loop,
                        map<std::string, MIAtom*> &atom_map,
                        vector<ANGLE> &angles,
                        Residue *current_res)
{
    CifTokenizer toker(loop._values);
    int nCols = loop._names.size();

    AngleKeyIndices index(loop._names);

    std::string token, restype, name1, name2;
    int col = 0;
    ANGLE angle;
    while (toker.GetToken(token))
    {
        if (col == 0)                       //Initialize the angle each
        {
            angle.Clear();
        }

        if (token == "?" || token == ".") //Skip values that are unknown
        {
            col++;                          //or undefined
            continue;
        }

        if (col == index.xRes)
        {
            restype = token;
        }
        else if (col == index.xAtom1)
        {
            angle.setAtom1(atom_map[token]);
            if (angle.getAtom1() == NULL)
            {
                return false;
            }
        }
        else if (col == index.xAtom2)
        {
            angle.setAtom2(atom_map[token]);
            if (angle.getAtom2() == NULL)
            {
                return false;
            }
        }
        else if (col == index.xAtom3)
        {
            if ((angle.atom3 = atom_map[token]) == NULL)
            {
                return false;
            }
        }
        else if (col == index.xAngle)
        {
            angle.ideal_angle = (float)atof(token.c_str());
        }
        else if (col == index.xTolerance)
        {
            angle.tolerance = (float)atof(token.c_str());
        }

        col++;
        if (col == nCols)
        {
            if (restype == current_res->type())
            {
                angle.res = current_res;
                angles.push_back(angle);
            }
            col = 0;
        }
    }
    return true;
}

bool mmCIF::SlurpTorsions(CifLoop &loop, map<std::string, TORSDICT> &torsion_map)
{
    //bool mmCIF::SlurpTorsions(CifLoop &loop, vector<TORSDICT> &torsions) {
    CifTokenizer toker(loop._values);
    int nCols = loop._names.size();

    TorsionKeyIndices index(loop._names);

    std::string token, res, name1, name2, torsion_name;
    int col = 0;
    int period;
    bool found_angle, found_period, found_tolerance, found_name;
    float angle, tolerance;
    TORSDICT torsion;
    while (toker.GetToken(token))
    {
        if (col == 0)                       //Initialize the torsion each
        {
            torsion.Clear();
            found_angle = false;
            found_period = false;
            found_tolerance = false;
            found_name = false;
        }

        if (token == "?" || token == ".") //Skip values that are unknown
        {
            col++;                          //or undefined
            continue;
        }

        if (col == index.xRes)
        {
            strncpy(torsion.restype, token.c_str(), MAXNAME);
        }
        else if (col == index.xAtom1)
        {
            strncpy(torsion.name[0], token.c_str(), chemlib::MAXATOMNAME);
        }
        else if (col == index.xAtom2)
        {
            strncpy(torsion.name[1], token.c_str(), chemlib::MAXATOMNAME);
        }
        else if (col == index.xAtom3)
        {
            strncpy(torsion.name[2], token.c_str(), chemlib::MAXATOMNAME);
        }
        else if (col == index.xAtom4)
        {
            strncpy(torsion.name[3], token.c_str(), chemlib::MAXATOMNAME);
        }
        else if (col == index.xName)
        {
            torsion_name = token;
            strncpy(torsion.type, token.c_str(), 11);
            found_name = true;
        }
        else if (col == index.xAngle)
        {
            angle = (float)atof(token.c_str());
            found_angle = true;
        }
        else if (col == index.xTolerance)
        {
            tolerance = (float)atof(token.c_str());
            found_tolerance = true;
        }
        else if (col == index.xPeriod)
        {
            period = atoi(token.c_str());
            found_period = true;
        }

        col++;
        if (col == nCols)
        {
            if (found_angle && ((found_tolerance && tolerance == 0) //Is it an improper?
                                || (!found_period || period == 0)))
            {
                if (found_period && period == 2)                //Improper with 2 possibles
                {
                    torsion.ideal[0] = angle;
                    torsion.ideal[1] = (angle > 0) ? angle - 180.0F : angle + 180.0F;
                    torsion.nideal = 2;
                }
                else if (found_period && period == 3) //Improper with 3 possibles
                {
                    torsion.ideal[0] = angle;
                    torsion.ideal[1] = (angle <= 60) ? angle + 120.0F : angle - 240.0F;
                    torsion.ideal[2] = (angle > -60) ? angle - 120.0F : angle + 240.0F;
                    torsion.nideal = 3;
                }
                else
                {
                    torsion.ideal[0] = angle; //The usual case for an improper: 1 possible
                    torsion.nideal = 1;
                }
            }
            else
            {
                torsion.nideal = 0;
            }
            //				torsions.push_back(torsion);
            if (found_name)
            {
                torsion_map[torsion_name] = torsion;
            }
            col = 0;
        }
    }
    return true;
}

bool mmCIF::SlurpTorsionValues(CifLoop &loop,
                               map<std::string, TORSDICT> &torsion_map)
{
    CifTokenizer toker(loop._values);
    int nCols = loop._names.size();

    TorsionValueKeyIndices index(loop._names);

    map<std::string, TORSDICT>::iterator t;

    int col = 0;
    std::string token, res_name, tors_name;
    bool angle_found = false;
    float angle;
    while (toker.GetToken(token))
    {
        if (col == 0)
        {
            res_name.clear();
            tors_name.clear();
            angle_found = false;
        }
        if (col == index.xRes)
        {
            res_name = token;
        }
        else if (col == index.xName)
        {
            tors_name = token;
        }
        else if (col == index.xAngle)
        {
            angle = (float)atof(token.c_str());
        }

        col++;
        if (col == nCols)
        {
            t = torsion_map.find(tors_name);
            if (t != torsion_map.end()
                && res_name == t->second.restype)
            {
                t->second.AddAngle(angle);
            }
            col = 0;
        }
    }
    return true;
}

bool mmCIF::SlurpPlanes(CifLoop &loop,
                        map<std::string, PLANEDICT> &plane_map)
{

    CifTokenizer toker(loop._values);
    int nCols = loop._names.size();

    PlaneKeyIndices index(loop._names);

    PLANEDICT plane;
    //	map<std::string, PLANEDICT> plane_map;
    map<std::string, PLANEDICT>::iterator p;

    std::string token, res;
    std::string plane_name, atom_name, res_name;
    int col = 0;                                //Track where in the index rotation we are

    while (toker.GetToken(token))
    {
        if (col == 0)                       //Initialize the plane each
        {
            plane.Clear();
            plane_name.clear();
            atom_name.clear();
            res_name.clear();
        }

        if (token == "?" || token == ".") //Skip values that are unknown
        {
            col++;                          //or undefined
            continue;
        }

        if (col == index.xRes)
        {
            res_name = token;
        }
        else if (col == index.xName)
        {
            plane_name = token;
        }
        else if (col == index.xAtom)
        {
            atom_name = token;
        }
        else if (col == index.xAtomCount)
        {
            plane.natoms = atoi(token.c_str());
        }
        else if (col == index.xHvyAtomCount)        //Don't need to store num heavy atoms
        {
        }
        else if (col == index.xTolerance)           //We don't store tolerances yet
        {
        }

        col++;
        if (col == nCols)
        {
            if (plane_map.find(plane_name) != plane_map.end()
                && !atom_name.empty())
            {
                p = plane_map.find(plane_name);
                if (!p->second.AddAtom(atom_name.c_str()))
                {
                    return false;
                }
            }
            else if (plane_map.find(plane_name) == plane_map.end()
                     && !atom_name.empty())
            {
                plane.natoms = 1;
                strncpy(plane.restype, res_name.c_str(), MAXNAME);
                strncpy(plane.name[0], atom_name.c_str(), MAXNAME);
                plane_map[plane_name] = plane;
            }
            else if (plane_map.find(plane_name) == plane_map.end()
                     && atom_name.empty())
            {
                plane.natoms = 0;
                strncpy(plane.restype, res_name.c_str(), MAXNAME);
                plane_map[plane_name] = plane;
            }
            col = 0;
        }
    }

    //	map<std::string, PLANEDICT>::iterator pe = plane_map.end();
    //	p = plane_map.begin();
    //	while ( p != pe) {
    //		planes.push_back(p->second);
    //		++p;
    //	}
    return true;
}

bool mmCIF::SlurpPlaneAtoms(CifLoop &loop,
                            map<std::string, PLANEDICT> &plane_map)
{

    CifTokenizer toker(loop._values);
    int nCols = loop._names.size();

    PlaneAtomKeyIndices index(loop._names);

    map<std::string, PLANEDICT>::iterator p;
    std::string token, res;
    std::string plane_name, atom_name, res_name;
    int col = 0;                                //Track where in the index rotation we are

    while (toker.GetToken(token))
    {
        if (col == 0)
        {
            plane_name.clear();
            atom_name.clear();
            res_name.clear();
        }

        if (col == index.xRes)
        {
            res_name = token;
        }
        else if (col == index.xName)
        {
            plane_name = token;
        }
        else if (col == index.xAtom)
        {
            atom_name = token;
        }

        col++;
        if (col == nCols)
        {
            p = plane_map.find(plane_name);
            if (p == plane_map.end()
                && !res_name.empty()
                && !plane_name.empty()
                && !atom_name.empty())
            {
                PLANEDICT plane;
                plane.natoms = 1;
                strncpy(plane.restype, res_name.c_str(), MAXNAME);
                strncpy(plane.name[0], atom_name.c_str(), MAXNAME);
                plane_map[plane_name] = plane;
            }
            else if (!res_name.empty()
                     && !atom_name.empty())
            {
                if (!p->second.AddAtom(atom_name.c_str()))
                {
                    return false;
                }
            }
            col = 0;
        }
    }
    return true;
}

bool mmCIF::SlurpChirals(CifLoop &loop, map<std::string, CHIRALDICT> &chiral_map)
{

    CifTokenizer toker(loop._values);
    int nCols = loop._names.size();

    ChiralKeyIndices index(loop._names);
    CHIRALDICT chiral;
    int col = 0;                                //Track where in the index rotation we are

    std::string token, chiral_name;
    while (toker.GetToken(token))
    {
        if (col == 0)
        {
            chiral.Clear();
            chiral_name.clear();
        }

        if (col == index.xRes)
        {
            strncpy(chiral.restype, token.c_str(), MAXNAME);
        }
        else if (col == index.xName)
        {
            chiral_name = token;
        }
        else if (col == index.xCenter)
        {
            strncpy(chiral.center, token.c_str(), chemlib::MAXATOMNAME);
        }
        else if (col == index.xAtom1)
        {
            strncpy(chiral.name[0], token.c_str(), chemlib::MAXATOMNAME);
        }
        else if (col == index.xAtom2)
        {
            strncpy(chiral.name[1], token.c_str(), chemlib::MAXATOMNAME);
        }
        else if (col == index.xAtom3)
        {
            strncpy(chiral.name[2], token.c_str(), chemlib::MAXATOMNAME);
        }
        else if (col == index.xVolume)
        {
            float volume = (float)atof(token.c_str());
            if (token == "positive" || token == "positiv" || volume > 0)
            {
                chiral.order = CLOCKWISE;
            }
            else if (token == "negative" || token == "negativ" || volume < 0)
            {
                chiral.order = COUNTERCLOCKWISE;
            }
            else
            {
                chiral.order = 0;
            }
        }
        else if (col == index.xDegree)          //Not currently using these
        {
        }
        else if (col == index.xHvyDegree)
        {
        }
        else if (col == index.xConfig)          //Would allow the user to specify chirality
        {
        }                                       //as "R" or "S"

        col++;
        if (col == nCols)
        {
            chiral_map[chiral_name] = chiral;
            col = 0;
        }
    }
    return true;
}

bool mmCIF::SlurpChiralAtoms(CifLoop &loop, map<std::string, CHIRALDICT> &chiral_map)
{
    CifTokenizer toker(loop._values);
    int nCols = 0;

    ChiralAtomKeyIndices index(loop._names);

    map<std::string, CHIRALDICT>::iterator c;
    map<std::string, int> n_subs;

    int col = 0;
    std::string token, chiral_name, res_name, atom_name;
    while (toker.GetToken(token))
    {
        if (col == 0)
        {
            chiral_name.clear();
            res_name.clear();
        }

        if (col == index.xRes)
        {
            res_name = token;
        }
        else if (col == index.xName)
        {
            chiral_name = token;
        }
        else if (col == index.xAtom)
        {
            atom_name = token;
        }

        col++;
        if (col == nCols)
        {
            c = chiral_map.find(chiral_name);
            if (c != chiral_map.end()
                && res_name == c->second.restype)
            {
                c->second.AddAtom(atom_name.c_str(), n_subs[chiral_name]);
                n_subs[chiral_name]++;
            }
            col = 0;
        }
    }
    return true;
}

bool mmCIF::SlurpHeader(CifLoop &loop, Residue *res)
{

    CifTokenizer toker(loop._values);
    HeaderKeyIndices index(loop._names);

    int col = 0;
    std::string token;
    while (toker.GetToken(token))
    {
        if (col == index.xDescLevel)
        {
            res->setName1(token[0]);
            return true;
        }
        col++;
    }
    return false;
}

}
