#include <vector>
#include <fstream>
#include <chemlib/Monomer.h>
#include "Chiral.h"
#include "math_util.h"
#include "sequence_util.h"
#include "atom_util.h"
#include "Ligand.h"

namespace chemlib
{

/////////////////////////////////////////////////////////////////////////////
// Function:    SetCenter
// Purpose:		Specify the atom that represents the tetrahedral center
// Input:       Ptr to the chiral atom
// Output:      None
// Requires:
/////////////////////////////////////////////////////////////////////////////
void Chiral::SetCenter(MIAtom *atom)
{
    _center = atom;
}

/////////////////////////////////////////////////////////////////////////////
// Function:    SetOrder
// Purpose:		Specify the chirality--i.e. the order of substituents
// Input:       1 or 2
// Output:      None
// Requires:	Use of convention defined in Xguicryst.h:
//				1=COUNTERCLOCKWISE 2=CLOCKWISE
/////////////////////////////////////////////////////////////////////////////

void Chiral::SetOrder(int order)
{
    _order = order;
}

/////////////////////////////////////////////////////////////////////////////
// Function:    AddSub
// Purpose:		Allow clients to build the chiral definition by specifying the
//				sequence of substituent atoms
// Input:       Ptr to the next substituent atom in sequence
// Output:      None
// Requires:
/////////////////////////////////////////////////////////////////////////////

void Chiral::AddSub(MIAtom *atom)
{
    _subs.push_back(atom);
}

/////////////////////////////////////////////////////////////////////////////
// Function:    GetCenter
// Purpose:		Allow clients access to the identity of the chiral atom
// Input:       None
// Output:      Ptr to the chiral atom
// Requires:
/////////////////////////////////////////////////////////////////////////////
MIAtom*Chiral::GetCenter()
{
    return _center;
}

/////////////////////////////////////////////////////////////////////////////
// Function:    GetOrder
// Purpose:		Reports the chiral order (of substituents)
// Input:       None
// Output:      1 for counterclockwise, 2 for clockwise
// Requires:
/////////////////////////////////////////////////////////////////////////////
int Chiral::GetOrder()
{
    return _order;
}

/////////////////////////////////////////////////////////////////////////////
// Function:    GetSub
// Purpose:		Allow clients access to the identity of the chiral atom
// Input:       None
// Output:      Ptr to the chiral atom
// Requires:
/////////////////////////////////////////////////////////////////////////////

MIAtom*Chiral::GetSub(int i)
{
    return _subs[i];
}

/////////////////////////////////////////////////////////////////////////////
// Function:    Clear
// Purpose:		Reinitialize the object to having no atoms or order specified
// Input:       None
// Output:      None
// Requires:
/////////////////////////////////////////////////////////////////////////////

void Chiral::Clear()
{
    _order = 0;
    _center = 0;
    _subs.clear();
}

/////////////////////////////////////////////////////////////////////////////
// Function:    Card
// Purpose:		Writes the chiral data to a string
// Input:       None
// Output:      A string with a single line, having a sequence of atom names
//				starting with the chiral center, then the order
// Requires:
/////////////////////////////////////////////////////////////////////////////

void Chiral::Card(std::string &s)
{
    s.append(_center->name());
    s.append(" ");

    MIAtom_iter atm;
    for (atm = _subs.begin(); atm != _subs.end(); ++atm)
    {
        s.append((*atm)->name());
        s.append(" ");
    }

    char buf[5];
    sprintf(buf, "%4d", _order);
    s.append(buf);
}

/////////////////////////////////////////////////////////////////////////////
// Function:    Measure
// Purpose:		Reports the chirality of a tetrahedral center, given the 3D coordinates
// Input:       Reads coordinates from the atom objects contained in this object
// Output:      1 for counterclockwise, 2 for clockwise
// Requires:
/////////////////////////////////////////////////////////////////////////////

int Chiral::Measure()
{
    float ref_axis[3];
    float v1[3];
    float v2[3];

    ref_axis[0] = _subs[0]->x() - _center->x();
    ref_axis[1] = _subs[0]->y() - _center->y();
    ref_axis[2] = _subs[0]->z() - _center->z();

    v1[0] = _subs[1]->x() - _center->x();
    v1[1] = _subs[1]->y() - _center->y();
    v1[2] = _subs[1]->z() - _center->z();

    v2[0] = _subs[2]->x() - _center->x();
    v2[1] = _subs[2]->y() - _center->y();
    v2[2] = _subs[2]->z() - _center->z();

    ideal_volume = Cross_and_dot_3D(ref_axis, v1, v2) / 6.0;

    if (ideal_volume > 0)
    {
        return COUNTERCLOCKWISE;
    }
    else
    {
        return CLOCKWISE;
    }
}

/////////////////////////////////////////////////////////////////////////////
// Function:    ToCHIRAL
// Purpose:		Produce Chiral data in a format amenable to MIFit
// Input:       A vector of MIAtom objects which include the constituent atoms
//				this chiral object
//				A vector of (MIFit) MIAtom structs so that this ftn can locate
//				the new versions of the constituent atoms
// Output:      An MIFit CHIRAL struct
// Requires:	The general use case is to first convert a vector of MIAtom objects
//				to a vector of MIAtom structs for the molecule, then use those
//				two vectors as input to this function.
/////////////////////////////////////////////////////////////////////////////

CHIRAL Chiral::ToCHIRAL(const MIAtomList &new_atoms,
                        const std::vector<const MIAtom*> &old_atoms)
{
    CHIRAL chrl;
    std::string error_message;
    int i = GetIndex((const MIAtom*) _center, old_atoms);

    if (i < 0)
    {
        error_message = "Corrupted chirality information.";
        throw error_message;
    }
    chrl.center = new_atoms[i];

    i = GetIndex(const_cast<const MIAtom*> (_subs[0]), old_atoms);
    if (i < 0)
    {
        error_message = "Corrupted chirality information.";
        throw error_message;
    }
    chrl.setAtom1(new_atoms[i]);

    i = GetIndex(const_cast<const MIAtom*> (_subs[1]), old_atoms);
    if (i < 0)
    {
        error_message = "Corrupted chirality information.";
        throw error_message;
    }
    chrl.setAtom2(new_atoms[i]);

    i = GetIndex(const_cast<const MIAtom*> (_subs[2]), old_atoms);
    if (i < 0)
    {
        error_message = "Corrupted chirality information.";
        throw error_message;
    }
    chrl.atom3 = new_atoms[i];

    chrl.order = _order;

    return chrl;
}

/*
    void GraphPotentials(const Ligand &lig, std::vector<double> &gp) {
        double det;
        int n = lig.GetNumAtoms();
        gp.reserve(n);

        TNT::Matrix<double> g, c, graph_potentials;
        graph_potentials.newsize(n,1);

        construct_g_matrix(lig, g);
        TNT::InvertMatrix(g, det);
        construct_c_matrix(lig, c);
        graph_potentials = g*c;

        for(int i=0; i<n; ++i) {
            gp[i] = graph_potentials[i][0];
        }
    }
 */


std::string FindChiralCenters(const Monomer &res, const std::vector<Bond> &bonds)
{
    std::string list;

    std::vector<double>  gp;
    GraphPotentials(res, bonds, gp);

    std::vector<double> nabor_gplist;

    unsigned int i;
    MIAtom_const_iter ab = res.atoms().begin();
    MIAtom_const_iter ae = res.atoms().end();
    MIAtom_const_iter atm;
    MIAtom_const_iter nb;
    MIAtom_const_iter ne;
    MIAtom_const_iter nabor;
    std::vector<double>::iterator ngp_begin;
    std::vector<double>::iterator ngp_end;
    std::vector<double>::iterator nabor_gp;

    for (atm = ab; atm != ae; ++atm)
    {

        if ((*atm)->hybrid() != 3                     //Only consider sp3 atoms that
            || (*atm)->atomicnumber() == 7    //are not nitrogens, and have
            || HeavyDegree(*atm) < 3)               //at least 3 non-hydrogen neighbors
        {
            continue;
        }

        bool is_tetra_chiral = true;
        nb = (*atm)->nabors().begin();
        ne = (*atm)->nabors().end();
        nabor_gplist.clear();
        for (nabor = nb; nabor != ne; ++nabor)
        {
            ngp_begin = nabor_gplist.begin();
            ngp_end = nabor_gplist.end();

            //int nbr_index = *nabor - res.atoms().begin();
            bool found = false;
            int nbr_index = 0;
            //				ab = res.atoms().begin();
            //				ae = res.atoms().end();
            //				for(; ab != ae; ++ab) {
            for (i = 0; i < res.atoms().size(); ++i)
            {
                if (*nabor == res.atom(i))
                {
                    found = true;
                    nbr_index = i;
                    break;
                }
            }
            if (!found)
            {
                continue;
            }

            /* As such this section here should not be needed
               if (nbr_index < 0 ||
                nbr_index >= res.atoms().size()) {
                continue;
               }
             */
            for (nabor_gp = ngp_begin; nabor_gp != ngp_end; ++nabor_gp)
            {
                if (fabs(gp[nbr_index] - *nabor_gp) < 0.00001)
                {
                    is_tetra_chiral = false;
                }
            }

            if (is_tetra_chiral)
            {
                nabor_gplist.push_back(gp[nbr_index]);
            }
            else
            {
                break;
            }
        }

        if (is_tetra_chiral)
        {
            (*atm)->chiral_class(CH_TETRAHEDRAL);
            list += " ";
            list += (*atm)->name();
        }
        else
        {
            (*atm)->chiral_class(CH_NONE);
        }
    }
    return list;
}

std::string FindChiralCenters(const Residue *res, std::vector<Bond> &bonds, bool copyChiralClasses)
{
    Ligand mol(res, bonds);
    std::string result = FindChiralCenters(*(mol.residues.front()), mol.bonds);

    if (copyChiralClasses)
    {
        MIAtom *a1_old;
        MIAtom *a1_new;
        for (int i = 0; i < res->atomCount(); ++i)
        {
            a1_old = res->atom(i);
            a1_new = mol.residues.front()->atom(i);
            a1_old->chiral_class(a1_new->chiral_class());
        }
    }
    return result;
}

void GraphPotentials(const Monomer &res,
                     const std::vector<Bond> &bonds,
                     std::vector<double> &gp)
{

    double det;
    int n = res.atoms().size();
    gp.resize(n);

    TNT::Matrix<double> g, c, graph_potentials;
    graph_potentials.newsize(n, 1);

    construct_g_matrix(res, bonds, g);
    TNT::InvertMatrix(g, det);
    construct_c_matrix(res, c);
    graph_potentials = g*c;

    for (int i = 0; i < n; ++i)
    {
        gp[i] = graph_potentials[i][0];
    }
}

void GraphPotentials(const Residue &res,
                     const std::vector<Bond> &bonds,
                     std::vector<double> &gp)
{

    double det;
    gp.resize(res.atomCount() );

    TNT::Matrix<double> g, c, graph_potentials;
    graph_potentials.newsize(res.atomCount(), 1);

    construct_g_matrix(res, bonds, g);
    TNT::InvertMatrix(g, det);
    construct_c_matrix(res, bonds, c);
    graph_potentials = g*c;

    for (int i = 0; i < res.atomCount(); ++i)
    {
        gp[i] = graph_potentials[i][0];
    }
}

void construct_g_matrix(const Monomer &res,
                        const std::vector<Bond> &bonds,
                        TNT::Matrix<double> &m)
{
    m.newsize(res.atoms().size(), res.atoms().size());
    m = 0;

    //Encode atom information on the diagonal
    for (int i = 0; i < (int)res.atoms().size(); ++i)
    {
        //			for(j=0; j<res.atoms().size(); ++j) {
        //				if (i == j) {
        m[i][i] = Degree(res.atom(i)) + 1;
        m[i][i] += static_cast<double> (res.atom(i)->atomicnumber()) / 10.0;
        m[i][i] += static_cast<double> (res.atom(i)->hybrid()) / 100.0;
        //				}
        //				else {
        //					m[i][j] = 0;
        //				}
        //			}
    }

    //Encode a connection table in the off-diagonal elements, using 0
    //for non-bonded atom pairs and -1 for bonded atom pairs
    std::vector<Bond>::const_iterator bnd, bb, be;
    bb = bonds.begin();
    be = bonds.end();
    int xatom1, xatom2;
    for (bnd = bb; bnd != be; ++bnd)
    {
        xatom1 = res.indexOfAtom(bnd->getAtom1());
        xatom2 = res.indexOfAtom(bnd->getAtom2());

        if (xatom1 != -1 && xatom2 != -1)
        {
            m[xatom1][xatom2] = -1;
            m[xatom2][xatom1] = -1;
        }
    }
}

void construct_g_matrix(const Residue &res,
                        const std::vector<Bond> &bonds,
                        TNT::Matrix<double> &m)
{
    m.newsize(res.atomCount(), res.atomCount());

    for (int i = 0; i < res.atomCount(); ++i)
    {
        m[i][i] = Degree(res.atom(i), bonds);
        m[i][i] += static_cast<double> (res.atom(i)->atomicnumber()) /10.0;
        m[i][i] += static_cast<double> (res.atom(i)->hybrid()) / 100.0;
    }

    std::vector<Bond>::const_iterator bnd, be = bonds.end();
    int xatom1, xatom2;
    for (bnd = bonds.begin(); bnd != be; ++bnd)
    {
        xatom1 = GetAtomIndex(bnd->getAtom1(), res);
        xatom2 = GetAtomIndex(bnd->getAtom2(), res);

        if (xatom1 != -1 && xatom2 != -1)
        {
            m[xatom1][xatom2] = -1;
            m[xatom2][xatom1] = -1;
        }
    }
}

void construct_c_matrix(const Ligand &lig, TNT::Matrix<double> &m)
{
    m.newsize(lig.GetNumAtoms(), 1);

    int n = 0;
    for (unsigned int i = 0; i < lig.residues.size(); ++i)
    {
        for (unsigned int j = 0; j < lig.residues[i]->atoms().size(); ++j)
        {
            m[n++][0] = Degree(lig.residues[i]->atom(j));
        }
    }
}

void construct_c_matrix(const Monomer &res, TNT::Matrix<double> &m)
{
    m.newsize(res.atoms().size(), 1);

    for (int i = 0; i < (int)res.atoms().size(); ++i)
    {
        m[i][0] = Degree(res.atom(i));
    }
}

void construct_c_matrix(const Residue &res,
                        const std::vector<Bond> &bonds,
                        TNT::Matrix<double> &m)
{
    m.newsize(res.atomCount(), 1);
    for (int i = 0; i < res.atomCount(); ++i)
    {
        m[i][0] = Degree(res.atom(i), bonds);
    }
}

} //namespace chemlib
