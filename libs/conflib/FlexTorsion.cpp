#include <math/mathlib.h>
#include <chemlib/chemlib.h>
#include <chemlib/RESIDUE_.h>

#include "FlexTorsion.h"

using namespace chemlib;
using namespace conflib;
using namespace std;

//Boost graph library
#include <boost/graph/graph_traits.hpp>
#include <boost/graph/adjacency_list.hpp>
#include <boost/graph/breadth_first_search.hpp>
#include <boost/graph/visitors.hpp>
#include <boost/property_map.hpp>
using namespace boost;




//Typedefs for boost graphs
//struct atom_address_t {
//	typedef vertex_property_tag kind;
//};
//typedef property<atom_address_t, MIAtom*> AtomAddressProperty;
typedef adjacency_list<vecS, vecS, undirectedS> AtomGraph;
typedef graph_traits<AtomGraph>::vertices_size_type AtomGraphIndex; //effectively an int
//typedef typename property_map<AtomGraph, atom_address_t>::type AtomAddressMap;


FlexTorsion::FlexTorsion(const TORSION &t, const vector<Bond> &bonds, bool set)
{
    //Look up the central bond of this torsion from the sequence of bonds, and
    //store the order for use in determining ideal angle values
    _bond = GetBond(t.getAtom2(), t.atom3, bonds);

    _a2 = SetUpTorsion(t.res, bonds, t.getAtom2(), t.atom3, _flex_atoms);

    if (_a2 == t.atom3)
    {
        _a1 = t.atom4;
        //		_a2 = t.atom3;
        _a3 = t.getAtom2();
        _a4 = t.getAtom1();
    }
    else
    {
        _a1 = t.getAtom1();
        _a2 = t.getAtom2();
        _a3 = t.atom3;
        _a4 = t.atom4;
    }

    if (t.nideal == 0)
    {
        GuessIdeal();
    }
    else
    {
        SetIdeal(t.nideal, t.ideal);
    }

    _dihedral_index = _dihedrals.begin();

    if (set)
    {
        Set();
    }
}

FlexTorsion::FlexTorsion(const FlexTorsion &rhs)
    : _a1(rhs._a1),
      _a2(rhs._a2),
      _a3(rhs._a3),
      _a4(rhs._a4),
      _dihedrals(rhs._dihedrals),
      _flex_atoms(rhs._flex_atoms),
      _bond(_bond)
{
    _dihedral_index = _dihedrals.begin();
}

FlexTorsion&FlexTorsion::operator=(const FlexTorsion &rhs)
{
    if (this == &rhs)                       //Check for self-assignment by comparing addresses
    {
        return *this;
    }

    _a1 = rhs._a1;
    _a2 = rhs._a2;
    _a3 = rhs._a3;
    _a4 = rhs._a4;
    _dihedrals = rhs._dihedrals;
    _dihedral_index = _dihedrals.begin();   //Don't copy this, point to local container
    _flex_atoms = rhs._flex_atoms;          //New copy points to the same atoms
    _bond = _bond;

    //Return a reference to the newly assigned object
    return *this;
}

void FlexTorsion::Set()
{
    float rot_angle = *_dihedral_index - (float)CalcAtomTorsion(_a1, _a2, _a3, _a4);

    float mat[4][3];
    initrotate(_a2->x(), _a2->y(), _a2->z(), _a3->x()-_a2->x(), _a3->y()-_a2->y(), _a3->z()-_a2->z(),
               rot_angle,
               mat);

    vector<MIAtom*>::iterator i = _flex_atoms.begin();
    vector<MIAtom*>::iterator e = _flex_atoms.end();
    float x, y, z;
    while (i != e)
    {
        x = (*i)->x();
        y = (*i)->y();
        z = (*i)->z();
        rotate(&x, &y, &z, mat);
        (*i)->setPosition(x, y, z);
        (*i)->set_search_flag(1);
        i++;
    }
}

float FlexTorsion::Measure()
{
    return (float)CalcAtomTorsion(_a1, _a2, _a3, _a4);
}

void FlexTorsion::Pick(int index)
{
    if (index < 0 || index >= (int)_dihedrals.size())
    {
        return;
    }

    _dihedral_index = _dihedrals.begin() + index;

    Set();
}

void FlexTorsion::SetIdeal(int n, const float *angles)
{
    _dihedrals.clear();

    for (int i = 0; i < n; ++i)
    {
        _dihedrals.push_back(angles[i]);
    }

    _dihedral_index = _dihedrals.begin();
}

void FlexTorsion::GuessIdeal()
{
    _dihedrals.clear();

    /* Can't quite use this strategy yet, largely because not all partial double bonds
       (amides, esters, anilines, etc.) are labeled as such, so they would end up treated
       like single bonds
       //Handle higher-order bonds here
       if (_bond->getOrder() == DOUBLEBOND) {
        _dihedrals.push_back(0.0F);
        _dihedrals.push_back(180.0F);
       }
       else if (_bond->getOrder() == PARTIALDOUBLEBOND) {
        _dihedrals.push_back(0.0F);
        _dihedrals.push_back(180.0F);
       }
       else if (_bond->getOrder() == TRIPLEBOND) {			//Triple bonds shoudn't ever be flexible,
        _dihedrals.push_back(180.0F);				//but just for completeness
       }
     */

    if (_a2->iscyclic() && _a3->iscyclic()              //This is the special case of biphenyl linkage
        && _a2->hybrid() == 2 && _a3->hybrid() == 2
        && _bond->iscyclic == 0
        && _bond->getOrder() == SINGLEBOND)
    {
        _dihedrals.push_back(-150.0F);
        _dihedrals.push_back(-90.0F);
        _dihedrals.push_back(-30.0F);
        _dihedrals.push_back(30.0F);
        _dihedrals.push_back(90.0F);
        _dihedrals.push_back(150.0F);
    }
    else if (_a2->hybrid() == 1
             && _a3->hybrid() == 1)
    {
        _dihedrals.push_back(0.0F);
    }
    else if (_a2->hybrid() == 1
             && _a3->hybrid() == 2)
    {
        _dihedrals.push_back(0.0F);
    }
    else if (_a2->hybrid() == 1
             && _a3->hybrid() == 3)
    {
        _dihedrals.push_back(0.0F);
    }
    else if (_a2->hybrid() == 2
             && _a3->hybrid() == 1)
    {
        _dihedrals.push_back(0.0F);
    }
    else if (_a2->hybrid() == 2
             && _a3->hybrid() == 2)
    {
        _dihedrals.push_back(0.0F);
        _dihedrals.push_back(90.0F);   //Added by JB to generate a correct conformation for ASC member ligand
        _dihedrals.push_back(180.0F);
    }
    else if (_a2->hybrid() == 2
             && _a3->hybrid() == 3)
    {
        _dihedrals.push_back(-150.0F);
        _dihedrals.push_back(-90.0F);
        _dihedrals.push_back(-30.0F);
        _dihedrals.push_back(30.0F);
        _dihedrals.push_back(90.0F);
        _dihedrals.push_back(150.0F);
    }
    else if (_a2->hybrid() == 3
             && _a3->hybrid() == 1)
    {
        _dihedrals.push_back(0.0F);
    }
    else if (_a2->hybrid() == 3
             && _a3->hybrid() == 2)
    {
        _dihedrals.push_back(-150.0F);
        _dihedrals.push_back(-90.0F);
        _dihedrals.push_back(-30.0F);
        _dihedrals.push_back(30.0F);
        _dihedrals.push_back(90.0F);
        _dihedrals.push_back(150.0F);
    }
    else if (_a2->hybrid() == 3
             && _a3->hybrid() == 3)
    {
        _dihedrals.push_back(-60.0F);
        _dihedrals.push_back(60.0F);
        _dihedrals.push_back(180.0F);
    }
    else
    {
        _dihedrals.push_back(-150.0F);
        _dihedrals.push_back(-90.0F);
        _dihedrals.push_back(-30.0F);
        _dihedrals.push_back(30.0F);
        _dihedrals.push_back(90.0F);
        _dihedrals.push_back(150.0F);
    }

    _dihedral_index = _dihedrals.begin();
}

bool FlexTorsion::Advance()
{

    bool rolled_over = false;
    _dihedral_index++;
    if (_dihedral_index == _dihedrals.end())
    {
        _dihedral_index = _dihedrals.begin();
        rolled_over = true;
    }


    Set();          //Update coordinates here


    return !rolled_over;        //Return true if advance went normally, false if we reset
                                //back to the start of the list of ideal angles
}

namespace conflib
{
int SumDistances(const AtomGraph &graph, AtomGraphIndex *distances, AtomGraphIndex source)
{

    //Initialize all the distances to zero
    std::fill_n(distances, num_vertices(graph), 0);

    //Search the graph, recording distances in the distances array
    breadth_first_search(graph,
                         vertex(source, graph),
                         visitor(make_bfs_visitor(record_distances(distances, on_tree_edge()))) );

    //Count up the distances
    return std::accumulate(distances, distances + num_vertices(graph), 0);
}

void GatherFlexAtoms(const AtomGraphIndex *distances, const RESIDUE *res, vector<MIAtom*> &atoms)
{
    int i;
    for (i = 0; i < res->atomCount(); ++i)
    {
        if (distances[i] > 0)
        {
            atoms.push_back(res->atom(i));
        }
    }
}

MIAtom *SetUpTorsion(const RESIDUE *res,
                     const std::vector<Bond> &bonds,
                     MIAtom *atom1,
                     MIAtom *atom2,
                     vector<MIAtom*> &flex_atoms)
{
    AtomGraph ag(res->atomCount());

    int xAtom1, xAtom2;
    for (unsigned int i = 0; i < bonds.size(); ++i)
    {
        if ((xAtom1 = GetAtomIndex(bonds[i].getAtom1(), *res)) < 0)
        {
            continue;
        }
        if ((xAtom2 = GetAtomIndex(bonds[i].getAtom2(), *res)) < 0)
        {
            continue;
        }


        //Don't put in the central bond of the torsion...removing
        //this bond splits the graph into two parts, which allows
        //us to search the subgraph on each side of the torsion
        if ((bonds[i].getAtom1() == atom1 && bonds[i].getAtom2() == atom2)
            || (bonds[i].getAtom1() == atom2 && bonds[i].getAtom2() == atom1))
        {
            continue;
        }

        add_edge(xAtom1, xAtom2, ag);
    }

    AtomGraphIndex *d2 = new AtomGraphIndex[res->atomCount()];
    AtomGraphIndex *d3 = new AtomGraphIndex[res->atomCount()];

    MIAtom *root;
    if (SumDistances(ag, d2, GetAtomIndex(atom1, *res)) <
        SumDistances(ag, d3, GetAtomIndex(atom2, *res)))
    {
        root = atom2;
        GatherFlexAtoms(d2, res, flex_atoms);
    }
    else
    {
        root = atom1;
        GatherFlexAtoms(d3, res, flex_atoms);
    }

    delete[] d2;
    delete[] d3;

    return root;
}

} //namespace conflib



