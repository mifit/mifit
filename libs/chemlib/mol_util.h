#ifndef MOL_UTIL_H
#define MOL_UTIL_H

#include "MIAtom.h"
#include "Bond.h"
#include "Residue.h"
#include "Ligand.h"
#include "RingSystem.h"
#include "math_util.h"
#include "Constraint.h"
#include "PLANE.h"

#include <set>

namespace chemlib
{

    inline void BondVector(const MIAtom *source, const MIAtom *target, double v[])
    {
        v[0] = target->x() - source->x();
        v[1] = target->y() - source->y();
        v[2] = target->z() - source->z();
    }

    inline void BondVector(const MIAtom *source, const MIAtom *target, std::vector<double> &v)
    {
        v[0] = target->x() - source->x();
        v[1] = target->y() - source->y();
        v[2] = target->z() - source->z();
    }

    const char *Atomic_Name(int atomicnumber);

    const char *Left_Atomic_Name(int atomic_number);

    int Atomic_Number(const char *symbol);
    int Atomic_Number_Nformat(const std::string &atname);

    float CovalentRadius(int atomicnumber);

    int Electroneg(int atomic_number);

    int MaxValence(int atomicnumber);

    double BondLimit(const char *s);

//@{
// Copy the coordinates from one residue to another, using atom names
// to match corresponding atoms.  OK if atoms are ordered differently
//@}
    bool MICopyCoords(RESIDUE *target, const RESIDUE *source);

    inline MIAtomList::size_type
    NumAtomsSum(MIAtomList::size_type sumSoFar, const Residue *res)
    {
        return sumSoFar + res->atoms().size();
    }

    inline int NumAromaticsSum(int sumSoFar, const RingSystem &rs)
    {
        return sumSoFar + rs.NumAromatics();
    }

    inline std::string AppendBond(std::string stringSoFar, const Bond &bond)
    {
        return stringSoFar + bond.Print();
    }

    inline std::string AppendRingSystem(std::string stringSoFar, const RingSystem &rs)
    {
        return stringSoFar + rs.Print();
    }

    inline std::string AppendAromatics(std::string stringSoFar, const RingSystem &rs)
    {
        return stringSoFar + rs.PrintAromatics();
    }

    inline double SignedVolume(MIAtom *center, MIAtom *a1, MIAtom *a2, MIAtom *a3)
    {
        double u[3], v[3], w[3];
        BondVector(center, a1, u);
        BondVector(center, a2, v);
        BondVector(center, a3, w);
        return Cross_and_dot_3D(u, v, w);
    }

    struct IsBondRotatable
        : std::unary_function<Bond, bool>
    {
        bool operator()(const Bond &bond) const;
    };

//@{
// Returns the distance (in Angstroms) between two atoms given their pointers.
//@}
    inline double AtomDist(const MIAtom &a, const MIAtom &b)
    {
        return sqrt((a.x() - b.x()) * (a.x() - b.x())
                    +(a.y() - b.y()) * (a.y() - b.y())
                    +(a.z() - b.z()) * (a.z() - b.z()));
    }

//@{
// Returns the angle (in angle) between three atoms given their pointers.
// The second atom is the vertex of the angle.
//@}
    double CalcAtomAngle(const chemlib::MIAtom&, const chemlib::MIAtom&, const chemlib::MIAtom&);

    inline double SquaredAtomDist(const MIAtom &a, const MIAtom &b)
    {
        return (a.x() - b.x()) * (a.x() - b.x())
               +(a.y() - b.y()) * (a.y() - b.y())
               +(a.z() - b.z()) * (a.z() - b.z());
    }

    inline void AtomStep(MIAtom *atom, const float *step_direction, float step_size)
    {
        atom->translate(step_size * step_direction[0],
                        step_size * step_direction[1],
                        step_size * step_direction[2]);
    }

    inline void AtomStep(MIAtom *atom, const double *step_direction, double step_size)
    {
        atom->translate((float)(step_size * step_direction[0]),
                        (float)(step_size * step_direction[1]),
                        (float)(step_size * step_direction[2]));
    }

    inline void MoveIntoPlane(MIAtom *atom,
                              const double p_normal[],
                              const double p_displace)
    {
        double dev;
        dev = -DotVect(atom->x(), atom->y(), atom->z(),
                       p_normal[0], p_normal[1], p_normal[2]);
        dev += p_displace;

        AtomStep(atom, p_normal, dev);
    }

    void GatherAtmPtrs(MIAtomList &target, const std::vector<Residue*> &source);

    void EnumerateTorsions(const Bond *bond, std::vector< MIAtomList > &torsions);

//	int GetIndex(const MIAtom *query, const std::vector<MIAtom *> &domain);


    double CalcAtomTorsion(MIAtom *na, MIAtom *nb, MIAtom *nc, MIAtom *nd);
    double CalcAtomTorsion(const MIAtom *na, const MIAtom *nb, const MIAtom *nc, const MIAtom *nd);

    void lsqplane(PLANE *plane);
    void lsqplane(Plane &plane, float *normal, float *displace);

    void DepthFirstSearch(MIAtom *root,
                          const MIAtom *prev,
                          const MIAtom *block,
                          std::vector <MIAtom*> &aggregate);

    void DepthFirstSearch(MIAtom *root,
                          const MIAtom *prev,
                          const MIAtom *block,
                          const std::vector<Bond> &bonds,
                          MIAtomList &aggregate);

    void ExclDepthFirstSearch(MIAtom *root,
                              const MIAtom *block,
                              const std::vector<Bond> &bonds,
                              MIAtomList &aggregate);


    void TrimBonds(std::vector<Bond> &trimmed_bonds,
                   const std::vector<Bond> &orig_bonds,
                   const MIAtomList &atoms);

    int CountResBonds(const RESIDUE &res, const std::vector<Bond> &bonds);

    int Build(Ligand &lig);

    bool IsPeptide(const Residue &res);
    bool IsNucleic(Residue *res);
    int IsDna(const Residue &res);
    bool IsNucleic(RESIDUE *res);
    int IsDna(const RESIDUE &res);
    bool IsWater(const RESIDUE *res);
    bool IsPolarH(const MIAtom *atom, const RESIDUE *res);

//@{
//
//@}
    void AminoOrNucleic(RESIDUE *res1, int reset);
//@{
// Return true if the residue is an amino acid.
//@}
    bool IsAminoAcid(RESIDUE *res);



//@{
// Given a residue and an atom name, return a pointer to the atom that matches, if any.
// Returns null if no match is found.
//@}
    MIAtom *atom_from_name(const char *name, const Residue &residue);

//@{
// Returns the residue that matches the name and chainid given a residue list.
// Always returns the first one found that matches the name so that if there
// were a later residue in the list with the same name, it will e hidden.
// @param res a linked list of residues.
// @param name a string containing the name such as "104".
// @param chain the chain id of the residue such as "A".  Use "*" to match a space.
//@}
    RESIDUE *residue_from_name(RESIDUE *res, const char *name, const char chain);

//@{
// Given a linked list of residues and an atom, returns the residue containing that atom.
// If the atom is not found, the function returns NULL.
//@}
    RESIDUE *residue_from_atom(RESIDUE*, MIAtom*);

    void BisectAtom(const MIAtom *source, double *v);

    std::vector<Residue*>::iterator
    ResSearch(const MIAtom *query,
              std::vector<Residue*>::iterator res_begin,
              std::vector<Residue*>::iterator res_end);

    int MaxNumBonds(const MIAtom *atom);
    int CurrentValence(const MIAtom &atom,
                       const std::vector<Bond> &bonds);

    int UnusedValences(const MIAtom &atom,
                       const std::vector<Bond> &bonds);


    int CountAromaticBonds(const MIAtom &atom, const std::vector<Bond> &bonds);

    bool AtSixSixFusion(const MIAtom &atom, const std::vector<Bond> &bonds);
    bool AtFiveSixFusion(const MIAtom &atom, const std::vector<Bond> &bonds);
    bool AtFiveFiveFusion(const MIAtom &atom, const std::vector<Bond> &bonds);

    void GuessBondOrders(RESIDUE *res, std::vector<Bond> &bonds);

    void HybridizeTerminalAtom(MIAtom &atom, std::vector<Bond> &bonds);
    float AverageBondAngle(MIAtom &atom);
    void HybridFromGeom(MIAtom &atom);
    int PredictValenceFrom3D(MIAtom &atom);

//@{
// Returns a 1-letter name given the 3-letter name for a residue.
//@}
    char singleletter(const char*);
//@{
// Returns a 3-letter name given the 1-letter name for a residue.
//@}
    const char *tripleletter(const char t);

//@{
// From a list of residues builds a vector of links (bonds/edges).
//@}
    int BuildCALinks(std::vector<Bond> &bonds, const RESIDUE *res);

    const char *Atomic_Name(int atomic_number);

} //namespace chemlib

namespace MolFrom3D
{
    const float PLANE_TOLERANCE = 0.05F;
}
#endif //MOL_UTIL_H
