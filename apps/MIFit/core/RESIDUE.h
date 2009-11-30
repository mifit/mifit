#ifndef MIFIT_MODEL_RESIDUE_H_
#define MIFIT_MODEL_RESIDUE_H_

#include <string>
#include <chemlib/chemlib.h>
#include "corelib.h"


class CMapHeaderBase;

// residue.linkage_type definitions
#define NTERMINUS ((unsigned short)1)
#define FIRST ((unsigned short)1)
#define MIDDLE ((unsigned short)2)
#define LAST ((unsigned short)4)
#define CTERMINUS ((unsigned short)4)
#define BRANCH ((unsigned short)8)
#define LINKAGEMASK ((unsigned short)15)
#define PEPTIDE ((unsigned short)16)
#define NUCLEIC ((unsigned short)32)


//@{
// Synthesizes a string representing the residue from the name, type, and chainid.
// Used for printing to the screen.
//@}
std::string resid(const chemlib::RESIDUE*);
//@{
// convert one letter residue names to single letter charge type
//@}
char chargetype(char);


//@{
// Returns the phi angle between two residues.
//@}
float phi(const chemlib::RESIDUE*, const chemlib::RESIDUE*);
//@{
// Returns the psi angle between two residues.
//@}
float psi(const chemlib::RESIDUE*, const chemlib::RESIDUE*);

//@{
// Read the colors into a residue from a buffer for loading.
// This function is defunct with the new XML format.
// It is needed for backwards compatibility.
//@}
int read_colors(const chemlib::RESIDUE *res, char *buf, int nbuf);
//@{
// Read the radii into a residue from a buffer for loading.
// This function is defunct with the new XML format.
// It is needed for backwards compatibility.
//@}
int read_radii(const chemlib::RESIDUE *res, char *buf, int nbuf);

//@{
// Makes a string for listing a chain in the tree control
//@}
const std::string chainstring(const chemlib::RESIDUE *res);

//@{
// Move the residue(s) in source onto target.
// @param source the residue(s) to be moved.
// @param target the residue(s) to be targeted (stays still).
// @param nres the number of residues to move starting at source and target in the linked lists.
//@}
bool MoveOnto(const chemlib::RESIDUE *source, chemlib::RESIDUE *target, int nres = 1);

//@{
// Calculates the symmetry residues around a center point.
//@}
chemlib::RESIDUE*SymmResidue(const chemlib::RESIDUE *Model, CMapHeaderBase * mh, float center[3], float r = 10.0F, int color = Colors::MAGENTA);

//@{
// Function for the backbone builder.
//@}
chemlib::RESIDUE *matchvec(const std::vector<chemlib::MIAtom*> &CA, const std::vector<chemlib::MIAtom*> &CB, std::string &pentdir, FILE*);
//@{
// Function for the backbone builder.
//@}
chemlib::RESIDUE *pdbvec(std::string &pentdir, const std::vector<chemlib::MIAtom*> &CA, const std::vector<chemlib::MIAtom*> &CB);

//@{
// Function for the backbone builder.
//@}
int atomvector(const chemlib::MIAtom *atom1, const chemlib::MIAtom *atom2);

//@{
// build a residue from a list of atoms.
// the other residue still need to be filled in.
//@}
chemlib::RESIDUE *make_res(const std::vector<chemlib::MIAtom*> &atoms);
//@{
// Build a chain.
//@}
void getchain(unsigned short chainid, chemlib::RESIDUE *reslist, chemlib::RESIDUE* &nter, chemlib::RESIDUE* &cter);
//@{
// Figures out the correct order for two residues providing both are in the list.
// Very useful for ordering two picks by the user into a useful range.
//@}
int order_ends(chemlib::RESIDUE* &nter, chemlib::RESIDUE* &cter, chemlib::RESIDUE *reslist);

//@{
// Returns the color for atoms being fit or refined.
// If the atom is being both fit and refined, the color will be for fitting.
//@}
inline int getColor(const chemlib::MIAtom *a)
{
    if (a->type() & chemlib::AtomType::FITATOM)
    {
        return Colors::GREEN;
    }
    if (a->type() & chemlib::AtomType::REFIATOM)
    {
        return Colors::CYAN;
    }
    return (int)a->color();
}

#endif /*MIFIT_MODEL_RESIDUE_H_*/

