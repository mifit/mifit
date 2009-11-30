#ifndef MI_MOL_UTIL_PRIVATE
#define MI_MOL_UTIL_PRIVATE

//@{
// Given two atoms and the two residues containing them, return true if they are h-bonded.
// Used in the case of two protein or nucleic atoms - not neccesisarily valid for a ligand.
//@}
namespace chemlib
{
    class MIAtom;
    class RESIDUE;
    bool hbondable(const MIAtom &a1, const MIAtom &a2, const Residue &r1, const Residue &r2);


    void clear_residue_from_atom_cache(const Residue *res);

}
#endif // ifndef MI_MOL_UTIL_PRIVATE
