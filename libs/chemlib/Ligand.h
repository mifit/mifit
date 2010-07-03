#ifndef LIGAND_H
#define LIGAND_H

#include <vector>

#include "MIAtom_fwd.h"
#include "Bond.h"
#include "Residue.h"
#include "RingSystem.h"
#include "Constraint.h"
#include "Substituent.h"
#include "util.h"
#include "ANGLE.h"

namespace chemlib
{

    class Ligand
    {

    public:
        Ligand(int type = Ligand::Other);
        Ligand(const Residue &res, const std::vector<Bond> &bonds);
        Ligand(const Residue *res, const std::vector<Bond> &bonds);
        Ligand(const std::vector<Monomer*> &orig_res,
               const std::vector<Bond> &orig_bonds);
        ~Ligand();

        std::vector<Bond> bonds; // array of bonds
        std::vector<Bond> connects;
        std::vector<Monomer*> residues; //Residue sequence
        std::vector<RingSystem> ringsystems;
        ConstraintList geometry;
        //		Residue * SymmResidues;
        int xmin, ymin, zmin, xmax, ymax, zmax;
        char link_here[MAXNAME], link_next[MAXNAME]; //MAXNAME is in define.h

        // Export internal data to these vectors
        void Export(std::vector<Residue*> &rdues, std::vector<Bond> &bnds);

        void GetConstraints(const MIAtomList &mifit_atms,
                            const std::vector<Residue*> &mifit_residues,
                            std::vector<Bond> &bond_lengths,
                            std::vector<ANGLE> &angles,
                            std::vector<TORSION> &torsions,
                            std::vector<TORSION> &impropers,
                            std::vector<PLANE> &planes,
                            std::vector<CHIRAL> &chirals);

        Monomer *AddRes(const std::string&, const std::string&, unsigned short linkage_type = 0,
                        unsigned short chain_id = ' ', char = 'U');
        Monomer *AddRes(const Residue&, const std::vector<Bond> &orig_bonds);
        int GetNumAtoms() const;

        void AddBond(Bond &bond);
        void AddBond(MIAtom*, MIAtom*);
        void AddBond(MIAtom*, MIAtom*, unsigned char, char);
        void ClearBondOrders();
        void GuessBondOrders();             //Uses 3D-coords for hybridization and orders
        void Connect(const Bond &connect);
        void FixAtomicNumbers();
        //		void FindCycles();
        //		void AddCycle(int);
        int FindRingSystems();
        void InitRingData();
        void SetRingFlags();
        void ResetSearchFlags();
        void FreeBonds();
        MIAtom *MoreCentralAtm(MIAtom*, MIAtom*);
        void DepthFirstSearch(MIAtom*, MIAtom*, std::vector <MIAtom*>&);
        MIAtom *GetNabor(const MIAtom*, int);
        MIAtom *GetNewNabor(MIAtom*, MIAtom*);
        Bond *GetBond(int);
        const Bond *GetBond(int) const;
        Bond *GetBond(const MIAtom*, const MIAtom*);
        bool AlreadyBonded(const MIAtom*, const MIAtom*) const;
        void GetResidueBonds(const Monomer*, std::vector <Bond*>&);
        void GetResidueIntBonds(const Monomer*, std::vector <Bond>&);
        enum
        {
            PDB,
            mmCIF,
            XML,
            New,
            SHELX,
            Pentamer,
            Smiles,
            Other
        };
        int ModelType;
        int GetIndex(MIAtom*, MIAtom**, int);
        void TranslateFragment(MIAtom*, MIAtom*, double*);      //Translate all atoms on one side of a bond
        void FlipBond(Bond&);                       //Rotate 180 deg the coordinates on one side of a bond
        void Flatten();
        void CleanUp2D();
        void FlattenRings();
        void Unroll();
        void UnrollAtom(MIAtom*, MIAtom*);
        void GroupSubstituents(MIAtom *origin, MIAtom *omit, std::vector<Substituent> &subs) const;
        void Translate(double xstep, double ystep, double zstep);
        void Scale(double scale_factor);
        void LSqrPlane(double[], double*) const;        //Calc a best fit plane through all the atoms
        void RandomizeCoords(std::pair<double, double> xbounds,
                             std::pair<double, double> ybounds,
                             std::pair<double, double> zbounds);

        void MapAtomPtrs(std::map<MIAtom*, int> &atom_map);
        void MapAtomPtrs(std::map<const MIAtom*, int> &atom_map) const;
    };

} //namespace chemlib

#endif //LIGAND_H
