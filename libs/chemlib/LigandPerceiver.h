#ifndef LIGAND_PERCEIVER_H
#define LIGAND_PERCEIVER_H

#include "MIAtom_fwd.h"
#include "Residue.h"
#include "Ligand.h"

namespace chemlib
{
    class LigandPerceiver
    {
    public:
        void AssignImpHydrogens(Ligand *lig);
        void AssignImpHydrogens(Monomer &res, const std::vector<Bond> bonds);
        void AssignImpHydrogens(MIAtom &atom, const std::vector<Bond> bonds);

        void AssignChirality(Ligand*);
        void AssignChirality(Monomer&, Ligand*);
        int DefaultChiralClass(MIAtom&);

        void AssignHybridization(Ligand*);
        void AssignHybridization(Monomer&, Ligand*);
        void AssignHybridization(MIAtom*, Ligand*);
        void AdjustHybridization(MIAtom*, Ligand*);
        void AssignAtomGeom(Ligand*);
        void AssignAtomGeom(Monomer&);
        void AssignAtomGeom(MIAtom&);

    };

    bool Is_Onium(MIAtom &atom, std::vector<Bond> bonds);

//	void PrepPolarAtom(MIAtom &atom);
//	void PrepPolarAtoms(Residue &res);

}   //namespace chemlib
#endif //LIGAND_PERCEIVER_H
