#ifndef COV_GEOM_H
#define COV_GEOM_H

#include "MIAtom_fwd.h"
#include "Bond.h"
#include "Residue.h"
#include "Constraint.h"
#include "Ligand.h"
#include "mol_util.h"

#include <vector>

namespace chemlib
{

    class CovalentGeometry
    {
    public:
        CovalentGeometry(Ligand*, const Residue*);
        void AssignResidue();
        void InitBondAngle(MIAtom*, MIAtom*, MIAtom*, Angle&);
        //	void AssignBondAngle(MIAtom *, MIAtom *, MIAtom *);
        //	void AssignBondAngle(MIAtom *, MIAtom *, MIAtom *, float value);
        float AssignBondLength(Bond*);
        void AssignBump(MIAtom*, MIAtom*);
        void AngleTolerance(Angle&);
        void LengthTolerance(BondLength&);
        void BumpTolerance(Bump&);
        //		void AngleValue(Angle &, BondLength &, BondLength &);
        void ConvertToDistance(Angle&, float, float);
        //	void OneThreeDistance(Angle &, BondLength &, BondLength &, float value);
        float AngleFromGeom(int);
        void SetRingData(Angle&);
        void CalcStrainAngle(Angle*);
        void LengthValue(BondLength&);

        void AngleValue(Angle&);
        void BumpValue(Bump&);
        //	void Print();
    private:
        void AddLength(const char*, const char*, int, int, unsigned char, float, BondLength&);
        void AddAngle(const char*, const char*, const char*, int, int, int,
                      unsigned char, unsigned char, float, Angle&, const Bond&, const Bond&);
        const Residue *_res;
        Ligand *_lig;
        ConstraintList _cl;
        float _sigmaangle;
        float _sigmabond;
        float _sigmabump;
        Bond *_bond;                                        //Used when assigning a bond length
    };


    float TrigAngleRemainder(std::vector <Angle> &fixed);
    float TetraAngleRemainder(std::vector <Angle> &fixed);
    float CheckAngle(MIAtom *atom1, MIAtom *atom2, std::vector<Angle> &angles);

    float IdealBondLength(const Bond &bond);
    float IdealBondLength(const Bond &bond);
    float IdealBondLength(int element1, int element2, unsigned char bnd_order);

    unsigned char PredictBondOrder(int element1, int element2, float distance);

} //namespace chemlib

#endif //COV_GEOM_H
