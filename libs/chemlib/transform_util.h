#ifndef TRANSFORM_UTIL_H
#define TRANSFORM_UTIL_H

#include "MIAtom_fwd.h"
#include "Ligand.h"

namespace chemlib
{

    inline double DistToOrigin(MIAtom *atom)
    {
        return sqrt(atom->x() * atom->x()
                    +atom->y() * atom->y()
                    +atom->z() * atom->z());
    }

    inline void DoReflectY(MIAtom *atom)
    {
        atom->setPosition(atom->x(), -atom->y(), atom->z());
    }

    struct DoTransformAtom : public std::binary_function<MIAtom*, const double*, bool>
    {
        bool operator()(MIAtom *atom, const double *mat) const;
    };


    inline bool DoTransformAtom::operator ()(MIAtom *atom, const double *mat) const
    {
        float x = atom->x();
        float y = atom->y();
        float z = atom->z();

        atom->setPosition((float)(mat[0] * x
                                  +mat[1] * y
                                  +mat[2] * z),
                          (float)(mat[3] * x
                                  +mat[4] * y
                                  +mat[5] * z),
                          (float)(mat[6] * x
                                  +mat[7] * y
                                  +mat[8] * z));
        return true;
    }

    struct DoTranslateAtom : public std::binary_function<MIAtom*, const double*, bool>
    {
        bool operator()(MIAtom *atom, const double *v) const;
    };

    inline bool DoTranslateAtom::operator ()(MIAtom *atom, const double *v) const
    {

        atom->translate((float)v[0], (float)v[1], (float)v[2]);

        return true;
    }

    inline void ReflectY(MIAtomList &atoms)
    {
        std::for_each(atoms.begin(), atoms.end(), DoReflectY);
    }

    inline void TransformAtoms3D(const double *mat, MIAtomList &atoms)
    {
        if (mat == 0)
        {
            return;
        }

        std::for_each(atoms.begin(), atoms.end(), std::bind2nd(DoTransformAtom(), mat));
    }

    inline void TranslateAtoms3D(const double *v, MIAtomList &atoms)
    {
        if (v == 0)
        {
            return;
        }

        std::for_each(atoms.begin(), atoms.end(), std::bind2nd(DoTranslateAtom(), v));
    }

    void RotateIntoXYPlane(MIAtom *ref, MIAtomList &branch);
    void RotateIntoXYPlane(MIAtomList &ref, MIAtomList &branch);
    void RotateIntoXYPlane(std::vector<double> &normal, MIAtomList &branch);
    void RotateIntoXYPlane(double *normal, MIAtom **atoms, int natoms);

    void FitToXYPlane(Ligand *mol);

    void TranslateToOrigin(MIAtom *ref, Ligand *mol);
    void RotateToXAxis(MIAtom *ref, Ligand *mol);

    double AngleToXYPlane(MIAtom *atom);
    void XAxisRotate(double theta, MIAtomList &domain);
    void CalcNormal(MIAtomList &ref, std::vector<double> &normal);

    void CalcRotationToZAxis(std::vector<double> &normal,
                             std::vector<double> &axis,
                             std::vector<double> &cos_sin);

    void CalcRotMatrix(std::vector<double> rot_axis,
                       double cos_alpha,
                       double sin_alpha,
                       double mat[3][3]);

    void CalcRotMatrix(double vx,
                       double vy,
                       double vz,
                       double cos_alpha,
                       double sin_alpha,
                       double mat[3][3]);

    void CalcRotMatrix(double *rot_axis,
                       double cos_alpha,
                       double sin_alpha,
                       double mat[3][3]);

    void RotateAtoms(double mat[3][3], MIAtomList rot_domain);
    void RotateAtoms(double mat[3][3], MIAtom **atoms, int natoms);

    void RotateAtom(const MIAtom *atom1, const MIAtom *atom2, MIAtom *atom3, float alpha);

    int dTorsion(MIAtom *a1, MIAtom *a2, MIAtom *a3, float dangle, float *dx, float *dy, float *dz);

    void FlipAtoms(Bond&, MIAtomList&);

    bool InvertChiralCenter(MIAtom *center, const std::vector<Bond> &bonds, std::string &error);

    void ReflectAtoms(MIAtomList &atoms);

    void CartesianToRelative(const MIAtom *ref1,
                             const MIAtom *ref2,
                             const MIAtom *ref3,
                             MIAtomList &atoms);

    void RelativeToCartesian(const MIAtom *ref1,
                             const MIAtom *ref2,
                             const MIAtom *ref3,
                             MIAtomList &atoms);
} //namespace chemlib


#endif //TRANSFORM_UTIL_H
