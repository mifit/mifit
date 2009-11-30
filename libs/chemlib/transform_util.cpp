#include <cmath>
#include <vector>
#include <sstream>
#include <fstream>
#include <functional>

#include <math/mathlib.h>


#include "MIAtom_fwd.h"
#include "Residue.h"
#include "mol_util.h"
#include "atom_util.h"
#include "util.h"
#include "transform_util.h"
#include "math_util.h"

using namespace std;

#define X 0
#define Y 1
#define Z 2

#if defined(_WIN32) || defined(__APPLE__)
#define TINY 1.0E-7
#else
#define TINY 1.0E-12
#endif

namespace chemlib
{
double AngleToXYPlane(MIAtom *atom)
{

    if (fabs(atom->z()) < TINY)
    {
        return 0.0;
    }

    double theta = atan2(atom->z(), atom->y());         //Measure the necessary rotation angle

    if (theta > PI / 2.0)                                   //Choose the smallest rotation
    {
        theta -= PI;                                        //that aligns with the XY-plane
    }
    if (theta < -PI / 2.0)
    {
        theta += PI;
    }

    return -theta;
}

void XAxisRotate(double theta, MIAtomList &domain)
{

    double cos_theta = cos(theta);                          //Store the rotation paramters
    double sin_theta = sin(theta);
    double y;                                               //Only need to store a temp for y,
    //since z can be modified in-place
    MIAtom_iter atm;
    for (atm = domain.begin(); atm != domain.end(); ++atm)
    {
        y = (*atm)->y();
        (*atm)->setPosition((*atm)->x(),
                            (float)(cos_theta * y - sin_theta * (*atm)->z()),
                            (float)(sin_theta * y + cos_theta * (*atm)->z()));
    }
}

/////////////////////////////////////////////////////////////////////////////
// Function:    RotateIntoXYPlane
// Purpose:		Rotates a std::vector of atoms around the x-axis, aligning a given
//				atom with the XY-plane (i.e. setting z=0)
// Input:       A pointer to the atom to align, and a std::vector of ptrs to atoms to rotate
// Output:      New coordinates written to x,y, and z members of ATOM structs
// Requires:
/////////////////////////////////////////////////////////////////////////////
void RotateIntoXYPlane(MIAtom *ref, MIAtomList &branch)
{
    if (fabs(ref->z()) < TINY)
    {
        return;
    }

    double theta = atan2(ref->z(), ref->y());           //Measure the necessary rotation angle

    if (theta > PI / 2.0)                                   //Choose the smallest rotation
    {
        theta -= PI;                                        //that aligns with the XY-plane
    }
    if (theta < -PI / 2.0)
    {
        theta += PI;
    }

    double cos_theta = cos(-theta);                         //Store the rotation paramters
    double sin_theta = sin(-theta);
    double y;                                               //Only need to store a temp for y,
    //since z can be modified in-place

    MIAtom_iter atm;
    for (atm = branch.begin(); atm != branch.end(); ++atm)
    {
        y = (*atm)->y();
        (*atm)->setPosition((*atm)->x(),
                            (float)(cos_theta * y - sin_theta * (*atm)->z()),
                            (float)(sin_theta * y + cos_theta * (*atm)->z()));
    }
}

void CalcNormal(MIAtomList &ref, std::vector<double> &normal)
{
    double v1[3];
    double v2[3];

    v1[0] = ref[1]->x() - ref[0]->x();                              //Calc the first atom-atom std::vector
    v1[1] = ref[1]->y() - ref[0]->y();
    v1[2] = ref[1]->z() - ref[0]->z();

    v2[0] = ref[2]->x() - ref[0]->x();                              //Calc the second atom-atom std::vector
    v2[1] = ref[2]->y() - ref[0]->y();
    v2[2] = ref[2]->z() - ref[0]->z();

    normal.clear();
    normal.push_back(v1[1] * v2[2] - v1[2] * v2[1]);            //v1 X v2 is the normal std::vector
    normal.push_back(v1[2] * v2[0] - v1[0] * v2[2]);            //to the three atoms
    normal.push_back(v1[0] * v2[1] - v1[1] * v2[0]);

    if (VectLength(normal) < TINY)
    {
        normal.clear();                                         //just use the X-axis to rotate
        normal.push_back(1);                                    //the atoms into the XY-plane
        normal.push_back(0);
        normal.push_back(0);
    }
}

void CalcRotationToZAxis(std::vector<double> &normal,
                         std::vector<double> &axis,
                         std::vector<double> &cos_sin)
{
    double x = normal[0];
    double y = normal[1];

    if (fabs(x) < TINY && fabs(y) < TINY)
    {
        axis.push_back(1);
        axis.push_back(0);
        axis.push_back(0);
        cos_sin.push_back(1);
        cos_sin.push_back(0);
        return;
    }

    double n_length = VectLength(normal);

    if (n_length < TINY)
    {
        axis.push_back(1);
        axis.push_back(0);
        axis.push_back(0);
        cos_sin.push_back(1);
        cos_sin.push_back(0);
        return;
    }

    //	if (cal < 0.0) {
    //		ScaleVect(normal, -1.0);
    //	}

    //	cal *= -1.0;
    //	double sal = -sqrt(1 - cal*cal);
    //
    //	double vx, vy;
    double alength = sqrt(x * x                     //Roll-your-own normalization
                          +y * y); //of v to unit length
    axis.clear();
    axis.reserve(3);
    axis.push_back(y / alength);                    //Rotation axis is the cross product of the normal
    axis.push_back(-x / alength);                   //std::vector with the positive z-axis
    axis.push_back(0);                              //i.e. (nx, ny, nz) X (0, 0, 1) = (ny, -nx, 0)

    cos_sin.clear();
    double cal = -normal[2] / n_length;
    cos_sin.push_back(cal);
    cos_sin.push_back(-sqrt(1 - cal*cal));
}

void CalcRotMatrix(std::vector<double> rot_axis,
                   double cos_alpha,
                   double sin_alpha,
                   double mat[3][3])
{

    double vx = rot_axis[0];
    double vy = rot_axis[1];
    //double vz = rot_axis[2];

    mat[X][X] = cos_alpha + (1-cos_alpha)*vx*vx;                    //cos_alphac the rotation matrix
    mat[X][Y] = (1-cos_alpha)*vx*vy;                                //that rotates through an
    mat[X][Z] = sin_alpha*vy;                                       //angle with sin and cos of
    mat[Y][X] = (1-cos_alpha)*vx*vy;                                //"sin_alpha" and "cos_alpha" about an
    mat[Y][Y] = cos_alpha + (1-cos_alpha)*vy*vy;                    //axis pointing in the
    mat[Y][Z] = -sin_alpha*vx;                                      //direction of the unit std::vector
    mat[Z][X] = -sin_alpha*vy;                                      //(vx,vy,0)
    mat[Z][Y] = sin_alpha*vx;
    mat[Z][Z] = cos_alpha;
}

void CalcRotMatrix(double vx,
                   double vy,
                   double /* vz */,
                   double cos_alpha,
                   double sin_alpha,
                   double mat[3][3])
{


    mat[X][X] = cos_alpha + (1-cos_alpha)*vx*vx;                    //cos_alphac the rotation matrix
    mat[X][Y] = (1-cos_alpha)*vx*vy;                                //that rotates through an
    mat[X][Z] = sin_alpha*vy;                                       //angle with sin and cos of
    mat[Y][X] = (1-cos_alpha)*vx*vy;                                //"sin_alpha" and "cos_alpha" about an
    mat[Y][Y] = cos_alpha + (1-cos_alpha)*vy*vy;                    //axis pointing in the
    mat[Y][Z] = -sin_alpha*vx;                                      //direction of the unit std::vector
    mat[Z][X] = -sin_alpha*vy;                                      //(vx,vy,0)
    mat[Z][Y] = sin_alpha*vx;
    mat[Z][Z] = cos_alpha;
}

void CalcRotMatrix(double *rot_axis,
                   double cos_alpha,
                   double sin_alpha,
                   double mat[3][3])
{

    double vx = rot_axis[0];
    double vy = rot_axis[1];
    //double vz = rot_axis[2];

    mat[X][X] = cos_alpha + (1-cos_alpha)*vx*vx;                    //cos_alphac the rotation matrix
    mat[X][Y] = (1-cos_alpha)*vx*vy;                                //that rotates through an
    mat[X][Z] = sin_alpha*vy;                                       //angle with sin and cos of
    mat[Y][X] = (1-cos_alpha)*vx*vy;                                //"sin_alpha" and "cos_alpha" about an
    mat[Y][Y] = cos_alpha + (1-cos_alpha)*vy*vy;                    //axis pointing in the
    mat[Y][Z] = -sin_alpha*vx;                                      //direction of the unit std::vector
    mat[Z][X] = -sin_alpha*vy;                                      //(vx,vy,0)
    mat[Z][Y] = sin_alpha*vx;
    mat[Z][Z] = cos_alpha;
}

void RotateAtoms(double mat[3][3],
                 MIAtomList rot_domain)
{
    double x, y, z;
    MIAtom_iter atm;
    for (atm = rot_domain.begin(); atm != rot_domain.end(); ++atm)
    {
        x = (*atm)->x();
        y = (*atm)->y();
        z = (*atm)->z();

        (*atm)->setPosition((float)(x*mat[X][X] + y*mat[X][Y] + z*mat[X][Z]),               //Rotate the atom
                            (float)(x*mat[Y][X] + y*mat[Y][Y] + z*mat[Y][Z]),
                            (float)(x*mat[Z][X] + y*mat[Z][Y] + z*mat[Z][Z]));
    }
}

void RotateAtoms(double mat[3][3],
                 MIAtom **atoms,
                 int natoms)
{
    double x, y, z;
    for (int i = 0; i < natoms; ++i)
    {
        x = atoms[i]->x();
        y = atoms[i]->y();
        z = atoms[i]->z();

        atoms[i]->setPosition((float)(x*mat[X][X] + y*mat[X][Y] + z*mat[X][Z]),             //Rotate the atom
                              (float)(x*mat[Y][X] + y*mat[Y][Y] + z*mat[Y][Z]),
                              (float)(x*mat[Z][X] + y*mat[Z][Y] + z*mat[Z][Z]));
    }
}

/////////////////////////////////////////////////////////////////////////////
// Function:    RotateIntoXYPlane
// Purpose:		Rotates a set of atoms, aligning a given std::vector with the Z-axis
// Input:       A std::vector of doubles to aligne, and a std::vector of ptrs to atoms to rotate
// Output:      New coordinates written to x,y, and z members of ATOM structs
// Requires:
/////////////////////////////////////////////////////////////////////////////
void RotateIntoXYPlane(std::vector<double> &normal, MIAtomList &branch)
{
    if (normal[0] == 0.0 && normal[1] == 0.0)
    {
        return;
    }
    double n_length = VectLength(normal);
    if (n_length == 0.0)
    {
        return;
    }

    double cal = -normal[2] / n_length;
    //	if (cal < 0.0) {
    //		ScaleVect(normal, -1.0);
    //	}

    //	cal *= -1.0;
    double sal = -sqrt(1 - cal*cal);

    double vx, vy;
    vx = normal[1];                         //Rotation axis v is the cross product of the normal
    vy = -normal[0];                                //std::vector with the positive z-axis
    //i.e. (nx, ny, nz) X (0, 0, 1) = (ny, -nx, 0)


    double length;
    length = sqrt(vy * vy                   //Roll-your-own normalization
                  +vx * vx);        //of v to unit length

    if (length == 0)
    {
        return;
    }

    vx /= length;
    vy /= length;

    double mat[3][3];
    mat[X][X] = cal + (1-cal)*vx*vx;                        //Calc the rotation matrix
    mat[X][Y] = (1-cal)*vx*vy;                              //that rotates through an
    mat[X][Z] = sal*vy;                                     //angle with sin and cos of
    mat[Y][X] = (1-cal)*vx*vy;                              //"sal" and "cal" about an
    mat[Y][Y] = cal + (1-cal)*vy*vy;                        //axis pointing in the
    mat[Y][Z] = -sal*vx;                                    //direction of the unit std::vector
    mat[Z][X] = -sal*vy;                                    //(vx,vy,0)
    mat[Z][Y] = sal*vx;
    mat[Z][Z] = cal;

    double x, y, z;
    MIAtom_iter atm;
    for (atm = branch.begin(); atm != branch.end(); ++atm)
    {
        x = (*atm)->x();
        y = (*atm)->y();
        z = (*atm)->z();

        (*atm)->setPosition((float)(x*mat[X][X] + y*mat[X][Y] + z*mat[X][Z]),               //Rotate the atom
                            (float)(x*mat[Y][X] + y*mat[Y][Y] + z*mat[Y][Z]),
                            (float)(x*mat[Z][X] + y*mat[Z][Y] + z*mat[Z][Z]));
    }
}

/////////////////////////////////////////////////////////////////////////////
// Function:    RotateIntoXYPlane
// Purpose:		Rotates a set of atoms, aligning a given std::vector of three
//				atoms with the XY-plane (i.e. setting z=0)
// Input:       A std::vector of ptr to atoms to align, and a std::vector of ptrs to atoms to rotate
// Output:      New coordinates written to x,y, and z members of ATOM structs
// Requires:
/////////////////////////////////////////////////////////////////////////////
void RotateIntoXYPlane(MIAtomList &ref, MIAtomList &branch)
{
    std::vector<double> normal;
    double v1[3];
    double v2[3];

    v1[0] = ref[1]->x() - ref[0]->x();                              //Calc the first atom-atom std::vector
    v1[1] = ref[1]->y() - ref[0]->y();
    v1[2] = ref[1]->z() - ref[0]->z();

    v2[0] = ref[2]->x() - ref[0]->x();                              //Calc the second atom-atom std::vector
    v2[1] = ref[2]->y() - ref[0]->y();
    v2[2] = ref[2]->z() - ref[0]->z();

    normal.push_back(v1[1] * v2[2] - v1[2] * v2[1]);            //v1 X v2 is the normal std::vector
    normal.push_back(v1[2] * v2[0] - v1[0] * v2[2]);            //to the three atoms
    normal.push_back(v1[0] * v2[1] - v1[1] * v2[0]);

    if (VectLength(normal) > TINY)
    {
        RotateIntoXYPlane(normal, branch);
    }
    else                                                        //If the atoms are co-linear,
    {
        normal.clear();                                         //just use the X-axis to rotate
        normal.push_back(1);                                    //the atoms into the XY-plane
        normal.push_back(0);
        normal.push_back(0);
        RotateIntoXYPlane(normal, branch);
    }
}

/////////////////////////////////////////////////////////////////////////////
// Function:    RotateIntoXYPlane
// Purpose:		Rotates a set of atoms, aligning a given std::vector with the Z-axis
// Input:       A std::vector of doubles to align, and an array of ptrs to atoms
// Output:      New coordinates written to x,y, and z members of ATOM structs
// Requires:
/////////////////////////////////////////////////////////////////////////////
void RotateIntoXYPlane(double *normal, MIAtom **atoms, int natoms)
{
    if (normal[0] == 0.0 && normal[1] == 0.0)
    {
        return;
    }
    double n_length = VectLength(normal);
    if (n_length == 0.0)
    {
        return;
    }

    double cal = -normal[2] / n_length;
    //	if (cal < 0.0) {
    //		ScaleVect(normal, -1.0);
    //	}

    //	cal *= -1.0;
    double sal = -sqrt(1 - cal*cal);

    double vx, vy;
    vx = normal[1];                                 //Rotation axis v is the cross product of the normal
    vy = -normal[0];                                //std::vector with the positive z-axis
    //i.e. (nx, ny, nz) X (0, 0, 1) = (ny, -nx, 0)

    double length;
    length = sqrt(vy * vy                   //Roll-your-own normalization
                  +vx * vx);                //of v to unit length

    if (length == 0)
    {
        return;
    }

    vx /= length;
    vy /= length;

    double mat[3][3];
    mat[X][X] = cal + (1-cal)*vx*vx;                        //Calc the rotation matrix
    mat[X][Y] = (1-cal)*vx*vy;                              //that rotates through an
    mat[X][Z] = sal*vy;                                     //angle with sin and cos of
    mat[Y][X] = (1-cal)*vx*vy;                              //"sal" and "cal" about an
    mat[Y][Y] = cal + (1-cal)*vy*vy;                        //axis pointing in the
    mat[Y][Z] = -sal*vx;                                    //direction of the unit std::vector
    mat[Z][X] = -sal*vy;                                    //(vx,vy,0)
    mat[Z][Y] = sal*vx;
    mat[Z][Z] = cal;

    double x, y, z;
    for (int i = 0; i < natoms; ++i)
    {
        x = atoms[i]->x();
        y = atoms[i]->y();
        z = atoms[i]->z();

        atoms[i]->setPosition((float)(x*mat[X][X] + y*mat[X][Y] + z*mat[X][Z]),             //Rotate the atom
                              (float)(x*mat[Y][X] + y*mat[Y][Y] + z*mat[Y][Z]),
                              (float)(x*mat[Z][X] + y*mat[Z][Y] + z*mat[Z][Z]));
    }
}

void FitToXYPlane(Ligand *mol)
{
    double p_normal[3];             //vector normal to best-fit plane
    double p_displace;              //displacement from origin of best-fit plane

    mol->LSqrPlane(p_normal, &p_displace);

    MIAtom **atoms = new MIAtom *[mol->GetNumAtoms()];

    int i, atm_count = 0;
    std::vector<Residue*>::iterator ri;
    for (ri = mol->residues.begin(); ri != mol->residues.end(); ++ri)
    {
        Residue *res = *ri;
        for (unsigned int i = 0; i < res->atoms().size(); ++i)
        {
            atoms[atm_count] = res->atom(i);
            atm_count++;
        }
    }

    RotateIntoXYPlane(p_normal, atoms, atm_count);

    for (i = 0; i < atm_count; ++i)             //Translate to z=0
    {
        atoms[i]->setPosition(atoms[i]->x(), atoms[i]->y(), 0.0f);
    }
}

void TranslateToOrigin(MIAtom *ref, Ligand *mol)
{
    double dx = ref->x();
    double dy = ref->y();
    double dz = ref->z();

    std::vector<Residue*>::iterator ri;
    for (ri = mol->residues.begin(); ri != mol->residues.end(); ++ri)
    {
        Residue *res = *ri;
        for (unsigned int i = 0; i < res->atoms().size(); i++)
        {
            res->atom(i)->translate(-(float)dx, -(float)dy, -(float)dz);
        }
    }

    ref->setPosition(0.0f, 0.0f, 0.0f);
}

void RotateToXAxis(MIAtom *ref, Ligand *mol)
{

    double bond_length = DistToOrigin(ref);         //Store the distance from this atom to the
    if (bond_length == 0.0)                         //origin
    {
        return;
    }

    double vy, vz;
    vy = -ref->z();                               //Rotation axis v is the cross product of the
    vz = ref->y();                                //atom coords with the positive x-axis
    //i.e. (x, y, z) X (1, 0, 0) = (0, -z, y)

    double length;
    length = sqrt(vy * vy                           //Roll-your-own normalization
                  +vz * vz);                //of v to unit length

    double cal, sal;
    if (length < TINY)
    {
        vy = 0.0;
        vz = 1.0;
        cal = -1.0;
        sal = 0.0;
    }
    else
    {
        vy /= length;
        vz /= length;
        //Determine the length of the rotation by
        cal = ref->x() / bond_length;                 //calculating the cosine and sine of the
        sal = -sqrt(1 - cal*cal);                   //angle with the positive x-axis (i.e. (1,0,0))
    }

    double mat[3][3];
    mat[X][X] = cal;                                            //Calc the rotation matrix
    mat[X][Y] = -sal*vz;                                        //that rotates through an
    mat[X][Z] = sal*vy;                                         //angle with sin and cos of
    mat[Y][X] = sal*vz;                                         //"sal" and "cal" about an
    mat[Y][Y] = cal + (1-cal)*vy*vy;                            //axis pointing in the
    mat[Y][Z] = (1-cal)*vy*vz;                                  //direction of the unit std::vector
    mat[Z][X] = -sal*vy;                                        //(vx,vy,vz)
    mat[Z][Y] = (1-cal)*vy*vz;
    mat[Z][Z] = cal + (1-cal)*vz*vz;

    double x, y, z;
    std::vector<Residue*>::iterator ri;
    for (ri = mol->residues.begin(); ri != mol->residues.end(); ++ri)
    {
        Residue *res = *ri;
        for (unsigned int i = 0; i < res->atoms().size(); i++)
        {

            x = res->atom(i)->x();
            y = res->atom(i)->y();
            z = res->atom(i)->z();

            res->atom(i)->setPosition((float)(x*mat[X][X] + y*mat[X][Y] + z*mat[X][Z]),            //Rotate the atom
                                      (float)(x*mat[Y][X] + y*mat[Y][Y] + z*mat[Y][Z]),
                                      (float)(x*mat[Z][X] + y*mat[Z][Y] + z*mat[Z][Z]));
        }
    }
}

//***RotateAtom***//
//Rotates MIAtom atom3 by alpha degrees around the axis formed by atom1 and atom2
//The direction of rotation is counterclockwise, when looking from atom2 to atom3
//The new coordinates are written to the x,y,z members of atom3
void RotateAtom(const MIAtom *atom1, const MIAtom *atom2, MIAtom *atom3, float alpha)
{
    double cal, sal;                                                //Temps to hold cos and sin
    double vx, vy, vz;                                              //Unit vector for axis of rotat
    double normfactor;                                              //Length of atom1->getAtom2() bond
    double mat[3][3];

    cal = cos(alpha * RAD2DEG);                                     //Store cos and sin of alpha
    sal = sin(alpha * RAD2DEG);

    normfactor = AtomDist(*atom1, *atom2);                            //Calc the bond length

    vx = (atom2->x() - atom1->x()) / normfactor;                //Get unit vector in direction
    vy = (atom2->y() - atom1->y()) / normfactor;                //from atom1 to atom2
    vz = (atom2->z() - atom1->z()) / normfactor;

    mat[X][X] = cal + (1-cal)*vx*vx;                                //Calc the rotation matrix
    mat[X][Y] = (1-cal)*vx*vy - sal*vz;                             //that rotates through an
    mat[X][Z] = (1-cal)*vx*vz + sal*vy;                             //angle with sin and cos of
    mat[Y][X] = (1-cal)*vx*vy + sal*vz;                             //"sal" and "cal" about an
    mat[Y][Y] = cal + (1-cal)*vy*vy;                                //axis pointing in the
    mat[Y][Z] = (1-cal)*vy*vz - sal*vx;                             //direction of the unit vector
    mat[Z][X] = (1-cal)*vx*vz - sal*vy;                             //(vx,vy,vz)
    mat[Z][Y] = (1-cal)*vy*vz + sal*vx;
    mat[Z][Z] = cal + (1-cal)*vz*vz;

    double x = atom3->x() - atom1->x();                             //Translate so that atom1, the
    double y = atom3->y() - atom1->y();                             //base of the rotation vector,
    double z = atom3->z() - atom1->z();                             //is at the origin

    atom3->setPosition((float)(x*mat[X][X] + y*mat[X][Y] + z*mat[X][Z]),                //Rotate the atom
                       (float)(x*mat[Y][X] + y*mat[Y][Y] + z*mat[Y][Z]),
                       (float)(x*mat[Z][X] + y*mat[Z][Y] + z*mat[Z][Z]));

    atom3->translate(atom1->x(), atom1->y(), atom1->z());                  //Translate back to the
    //original reference frame
}

int dTorsion(MIAtom *a1,
             MIAtom *a2,
             MIAtom *a3,
             float dangle,
             float *dx,
             float *dy,
             float *dz)
{

    float orig_a3[3] = {a3->x(), a3->y(), a3->z()};       //Store the original coords

    RotateAtom(a1, a2, a3, dangle);                             //Rotate atom3 to the desired angle

    *dx = a3->x() - orig_a3[0];                               //Measure the translation (dx,dy,dz)
    *dy = a3->y() - orig_a3[1];                               //that equals this rotation
    *dz = a3->z() - orig_a3[2];

    a3->setPosition(orig_a3[0],                                     //Restore the original coordinates
                    orig_a3[1],
                    orig_a3[2]);

    return 1;
}

void FlipAtoms(Bond &bond, MIAtomList &atoms)
{
    double refl[2][2];

    double bv[2];
    bv[0] = bond.getAtom2()->x() - bond.getAtom1()->x();
    bv[1] = bond.getAtom2()->y() - bond.getAtom1()->y();

    //Shortcut the normalization of bv to a unit vector, do it on the fly instead.
    //only need the square of the bond_length, since we'll always be
    //multiplying two components of the bond vector togther.

    double sq_bond_length = bv[0] * bv[0] + bv[1] * bv[1];

    //As a shortcut, only calc one of the cross terms, since the
    //matrix for a reflection is symmetric.
    refl[0][0] = (2 * bv[0] * bv[0] / sq_bond_length) - 1;
    refl[0][1] =  2 * bv[0] * bv[1] / sq_bond_length;
    refl[1][1] = (2 * bv[1] * bv[1] / sq_bond_length) - 1;


    double x, y;
    MIAtom_iter atm;
    for (atm = atoms.begin(); atm != atoms.end(); ++atm)
    {
        x = (*atm)->x() - bond.getAtom1()->x();         //Translate so that the bond is based
        y = (*atm)->y() - bond.getAtom1()->y();         //at the origin

        (*atm)->setPosition((float)(refl[0][0] * x + refl[0][1] * y + bond.getAtom1()->x()),          //Reflect and
                            (float)(refl[0][1] * x + refl[1][1] * y + bond.getAtom1()->y()), //translate back
                            0.0f);
    }
}

bool InvertChiralCenter(MIAtom *center, const vector<Bond> &bonds, std::string &error)
{
    //Check if there are an appropriate # of neighbors
    //(Currently only 4 is supported)
    int degree = Degree(center, bonds);

    vector<MIAtom*> core;
    core.push_back(center);
    GetNabors(center, bonds, core);

    vector<MIAtom*> sub1;
    vector<MIAtom*> sub2;

    if (degree == 3)
    {
        double v[3];
        DirectNextBond(center, bonds, v);
        if (VectLength(v) < 0.01)
        {
            error = "Trigonal planar atoms cannot be inverted";    //Eventually replace. Construct perpendicular,
            return false;                                            //add dummy, and go ahead and invert
        }
        MIAtom dummy;
        dummy.setPosition((float)(center->x() + v[0]),
                          (float)(center->y() + v[1]),
                          (float)(center->z() + v[2]));
        core.push_back(&dummy);

        ExclDepthFirstSearch(core[1], center, bonds, sub1);
        ExclDepthFirstSearch(core[2], center, bonds, sub1);
        ExclDepthFirstSearch(core[3], center, bonds, sub2);
    }
    else if (degree == 4)
    {

        ExclDepthFirstSearch(core[1], center, bonds, sub1);
        ExclDepthFirstSearch(core[2], center, bonds, sub1);
        ExclDepthFirstSearch(core[3], center, bonds, sub2);
        ExclDepthFirstSearch(core[4], center, bonds, sub2);
    }
    else
    {
        error = "Chirality inversion is only supported for atoms with 3 or 4 neighbors";
        return false;
    }


    sort(sub1.begin(), sub1.end());
    sort(sub2.begin(), sub2.end());

    if (Has_Intersection(sub1.begin(), sub1.end(),
                         sub2.begin(), sub2.end()))             //Abort if sub1 and sub2 are not
    {
        error = "Could not invert ring atom";                  //disjoint, meaning there is
        return false;                                           //ring bridging the atoms that
    }                                                           //we had hoped to reorient separately

    //Now convert the coordinates for each substituent to
    //coordinates relative to a separate atom triplet
    CartesianToRelative(core[1], core[0], core[2], sub1);
    CartesianToRelative(core[3], core[0], core[4], sub2);

    //Convert the core atoms to their mirror image
    ReflectAtoms(core);

    //Convert the relative coordinates back to cartesian,
    //rebuilding the molecule around the new configuration
    RelativeToCartesian(core[1], core[0], core[2], sub1);
    RelativeToCartesian(core[3], core[0], core[4], sub2);

    return true;
}

void ReflectAtoms(vector<MIAtom*> &atoms)
{
    //Translate so that the first atom is at the origin
    double trans_vect[] = {-atoms[0]->x(), -atoms[0]->y(), -atoms[0]->z()};
    TranslateAtoms3D(trans_vect, atoms);

    //Reflect one axis (happened to choose Y for appearances)
    ReflectY(atoms);

    //Translate back
    trans_vect[0] = -trans_vect[0];
    trans_vect[1] = -trans_vect[1];
    trans_vect[2] = -trans_vect[2];
    TranslateAtoms3D(trans_vect, atoms);
}

void CartesianToRelative(const MIAtom *ref1,
                         const MIAtom *ref2,
                         const MIAtom *ref3,
                         vector<MIAtom*> &atoms)
{
    double b1[3] = {ref1->x() - ref2->x(),
                    ref1->y() - ref2->y(),
                    ref1->z() - ref2->z()};
    double b2[3] = {ref3->x() - ref2->x(),
                    ref3->y() - ref2->y(),
                    ref3->z() - ref2->z()};

    double x_prime[3] = {b1[0], b1[1], b1[2]};
    double y_prime[3];
    double z_prime[3];

    CrossVects(x_prime, b2, y_prime);
    CrossVects(x_prime, y_prime, z_prime);

    NormVect(x_prime);
    NormVect(y_prime);
    NormVect(z_prime);

    double coord_transform[9] = {x_prime[0], x_prime[1], x_prime[2],
                                 y_prime[0], y_prime[1], y_prime[2],
                                 z_prime[0], z_prime[1], z_prime[2]};
    double tv[3] = {-ref2->x(), -ref2->y(), -ref2->z()};

    TranslateAtoms3D(tv, atoms);
    TransformAtoms3D(coord_transform, atoms);

}

void RelativeToCartesian(const MIAtom *ref1,
                         const MIAtom *ref2,
                         const MIAtom *ref3,
                         vector<MIAtom*> &atoms)
{
    double b1[3] = {ref1->x() - ref2->x(),
                    ref1->y() - ref2->y(),
                    ref1->z() - ref2->z()};
    double b2[3] = {ref3->x() - ref2->x(),
                    ref3->y() - ref2->y(),
                    ref3->z() - ref2->z()};

    double x_prime[3] = {b1[0], b1[1], b1[2]};
    double y_prime[3];
    double z_prime[3];

    CrossVects(x_prime, b2, y_prime);
    CrossVects(x_prime, y_prime, z_prime);

    NormVect(x_prime);
    NormVect(y_prime);
    NormVect(z_prime);

    double back_transform[9] = {x_prime[0], y_prime[0], z_prime[0],
                                x_prime[1], y_prime[1], z_prime[1],
                                x_prime[2], y_prime[2], z_prime[2]};
    double tv[3] = {ref2->x(), ref2->y(), ref2->z()};

    //		double back_transform[9];
    //		InvertMatrix3D(coord_transform, back_transform);

    TransformAtoms3D(back_transform, atoms);
    TranslateAtoms3D(tv, atoms);

}

} //namespace chemlib
