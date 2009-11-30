#include <vector>
#include <sstream>
#include <fstream>

#include <math/mathlib.h>
#include "Substituent.h"
#include "mol_util.h"
#include "Ligand.h"
#include "RingSystem.h"
#include "transform_util.h"

namespace chemlib
{

void Substituent::Clear()
{
    _origin = 0;
    _atoms.clear();
    _first_shell_atoms.clear();
    _direction.clear();
    _issimple = true;
    _iscyclic = false;
}

void Substituent::AddBranch(MIAtom *root)
{
    DepthFirstSearch(root, _origin, _origin, _atoms);
    UpdateFirstShell();
    UpdateFlags();
}

void Substituent::CalcDirection()
{
    _direction.clear();

    if (_issimple)
    {
        _direction.push_back(_first_shell_atoms[0]->x());
        _direction.push_back(_first_shell_atoms[0]->y());
        _direction.push_back(_first_shell_atoms[0]->z());
    }
    else
    {
        _direction.push_back(0.0);
        _direction.push_back(0.0);
        _direction.push_back(0.0);

        double scale;

        MIAtom_iter atm;
        for (atm = _first_shell_atoms.begin();
             atm != _first_shell_atoms.end();
             ++atm)
        {
            scale =  1 / sqrt((*atm)->x() * (*atm)->x()
                              +(*atm)->y() * (*atm)->y()
                              +(*atm)->z() * (*atm)->z());

            _direction[0] += (*atm)->x() * scale;
            _direction[1] += (*atm)->y() * scale;
            _direction[2] += (*atm)->z() * scale;
        }

        _direction[0] /= _first_shell_atoms.size();
        _direction[1] /= _first_shell_atoms.size();
        _direction[2] /= _first_shell_atoms.size();
    }

    CalcTheta();

}

void Substituent::CalcTheta()
{
    _theta = atan2(_direction[1], _direction[0]);   //Returns a value between -PI and PI

    if (_theta < 0)                                 //Adjust theta to run between 0 and 2*PI
    {
        _theta += 2 * PI;
    }
}

void Substituent::Squash()
{
    if (!_issimple && _iscyclic)
    {
        CycleSquash();
    }
    else if (!_issimple && !_iscyclic)
    {
        DoubleSquash();
    }
    else
    {
        SingleSquash();
    }

    CalcDirection();
}

void Substituent::SingleSquash()
{

    XAxisRotate(
        AngleToXYPlane(_first_shell_atoms.front()),
        _atoms
        );
}

void Substituent::DoubleSquash()
{
    MIAtomList planar_atoms;                 //Create vector with the three atoms to
    planar_atoms.push_back(_origin);                    //be rotated into the XY-plane, the origin
    planar_atoms.push_back(_first_shell_atoms[0]);      //atom and the first two neighbors
    planar_atoms.push_back(_first_shell_atoms[1]);

    std::vector<double> normal;
    std::vector<double> rot_axis;
    std::vector<double> cos_sin;
    double rot_mat[3][3];

    CalcNormal(planar_atoms, normal);
    CalcRotationToZAxis(normal, rot_axis, cos_sin);
    CalcRotMatrix(rot_axis, cos_sin[0], cos_sin[1], rot_mat);
    RotateAtoms(rot_mat, _atoms);

    //	RotateIntoXYPlane(planar_atoms, _atoms);		//Rotates these three atoms into XY-plane
}

void Substituent::CycleSquash()
{

    const RingSystem &rs = _lig->ringsystems[_origin->ring_system()]; //Fetch the ringsystem object

    std::vector<double> normal;
    rs.GetNormal(normal);                           //Get a normal to the plane of the ring atoms

    std::vector<double> rot_axis;
    std::vector<double> cos_sin;
    double rot_mat[3][3];

    CalcRotationToZAxis(normal, rot_axis, cos_sin);
    CalcRotMatrix(rot_axis, cos_sin[0], cos_sin[1], rot_mat);
    RotateAtoms(rot_mat, _atoms);

    //	RotateIntoXYPlane(normal, _atoms);				//Rotate all the ring atoms into the XY-plane
}

void Substituent::UpdateFlags()
{
    if (_first_shell_atoms.size() < 2)
    {
        _issimple = true;
    }
    else
    {
        _issimple = false;
    }

    if (_origin->iscyclic() && _first_shell_atoms[0]->iscyclic()
        && _origin->ring_system() == _first_shell_atoms[0]->ring_system())
    {
        _iscyclic = true;
    }
    else
    {
        _iscyclic = false;
    }
}

void Substituent::UpdateFirstShell()
{

    _first_shell_atoms.clear();                 //For simplicity, start with a clean slate

    MIAtom_iter atm;
    for (atm = _atoms.begin(); atm != _atoms.end(); ++atm)       //Loop over atoms in substituent
    {
        if (std::find(_origin->nabors().begin(),
                      _origin->nabors().end(),                //Check if this atom is a neighbor
                      *atm) != _origin->nabors().end())       //of the origin atom

        {
            _first_shell_atoms.push_back(*atm);         //Store as a first-shell atom
        }
    }
}

void Substituent::ZAxisRotate(double theta)
{

    double cos_theta = cos(theta);              //Store the rotation paramters
    double sin_theta = sin(theta);

    double x;
    MIAtom_iter atm;

    for (atm = _atoms.begin(); atm != _atoms.end(); ++atm)
    {
        x = (*atm)->x();
        (*atm)->setPosition((float)(cos_theta * x - sin_theta * (*atm)->y()),
                            (float)(sin_theta * x + cos_theta * (*atm)->y()),
                            0.0f);
    }
}

}
