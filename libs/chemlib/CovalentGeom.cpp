#include <vector>
#include <sstream>
#include <fstream>

#include "mathlib.h"

#include "CovalentGeom.h"
#include "RingSystem.h"
#include "Dictionary.h"
#include "LigandPerceiver.h"
#include "math_util.h"
#include "mol_util.h"


namespace chemlib {

//Constructor
CovalentGeometry::CovalentGeometry(Ligand* mol, const Residue* res) {
  _lig = mol;                                               //Store pointers to the
  _res = res;                                               //molecule and residue to process
  _cl = mol->geometry;
  _sigmaangle = 0.03F;
  _sigmabond = 0.005F;
  _sigmabump = 1.0F;
}

void CovalentGeometry::AssignResidue() {
  unsigned int i, j;
  MIAtom* atm;
  Bond* bond1, * bond2;
  MIAtom_const_iter nabor1;
  MIAtom_const_iter nabor2;
  std::vector<Bond*> bonds;
  std::vector<Bond*>::iterator bnd;

  _lig->GetResidueBonds(_res, bonds);
  for (bnd = bonds.begin(); bnd != bonds.end(); ++bnd) {
    AssignBondLength(*bnd);
    //		bnd_length = AssignBondLength(*bnd);
    //		bond_lengths[*bnd] = bnd_length;
  }

  double default_angle;
  Angle angle;
  std::vector<Angle> strained;
  std::vector<Angle> angles;
  std::vector<Angle>::iterator ang;

  for (i = 0; i < _res->atoms().size(); ++i) {
    atm = (_res->atom(i));
    angles.clear();
    strained.clear();

    for (nabor1 = atm->nabors().begin(); nabor1 != atm->nabors().end(); ++nabor1) {
      for (nabor2 = nabor1; nabor2 != atm->nabors().end(); ++nabor2) {
        if (*nabor1 == *nabor2) {
          continue;
        }
        InitBondAngle(*nabor1, atm, *nabor2, angle);
        angle.ideal_angle = 0;                      //Marks the angles as yet undefined

        angles.push_back(angle);
      }
    }

    for (ang = angles.begin(); ang != angles.end(); ++ang) {
      CalcStrainAngle(&(*ang));
      if (ang->ideal_angle != 0) {
        strained.push_back(*ang);
      }
    }

    if (strained.size() > 0) {
      if (atm->geom() == TETRAHEDRAL) {
        default_angle = TetraAngleRemainder(strained);              //Get the value in degrees
      } else if (atm->geom() == TRIGONAL_PLANAR) {
        default_angle = TrigAngleRemainder(strained);
      } else {
        default_angle = AngleFromGeom(atm->geom());
      }
    } else {
      default_angle = AngleFromGeom(atm->geom());
    }


    for (ang = angles.begin(); ang != angles.end(); ++ang) {
      if (ang->ideal_angle == 0) {
        ang->ideal_angle = (float)default_angle;
        AngleValue(*ang);
      }
      if (!_lig->GetBond(ang->getAtom1(), ang->atom3)) {
        bond1 = _lig->GetBond(ang->getAtom1(), ang->getAtom2());
        bond2 = _lig->GetBond(ang->getAtom2(), ang->atom3);

        ConvertToDistance(*ang,
          bond1->ideal_length,
          bond2->ideal_length);
        //								  bond_lengths[bond1],
        //								  bond_lengths[bond2]);

        _lig->geometry.AddAngle(*ang);
      }
    }
  }             //Continue loop over atoms

  //Assign the Bumps

  for (i = 0; i < _res->atoms().size(); ++i) {
    for (j = i+1; j < _res->atoms().size(); ++j) {
      AssignBump((_res->atom(i)), (_res->atom(j)));
    }
  }
}

void CovalentGeometry::InitBondAngle(MIAtom* atom1,
                                     MIAtom* atom2,
                                     MIAtom* atom3,
                                     Angle& angle) {
  angle.setAtom1(atom1);
  angle.setAtom2(atom2);
  angle.atom3 = atom3;
  angle.res = _res;

  AngleTolerance(angle);

  SetRingData(angle);
}

void CovalentGeometry::SetRingData(Angle& angle) {
  RingSystem* rs;
  Aromatic* arom;

  //The angle is cyclic only if all of its member atoms are in
  //the same ring system

  if (angle.getAtom1()->iscyclic()
      && angle.getAtom2()->iscyclic()
      && angle.atom3->iscyclic()
      && (angle.getAtom1()->ring_system() == angle.getAtom2()->ring_system())
      && (angle.getAtom2()->ring_system() == angle.atom3->ring_system())) {

    angle.iscyclic = true;
    angle.ring_system = angle.getAtom1()->ring_system();

    rs = &_lig->ringsystems[angle.ring_system];
    angle.smallest_ring_size =  rs->SmallestRing(angle.getAtom1(),
                                 angle.getAtom2(),
                                 angle.atom3);

    angle.isaromatic = _lig->GetBond(angle.getAtom1(), angle.getAtom2() )->isaromatic
                       && _lig->GetBond(angle.getAtom2(), angle.atom3)->isaromatic;

    if (angle.isaromatic) {
      arom = rs->GetAromaticSys(angle.getAtom1(), angle.getAtom2(), angle.atom3);

      angle.smallest_aromatic_ring = arom->SmallestRing(angle.getAtom1(),
                                       angle.getAtom2(),
                                       angle.atom3);
    } else {
      angle.smallest_aromatic_ring = -1;
    }
  } else {
    angle.iscyclic = false;
  }
}

float CovalentGeometry::AssignBondLength(Bond* bond) {
  BondLength bondlength;

  bondlength.setAtom1(bond->getAtom1());
  bondlength.setAtom2(bond->getAtom2());

  LengthTolerance(bondlength);
  LengthValue(bondlength);

  bondlength.dict_include = true;

  _lig->geometry.AddBond(bondlength);

  bond->ideal_length = bondlength.ideal_dist;
  return bondlength.ideal_dist;
}

void CovalentGeometry::AngleTolerance(Angle& angle) {
  angle.tolerance = _sigmaangle;
}

void CovalentGeometry::LengthTolerance(BondLength& bondlength) {
  bondlength.tolerance = _sigmabond;
}

void CovalentGeometry::BumpTolerance(Bump& bump) {
  bump.tolerance = _sigmabump;
}

void CovalentGeometry::ConvertToDistance(Angle& angle, float d1, float d2) {

  angle.ideal_angle = (float)sqrt(d1*d1 + d2*d2 - 2*cos(DEG2RAD*angle.ideal_angle)*d1*d2);

}

//Given a geometry of the central, returns the value the theoretically ideal bond angle
//This disregards ring strain and a host of other factors.  Other strategies will be
//needed for highly accurate values.

float CovalentGeometry::AngleFromGeom(int geometry) {
  if (geometry == LINEAR) {
    return 180.0F;
  } else if (geometry == TRIGONAL_PLANAR) {
    return 120.0F;
  } else if (geometry == TETRAHEDRAL) {
    return 109.47F;
  } else if (geometry == SQUARE_PLANAR) {
    return 90.0F;
  } else if (geometry == TRIGONAL_BIPYRAMIDAL) {
    return 109.47F;
  } else if (geometry == OCTAHEDRAL) {
    return 90.0F;
  }
  //  cout << "Warning: unrecognized geometry!" << endl;
  return 0;
}

void CovalentGeometry::CalcStrainAngle(Angle* angle) {
  if (angle->iscyclic == false) {
    return;
  }

  if (angle->smallest_ring_size == 3) {
    angle->ideal_angle = 60.0;
  } else if (angle->smallest_ring_size == 4) {
    angle->ideal_angle = 90.0;
  } else if (angle->isaromatic
             && angle->smallest_aromatic_ring == 5) {
    angle->ideal_angle = 108.0;
  } else if (angle->smallest_ring_size == 5) {            //Added to make cyclopentanes non-planar
    angle->ideal_angle = 103.5;
  }
}

void CovalentGeometry::AngleValue(Angle& angle) {
  const Bond& bond1 = *(_lig->GetBond(angle.getAtom1(), angle.getAtom2()));
  const Bond& bond2 = *(_lig->GetBond(angle.getAtom2(), angle.atom3));

  AddAngle("C", "N", "C", 2, 2, 2, SINGLEBOND, SINGLEBOND, 128.0F, angle, bond1, bond2);

}

void CovalentGeometry::LengthValue(BondLength& bondlength) {
  _bond = _lig->GetBond(bondlength.getAtom1(), bondlength.getAtom2());

  bondlength.ideal_dist = 0.0F;

  //***Carbon-Carbon Bonds***//
  AddLength("C", "C", 3, 3, SINGLEBOND, 1.53F, bondlength);
  AddLength("C", "C", 3, 2, SINGLEBOND, 1.51F, bondlength);
  AddLength("C", "C", 3, 1, SINGLEBOND, 1.47F, bondlength);
  AddLength("C", "C", 2, 2, DOUBLEBOND, 1.32F, bondlength);
  AddLength("C", "C", 2, 2, PARTIALDOUBLEBOND, 1.39F, bondlength);
  AddLength("C", "C", 2, 2, SINGLEBOND, 1.49F, bondlength);
  AddLength("C", "C", 2, 1, SINGLEBOND, 1.43F, bondlength);
  AddLength("C", "C", 1, 1, TRIPLEBOND, 1.18F, bondlength);

  //***Carbon-Nitrogen Bonds***//
  AddLength("C", "N", 3, 3, SINGLEBOND, 1.49F, bondlength);
  AddLength("C", "N", 3, 2, SINGLEBOND, 1.46F, bondlength);
  //	AddLength("C", "N", 3, 1, SINGLEBOND, 1.41F, bondlength);			//Occurs in isonitrile attchmnts
  AddLength("C", "N", 2, 3, SINGLEBOND, 1.47F, bondlength);
  AddLength("C", "N", 2, 2, DOUBLEBOND, 1.28F, bondlength);
  AddLength("C", "N", 2, 2, PARTIALDOUBLEBOND, 1.35F, bondlength);
  AddLength("C", "N", 2, 2, SINGLEBOND, 1.34F, bondlength);         //Occurs in amides
  //	AddLength("C", "N", 2, 1, SINGLEBOND, , bondlength);				//Occurs in isonitrile attchmnts
  //	AddLength("C", "N", 1, 3, SINGLEBOND, , bondlength);				//Occurs in cyanoamines
  //	AddLength("C", "N", 1, 2, DOUBLEBOND, , bondlength);				//Occurs in isocyanates
  AddLength("C", "N", 1, 1, TRIPLEBOND, 1.14F, bondlength);         //Occurs in nitriles and
  //isonitriles
  //***Carbon-Oxygen Bonds***//
  AddLength("C", "O", 3, 3, SINGLEBOND, 1.43F, bondlength);
  AddLength("C", "O", 3, 2, SINGLEBOND, 1.45F, bondlength);             //Occurs in esters
  AddLength("C", "O", 2, 3, SINGLEBOND, 1.36F, bondlength);
  AddLength("C", "O", 2, 2, DOUBLEBOND, 1.22F, bondlength);
  AddLength("C", "O", 2, 2, PARTIALDOUBLEBOND, 1.30F, bondlength);
  AddLength("C", "O", 2, 2, SINGLEBOND, 1.33F, bondlength);
  AddLength("C", "O", 1, 3, SINGLEBOND, 1.21F, bondlength);
  AddLength("C", "O", 1, 2, DOUBLEBOND, 1.17F, bondlength);             //Occurs in isocyanates
  AddLength("C", "O", 1, 1, TRIPLEBOND, 1.00F, bondlength);


  //***Carbon-Sulfur bonds***//
  AddLength("C", "S", 3, 3, SINGLEBOND, 1.79F, bondlength);
  AddLength("C", "S", 3, 2, SINGLEBOND, 1.81F, bondlength);
  AddLength("C", "S", 2, 3, SINGLEBOND, 1.76F, bondlength);
  AddLength("C", "S", 2, 2, SINGLEBOND, 1.77F, bondlength);
  AddLength("C", "S", 2, 2, DOUBLEBOND, 1.68F, bondlength);
  AddLength("C", "S", 2, 2, PARTIALDOUBLEBOND, 1.72F, bondlength);
  AddLength("C", "S", 1, 3, SINGLEBOND, 1.64F, bondlength);         //Occurs in Thiocyanates
  AddLength("C", "S", 1, 3, PARTIALDOUBLEBOND, 1.64F, bondlength);  //Occurs in Thiocyanates
  AddLength("C", "S", 1, 2, DOUBLEBOND, 1.56F, bondlength);         //Occurs in Isothiocyanates

  //***Carbon-Halogen bonds***//
  AddLength("C", "F", 3, 3, SINGLEBOND, 1.35F, bondlength);
  AddLength("C", "F", 2, 3, SINGLEBOND, 1.34F, bondlength);
  //	AddLength("C", "F", 1, 3, SINGLEBOND, 0.0F, bondlength);			//No data here
  AddLength("C", "Cl", 3, 3, SINGLEBOND, 1.78F, bondlength);
  AddLength("C", "Cl", 2, 3, SINGLEBOND, 1.74F, bondlength);
  AddLength("C", "Cl", 1, 3, SINGLEBOND, 1.64F, bondlength);
  AddLength("C", "Br", 3, 3, SINGLEBOND, 1.94F, bondlength);
  AddLength("C", "Br", 2, 3, SINGLEBOND, 1.85F, bondlength);
  AddLength("C", "Br", 1, 3, SINGLEBOND, 1.79F, bondlength);
  AddLength("C", "I", 3, 3, SINGLEBOND, 2.16F, bondlength);
  AddLength("C", "I", 2, 3, SINGLEBOND, 2.10F, bondlength);
  AddLength("C", "I", 1, 3, SINGLEBOND, 1.99F, bondlength);

  //***Carbon-Hydrogen bonds***//
  AddLength("C", "H", 3, 3, SINGLEBOND, 1.09F, bondlength);
  AddLength("C", "H", 2, 3, SINGLEBOND, 1.08F, bondlength);
  AddLength("C", "H", 1, 3, SINGLEBOND, 1.06F, bondlength);

  //***Carbon-Phosporous bonds**//
  AddLength("C", "P", 3, 3, SINGLEBOND, 1.84F, bondlength);
  AddLength("C", "P", 2, 3, SINGLEBOND, 1.82F, bondlength);

  //***Carbon-"Other" bonds***//
  AddLength("C", "B", 3, 3, SINGLEBOND, 1.56F, bondlength);
  AddLength("C", "As", 3, 3, SINGLEBOND, 1.97F, bondlength);
  AddLength("C", "As", 2, 3, SINGLEBOND, 1.96F, bondlength);
  AddLength("C", "Si", 3, 3, SINGLEBOND, 1.87F, bondlength);
  AddLength("C", "Si", 2, 3, SINGLEBOND, 1.84F, bondlength);
  AddLength("C", "Se", 3, 3, SINGLEBOND, 1.97F, bondlength);
  AddLength("C", "Se", 2, 3, SINGLEBOND, 1.93F, bondlength);
  AddLength("C", "Se", 2, 2, PARTIALDOUBLEBOND, 1.86F, bondlength);
  AddLength("C", "Te", 3, 3, SINGLEBOND, 2.16F, bondlength);
  AddLength("C", "Te", 2, 2, PARTIALDOUBLEBOND, 2.12F, bondlength);
  AddLength("C", "Te", 2, 2, DOUBLEBOND, 2.04F, bondlength);
  AddLength("C", "Sb", 2, 3, SINGLEBOND, 2.19F, bondlength);
  AddLength("C", "Au", 1, 3, SINGLEBOND, 2.15F, bondlength);
  AddLength("C", "Pb", 3, 3, SINGLEBOND, 2.31F, bondlength);

  //***Nitrogen-Nitrogen bonds***//
  AddLength("N", "N", 3, 3, SINGLEBOND, 1.43F, bondlength);
  AddLength("N", "N", 3, 2, SINGLEBOND, 1.41F, bondlength);         //Occurs in nitrosamines
  AddLength("N", "N", 2, 2, DOUBLEBOND, 1.24F, bondlength);
  AddLength("N", "N", 2, 2, PARTIALDOUBLEBOND, 1.33F, bondlength);
  AddLength("N", "N", 2, 2, SINGLEBOND, 1.40F, bondlength);             //Occurs in nitrosimines
  AddLength("N", "N", 2, 1, DOUBLEBOND, 1.17F, bondlength);
  AddLength("N", "N", 1, 1, TRIPLEBOND, 1.10F, bondlength);

  //***Nitrogen-Oxygen bonds***//
  AddLength("N", "O", 3, 3, SINGLEBOND, 1.46F, bondlength);             //Occurs in hydroxamines
  AddLength("N", "O", 2, 3, SINGLEBOND, 1.39F, bondlength);             //Occurs in isonitroso
  AddLength("N", "O", 2, 2, DOUBLEBOND, 1.22F, bondlength);
  AddLength("N", "O", 2, 2, PARTIALDOUBLEBOND, 1.41F, bondlength);
  AddLength("N", "O", 2, 2, SINGLEBOND, 1.39F, bondlength);
  AddLength("N", "O", 1, 2, DOUBLEBOND, 1.19F, bondlength);

  //***Nitrogen-Sulfur bonds***//
  AddLength("N", "S", 3, 3, SINGLEBOND, 1.63F, bondlength);
  AddLength("N", "S", 2, 3, SINGLEBOND, 1.68F, bondlength);
  AddLength("N", "S", 2, 2, SINGLEBOND, 1.66F, bondlength);
  AddLength("N", "S", 2, 2, DOUBLEBOND, 1.54F, bondlength);
  AddLength("N", "S", 2, 2, PARTIALDOUBLEBOND, 1.73F, bondlength);

  //***Nitrogen-Halogen bonds***//
  AddLength("N", "F", 3, 3, SINGLEBOND, 1.39F, bondlength);
  AddLength("N", "F", 2, 3, SINGLEBOND, 1.48F, bondlength);
  AddLength("N", "Cl", 3, 3, SINGLEBOND, 1.75F, bondlength);
  AddLength("N", "Cl", 2, 3, SINGLEBOND, 1.9F, bondlength);
  AddLength("N", "Br", 3, 3, SINGLEBOND, 1.84F, bondlength);
  AddLength("N", "Br", 2, 3, SINGLEBOND, 2.14F, bondlength);
  AddLength("N", "I", 3, 3, SINGLEBOND, 2.26F, bondlength);

  //***Nitrogen-Phosphorus bonds***//
  AddLength("N", "P", 3, 3, SINGLEBOND, 1.68F, bondlength);
  AddLength("N", "P", 2, 3, SINGLEBOND, 1.65F, bondlength);
  AddLength("N", "P", 2, 2, SINGLEBOND, 1.72F, bondlength);
  AddLength("N", "P", 2, 3, DOUBLEBOND, 1.59F, bondlength);
  AddLength("N", "P", 2, 2, PARTIALDOUBLEBOND, 1.59F, bondlength);

  //***Nitrogen-Hydrogen bonds***//
  AddLength("N", "H", 3, 3, SINGLEBOND, 1.02F, bondlength);
  AddLength("N", "H", 2, 3, SINGLEBOND, 1.01F, bondlength);

  //***Nitrogen-"Other" bonds***//
  AddLength("N", "Si", 3, 3, SINGLEBOND, 1.75F, bondlength);
  AddLength("N", "Se", 3, 3, SINGLEBOND, 1.85F, bondlength);
  AddLength("N", "Se", 2, 3, SINGLEBOND, 1.81F, bondlength);
  AddLength("N", "Se", 2, 2, DOUBLEBOND, 1.79F, bondlength);
  AddLength("N", "Te", 3, 3, SINGLEBOND, 2.00F, bondlength);
  AddLength("N", "Cu", 2, 3, SINGLEBOND, 1.58F, bondlength);
  AddLength("N", "Cu", 2, 5, SINGLEBOND, 1.58F, bondlength);
  AddLength("N", "Mg", 2, 3, SINGLEBOND, 1.90F, bondlength);
  AddLength("N", "Ru", 2, 3, SINGLEBOND, 2.08F, bondlength);
  AddLength("N", "Zn", 2, 3, SINGLEBOND, 1.78F, bondlength);
  AddLength("N", "Os", 2, 3, SINGLEBOND, 2.09F, bondlength);
  AddLength("N", "Pt", 2, 3, SINGLEBOND, 2.08F, bondlength);
  AddLength("N", "Pt", 3, 3, SINGLEBOND, 2.09F, bondlength);
  AddLength("N", "Rh", 2, 3, SINGLEBOND, 1.90F, bondlength);
  AddLength("N", "Rh", 3, 3, SINGLEBOND, 1.90F, bondlength);
  AddLength("N", "Re", 2, 3, SINGLEBOND, 2.09F, bondlength);
  AddLength("N", "Ir", 3, 3, SINGLEBOND, 2.09F, bondlength);

  //***Oxygen-Oxygen bonds***//
  AddLength("O", "O", 3, 3, SINGLEBOND, 1.47F, bondlength);         //Occurs in peroxides
  AddLength("O", "O", 2, 2, DOUBLEBOND, 1.21F, bondlength);         //Only occurs in molecular O2?
  AddLength("O", "O", 2, 2, PARTIALDOUBLEBOND, 1.28F, bondlength);  //Only occurs in ozone?

  //***Oxygen-Sulfur bonds***//
  AddLength("S", "O", 3, 2, SINGLEBOND, 1.57F, bondlength);
  AddLength("S", "O", 3, 2, PARTIALDOUBLEBOND, 1.47F, bondlength);  //Occurs in sulfoxides
  AddLength("S", "O", 3, 2, DOUBLEBOND, 1.42F, bondlength);
  AddLength("S", "O", 2, 2, DOUBLEBOND, 1.49F, bondlength);         //Occurs in sulfoxides


  //***Oxygen-Halogen bonds***//
  AddLength("O", "F", 3, 3, SINGLEBOND, 1.42F, bondlength);
  AddLength("O", "Cl", 3, 3, SINGLEBOND, 1.70F, bondlength);
  //	AddLength("O", "Cl", 2, 3, DOUBLEBOND, , bondlength);
  //	AddLength("O", "Cl", 2, 3, SINGLEBOND, , bondlength);
  AddLength("O", "Cl", 2, 3, PARTIALDOUBLEBOND, 1.99F, bondlength);
  //	AddLength("O", "Cl", 2, 3, DOUBLEBOND, , bondlength);
  AddLength("O", "I", 3, 3, SINGLEBOND, 2.14F, bondlength);
  //	AddLength("O", "I", 2, 2, DOUBLEBOND, , bondlength);

  //***Oxygen-Hydrogen bonds***//
  AddLength("O", "H", 3, 3, SINGLEBOND, 0.97F, bondlength);
  AddLength("O", "H", 2, 3, SINGLEBOND, 0.97F, bondlength);

  //***Oxygen-Phosphorus bonds***//
  AddLength("P", "O", 3, 3, SINGLEBOND, 1.65F, bondlength);
  AddLength("P", "O", 3, 2, SINGLEBOND, 1.50F, bondlength);
  AddLength("P", "O", 3, 2, PARTIALDOUBLEBOND, 1.50F, bondlength);
  AddLength("P", "O", 3, 2, DOUBLEBOND, 1.47F, bondlength);

  //***Oxygen-"Other" bonds***//
  AddLength("O", "As", 2, 2, DOUBLEBOND, 1.66F, bondlength);
  AddLength("O", "As", 2, 2, SINGLEBOND, 1.71F, bondlength);
  AddLength("O", "B", 3, 3, SINGLEBOND, 1.42F, bondlength);     //Occurs in borates
  //	AddLength("O", "B", 3, 2, SINGLEBOND, , bondlength);			//Occurs in boronic & boric acids
  AddLength("O", "B", 2, 1, DOUBLEBOND, 1.20F, bondlength);
  AddLength("O", "B", 2, 2, DOUBLEBOND, 1.19F, bondlength);
  AddLength("O", "Si", 3, 3, SINGLEBOND, 1.63F, bondlength);
  AddLength("O", "Si", 2, 2, DOUBLEBOND, 1.50F, bondlength);
  AddLength("O", "Se", 3, 3, SINGLEBOND, 1.60F, bondlength);
  AddLength("O", "Se", 2, 2, DOUBLEBOND, 1.97F, bondlength);
  AddLength("O", "Te", 3, 3, SINGLEBOND, 2.13F, bondlength);
  AddLength("O", "Cd", 3, 3, SINGLEBOND, 1.91F, bondlength);
  AddLength("O", "K", 3, 3, SINGLEBOND, 2.36F, bondlength);
  AddLength("O", "Ca", 3, 3, SINGLEBOND, 2.13F, bondlength);
  AddLength("O", "Na", 3, 3, SINGLEBOND, 1.99F, bondlength);
  AddLength("O", "Mg", 3, 3, SINGLEBOND, 1.84F, bondlength);
  AddLength("O", "V", 3, 3, SINGLEBOND, 1.84F, bondlength);
  AddLength("O", "Mn", 3, 3, SINGLEBOND, 1.80F, bondlength);
  AddLength("O", "Co", 3, 3, SINGLEBOND, 1.72F, bondlength);
  AddLength("O", "Cu", 3, 3, SINGLEBOND, 1.50F, bondlength);
  AddLength("O", "Zn", 3, 3, SINGLEBOND, 1.70F, bondlength);
  AddLength("O", "Sb", 3, 3, SINGLEBOND, 2.04F, bondlength);
  AddLength("O", "Re", 3, 3, SINGLEBOND, 1.74F, bondlength);
  AddLength("O", "Re", 2, 3, SINGLEBOND, 2.15F, bondlength);
  AddLength("O", "Re", 2, 2, PARTIALDOUBLEBOND, 2.00F, bondlength);
  AddLength("O", "U", 3, 3, SINGLEBOND, 2.22F, bondlength);
  AddLength("O", "Cu", 2, 2, DOUBLEBOND, 1.24F, bondlength);
  AddLength("O", "W", 2, 2, DOUBLEBOND, 1.75F, bondlength);
  AddLength("O", "W", 2, 2, PARTIALDOUBLEBOND, 2.01F, bondlength);
  AddLength("O", "B", 2, 2, DOUBLEBOND, 1.19F, bondlength);

  //***Sulfur-Sulfur bonds***//
  AddLength("S", "S", 3, 3, SINGLEBOND, 2.04F, bondlength);
  AddLength("S", "S", 3, 2, SINGLEBOND, 1.90F, bondlength);

  //***Sulfur-Halogen bonds***//
  AddLength("S", "F", 3, 3, SINGLEBOND, 1.53F, bondlength);
  AddLength("S", "F", 2, 3, SINGLEBOND, 1.53F, bondlength);
  AddLength("S", "Cl", 3, 3, SINGLEBOND, 2.08F, bondlength);
  AddLength("S", "Cl", 2, 3, SINGLEBOND, 2.08F, bondlength);
  AddLength("S", "Br", 3, 3, SINGLEBOND, 2.27F, bondlength);
  AddLength("S", "Br", 2, 3, SINGLEBOND, 2.27F, bondlength);
  AddLength("S", "I", 2, 3, SINGLEBOND, 2.66F, bondlength);

  //***Sulfur-Hydrogen bonds***//
  AddLength("S", "H", 3, 3, SINGLEBOND, 1.33F, bondlength);

  //***Sulfur-Phosphorus bonds***//
  AddLength("S", "P", 2, 4, DOUBLEBOND, 1.93F, bondlength);

  //***Sulfur-"Other" bonds***//
  AddLength("S", "W", 3, 3, SINGLEBOND, 2.41F, bondlength);
  AddLength("S", "Pb", 3, 3, SINGLEBOND, 2.56F, bondlength);


  //***Phosphorus-Phosphorus bonds***//
  AddLength("P", "P", 3, 3, SINGLEBOND, 2.21F, bondlength);
  AddLength("P", "P", 2, 2, DOUBLEBOND, 2.03F, bondlength);

  //***Halogen-Halogen bonds***//
  AddLength("F", "F", 3, 3, SINGLEBOND, 1.41F, bondlength);
  AddLength("Cl", "Br", 3, 3, SINGLEBOND, 1.99F, bondlength);
  AddLength("Br", "Br", 3, 3, SINGLEBOND, 2.54F, bondlength);
  AddLength("I", "I", 3, 3, SINGLEBOND, 2.92F, bondlength);

  //***Other bonds***//
  AddLength("As", "As", 3, 3, SINGLEBOND, 2.46F, bondlength);
  //	AddLength("P", "B", DOUBLEBOND, 1.55F, bondlength);
  AddLength("P", "F", 3, 3, SINGLEBOND, 1.55F, bondlength);
  AddLength("P", "F", 3, 3, SINGLEBOND, 1.55F, bondlength);
  AddLength("P", "I", 2, 3, SINGLEBOND, 2.49F, bondlength);
  AddLength("P", "Se", 4, 2, SINGLEBOND, 2.09F, bondlength);
  AddLength("P", "Te", 4, 2, SINGLEBOND, 2.34F, bondlength);
  AddLength("P", "Si", 3, 3, SINGLEBOND, 2.26F, bondlength);
  AddLength("Se", "Se", 3, 3, SINGLEBOND, 2.34F, bondlength);
  AddLength("Si", "F", 3, 3, SINGLEBOND, 1.59F, bondlength);
  AddLength("Si", "Si", 3, 3, SINGLEBOND, 2.36F, bondlength);
  AddLength("Te", "F", 3, 3, SINGLEBOND, 1.97F, bondlength);
  AddLength("Te", "Te", 3, 3, SINGLEBOND, 2.73F, bondlength);
  AddLength("Cu", "Cu", 3, 3, SINGLEBOND, 1.74F, bondlength);
  AddLength("Cu", "Cl", 3, 3, SINGLEBOND, 1.86F, bondlength);
  AddLength("Pt", "Cl", 3, 3, SINGLEBOND, 2.36F, bondlength);
  AddLength("Be", "F", 3, 3, SINGLEBOND, 1.39F, bondlength);
  AddLength("Al", "F", 3, 3, SINGLEBOND, 1.71F, bondlength);


  if (bondlength.ideal_dist == 0.0F) {
    bondlength.ideal_dist = CovalentRadius(_bond->getAtom1()->atomicnumber()) +
                            CovalentRadius(_bond->getAtom2()->atomicnumber());
  }

  _bond->ideal_length = bondlength.ideal_dist;

}

void CovalentGeometry::AddLength(const char* a1name,
                                 const char* a2name,
                                 int a1hyb,
                                 int a2hyb,
                                 unsigned char bondorder,
                                 float length,
                                 BondLength& bondlength) {
  char a1symbol[10];
  char a2symbol[10];

  if (strlen(a1name) == 1) {                             //Convert the atom symybols
    a1symbol[0] = ' ';                                 //to a standard form
    a1symbol[1] = toupper(a1name[0]);                   //(All uppercase, right-justified)
    a1symbol[2] = '\0';
  } else {
    a1symbol[0] = toupper(a1name[0]);
    a1symbol[1] = toupper(a1name[1]);
    a1symbol[2] = '\0';
  }

  if (strlen(a2name) == 1) {
    a2symbol[0] = ' ';
    a2symbol[1] = toupper(a2name[0]);
    a2symbol[2] = '\0';
  } else {
    a2symbol[0] = toupper(a2name[0]);
    a2symbol[1] = toupper(a2name[1]);
    a2symbol[2] = '\0';
  }

  int a1 = Atomic_Number(a1symbol);
  int a2 = Atomic_Number(a2symbol);

  if  (_bond->getAtom1()->atomicnumber() == a1                //See if this bondtype matches
       && _bond->getAtom2()->atomicnumber() == a2                 //the current bond
       && _bond->getAtom1()->hybrid() == a1hyb
       && _bond->getAtom2()->hybrid() == a2hyb
       && _bond->getOrder() == bondorder) {
    bondlength.ideal_dist = length;
  }

  if (a1 == a2                                          //Return early if we have a
      && a1hyb == a2hyb) {                              //symmetric bond
    return;
  }

  if  (_bond->getAtom1()->atomicnumber() == a2                //Try switching the direction of
       && _bond->getAtom2()->atomicnumber() == a1                 //the bond.
       && _bond->getAtom1()->hybrid() == a2hyb
       && _bond->getAtom2()->hybrid() == a1hyb
       && _bond->getOrder() == bondorder) {
    bondlength.ideal_dist = length;
  }

}

void CovalentGeometry::AddAngle(const char* a1name,
                                const char* a2name,
                                const char* a3name,
                                int a1hyb,
                                int a2hyb,
                                int a3hyb,
                                unsigned char bondorder1,
                                unsigned char bondorder2,
                                float degrees,
                                Angle& angle,
                                const Bond& bond1,
                                const Bond& bond2) {
  char a1symbol[10];
  char a2symbol[10];
  char a3symbol[10];

  if (strlen(a1name) == 1) {                             //Convert the atom symybols
    a1symbol[0] = ' ';                                 //to a standard form
    a1symbol[1] = toupper(a1name[0]);                   //(All uppercase, right-justified)
    a1symbol[2] = '\0';
  } else {
    a1symbol[0] = toupper(a1name[0]);
    a1symbol[1] = toupper(a1name[1]);
    a1symbol[2] = '\0';
  }

  if (strlen(a2name) == 1) {
    a2symbol[0] = ' ';
    a2symbol[1] = toupper(a2name[0]);
    a2symbol[2] = '\0';
  } else {
    a2symbol[0] = toupper(a2name[0]);
    a2symbol[1] = toupper(a2name[1]);
    a2symbol[2] = '\0';
  }

  if (strlen(a3name) == 1) {
    a3symbol[0] = ' ';
    a3symbol[1] = toupper(a3name[0]);
    a3symbol[2] = '\0';
  } else {
    a3symbol[0] = toupper(a3name[0]);
    a3symbol[1] = toupper(a3name[1]);
    a3symbol[2] = '\0';
  }

  int a1 = Atomic_Number(a1symbol);
  int a2 = Atomic_Number(a2symbol);
  int a3 = Atomic_Number(a3symbol);

  if  (angle.getAtom1()->atomicnumber() == a1                 //See if this bondtype matches
       && angle.getAtom2()->atomicnumber() == a2              //the current bond
       && angle.atom3->atomicnumber() == a3
       && angle.getAtom1()->hybrid() == a1hyb
       && angle.getAtom2()->hybrid() == a2hyb
       && angle.atom3->hybrid() == a3hyb
       && bond1.getOrder() == bondorder1
       && bond2.getOrder() == bondorder2) {
    angle.ideal_angle = degrees;
  }

  if  (angle.getAtom1()->atomicnumber() == a3                 //Try switching the direction of
       && angle.getAtom2()->atomicnumber() == a2              //the bond
       && angle.atom3->atomicnumber() == a1
       && angle.getAtom1()->hybrid() == a3hyb
       && angle.getAtom2()->hybrid() == a2hyb
       && angle.atom3->hybrid() == a1hyb
       && bond1.getOrder() == bondorder2
       && bond2.getOrder() == bondorder1) {
    angle.ideal_angle = degrees;
  }

}

void CovalentGeometry::AssignBump(MIAtom* atom1, MIAtom* atom2) {

  if (_lig->GetBond(atom1, atom2) != 0) {                    //Check if atoms are bonded
    return;
  }

  MIAtom_const_iter nbr1;                     //Check if atoms share a neighbor
  MIAtom_const_iter nbr2;                     //(1-3 interactions)
  for (nbr1 = atom1->nabors().begin(); nbr1 != atom1->nabors().end(); ++nbr1) {
    for (nbr2 = atom2->nabors().begin(); nbr2 != atom2->nabors().end(); ++nbr2) {
      if (*nbr1 == *nbr2) {
        return;
      }
    }
  }

  if (atom1->iscyclic()                                        //Don't track bumps between
      && atom2->iscyclic()                                    //atoms in the same ring system
      && (atom1->ring_system() == atom2->ring_system())) {
    return;
  }

  Bump bump;                                                //Create the bump constraint
  bump.setAtom1(atom1);
  bump.setAtom2(atom2);
  BumpValue(bump);
  BumpTolerance(bump);

  _lig->geometry.AddBump(bump);                             //Add bump to molecule geometry
}

void CovalentGeometry::BumpValue(Bump& bump) {
  bump.min_d = 1.5f * (CovalentRadius(bump.getAtom1()->atomicnumber())
                       + CovalentRadius(bump.getAtom2()->atomicnumber()));
}

//Given a set of fixed angles (in degrees) between bond vectors emanating from a
//tetrahedral center, determine a value (theta, in degrees) for the remaining
//angles so that all the remaining angles can be set with equal value

float TrigAngleRemainder(std::vector<Angle>& fixed) {
  float x;

  if (fixed.size() == 1) {
    x = (360.0f - fixed[0].ideal_angle) / 2.0f;
  } else if (fixed.size() == 2) {
    x = (360.0f - fixed[0].ideal_angle - fixed[1].ideal_angle);
  } else {
    return 120.0f;
  }

  if (x > 180.0) {
    x -= 180.0;
  }

  return x;
}

#define REM_MAX_ITERATIONS 100
#define REM_TOLERANCE 0.000001


//Given a set of fixed angles (in degrees) between bond vectors emanating from a
//tetrahedral center, determine a value (theta, in degrees) for the remaining
//angles so that all the remaining angles can be set with equal value

float TetraAngleRemainder(std::vector<Angle>& fixed) {
  int i, j;

  //Handle the special case of only one specified angle here
  double cos_half;
  if (fixed.size() == 1) {
    cos_half = cos(DEG2RAD * 0.5 * fixed[0].ideal_angle);
    return (float)(2 * RAD2DEG * acos((sqrt(cos_half*cos_half + 8) - cos_half) / 4));
  }

  //Cases with more than 5 angles are either completely specified, or not tetrahedral
  if (fixed.size() > 5) {
    return 109.5;
  }

  //Create a vector of pointers to the substituent atoms covered by the input angles
  //and check that the angles values are valid

  std::vector<Angle>::const_iterator ang;
  MIAtomList sub_atoms;
  for (ang = fixed.begin(); ang != fixed.end(); ++ang) {
    (sub_atoms.push_back(ang->getAtom1()));
    (sub_atoms.push_back(ang->atom3));
    if (!(ang->ideal_angle > 0)
        || !(ang->ideal_angle < 180)) {
      return 109.5;
    }
  }



  //Check that the angles meet certain criteria necessary for the remaining angles
  //to be equal to one another.  If this is impossible, return 109.5.
  std::vector<Angle>::const_iterator ang1, ang2;
  double lower_bound, upper_bound;
  double min_angle = 0.0;
  double max_angle = 180.0;

  for (ang1 = fixed.begin(); ang1 != fixed.end(); ++ang1) {
    lower_bound = ang1->ideal_angle / 2;
    upper_bound = 180 - lower_bound;

    if (lower_bound > min_angle) {
      min_angle = lower_bound;
    }
    if (upper_bound < max_angle) {
      max_angle = upper_bound;
    }
  }

  for (ang1 = fixed.begin(); ang1 != fixed.end(); ++ang1) {
    for (ang2 = fixed.begin(); ang2 != fixed.end(); ++ang2) {
      lower_bound = fabs(ang1->ideal_angle - ang2->ideal_angle);
      upper_bound = ang1->ideal_angle + ang2->ideal_angle;

      if (lower_bound > min_angle) {
        min_angle = lower_bound;
      }
      if (upper_bound < max_angle) {
        max_angle = upper_bound;
      }
    }
  }


  //Remove dupes from the vector of atoms
  MIAtom_iter new_end;

  std::sort(sub_atoms.begin(), sub_atoms.end());
  new_end = std::unique(sub_atoms.begin(), sub_atoms.end());
  sub_atoms.erase(new_end, sub_atoms.end());


  //Return if there aren't between two and four substituent atoms
  if (sub_atoms.size() < 2 || sub_atoms.size() > 4) {
    return 109.5F;
  }

  //Handle the special case of "spiro" centers analytically

  if (fixed.size() == 2 && sub_atoms.size() == 4) {
    return (float)(2 * RAD2DEG * acos(cos(DEG2RAD * 0.5 * fixed[0].ideal_angle) *
                     cos(DEG2RAD * 0.5 * fixed[1].ideal_angle)));
  }

  //Handle the special case where two angles to (one) common bond are fixed
  //		if (fixed.size() == 2 && sub_atoms.size() == 3) {
  //			return Remainder

  //Handle one more special case--three angles spanning only three atoms
  if (fixed.size() == 3 && sub_atoms.size() == 3
      && fixed[0].ideal_angle == fixed[1].ideal_angle
      && fixed[0].ideal_angle == fixed[2].ideal_angle) {

    if (fixed[0].ideal_angle >= 120.0) {
      return 109.5;                                 //To guard sqrt() from negative arguments
    } else {
      return (float)(RAD2DEG *
                     acos(-sqrt((2 * cos(DEG2RAD * fixed[0].ideal_angle) + 1) / 3.0)));
    }
  }


  //If necessary, reorder to assure that the last two atoms do not have a fixed angle
  if (sub_atoms.size() == 4) {
    while (CheckAngle(sub_atoms[2], sub_atoms[3], fixed) > 0) {
      std::next_permutation(sub_atoms.begin(), sub_atoms.end());
    }
  }

  //Create a matrix of pointers to the angle values, and init the upper half
  //to nulls (Don't need last one([2][3]), since we calculate that from the others.)
  double* cos_angle[3][4];

  for (i = 0; i < 2; ++i) {
    for (j = i+1; j < 4; ++j) {
      cos_angle[i][j] = 0;
    }
  }

  //Set the fixed values of the angle
  double cos_fixed_angles[6];
  int nfixed = 0, index1, index3;

  for (ang = fixed.begin(); ang != fixed.end(); ++ang) {
    index1 = std::distance(sub_atoms.begin(),
               std::find(sub_atoms.begin(), sub_atoms.end(), ang->getAtom1()));
    index3 = std::distance(sub_atoms.begin(),
               std::find(sub_atoms.begin(), sub_atoms.end(), ang->atom3));
    if (index1 > 3 || index3 > 3 || index1 == index3
        || (index1 == 2 && index3 == 3)
        || (index3 == 2 && index1 == 2)) {
      //throw error
    }
    if (index1 < index3) {
      cos_fixed_angles[nfixed] = cos(DEG2RAD * ang->ideal_angle);
      cos_angle[index1][index3] = &cos_fixed_angles[nfixed];
      nfixed++;
    } else if (index3 < index1) {
      cos_fixed_angles[nfixed] = cos(DEG2RAD * ang->ideal_angle);
      cos_angle[index3][index1] = &cos_fixed_angles[nfixed];
      nfixed++;
    } else {
      //Throw error
    }
  }

  double cos_theta[3] = {cos(DEG2RAD*upper_bound),
                         cos(DEG2RAD*(upper_bound+lower_bound)/2),
                         cos(DEG2RAD*lower_bound)};

  nfixed = 0;


  for (i = 0; i < 2; ++i) {
    for (j = i+1; j < 4; ++j) {
      if (cos_angle[i][j] == 0) {
        cos_angle[i][j] = &(cos_theta[1]);
      } else {
        nfixed++;
      }
    }
  }


  double v1[3], v2[3], v3[3], v4[3];
  double cos_last_angle;
  double sin_v1_v2 = sqrt(1.0 - *cos_angle[0][1] * *cos_angle[0][1]);
  if (nfixed >= 0 && nfixed <= 5) {
    for (i = 0; i < REM_MAX_ITERATIONS; ++i) {

      v1[0] = 1.0;
      v1[1] = 0.0;
      v1[2] = 0.0;

      v2[0] = *cos_angle[0][1];
      v2[1] = sin_v1_v2;
      v2[2] = 0.0;

      v3[0] = *cos_angle[0][2];
      v3[1] = (*cos_angle[1][2] - (*cos_angle[0][1] * *cos_angle[0][2]))/sin_v1_v2;
      v3[2] = (1 - v3[0]*v3[0] - v3[1]*v3[1] < 0) ?
              0 :
              (sqrt(1 - v3[0]*v3[0] - v3[1]*v3[1]));
      //				if (1 - v3[0]*v3[0] - v3[1]*v3[1] < 0) {
      //					v3[2] = 100.0;
      //				}
      //				else {
      //					v3[2] = -sqrt(1 - v3[0]*v3[0] - v3[1]*v3[1]);
      //				}

      v4[0] = *cos_angle[0][3];
      v4[1] = (*cos_angle[1][3] - (*cos_angle[0][1] * *cos_angle[0][3]))/sin_v1_v2;
      v4[2] = (1 - v4[0]*v4[0] - v4[1]*v4[1] < 0) ?
              0 :
              (-sqrt(1 - v4[0]*v4[0] - v4[1]*v4[1]));
      //				if (1 - v4[0]*v4[0] - v4[1]*v4[1] < 0) {
      //					v4[2] = 100.0;
      //				}
      //				else {
      //					v4[2] = -sqrt(1 - v4[0]*v4[0] - v4[1]*v4[1]);
      //				}

      //				cos_last_angle = (v4[2] < 10.0) ? CosVectorAngle(v3,v4) : 100.0;
      cos_last_angle = CosVectorAngle(v3, v4);

      if (fabs(cos_last_angle - cos_theta[1]) < REM_TOLERANCE) {
        return (float)(RAD2DEG * acos((cos_last_angle + cos_theta[1]) / 2.0));
      }

      if (cos_last_angle > cos_theta[1]) {
        cos_theta[0] = cos_theta[1];
        cos_theta[1] = (cos_theta[0] + cos_theta[2]) / 2;
      } else if (cos_last_angle < cos_theta[1]) {
        cos_theta[2] = cos_theta[1];
        cos_theta[1] = (cos_theta[0] + cos_theta[2]) / 2;
      }

      if (cos_angle[0][1] == &(cos_theta[1])) {
        sin_v1_v2 = sqrt(1.0 - *cos_angle[0][1] * *cos_angle[0][1]);
      }

    }         //End for loop over iterations
    return (float)(RAD2DEG * acos((cos_last_angle + cos_theta[1]) / 2));
  } else if (nfixed == 6) {
    return 109.5F;
  }
  return 109.5F;
}

//Search thru a vector of angles, returning the ideal_angle of the first angle whose
//substituent atoms match the two input atoms
float CheckAngle(MIAtom* atom1, MIAtom* atom2, std::vector<Angle>& angles) {
  std::vector<Angle>::const_iterator ang;
  for (ang = angles.begin(); ang != angles.end(); ++ang) {
    //Check both permutations
    if (ang->getAtom1() == atom1 && ang->atom3 == atom2) {
      return ang->ideal_angle;
    } else if (ang->getAtom1() == atom2 && ang->atom3 == atom1) {
      return ang->ideal_angle;
    }
  }
  return -1;
}

float IdealBondLength(const Bond& bond) {
  return IdealBondLength(bond.getAtom1()->atomicnumber(),
           bond.getAtom2()->atomicnumber(),
           bond.getOrder());
}

float IdealBondLength(int element1, int element2, unsigned char bnd_order) {
  //	MIAtom *atom1 = bond.getAtom1();
  //	MIAtom *atom2 = bond.getAtom2();
  //	unsigned char bnd_order = bond.order;

  //
  //***Carbon-Carbon Bonds***//
  //
  if (element1 == 6 && element2 == 6) {
    switch (bnd_order) {
      case SINGLEBOND:    return 1.51F;
      case PARTIALDOUBLEBOND: return 1.39F;
      case DOUBLEBOND:    return 1.32F;
      case TRIPLEBOND:    return 1.18F;
    }
  }

  //
  //***Carbon-Nitrogen Bonds***//
  //
  if (element1 == 6 && element2 == 7
      || element1 == 7 && element2 == 6) {
    switch (bnd_order) {
      case SINGLEBOND:    return 1.47F;
      case PARTIALDOUBLEBOND: return 1.35F;
      case DOUBLEBOND:    return 1.28F;
      case TRIPLEBOND:    return 1.14F;
    }
  }

  //
  //***Carbon-Oxygen Bonds***//
  //
  if (element1 == 6 && element2 == 8
      || element1 == 8 && element2 == 6) {
    switch (bnd_order) {
      case SINGLEBOND:    return 1.40F;
      case PARTIALDOUBLEBOND: return 1.30F;
      case DOUBLEBOND:    return 1.22F;
      case TRIPLEBOND:    return 1.00F;
    }
  }
  //
  //***Carbon-Sulfur Bonds***//
  //
  if (element1 == 6 && element2 == 16
      || element1 == 16 && element2 == 6) {
    switch (bnd_order) {
      case SINGLEBOND:    return 1.79F;
      case PARTIALDOUBLEBOND: return 1.72F;
      case DOUBLEBOND:    return 1.68F;
    }
  }

  //
  //***Carbon-Fluorine Bonds***//
  //
  if (element1 == 6 && element2 == 9
      || element1 == 9 && element2 == 6) {
    switch (bnd_order) {
      case SINGLEBOND:    return 1.34F;
    }
  }

  //
  //***Carbon-Chlorine Bonds***//
  //
  if (element1 == 6 && element2 == 17
      || element1 == 17 && element2 == 6) {
    switch (bnd_order) {
      case SINGLEBOND:    return 1.76F;
    }
  }

  //
  //***Carbon-Bromine Bonds***//
  //
  if (element1 == 6 && element2 == 35
      || element1 == 35 && element2 == 6) {
    switch (bnd_order) {
      case SINGLEBOND:    return 1.89F;
    }
  }

  //
  //***Carbon-Iodine Bonds***//
  //
  if (element1 == 6 && element2 == 53
      || element1 == 53 && element2 == 6) {
    switch (bnd_order) {
      case SINGLEBOND:    return 2.13F;
    }
  }

  //
  //***Carbon-Hydrogen bonds***//
  //
  if (element1 == 6 && element2 == 1
      || element1 == 1 && element2 == 6) {
    switch (bnd_order) {
      case SINGLEBOND:    return 1.09F;
    }
  }

  //
  //***Carbon-Phosphorous bonds**//
  //
  if (element1 == 6 && element2 == 15
      || element1 == 15 && element2 == 6) {
    switch (bnd_order) {
      case SINGLEBOND:    return 1.83F;
    }
  }

  //
  //***Nitrogen-Nitrogen bonds***//
  //
  if (element1 == 7 && element2 == 7) {
    switch (bnd_order) {
      case SINGLEBOND:    return 1.41F;
      case PARTIALDOUBLEBOND: return 1.33F;
      case DOUBLEBOND:    return 1.24F;
      case TRIPLEBOND:    return 1.10F;
    }
  }

  //
  //***Nitrogen-Oxygen Bonds***//
  //
  if (element1 == 7 && element2 == 8
      || element1 == 8 && element2 == 7) {
    switch (bnd_order) {
      case SINGLEBOND:    return 1.42F;
      case PARTIALDOUBLEBOND: return 1.41F;
      case DOUBLEBOND:    return 1.22F;
    }
  }

  //
  //***Nitrogen-Sulfur Bonds***//
  //
  if (element1 == 7 && element2 == 16
      || element1 == 16 && element2 == 7) {
    switch (bnd_order) {
      case SINGLEBOND:    return 1.65F;
      case PARTIALDOUBLEBOND: return 1.73F;
      case DOUBLEBOND:    return 1.54F;
    }
  }

  //
  //***Nitrogen-Phosphorus Bonds***//
  //
  if (element1 == 7 && element2 == 15
      || element1 == 15 && element2 == 7) {
    switch (bnd_order) {
      case SINGLEBOND:    return 1.68F;
      case PARTIALDOUBLEBOND: return 1.59F;
      case DOUBLEBOND:    return 1.59F;
    }
  }

  //
  //***Nitrogen-Hydrogen bonds***//
  //
  if (element1 == 7 && element2 == 1
      || element1 == 1 && element2 == 7) {
    switch (bnd_order) {
      case SINGLEBOND:    return 1.02F;
    }
  }
  //
  //***Oxygen-Oxygen bonds***//
  //
  if (element1 == 8 && element2 == 8) {
    switch (bnd_order) {
      case SINGLEBOND:    return 1.47F;
      case PARTIALDOUBLEBOND: return 1.28F;
      case DOUBLEBOND:    return 1.21F;
    }
  }
  //
  //***Oxygen-Sulfur bonds***//
  //
  if (element1 == 8 && element2 == 16
      || element1 == 16 && element2 == 8) {
    switch (bnd_order) {
      case SINGLEBOND:    return 1.57F;
      case PARTIALDOUBLEBOND: return 1.47F;
      case DOUBLEBOND:    return 1.44F;
    }
  }


  //
  //***Oxygen-Hydrogen bonds***//
  //
  if (element1 == 8 && element2 == 1
      || element1 == 1 && element2 == 8) {
    switch (bnd_order) {
      case SINGLEBOND:    return 0.97F;
    }
  }

  //
  //***Oxygen-Phosphorus bonds***//
  //
  if (element1 == 8 && element2 == 15
      || element1 == 15 && element2 == 8) {
    switch (bnd_order) {
      case SINGLEBOND:    return 1.50F;
      case PARTIALDOUBLEBOND: return 1.50F;
      case DOUBLEBOND:    return 1.47F;
    }
  }

  //
  //***Sulfur-Sulfur bonds***//
  //
  if (element1 == 16 && element2 == 16) {
    switch (bnd_order) {
      case SINGLEBOND:    return 1.97F;
    }
  }


  //
  //***Sulfur-Hydrogen bonds***//
  //
  if (element1 == 16 && element2 == 1
      || element1 == 1 && element2 == 16) {
    switch (bnd_order) {
      case SINGLEBOND:    return 1.33F;
    }
  }

  //
  //***Sulfur-Phosphorus bonds***//
  //
  if (element1 == 16 && element2 == 15
      || element1 == 15 && element2 == 16) {
    switch (bnd_order) {
      case DOUBLEBOND:    return 1.93F;
    }
  }

  //
  //***Phosphorus-Phosphorus bonds***//
  //
  if (element1 == 15 && element2 == 15) {
    switch (bnd_order) {
      case SINGLEBOND:    return 2.21F;
      case DOUBLEBOND:    return 2.03F;
    }
  }

  return CovalentRadius(element1) +
         CovalentRadius(element2);
}

unsigned char PredictBondOrder(int element1, int element2, float distance) {

  //
  //***Carbon-Carbon Bonds***//
  //
  if (element1 == 6 && element2 == 6) {
    if (distance < 1.23F) {
      return TRIPLEBOND;
    }
    if (distance < 1.32F) {
      return DOUBLEBOND;
    }
    if (distance < 1.42F) {
      return PARTIALDOUBLEBOND;
    }
    return SINGLEBOND;
  }

  //
  //***Carbon-Nitrogen Bonds***//
  //
  if (element1 == 6 && element2 == 7
      || element1 == 7 && element2 == 6) {
    if (distance < 1.21F) {
      return TRIPLEBOND;
    }
    if (distance < 1.28F) {
      return DOUBLEBOND;
    }
    if (distance < 1.41F) {
      return PARTIALDOUBLEBOND;                             //amides unfortunately included here
    }
    return SINGLEBOND;
  }

  //
  //***Carbon-Oxygen Bonds***//
  //
  if (element1 == 6 && element2 == 8
      || element1 == 8 && element2 == 6) {
    if (distance < 1.27F) {
      return DOUBLEBOND;
    }
    if (distance < 1.34F) {
      return PARTIALDOUBLEBOND;
    }
    return SINGLEBOND;
  }
  //
  //***Carbon-Phosphorous bonds**//
  //
  if (element1 == 6 && element2 == 15
      || element1 == 15 && element2 == 6) {
    return SINGLEBOND;
  }

  //
  //***Carbon-Sulfur Bonds***//
  //
  if (element1 == 6 && element2 == 16
      || element1 == 16 && element2 == 6) {
    if (distance < 1.70F) {
      return DOUBLEBOND;
    }
    return SINGLEBOND;
  }

  //
  //***Nitrogen-Nitrogen bonds***//
  //
  if (element1 == 7 && element2 == 7) {
    if (distance < 1.08F) {
      return TRIPLEBOND;
    }
    if (distance < 1.29F) {
      return DOUBLEBOND;
    }
    if (distance < 1.35F) {
      return PARTIALDOUBLEBOND;
    }
    return SINGLEBOND;
  }

  //
  //***Nitrogen-Oxygen Bonds***//
  //
  if (element1 == 7 && element2 == 8
      || element1 == 8 && element2 == 7) {
    if (distance < 1.30F) {
      return DOUBLEBOND;
    }
    if (distance < 1.34F) {
      return PARTIALDOUBLEBOND;
    }
    return SINGLEBOND;
  }

  //
  //***Nitrogen-Phosphorus Bonds***//
  //
  if (element1 == 7 && element2 == 15
      || element1 == 15 && element2 == 7) {
    if (distance < 1.59F) {
      return DOUBLEBOND;
    }
    return SINGLEBOND;
  }

  //
  //***Nitrogen-Sulfur Bonds***//
  //
  if (element1 == 7 && element2 == 16
      || element1 == 16 && element2 == 7) {
    if (distance < 1.57F) {
      return DOUBLEBOND;
    }
    return SINGLEBOND;
  }
  //
  //***Oxygen-Oxygen bonds***//
  //
  if (element1 == 8 && element2 == 8) {
    if (distance < 1.27F) {
      return DOUBLEBOND;                                    //molecular oxygen
    }
    if (distance < 1.34F) {
      return PARTIALDOUBLEBOND;                             //ozone
    }
    return SINGLEBOND;
  }

  //
  //***Oxygen-Phosphorus bonds***//
  //
  if (element1 == 8 && element2 == 15
      || element1 == 15 && element2 == 8) {
    if (distance < 1.52F) {
      return DOUBLEBOND;
    }
    return SINGLEBOND;
  }

  //
  //***Oxygen-Sulfur bonds***//
  //
  if (element1 == 8 && element2 == 16
      || element1 == 16 && element2 == 8) {
    if (distance < 1.53F) {
      return DOUBLEBOND;
    }
    return SINGLEBOND;
  }

  //
  //***Sulfur-Sulfur bonds***//
  //
  if (element1 == 16 && element2 == 16) {
    if (distance < 1.92F) {
      return DOUBLEBOND;                                //This has a problem with -N=SS groups
    }
    return SINGLEBOND;                                  //but they are rare  -kwb
  }

  //
  //***Sulfur-Phosphorus bonds***//
  //
  if (element1 == 16 && element2 == 15
      || element1 == 15 && element2 == 16) {
    if (distance < 1.98F) {
      return DOUBLEBOND;
    }
    return SINGLEBOND;
  }

  //
  //***Phosphorus-Phosphorus bonds***//
  //
  if (element1 == 15 && element2 == 15) {
    if (distance < 2.04F) {
      return DOUBLEBOND;
    }
    return SINGLEBOND;
  }

  //
  //***Carbon-Selenium bonds***//
  //
  if (element1 == 15 && element2 == 15) {
    if (distance < 2.04F) {
      return DOUBLEBOND;
    }
    return SINGLEBOND;
  }
  //

  //***Carbon-Tellurium bonds***//
  //
  if (element1 == 15 && element2 == 15) {
    if (distance < 2.05F) {
      return DOUBLEBOND;
    }
    if (distance < 2.07F) {
      return PARTIALDOUBLEBOND;
    }
    return SINGLEBOND;
  }

  //
  //***Oxygen-Chlorine bonds***//
  //
  if (element1 == 15 && element2 == 15) {
    if (distance < 2.04F) {
      return DOUBLEBOND;
    }
    return SINGLEBOND;
  }

  return SINGLEBOND;                    //Single bonds are the default
}

} //namespace conflib
