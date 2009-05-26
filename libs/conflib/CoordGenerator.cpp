#include "mathlib.h"
#include "chemlib.h"

#include "CoordGenerator.h"
#include "LigRefiner.h"
#include "AtomRefiner.h"
#include "sdg.h"
#include "sdg_parms.h"

using namespace std;
using namespace chemlib;
using namespace conflib;

namespace conflib {

void PositionAtom(MIAtom* atom, const MIAtom* root, double length, const std::vector<Bond>& /* bonds */) {
  switch (root->nabors().size()) {
    case 0:
      atom->setPosition(root->x() + (float)length, root->y(), root->z());
      //		case 1:
      //			PositionAtom( atom, root, &root->nabors[0], length, AngleFromGeom(root->geom) );
      //			if (root->nabors.front().nabors.size() > 0 &&
      //				root->geom != LINEAR)  {
      //				SetExtended( atom, root, &root->nabors[0] );
      //			}
  }
}

/////////////////////////////////////////////////////////////////////////////
// Function:    PositionAtom
// Purpose:		Given an atom object, sets the coordinates to be consistent
//			    with the given internal coordinates
// Input:       Three atoms, generally the first atom is bonded to the new atom
//				The distance in angstroms between the new atom and "root1"
//				The angle in degrees between the new atom, root1, and root2
//				The dihedral angle between all four atoms, around the root1-root2 bond
// Output:      Writes new coordinates to the atom object
// Requires:	Note that if root1, root2, and root3 are collinear atoms, the
//				the atom position is underdetermined, because the dihedral angle
//				is meaningless
/////////////////////////////////////////////////////////////////////////////
void PositionAtom(MIAtom* atom,
                  const MIAtom* root1,
                  const MIAtom* root2,
                  const MIAtom* root3,
                  float length,
                  float angle,
                  float torsion) {

  atom->setPosition(root1->x(), root1->y(), root1->z());                           //Start new atom on top of the previous atom

  double cosba = cos(DEG2RAD * angle);
  double sinba = sin(DEG2RAD * angle);

  double v[3], temp[3];                                 //set v to be the unit vector from the
  BondVector(root1, root2, v);                          //attachment pt to a connected atom
  NormVect(v, v);

  atom->translate((float)(v[0] * cosba * length),                //Move to an arbitrary point
      (float)(v[1] * cosba * length),                //with the correct bond length
      (float)(v[2] * cosba * length));                //and angle.

  temp[0] = -v[1];
  temp[1] = v[0];
  temp[2] = 0;
  NormVect(temp, temp);

  atom->translate((float)(temp[0] * sinba * length),
      (float)(temp[1] * sinba * length),
      0.0f);		//Have chosen delta(z)=0

  double chi;
  chi = CalcAtomTorsion(root3, root2, root1, atom);         //Get the current torsion value

  RotateAtom(root2, root1, atom, torsion-(float)chi);
}

/////////////////////////////////////////////////////////////////////////////
// Function:	PositionAtom
// Purpose:		Given an atom object, sets the coordinates to be consistent with
//				the given bond length and angle
// Input:       Two atoms, with the first atom is bonded to the new atom and the
//				second atom bonded to the first.
//				The distance in angstroms between the new atom and the first atom
//				The angle in degrees between the new atom, the first atom, and the
//				second atom
// Output:		Writes new coordinates to the atom object
// Requires:
/////////////////////////////////////////////////////////////////////////////
void PositionAtom(MIAtom* atom,
                  const MIAtom* root1,
                  const MIAtom* root2,
                  float length,
                  float angle) {
  double cosba = cos(DEG2RAD * angle);
  double sinba = sin(DEG2RAD * angle);

  //Calculate the bond vector: the unit vector along the bond from root1 to root2
  double bv[3];
  BondVector(root1, root2, bv);
  NormVect(bv, bv);

  //Construct an orthogonal vector: a unit vector that is orthogonal to the bond vector
  double ov[3];
  ov[0] = -bv[1];
  ov[1] = bv[0];
  ov[2] = 0;
  NormVect(ov, ov);

  //Position the atom along the vector (cosba * bv + sinba * ov)
  atom->setPosition((float)(root1->x() + length * (bv[0] * cosba + ov[0] * sinba)),
      (float)(root1->y() + length * (bv[1] * cosba + ov[1] * sinba)),
      (float)(root1->z() + length * (bv[2] * cosba)));                     //note that ov[2] = 0
}

}

void CoordGenerator::sdgGenMolStruct() {
  //	GenMolStruct();

  std::vector<sdgDistance> dists;
  std::vector<sdgVolume> vols;
  std::vector<BondLength>::iterator bnd;
  std::vector<Angle>::iterator ang;
  std::vector<Improper>::iterator tor;
  std::vector<Chiral>::iterator chir;
  //	std::vector<MIAtom>::iterator atm;

  for (bnd = _lig->geometry.Bonds.begin(); bnd != _lig->geometry.Bonds.end(); ++bnd) {
    sdgDistance dist(*(bnd->getAtom1()), *(bnd->getAtom2()),
                     bnd->ideal_dist, false);
    dists.push_back(dist);
  }
  for (ang = _lig->geometry.Angles.begin(); ang != _lig->geometry.Angles.end(); ++ang) {
    sdgDistance dist(*(ang->getAtom1()), *(ang->atom3), ang->ideal_angle, false);
    dists.push_back(dist);
  }
  //	for(atm=_res->atoms().begin(); atm!=_res->atoms().end(); ++atm) {
  //		if (atm->hybrid == 2 && atm->nabors.size() >= 3) {
  //			sdgVolume vol(*atm, *(atm->nabors[0]), *(atm->nabors[1]), *(atm->nabors[2]), 0.0);
  //			vols.push_back(vol);
  //		}
  //	}

  for (tor = _lig->geometry.Impropers.begin(); tor != _lig->geometry.Impropers.end(); ++tor) {
    double d = tor->ConvertToDistance(_lig->geometry);
    sdgDistance dist(*(tor->GetAtom(0)), *(tor->GetAtom(3)), d+0.3F, d-0.3, true);
    if (find(dists.begin(), dists.end(), dist) == dists.end()) {
      dists.push_back(dist);
    }
  }

  for (chir = _lig->geometry.Chirals.begin(); chir != _lig->geometry.Chirals.end(); ++chir) {
    GetTargetVolume(*chir, dists);
    if (chir->GetOrder() == COUNTERCLOCKWISE) {
      sdgVolume vol(*(chir->GetCenter()), *(chir->GetSub(0)), *(chir->GetSub(1)),
                    *(chir->GetSub(2)), 0.0, sdg_params::SDG_OPENEND_HIGH);
      vols.push_back(vol);
    }
    if (chir->GetOrder() == CLOCKWISE) {
      sdgVolume vol(*(chir->GetCenter()), *(chir->GetSub(0)), *(chir->GetSub(1)),
                    *(chir->GetSub(2)), 0.0, sdg_params::SDG_OPENEND_LOW);
      vols.push_back(vol);
    }
  }

  sdgEngine engine(dists, vols, _res->atoms());
  engine.DoOptimize();

}

double CoordGenerator::sdgRefine() {

  std::vector<sdgDistance> dists;
  std::vector<sdgVolume> vols;
  MIAtom_const_iter atm, atm2;
  std::vector<BondLength>::iterator bnd;
  std::vector<Angle>::iterator ang;
  //	std::vector<Bump>::iterator bmp;
  std::vector<Improper>::iterator tor;
  std::vector<Chiral>::iterator chir;

  for (bnd = _lig->geometry.Bonds.begin(); bnd != _lig->geometry.Bonds.end(); ++bnd) {
    sdgDistance dist(*(bnd->getAtom1()), *(bnd->getAtom2()), bnd->ideal_dist-0.002, bnd->ideal_dist+0.002, true);
    dists.push_back(dist);
  }
  for (ang = _lig->geometry.Angles.begin(); ang != _lig->geometry.Angles.end(); ++ang) {
    sdgDistance dist(*(ang->getAtom1()), *(ang->atom3), ang->ideal_angle-0.03, ang->ideal_angle+0.03, true);
    dists.push_back(dist);
  }
  /*
      for(bnd=_lig->geometry.Bonds.begin(); bnd!=_lig->geometry.Bonds.end(); ++bnd) {
          sdgDistance dist(*(bnd->getAtom1()), *(bnd->getAtom2()), bnd->ideal_dist - bnd->tolerance,
              bnd->ideal_dist + bnd->tolerance, true);
          dists.push_back(dist);
      }
      for(ang=_lig->geometry.Angles.begin(); ang!=_lig->geometry.Angles.end(); ++ang) {
          sdgDistance dist(*(ang->getAtom1()), *(ang->atom3), ang->ideal_angle - ang->tolerance,
              ang->ideal_angle + ang->tolerance, true);
          dists.push_back(dist);
      }

      for(bmp=_lig->geometry.Bumps.begin(); bmp!=_lig->geometry.Bumps.end(); ++bmp) {
          sdgDistance dist(*(bmp->getAtom1()), *(bmp->getAtom2()), bmp->min_d, true);
          if (find(dists.begin(), dists.end(), dist) == dists.end()) {
              dists.push_back(dist);
          }
      }

   */
  for (atm = _res->atoms().begin(); atm != _res->atoms().end(); ++atm) {
    if ((*atm)->hybrid() == 2 && (*atm)->nabors().size() > 2) {
      sdgVolume vol(**atm, *(*atm)->nabors()[0], *(*atm)->nabors()[1], *(*atm)->nabors()[2],
                    -0.01, 0.01, sdg_params::SDG_OPENEND_NONE);
      vols.push_back(vol);
    }
  }

  for (atm = _res->atoms().begin(); atm != _res->atoms().end(); ++atm) {
    for (atm2 = atm+1; atm2 != _res->atoms().end(); ++atm2) {
      sdgDistance dist(**atm, **atm2, 2.75, true);
      if (find(dists.begin(), dists.end(), dist) == dists.end()) {
        dists.push_back(dist);
      }
    }
  }
  /*

      for(tor=_lig->geometry.Impropers.begin(); tor!=_lig->geometry.Impropers.end(); ++tor) {
          double d = tor->ConvertToDistance(_lig->geometry);
          sdgDistance dist(*(tor->GetAtom(0)), *(tor->GetAtom(3)), d-0.2, d+0.2, true);
          if (find(dists.begin(), dists.end(), dist) == dists.end()) {
              dists.push_back(dist);
          }
      }
   */
  for (tor = _lig->geometry.Impropers.begin(); tor != _lig->geometry.Impropers.end(); ++tor) {
    if (tor->GetAngle(0) != 0.0 && tor->GetAngle(0) != 180.0) {
      continue;
    }

    sdgVolume vol(*(tor->GetAtom(2)), *(tor->GetAtom(0)),
                  *(tor->GetAtom(1)), *(tor->GetAtom(3)), -0.01, 0.01, sdg_params::SDG_OPENEND_NONE);
    if (find(vols.begin(), vols.end(), vol) == vols.end() ) {
      vols.push_back(vol);
    }
  }

  double target_vol;
  for (chir = _lig->geometry.Chirals.begin(); chir != _lig->geometry.Chirals.end(); ++chir) {
    target_vol = GetTargetVolume(*chir, dists);
    if (chir->GetOrder() == CLOCKWISE) {
      sdgVolume vol(*(chir->GetCenter()), *(chir->GetSub(0)), *(chir->GetSub(1)),
                    *(chir->GetSub(2)), 0.5 * target_vol, 2.0 * target_vol, sdg_params::SDG_OPENEND_NONE);
      vols.push_back(vol);
    }
    if (chir->GetOrder() == COUNTERCLOCKWISE) {
      sdgVolume vol(*(chir->GetCenter()), *(chir->GetSub(0)), *(chir->GetSub(1)),
                    *(chir->GetSub(2)), -2.0 * target_vol, -0.5 * target_vol, sdg_params::SDG_OPENEND_NONE);
      vols.push_back(vol);
    }
  }

  sdgEngine engine(dists, vols, _res->atoms());
  //	for(int i=0; i<10, ++i) {
  //		srand(87531);
  pair<double, double> interval(-20, 20);
  double score = 1.0;
  int nTrials = 0;
  while (score > 0.2 && nTrials < 50) {
    _lig->RandomizeCoords(interval, interval, interval);
    score = engine.DoOptimize();
    nTrials++;
  }
  //	}
  return score;
  //	LigRefiner lr;
  //	lr.SetModel(_lig);
  //	lr.AddRefiRes(_res);
  //	return lr.Refine();			//Now optimize all the atom positions at once
}

double CoordGenerator::GenMolStruct(std::string& log) {

  MIAtom* atom;
  MIAtom_const_iter nabor;

  _lig->ResetSearchFlags();

  for (unsigned int i = 0; i < _res->atoms().size(); ++i) {
    atom = _res->atom(i);

    if (atom->search_flag() == 0) {
      PlaceAtom(atom);
      char buf[1024];
      sprintf(buf, "Placed atom %s at %f %f %f\n", atom->name(), atom->x(), atom->y(), atom->z());
      log += buf;
    }



    if (atom->iscyclic()) {
      for (nabor = atom->nabors().begin(); nabor != atom->nabors().end(); ++nabor) {
        if ((*nabor)->iscyclic()) {                               //Place all cyclic neighbors now
          if ((*nabor)->search_flag() == 0) {

            PlaceAtom(*nabor);
            char buf[1024];
            sprintf(buf, "Placed atom %s at %f %f %f\n", (*nabor)->name(), (*nabor)->x(), (*nabor)->y(), (*nabor)->z());
            log += buf;
          }
        }
      }
    }
  }

  //	sdgGenMolStruct();
  //	sdgRefine();
  LigRefiner lr;
  lr.AddRefiRes(_res);
  lr.SetModel(_lig);
  double score = lr.Refine();           //Now optimize all the atom positions at once

  _lig->ResetSearchFlags();
  return score;
}

void CoordGenerator::SumRootVectors() {
  MIAtom_const_iter nabor;
  _newbond[0] = _newbond[1] = _newbond[2] = 0.0F;
  double norm;
  int nvects = 0;
  std::vector< std::vector<double> > bond_vectors;
  std::vector<double> zero;
  zero.push_back(0.0);
  zero.push_back(0.0);
  zero.push_back(0.0);

  for (nabor = _root->nabors().begin(); nabor != _root->nabors().end(); ++nabor) {
    if ((*nabor)->search_flag() != 0) {
      bond_vectors.push_back(zero);
      BondVector(_root, (*nabor), bond_vectors[nvects]);
      NormVect(bond_vectors[nvects]);
      _newbond[0] -= (float)bond_vectors[nvects][0];
      _newbond[1] -= (float)bond_vectors[nvects][1];
      _newbond[2] -= (float)bond_vectors[nvects][2];
      nvects++;
      //			_newbond[0] -= (*nabor)->x - _root->x;
      //			_newbond[1] -= (*nabor)->y - _root->y;
      //			_newbond[2] -= (*nabor)->z - _root->z;
    }
  }

  norm = sqrt(_newbond[0] * _newbond[0] +
           _newbond[1] * _newbond[1] +
           _newbond[2] * _newbond[2]);

  for (int i = 0; i <= 2; ++i) {
    _newbond[i] /= (float)norm;
  }
}

void CoordGenerator::GetRootAtom(MIAtom* atom) {
  MIAtom_const_iter nabor;
  for (nabor = atom->nabors().begin(); nabor != atom->nabors().end(); ++nabor) {
    if ((*nabor)->search_flag() != 0) {
      _root = *nabor;
      return;
    }
  }
  _root = 0;                                //Reach here only if there are no neihgbors placed
}

int CoordGenerator::GetRootNeighbor() {
  if (_root == 0) {
    _rnabor = 0;
    return 0;
  }
  int count = 0;
  MIAtom_const_iter nabor;
  for (nabor = _root->nabors().begin(); nabor != _root->nabors().end(); ++nabor) {
    if ((*nabor)->search_flag() != 0) {
      _rnabor = *nabor;
      count++;
    }
  }
  if (count == 0) {
    _rnabor = 0;
  }

  return count;
}

//Gets the std::vector from the "RootNabor" atom to the root atom
void CoordGenerator::GetRootVector() {
  if (_root != 0 && _rnabor != 0) {
    _rvect[0] = _root->x() - _rnabor->x();
    _rvect[1] = _root->y() - _rnabor->y();
    _rvect[2] = _root->z() - _rnabor->z();
  } else {
    _rvect[0] = 1; _rvect[1] = 0; _rvect[2] = 0;
  }
}

void CoordGenerator::PlaceAtom(MIAtom* atom) {

  if (atom->search_flag() != 0) {                 //Skip if this atom has been
    return;                                         //placed already!
  }

  GatherAtomConstraints(atom);

  CalcInitialCoordinates(atom);


  if (_current.Angles.size() == 0                   //Exit if there are no
      && _current.Bonds.size() == 0                 //constraints
      && _current.Planes.size() == 0
      && _current.Impropers.size() == 0) {
    atom->set_search_flag(1);
    return;
  }

  AtomRefiner ar(this);                             //Instantiate a refiner//
  ar.SetConstraintList(&_current);                  //Load the geometries gathered for this atom//
  ar.RefineAtom(atom, _root, _rnabor);              //Refinement//
  atom->set_search_flag(1);                        //Record that this atom has been placed//
}

void CoordGenerator::CalcInitialCoordinates(MIAtom* atom) {
  if (CountPlaced() == 0) {
    atom->setPosition(0.0f, 0.0f, 0.0f);
    _root = 0;
    _rnabor = 0;
    return;
  }

  //Get a starting point for this atom (i.e. the previous atom)
  GetRootAtom(atom);

  if (_root == 0) {
    atom->setPosition(20.0f, 20.0f, 20.0f);
    return;
  }

  double cosba = 0.0, sinba = 0.0;

  if (_root->geom() == TRIGONAL_PLANAR) {
    cosba =  0.5;
    sinba = 0.866;                                  //Magic numbers come from
  }                                                 //cosines and sines of the theoretically
  else if (_root->geom() == TETRAHEDRAL) {        //ideal bond angles for
    cosba = 0.3338;                                 //atoms of each geometry.
    sinba = 0.9246;
  } else if (_root->geom() == LINEAR) {
    cosba = 1;
    sinba = 0;
  }

  float length1_2 = _lig->GetBond(atom, _root)->ideal_length;

  if (CountPlaced() == 1) {
    atom->setPosition(length1_2, 0.0f, 0.0f);
    _rnabor = 0;
    return;
  }

  int n_nabors = GetRootNeighbor();

  if (_rnabor == 0) {
    atom->setPosition(_root->x() + length1_2, _root->y(), _root->z());
    return;
  }

  float length2_3 = _lig->GetBond(_root, _rnabor)->ideal_length;
  float length1_3;

  if (CountPlaced() == 2 && _current.Bonds.size() == 1 && _current.Angles.size() == 1) {
    length1_3 = _current.Angles.front().ideal_angle;
    double angle = Nonstringent_acos((length1_2 * length1_2 +
                                      length2_3 * length2_3 -
                                      length1_3 * length1_3) / (2 * length1_2 * length2_3));

    PositionAtom(atom, _root, _rnabor, length1_2, (float)(RAD2DEG * angle));
    return;
  }


  Bond* bond = _lig->GetBond(atom, _root);

  double temp[3];
  atom->setPosition(_root->x(),                       //Start new atom on top of the previous atom
      _root->y(),
      _root->z());

  if (n_nabors == 1) {
    GetRootVector();

    NormVect(_rvect, _rvect);
    if (cosba != 0) {
      atom->translate((float)(_rvect[0] * cosba * bond->ideal_length),       //Move to an arbitrary point
          (float)(_rvect[1] * cosba * bond->ideal_length),       //with the correct bond length
          (float)(_rvect[2] * cosba * bond->ideal_length));       //and angle.
    }

    temp[0] = -_rvect[1];
    temp[1] = _rvect[0];
    temp[2] = 0;
    NormVect(temp, temp);
    if (sinba != 0) {
      atom->translate((float)(temp[0] * sinba * bond->ideal_length),
          (float)(temp[1] * sinba * bond->ideal_length),
          0.0f);  //Have chosen delta(z)=0
    }
  } else if (n_nabors > 1) {
    SumRootVectors();
    AtomStep(atom, _newbond, bond->ideal_length);
  }

}

void CoordGenerator::SetExtended(MIAtom* atom) {
  if (atom->iscyclic()) {
    return;
  }

  MIAtom* atom1;
  bool found = false;
  MIAtom_const_iter atm;
  for (atm = _rnabor->nabors().begin(); atm != _rnabor->nabors().end(); ++atm) {
    if (*atm == _root) {                                            //We need a different atom
      continue;                                                     //than _root
    }
    if ((*atm)->search_flag() != 0) {                     //Look until we find an atom
      atom1 = *atm;                                         //that has been placed
      found = true;
      break;
    }
  }

  if (found == false) {                                             //Can't set the torsion
    return;                                                         //if we don't have a 4th atom
  }                                                                 //placed

  float chi;
  chi = (float)CalcAtomTorsion(atom1, _rnabor, _root, atom);        //Get the current torsion value

  RotateAtom(_rnabor, _root, atom, 180-chi);
}

void CoordGenerator::GatherAtomConstraints(MIAtom* atom) {
  std::vector<BondLength>::iterator bnd;
  std::vector<Angle>::iterator ang;
  std::vector<Plane>::iterator pln;
  std::vector<Improper>::iterator imp;
  std::vector<Bump>::iterator bmp;
  std::vector<Chiral>::iterator chrl;


  _current.Clear();                         //There's also a delete[] call in this ftn

  for (bnd = _lig->geometry.Bonds.begin(); bnd != _lig->geometry.Bonds.end(); ++bnd) {

    if (atom == bnd->getAtom1() && bnd->getAtom2()->search_flag() != 0) {
      _current.Bonds.push_back(*bnd);
    }
    if (atom == bnd->getAtom2() && bnd->getAtom1()->search_flag() != 0) {
      _current.Bonds.push_back(*bnd);
    }
  }

  for (ang = _lig->geometry.Angles.begin(); ang != _lig->geometry.Angles.end(); ++ang) {
    if (atom == ang->getAtom1() && ang->atom3->search_flag() != 0) {
      _current.Angles.push_back(*ang);
    }
    if (atom == ang->atom3 && ang->getAtom1()->search_flag() != 0) {
      _current.Angles.push_back(*ang);
    }
  }

  bool in_plane;
  int i, n = 0;
  std::vector<Plane> p(100);   //A partial plane containing only atoms placed so far
  for (pln = _lig->geometry.Planes.begin(); pln != _lig->geometry.Planes.end(); ++pln) {
    in_plane = false;
    for (i = 0; i < pln->NumAtoms(); ++i) {
      if (pln->GetAtom(i) == atom) {
        in_plane = true;
      }
    }

    if (!in_plane) {
      continue;
    }

    if (CountPlaced(&*pln) >= 3) {
      //			ClearPlane(&p);
      PartialPlane(&*pln, &p[n]);
      //			_current.AddPlane(p[n]);
      _current.Planes.push_back(p[n]);
      n++;
    }
  }

  for (imp = _lig->geometry.Impropers.begin(); imp != _lig->geometry.Impropers.end(); ++imp) {
    if (atom == imp->GetAtom(0)
        && imp->GetAtom(1)->search_flag() != 0
        && imp->GetAtom(2)->search_flag() != 0
        && imp->GetAtom(3)->search_flag() != 0) {
      _current.Impropers.push_back(*imp);
    }
    if (atom == imp->GetAtom(3)
        && imp->GetAtom(0)->search_flag() != 0
        && imp->GetAtom(1)->search_flag() != 0
        && imp->GetAtom(2)->search_flag() != 0) {
      _current.Impropers.push_back(*imp);
    }
  }

  for (chrl = _lig->geometry.Chirals.begin(); chrl != _lig->geometry.Chirals.end(); ++chrl) {
    if (atom == chrl->GetSub(0)
        && chrl->GetCenter()->search_flag() != 0
        && chrl->GetSub(1)->search_flag() != 0
        && chrl->GetSub(2)->search_flag() != 0) {
      _current.Chirals.push_back(*chrl);
    }
    if (atom == chrl->GetSub(1)
        && chrl->GetCenter()->search_flag() != 0
        && chrl->GetSub(0)->search_flag() != 0
        && chrl->GetSub(2)->search_flag() != 0) {
      _current.Chirals.push_back(*chrl);
    }
    if (atom == chrl->GetSub(2)
        && chrl->GetCenter()->search_flag() != 0
        && chrl->GetSub(0)->search_flag() != 0
        && chrl->GetSub(1)->search_flag() != 0) {
      _current.Chirals.push_back(*chrl);
    }
  }

  for (bmp = _lig->geometry.Bumps.begin(); bmp != _lig->geometry.Bumps.end(); ++bmp) {

    if (atom == bmp->getAtom1() && bmp->getAtom2()->search_flag() != 0) {
      _current.Bumps.push_back(*bmp);
    }
    if (atom == bmp->getAtom2() && bmp->getAtom1()->search_flag() != 0) {
      _current.Bumps.push_back(*bmp);
    }
  }

}

int CoordGenerator::CountPlaced() {
  unsigned int count = 0;
  for (unsigned int i = 0; i < _res->atoms().size(); ++i) {
    if (_res->atom(i)->search_flag() != 0) {
      ++count;
    }
  }
  return count;
}

int CoordGenerator::CountPlaced(Plane* plane) {
  int count = 0;
  for (int i = 0; i < plane->NumAtoms(); ++i) {
    if (plane->GetAtom(i)->search_flag() != 0) {
      ++count;
    }
  }
  return count;
}

void CoordGenerator::PartialPlane(Plane* source, Plane* target) {

  for (int i = 0; i < source->NumAtoms(); ++i) {
    if (source->GetAtom(i)->search_flag() != 0) {
      target->AddAtom(source->GetAtom(i));
    }
  }
  if (source->GetTolerance() != 0) {
    target->SetTolerance(source->GetTolerance());
  }
}

namespace conflib {
double GetTargetVolume(Chiral& chiral, std::vector<sdgDistance>& dists) {

  vector<double> ideal_dists;

  ideal_dists.push_back(FindDistance(chiral.GetCenter(), chiral.GetSub(0), dists).GetIdeal() );
  ideal_dists.push_back(FindDistance(chiral.GetCenter(), chiral.GetSub(1), dists).GetIdeal() );
  ideal_dists.push_back(FindDistance(chiral.GetCenter(), chiral.GetSub(2), dists).GetIdeal() );
  ideal_dists.push_back(FindDistance(chiral.GetSub(0), chiral.GetSub(1), dists).GetIdeal() );
  ideal_dists.push_back(FindDistance(chiral.GetSub(0), chiral.GetSub(2), dists).GetIdeal() );
  ideal_dists.push_back(FindDistance(chiral.GetSub(1), chiral.GetSub(2), dists).GetIdeal() );

  return CalcTetraVolume(ideal_dists);
}

sdgDistance FindDistance(MIAtom* a1, MIAtom* a2, std::vector<sdgDistance>& dists) {

  sdgDistance d(*a1, *a2, 2.5, false);
  std::vector<sdgDistance>::iterator i = find(dists.begin(), dists.end(), d);

  if (i == dists.end()) {
    return d;
  } else {
    return *i;
  }
}

} //namespace conflib
