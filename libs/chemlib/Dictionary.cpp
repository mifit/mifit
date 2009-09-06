#include "Dictionary.h"
#include "RingSystem.h"
#include "mol_util.h"

#include <vector>

using namespace std;

namespace chemlib {

LigDictionary::LigDictionary(Ligand* mol, const Residue* res) {
  lig = mol;                                                //Store pointers to the
  _res = res;                                               //molecule and residue to process
  _sigmaplane = 0.03F;
}

void LigDictionary::AddPlane(Plane& plane) {
  plane.SetResidue(_res);
  if (find(planes.begin(), planes.end(), plane) == planes.end()) {
    planes.push_back(plane);
  }
}

void LigDictionary::AddTorsion(Torsion& torsion) {
  torsion.SetResidue(_res);
  if (find(torsions.begin(), torsions.end(), torsion) == torsions.end()) {
    torsions.push_back(torsion);
  }
}

/////////////////////////////////////////////////////////////////////////////
// Function:    AddImproper
// Purpose:		Adds a new Improper constrain to the dictionary and assigns it
//				a number.
// Input:       A reference to an Improper object and an (optional) flag as to
//				whether to include it in the dictionary for ligand fitting
//				(Default value of the flag is true)
// Output:      Adds impropers to the dictionary object
// Requires:
/////////////////////////////////////////////////////////////////////////////
void LigDictionary::AddImproper(Improper& imp, bool report) {
  imp.SetResidue(_res);

  if (find(impropers.begin(), impropers.end(), imp) == impropers.end()) {
    if (report) {
      imp.IncludeInDict();
    } else {
      imp.ExcludeFromDict();
    }

    imp.SetIndex(impropers.size() + 1);
    impropers.push_back(imp);
  }
}

/////////////////////////////////////////////////////////////////////////////
// Function:    UnspecDoubleBond
// Input:       A pointer to a double bond
// Output:      Adds impropers to the dictionary with two acceptable angles,
//				0 and 180 degrees
// Requires:
/////////////////////////////////////////////////////////////////////////////
void LigDictionary::UnspecDoubleBond(Bond* bond) {

  std::vector<double> angles;                   //Initialize a vector with the ideal angles
  angles.push_back(0);
  angles.push_back(180);

  std::vector< MIAtomList > torsions;    //Create vectors with the 4-atom sequences
  EnumerateTorsions(bond, torsions);            //of all the torsions around this bond

  std::vector< MIAtomList >::iterator atms;
  Improper imp;

  for (atms = torsions.begin(); atms != torsions.end(); ++atms) {
    imp.ReInit(*atms, angles);
  }

  if (bond->iscyclic) {
    AddImproper(imp, true);                     //Report this improper in the dictionary only
  }                                             //if it is in a ring...otherwise its only for
  else {                                        //internal use in coordinate generation
    AddImproper(imp, false);
  }
}

/////////////////////////////////////////////////////////////////////////////
// Function:    AnalyzeRings
// Purpose:		Generates Plane and Improper entries associated with rings
// Input:       Data lookups from parent molecule
// Output:      Stores entries in the internal data of the LigDictionary object
// Requires:	FindRingSystems() has been called for the parent molecule
/////////////////////////////////////////////////////////////////////////////

void LigDictionary::AnalyzeRings() {
  std::vector<RingSystem>::iterator rs;
  for (rs = lig->ringsystems.begin(); rs != lig->ringsystems.end(); ++rs) {
    rs->GeneratePlanes(*this);
    rs->GenerateImpropers(*this);
    rs->CyclohexImpropers(*this);
  }
}

/*
   //Infer the required angle for a ring improper from whether the end atoms are in
   //the ring.  (The two middle atoms must be in the ring.)  If both are in or both
   //are out, the angle shoud be zero.  If only atom is outside the ring, angle=180.
   double LigDictionary::ImproperAngle(MIAtom *atom1, MIAtom *atom4, Cycle &rng) {
    bool in_ring1 = rng.Contains(atom1);
    bool in_ring4 = rng.Contains(atom4);

    if (in_ring1 && in_ring4) {
        return 0.0;							//All atoms are in the ring
    }
    else if (in_ring1 || in_ring4) {
        return 180.0;						//Only one atom is outside the ring
    }
    else {
        return 0.0;							//Both end atoms are out of the ring
    }
   }
 */
void LigDictionary::FindAcycPlanes() {
  MIAtom* ctr;
  vector<MIAtom*>::const_iterator atm;
  Plane pln(_sigmaplane);

  unsigned int i;

  for (i = 0; i < _res->atoms().size(); ++i) {            //Loop over all atoms
    ctr = _res->atom(i);

    if (ctr->iscyclic()) {                        //Skip atoms in rings
      continue;
    }

    if (ctr->hybrid() == 2
        && ctr->nabors().size() == 3) {
      pln.SetResidue(_res);
      pln.AddAtom(ctr);                         //Add central atom to plane
      for (atm = ctr->nabors().begin();           //Find neighboring atoms
           atm != ctr->nabors().end();
           ++atm) {
        pln.AddAtom(*atm);                      //Add neighbors to the plane
      }
      AddPlane(pln);
      pln.Clear();
    }
  }   //Loop atoms (i)
}

void LigDictionary::GenDoubleBondImpropers() {
  vector<Bond*> bonds;
  vector<Bond*>::iterator bnd;

  lig->GetResidueBonds(_res, bonds);            //Assume we just want a single residue

  MIAtom* atom2;
  MIAtom* atom3;
  Bond* bnd1;
  Bond* bnd2;
  vector<MIAtom*>::const_iterator atm1, atm4;
  Improper imp;


  for (bnd = bonds.begin(); bnd != bonds.end(); ++bnd) {
    if ((*bnd)->getOrder() != 2) {                      //Treat only double bonds
      continue;
    }

    if ((*bnd)->getAtom1()->geom() != TRIGONAL_PLANAR         //Probably unnecessary, but check
        || (*bnd)->getAtom2()->geom() != TRIGONAL_PLANAR) {   //that both geometries are trigonal
      continue;
    }

    if ((*bnd)->getAtom1()->nabors().size() == 1              //Skip bonds with terminal atoms
        || (*bnd)->getAtom2()->nabors().size() == 1) {
      continue;
    }

    atom2 = (*bnd)->getAtom1();
    atom3 = (*bnd)->getAtom2();

    //Check flags for consistency of bond directions (e.g. check that the atom
    //doesn't have two atoms in the "up" direction)

    if (atom2->nabors().size() >= 3 && !CheckConsistency(atom2->bondnumbers())) {
      continue;
    }

    if (atom3->nabors().size() >= 3 && !CheckConsistency(atom3->bondnumbers())) {
      continue;
    }


    for (atm1 = atom2->nabors().begin(); atm1 != atom2->nabors().end(); ++atm1) {
      if (*atm1 == atom3) {
        continue;
      }
      for (atm4 = atom3->nabors().begin(); atm4 != atom3->nabors().end(); ++atm4) {
        if (*atm4 == atom2) {
          continue;
        }

        bnd1 = lig->GetBond(*atm1, atom2);
        bnd2 = lig->GetBond(*atm4, atom3);

        if (*atm1 == *atm4) {
          continue;                                         //fix for cyclopropenes

        }
        if ((bnd1->stereo == STEREO_UP
            && bnd2->stereo == STEREO_UP)
            || (bnd1->stereo == STEREO_DOWN
            && bnd2->stereo == STEREO_DOWN)) {
          imp.AddAtom(*atm1);
          imp.AddAtom(atom2);
          imp.AddAtom(atom3);
          imp.AddAtom(*atm4);
          imp.AddAngle(180.0F);
          imp.SetIndex(impropers.size()+1);
          AddImproper(imp, false);
          imp.Clear();
        } else if ((bnd1->stereo == STEREO_DOWN
                   && bnd2->stereo == STEREO_UP)
                   || (bnd1->stereo == STEREO_UP
                   && bnd2->stereo == STEREO_DOWN)) {
          imp.AddAtom(*atm1);
          imp.AddAtom(atom2);
          imp.AddAtom(atom3);
          imp.AddAtom(*atm4);
          imp.AddAngle(0.0F);
          imp.SetIndex(impropers.size()+1);
          AddImproper(imp, false);
          imp.Clear();
        } else if (bnd1->stereo == 0                            //Double bond with unspecified
                   || bnd2->stereo == 0) {                      //stereochemistry
          imp.AddAtom(*atm1);
          imp.AddAtom(atom2);
          imp.AddAtom(atom3);
          imp.AddAtom(*atm4);
          imp.AddAngle(0.0F);
          imp.AddAngle(180.0F);
          imp.SetIndex(impropers.size()+1);
          AddImproper(imp, false);
          imp.Clear();
        }
      }
    }
  }
}

//Generate torsions for all acyclic bonds between non-terminal atoms
void LigDictionary::GenTorsions() {
  vector<Bond*> bonds;
  vector<Bond*>::iterator bnd;
  Torsion tors;
  MIAtom* atom1, * atom2, * atom3, * atom4;
  lig->GetResidueBonds(_res, bonds);            //Assume we just want a single residue

  for (bnd = bonds.begin(); bnd != bonds.end(); ++bnd) {
    if ((*bnd)->getAtom1()->nabors().size() == 1              //Skip bonds with terminal atoms
        || (*bnd)->getAtom2()->nabors().size() == 1) {
      continue;
    }

    if ((*bnd)->iscyclic) {                                     //Skip bonds in rings
      continue;
    }

    atom2 = lig->MoreCentralAtm((*bnd)->getAtom1(), (*bnd)->getAtom2());

    if (atom2 == (*bnd)->getAtom1()) {
      atom3 = (*bnd)->getAtom2();
    } else {
      atom3 = (*bnd)->getAtom1();
    }


    atom1 = lig->GetNewNabor(atom2, atom3);                 //Pick arbitrary neighbor atoms
    atom4 = lig->GetNewNabor(atom3, atom2);                 //for the ends

    tors.AddAtom(atom1);
    tors.AddAtom(atom2);
    tors.AddAtom(atom3);
    tors.AddAtom(atom4);
    tors.SetIndex(torsions.size()+1);
    AddTorsion(tors);
    //		torsions.push_back(tors);
    tors.Clear();
  }

}

bool LigDictionary::CheckConsistency(const std::vector<int>& bond_indices) {
  int n_up = 0;
  int n_down = 0;
  std::vector<int>::const_iterator x_bond;

  for (x_bond = bond_indices.begin(); x_bond != bond_indices.end(); ++x_bond) {
    if ((lig->GetBond(*x_bond))->stereo == STEREO_UP) {
      ++n_up;
    } else if ((lig->GetBond(*x_bond))->stereo == STEREO_DOWN) {
      ++n_down;
    }
  }

  if (n_up > 1) {
    return false;
  } else if (n_down > 1) {
    return false;
  } else {
    return true;
  }
}

void LigDictionary::GenChirals() {
  const MIAtom* ctr;
  unsigned int i;

  for (i = 0; i < _res->atoms().size(); ++i) {                //Loop over all atoms
    ctr = _res->atom(i);

    if (ctr->chiral_class() == CH_NONE) {
      continue;
    }

    if (ctr->chiral_class() == CH_TETRAHEDRAL && ctr->bondnumbers().size() == 4) {

      AddChiral(ctr, 0, 1, 2, ctr->chiral_order());
      AddChiral(ctr, 0, 2, 3, ctr->chiral_order());
      AddChiral(ctr, 0, 3, 1, ctr->chiral_order());
      AddChiral(ctr, 3, 2, 1, ctr->chiral_order());
    } else if (ctr->chiral_class() == CH_TETRAHEDRAL && ctr->bondnumbers().size() == 3) {
      AddChiral(ctr, 0, 1, 2, ctr->chiral_order());
    }
  }
}

void LigDictionary::AddChiral(const MIAtom* ctr, int n1, int n2, int n3, int order) {
  Chiral chiral;

  chiral.SetCenter(const_cast<MIAtom*> (ctr));
  chiral.AddSub(lig->GetNabor(ctr, n1));
  chiral.AddSub(lig->GetNabor(ctr, n2));
  chiral.AddSub(lig->GetNabor(ctr, n3));
  chiral.SetOrder(order);

  chirals.push_back(chiral);
}

/*
   void LigDictionary::Write(ofstream &test_out) {
    vector<Torsion>::iterator tor;
    vector<Plane>::iterator pln;
    vector<Improper>::iterator imp;

    std::string s;
   //	ofstream test_out("test3.pdb");

    for (tor=torsions.begin();tor!=torsions.end();++tor) {
        s.assign("TORSION ");
        s.append(_res->type());
        s.append(" CHI");
        tor->Card(s);
        test_out << s.c_str() << endl;
    }

    for (pln=planes.begin();pln!=planes.end();++pln) {
        s.assign("PLANE ");
        s.append(_res->type());
        pln->Card(s);
        test_out << s.c_str() << endl;
    }

    for (imp=impropers.begin();imp!=impropers.end();++imp) {
        if (imp->DictTest() == false) {
            continue;
        }
        s.assign("TORSION ");
        s.append(_res->type());
        s.append(" IMP");
        imp->Card(s);
        test_out << s.c_str() << endl;
    }
   }
 */
/*
   void LigDictionary::Print() {
    vector<Torsion>::iterator tor;
    vector<Plane>::iterator pln;
    vector<Improper>::iterator imp;
    vector<Chiral>::iterator chrl;

    std::cout << "Printing Torsions... ";
    if (torsions.size() == 0) {
        std::cout << "\t" << "no torsions." << endl;
    }
    else {
        std::cout << torsions.size() << " torsions:" << endl;
    }
    std::string s;
    for (tor=torsions.begin();tor!=torsions.end();++tor) {
        s.assign("TORSION ");
        s.append(_res->type());
        s.append(" CHI");
        tor->Card(s);
        std::cout << s.c_str() << endl;
    }

    std::cout << "Printing Planes... ";
    if (planes.size() == 0) {
        std::cout << "\t" << "no planes." << endl;
    }
    else {
        std::cout << planes.size() << " planes:" << endl;
    }
    for (pln=planes.begin();pln!=planes.end();++pln) {
        s.assign("PLANE ");
        s.append(_res->type());
        pln->Card(s);
        std::cout << s.c_str() << endl;
    }

    std::cout << "Printing Improper List... ";
    if (impropers.size() == 0) {
        std::cout << "\t" << "no impropers." << endl;
    }
    else {
        std::cout << impropers.size() << " impropers:" << endl;
    }
    for (imp=impropers.begin();imp!=impropers.end();++imp) {
        s.assign("TORSION ");
        s.append(_res->type());
        s.append(" IMP");
        imp->Card(s);
        std::cout << s.c_str() << endl;
    }

    std::cout << "Printing Chiral List... ";
    if (chirals.size() == 0) {
        std::cout << "\t" << "no chirals." << endl;
    }
    else {
        std::cout << chirals.size() << " chirals:" << endl;
    }
    for (chrl=chirals.begin();chrl!=chirals.end();++chrl) {
        s.assign("CHIRAL ");
        s.append(_res->type());
        chrl->Card(s);
        std::cout << s.c_str() << endl;
    }
   }
 */
void LigDictionary::SetMolGeometry() {
  vector<Improper>::iterator imp;
  vector<Plane>::iterator pln;
  vector<Chiral>::iterator chrl;
  vector<Torsion>::iterator tors;

  lig->geometry.Impropers.clear();
  lig->geometry.Planes.clear();
  lig->geometry.Chirals.clear();
  lig->geometry.Torsions.clear();

  for (imp = impropers.begin(); imp != impropers.end(); ++imp) {        //Loop over impropers
    lig->geometry.AddImproper(*imp);
  }

  for (pln = planes.begin(); pln != planes.end(); ++pln) {
    lig->geometry.AddPlane(*pln);                           //Add the temp plane to the list
  }                                                         //of geometric features for the
                                                            //the molecule.
  for (chrl = chirals.begin(); chrl != chirals.end(); ++chrl) {
    lig->geometry.AddChiral(*chrl);
  }

  for (tors = torsions.begin(); tors != torsions.end(); ++tors) {
    lig->geometry.AddTorsion(*tors);
  }
}

}
