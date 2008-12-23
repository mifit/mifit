#include "FigRefiner.h"
#include "chemlib.h"

using namespace chemlib;
using namespace moldraw;

namespace moldraw {

void FigOpt(Ligand& site, const std::vector<HBond>& hbonds) {
  double fitness = CalcFigFitness(site, hbonds);
  double new_fitness;
  bool no_flips = false;
  int ncycles = 0;

  std::vector<Bond>::iterator bnd;
  while (no_flips == false && ncycles <= MAX_FLIP_CYCLES) {
    no_flips = true;
    ncycles++;

    for (bnd = site.bonds.begin(); bnd != site.bonds.end(); ++bnd) {
      if (!bnd->IsRotatable()) {
        continue;
      }

      site.FlipBond(*bnd);
      new_fitness = CalcFigFitness(site, hbonds);
      if (new_fitness < fitness) {
        fitness = new_fitness;
        no_flips = false;
      } else {
        site.FlipBond(*bnd);
      }
    }
  }
}

double CalcFigFitness(const Ligand& site, const std::vector<HBond>& hbonds) {
  return AtomClashTotal(site) + BOND_CLASH_WEIGHT* BondClashTotal(site, hbonds);
}

double AtomClashTotal(const Ligand& site) {
  double total = 0;
  unsigned int xRes1, xRes2;            //Residue indices
  unsigned int xAtom1, xAtom2;          //Atom indices

  for (xRes1 = 0; xRes1 < site.residues.size(); ++xRes1) {
    for (xAtom1 = 0; xAtom1 < site.residues[xRes1]->atoms().size(); ++xAtom1) {
      for (xRes2 = xRes1; xRes2 < site.residues.size(); ++xRes2) {
        for (xAtom2 = 0; xAtom2 < site.residues[xRes2]->atoms().size(); ++xAtom2) {

          if (xRes1 == xRes2 && xAtom1 >= xAtom2) {
            continue;
          }
          if (site.AlreadyBonded(site.residues[xRes1]->atom(xAtom1),
                site.residues[xRes2]->atom(xAtom2))) {
            continue;
          }
          total += AtomAtomClashScore(*site.residues[xRes1]->atom(xAtom1),
                     *site.residues[xRes2]->atom(xAtom2));
        }
      }
    }
  }

  return total;
}

double BondClashTotal(const Ligand& site, const std::vector<HBond>& hbonds) {
  double total = 0;

  //Step 1: Calculate the penalties for intersecting HBonds
  std::vector<HBond>::const_iterator hb1;
  std::vector<HBond>::const_iterator hb2;
  for (hb1 = hbonds.begin(); hb1 != hbonds.end(); ++hb1) {
    for (hb2 = hb1+1; hb2 != hbonds.end(); ++hb2) {
      if (hb1->donor == hb2->donor                      //Be sure that there
          || hb1->donor == hb2->acceptor                    //are four separate atoms here, not
          || hb1->acceptor == hb2->donor                    //two hbonds to the same atom
          || hb1->acceptor == hb2->acceptor) {
        continue;
      }

      total += BondBondClashScore(*hb1, *hb2);
    }
  }

  //Step 2: Calculate the scores for atoms clashing with HBonds
  std::vector<Residue*>::const_iterator ri;
  unsigned int xAtom;
  for (hb1 = hbonds.begin(); hb1 != hbonds.end(); ++hb1) {
    for (ri = site.residues.begin(); ri != site.residues.end(); ++ri) {
      Residue *res=*ri;
      for (xAtom = 0; xAtom < res->atoms().size(); ++xAtom) {
        if (hb1->donor == res->atom(xAtom)
            || hb1->acceptor == res->atom(xAtom)) {
          continue;
        }

        total += BondAtomClashScore(*hb1, *res->atom(xAtom));
      }
    }
  }

  return total;

}

double AtomAtomClashScore(const MIAtom& atom1, const MIAtom& atom2) {
  double sq_dist = (atom2.x() - atom1.x()) * (atom2.x() - atom1.x()) +
                   (atom2.y() - atom1.y()) * (atom2.y() - atom1.y()) +
                   (atom2.z() - atom1.z()) * (atom2.z() - atom1.z());

  if (sq_dist > SQ_ATOM_ATOM_THRESHOLD) {
    return 0;
  } else if (sq_dist > SMALL_DISTANCE) {
    return 1 / sq_dist;
  } else {
    return 1 / SMALL_DISTANCE;                  //Clamp to prevent infinite score
  }
}

double BondAtomClashScore(const HBond& bond, const MIAtom& atom) {

  //Calculate the bond vector
  double bond_vect[3];
  bond_vect[0] = bond.acceptor->x() - bond.donor->x();
  bond_vect[1] = bond.acceptor->y() - bond.donor->y();
  bond_vect[2] = bond.acceptor->z() - bond.donor->z();

  //Calculate the vector from bond->donor to the query atom
  double atom_vect[3];
  atom_vect[0] = atom.x() - bond.donor->x();
  atom_vect[1] = atom.y() - bond.donor->y();
  atom_vect[2] = atom.z() - bond.donor->z();

  //Calculate the the dot product (angle) between these two vectors.	If this angle
  //is greater than 90 degrees, do not consider further
  double dot = DotVect(bond_vect, atom_vect);

  if (dot < 0) {
    return 0;
  }

  //Now calculate the distance from a point to a line.  This is an application
  //of the Pythagorean theorem, where atom_vect is the hypotoneuse, a projection
  //onto the bond_vect is one leg, and the other leg is the desired distance
  double point_to_line = SquaredVectLength(atom_vect) - dot * dot / (SquaredVectLength(bond_vect));

  if (point_to_line > SQ_BOND_ATOM_THRESHOLD) {
    return 0;
  }

  //Now check if the perpendicular from the atom to the bond vector lies between
  //the two atoms of the bond.  If so, calculate the clash score
  if (SquaredVectLength(atom_vect) - point_to_line <= SquaredVectLength(bond_vect)) {
    return BondClashScore(point_to_line);
  } else {
    return 0;
  }
}

double BondBondClashScore(const HBond& bond1, const HBond& bond2) {
  double midpt_connector[3];            //A vector between the midpoints of the two bonds

  midpt_connector[0] = 0.5 *(bond2.donor->x() + bond2.acceptor->x() -
                             bond1.donor->x() - bond1.acceptor->x());
  midpt_connector[1] = 0.5 *(bond2.donor->y() + bond2.acceptor->y() -
                             (bond1.donor->y() + bond1.acceptor->y()));
  midpt_connector[2] = 0.5 *(bond2.donor->z() + bond2.acceptor->z() -
                             (bond1.donor->z() + bond1.acceptor->z()));

  if (SquaredVectLength(midpt_connector) < SQ_BOND_BOND_THRESHOLD) {
    return BOND_INTERSECT_PENALTY;
  } else {
    return 0;
  }
}

} //namespace moldraw
