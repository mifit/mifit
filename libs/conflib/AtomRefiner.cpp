#include <nongui/nonguilib.h>
#include <vector>
#include <fstream>
#include <sstream>

#include "AtomRefiner.h"
#include <chemlib/chemlib.h>
#include "CoordGenerator.h"
#include <cmath>

using namespace chemlib;
using namespace conflib;


AtomRefiner::AtomRefiner(CoordGenerator* xyzGen) {
  srand(22);
  //	_mol = mol;
  BondWeight = AngleWeight = PlaneWeight = TorsionWeight = 1;
  _nTwists = 360;                                               //Scan torsions every 3 degrees
  _nCycles = 5000;
  _xyzGen = xyzGen;
}

AtomRefiner::~AtomRefiner() {
}

void AtomRefiner::RefineAtom(MIAtom* atom, MIAtom* root, MIAtom* rnabor) {
  float step[3], step_size;
  double best, score, worst, orig;
  float best_coord[3], orig_coord[3];
  //	float rot_mat[4][3];
  double twist_size = 360.0 / _nTwists;

  best_coord[0] = orig_coord[0] = atom->x();
  best_coord[1] = orig_coord[1] = atom->y();
  best_coord[2] = orig_coord[2] = atom->z();
  best = ScoreAtom(atom);
  worst = best;

  int i;
  if (root != 0 && rnabor != 0) {
    for (i = 1; i < _nTwists; i++) {
      RotateAtom(rnabor, root, atom, (float)(twist_size * i));

      score = ScoreAtom(atom);

      if (score < best) {
        best = score;
        best_coord[0] = atom->x();
        best_coord[1] = atom->y();
        best_coord[2] = atom->z();
      } else if (score > worst) {
        worst = score;
      }

      atom->setPosition(orig_coord[0], orig_coord[1], orig_coord[2]);
    }

    if (worst - best < 0.001) {
      _xyzGen->SetExtended(atom);
      best = ScoreAtom(atom);
      best_coord[0] = atom->x();
      best_coord[1] = atom->y();
      best_coord[2] = atom->z();
    } else {
      atom->setPosition(best_coord[0], best_coord[1], best_coord[2]);
    }

    orig_coord[0] = best_coord[0];
    orig_coord[1] = best_coord[1];
    orig_coord[2] = best_coord[2];
  }

  //Store some stuff so that we can reset if the optimization made
  //no progress (helps the program be more consistent & reproducible)
  orig = best;
  for (i = 0; i < _nCycles; i++) {
    step[0] = -1 + 2 * (rand() / (float)RAND_MAX);
    step[1] = -1 + 2 * (rand() / (float)RAND_MAX);
    step[2] = -1 + 2 * (rand() / (float)RAND_MAX);
    step_size = 0.01f * (rand() / (float)RAND_MAX);

    AtomStep(atom, step, step_size);
    score = ScoreAtom(atom);

    if (score < best) {
      best = score;
      best_coord[0] = atom->x();
      best_coord[1] = atom->y();
      best_coord[2] = atom->z();
    } else {
      atom->setPosition(best_coord[0], best_coord[1], best_coord[2]);
    }

  }

  if (orig - best < 0.001) {
    atom->setPosition(orig_coord[0], orig_coord[1], orig_coord[2]);
  }
}

void AtomRefiner::Min1D(MIAtom* atom, float* grad, float& min) {
  float a[3];           //First ("lower" bound)--also the starting point
  float b[3];       //First of two search points within the bracket
  float c[3];       //Second search point
  float d[3];           //"Upper" bound
  float fb, fc;         //Scores of the inner points
  float dxb, dyb, dzb;      //Estimated derivatives (differentials) of the inner points
  float dxc, dyc, dzc;

  a[0] = atom->x();
  a[1] = atom->y();                           //Store the original position--this is one
  a[2] = atom->z();                           //end of our bracket.

  Bracket3DMin(atom, grad, min, b, d);      //Get points b and d, along vector grad,
                                            //where score(b) is less than
                                            //both score(d) and score(original position)
  c[0] = (b[0] + d[0]) / 2;
  c[1] = (b[1] + d[1]) / 2;                 //Create another point in the bracket [b,d]
  c[2] = (b[2] + d[2]) / 2;

  MoveAtom(atom, b);                        //Evaluate the score and differentials at point b
  fb = ScoreAtom(atom);
  dxb = atom->dx();
  dyb = atom->dy();
  dzb = atom->dz();

  MoveAtom(atom, c);                        //Evaluate the score and differential at point b
  fc = ScoreAtom(atom);
  dxc = atom->dx();
  dyc = atom->dy();
  dzc = atom->dz();

  while ((b[0]-c[0])*(b[0]-c[0]) +
         (b[1]-c[1])*(b[1]-c[1]) +
         (b[2]-c[2])*(b[2]-c[2]) > 0.0000001) {
    if (fb < fc) {
      d[0] = c[0];                              //Constrict the boundary
      d[1] = c[1];
      d[2] = c[2];
      c[0] = b[0];                              //Replace the worse point with the better pt
      c[1] = b[1];
      c[2] = b[2];
      fc = fb;
      b[0] = (b[0] + a[0]) / 2;                 //Try to improve our better pt
      b[1] = (b[1] + a[1]) / 2;
      b[2] = (b[2] + a[2]) / 2;
      MoveAtom(atom, b);
      fb = ScoreAtom(atom);
      dxb = atom->dx();                           //Differentials just for bookkeeping
      dyb = atom->dy();
      dzb = atom->dz();
    } else {
      a[0] = b[0];                              //Do the reverse of the prev section
      a[1] = b[1];
      a[2] = b[2];
      b[0] = c[0];
      b[1] = c[1];
      b[2] = c[2];
      fb = fc;
      c[0] = (c[0] + d[0]) / 2;
      c[1] = (c[1] + d[1]) / 2;
      c[2] = (c[2] + d[2]) / 2;
      MoveAtom(atom, c);
      fc = ScoreAtom(atom);
      dxc = atom->dx();
      dyc = atom->dy();
      dzc = atom->dz();
    }
  }
  if (fb < fc) {                            //Pick the better point
    MoveAtom(atom, b);
    atom->resetDelta();
    atom->addDelta(dxb, dyb, dzb);
    min = fb;                               //return the best score
    return;
  } else {
    MoveAtom(atom, c);
    atom->resetDelta();
    atom->addDelta(dxc, dyc, dzc);
    min = fc;
    return;
  }
}

void AtomRefiner::Bracket3DMin(MIAtom* atom, float* grad, float& min, float* b, float* d) {
  float tiny_step, large_step;
  float orig_coord[3];
  float inner_score, outer_score;
  float gg;

  gg = grad[0] * grad[0] +
       grad[1] * grad[1] +
       grad[2] * grad[2];

  if (gg < 0.00000001) {
    b[0] = d[0] = atom->x();
    b[1] = d[1] = atom->y();
    b[2] = d[2] = atom->z();
    return;
  }

  orig_coord[0] = atom->x();
  orig_coord[1] = atom->y();
  orig_coord[2] = atom->z();

  tiny_step = 1.0F;
  AtomStep(atom, grad, tiny_step);
  inner_score = ScoreAtom(atom);
  MoveAtom(atom, orig_coord);
  while (inner_score > min && tiny_step > 0.00001) {
    tiny_step /= 2;
    AtomStep(atom, grad, tiny_step);
    inner_score = ScoreAtom(atom);
    MoveAtom(atom, orig_coord);
  }
  if (inner_score > min) {
    tiny_step = -1.0F;
    AtomStep(atom, grad, tiny_step);
    inner_score = ScoreAtom(atom);
    MoveAtom(atom, orig_coord);
    while (inner_score > min && tiny_step < -0.00001) {
      tiny_step /= 2;
      AtomStep(atom, grad, tiny_step);
      inner_score = ScoreAtom(atom);
      MoveAtom(atom, orig_coord);
    }

    if (inner_score > min) {
      b[0] = d[0] = orig_coord[0];
      b[1] = d[1] = orig_coord[1];
      b[2] = d[2] = orig_coord[2];
      return;
    } else {
      b[0] = orig_coord[0] + tiny_step * grad[0];
      b[1] = orig_coord[1] + tiny_step * grad[1];
      b[2] = orig_coord[2] + tiny_step * grad[2];
    }
  } else {
    b[0] = orig_coord[0] + tiny_step * grad[0];
    b[1] = orig_coord[1] + tiny_step * grad[1];
    b[2] = orig_coord[2] + tiny_step * grad[2];
  }

  large_step = (tiny_step > 0) ? 2.0F : -2.0F;
  AtomStep(atom, grad, large_step);
  outer_score = ScoreAtom(atom);
  MoveAtom(atom, orig_coord);
  while (outer_score <= inner_score) {
    large_step *= 2;
    AtomStep(atom, grad, large_step);
    outer_score = ScoreAtom(atom);
    MoveAtom(atom, orig_coord);
  }
  d[0] = orig_coord[0] + large_step * grad[0];
  d[1] = orig_coord[1] + large_step * grad[1];
  d[2] = orig_coord[2] + large_step * grad[2];
}

void AtomRefiner::MoveAtom(MIAtom* atom, float* v) {
  atom->setPosition(v[0], v[1], v[2]);
}

float AtomRefiner::ScoreAtom(MIAtom* atom) {
  float score = 0.0F;
  atom->resetDelta();
  score += minimize_bonds(atom);
  score += minimize_angles(atom);
  score += minimize_planes(atom);
  score += minimize_torsions(atom);
  score += ScoreBumps(atom);
  score += ScoreChirals(atom);

  return score;
}

float AtomRefiner::minimize_bonds(MIAtom* refi_atom) {
  unsigned int i;
  MIAtom* a1, * a2;
  float r, d;
  double adjust;                //Factor by which the length will be adjusted
  //	float bweight = 1.0;
  //	float weight;
  float dsum = 0.0F;
  float dd;
  //	char buf[100];
  if (BondWeight == 0 || cons->Bonds.size() == 0) {
    return 0.0F;
  }
  //	bweight = BondWeight;

  /* find deviations */
  for (i = 0; i < cons->Bonds.size(); i++) {
    a1 = cons->Bonds[i].getAtom1();
    a2 = cons->Bonds[i].getAtom2();
    /* build differentials */
    r = (float)AtomDist(*a1, *a2);
    d = r - cons->Bonds[i].ideal_dist;
    dd = d*d;

    dsum += dd;

    adjust = 2 * d / r;
    if (refi_atom == a1) {
      a1->addDelta((float)(adjust*(a2->x() - a1->x())),
          (float)(adjust*(a2->y() - a1->y())),
          (float)(adjust*(a2->z() - a1->z())));
    } else if (refi_atom == a2) {
      a2->addDelta(-(float)(adjust*(a2->x() - a1->x())),
          -(float)(adjust*(a2->y() - a1->y())),
          -(float)(adjust*(a2->z() - a1->z())));
    }
  }

  return dsum;
}

float AtomRefiner::minimize_angles(MIAtom* refi_atom) {
  unsigned int i;
  MIAtom* a1, * a2;
  float r, d;
  double adjust;
  //	float aweight = 1.0;
  //	float weight;
  float dsum = 0.0;
  float dd;
  //	char buf[100];

  if (AngleWeight == 0 || cons->Angles.size() == 0) {
    return 0.0;
  }
  //	aweight = AngleWeight;

  /* find deviations */
  //	for(i=0;i<RefiAngles.size();i++){
  for (i = 0; i < cons->Angles.size(); i++) {
    /* angles are constrained by 1-3 distance rather
     * than by angle value
     */
    a1 = cons->Angles[i].getAtom1();
    a2 = cons->Angles[i].atom3;
    /* build derivatives */
    r = (float)AtomDist(*a1, *a2);
    d = r - cons->Angles[i].ideal_angle;
    dd = d*d;
    //	    ss = cons->Angles[i].tolerance*cons->Angles[i].tolerance;
    dsum += dd;
    //	    if ( dd/ss < 0.0001)continue;/* don't waste time on very small shifts */
    //	    if ((weight=dd/ss) > 25.0 ){
    //	        if(RefiVerbose)
    //	        printf("Very bad angle: %s %s: %s -> %s: ideal: %6.2f actual: %6.2f\n",angles[i].res->type(),angles[i].res->name(),angles[i].getAtom1()->name,angles[i].atom3->name,angles[i].ideal_angle, r);
    //	        weight = 25.0; /* clamp weight */
    //	    }
    //	    weight *= AngleWeight; /* times overall angle weight */
    adjust =  2 * d / r;
    if (refi_atom == a1) {
      a1->addDelta((float)( adjust*(a2->x() - a1->x())),
          (float)( adjust*(a2->y() - a1->y())),
          (float)( adjust*(a2->z() - a1->z())));
    } else if (refi_atom == a2) {
      a2->addDelta(-(float)(adjust*(a2->x() - a1->x())),
          -(float)(adjust*(a2->y() - a1->y())),
          -(float)(adjust*(a2->z() - a1->z())));
    }
  }
  //	if(RefiVerbose){
  //	Logger::log("RMS error in angle lengths: %0.3f (sigma =%0.3f)",sqrt(dsum/(float)nsum),sqrt(ssum/(float)nsum));
  //	}
  return dsum;
}

float AtomRefiner::minimize_planes(MIAtom* refi_atom) {
  unsigned int i;
  float d, dx, dy, dz;
  float dsum = 0.0;

  if (PlaneWeight == 0 || cons->Planes.size() == 0) {
    return 0.0;
  }

  float pln_normal[3];
  float pln_displace;

  /* find deviations */
  for (i = 0; i < cons->Planes.size(); i++) {
    /* build derivatives */
    lsqplane(cons->Planes[i], pln_normal, &pln_displace);
    /*
       printf("plane %d = %f %f %f = %f\n",i,planes[i].vm[0],planes[i].vm[1],planes[i].vm[2],planes[i].d);
     */
    d =  refi_atom->x() * pln_normal[0];
    d += refi_atom->y() * pln_normal[1];
    d += refi_atom->z() * pln_normal[2];
    d -= pln_displace;
    //	   dd = d*d;
    //	   dsum += dd;
    dsum += (int) abs((int) d);
    //	   ss = cons->Planes[i].tolerance*cons->Planes[i].tolerance;
    //	   dsum += dd; ssum += ss; nsum++;
    //	   if ( dd/ss < 0.05)continue;/* don't waste time on very small shifts */
    //	   if ((weight=dd/ss) > 5.0 ){ /* 5.0 in unsquared units */
    //	        if(RefiVerbose)printf("Very bad out of plane: %s %s: %s: %6.2f\n",planes[i].res->type(),planes[i].res->name(),planes[i].atoms[j]->name, d);
    //	        printf("Very bad out of plane: %s %s: %s: %6.2f\n",planes[i].res->type(),planes[i].res->name(),planes[i].atoms[j]->name, d);
    //		        weight = 5.0; /* clamp weight */
    //	   }
    //	   weight *= PlaneWeight;
    dx = 2 * d * pln_normal[0];
    dy = 2 * d * pln_normal[1];
    dz = 2 * d * pln_normal[2];
    //	   dx *= weight;
    //	   dy *= weight;
    //	   dz *= weight;
    refi_atom->addDelta(-dx, -dy, -dz);

    //	   refi_atom->weight += weight;
  }
  //	if(RefiVerbose){
  //	Logger::log("RMS error in planes: %0.3f (sigma =%0.3f)",sqrt(dsum/(float)nsum),sqrt(ssum/(float)nsum));
  //	}
  return dsum;
}

float AtomRefiner::minimize_torsions(MIAtom* refi_atom) {
  float dx, dy, dz;
  float chi, dchi, d;
  int i, j;
  int n = cons->Impropers.size();
  Improper* T;
  float dsum = 0.0;

  for (i = 0; i < n; i++) {
    T = &cons->Impropers[i];
    chi = (float)CalcAtomTorsion(T->GetAtom(0), T->GetAtom(1), T->GetAtom(2), T->GetAtom(3));
    if (chi < 0.0) {
      chi += 360.0;
    }
    dchi = 0.0;
    for (j = 0; j < T->NumAngles(); j++) {
      d = (float)( T->GetAngle(j) - chi);
      if (d < -180.0) {
        d += 360.0;
      }
      if (d >  180.0) {
        d -= 360.0;
      }
      if (fabs(d) < fabs(dchi) || j == 0) {
        dchi = d;
      }
    }
    /* clamp to avoid excessive movements in a single cycle
     * and to preserve small angle approximation of dx,dy,dz */
    //		if(dchi > 5.0) dchi = 5.0;
    //		if(dchi < -5.0) dchi = -5.0;

    dsum += (int) abs((int) dchi);

    //refi_atom can be at either end, but never in the middle of the torsion
    if (refi_atom == T->GetAtom(3)) {
      dTorsion(T->GetAtom(1), T->GetAtom(2), refi_atom, dchi, &dx, &dy, &dz);
      refi_atom->addDelta(dx, dy, dz);
    }
    if (refi_atom == T->GetAtom(0)) {
      dTorsion(T->GetAtom(2), T->GetAtom(1), refi_atom, dchi, &dx, &dy, &dz);
      refi_atom->addDelta(dx, dy, dz);
    }
  }
  return dsum;
}

float AtomRefiner::ScoreBumps(MIAtom*) {
  float d;
  float dsum = 0;
  //	if(BondWeight == 0 || cons->Bonds.size() == 0) return 0.0F;
  /* find clashes */
  for (unsigned int i = 0; i < cons->Bumps.size(); i++) {
    d = (float)AtomDist(*cons->Bumps[i].getAtom1(),
          *cons->Bumps[i].getAtom2());

    if (d < cons->Bumps[i].min_d) {
      dsum += 1000000 * (cons->Bumps[i].min_d - d);
    }
  }

  return dsum;
}

float AtomRefiner::ScoreChirals(MIAtom*) {
  float csum = 0;

  /* find chirality violations */
  for (unsigned int i = 0; i < cons->Chirals.size(); i++) {
    if (cons->Chirals[i].GetOrder() != cons->Chirals[i].Measure()) {
      csum += 1000000;
    }
  }

  return csum;
}

