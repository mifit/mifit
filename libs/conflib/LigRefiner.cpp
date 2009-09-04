#include <nongui/nonguilib.h>
#include <vector>
#include <fstream>
#include <sstream>

#include <chemlib/chemlib.h>

#include "LigRefiner.h"

//Adapted from Duncan's (MIFit) GeomRefiner Class

using namespace chemlib;
using namespace conflib;

LigRefiner::LigRefiner() {
  BondWeight = AngleWeight = PlaneWeight = BumpWeight = TorsionWeight = 1.0F;
  nCycles = 500;
}

double LigRefiner::Refine() {
  //	if(!IsRefining()) return;
  // does one cycle of refinement
  double score = 0.0;

  for (int i = 0; i < nCycles; i++) {
    resetderivatives();

    score = 0;
    score += minimize_bonds(&RefiBonds[0], RefiBonds.size());
    score += minimize_angles(&RefiAngles[0], RefiAngles.size());
    score += minimize_planes(&RefiPlanes[0], RefiPlanes.size());
    score += minimize_torsions();
    score += minimize_bumps(&RefiBumps[0], RefiBumps.size());

    applyderivatives();
  }
  return score;
}

double LigRefiner::minimize_bonds(BondLength* bonds, unsigned int nbonds) {
  unsigned int i;
  MIAtom* a1, * a2;
  float r, d, dx, dy, dz;
  float bweight = 1.0, weight;
  float dsum = 0.0, ssum = 0.0;
  float dd, ss;
  int nsum = 0;
  //	char buf[100];
  if (BondWeight == 0 || nbonds == 0) {
    return (0);
  }
  bweight = BondWeight;

  /* find deviations */
  for (i = 0; i < nbonds; i++) {
    a1 = bonds[i].getAtom1();
    a2 = bonds[i].getAtom2();
    /* build derivatives */
    r = (float)AtomDist(*a1, *a2);
    d = r - bonds[i].ideal_dist;
    dd = d*d;
    ss = bonds[i].tolerance*bonds[i].tolerance;
    dsum += dd; ssum += ss; nsum++;
    if (dd/ss < 0.05) {
      continue;                    /* don't waste time on very small shifts */
    }
    if ((weight = dd/ss) > 25.0) {   /* 5.0 in unsquared units */
      weight = 25.0;       /* clamp weight */
    }
    weight *= bweight;
    dx = d*(a2->x() - a1->x())/r/2.0f;     /*move each atom half the distance*/
    dy = d*(a2->y() - a1->y())/r/2.0f;     /*so that the total movement is d */
    dz = d*(a2->z() - a1->z())/r/2.0f;
    dx *= weight;
    dy *= weight;
    dz *= weight;
    a1->addDelta(dx, dy, dz);
    a2->addDelta(-dx, -dy, -dz);
    a1->addWeight(weight);
    a2->addWeight(weight);
  }
  //	if(RefiVerbose){
  //	Logger::log("RMS error in bond lengths: %0.3f (sigma =%0.3f)",sqrt(dsum/(float)nsum),sqrt(ssum/(float)nsum));
  //	}
  //	return(RefiBonds.size());
  return dsum;
}

double LigRefiner::minimize_angles(Angle* angles, unsigned int nangles) {
  unsigned int i;
  MIAtom* a1, * a2;
  float r, d, dx, dy, dz;
  float aweight = 1.0;
  float weight;
  float dsum = 0.0, ssum = 0.0;
  float dd, ss;
  int nsum = 0;
  //	char buf[100];

  if (AngleWeight == 0 || nangles == 0) {
    return (0);
  }
  aweight = AngleWeight;

  /* find deviations */
  //	for(i=0;i<RefiAngles.size();i++){
  for (i = 0; i < nangles; i++) {
    /* angles are constrained by 1-3 distance rather
     * than by angle value
     */
    a1 = angles[i].getAtom1();
    a2 = angles[i].atom3;
    /* build derivatives */
    r = (float)AtomDist(*a1, *a2);
    d = r - angles[i].ideal_angle;
    dd = d*d;
    ss = angles[i].tolerance*angles[i].tolerance;
    dsum += dd; ssum += ss; nsum++;
    if (dd/ss < 0.05) {
      continue;                    /* don't waste time on very small shifts */
    }
    if ((weight = dd/ss) > 25.0) {
      //	        if(RefiVerbose)
      //	        printf("Very bad angle: %s %s: %s -> %s: ideal: %6.2f actual: %6.2f\n",angles[i].res->type(),angles[i].res->name(),angles[i].getAtom1()->name,angles[i].atom3->name,angles[i].ideal_angle, r);
      weight = 25.0;       /* clamp weight */
    }
    weight *= aweight;     /* times overall angle weight */
    dx = d*(a2->x() - a1->x())/r/2.0f;
    dy = d*(a2->y() - a1->y())/r/2.0f;
    dz = d*(a2->z() - a1->z())/r/2.0f;
    dx *= weight;
    dy *= weight;
    dz *= weight;
    a1->addDelta(dx, dy, dz);
    a2->addDelta(-dx, -dy, -dz);
    a1->addWeight(weight);
    a2->addWeight(weight);
  }
  //	if(RefiVerbose){
  //	Logger::log("RMS error in angle lengths: %0.3f (sigma =%0.3f)",sqrt(dsum/(float)nsum),sqrt(ssum/(float)nsum));
  //	}
  //	return(RefiAngles.size());
  return dsum;
}

double LigRefiner::minimize_planes(Plane* planes, unsigned int nplanes) {
  unsigned int i, j;
  MIAtom* a1;
  float d, dx, dy, dz;
  float bweight = 1.0, weight;
  float dsum = 0.0, ssum = 0.0;
  float dd, ss;
  int nsum = 0;
  //	char buf[130];

  float pln_normal[3];
  float pln_displace;

  if (PlaneWeight == 0 || nplanes == 0) {
    return (0);
  }
  bweight = PlaneWeight;

  /* find deviations */
  for (i = 0; i < nplanes; i++) {
    /* build derivatives */
    lsqplane(planes[i], pln_normal, &pln_displace);
    /*
       printf("plane %d = %f %f %f = %f\n",i,planes[i].vm[0],planes[i].vm[1],planes[i].vm[2],planes[i].d);
     */
    for (j = 0; j < (unsigned int)planes[i].NumAtoms(); j++) {
      a1 = planes[i].GetAtom(j);
      d =  a1->x() * pln_normal[0];
      d += a1->y() * pln_normal[1];
      d += a1->z() * pln_normal[2];
      d -= pln_displace;
      dd = d*d;
      ss = planes[i].GetTolerance()*planes[i].GetTolerance();
      dsum += dd; ssum += ss; nsum++;
      if (dd/ss < 0.05) {
        continue;                  /* don't waste time on very small shifts */
      }
      if ((weight = dd/ss) > 5.0) { /* 5.0 in unsquared units */
        //	        if(RefiVerbose)printf("Very bad out of plane: %s %s: %s: %6.2f\n",planes[i].res->type(),planes[i].res->name(),planes[i].atoms[j]->name, d);
        //	        printf("Very bad out of plane: %s %s: %s: %6.2f\n",planes[i].res->type(),planes[i].res->name(),planes[i].atoms[j]->name, d);
        weight = 3.0;         /* clamp weight */
      }
      weight *= bweight;
      dx = d * pln_normal[0];
      dy = d * pln_normal[1];
      dz = d * pln_normal[2];
      dx *= weight;
      dy *= weight;
      dz *= weight;
      a1->addDelta(-dx, -dy, -dz);
      a1->addWeight(weight);
    }
  }
  //	if(RefiVerbose){
  //	Logger::log("RMS error in planes: %0.3f (sigma =%0.3f)",sqrt(dsum/(float)nsum),sqrt(ssum/(float)nsum));
  //	}
  //	return(nplanes);
  return dsum;
}

double LigRefiner::minimize_bumps(Bump* bumps, unsigned int nbumps) {
  unsigned int i;
  MIAtom* a1, * a2;
  float r, d, dx, dy, dz;
  float bweight = 1.0, weight;
  int sliderweight = 10;
  float dsum = 0.0, ssum = 0.0;
  float dd, ss;
  int nsum = 0;
  //	int color = WHITE;
  float mx, my, mz;
  if (sliderweight == 0 || nbumps == 0) {
    return (0);
  }
  bweight = (float)sliderweight/10.0f;

  //	Vu.clear();
  /* find deviations */
  for (i = 0; i < nbumps; i++) {
    a1 = bumps[i].getAtom1();
    a2 = bumps[i].getAtom2();
    /* build derivatives */
    r = (float)AtomDist(*a1, *a2);
    d = r - bumps[i].min_d;
    if (d > 0.0) {
      continue;
    }
    dd = d*d;
    ss = bumps[i].tolerance*bumps[i].tolerance;
    if (d < -2.0*bumps[i].tolerance) {
      weight = 2*dd/ss;
      //			color = RED;
    } else {
      weight = dd/ss;
      //			color = WHITE;
    }
    mx = (a2->x()+a1->x())/2.0f;
    my = (a2->y()+a1->y())/2.0f;
    mz = (a2->z()+a1->z())/2.0f;
    dsum += dd; ssum += ss; nsum++;
    /* don't waste time on very small shifts */
    if (dd/ss < 0.05) {
      continue;
    }
    weight *= bweight;
    dx = d*(a2->x() - a1->x())/r/2.0f;     /*move each atom half the distance*/
    dy = d*(a2->y() - a1->y())/r/2.0f;     /*so that the total movement is d */
    dz = d*(a2->z() - a1->z())/r/2.0f;
    dx *= weight;
    dy *= weight;
    dz *= weight;
    a1->addDelta(dx, dy, dz);
    a2->addDelta(-dx, -dy, -dz);
    a1->addWeight(weight);
    a2->addWeight(weight);
  }
  //	return(nbumps);
  return dsum;
}

double LigRefiner::minimize_torsions() {
  float w;
  float dx, dy, dz;
  float chi, dchi, d;
  int sweight = 10;
  float bweight;
  float dsum = 0.0;
  int i, j;
  int n = RefiTorsions.size();
  bweight = (float)sweight/10.0f;
  for (i = 0; i < n; i++) {
    chi = (float)CalcAtomTorsion(RefiTorsions[i].GetAtom(0),
            RefiTorsions[i].GetAtom(1),
            RefiTorsions[i].GetAtom(2),
            RefiTorsions[i].GetAtom(3));
    if (chi < 0.0) {
      chi += 360.0;
    }
    dchi = 0.0;
    for (j = 0; j < RefiTorsions[i].NumAngles(); j++) {
      d = (float)RefiTorsions[i].GetAngle(j) - chi;
      if (d < -180.0f) {
        d += 360.0f;
      }
      if (d >  180.0f) {
        d -= 360.0f;
      }
      if (fabs(d) < fabs(dchi) || j == 0) {
        dchi = d;
      }
    }
    if (fabs(dchi) > 3.0) {
      /* clamp to avoid excessive movements in a single cycle
       * and to preserve small angle approximation of dx,dy,dz */
      if (dchi > 5.0) {
        dchi = 5.0;
      }
      if (dchi < -5.0) {
        dchi = -5.0;
      }
      w  = dchi*dchi/100.0f*bweight;
      dTorsion(RefiTorsions[i].GetAtom(1),
        RefiTorsions[i].GetAtom(2),
        RefiTorsions[i].GetAtom(3),
        dchi,
        &dx, &dy, &dz);
      dsum += dx*dx + dy*dy + dz*dz;
      dx *= w;
      dy *= w;
      dz *= w;
      RefiTorsions[i].GetAtom(3)->addDelta(dx, dy, dz);
      RefiTorsions[i].GetAtom(3)->addWeight(w);
    }
  }
  //	return n;
  return 0.0;
}

int LigRefiner::resetderivatives() {
  if (RefiRes.size() == 0) {
    return 0;
  }
  Residue* res;
  unsigned int i;
  unsigned int n = 0;
  while (n < RefiRes.size()) {
    res = RefiRes[n];
    for (i = 0; i < res->atoms().size(); i++) {
      res->atom(i)->resetDelta();
      res->atom(i)->resetWeight();
    }
    n++;
  }
  return 1;
}

int LigRefiner::applyderivatives() {
  if (RefiRes.size() == 0) {
    return 0;
  }
  Residue* res;
  unsigned int i;
  MIAtom* a;
  float dx, dy, dz;
  float maxd = 0.3F;
  unsigned int n = 0;
  while (n < RefiRes.size()) {
    res = RefiRes[n];
    for (i = 0; i < res->atoms().size(); i++) {
      a = res->atom(i);
      if (a->weight() < .05) {
        a->resetDelta();
        a->resetWeight();
        a->addWeight(1.0f);
        continue;
      }
      dx = a->dx()/a->weight();
      if (dx > maxd) {
        dx = maxd;
      }
      if (dx < -maxd) {
        dx = -maxd;
      }
      dy = a->dy()/a->weight();
      if (dy > maxd) {
        dy = maxd;
      }
      if (dy < -maxd) {
        dy = -maxd;
      }
      dz = a->dz()/a->weight();
      if (dz > maxd) {
        dz = maxd;
      }
      if (dz < -maxd) {
        dz = -maxd;
      }
      a->translate(dx, dy, dz);
      a->resetDelta();
      a->addDelta(dx, dy, dz);
      a->resetWeight();
      a->addWeight(1.0f);
    }
    n++;
  }
  return 1;
}

void LigRefiner::AddRefiRes(Residue* res) {
  RefiRes.push_back(res);
}

long LigRefiner::SetRefiRes(std::vector<Residue*>::iterator res1,
                            std::vector<Residue*>::iterator res2) {
  RefiRes.assign(res1, res2);
  return RefiRes.size();
}

void LigRefiner::SetModel(Ligand* model) {
  CurrentModel = model;
  RefiBonds = CurrentModel->geometry.Bonds;
  RefiAngles = CurrentModel->geometry.Angles;
  RefiPlanes = CurrentModel->geometry.Planes;
  RefiTorsions = CurrentModel->geometry.Impropers;
  BuildBumps();
}

inline bool LigandRefiner_IsHydrogen(const MIAtom* atom) {
  return (atom->atomicnumber() == 1);
}

int LigRefiner::BuildBumps() {
  Residue* res, * res2;
  MIAtom* a1, * a2;
  Bump bond;
  unsigned int n2;
  unsigned int n;
  unsigned int i, j;
  int k;
  int found;
  float d;
  //	int nres = RefiRes;
  //	Residue * reslist = RefiRes;
  //	res = reslist;
  /* loop thru every atom in the reslist and add a bump if
   * reasonably close */
  n = 0;
  RefiBumps.clear();
  while (n < RefiRes.size()) {
    res = RefiRes[n];
    for (i = 0; i < res->atoms().size(); i++) {
      n2 = 0;
      a1 = res->atom(i);
      while (n2 < RefiRes.size()) {
        res2 = RefiRes[n2];
        for (j = 0; j < res2->atoms().size(); j++) {
          a2 = res2->atom(j);
          if (a1 < a2
              && !LigandRefiner_IsHydrogen(a1) && !LigandRefiner_IsHydrogen(a2)
              && (d = (float)AtomDist(*a1, *res2->atom(j))) < 4.3f) {
            /* are they bonded ? */
            found = 0;
            //						if(res==res2 && d<2.95 && !MIIsMainChainAtom(a1) && !MIIsMainChainAtom(a2)) found =1;
            //Will need				if(res->next())
            //to fix for				if(res->next()==res2 && d<3.5){
            //multi-res ligs				if(MIIsMainChainAtom(a1)&&MIIsMainChainAtom(a2))found=1;
            //						}
            //						if(found ==0)
            for (k = 0; (unsigned int)k < RefiBonds.size(); k++) {
              if (a1 == RefiBonds[k].getAtom1()) {
                if (a2 == RefiBonds[k].getAtom2()) {
                  found = 1;
                  break;
                }
              }
              if (a1 == RefiBonds[k].getAtom2()) {
                if (a2 == RefiBonds[k].getAtom1()) {
                  found = 1;
                  break;
                }
              }
            }
            if (found == 0) {
              for (k = 0; (unsigned int)k < RefiAngles.size(); k++) {
                if (a1 == RefiAngles[k].getAtom1()) {
                  if (a2 == RefiAngles[k].atom3) {
                    found = 1;
                    break;
                  }
                }
                if (a1 == RefiAngles[k].atom3) {
                  if (a2 == RefiAngles[k].getAtom1()) {
                    found = 1;
                    break;
                  }
                }
              }
            }
            if (found == 0) {
              bond.setAtom1(a1);
              bond.setAtom2(a2);
              if ((a1->name()[0] == 'O' && a2->name()[0] == 'N')
                  || (a2->name()[0] == 'O' && a1->name()[0] == 'N')) {
                bond.min_d = 2.75F;
                bond.tolerance = 0.3F;
              } else {
                bond.min_d = 3.1F;
                bond.tolerance = 0.1F;
              }
              RefiBumps.push_back(bond);
            }
          }
        }
        n2++;
      }
    }
    n++;
  }
  return 1;
}

void LigRefiner::Clear() {
  RefiAngles.clear();
  RefiBonds.clear();
  RefiPlanes.clear();
  RefiTorsions.clear();
  RefiBumps.clear();
}

