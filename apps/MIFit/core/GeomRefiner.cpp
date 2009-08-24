#include <vector>
#include <algorithm>

#include "nonguilib.h"
#include "mathlib.h"
#include "chemlib.h"
#include "RESIDUE_.h"
#include "conflib.h"
#include "maplib.h"


#include "GeomRefiner.h"
#include "RESIDUE.h"
#include "Molecule.h"

using namespace std;
using namespace chemlib;


bool GeomRefiner::EditEntry(const char* type) {
  if (IsRefining()) {
    Logger::message("Can not edit dictionary entry while refining\nCancel Refine and start again");
    return false;
  }
  RESIDUE* dictres = dict.GetDictResidue(type, 0);
  //RESIDUE * reslist = new RESIDUE;
  RESIDUE* reslist = new RESIDUE(*dictres);
  Molecule* model = new Molecule(reslist, "Dictionary", NULL, NULL, 0, MoleculeType::Other);

  //Store state of constraint prefs to restore and
  bool tmp_ca = dict.GetConstrainCA();
  bool tmp_ends = dict.GetConstrainEnds();

  //Don't add constraints to refinement in the editor
  dict.SetConstrainCA(false);
  dict.SetConstrainEnds(false);
  SetRefiRes(reslist, reslist, model);

  //Restore state
  dict.SetConstrainCA(tmp_ca);
  dict.SetConstrainEnds(tmp_ends);

  RefiAllTorsions(reslist);
  CurrentModel->compound = "Dictionary";
  return true;
}

unsigned long GeomRefiner::FindGeomErrors(Molecule* model, float error_threshold) {
  if (dict.EmptyDictCheck() == false) {
    return 0;
  }
  if (!model) {
    return 0;
  }
  if (model->getResidues() == NULL) {
    return 0;
  }
  if (error_threshold <= 0.0) {
    error_threshold = 5.0;
  }
  bool found_geom = false;
  if (!IsRefining()) {
    MIIter<RESIDUE> ri = model->GetResidues();
    RESIDUE* firstres = ri;
    ri.Last();
    RESIDUE* last = ri;
    if (!Residue::isValid(firstres) || !Residue::isValid(last)) {
      return 0;
    }
    Logger::footer("Finding all geometry...");
    SetRefiRes(firstres, last, model);
    found_geom = true;
  }
  unsigned int i, j;
  unsigned long nadd = 0;
  std::string s;
  double d, ideal = 0, diff;
  double rmsangle, rmsdist, rmsplane;
  double sumangle = 0, sumdist = 0, sumplane = 0;
  // find bad bonds and list
  Logger::footer("Analyzing bonds...");
  for (i = 0; i < dict.RefiBonds.size(); i++) {
    d = AtomDist(*dict.RefiBonds[i].getAtom1(), *dict.RefiBonds[i].getAtom2());
    ideal = dict.RefiBonds[i].ideal_length;
    diff = fabs(ideal-d);
    sumdist += diff*diff;
    if (diff >= error_threshold* dict.RefiBonds[i].tolerance) {
      Annotation* ann = new Annotation;
      s= format("Bad bond %s-%s %0.2f(%0.2f)", dict.RefiBonds[i].getAtom1()->name(), dict.RefiBonds[i].getAtom2()->name(),
        d, ideal);
      ann->m_x = (dict.RefiBonds[i].getAtom1()->x() + dict.RefiBonds[i].getAtom2()->x())/2.0F;
      ann->m_y = (dict.RefiBonds[i].getAtom1()->y() + dict.RefiBonds[i].getAtom2()->y())/2.0F;
      ann->m_z = (dict.RefiBonds[i].getAtom1()->z() + dict.RefiBonds[i].getAtom2()->z())/2.0F;
      ann->m_text = s;
      ann->m_type = Annotation::Geom_error;
      model->addAnnotation(ann);
      nadd++;
    }
  }
  // find bad angles and list
  Logger::footer("Analyzing angles...");
  for (i = 0; i < dict.RefiAngles.size(); i++) {
    d = AtomDist(*dict.RefiAngles[i].getAtom1(), *dict.RefiAngles[i].atom3);
    ideal = dict.RefiAngles[i].ideal_angle;
    diff = fabs(ideal-d);
    sumangle += diff*diff;
    if (diff >= error_threshold*dict.RefiAngles[i].tolerance) {
      Annotation* ann = new Annotation;
      s=format("Bad angle %s-%s-%s %0.2f(%0.2f)", dict.RefiAngles[i].getAtom1()->name(), dict.RefiAngles[i].getAtom2()->name(), dict.RefiAngles[i].atom3->name(),
        d, ideal);
      ann->m_x = (dict.RefiAngles[i].getAtom1()->x() + dict.RefiAngles[i].getAtom2()->x() + dict.RefiAngles[i].atom3->x())/3.0F;
      ann->m_y = (dict.RefiAngles[i].getAtom1()->y() + dict.RefiAngles[i].getAtom2()->y() + dict.RefiAngles[i].atom3->y())/3.0F;
      ann->m_z = (dict.RefiAngles[i].getAtom1()->z() + dict.RefiAngles[i].getAtom2()->z() + dict.RefiAngles[i].atom3->z())/3.0F;
      ann->m_text = s;
      ann->m_type = Annotation::Geom_error;
      model->addAnnotation(ann);
      nadd++;
    }
  }
  // find bad planes and list
  Logger::footer("Analyzing planes...");
  MIAtom* a1;
  for (i = 0; i < dict.RefiPlanes.size(); i++) {
    /* build derivatives */
    lsqplane(dict.RefiPlanes[i]);
    int nsum = 0;
    float dsum = 0;
    float x = 0, y = 0, z = 0;
    for (int j = 0; j < (int)dict.RefiPlanes[i].natoms; j++) {
      a1 = dict.RefiPlanes[i].atoms[j];
      d =  a1->x() * dict.RefiPlanes[i].vm[0];
      d += a1->y() * dict.RefiPlanes[i].vm[1];
      d += a1->z() * dict.RefiPlanes[i].vm[2];
      d -= dict.RefiPlanes[i].d;
      d = d*d;
      dsum += d;
      nsum++;
      x += a1->x();
      y += a1->y();
      z += a1->z();
    }
    diff = sqrt(dsum/(double)nsum);
    sumplane += diff*diff;
    if (diff >= error_threshold*dict.RefiPlanes[i].tolerance) {
      Annotation* ann = new Annotation;
      s=format("Bad plane %0.2f", diff);
      ann->m_x = x/(float)dict.RefiPlanes[i].natoms;
      ann->m_y = y/(float)dict.RefiPlanes[i].natoms;
      ann->m_z = z/(float)dict.RefiPlanes[i].natoms;
      ann->m_text = s;
      ann->m_type = Annotation::Geom_error;
      model->addAnnotation(ann);
      nadd++;
    }
  }
  // analyze torsions
  Logger::footer("Analyzing torsions...");
  float chi, dchi;
  for (i = 0; i < dict.RefiTorsions.size(); i++) {
    chi = CalcAtomTorsion(dict.RefiTorsions[i].getAtom1(), dict.RefiTorsions[i].getAtom2(), dict.RefiTorsions[i].atom3, dict.RefiTorsions[i].atom4);
    if (chi < 0.0) {
      chi += 360.0;
    }
    dchi = 0.0;
    for (j = 0; j < (unsigned int)dict.RefiTorsions[i].nideal; j++) {
      d = dict.RefiTorsions[i].ideal[j] - chi;
      if (d < -180.0) {
        d += 360.0;
      }
      if (d >  180.0) {
        d -= 360.0;
      }
      if (fabs(d) < fabs(dchi) || j == 0) {
        dchi = d;
        ideal = dict.RefiTorsions[i].ideal[j];
      }
    }
    if (fabs(dchi) >= error_threshold*dict.GetSigmaTorsion()) {
      Annotation* ann = new Annotation;
      s=format("Bad torsion %s %s-%s-%s-%s %0.2f(%0.2f)",
        dict.RefiTorsions[i].type, dict.RefiTorsions[i].getAtom1()->name(), dict.RefiTorsions[i].getAtom2()->name(),
        dict.RefiTorsions[i].atom3->name(), dict.RefiTorsions[i].atom4->name(), chi, ideal);
      ann->m_x = (dict.RefiTorsions[i].getAtom2()->x() + dict.RefiTorsions[i].atom3->x())/2.0F;
      ann->m_y = (dict.RefiTorsions[i].getAtom2()->y() + dict.RefiTorsions[i].atom3->y())/2.0F;
      ann->m_z = (dict.RefiTorsions[i].getAtom2()->z() + dict.RefiTorsions[i].atom3->z())/2.0F;
      ann->m_text = s;
      ann->m_type = Annotation::Geom_error;
      model->addAnnotation(ann);
      nadd++;
    }
  }
  float phi, psi, omega, omegap;
  for (i = 0; i < dict.RefiPhiPsis.size(); i += 4) {
    chi = CalcAtomTorsion(dict.RefiPhiPsis[i].getAtom1(), dict.RefiPhiPsis[i].getAtom2(), dict.RefiPhiPsis[i].atom3, dict.RefiPhiPsis[i].atom4);
    if (chi < 0.0) {
      chi += 360.0;
    }
    dchi = 0.0;
    phi = CalcAtomTorsion(dict.RefiPhiPsis[i].getAtom1(), dict.RefiPhiPsis[i].getAtom2(), dict.RefiPhiPsis[i].atom3, dict.RefiPhiPsis[i].atom4);
    psi = CalcAtomTorsion(dict.RefiPhiPsis[i+1].getAtom1(), dict.RefiPhiPsis[i+1].getAtom2(), dict.RefiPhiPsis[i+1].atom3, dict.RefiPhiPsis[i+1].atom4);
    omega = CalcAtomTorsion(dict.RefiPhiPsis[i+2].getAtom1(), dict.RefiPhiPsis[i+2].getAtom2(), dict.RefiPhiPsis[i+2].atom3, dict.RefiPhiPsis[i+2].atom4);
    omegap = CalcAtomTorsion(dict.RefiPhiPsis[i+3].getAtom1(), dict.RefiPhiPsis[i+3].getAtom2(), dict.RefiPhiPsis[i+3].atom3, dict.RefiPhiPsis[i+3].atom4);

    d = phipsi_energy(phi, psi);
    if (d <= 3.0F && strcmp(dict.RefiPhiPsis[i].res->type().c_str(), "GLY") != 0) {
      Annotation* ann = new Annotation;
      s=format("Bad phi-psi %s %0.1f, %0.1f", resid(dict.RefiPhiPsis[i].res).c_str(), phi, psi);
      ann->m_x = dict.RefiPhiPsis[i].atom3->x();
      ann->m_y = dict.RefiPhiPsis[i].atom3->y();
      ann->m_z = dict.RefiPhiPsis[i].atom3->z();
      ann->m_text = s;
      ann->m_type = Annotation::Geom_error;
      model->addAnnotation(ann);
      nadd++;
    }
    for (j = 0; j < (unsigned int)dict.RefiPhiPsis[i+2].nideal; j++) {
      d = dict.RefiPhiPsis[i+2].ideal[j] - omega;
      if (d < -180.0) {
        d += 360.0;
      }
      if (d >  180.0) {
        d -= 360.0;
      }
      if (fabs(d) < fabs(dchi) || j == 0) {
        dchi = d;
        ideal = dict.RefiPhiPsis[i+2].ideal[j];
      }
    }
    if (fabs(dchi) >= error_threshold*dict.GetSigmaTorsion()) {
      Annotation* ann = new Annotation;
      s=format("Bad phi-psi %s %0.2f(180.0)", dict.RefiPhiPsis[i+2].type, dchi);
      ann->m_x = (dict.RefiPhiPsis[i+2].getAtom2()->x() + dict.RefiPhiPsis[i+2].atom3->x())/2.0F;
      ann->m_y = (dict.RefiPhiPsis[i+2].getAtom2()->y() + dict.RefiPhiPsis[i+2].atom3->y())/2.0F;
      ann->m_z = (dict.RefiPhiPsis[i+2].getAtom2()->z() + dict.RefiPhiPsis[i+2].atom3->z())/2.0F;
      ann->m_text = s;
      ann->m_type = Annotation::Geom_error;
      model->addAnnotation(ann);
      nadd++;
    }
    for (j = 0; j < (unsigned int)dict.RefiPhiPsis[i+3].nideal; j++) {
      d = dict.RefiPhiPsis[i+3].ideal[j] - omegap;
      if (d < -180.0) {
        d += 360.0;
      }
      if (d >  180.0) {
        d -= 360.0;
      }
      if (fabs(d) < fabs(dchi) || j == 0) {
        dchi = d;
        ideal = dict.RefiPhiPsis[i+3].ideal[j];
      }
    }
    if (fabs(dchi) >= error_threshold*dict.GetSigmaTorsion()) {
      Annotation* ann = new Annotation;
      s=format("Bad phi-psi %s %0.2f(0.0)", dict.RefiPhiPsis[i+3].type, dchi);
      ann->m_x = (dict.RefiPhiPsis[i+3].getAtom2()->x() + dict.RefiPhiPsis[i+3].atom3->x())/2.0F;
      ann->m_y = (dict.RefiPhiPsis[i+3].getAtom2()->y() + dict.RefiPhiPsis[i+3].atom3->y())/2.0F;
      ann->m_z = (dict.RefiPhiPsis[i+3].getAtom2()->z() + dict.RefiPhiPsis[i+3].atom3->z())/2.0F;
      ann->m_text = s;
      ann->m_type = Annotation::Geom_error;
      model->addAnnotation(ann);
      nadd++;
    }
  }
  // find bad bumps and list
  for (i = 0; i < dict.RefiBumps.size(); i++) {
    MIAtom* atom1 = dict.RefiBumps[i].getAtom1();
    MIAtom* atom2 = dict.RefiBumps[i].getAtom2();
    d = AtomDist(*atom1, *atom2);
    ideal = dict.RefiBumps[i].ideal_length;
    diff = d - ideal;
    if (diff <= -error_threshold* dict.RefiBumps[i].tolerance) {
      if (atom1->altloc() != ' ' && atom2->altloc() != ' ') {
        continue;
      }
      if (atom1->altloc() != ' ' || atom2->altloc() != ' ') {
        MIIter<RESIDUE> ri = model->GetResidues();
        RESIDUE* res1 = residue_from_atom(ri, atom1);
        RESIDUE* res2 = residue_from_atom(ri, atom2);
        if (res1 == res2) {
          continue;
        }
      }
      Annotation* ann = new Annotation;
      s=format("Bad bump %s-%s %0.2f(%0.2f)", atom1->name(), atom2->name(), d, ideal);
      ann->m_x = (atom1->x() + atom2->x())/2.0F;
      ann->m_y = (atom1->y() + atom2->y())/2.0F;
      ann->m_z = (atom1->z() + atom2->z())/2.0F;
      ann->m_text = s;
      ann->m_type = Annotation::Geom_error;
      model->addAnnotation(ann);
      nadd++;
    }
  }

  if (dict.RefiBonds.size() > 0) {
    rmsdist = sqrt(sumdist/(double)dict.RefiBonds.size());
    s=format("Rms deviation of bonds from ideal: %0.3f A", rmsdist);
    Logger::log(s);
  }
  if (dict.RefiAngles.size() > 0) {
    rmsangle = sqrt(sumangle/(double)dict.RefiAngles.size());
    s=format("Rms deviation of angles (1-3 distance) from ideal: %0.3f A", rmsangle);
    Logger::log(s);
  }
  if (dict.RefiPlanes.size() > 0) {
    rmsplane = sqrt(sumplane/(double)dict.RefiPlanes.size());
    s=format("Rms deviation of planes from ideal: %0.3f A", rmsplane);
    Logger::log(s);
  }

  if (found_geom) {
    Cancel();
  }

  return nadd;
}

void GeomRefiner::EditEntryCleanup(bool /*isok*/) {
  MIMoleculeBase* model = CurrentModel;   // Just incase the model was rebuilt in the editor
  delete model;
  SaveToken--;

  unlockRefineTarget();
  clearRefineTarget();
}

vector<std::string> GetDictResList(GeomRefiner* gr) {
  std::vector<std::string> reslist;
  std::vector<std::string> result;
  reslist = gr->dict.GetDictResList();
  for (unsigned int i = 0; i < reslist.size(); ++i) {
    result.push_back(reslist[i].c_str());
  }
  return result;
}

