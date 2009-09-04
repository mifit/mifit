#include <chemlib/chemlib.h>

#include "Atom.h"
#include "DoubleBond.h"
#include "Drawing.h"
#include "FigRefiner.h"
#include "HBond.h"
#include "HPhobe.h"
#include "HydrogenBond.h"
#include "Label.h"
#include "Shape.h"
#include "SingleBond.h"
#include "Spokes.h"
#include "UnifiedRes.h"

#include "figure.h"


#include <numeric>
#include <functional>

#define SQR(a) ((a)*(a))



//Parameters for placing whole residues (without showing individual atoms), which
//are placed for residues which have hydrophobic contacts to the ligand but no
//hydrogen bonds to the ligand
#define SUN_ATOM_DIST 3.0F
#define SQ_SUN_ATOM_DIST (SUN_ATOM_DIST*SUN_ATOM_DIST)
#define SUN_SUN_DIST 4.0F
#define SQ_SUN_SUN_DIST (SUN_SUN_DIST*SUN_SUN_DIST)

#define RESLABEL_ATOM_DIST 0.30F
#define SQ_RESLABEL_ATOM_DIST (RESLABEL_ATOM_DIST*RESLABEL_ATOM_DIST)
#define RESLABEL_SUN_DIST 3.0F
#define SQ_RESLABEL_SUN_DIST (RESLABEL_ATOM_DIST*RESLABEL_ATOM_DIST)


//using namespace chemlib;



namespace moldraw {

bool DrawSiteASP(Drawing *dp,
                 std::vector<chemlib::Residue*>& lig,
                 std::vector<chemlib::Residue*>& rec,
                 std::vector<chemlib::Bond>& bonds,
                 std::vector<HBond>& hbonds,
                 std::vector<HPhobe>& hphobes,
                 std::vector<chemlib::Residue*>& suns);

void HBondsToLigand(std::vector<chemlib::Residue*>& lig,
                    std::vector<chemlib::Residue*>& rec,
                    std::vector<chemlib::Bond>& bonds,
                    std::vector<HBond>& hbonds,
                    std::vector<chemlib::Residue*>& site);

void HphobesToLigand(std::vector<chemlib::Residue*>& lig,
                     std::vector<chemlib::Residue*>& rec,
                     std::vector<HPhobe>& hphobe,
                     std::vector<chemlib::Residue*>& site);

void Draw(std::vector<chemlib::Residue*>& residues,
          std::vector<chemlib::Bond>& bonds,
          const std::vector<UnifiedRes>& suns,
          const std::vector<HPhobe>& hphobes,
          const std::vector<std::string>& lig_names,
          Drawing* dp);

void DrawAtomLabels(const chemlib::Residue& res, Drawing* dp);

void DrawResidueLabel(const chemlib::Residue& res,
                      chemlib::Ligand& site,
                      std::vector<UnifiedRes>& suns,
                      Drawing* dp);

void PlaceSuns(std::vector<UnifiedRes>& placed_suns,
               std::vector<HPhobe>& hphobes,
               std::vector<chemlib::Residue*>& suns,
               chemlib::Ligand& site);

double ScoreSun(double* pos,
                std::vector<UnifiedRes>& placed_suns,
                chemlib::Ligand& site);

float ScoreResLabel(float* pos,
                    std::vector<UnifiedRes>& suns,
                    chemlib::Ligand& site);

PaletteColor GetAtomPlotColor(const chemlib::MIAtom& atom);



bool GenerateFigureASP(Drawing *dp,
                       chemlib::MIMoleculeBase* mol,
                       const char* lig_resnumber,
                       const unsigned short chain_id) {
  const std::vector<chemlib::Bond>& bonds = mol->getBonds();


  chemlib::Ligand tmp_mol;


  int nresidues = 0;
  if (mol) {
    nresidues = mol->getnresidues();
  }
  if (nresidues == 0) {
    return false;
  }
  tmp_mol.residues.reserve(nresidues);
  for (MIIter<chemlib::RESIDUE> r = mol->GetResidues(); r; ++r) {
    tmp_mol.AddRes(r, bonds);
  }

  std::vector<chemlib::Bond>::iterator bnd;
  for (bnd = tmp_mol.bonds.begin(); bnd != tmp_mol.bonds.end(); ++bnd) {
    bnd->setOrder(1);
  }

  std::vector<chemlib::Residue*> lig;
  std::vector<chemlib::Residue*> rec;
  unsigned int i;
  for (i = 0; i < tmp_mol.residues.size(); ++i) {
    if ((strcmp(tmp_mol.residues[i]->name().c_str(), lig_resnumber) == 0)
        && (chain_id == tmp_mol.residues[i]->chain_id())) {
      lig.push_back(tmp_mol.residues[i]);
    } else {
      rec.push_back(tmp_mol.residues[i]);
    }
  }

  std::vector<chemlib::Residue*> site;                  //List of residues that hbond to the ligand
  std::vector<chemlib::Residue*> suns;                  //List of residues that (only) hydrophobically
                                                        //interact with the ligand
  //Gather hydrogen-bonded residues
  std::vector<HBond> hbonds;
  HBondsToLigand(lig,
                 rec,
                 tmp_mol.bonds,
                 hbonds,
                 site);

  //Gather residues hydrophobically "bonded" to the ligand
  std::vector<HPhobe> hphobes;
  HphobesToLigand(lig,
                  rec,
                  hphobes,
                  suns);
  /*
    std::vector<chemlib::Residue *> tmp_site;
    for(i = 0; i<tmp_mol.residues.size(); ++i) {
    if((tmp_mol.residues[i].name == "46") &&
    (tmp_mol.residues[i].type == "HOH")) {
    tmp_site.push_back(&tmp_mol.residues[i]);
    }
    }
    std::vector<chemlib::Residue *> tmp_site2;
    for(i = 0; i<tmp_mol.residues.size(); ++i) {
    if((tmp_mol.residues[i].name == "196") &&
    (tmp_mol.residues[i].type == "PHE")) {
    tmp_site2.push_back(&tmp_mol.residues[i]);
    }
    }

    HBondsToLigand(tmp_site,
    tmp_site2,
    tmp_mol.bonds,
    hbonds,
    site);
  */
  return DrawSiteASP(dp, lig, site, tmp_mol.bonds, hbonds, hphobes, suns);
}


bool DrawSiteASP(Drawing *dp,
                 std::vector<chemlib::Residue*>& lig,
                 std::vector<chemlib::Residue*>& rec,
                 std::vector<chemlib::Bond>& bonds,
                 std::vector<HBond> &hbonds,
                 std::vector<HPhobe> &hphobes,
                 std::vector<chemlib::Residue*>& suns) {

  //Create a vector of covalent & hydrogen bonds
  std::vector<chemlib::Bond> all_bonds = bonds;             //covalent & hydrogen bonds
  chemlib::Bond* tmp_bond;
  std::vector<HBond>::iterator hbnd;
  for (hbnd = hbonds.begin(); hbnd != hbonds.end(); ++hbnd) {
    tmp_bond = new chemlib::Bond(hbnd->donor, hbnd->acceptor, HYDROGENBOND, 0);
    tmp_bond->ideal_length = AtomDist(*hbnd->donor, *hbnd->acceptor);
    all_bonds.push_back(*tmp_bond);
    delete tmp_bond;
  }

  //Create a vector of all the residues to
  //include in the molecule
  std::vector<chemlib::Residue*> site_res;
  site_res.insert(site_res.begin(), lig.begin(), lig.end());
  site_res.insert(site_res.begin(), rec.begin(), rec.end());

  std::vector<std::string> lig_names;
  unsigned int i;
  for (i = 0; i < lig.size(); ++i) {
    lig_names.push_back(lig[i]->name().c_str());
  }

  //Create a new molecule
  chemlib::Ligand site(site_res, all_bonds);

  //Extract hydrophobic bonds from the ligand

  std::vector<HPhobe>::iterator hp;
  for (hp = hphobes.begin(); hp != hphobes.end(); ++hp) {
    RePointHphobe(*hp, site);
  }

  //Sever a subset of the hydrogen bonds (those that create loops)

  //Flatten the molecule
  site.FindRingSystems();
  site.Flatten();



  //Optimize the figure
  FigOpt(site, hbonds);

  //
  //Restore the hydrogen bonds
  std::vector<UnifiedRes> placed_suns;

  PlaceSuns(placed_suns, hphobes, suns, site);

  Draw(site.residues,
       site.bonds,
       placed_suns,
       hphobes,
       lig_names,
       dp);

  //Add unified residues (suns)

  std::vector<UnifiedRes>::iterator sun;
  for (sun = placed_suns.begin(); sun != placed_suns.end(); ++sun) {
    sun->Draw(dp);
  }

  //Add labels
  std::vector<chemlib::Residue*>::iterator res;
  for (res = site.residues.begin(); res != site.residues.end(); ++res) {
    DrawAtomLabels(**res, dp);
    DrawResidueLabel(**res, site, placed_suns, dp);
  }
  //Label(site, hbonds, hphobes);

  dp->Finish();

  return true;
}

void HBondsToLigand(std::vector<chemlib::Residue*>& lig,
                    std::vector<chemlib::Residue*>& rec,
                    std::vector<chemlib::Bond>& bonds,
                    std::vector<HBond>& hbonds,
                    std::vector<chemlib::Residue*>& site) {

  std::vector<chemlib::MIAtom*> lig_donors;         //Atoms in the ligand than can be h-bond donors
  std::vector<chemlib::MIAtom*> lig_acceptors;
  std::vector<chemlib::MIAtom*> lig_polar;              //Union of donors and acceptors

  HBond hbond;

  chemlib::LigandPerceiver lp;
  std::vector<chemlib::Residue*>::iterator lig_res;
  for (lig_res = lig.begin(); lig_res != lig.end(); ++lig_res) {
    lp.AssignImpHydrogens(**lig_res, bonds);
  }

  GetDonors(lig, lig_donors);
  GetAcceptors(lig, lig_acceptors);

  if (lig_donors.empty() && lig_acceptors.empty()) {
    return;
  }

  std::sort(lig_donors.begin(), lig_donors.end());
  std::sort(lig_acceptors.begin(), lig_acceptors.end());        //Create the list of all atoms that
  std::set_union(lig_donors.begin(),                            //are either donors or acceptors or
                 lig_donors.end(),                                           //both.
                 lig_acceptors.begin(),
                 lig_acceptors.end(),
                 std::back_inserter(lig_polar));

  chemlib::Box lig_box(lig, MAX_HBOND_LENGTH);                          //Construct a box encompassing the
  bool retain_res;                                          //area where hydrogen bond partners
  unsigned int i;                                                   //with the ligand may be


  std::vector<chemlib::Residue*>::iterator res;
  std::vector<chemlib::MIAtom*>::iterator ligatm;

  for (res = rec.begin(); res != rec.end(); ++res) {        //Loop over receptor residues
    retain_res = false;

    for (i = 0; i < (*res)->atoms().size(); ++i) {            //Loop over atoms in the residue

      if (!lig_box.Contains(*(*res)->atom(i))) {           //Check if this atom is close enough
        continue;                                           //to the ligand to merit consideration
      }

      if (IsDonor(*(*res)->atom(i), **res)                     //Handle receptor atoms that are both
          && IsAcceptor(*(*res)->atom(i), **res)) {   //donors and acceptors here

        for (ligatm = lig_polar.begin();
             ligatm != lig_polar.end();
             ++ligatm) {

          if (CheckHBond((*res)->atom(i), *ligatm, hbond)) {
            retain_res = true;
            hbonds.push_back(hbond);
          } else if (CheckHBond(*ligatm, (*res)->atom(i), hbond)) {
            retain_res = true;
            hbonds.push_back(hbond);
          }
        }
      } else if (IsDonor(*(*res)->atom(i), **res)) {           //Handle receptor atoms that are only
        for (ligatm = lig_acceptors.begin();                    //donors here
             ligatm != lig_acceptors.end();
             ++ligatm) {
          if (CheckHBond((*res)->atom(i), *ligatm, hbond)) {
            retain_res = true;
            hbonds.push_back(hbond);
          }
        }
      } else if (IsAcceptor(*(*res)->atom(i), **res)) {        //Handle receptor atoms that are only
        for (ligatm = lig_donors.begin();                       //acceptors here
             ligatm != lig_donors.end();
             ++ligatm) {
          if (CheckHBond(*ligatm, (*res)->atom(i), hbond)) {
            retain_res = true;
            hbonds.push_back(hbond);
          }
        }
      }
    }           //Loop over residue atoms

    if (retain_res) {                                       //If we've found any hydrogen bonds
      site.push_back(*res);                                 //to ligand, add this residue to the site
    }

  }         //Loop over residues
}

void HphobesToLigand(std::vector<chemlib::Residue*>& lig,
                     std::vector<chemlib::Residue*>& rec,
                     std::vector<HPhobe>& hphobes,
                     std::vector<chemlib::Residue*>& site) {

  std::vector<chemlib::MIAtom*> lig_nonpolars;          //Atoms in the ligand than can form hydrophobic
  //interactions
  unsigned int i;
  HPhobe hphobe;

  GetNonPolars(lig, lig_nonpolars);

  if (lig_nonpolars.empty()) {
    return;
  }

  chemlib::Box lig_box(lig, MAX_HPHOBE_LENGTH);                     //Construct a box encompassing any
  bool retain_res;                                          //prospective hydrogen bond partners

  std::vector<chemlib::Residue*>::iterator res;
  std::vector<chemlib::MIAtom*>::iterator ligatm;

  for (res = rec.begin(); res != rec.end(); ++res) {        //Loop over receptor residues
    retain_res = false;

    for (i = 0; i < (*res)->atoms().size(); ++i) {                    //Loop over atoms in the residue


      if (!lig_box.Contains(*(*res)->atom(i))) {           //Check if this atom is close enough
        continue;                                           //to the ligand to merit consideration
      }

      if (!IsNonPolar(*(*res)->atom(i))) {                 //Skip polar atoms
        continue;
      }

      for (ligatm = lig_nonpolars.begin();
           ligatm != lig_nonpolars.end();
           ++ligatm) {
        if (CheckHPhobe((*res)->atom(i), *ligatm, hphobe)) {
          hphobe.res1 = *res;
          retain_res = true;
          hphobes.push_back(hphobe);
        }
      }
    }

    if (retain_res) {                                       //If we've found any hphobes
      site.push_back(*res);                                 //to ligand, add this residue to the site
    }
  }
}

#define SQRT_OF_HALF 0.70710678118654752440084436210485F

void DrawAtomLabels(const chemlib::Residue& res,
                    Drawing* dp) {
  std::vector<chemlib::MIAtom*>::const_iterator atm;
  Label* l;
  std::string* label_txt;
  for (atm = res.atoms().begin(); atm != res.atoms().end(); ++atm) {
    label_txt = new std::string((*atm)->name());

    l = new moldraw::Label(dp,
                           (*atm)->x(),
                           (*atm)->y(),
                           Drawing::ATOM_LABEL_OFFSET * Drawing::ATOM_RADIUS,
                           SQRT_OF_HALF,
                           SQRT_OF_HALF,
                           Drawing::ATOM_LABEL_SIZE,
                           *label_txt);
    l->Draw();
    delete l;
    delete label_txt;
  }
}

void DrawResidueLabel(const chemlib::Residue& res,
                      chemlib::Ligand& site,
                      std::vector<UnifiedRes>& suns,
                      Drawing* dp) {

  float pos[2];
  double init_vect[2];

  chemlib::BisectAtom(res.atom(0), init_vect);

  pos[0] = res.atoms().front()->x() + init_vect[0] * 0.8;
  pos[1] = res.atoms().front()->y() + init_vect[1] * 0.8;

  float score = ScoreResLabel(pos, suns, site);
  float trial_score, step;

  int i = 0;
  while (i < Drawing::MAX_RESLABEL_CYCLES && score > 0) {
    step = -0.1F + 0.2F * ((float)rand() / (float)RAND_MAX);
    pos[0] += step;

    if ((trial_score = ScoreResLabel(pos, suns, site)) > score) {
      pos[0] -= step;
    } else {
      score = trial_score;
    }

    pos[1] += step;
    if ((trial_score = ScoreResLabel(pos, suns, site)) > score) {
      pos[1] -= step;
    } else {
      score = trial_score;
    }
    ++i;
  }

  std::string label_text(res.type().c_str());
  label_text += " ";
  label_text += res.name().c_str();
  if (res.chain_id() % 256 != ' ') {
    label_text += "(";
    label_text += res.chain_id() % 256;
    label_text += ")";
  }

  moldraw::Label l(dp,
                   pos[0],
                   pos[1],
                   0.0,
                   1.0,
                   0.0,
                   Drawing::RESIDUE_LABEL_SIZE,
                   label_text);
  l.Draw();
}

/*
  void DrawSuns(std::vector<chemlib::Residue> &suns,
  std::vector<chemlib::Residue> &residues,
  std::vector<HPhobe> &hphobes,
  Drawing *dp) {
  Sun *s;

  double min_angle, max_angle, angle;
  double dir[3];
  std::vector<double> angles;

  std::vector<HPhobe>::iterator hp;

  std::vector<chemlib::Residue>::iterator sun;

  for(sun=suns.begin(); sun!=suns.end(); ++sun){
  for(hp=hphobes.begin(); hp!=hphobes.end(); ++hp) {
  if(hp->res1 == *sun) {
  BondVector(&(*sun)->atom(0), hp->getAtom2(), dir);
  angle = RAD2DEG * atan2(dir[1], dir[0]);
  angles.push_back(angle);
  }
  if(hp->res2 == *sun) {
  BondVector(&(*sun)->atom(0), hp->getAtom1(), dir);
  angle = RAD2DEG * atan2(dir[1], dir[0]);
  angles.push_back(angle);
  }
  }
  if(angles.empty()) {
  continue;
  }

  AngleRange(angles, &min_angle, &max_angle);

  s = new moldraw::Sun(dp,
  sun->atoms.front().x,
  sun->atoms.front().y,
  SUN_RADIUS,
  sin(DEG2RAD * (min_angle + max_angle) / 2.0),
  cos(DEG2RAD * (min_angle + max_angle) / 2.0),
  max_angle - min_angle,
  SUN_BORDER_WIDTH,
  "Test 22");
  s->Draw();
  delete s;
  }
  }
*/


void Draw(std::vector<chemlib::Residue*>& residues,
          std::vector<chemlib::Bond>& bonds,
          const std::vector<UnifiedRes>& suns,
          const std::vector<HPhobe>& hphobes,
          const std::vector<std::string>& lig_names,
          Drawing* dp) {

  Atom* a;
  SingleBond* b;
  HydrogenBond* hb;
  //		float atom_radius = 6.0F;
  //		float residue_radius = 13.5F;
  //		float margin = 1.0F;

  //	dp->AddMargin(1.0);

  chemlib::Box molbox(residues, (float)(5.0f*Drawing::ATOM_RADIUS));

  std::vector<UnifiedRes>::const_iterator sun;
  chemlib::MIAtom atom;
  for (sun = suns.begin(); sun != suns.end(); ++sun) {
    atom.setPosition(sun->x, sun->y, 0.0f);
    molbox.Expand(atom, (float)((Drawing::UNIRES_RADIUS * 1.5F)));
  }


  dp->FitToPage(molbox.Origin(), molbox.Dimensions());
  chemlib::MIAtom* a1, * a2;
  std::vector<chemlib::Residue*>::iterator r1, r2;
  std::vector<std::string>::const_iterator res_name;
  bool is_lig;

  std::vector<chemlib::Bond>::const_iterator bnd;
  for (bnd = bonds.begin();
       bnd != bonds.end();
       ++bnd) {
    is_lig = false;
    a1 = bnd->getAtom1();
    a2 = bnd->getAtom2();
    r1 = ResSearch(a1, residues.begin(), residues.end());
    r2 = ResSearch(a2, residues.begin(), residues.end());
    if (r1 != residues.end() && r1 != residues.end() && r1 == r2) {
      res_name = std::find(lig_names.begin(), lig_names.end(), (*r1)->name().c_str());
      is_lig = (res_name != lig_names.end());
    }


    if (bnd->getOrder() != 5) {
      b = new moldraw::SingleBond(dp,
                                  bnd->getAtom1()->x(),
                                  bnd->getAtom1()->y(),
                                  bnd->getAtom2()->x(),
                                  bnd->getAtom2()->y(),
                                  is_lig ? Drawing::LIG_BOND_WIDTH : Drawing::STD_BOND_WIDTH);
      b->Draw();
      delete b;
    }
  }

  //Kind of dumb to have this in a separate loop?? -KWB 1/19/5
  for (bnd = bonds.begin();
       bnd != bonds.end();
       ++bnd) {
    if (bnd->getOrder() == 5) {
      hb = new moldraw::HydrogenBond(dp,
                                     bnd->getAtom1()->x(),
                                     bnd->getAtom1()->y(),
                                     bnd->getAtom2()->x(),
                                     bnd->getAtom2()->y(),
                                     Drawing::STD_BOND_WIDTH,
                                     Drawing::HBOND_LABEL_SIZE,
                                     bnd->ideal_length);
      hb->Draw();
      delete hb;

      std::string label_text;
      label_text=format("%.2f", bnd->ideal_length);
      //				int height, width;
      //				measure.GetTextExtent(label_text, &height, &width);

      moldraw::Label l(dp,
                       (bnd->getAtom1()->x() + bnd->getAtom2()->x()) / 2.0,
                       (bnd->getAtom1()->y() + bnd->getAtom2()->y()) / 2.0,
                       0.0,
                       1.0,
                       0.0,
                       Drawing::HBOND_LABEL_SIZE,
                       label_text);
      l.Draw();

    }
  }

  moldraw::Spokes* spks;
  std::vector<HPhobe>::const_iterator hp;
  std::set<chemlib::Residue*> done_residues;
  std::set<chemlib::MIAtom*> done_atoms;


  for (hp = hphobes.begin(); hp != hphobes.end(); ++hp) {
    if (done_residues.find(hp->res1) == done_residues.end() ||
        done_atoms.find(hp->getAtom2()) == done_atoms.end()) {
      spks = new moldraw::Spokes(dp,
                                 hp->res1->x(),
                                 hp->res1->y(),
                                 Drawing::UNIRES_RADIUS,
                                 hp->getAtom2()->x(),
                                 hp->getAtom2()->y(),
                                 Drawing::ATOM_RADIUS,
                                 Drawing::SPOKE_EXTENT,
                                 Drawing::ATOM_SPOKE_WIDTH,
                                 Drawing::SUN_BORDER_WIDTH);
      spks->Draw();
      delete spks;
      done_residues.insert(hp->res1);
      done_atoms.insert(hp->getAtom2());
    }
  }

  std::vector<chemlib::Residue*>::iterator ri;
  unsigned int i;
  PaletteColor color;

  for (ri = residues.begin(); ri != residues.end(); ++ri) {
    chemlib::Residue *res=*ri;
    for (i = 0; i < res->atoms().size(); i++) {
      color = GetAtomPlotColor(*res->atom(i));
      a = new moldraw::Atom(dp,
                            res->atom(i)->x(),
                            res->atom(i)->y(),
                            Drawing::ATOM_RADIUS,
                            color);
      a->Draw();
      delete a;
    }
  }
}

void PlaceSuns(std::vector<UnifiedRes>& placed_suns,
               std::vector<HPhobe>& hphobes,
               std::vector<chemlib::Residue*>& suns,
               chemlib::Ligand& site) {
  std::vector<chemlib::Residue*>::iterator sun;
  std::vector<HPhobe>::iterator hp;

  double init_vect[2];
  double pos[2];

  double step;
  double score, trial_score;
  int i;
  std::vector<chemlib::MIAtom*> target_atoms;
  chemlib::MIAtom* start_atom;

  UnifiedRes new_res;

  for (sun = suns.begin(); sun != suns.end(); ++sun) {
    target_atoms.clear();

    //Find a hydrophobic interaction to this residue

    for (hp = hphobes.begin(); hp != hphobes.end(); ++hp) {
      if (hp->res1 == *sun) {
        target_atoms.push_back(hp->getAtom2());
        //					start_atom = hp->getAtom2();
        //					break;
      }
      if (hp->res2 == *sun) {
        target_atoms.push_back(hp->getAtom1());
        //					start_atom = hp->getAtom1();
        //					break;
      }
    }

    if (target_atoms.empty()) {                 //Skip if we couldn't find an interacting atom
      continue;
    }

    start_atom = target_atoms.front();

    chemlib::BisectAtom(start_atom, init_vect);

    pos[0] = start_atom->x() + init_vect[0] + 2.5;
    pos[1] = start_atom->y() + init_vect[1] + 2.5;

    score = ScoreSun(pos, placed_suns, site);

    i = 0;
    while (i < Drawing::MAX_SUN_CYCLES && score > 0) {
      step = -0.2 + 0.4 * (rand() / (float)RAND_MAX);
      pos[0] += step;
      if ((trial_score = ScoreSun(pos, placed_suns, site)) > score) {
        pos[0] -= step;
      } else {
        score = trial_score;
      }

      pos[1] += step;
      if ((trial_score = ScoreSun(pos, placed_suns, site)) > score) {
        pos[1] -= step;
      } else {
        score = trial_score;
      }
      ++i;
    }

    //Place residue
    (*sun)->setPosition(pos[0], pos[1]);
    new_res.type = (*sun)->type().c_str();
    new_res.name = (*sun)->name().c_str();
    new_res.chain_id = (*sun)->chain_id();
    new_res.x = pos[0];
    new_res.y = pos[1];
    new_res.partners = target_atoms;
    placed_suns.push_back(new_res);
  }
}

double ScoreSun(double* pos,
                std::vector<UnifiedRes>& placed_suns,
                chemlib::Ligand& site) {
  double score = 0;
  double dx, dy;

  std::vector<chemlib::Residue*>::iterator res;
  std::vector<chemlib::MIAtom*>::const_iterator atm;

  for (res = site.residues.begin(); res != site.residues.end(); ++res) {
    for (atm = (*res)->atoms().begin(); atm != (*res)->atoms().end(); ++atm) {
      dx = pos[0] - (*atm)->x();
      if (dx > SUN_ATOM_DIST || dx < -SUN_ATOM_DIST) {
        continue;
      }

      dy = pos[1] - (*atm)->y();
      if (dy > SUN_ATOM_DIST || dy < -SUN_ATOM_DIST) {
        continue;
      }

      if (dx*dx + dy*dy < SQ_SUN_ATOM_DIST) {
        score += SQR(SQ_SUN_ATOM_DIST - (dx*dx + dy*dy));
      }
    }
  }

  std::vector<UnifiedRes>::iterator unires;

  for (unires = placed_suns.begin(); unires != placed_suns.end(); ++unires) {
    dx = pos[0] - unires->x;
    if (dx > SUN_SUN_DIST || dx < -SUN_SUN_DIST) {
      continue;
    }

    dy = pos[1] - unires->y;
    if (dy > SUN_SUN_DIST || dy < -SUN_SUN_DIST) {
      continue;
    }

    if (dx*dx + dy*dy < SQ_SUN_SUN_DIST) {
      score += SQR(SQ_SUN_SUN_DIST - (dx*dx + dy*dy));
    }
  }

  return score;
}

float ScoreResLabel(float* pos,
                    std::vector<UnifiedRes>& suns,
                    chemlib::Ligand& site) {
  float score = 0;
  float dx, dy;

  std::vector<chemlib::Residue*>::iterator res;
  std::vector<chemlib::MIAtom*>::const_iterator atm;

  for (res = site.residues.begin(); res != site.residues.end(); ++res) {
    for (atm = (*res)->atoms().begin(); atm != (*res)->atoms().end(); ++atm) {
      dx = pos[0] - (*atm)->x();
      if (dx > RESLABEL_ATOM_DIST || dx < -RESLABEL_ATOM_DIST) {
        continue;
      }

      dy = pos[1] - (*atm)->y();
      if (dy > RESLABEL_ATOM_DIST || dy < -RESLABEL_ATOM_DIST) {
        continue;
      }

      if (dx*dx + dy*dy < SQ_RESLABEL_ATOM_DIST) {
        //					score += 7;
        score += SQR(SQ_RESLABEL_ATOM_DIST - (dx*dx + dy*dy));
      }


    }
  }

  std::vector<UnifiedRes>::iterator unires;

  for (unires = suns.begin(); unires != suns.end(); ++unires) {
    dx = pos[0] - unires->x;
    if (dx > RESLABEL_SUN_DIST || dx < -RESLABEL_SUN_DIST) {
      continue;
    }

    dy = pos[1] - unires->y;
    if (dy > RESLABEL_SUN_DIST || dy < -RESLABEL_SUN_DIST) {
      continue;
    }

    if (dx*dx + dy*dy < SQ_RESLABEL_SUN_DIST) {
      //				score += 10;
      score += SQR(SQ_RESLABEL_SUN_DIST - (dx*dx + dy*dy));
    }
  }

  return score;
}

PaletteColor GetAtomPlotColor(const chemlib::MIAtom& atom) {

  int ci=color_by_name(atom.name());
  if (atom.atomicnumber() == 6) {  //Carbon
    ci=Colors::BLACK;
  }
  if (ci==Colors::WHITE) {
    ci=Colors::GREEN;
  }
  int c = PaletteIndex(ci);
  return PaletteColor(Colors::RPallette[c],Colors::GPallette[c],Colors::BPallette[c]);
}

} //namespace moldraw

