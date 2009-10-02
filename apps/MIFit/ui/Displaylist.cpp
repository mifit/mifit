#include <vector>
#include <algorithm>
#include <cmath>
#include <cstring>

#include <nongui/nonguilib.h>
#include <chemlib/chemlib.h>
#include <chemlib/RESIDUE_.h>
#include <map/maplib.h>
#include "core/corelib.h"
#include "ui/MIDialog.h"

#include "EMap.h"
#include "id.h"
#include "Displaylist.h"
#include "molw.h"

using namespace std;
using namespace chemlib;

Displaylist::DisplaylistRefCountMap Displaylist::refCounts;

bool Displaylist::isValid(Displaylist* list) {
  bool result = false;
  if (list != NULL && refCounts.find(list) != refCounts.end()) {
    result = refCounts[list] > 0;
  }
  return result;
}


Displaylist::Displaylist() { // create a NULL list
  ++refCounts[this];
  PickedAtom = NULL;
  PickedResidue = NULL;
  PickedMolecule = NULL;
  current = NULL;
  m_currentmap = NULL;
  delete_level = 0;
}

Displaylist::~Displaylist() {

  // delete all items in list
  std::list<Molecule*>::iterator node;
  for (node = Models.begin(); node != Models.end(); node++) {
    Molecule* m = *node;
    m->disconnect(this);
    delete m;
  }
  for (unsigned int i = 0; i < Maps.size(); i++) {
    delete Maps[i];
  }
  --refCounts[this];
  if (refCounts[this] == 0) {
    refCounts.erase(this);
  }
}

Molecule* Displaylist::operator [](int elem) {
  if (elem < 0 || elem >= (int) Models.size()) {
    return NULL;
  }
  std::list<Molecule*>::iterator mp = Models.begin();
  int m = 0;
  while (m < elem) {
    mp++; m++;
  }
  return *mp;
}

int Displaylist::AddItem(RESIDUE* reslist, std::string cmpd, FILE* fp, std::vector<Bond> * connects, int type) {
  Bond* edges = NULL;
  if (connects->size() > 0) {
    edges = &((*connects)[0]);
  }
  Molecule* mol = new Molecule(reslist, cmpd, fp, edges, connects->size(), type);
  return AddItem(mol);
}


void Displaylist::symmetryToBeCleared(MIMoleculeBase* mol) {
  // treat as if all symmetry residues are being deleted
  std::vector<RESIDUE*> deaders;
  for (MIIter<RESIDUE> res=mol->GetSymmResidues(); Residue::isValid(res); ++res) {
    deaders.push_back(res);
  }
  if (deaders.size()) {
    delete_level++;
    residuesToBeDeleted(mol, deaders);
    delete_level--;
  }
}

void Displaylist::moleculeToBeDeleted(MIMoleculeBase* mol) {
  if (!mol)  //is this even possible?
    return;

 std::vector<RESIDUE*> residues;
  for (MIIter<RESIDUE> res=mol->GetResidues(); res; ++res) {
    residues.push_back(res);
  }
  for (MIIter<RESIDUE> res=mol->GetSymmResidues(); Residue::isValid(res); ++res) {
    residues.push_back(res);
  }
  delete_level++;
  residuesToBeDeleted(mol,residues);
  delete_level--;

  //remove the molecule from the list, and reset current, if necessary
  //note: this is what the old DeleteItem function used to do
  std::list<Molecule*>::iterator node= std::find(Models.begin(), Models.end(), mol);
  if (node == Models.end()) {
    return; 
  }
  if (current == mol) {
    SetCurrent(NULL);
  }
  Models.remove((Molecule*)mol);
  SetCurrent(FirstItem());

  SetPicked(NULL, NULL, NULL);
}

void Displaylist::residuesToBeDeleted(MIMoleculeBase *mol, std::vector<RESIDUE*> &residues) {

  MIAtomList atoms;
  for (size_t i=0; i < residues.size(); ++i) {
    atoms.insert(atoms.end(), residues[i]->atoms().begin(), residues[i]->atoms().end());
    if (PickedResidue == residues[i])
      PickedResidue=NULL;
  }
  delete_level++;
  atomsToBeDeleted(mol,atoms);
  delete_level--;
  
  if (delete_level==0)
    SetPicked(NULL, NULL, NULL);
}


static std::set<MIAtom*> *ATOM_SET=0;
bool Contains(const CONTACT &c)
{
  MIAtom* a1=c.line.getAtom1();
  MIAtom* a2=c.line.getAtom2();
  return (std::find(ATOM_SET->begin(), ATOM_SET->end(), a1)!=ATOM_SET->end() ||
          std::find(ATOM_SET->begin(), ATOM_SET->end(), a2)!=ATOM_SET->end());
}

void Displaylist::atomsToBeDeleted(MIMoleculeBase *, const MIAtomList &atoms) {

  //Prune contacts list:
  std::set<MIAtom*> atom_set;
  atom_set.insert(atoms.begin(), atoms.end());
  ATOM_SET=&atom_set;
  std::vector<CONTACT>::iterator new_end=
    std::remove_if(Contacts.begin(), Contacts.end(),Contains);
  Contacts.erase(new_end, Contacts.end());


  if (delete_level==0)
    SetPicked(NULL, NULL, NULL);
}

static void AddAnnotations(Molecule *model)
{
  if (model == NULL) {
    return;
  }

  std::map<RESIDUE*,std::string> annotations;
  std::string err_type;
  std::vector<std::string> headers=*model->GetFileHead();

  for (size_t i=0; i < headers.size(); ++i)
  {
    // get class of error
    if (strncmp("REMARK 500 SUBTOPIC:",headers[i].c_str(),20)==0 ||
        strncmp("REMARK 501 SUBTOPIC:",headers[i].c_str(),20)==0) {
      err_type=""; // clear type on new subtopic
      if (headers[i].size() > 22) {
        if (headers[i].find("BOND LENGTHS",21)!=std::string::npos) {
          err_type=std::string("Geometry (bond length)");
        } else if (headers[i].find("BOND ANGLES",21)!=std::string::npos) {
          err_type=std::string("Geometry (bond angle)");
        } else if (headers[i].find("CHIRAL",21)!=std::string::npos) {
          err_type=std::string("Geometry (chiral)");
        } else if (headers[i].find("NON-CIS",21)!=std::string::npos) {
          err_type=std::string("Omega");
        } else if (headers[i].find("CLOSE CONTACTS",21)!=std::string::npos) {
          err_type=std::string("Van der Waals");
        } else if (headers[i].find("TORSION ANGLES",21)!=std::string::npos) {
          err_type=std::string("Phi-psi");
        } else if (headers[i].find("ELECTRON DENSITY",21)!=std::string::npos) {
          err_type=std::string("Density");
        } else if (headers[i].find("ELECTRON DENSITY",21)!=std::string::npos) {
          err_type=std::string("Density");
        } else if (headers[i].find("ROTAMER",21)!=std::string::npos) {
          err_type=std::string("Rotamer chi-1");
        } else if (headers[i].find("CIS PEPTIDE",21)!=std::string::npos) {
          err_type=std::string("Cis peptide");
        }
      }
    }

    int dum;
    char rname[4],rnum[6];
    char chain;
    if (err_type != "" &&
        (strncmp("REMARK 500   ",headers[i].c_str(),13)==0 ||
         strncmp("REMARK 501   ",headers[i].c_str(),13)==0)) {
      if ((sscanf(headers[i].c_str(),
                  "REMARK 500   %d %s %c %s", &dum,rname,&chain,rnum)==4) ||
          (sscanf(headers[i].c_str(),
                  "REMARK 501   %d %s %c %s", &dum,rname,&chain,rnum)==4))
      {
        RESIDUE* res = residue_from_name(model->getResidues(), rnum, chain);
        if (res != NULL) {
          std::string tmp=annotations[res];
          if (tmp.size()==0)
            tmp=::format("Error in %c %s %s: %s", chain, rnum, rname, err_type.c_str());
          else
            tmp=::format("%s, %s",tmp.c_str(),err_type.c_str());
          annotations[res]=tmp;
        }
      }
    }
  }

  bool hidden=(MIConfig::Instance()->GetProfileInt("DisplayView", "autoShowError", 1) == 0);


  for (std::map<RESIDUE*,std::string>::iterator i=annotations.begin(); i!=annotations.end(); ++i) {
    RESIDUE *res=i->first;
    std::string &text=i->second;
    MIAtom* atom = atom_from_name("CA", *res);
    if (atom == NULL && res->atomCount() > 0) {
      atom = res->atom(0);
    }
    if (atom != NULL) {
      Annotation* annotation = new Annotation(text.c_str(), atom->x(), atom->y(), atom->z());
      annotation->setHidden(hidden);
      model->addAnnotation(annotation);
    }
  }
}


int Displaylist::AddItem(Molecule* node) { //  add a new node in displaylist
  if (node == NULL) {
    Logger::message("Additem: Unable to allocate new object");
    return (0);
  } else {

    connect(node,
            SIGNAL(atomsToBeDeleted(chemlib::MIMoleculeBase*,chemlib::MIAtomList)),
            this, SLOT(atomsToBeDeleted(chemlib::MIMoleculeBase*,chemlib::MIAtomList)));
    connect(node, SIGNAL(residuesToBeDeleted(chemlib::MIMoleculeBase*,std::vector<chemlib::RESIDUE*>&)),
            this, SLOT(residuesToBeDeleted(chemlib::MIMoleculeBase*,std::vector<chemlib::RESIDUE*>&)));
    connect(node, SIGNAL(moleculeToBeDeleted(chemlib::MIMoleculeBase*)),
            this, SLOT(moleculeToBeDeleted(chemlib::MIMoleculeBase*)));
    connect(node, SIGNAL(symmetryToBeCleared(chemlib::MIMoleculeBase*)),
            this, SLOT(symmetryToBeCleared(chemlib::MIMoleculeBase*)));

    Models.push_back(node);
    modelAdded(node);
    SetCurrent(node);
    AddAnnotations(node);
    return (Models.size());
  }
  return 0;
}


void Displaylist::LabelPick(bool toggle) {
  std::list<Molecule*>::iterator node = Models.begin();
  while (node != Models.end()) {
    Molecule* molecule = *node;
    if (molecule == PickedMolecule) {
      if (toggle) {
        if (molecule->isAtomLabeled(PickedAtom)) {
          molecule->unlabelAtom(PickedAtom);
        } else {
          molecule->labelAtom(PickedAtom, PickedResidue);
        }
      }
      break;
    }
    node++;
  }
}

void Displaylist::ClearLabels() {
  std::list<Molecule*>::iterator node = Models.begin();
  while (node != Models.end()) {
    (*node)->clearAtomLabels();
    node++;
  }
}

/*
   Molecule * Displaylist::LastItem()
   {
    std::list<Molecule*>::iterator node = Models.begin();
    if(!node) return(NULL);
    while( node->getnext() != NULL) node++;
        return(node);
   }*/

Molecule* Displaylist::SetCurrent(Molecule* mol) {
  if (current == mol) {
    return current;
  }
  Molecule* oldCurrent = current;
  if (mol == NULL) {
    Logger::debug("current empty");
    current = NULL;
  } else {
    // look for mol in list and set to current if found;
    std::list<Molecule*>::iterator node = std::find(Models.begin(), Models.end(), mol);
    if (node != Models.end()) {
      current = mol;
      Logger::debug(mol->compound.c_str());
    } else {
      Logger::log("Set current failed");
    }
  }
  currentMoleculeChanged(oldCurrent, current);
  return current;
}

int Displaylist::FindContacts(MIAtom* a, RESIDUE* r, float cutoff, ViewPoint* vp) {
  int n = 0, i;
  long lcut = (long)(cutoff*CRS);
  float d;
  MIAtom* a2;
  std::list<Molecule*>::iterator node = Models.begin();
  while (node != Models.end()) {
    for (MIIter<RESIDUE> res = (*node)->GetResidues(); res; ++res) {
      if (r != res) {
        for (i = 0; i < res->atomCount(); i++) {
          a2 = res->atom(i);
          if (IsOnScreen(a2, vp)) {
            if (abs((int)(a->x() - a2->x())) < lcut) {
              if (abs((int)(a->z() - a2->z())) < lcut) {
                if (abs((int)(a->y() - a2->y())) < lcut) {
                  if ((d = (float)AtomDist(*a, *a2)) <= cutoff) {
                    if (AddContact(a, a2, d)) {
                      n++;
                    }
                  }
                }
              }
            }
          }
        }
      }
    }
    node++;
  }
  return n;
}

bool Displaylist::AddContact(MIAtom* a1, MIAtom* a2, float d) {
  CONTACT Contact;
  Contact.line.setAtom1(a1);
  Contact.line.setAtom2(a2);
  Contact.color = Colors::GREEN;
  Contact.d = d;
  Contacts.push_back(Contact);
  return true;
}

bool Displaylist::AddLine(float x1, float y1, float z1, float x2, float y2, float z2, int color) {
  PLINE Vu;
  Vu.p1.x = x1;
  Vu.p1.y = y1;
  Vu.p1.z = z1;
  Vu.p2.x = x2;
  Vu.p2.y = y2;
  Vu.p2.z = z2;
  Vu.p1.color = color;
  Vu.p2.color = color;
  Vus.push_back(Vu);
  return true;
}

bool Displaylist::SetCurrentMap(EMap* emap) {
  EMap* oldMap = m_currentmap;
  /* map number is indexed from 1 */
  if (Maps.size() == 0) {
    return false;
  }
  for (unsigned int i = 0; i < Maps.size() ; i++) {
    if (Maps[i] == emap) {
      m_currentmap = Maps[i];
      currentMapChanged(oldMap, m_currentmap);
      return true;
    }
  }
  return false;
}

void Displaylist::AddMap(EMap* emap) {
  Maps.push_back(emap);
  m_currentmap = emap;
  mapAdded(emap);
}

void Displaylist::DeleteMap(EMap* emap) {
  MapList::iterator map = find(Maps.begin(), Maps.end(), emap);
  mapToBeDeleted(emap);
  Maps.erase(map);
  if (m_currentmap == emap) {
    if (Maps.size() > 0) {
      m_currentmap = *Maps.begin();
    } else {
      m_currentmap = NULL;
    }
  }
  delete emap;
}

void Displaylist::ChooseActiveMap() {
  std::vector<std::string> choices;
  for (int i = 0; i < MapCount(); i++) {
    if (i >= 100) {
      break;
    }
    EMap* emap = Maps[i];
    std::string choice = ftoa(i);
    choice += ": ";
    choice += emap->MapID().c_str();
    choices.push_back(choice);
  }
  int selection = MIGetSingleChoiceIndex("Choose Map to make active", "MIFit", choices);
  if (selection == -1) {
    return;
  }
  if (selection < MapCount()) {
    SetCurrentMap(Maps[selection]);
  }
}

//DEL void Displaylist::CheckCenter(float x, float y, float z, bool symm_link)
//DEL {
//DEL 	for(unsigned int i=0; i<Maps.size(); i++) Maps[i]->CheckCenter(x, y, z);
//DEL 	std::list<Molecule*>::iterator node;
//DEL 	for(node = Models.begin(); node != Models.end(); node++){
//DEL 		(*node)->CheckCenter(x,y,z, symm_link);
//DEL 	}
//DEL }

long Displaylist::SurfaceCurrent(MIAtomList atoms, float radius_mult) {
  unsigned int i;
  SURFDOT* adots;
  long nadots = 0;
  long maxadots = 0;
  void* hdots = NULL;
  float srfdotsper = 0.20F;
  int l;
  nadots = 0;

  std::vector<SURFDOT>().swap(CurrentDots); // was CurrentDots.clear();
  for (i = 0; i < atoms.size(); i++) {
    nadots = atomsurf(atoms[i], radius_mult, &atoms[0], atoms.size(), &adots, nadots, &maxadots, hdots, srfdotsper);
  }
  CurrentDots.reserve(nadots);
  for (l = 0; l < nadots; l++) {
    CurrentDots.push_back(adots[l]);
  }
  if (hdots != NULL) {
    free(hdots);
    hdots = NULL;
  }
  return CurrentDots.size();
}

void Displaylist::ClearCurrentSurface() {
  std::vector<SURFDOT>().swap(CurrentDots); // was CurrentDots.clear();
}

void Displaylist::ProbeSurface(RESIDUE* reslist, MIAtomList a) {
  vector<SURFDOT>::iterator p;
  MIAtom_iter pa;
  RESIDUE* res = reslist;
  float r1, dx, dy, dz, dr;
  int i;
  MIAtomList b;
  MIAtom* atom;

  while (res != NULL) {
    for (i = 0; i < res->atomCount(); i++) {
      atom = res->atom(i);
      for (pa = a.begin(); pa != a.end(); pa++) {
        if (atom == *pa) {
          continue;
        }
      }
      r1 = atom->getRadius();
      for (pa = a.begin(); pa != a.end(); pa++) {
        if (AtomDist(**pa, *atom) < (r1 + (*pa)->getRadius())) {
          b.push_back(atom);
          break;
        }
      }
    }
    res = res->next();
  }
  for (p = CurrentDots.begin(); p != CurrentDots.end(); p++) {
    dr = 0.0;
    for (pa = b.begin(); pa != b.end(); pa++) {
      r1 = (*pa)->getRadius();
      r1 = r1*r1;
      dx = (*pa)->x() - (*p).x;
      dx = dx*dx;
      if (dx < r1) {
        dy = (*pa)->y() - (*p).y;
        dy = dy*dy;
        if (dx+dy < r1) {
          dz = (*pa)->z() - (*p).z;
          dz = dz*dz;
          dr = (dx+dy+dz*dz)-r1;
        }
      }
    }
    if (dr < -0.4) {
      (*p).color = Colors::WHITE;
      if (dr < -1.2) {
        (*p).color = Colors::PINK;
        if (dr < -2.4) {
          (*p).color = Colors::RED;
        }
      }
    } else {
      (*p).color = -abs((*p).color);
    }
  }
}

void Displaylist::UpdateContacts() {
  for (size_t i = 0; i < Contacts.size(); i++) {
    Contacts[i].d = AtomDist(*Contacts[i].line.getAtom1(), *Contacts[i].line.getAtom2());
  }
}

bool Displaylist::IsModified() {
  bool changed = false;
  std::list<Molecule*>::iterator node = Models.begin();
  while (node != Models.end()) {
    changed |= (*node)->GetModified();
    node++;
  }
  return changed;
}

void Displaylist::SetModified(bool value) {
  std::list<Molecule*>::iterator node;
  for (node = Models.begin(); node != Models.end(); ++node) {
    (*node)->SetModified(value);
  }  
}

APOINT Displaylist::GetVuCenter() {
  APOINT p;
  p.x = 0;
  p.y = 0;
  p.z = 0;
  if (Vus.size() > 0) {
    for (size_t i = 0; i < Vus.size(); i++) {
      p.x += Vus[i].p1.x;
      p.x += Vus[i].p2.x;
      p.y += Vus[i].p1.y;
      p.y += Vus[i].p2.y;
      p.z += Vus[i].p1.z;
      p.z += Vus[i].p2.z;
    }
    p.x /= ((float)Vus.size()*2.0F);
    p.y /= ((float)Vus.size()*2.0F);
    p.z /= ((float)Vus.size()*2.0F);
  }
  return p;
}

void Displaylist::SetPicked(Molecule* mol, RESIDUE* res, MIAtom* atom) {
  if (PickedMolecule == mol && PickedResidue == res && PickedAtom == atom) {
    return;
  }
  PickedMolecule = mol;
  PickedResidue = res;
  PickedAtom = atom;
  selectionChanged(mol, res, atom);
}

