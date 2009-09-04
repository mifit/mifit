#ifndef MIFIT_MODEL_DISPLAYLIST_H_
#define MIFIT_MODEL_DISPLAYLIST_H_

class Displaylist;

#include <cmath>
#include <boost/signal.hpp>

#include "core/corelib.h"
#include "macafxwin.h"
#include "Xguicryst.h"

class EMap;

/**
 * List of objects to be drawn on the canvas.
 */
class Displaylist {

  typedef std::map<Displaylist*, size_t> DisplaylistRefCountMap;

  /**
   * Stores reference counts to objects of this class. Currently,
   * only increments and decrements in constructor and destructor.
   * Used to check if object has been deleted with isValid method.
   */
  static DisplaylistRefCountMap refCounts;

  Displaylist(const Displaylist& /* list */) {}
  Displaylist& operator=(const Displaylist& /* list */) { return *((Displaylist*)0); } // NOTE: broken implementation just to avoid compiler warning, do not use!

  std::vector<boost::signals::connection> connections;
  
public:

  /**
   * Returns whether the given displaylist is still valid (not been deleted).
   */
  static bool isValid(Displaylist* list);

  typedef std::list<Molecule*> ModelList;
  typedef std::vector<EMap*> MapList;

private:
  EMap* m_currentmap;
  ModelList Models;
  Molecule* current;
  chemlib::MIAtom* PickedAtom;
  chemlib::RESIDUE* PickedResidue;
  Molecule* PickedMolecule;
  std::vector<CONTACT> Contacts;

  MapList Maps;
  std::vector<PLINE> Vus;
  std::vector<SURFDOT> CurrentDots;

  // 
  unsigned int delete_level;

  

public:

  Displaylist();
  ~Displaylist();

  APOINT GetVuCenter();
  size_t VuSize();
  bool IsModified();
  void SetModified(bool value);

  /**
   * Returns the annotation closest to the mousepoint, x, y.
   */
  void UpdateContacts();


  /**
   * Makes a dot-surface representing the contacts around an atom list.
   * @param reslist the residues to be checked for contats with the atom list
   * @param atoms a vector of atoms around which contacts are to be found
   */
  void ProbeSurface(chemlib::RESIDUE* reslist, std::vector<chemlib::MIAtom*> atoms);

  /**
   * Clear the current surface.
   */
  void ClearCurrentSurface();

  /**
   * Surface the current atoms
   * @param atoms vector of atoms to be surfaced.
   * @param radius_mult radius multipier for creating an extended surface
   */
  long SurfaceCurrent(std::vector<chemlib::MIAtom*> atoms, float radius_mult);

  /**
   * Displays a list of maps for the user to choose from.
   */
  void ChooseActiveMap();

  /**
   * Returns the number of maps in the display list.
   */
  int MapCount();

  /**
   * Add an EMap to the display list.
   * @param map a pointer to the map.  Do not delete the map use DeleteMap.
   */
  void AddMap(EMap* map);

  /**
   * Remove a map from the list and destroy it.
   * @param map a pointer to the map.
   */
  void DeleteMap(EMap* map);

  /**
   * Return a pointer to the current EMap.
   */
  EMap* GetCurrentMap();

  /**
   * Return the map at position i
   * @param i zero-based index of the map in the list
   */
  EMap* GetMap(int i);

  /**
   * Make an EMap current.
   * @param map a pointer to the EMap to make current.
   */
  bool SetCurrentMap(EMap* map);

  /**
   * Return the number of models in the list.
   */
  int NumberItems();

  /**
   * Add a model to the list
   */
  int AddItem(Molecule*);

  /**
   * Add a new model from a residue list.
   */
  int AddItem(chemlib::RESIDUE *, std::string, FILE *, std::vector<chemlib::Bond> *, int = MoleculeType::PDB);

  /**
   * Delete a model from the list.
   */
  //int DeleteItem(Molecule*);

  /**
   * Return the current model.
   */
  Molecule* CurrentItem();

  /**
   * Return the current model. Same as CurrentItem.
   */
  Molecule* GetCurrentModel();

  /**
   * Return the first model in the list.  Returns null if the list is empty.
   */
  Molecule* FirstItem();

  /**
   * Make a model the current one.
   */
  Molecule* SetCurrent(Molecule* node);

  /**
   * Get the last item in the list.  Returns null if the list is empty.
   */
  Molecule* LastItem();

  /**
   * Return a node iterator to the beginning of the list
   */
  std::list<Molecule*>::iterator begin();

  /**
   * Return a node iterator to the end of the list
   */
  std::list<Molecule*>::iterator end();


  /**
   * Return a pointer to the atom picked.
   */
  chemlib::MIAtom* GetPickedAtom();

  void SetPicked(Molecule* mol, chemlib::RESIDUE* res, chemlib::MIAtom* atom);

  /**
   * Return a pointer to the residue picked.
   */
  chemlib::RESIDUE* GetPickedResidue();

  /**
   * Return a pointer to the model picked.
   */
  Molecule* GetPickedMolecule();

  /**
   * Pick the lable closest to a mouse click.
   */
  void LabelPick(bool toggle = 1);

  /**
   * Clear the labels list.
   */
  void ClearLabels();

  /**
   * Find the contacts to an atom in the context of a residue list.
   * Contacts are lines between adjacent atoms with a distance label.
   */
  int FindContacts(chemlib::MIAtom *, chemlib::RESIDUE *res, float, ViewPoint *vp);

  /**
   * Add a contact between two atoms.
   * Contacts are lines between adjacent atoms with a distance label.
   */
  bool AddContact(chemlib::MIAtom *, chemlib::MIAtom *, float);

  /**
   * Clear the contact list.
   */
  void ClearContacts();

  /**
   * Add a line to the Vu object
   */
  bool AddLine(float x1, float y1, float z1, float x2, float y2, float z2, int color);

  /**
   * Clear the vu list.
   */
  void ClearVus();

  /**
   * Index operator to access the model list.
   */
  Molecule* operator [](int elem);

  ModelList& getMolecules();
  std::vector<CONTACT>& getContacts();
  MapList& getMaps();
  std::vector<PLINE>& getLines();
  std::vector<SURFDOT>& getCurrentDots();

  boost::signal1<void, Molecule*> modelAdded;
  boost::signal2<void, Molecule* /*old*/, Molecule* /*new*/> currentMoleculeChanged;
  boost::signal1<void, EMap*> mapAdded;
  boost::signal1<void, EMap*> mapToBeDeleted;
  boost::signal2<void, EMap* /*old*/, EMap* /*new*/> currentMapChanged;
  boost::signal3<void, Molecule*, chemlib::RESIDUE*, chemlib::MIAtom*> selectionChanged;


  void atomsToBeDeleted(chemlib::MIMoleculeBase* model, const std::vector<chemlib::MIAtom*> &atoms);
  void residuesToBeDeleted(chemlib::MIMoleculeBase* model, std::vector<chemlib::RESIDUE*> &res);
  void moleculeToBeDeleted(chemlib::MIMoleculeBase *model);
  void symmetryToBeCleared(chemlib::MIMoleculeBase *model);
};

inline std::list<Molecule*>& Displaylist::getMolecules() {
  return Models;
}

inline std::vector<CONTACT>& Displaylist::getContacts() {
  return Contacts;
}

inline Displaylist::MapList& Displaylist::getMaps() {
  return Maps;
}

inline std::vector<PLINE>& Displaylist::getLines() {
  return Vus;
}

inline std::vector<SURFDOT>& Displaylist::getCurrentDots() {
  return CurrentDots;
}

inline size_t Displaylist::VuSize() {
  return Vus.size();
}

inline int Displaylist::MapCount() {
  return Maps.size();
}

inline EMap* Displaylist::GetCurrentMap() {
  return m_currentmap;
}

inline EMap* Displaylist::GetMap(int i) {
  return (EMap*)Maps[i];
}

inline int Displaylist::NumberItems() {
  return Models.size();
}

inline Molecule* Displaylist::CurrentItem() {
  return (current);
}

inline Molecule* Displaylist::GetCurrentModel() {
  return (current);
}

inline Molecule* Displaylist::FirstItem() {
  if (Models.size() > 0) {
    return (Models.front());
  } else {
    return NULL;
  }
}

inline Molecule* Displaylist::LastItem() {
  if (Models.size() > 0) {
    return Models.back();
  } else {
    return NULL;
  }
}

inline std::list<Molecule*>::iterator Displaylist::begin() {
  return Models.begin();
}

inline std::list<Molecule*>::iterator Displaylist::end() {
  return Models.end();
}

inline chemlib::MIAtom* Displaylist::GetPickedAtom() {
  return PickedAtom;
}

inline chemlib::RESIDUE* Displaylist::GetPickedResidue() {
  return PickedResidue;
}

inline Molecule* Displaylist::GetPickedMolecule() {
  return PickedMolecule;
}

inline void Displaylist::ClearContacts() {
  std::vector<CONTACT>().swap(Contacts); // was Contacts.clear();
}

inline void Displaylist::ClearVus() {
  std::vector<PLINE>().swap(Vus); // was Vus.clear();
}

#endif /*MIFIT_MODEL_DISPLAYLIST_H_*/
