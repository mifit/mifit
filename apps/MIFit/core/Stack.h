#ifndef mifit_model_Stack_h
#define mifit_model_Stack_h

class Stack;

#include <stack>
#include <deque>

#include "CRect.h"
#include "Molecule.h"

typedef struct {
  chemlib::MIAtom* atom;
  chemlib::RESIDUE* residue;
  Molecule* molecule;
} StackItem;


/**
 * Implements an atom stack for storing atom picks.
 */
class Stack {
public:
  typedef std::deque<StackItem> DataContainer;

private:
  DataContainer data;

  bool minimized;
  bool changed;
  CRect ClearBox;
  CRect HideBox;
  CRect PopBox;

  void dataPop();
  void Purge(chemlib::MIMoleculeBase* model);
  void Purge(chemlib::RESIDUE* res);
  void Purge(chemlib::MIAtom* atom);

public:

  Stack();
  ~Stack();




  int size();
  bool empty();
  StackItem top();
  const DataContainer& getData();
  void Pop(chemlib::MIAtom*& atom, chemlib::RESIDUE*& res);
  void Pop(chemlib::MIAtom*& atom, chemlib::RESIDUE*& res, Molecule*& m);
  chemlib::MIAtom* Pop();
  void Peek(chemlib::MIAtom*& atom, chemlib::RESIDUE*& res, Molecule*& m);
  void Push(chemlib::MIAtom* atom, chemlib::RESIDUE* res, Molecule* m);

  void Clear();
  void ExpandTopAllAtoms();
  void ExpandTop2AllAtoms();
  void ExpandTop2Range();
  bool InStack(chemlib::RESIDUE* res);

  bool StackChanged();
  void ClearChanged();

  CRect& getClearBox();
  CRect& getHideBox();
  CRect& getPopBox();
  bool PickClearBox(int sx, int sy);
  bool PickHideBox(int sx, int sy);
  bool PickPopBox(int sx, int sy);
  void Minimize();
  void Maximize();
  bool isMinimized();
  void ToggleMinMax();

  // note: in practice these slots are called by CMolwView, not directly by the signal
  void residuesToBeDeleted(chemlib::MIMoleculeBase* model, std::vector<chemlib::RESIDUE*>& residues);
  void atomsToBeDeleted(chemlib::MIMoleculeBase* model, const std::vector<chemlib::MIAtom*>& atoms);
  void moleculeToBeDeleted(chemlib::MIMoleculeBase* model);

  boost::signal1<void, bool> emptyChanged;
};

inline bool Stack::isMinimized() {
  return minimized;
}

inline int Stack::size() {
  return data.size();
}

inline bool Stack::empty() {
  return data.empty();
}

inline StackItem Stack::top() {
  return data.back();
}

inline const Stack::DataContainer& Stack::getData() {
  return data;
}

inline CRect& Stack::getClearBox() {
  return ClearBox;
}

inline CRect& Stack::getHideBox() {
  return HideBox;
}

inline CRect& Stack::getPopBox() {
  return PopBox;
}

inline void Stack::Clear() {
  DataContainer().swap(data); // was data.clear();
}

inline bool Stack::PickClearBox(int sx, int sy) {
  return ClearBox.Within(sx, sy);
}

inline bool Stack::PickHideBox(int sx, int sy) {
  return HideBox.Within(sx, sy);
}

inline bool Stack::PickPopBox(int sx, int sy) {
  return PopBox.Within(sx, sy);
}

inline void Stack::Minimize() {
  minimized = true;
}

inline void Stack::Maximize() {
  minimized = false;
}

inline void Stack::ToggleMinMax() {
  minimized = !minimized;
}

#endif

