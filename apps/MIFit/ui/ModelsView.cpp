#include <QTreeWidget>
#include <QApplication>
#include <QToolButton>
#include <QLineEdit>
#include <QVBoxLayout>
#include <QHBoxLayout>
#include <QResizeEvent>
#include <QSplitter>
#include <QPushButton>
#include <QSettings>
#include <QActionGroup>
#include <QStackedLayout>

#include <set>

#include <map/maplib.h>
#include "core/corelib.h"
#include <nongui/nonguilib.h>
#include <chemlib/RESIDUE_.h>

#include "ModelsView.h"
#include "EMap.h"
#include "id.h"
#include "Displaylist.h"
#include "TreeData.h"
#include "Application.h"

#include "MIMenu.h"
#include "MIGLWidget.h"
#include "MIEventHandler.h"
#include "MIMainWindow.h"
#include "MIQTreeWidget.h"

#include "ui/MIDialog.h"

#include "surf.h"

#include <images/model.xpm>
#include <images/modelSelected.xpm>
#include <images/map.xpm>
#include <images/mapSelected.xpm>
#include <images/crystal.xpm>
#include <images/chain.xpm>
#include <images/chainSelected.xpm>
#include <images/chainHidden.xpm>
#include <images/chainHiddenSelected.xpm>
#include <images/chainPartial.xpm>
#include <images/chainPartialSelected.xpm>
#include <images/residue.xpm>
#include <images/residueSelected.xpm>
#include <images/residueHidden.xpm>
#include <images/residueHiddenSelected.xpm>
#include <images/residuePartial.xpm>
#include <images/residuePartialSelected.xpm>
#include <images/atom.xpm>
#include <images/atomSelected.xpm>
#include <images/atomHidden.xpm>
#include <images/atomHiddenSelected.xpm>
#include <images/synced.xpm>

#include "MIEventHandlerMacros.h"

const int ID_MODELSVIEW_OFFSET = 10000;
const int ID_MODELSVIEW_MODELSTREE_SAVETOCRYSTALS = ID_MODELSVIEW_OFFSET+1;
const int ID_MODELSVIEW_SYNCEDTOOL = ID_MODELSVIEW_OFFSET+2;
const int ID_MODELSVIEW_RESIDUESTREE_INSERT = ID_MODELSVIEW_OFFSET+3;
const int ID_MODELSVIEW_MODELSTREE_COPY = ID_MODELSVIEW_OFFSET+4;
const int ID_MODELSVIEW_MODELSTREE_INSERT = ID_MODELSVIEW_OFFSET+5;
const int ID_MODELSVIEW_RESIDUESTREE_PASTE = ID_MODELSVIEW_OFFSET+6;

using namespace chemlib;

bool syncView;
typedef std::map<const char*, int> ImageIndexMap;

class AtomsTree : public MIQTreeWidget, public MIEventHandler {
  Q_OBJECT

  QTreeWidgetItem* rootId;
  MIGLWidget* view;
  Molecule* model;
  RESIDUE* residue;
  MIAtom* currentAtom;

  ImageIndexMap imageIndex;
  typedef std::map<MIAtom*, TreeData*> AtomToDataMap;
  AtomToDataMap atomToData;

  MIMenu* _menu;
  QMenu *solidSurfMenu;
  QActionGroup* solidSurfActionGroup;
  bool _working;

  virtual void contextMenuEvent(QContextMenuEvent* event);

public:

  AtomsTree(QWidget* parent);
  virtual ~AtomsTree();

  void setView(MIGLWidget* view);
  void setModel(Molecule* model);
  void setResidue(RESIDUE* residue);
  void setCurrentAtom(MIAtom* atom);
  MIAtom* getCurrentAtom();
  void stylizeAtom(MIAtom* atom);

private slots:
  void atomChanged(chemlib::MIMoleculeBase* model, chemlib::MIAtomList& atom);
  void atomsToBeDeleted(chemlib::MIMoleculeBase* model, const chemlib::MIAtomList& atoms);

  void OnItemClicked(QTreeWidgetItem *item, int column); // single click
  void OnItemActivated(QTreeWidgetItem *item, int column); // double click

  void ShowItem();
  void ColorItem();
  void DeleteItem();
  void EditItem();
  void solidSurfaceActionTriggered(QAction* action);
  void updateSolidSurfMenu();

};

AtomsTree::AtomsTree(QWidget* parent)
: MIQTreeWidget(parent),MIEventHandler(this),
  view(NULL), model(NULL), residue(NULL) {

  setSelectionMode(QAbstractItemView::ExtendedSelection);

  std::vector<QIcon> imageList;
  imageIndex["atom"] = imageList.size();
  QIcon atomImage=QIcon(QPixmap(atom_xpm));
  imageList.push_back(atomImage);
  imageIndex["atomSelected"] = imageList.size();
  QIcon atomSelectedImage=QIcon(QPixmap(atomSelected_xpm));
  imageList.push_back(atomSelectedImage);
  imageIndex["atomHidden"] = imageList.size();
  QIcon atomHiddenImage=QIcon(QPixmap(atomHidden_xpm));
  imageList.push_back(atomHiddenImage);
  imageIndex["atomHiddenSelected"] = imageList.size();
  QIcon atomHiddenSelectedImage=QIcon(QPixmap(atomHiddenSelected_xpm));
  imageList.push_back(atomHiddenSelectedImage);

  AssignImageList(imageList);

  std::string rootText = std::string("Atoms List");
  setHeaderLabel(rootText.c_str());
  rootId = invisibleRootItem();
  _menu = new MIMenu(*this);
  _menu->Append(ID_MODELSVIEW_ATOMSTREE_SHOW, "Show/Hide", "Show or hide this atom", false);
  _menu->Append(ID_MODELSVIEW_ATOMSTREE_DELETE, "Delete", "Delete this atom", false);
  _menu->Append(ID_MODELSVIEW_ATOMSTREE_COLOR, "Color", "Color this atom", false);
  _menu->Append(ID_MODELSVIEW_ATOMSTREE_EDIT, "Edit", "Edit the properties of this atom", false);

  // Solid surface menu created here, but filled when view set
  solidSurfMenu = new QMenu(this);
  solidSurfMenu->setTitle(tr("Solid Surface"));
  _menu->addMenu(solidSurfMenu);
  connect(solidSurfMenu, SIGNAL(aboutToShow()), this, SLOT(updateSolidSurfMenu()));

  _working = false;


  connect(this, SIGNAL(itemClicked(QTreeWidgetItem *, int)),
          this, SLOT(OnItemClicked(QTreeWidgetItem *, int)));

  connect(this, SIGNAL(itemActivated(QTreeWidgetItem *, int)),
          this, SLOT(OnItemActivated(QTreeWidgetItem *, int)));

BEGIN_EVENT_TABLE(this, NONE)
EVT_MENU(ID_MODELSVIEW_ATOMSTREE_DELETE, AtomsTree::DeleteItem)
EVT_MENU(ID_MODELSVIEW_ATOMSTREE_SHOW, AtomsTree::ShowItem)
EVT_MENU(ID_MODELSVIEW_ATOMSTREE_EDIT, AtomsTree::EditItem)
EVT_MENU(ID_MODELSVIEW_ATOMSTREE_COLOR, AtomsTree::ColorItem)

END_EVENT_TABLE()
}

AtomsTree::~AtomsTree() {
  setVisible(false);
  delete _menu;
}

static QAction* solidsurf_menu_action(QMenu* menu, QActionGroup* group, int actionId, const char* text) {
  QAction* a = menu->addAction(QObject::tr(text));
  a->setData(actionId);
  group->addAction(a);
  return a;
}

static void fillSolidSurfMenu(QWidget * /* parent */, MIGLWidget *view, QMenu *menu) {
  if (view) {
    foreach (QAction *a, view->solidSurfCommonActionGroup()->actions()) {
      menu->addAction(a);
    }
  }
}

void AtomsTree::setView(MIGLWidget* view) {
  this->view = view;
  solidSurfMenu->clear();
  if (view) {
    solidSurfActionGroup = new QActionGroup(this);
    connect(solidSurfActionGroup, SIGNAL(triggered(QAction*)), this, SLOT(solidSurfaceActionTriggered(QAction*)));

    solidsurf_menu_action(solidSurfMenu, solidSurfActionGroup, ID_SOLIDSURFACE_BUILD, "Build surface");
    solidsurf_menu_action(solidSurfMenu, solidSurfActionGroup, ID_SOLIDSURFACE_COLOR, "Color surface");
    solidsurf_menu_action(solidSurfMenu, solidSurfActionGroup, ID_SOLIDSURFACE_COLOR_BY_ATOM, "Color surface by atom type");
  }
  fillSolidSurfMenu(this, view, solidSurfMenu);
}

void AtomsTree::atomChanged(MIMoleculeBase*, MIAtomList& atoms) {
  MIAtom_iter iter = atoms.begin();
  while (iter != atoms.end()) {
    MIAtom* atom = *iter;
    ++iter;
    stylizeAtom(atom);
  }
}

void AtomsTree::atomsToBeDeleted(MIMoleculeBase*, const MIAtomList& atoms) {
  bool setCurrent=false;
  for (unsigned int i=0;i<atoms.size();++i) {
    MIAtom* atom=atoms[i];
    if (atomToData.find(atom) != atomToData.end()) {
      TreeData* data = atomToData[atom];
      if (data == NULL || !validTreeData(data)) {
        atomToData.erase(atom);
        continue;
      }
      QTreeWidgetItem* item = data->GetId();
      if (item) {
        if (currentAtom == atom) {
          setCurrent=true;
        }
        atomToData.erase(atom);
        Delete(item);
      }
    }
  }

  if (setCurrent)
    setCurrentAtom(NULL);
}

void AtomsTree::setModel(Molecule* model) {
  if (this->model != NULL) {
    disconnect(this->model, SIGNAL(atomChanged(chemlib::MIMoleculeBase*,chemlib::MIAtomList&)));
    disconnect(this->model, SIGNAL(atomsToBeDeleted(chemlib::MIMoleculeBase*,chemlib::MIAtomList)));
  }
  this->model = model;
  if (model != NULL) {
    connect(model, SIGNAL(atomChanged(chemlib::MIMoleculeBase*,chemlib::MIAtomList&)),
            this, SLOT(atomChanged(chemlib::MIMoleculeBase*,chemlib::MIAtomList&)));
    connect(model, SIGNAL(atomsToBeDeleted(chemlib::MIMoleculeBase*,chemlib::MIAtomList)),
            this, SLOT(atomsToBeDeleted(chemlib::MIMoleculeBase*,chemlib::MIAtomList)));
  }
}

void AtomsTree::setResidue(RESIDUE* residue) {
  currentAtom = NULL;
  atomToData.clear();
  this->residue = residue;
  setVisible(false);
  DeleteChildren(rootId);
  if (residue != NULL) {
    for (int ia = 0; ia < residue->atomCount(); ++ia) {
      MIAtom* atom = residue->atom(ia);
      TreeData* data = new TreeData;
      data->atom = atom;
      QTreeWidgetItem* item = appendItem(rootId, MIAtom::liststring(atom).c_str(), imageIndex["atom"], imageIndex["atomSelected"], data);
      if (!item) {
        Logger::log("Error adding atom to atoms list");
        continue;
      }
      atomToData[atom] = data;
      stylizeAtom(atom);
    }
    if (currentAtom == NULL) {
      setCurrentAtom(atom_from_name("CA", *residue));
    }
    if (currentAtom == NULL && residue->atomCount() > 0) {
      setCurrentAtom(residue->atom(0));
    }

  }
  setVisible(true);
  update();
}

void AtomsTree::setCurrentAtom(MIAtom* atom) {
  if (atom == currentAtom) {
    return;
  }
  MIAtom* oldAtom = currentAtom;
  currentAtom = atom;
  stylizeAtom(oldAtom);
  stylizeAtom(atom);
  if (atomToData.find(currentAtom) != atomToData.end()) {
    TreeData* data = atomToData[currentAtom];
    if (!validTreeData(data)) {
      atomToData.erase(currentAtom);
      return;
    }
    QTreeWidgetItem* item = data->GetId();
    if (item) {
      scrollToItem(item);
    }
  }
}

MIAtom* AtomsTree::getCurrentAtom() {
  return currentAtom;
}

void AtomsTree::stylizeAtom(MIAtom* atom) {
  if (atom == NULL) {
    return;
  }
  TreeData* data = atomToData[atom];
  if (!validTreeData(data)) {
    atomToData.erase(atom);
    return;
  }
  QTreeWidgetItem* item = data->GetId();
  if (!item) {
    return;
  }
  int image = imageIndex["atom"];
  if (atom == currentAtom) {
    if (atom->isHidden()) {
      image = imageIndex["atomHiddenSelected"];
    } else {
      image = imageIndex["atomSelected"];
    }
  } else if (atom->isHidden()) {
    image = imageIndex["atomHidden"];
  }
  item->setText(0, MIAtom::liststring(atom).c_str());
  item->setIcon(0,_imageList[ image]);
  item->setIcon(0,_imageList[ image ]);
}

void AtomsTree::OnItemClicked(QTreeWidgetItem *item, int) {
  if (_working) {
    return;
  }

  TreeData* data = (TreeData*) GetItemData(item);
  if (data == NULL || !validTreeData(data)) {
    return;
  }

  if (data->atom != NULL) {
    MIAtom* atom = data->atom;
    setCurrentAtom(atom);
    if (syncView && currentAtom != NULL) {
      view->select(model, residue, currentAtom);
    }
  }
}

void AtomsTree::OnItemActivated(QTreeWidgetItem *item, int) {
  TreeData* data = (TreeData*) GetItemData(item);
  if (data == NULL || !validTreeData(data)) {
    return;
  }
  if (data->atom != NULL) {
    MIAtom* atom = data->atom;
    view->moveTo(atom->x(), atom->y(), atom->z());
  }
}

void AtomsTree::contextMenuEvent(QContextMenuEvent* event) {
  QPoint pos = event->pos();
  QTreeWidgetItem* item = itemAt(pos);
  if (item && !selectedItems().contains(item)) {
    clearSelection();
    item->setSelected(true);
  }
  if (selectedItems().size()) {
    _menu->doExec(QCursor::pos());
  }
}

void AtomsTree::ShowItem() {
  QList<QTreeWidgetItem *> selected;
  GetSelections(selected);
  if (selected.size() == 0) {
    return;
  }
  MIAtomList atoms;
  for (int i = 0; i < selected.size(); ++i) {
    QTreeWidgetItem* item = selected[i];
    if (!item) {
      continue;
    }
    TreeData* data = (TreeData*) GetItemData(item);
    if (data == NULL || !validTreeData(data)) {
      continue;
    }
    if (data->atom != NULL) {
      MIAtom* atom = data->atom;
      atoms.push_back(atom);
    }
  }
  if (atoms.size())
    model->toggleAtomsHidden(atoms);
}

void AtomsTree::ColorItem() {
  QList<QTreeWidgetItem *> selected;
  GetSelections(selected);
  if (selected.size() == 0) {
    return;
  }
  MIAtomList atoms;
  for (int i = 0; i < selected.size(); ++i) {
    QTreeWidgetItem* item = selected[i];
    if (!item) {
      continue;
    }
    TreeData* data = (TreeData*) GetItemData(item);
    if (data == NULL || !validTreeData(data)) {
      continue;
    }
    if (MIAtom::isValid(data->atom)) {
      atoms.push_back(data->atom);
    }
  }
  if (atoms.size() == 0) {
    return;
  }
  int color = MIColorChooser(atoms[0]->color());
  model->setAtomsColor(atoms, color);
}

void AtomsTree::DeleteItem() {
  QList<QTreeWidgetItem *> selected;
  GetSelections(selected);
  if (selected.size() == 0) {
    return;
  }
  MIAtomList atoms;
  for (int i = 0; i < selected.size(); ++i) {
    QTreeWidgetItem* item = selected[i];
    if (!item) {
      continue;
    }
    TreeData* data = (TreeData*) GetItemData(item);
    if (data == NULL || !validTreeData(data)) {
      continue;
    }
    if (MIAtom::isValid(data->atom)) {
      atoms.push_back(data->atom);
    }
  }
  if (atoms.size() == 0) {
    return;
  }
  model->DeleteAtoms(atoms);
}

unsigned int CountAtoms(Molecule *model) {
  unsigned int count=0;
  for (MIIter<RESIDUE> currRes = model->GetResidues(); currRes; ++currRes) {
    count += currRes->atomCount();
  }
  return count;
}

void AtomsTree::updateSolidSurfMenu() {
  QList<QTreeWidgetItem *> selected;
  GetSelections(selected);
  if (selected.size() == 0) {
    return;
  }
  if (view) {
    view->updateSolidsurfMenu();

    foreach (QAction *action, solidSurfActionGroup->actions()) {
      switch (action->data().toInt()) {

        case ID_SOLIDSURFACE_COLOR:
        case ID_SOLIDSURFACE_COLOR_BY_ATOM:
          action->setEnabled(MISurfaceCount());
          break;

        case ID_SOLIDSURFACE_BUILD:
          action->setEnabled(true);
          break;

        default:
          break;
      }
    }
  }

}

void AtomsTree::solidSurfaceActionTriggered(QAction* action)
{
  QList<QTreeWidgetItem *> selected;
  GetSelections(selected);
  if (selected.size() == 0) {
    return;
  }

  std::set<MIAtom*> atoms;
  for (int i = 0; i < selected.size(); ++i) {
    QTreeWidgetItem* item = selected[i];
    if (!item) {
      continue;
    }
    TreeData* data = (TreeData*) GetItemData(item);
    if (data == NULL || !validTreeData(data)) {
      continue;
    }
    if (MIAtom::isValid(data->atom)) {
      atoms.insert(data->atom);
    }
  }
  if (atoms.size() == 0) {
    return;
  }

  std::vector<Molecule *> mols;
  std::vector<unsigned int> sel;

  mols.push_back(model);


  // get total atom count
  unsigned int count = CountAtoms(model);
  sel.resize(count);

  unsigned int atom_index = 0;
  for (MIIter<RESIDUE> currRes = model->GetResidues(); currRes; ++currRes) {
    for (int i = 0; i < currRes->atomCount(); ++i) {
      MIAtom* atom = currRes->atom(i);
      sel[atom_index] = (atoms.find(atom) != atoms.end());
      atom_index++;
    }
  }


  view->solidSurfaceCommand(action->data().toInt(), mols, sel);
}


void AtomsTree::EditItem() {
  QList<QTreeWidgetItem *> selected;
  GetSelections(selected);
  if (selected.size() == 0) {
    return;
  }

  MIAtomList atoms;
  for (int i = 0; i < selected.size(); ++i) {
    QTreeWidgetItem* item = selected[i];
    if (!item) {
      continue;
    }
    TreeData* data = (TreeData*) GetItemData(item);
    if (data == NULL || !validTreeData(data)) {
      continue;
    }
    if (MIAtom::isValid(data->atom)) {
      atoms.push_back(data->atom);
    }
  }
  if (atoms.size() == 0) {
    return;
  }
  std::string s;
  s=::format("%0.4f %0.4f", atoms[0]->BValue(), atoms[0]->occ());
  std::string str;
  float bvalue;
  float occ;
  MIGetStringDialog dlg(0, "Edit atom", "Edit B-Value and occupancy");
  if (dlg.GetValue(s.c_str(), str) && str.size()) {
    sscanf(str.c_str(), "%f%f", &bvalue, &occ);
  } else {
    return;
  }
  model->setAtomsBValueAndOccupancy(atoms, bvalue, occ);
}

class ResiduesTree : public MIQTreeWidget, public MIEventHandler {
  Q_OBJECT

  QTreeWidgetItem* rootId;
  AtomsTree* atomsTree;
  MIGLWidget* view;
  Molecule* model;
  RESIDUE* currentResidue;

  ImageIndexMap imageIndex;
  typedef std::map<RESIDUE*, TreeData*> ResidueToDataMap;
  ResidueToDataMap residueToData;

  MIMenu* _menu;
  QMenu *solidSurfMenu;
  QActionGroup* solidSurfActionGroup;
  bool _working;

  virtual void contextMenuEvent(QContextMenuEvent* event);

public:

  ResiduesTree(QWidget* parent);
  virtual ~ResiduesTree();

  void setView(MIGLWidget* view);
  void setAtomsTree(AtomsTree* tree);
  void setModel(Molecule* model);
  Molecule* getModel();
  void setChain(RESIDUE* chain);
  void setCurrentResidue(RESIDUE* residue);
  RESIDUE* getCurrentResidue();
  void stylizeResidue(RESIDUE* residue);

private slots:
  void residuesToBeDeleted(chemlib::MIMoleculeBase* model, std::vector<chemlib::RESIDUE*>& residues);
  void atomChanged(chemlib::MIMoleculeBase* model, chemlib::MIAtomList& atoms);

  void OnItemClicked(QTreeWidgetItem *item, int column); // single click
  void OnItemActivated(QTreeWidgetItem *item, int column); // double click

  void ShowItem();
  void ColorItem();
  void DeleteItem();
  void EditItem();
  void CopyItem();
  void InsertItem();
  void PasteItem();
  void EditItemAtoms();
  void solidSurfaceActionTriggered(QAction* action);
  void updateSolidSurfMenu();
};


ResiduesTree::ResiduesTree(QWidget* parent)
: MIQTreeWidget(parent),MIEventHandler(this),
  atomsTree(NULL), view(NULL), model(NULL), currentResidue(NULL) {

  setSelectionMode(QAbstractItemView::ExtendedSelection);

 std::vector<QIcon> imageList;
  imageIndex["residue"] = imageList.size();
  QIcon residueImage=QIcon(QPixmap(residue_xpm));
  imageList.push_back(residueImage);
  imageIndex["residueSelected"] = imageList.size();
  QIcon residueSelectedImage=QIcon(QPixmap(residueSelected_xpm));
  imageList.push_back(residueSelectedImage);
  imageIndex["residuePartial"] = imageList.size();
  QIcon residuePartialImage=QIcon(QPixmap(residuePartial_xpm));
  imageList.push_back(residuePartialImage);
  imageIndex["residuePartialSelected"] = imageList.size();
  QIcon residuePartialSelectedImage=QIcon(QPixmap(residuePartialSelected_xpm));
  imageList.push_back(residuePartialSelectedImage);
  imageIndex["residueHidden"] = imageList.size();
  QIcon residueHiddenImage=QIcon(QPixmap(residueHidden_xpm));
  imageList.push_back(residueHiddenImage);
  imageIndex["residueHiddenSelected"] = imageList.size();
  QIcon residueHiddenSelectedImage=QIcon(QPixmap(residueHiddenSelected_xpm));
  imageList.push_back(residueHiddenSelectedImage);

  AssignImageList(imageList);

  std::string rootText = std::string("Residues List");
  setHeaderLabel(rootText.c_str());
  rootId = invisibleRootItem();

  _menu = new MIMenu(*this);
  _menu->Append(ID_MODELSVIEW_RESIDUESTREE_SHOW, "Show/Hide", "Show or hide this residue", false);
  _menu->Append(ID_MODELSVIEW_RESIDUESTREE_DELETE, "Delete", "Delete this residue", false);
  _menu->Append(ID_MODELSVIEW_RESIDUESTREE_COLOR, "Color", "Color this residue", false);
  _menu->Append(ID_MODELSVIEW_RESIDUESTREE_EDIT, "Edit", "Edit the properties of this residue", false);
  _menu->Append(ID_MODELSVIEW_RESIDUESTREE_EDITATOMS, "Edit atoms", "Edit the properties of the atoms in this residue", false);
  _menu->Append(ID_MODELSVIEW_RESIDUESTREE_COPY, "Copy", "Copy residues", false);
  _menu->Append(ID_MODELSVIEW_RESIDUESTREE_INSERT, "Insert", "Insert residues", false);
  _menu->Append(ID_MODELSVIEW_RESIDUESTREE_PASTE, "Paste", "Paste residues", false);

  // Solid surface menu created here, but filled when view set
  solidSurfMenu = new QMenu(this);
  solidSurfMenu->setTitle(tr("Solid Surface"));
  _menu->addMenu(solidSurfMenu);
  connect(solidSurfMenu, SIGNAL(aboutToShow()), this, SLOT(updateSolidSurfMenu()));

  _working = false;


  connect(this, SIGNAL(itemClicked(QTreeWidgetItem *, int)),
          this, SLOT(OnItemClicked(QTreeWidgetItem *, int)));

  connect(this, SIGNAL(itemActivated(QTreeWidgetItem *, int)),
          this, SLOT(OnItemActivated(QTreeWidgetItem *, int)));

BEGIN_EVENT_TABLE(this, NONE)
EVT_MENU(ID_MODELSVIEW_RESIDUESTREE_DELETE, ResiduesTree::DeleteItem)
EVT_MENU(ID_MODELSVIEW_RESIDUESTREE_SHOW, ResiduesTree::ShowItem)
EVT_MENU(ID_MODELSVIEW_RESIDUESTREE_EDIT, ResiduesTree::EditItem)
EVT_MENU(ID_MODELSVIEW_RESIDUESTREE_COLOR, ResiduesTree::ColorItem)
EVT_MENU(ID_MODELSVIEW_RESIDUESTREE_COPY, ResiduesTree::CopyItem)
EVT_MENU(ID_MODELSVIEW_RESIDUESTREE_INSERT, ResiduesTree::InsertItem)
EVT_MENU(ID_MODELSVIEW_RESIDUESTREE_PASTE, ResiduesTree::PasteItem)
EVT_MENU(ID_MODELSVIEW_RESIDUESTREE_EDITATOMS, ResiduesTree::EditItemAtoms)

END_EVENT_TABLE()
}

ResiduesTree::~ResiduesTree() {
  setVisible(false);
  delete _menu;
}

void ResiduesTree::setView(MIGLWidget* view) {
  this->view = view;
  solidSurfMenu->clear();
  if (view) {
    solidSurfActionGroup = new QActionGroup(this);
    connect(solidSurfActionGroup, SIGNAL(triggered(QAction*)), this, SLOT(solidSurfaceActionTriggered(QAction*)));

    solidsurf_menu_action(solidSurfMenu, solidSurfActionGroup, ID_SOLIDSURFACE_BUILD, "Build surface");
    solidsurf_menu_action(solidSurfMenu, solidSurfActionGroup, ID_SOLIDSURFACE_COLOR, "Color surface");
    solidsurf_menu_action(solidSurfMenu, solidSurfActionGroup, ID_SOLIDSURFACE_COLOR_BY_ATOM, "Color surface by atom type");
  }
  fillSolidSurfMenu(this, view, solidSurfMenu);
}

void ResiduesTree::setAtomsTree(AtomsTree* tree) {
  atomsTree = tree;
}

void ResiduesTree::setModel(Molecule* model) {
  if (this->model != NULL) {
    disconnect(this->model, SIGNAL(residuesToBeDeleted(chemlib::MIMoleculeBase*,std::vector<chemlib::RESIDUE*>&)));
    disconnect(this->model, SIGNAL(atomChanged(chemlib::MIMoleculeBase*,chemlib::MIAtomList&)));
  }
  this->model = model;
  if (model != NULL) {
    connect(model, SIGNAL(residuesToBeDeleted(chemlib::MIMoleculeBase*,std::vector<chemlib::RESIDUE*>&)),
            this, SLOT(residuesToBeDeleted(chemlib::MIMoleculeBase*,std::vector<chemlib::RESIDUE*>&)));
    connect(model, SIGNAL(atomChanged(chemlib::MIMoleculeBase*,chemlib::MIAtomList&)),
            this, SLOT(atomChanged(chemlib::MIMoleculeBase*,chemlib::MIAtomList&)));
  }
}

Molecule* ResiduesTree::getModel() {
  return model;
}

void ResiduesTree::residuesToBeDeleted(MIMoleculeBase*, std::vector<RESIDUE*>& residues) {
  bool setCurrent=false;
  std::vector<RESIDUE*>::iterator iter;
  for (iter = residues.begin(); iter != residues.end(); ++iter) {
    RESIDUE* residue = *iter;
    if (residueToData.find(residue) != residueToData.end()) {
      TreeData* data = residueToData[residue];
      if (data == NULL || !validTreeData(data)) {
        residueToData.erase(residue);
        continue;
      }
      QTreeWidgetItem* item = data->GetId();
      if (item) {
        if (currentResidue == residue) {
          setCurrent=true;
        }
        residueToData.erase(residue);
        _working=true;
        Delete(item);
        _working=false;
	  }
    }
  }
  if (setCurrent)
    setCurrentResidue(NULL);
}


void ResiduesTree::setChain(RESIDUE* chain) {
  setVisible(false);
  setCurrentResidue(NULL);
  ResidueToDataMap::iterator iter = residueToData.begin();
  for (iter = residueToData.begin(); iter != residueToData.end(); ++iter) {
    invalidateTreeData(iter->second);
  }
  residueToData.clear();
  DeleteChildren(rootId);
  if (chain != NULL) {
    int chain_id = chain->chain_id();
    RESIDUE* res = chain;
    while ((res != NULL) && res->chain_id() == chain_id) {
      TreeData* data = new TreeData;
      data->residue = res;
      QTreeWidgetItem* item = appendItem(rootId, RESIDUE::liststring(res).c_str(), imageIndex["residue"], imageIndex["residueSelected"], data);
      if (!item) {
        continue;
      }
      residueToData[res] = data;
      stylizeResidue(res);
      res = res->next();
    }
  }
  setCurrentResidue(chain);
  setVisible(true);
  update();
}

void ResiduesTree::setCurrentResidue(RESIDUE* residue) {
  if (residue == currentResidue) {
    return;
  }
  RESIDUE* oldResidue = currentResidue;
  currentResidue = residue;
  stylizeResidue(oldResidue);
  stylizeResidue(residue);
  atomsTree->setModel(model);
  atomsTree->setResidue(residue);
  if (residueToData.find(currentResidue) != residueToData.end()) {
    TreeData* data = residueToData[currentResidue];
    if (!validTreeData(data)) {
      residueToData.erase(currentResidue);
      return;
    }
    QTreeWidgetItem* item = data->GetId();
    if (item) {
      scrollToItem(item);
    }
  }
}

RESIDUE* ResiduesTree::getCurrentResidue() {
  return currentResidue;
}

void ResiduesTree::atomChanged(MIMoleculeBase* model, MIAtomList& atoms) {
  std::set<RESIDUE*> residues;
  MIAtom_iter iter = atoms.begin();
  while (iter != atoms.end()) {
    MIAtom* atom = *iter;
    ++iter;
    RESIDUE* res = residue_from_atom(model->getResidues(), atom);
    residues.insert(res);
  }
  std::set<RESIDUE*>::iterator iter2 = residues.begin();
  while (iter2 != residues.end()) {
    RESIDUE* res = *iter2;
    ++iter2;
    stylizeResidue(res);
  }
}

void ResiduesTree::stylizeResidue(RESIDUE* residue) {
  if (residue == NULL) {
    return;
  }
  if (residueToData.find(residue) == residueToData.end()) {
    return;
  }
  TreeData* data = residueToData[residue];
  if (!validTreeData(data)) {
    residueToData.erase(residue);
    return;
  }
  QTreeWidgetItem* item = data->GetId();
  if (!item) {
    return;
  }

  bool partial = false;
  int shownAtoms = 0;
  int hiddenAtoms = 0;
  for (int ia = 0; ia < residue->atomCount(); ++ia) {
    MIAtom* atom = residue->atom(ia);
    if (atom->isHidden()) {
      ++hiddenAtoms;
      if (shownAtoms > 0) {
        partial = true;
        break;
      }
    } else {
      ++shownAtoms;
      if (hiddenAtoms > 0) {
        partial = true;
        break;
      }
    }
  }
  bool hidden = hiddenAtoms == residue->atomCount();

  int image = 0;
  if (residue == currentResidue) {
    if (hidden) {
      image = imageIndex["residueHiddenSelected"];
    } else if (partial) {
      image = imageIndex["residuePartialSelected"];
    } else {
      image = imageIndex["residueSelected"];
    }
  } else {
    if (hidden) {
      image = imageIndex["residueHidden"];
    } else if (partial) {
      image = imageIndex["residuePartial"];
    } else {
      image = imageIndex["residue"];
    }
  }
  item->setIcon(0,_imageList[ image]);
  item->setIcon(0,_imageList[ image ]);
  item->setText(0, RESIDUE::liststring(residue).c_str());
}

void ResiduesTree::OnItemClicked(QTreeWidgetItem *item, int) {
  if (_working) {
    return;
  }

  TreeData* data = (TreeData*) GetItemData(item);
  if (data == NULL || !validTreeData(data)) {
    return;
  }

  if (data->residue != NULL) {
    RESIDUE* residue = data->residue;
    setCurrentResidue(residue);
    if (syncView && currentResidue != NULL) {
      view->select(model, currentResidue, atomsTree->getCurrentAtom());
    }
  }
}

void ResiduesTree::OnItemActivated(QTreeWidgetItem *item, int) {
  if (view != NULL) {
    TreeData* data = (TreeData*) GetItemData(item);
    if (data == NULL || !validTreeData(data)) {
      return;
    }
    if (data->residue != NULL) {
      RESIDUE* residue = data->residue;
      view->setFocusResidue(residue);
    }
  }
}

void ResiduesTree::contextMenuEvent(QContextMenuEvent* event) {
  if (Application::instance()->GetResidueBuffer() == NULL) {
    _menu->Enable(ID_MODELSVIEW_RESIDUESTREE_INSERT, false);
    _menu->Enable(ID_MODELSVIEW_RESIDUESTREE_PASTE, false);
  } else {
    _menu->Enable(ID_MODELSVIEW_RESIDUESTREE_INSERT, true);
    _menu->Enable(ID_MODELSVIEW_RESIDUESTREE_PASTE, true);
  }

  QPoint pos = event->pos();
  QTreeWidgetItem* item = itemAt(pos);
  if (item && !selectedItems().contains(item)) {
    clearSelection();
    item->setSelected(true);
  }
  if (selectedItems().size()) {
    _menu->doExec(QCursor::pos());
  }
}

void ResiduesTree::ShowItem() {
  QList<QTreeWidgetItem *> selected;
  GetSelections(selected);
  if (selected.size() == 0) {
    return;
  }
  std::vector<RESIDUE*> residues;
  for (int i = 0; i < selected.size(); ++i) {
    QTreeWidgetItem* item = selected[i];
    if (!item) {
      continue;
    }
    TreeData* data = (TreeData*) GetItemData(item);
    if (data == NULL || !validTreeData(data)) {
      continue;
    }
    if (data->residue != NULL) {
      RESIDUE* residue = data->residue;
      residues.push_back(residue);
    }
  }
  if (residues.size())
    model->toggleResiduesHidden(residues);
}

void ResiduesTree::ColorItem() {
  QList<QTreeWidgetItem *> selected;
  GetSelections(selected);
  if (selected.size() == 0) {
    return;
  }
  std::vector<RESIDUE*> residues;
  for (int i = 0; i < selected.size(); ++i) {
    QTreeWidgetItem* item = selected[i];
    if (!item) {
      continue;
    }
    TreeData* data = (TreeData*) GetItemData(item);
    if (data == NULL || !validTreeData(data)) {
      continue;
    }
    if (data->residue != NULL) {
      RESIDUE* residue = data->residue;
      residues.push_back(residue);
    }
  }


  MIGenericDialog dlg(this, "Choose Color");
  MIData data;
  data["Color:"].color[0] = (unsigned char)view->WhenShownColor;
  data["Color:"].isColorIndex = true;

  data["Method:"].radio = view->WhenShownColorMethod;
  data["Method:"].radio_count = 7;
  data["Method:"].radio_labels.push_back("Carbon Only");
  data["Method:"].radio_labels.push_back("All Atoms");
  data["Method:"].radio_labels.push_back("Secondary Structure");
  data["Method:"].radio_labels.push_back("B-Value");
  data["Method:"].radio_labels.push_back("Atom Type");
  data["Method:"].radio_labels.push_back("Hydrophobicity");
  data["Method:"].radio_labels.push_back("Shapley");

  if (!dlg.GetResults(data)) {
    return;
  }
  int color = (int)data["Color:"].color[0];
  int colorMethod = (int)data["Method:"].radio;

  model->setResiduesColor(residues, color, colorMethod);
}

void ResiduesTree::DeleteItem() {
  QList<QTreeWidgetItem *> selected;
  GetSelections(selected);
  if (selected.size() == 0) {
    return;
  }
  std::vector<RESIDUE*> residues;
  for (int i = 0; i < selected.size(); ++i) {
    QTreeWidgetItem* item = selected[i];
    if (!item) {
      continue;
    }
    TreeData* data = (TreeData*) GetItemData(item);
    if (data == NULL || !validTreeData(data)) {
      continue;
    }
    if (data->residue != NULL) {
      RESIDUE* residue = data->residue;
      residues.push_back(residue);
    }
  }
  if (residues.size() == 1) {
    std::string mess;
    RESIDUE* residue = residues[0];
    mess=::format("Are you sure you want to delete residue %s?", residue->name().c_str());
    if (MIMessageBox(mess.c_str(), "Confirm Delete Residue", MIDIALOG_YES_NO) == MI_YES) {
      //REDUNDANT: view->Purge(residue);
      model->DeleteRes(residue);
    }
  } else if (residues.size() > 1) {
    std::string mess;
    mess=::format("Are you sure you want to delete the %d residues selected?", residues.size());
    if (MIMessageBox(mess.c_str(), "Confirm Delete Residues", MIDIALOG_YES_NO) == MI_YES) {
//REDUNDANT
//       std::vector<RESIDUE*>::iterator iter = residues.begin();
//       while (iter != residues.end()) {
//         RESIDUE* residue = *iter;
//         ++iter;
//         view->Purge(residue);
//       }
      model->DeleteResidues(residues);
    }
  }
}

void ResiduesTree::updateSolidSurfMenu() {
  QList<QTreeWidgetItem *> selected;
  GetSelections(selected);
  if (selected.size() == 0) {
    return;
  }
  if (view) {
    view->updateSolidsurfMenu();

    foreach (QAction *action, solidSurfActionGroup->actions()) {
      switch (action->data().toInt()) {

        case ID_SOLIDSURFACE_COLOR:
        case ID_SOLIDSURFACE_COLOR_BY_ATOM:
          action->setEnabled(MISurfaceCount() != 0);
          break;

        case ID_SOLIDSURFACE_BUILD:
          action->setEnabled(true);
          break;

        default:
          break;
      }
    }
  }
}

void ResiduesTree::solidSurfaceActionTriggered(QAction* action)
{
  QList<QTreeWidgetItem *> selected;
  GetSelections(selected);
  if (selected.size() == 0) {
    return;
  }

  std::set<RESIDUE*> residues;
  for (int i = 0; i < selected.size(); ++i) {
    QTreeWidgetItem* item = selected[i];
    if (!item) {
      continue;
    }
    TreeData* data = (TreeData*) GetItemData(item);
    if (data == NULL || !validTreeData(data)) {
      continue;
    }
    if (data->residue != NULL) {
      residues.insert(data->residue);
    }
  }
  if (residues.size() == 0) {
    return;
  }

  std::vector<Molecule *> mols;
  std::vector<unsigned int> sel;

  mols.push_back(model);

  // get total atom count
  unsigned int count = CountAtoms(model);
  sel.resize(count);

  unsigned int atom_index = 0;
  for (MIIter<RESIDUE> currRes = model->GetResidues(); currRes; ++currRes) {
    if (residues.find(currRes) != residues.end()) {
      for (int i = 0; i < currRes->atomCount(); ++i) {
        //MIAtom* atom = currRes->atom(i);
        sel[atom_index] = true;
        atom_index++;
      }
    }
  }


  view->solidSurfaceCommand(action->data().toInt(), mols, sel);
}


void ResiduesTree::EditItem() {
  QList<QTreeWidgetItem *> selected;
  GetSelections(selected);
  if (selected.size() == 0) {
    return;
  }

  std::string newname("");
  bool first=true;
  std::vector<RESIDUE *> residues;

  for (int i = 0; i < selected.size(); ++i) {
    QTreeWidgetItem* item = selected[i];
    if (!item) {
      continue;
    }
    TreeData* data = (TreeData*) GetItemData(item);
    if (data == NULL || !validTreeData(data)) {
      continue;
    }
    if (data->residue != NULL) {
      residues.push_back(data->residue);
      if (first) {
        first=false;
        RESIDUE* residue = data->residue;
        std::string s;
        s=::format("Enter new starting residue number (%d digits max)", (int)(MAXNAME -1));
        std::string str(residue->name());
        MIGetStringDialog dlg(0, "Edit residue number", s.c_str());
        if (!dlg.GetValue(str, newname) || !newname.size()) {
          return;
        }
      }
    }
  }

  if (residues.size()) {
    model->setResidueNames(residues, newname.c_str());
  }
}

void ResiduesTree::EditItemAtoms() {
  QList<QTreeWidgetItem *> selected;
  GetSelections(selected);
  if (selected.size() == 0) {
    return;
  }
  // Collect atoms from all residues
  MIAtomList atoms;
  for (int i = 0; i < selected.size(); ++i) {
    QTreeWidgetItem* item = selected[i];
    if (!item) {
      continue;
    }
    TreeData* data = (TreeData*) GetItemData(item);
    if (data == NULL || !validTreeData(data)) {
      continue;
    }
    RESIDUE* res = data->residue;
    if (Residue::isValid(res) && res->atomCount() > 0) {
      for (int j = 0; j < res->atomCount(); ++j) {
        atoms.push_back(res->atom(j));
      }
    }
  }
  if (atoms.size() == 0) {
    return;
  }
  // Prompt for new b-value and occupancy
  MIAtom* atom = atoms[0];
  std::string s;
  s=::format("%0.4f %0.4f", atom->BValue(), atom->occ());
  std::string str;
  float bvalue;
  float occ;
  MIGetStringDialog dlg(0, "Edit atom", "Edit B-Value and occupancy");
  if (dlg.GetValue(s.c_str(), str) && str.size()) {
    sscanf(str.c_str(), "%f%f", &bvalue, &occ);
  } else {
    return;
  }
  // Set new values on all atoms
  model->setAtomsBValueAndOccupancy(atoms, bvalue, occ);
}

void ResiduesTree::CopyItem() {
  QList<QTreeWidgetItem *> selected;
  GetSelections(selected);
  if (selected.size() == 0) {
    return;
  }
  Application::instance()->ClearResidueBuffer();
  for (int i = 0; i < selected.size(); ++i) {
    QTreeWidgetItem* item = selected[i];
    if (!item) {
      continue;
    }
    TreeData* data = (TreeData*) GetItemData(item);
    if (data == NULL || !validTreeData(data)) {
      continue;
    }
    if (data->residue != NULL) {
      RESIDUE* residue = data->residue;
      Application::instance()->CopyResidueBuffer(residue);
    }
  }
}

void ResiduesTree::InsertItem() {
  QList<QTreeWidgetItem *> selected;
  GetSelections(selected);
  if (selected.size() == 0) {
    return;
  }
  for (int i = 0; i < selected.size(); ++i) {
    QTreeWidgetItem* item = selected[i];
    if (!item) {
      continue;
    }
    TreeData* data = (TreeData*) GetItemData(item);
    if (data == NULL || !validTreeData(data)) {
      continue;
    }
    if (data->residue != NULL) {
      RESIDUE* residue = data->residue;
      unsigned short chain_id = residue->chain_id();
      RESIDUE* buffer = Application::instance()->GetResidueBuffer();
      model->InsertResidues(residue, buffer, 0, chain_id);
      model->SortChains();
      break;
    }
  }

  model->Build();

}

void ResiduesTree::PasteItem() {
  QList<QTreeWidgetItem *> selected;
  GetSelections(selected);
  if (selected.size() == 0) {
    return;
  }
  std::vector<RESIDUE*> toDelete;
  for (int i = 0; i < selected.size(); ++i) {
    QTreeWidgetItem* item = selected[i];
    if (!item) {
      continue;
    }
    TreeData* data = (TreeData*) GetItemData(item);
    if (data == NULL || !validTreeData(data)) {
      continue;
    }
    if (data->residue != NULL) {
      RESIDUE* residue = data->residue;
      toDelete.push_back(residue);
    }
  }
  std::vector<RESIDUE*> insertedResidues;
  for (int i = 0; i < selected.size(); ++i) {
    QTreeWidgetItem* item = selected[i];
    if (!item) {
      continue;
    }
    TreeData* data = (TreeData*) GetItemData(item);
    if (data == NULL || !validTreeData(data)) {
      continue;
    }
    if (data->residue != NULL) {
      RESIDUE* residue = data->residue;
      unsigned short chain_id = residue->chain_id();
      RESIDUE* buffer = Application::instance()->GetResidueBuffer();
      model->InsertResidues(residue, buffer, 0, chain_id);
      break;
    }
  }

  model->DeleteResidues(toDelete);

  model->Build();
  std::vector<RESIDUE*>::iterator iter = insertedResidues.begin();
  while (iter != insertedResidues.end()) {
    RESIDUE* res = *iter;
    if (iter == insertedResidues.begin()) {
      setCurrentResidue(res);
    }
    ++iter;
    TreeData* data = residueToData[res];
    if (data == NULL || !validTreeData(data)) {
      residueToData.erase(res);
      continue;
    }
    QTreeWidgetItem* item = data->GetId();
    if (!item) {
      continue;
    }
    item->setSelected(true);
  }
  setFocus(Qt::MouseFocusReason);
}


class ModelsTree : public MIQTreeWidget, public MIEventHandler {
  Q_OBJECT

  MIGLWidget* view;
  ResiduesTree* residuesTree;
  AtomsTree* atomsTree;
  Molecule* currentModel;
  RESIDUE* currentChain;

  bool colorDialog(int& color, int& colorMethod);
  void refreshChains(Molecule* model);
  RESIDUE* findChain(RESIDUE* residue);
  Molecule* findModel(RESIDUE* residue);

  int truncateWidth;

  MIMenu* _menu;
  QMenu *solidSurfMenu;
  QActionGroup* solidSurfActionGroup;

  std::map<int, QAction*> _menu_map;
  std::vector<int> _menu_item_ids;
  void ResetMenu(bool show_all = true);


  bool _working;

  virtual void contextMenuEvent(QContextMenuEvent* event);

public:

  ModelsTree(QWidget* parent);
  virtual ~ModelsTree();

  void setView(MIGLWidget* view);
  void setResiduesTree(ResiduesTree* tree);
  void setAtomsTree(AtomsTree* tree);
  Displaylist* displaylist();
  void addModels(Displaylist* displaylist);
  void addModel(Molecule* model);
  Molecule* getCurrentModel();
  void addModelCrystal(Molecule* model);
  void addMaps(Displaylist* displaylist);
  void addMap(EMap* map);
  void addMapCrystal(EMap* map);
  void addChains(Molecule* model);
  void setCurrentChain(RESIDUE* residue);
  RESIDUE* getCurrentChain();
  void stylizeModel(Molecule* model);
  void stylizeMap(EMap* map);
  void stylizeChain(RESIDUE* residue);
  void stylizeCrystal(CMapHeaderBase* mapHeader);
  void restylize();
  void setTruncateWidth(int width);
  std::string truncateLeft(const std::string& name, int shorter = 0);
  std::string truncateRight(const std::string& name, int shorter = 0);
  void goToResidue(const std::string& name);

private slots:
  void modelAdded(Molecule* model);
  void moleculeChanged(chemlib::MIMoleculeBase* model);
  void atomChanged(chemlib::MIMoleculeBase* model, chemlib::MIAtomList& atoms);
  void modelToBeDeleted(chemlib::MIMoleculeBase* model);
  void currentModelChanged(Molecule* oldModel, Molecule* newModel);
  void mapAdded(EMap* model);
  void mapToBeDeleted(EMap* model);
  void currentMapChanged(EMap* oldModel, EMap* newModel);
  void mapFftRecalculated(EMap* map);
  void residuesToBeDeleted(chemlib::MIMoleculeBase* model, std::vector<chemlib::RESIDUE*>& residues);
  void residuesDeleted(chemlib::MIMoleculeBase* model);
  void mapHeaderChanged(CMapHeaderBase* mapHeader);
  void focusResidueChanged(chemlib::RESIDUE* residue);
  void selectionChanged(Molecule* model, chemlib::RESIDUE* residue, chemlib::MIAtom* atom);

private:
  QTreeWidgetItem* rootId;
  QTreeWidgetItem* FindLastModelItem();

  ImageIndexMap imageIndex;

  typedef std::map<Molecule*, TreeData*> ModelToDataMap;
  ModelToDataMap modelToData;
  ModelToDataMap modelToCrystalData;

  typedef std::map<EMap*, TreeData*> MapToDataMap;
  MapToDataMap mapToData;
  MapToDataMap mapToCrystalData;

  typedef std::map<RESIDUE*, Molecule*> ResidueToModelMap;
  ResidueToModelMap chainToModel;
  typedef std::map<RESIDUE*, TreeData*> ResidueToDataMap;
  ResidueToDataMap residueToChainData;

  typedef std::map<CMapHeaderBase*, TreeData*> MapHeaderToDataMap;
  MapHeaderToDataMap mapHeaderToCrystalData;


private slots:
  void OnItemClicked(QTreeWidgetItem *item, int column); // single click
  void OnItemActivated(QTreeWidgetItem *item, int column); // double click

  void ShowItem();
  void ColorItem();
  void DeleteItem();
  void EditItem();
  void SortChains();
  void CopyItem();
  void InsertItem();
  void MapProperties();
  void MapExport();
  void ModelExport();
  void OnSaveToCrystals();
  void solidSurfaceActionTriggered(QAction* action);
  void updateSolidSurfMenu();

  void ForeignDispatcher(const MIActionEvent &evt);
};


ModelsTree::ModelsTree(QWidget* parent)
: MIQTreeWidget(parent),MIEventHandler(this),
  view(NULL), currentModel(NULL), currentChain(NULL) {

  setSelectionMode(QAbstractItemView::ExtendedSelection);


  std::vector<QIcon> imageList;
  QIcon modelImage=QIcon(QPixmap(model_xpm));
  imageList.push_back(modelImage);
  QIcon modelSelectedImage=QIcon(QPixmap(modelSelected_xpm));
  imageList.push_back(modelSelectedImage);
  QIcon mapImage=QIcon(QPixmap(map_xpm));
  imageList.push_back(mapImage);
  QIcon mapSelectedImage=QIcon(QPixmap(mapSelected_xpm));
  imageList.push_back(mapSelectedImage);
  QIcon crystalImage=QIcon(QPixmap(crystal_xpm));
  imageList.push_back(crystalImage);
  QIcon chainImage=QIcon(QPixmap(chain_xpm));
  imageList.push_back(chainImage);
  QIcon chainSelectedImage=QIcon(QPixmap(chainSelected_xpm));
  imageList.push_back(chainSelectedImage);
  QIcon chainPartialImage=QIcon(QPixmap(chainPartial_xpm));
  imageList.push_back(chainPartialImage);
  QIcon chainPartialSelectedImage=QIcon(QPixmap(chainPartialSelected_xpm));
  imageList.push_back(chainPartialSelectedImage);
  QIcon chainHiddenImage=QIcon(QPixmap(chainHidden_xpm));
  imageList.push_back(chainHiddenImage);
  QIcon chainHiddenSelectedImage=QIcon(QPixmap(chainHiddenSelected_xpm));
  imageList.push_back(chainHiddenSelectedImage);

  AssignImageList(imageList);

  std::string rootText = std::string("Models List");
  setHeaderLabel(rootText.c_str());
  rootId = invisibleRootItem();

  setTruncateWidth(width());

  _menu = new MIMenu(*this);
  _menu->Append(ID_MODELSVIEW_MODELSTREE_SHOW, "Show/Hide", "Show or hide this item", false);
  _menu->Append(ID_MODELSVIEW_MODELSTREE_DELETE, "Delete", "Delete this item", false);
  _menu->Append(ID_MODELSVIEW_MODELSTREE_COLOR, "Color", "Color for this item", false);
  _menu->Append(ID_MODELSVIEW_MODELSTREE_EDIT, "Edit", "Edit this item", false);
  _menu->Append(ID_MODELSVIEW_MODELSTREE_SORTCHAINS, "Sort Chains", "Sort the chains", false);
  _menu->Append(ID_MODELSVIEW_MODELSTREE_COPY, "Copy", "Copy", false);
  _menu->Append(ID_MODELSVIEW_MODELSTREE_INSERT, "Insert", "Insert", false);
  _menu->Append(ID_MODELSVIEW_MODELSTREE_MAPPROPERTIES, "Properties", "Properties", false);
  _menu->Append(ID_MAP_CONTOUR,"Contour","Contour this map", false);
  _menu->Append(ID_MAP_CONTOURLEVELS,"Contour options...","Contour options for this map", false);
  _menu->Append(ID_MAP_FFT, "FFT Phases...", "re-FFT the map from the phases to change resolution, coefficients, etc.", false);
  _menu->Append(ID_MAP_SFCALC, "Calculate Structure Factors...", "Calculate structure factors from the model", false);
  _menu->Append(ID_MAP_REINDEX, "Reindex Reflections", "Change to an alternate indexing of the data", false);
  _menu->Append(ID_MAP_ADDFREE, "Add Free R-flag", "Add A Free R-flag field to the data - do this only once for a given dataset", false);
  _menu->Append(ID_FIND_DENSITY, "Find Ligand Density", "Find unaccounted for density large enough to be a ligand", false);
  _menu->Append(ID_MODELSVIEW_MODELSTREE_MAPEXPORT, "Export", "Export", false);
  _menu->Append(ID_MODELSVIEW_MODELSTREE_MODELEXPORT, "Export PDB...");
  _menu->Append(ID_MODELSVIEW_MODELSTREE_SAVETOCRYSTALS, "Save to crystals", "Save to crystals", false);

  // Solid surface menu created here, but filled when view set
  solidSurfMenu = new QMenu(this);
  solidSurfMenu->setTitle(tr("Solid Surface"));
  _menu->addMenu(solidSurfMenu);
  connect(solidSurfMenu, SIGNAL(aboutToShow()), this, SLOT(updateSolidSurfMenu()));

  _menu_map[ID_MODELSVIEW_MODELSTREE_SHOW] = _menu->findAction(ID_MODELSVIEW_MODELSTREE_SHOW);
  _menu_map[ID_MODELSVIEW_MODELSTREE_DELETE] = _menu->findAction(ID_MODELSVIEW_MODELSTREE_DELETE);
  _menu_map[ID_MODELSVIEW_MODELSTREE_COLOR] = _menu->findAction(ID_MODELSVIEW_MODELSTREE_COLOR);
  _menu_map[ID_MODELSVIEW_MODELSTREE_EDIT] = _menu->findAction(ID_MODELSVIEW_MODELSTREE_EDIT);
  _menu_map[ID_MODELSVIEW_MODELSTREE_SORTCHAINS] = _menu->findAction(ID_MODELSVIEW_MODELSTREE_SORTCHAINS);
  _menu_map[ID_MODELSVIEW_MODELSTREE_COPY] = _menu->findAction(ID_MODELSVIEW_MODELSTREE_COPY);
  _menu_map[ID_MODELSVIEW_MODELSTREE_INSERT] = _menu->findAction(ID_MODELSVIEW_MODELSTREE_INSERT);
  _menu_map[ID_MODELSVIEW_MODELSTREE_MAPPROPERTIES] = _menu->findAction(ID_MODELSVIEW_MODELSTREE_MAPPROPERTIES);
  _menu_map[ID_MAP_CONTOUR] = _menu->findAction(ID_MAP_CONTOUR);
  _menu_map[ID_MAP_CONTOURLEVELS] = _menu->findAction(ID_MAP_CONTOURLEVELS);
  _menu_map[ID_MAP_FFT] = _menu->findAction(ID_MAP_FFT);
  _menu_map[ID_MAP_SFCALC] = _menu->findAction(ID_MAP_SFCALC);
  _menu_map[ID_MAP_REINDEX] = _menu->findAction(ID_MAP_REINDEX);
  _menu_map[ID_MAP_ADDFREE] = _menu->findAction(ID_MAP_ADDFREE);
  _menu_map[ID_FIND_DENSITY] = _menu->findAction(ID_FIND_DENSITY);
  _menu_map[ID_MODELSVIEW_MODELSTREE_MAPEXPORT] = _menu->findAction(ID_MODELSVIEW_MODELSTREE_MAPEXPORT);
  _menu_map[ID_MODELSVIEW_MODELSTREE_MODELEXPORT] = _menu->findAction(ID_MODELSVIEW_MODELSTREE_MODELEXPORT);
  _menu_map[ID_MODELSVIEW_MODELSTREE_SAVETOCRYSTALS] = _menu->findAction(ID_MODELSVIEW_MODELSTREE_SAVETOCRYSTALS);
  //_menu_map[ID_SOLIDSURFACE_MENU] = _menu->findAction(ID_SOLIDSURFACE_MENU);

  _menu_item_ids.push_back(ID_MODELSVIEW_MODELSTREE_SHOW);
  _menu_item_ids.push_back(ID_MODELSVIEW_MODELSTREE_DELETE);
  _menu_item_ids.push_back(ID_MODELSVIEW_MODELSTREE_COLOR);
  _menu_item_ids.push_back(ID_MODELSVIEW_MODELSTREE_EDIT);
  _menu_item_ids.push_back(ID_MODELSVIEW_MODELSTREE_SORTCHAINS);
  _menu_item_ids.push_back(ID_MODELSVIEW_MODELSTREE_COPY);
  _menu_item_ids.push_back(ID_MODELSVIEW_MODELSTREE_INSERT);
  _menu_item_ids.push_back(ID_MODELSVIEW_MODELSTREE_MAPPROPERTIES);
  _menu_item_ids.push_back(ID_MAP_CONTOUR);
  _menu_item_ids.push_back(ID_MAP_CONTOURLEVELS);
  _menu_item_ids.push_back(ID_MAP_FFT);
  _menu_item_ids.push_back(ID_MAP_SFCALC);
  _menu_item_ids.push_back(ID_MAP_REINDEX);
  _menu_item_ids.push_back(ID_MAP_ADDFREE);
  _menu_item_ids.push_back(ID_FIND_DENSITY);
  _menu_item_ids.push_back(ID_MODELSVIEW_MODELSTREE_MAPEXPORT);
  _menu_item_ids.push_back(ID_MODELSVIEW_MODELSTREE_MODELEXPORT);
  _menu_item_ids.push_back(ID_MODELSVIEW_MODELSTREE_SAVETOCRYSTALS);
  //_menu_item_ids.push_back(ID_SOLIDSURFACE_MENU);
  _working = false;


  connect(this, SIGNAL(itemClicked(QTreeWidgetItem *, int)),
          this, SLOT(OnItemClicked(QTreeWidgetItem *, int)));

  connect(this, SIGNAL(itemActivated(QTreeWidgetItem *, int)),
          this, SLOT(OnItemActivated(QTreeWidgetItem *, int)));

BEGIN_EVENT_TABLE(this, NONE)
EVT_MENU(ID_MODELSVIEW_MODELSTREE_DELETE, ModelsTree::DeleteItem)
EVT_MENU(ID_MODELSVIEW_MODELSTREE_SHOW, ModelsTree::ShowItem)
EVT_MENU(ID_MODELSVIEW_MODELSTREE_COLOR, ModelsTree::ColorItem)
EVT_MENU(ID_MODELSVIEW_MODELSTREE_EDIT, ModelsTree::EditItem)
EVT_MENU(ID_MODELSVIEW_MODELSTREE_SORTCHAINS, ModelsTree::SortChains)
EVT_MENU(ID_MODELSVIEW_MODELSTREE_MAPPROPERTIES, ModelsTree::MapProperties)
EVT_MENU(ID_MODELSVIEW_MODELSTREE_MAPEXPORT, ModelsTree::MapExport)
EVT_MENU(ID_MODELSVIEW_MODELSTREE_MODELEXPORT, ModelsTree::ModelExport)
EVT_MENU(ID_MODELSVIEW_MODELSTREE_SAVETOCRYSTALS, ModelsTree::OnSaveToCrystals)
EVT_MENU(ID_MODELSVIEW_MODELSTREE_COPY, ModelsTree::CopyItem)
EVT_MENU(ID_MODELSVIEW_MODELSTREE_INSERT, ModelsTree::InsertItem)

EVT_MENU(ID_MAP_CONTOUR, ModelsTree::ForeignDispatcher)
EVT_MENU(ID_MAP_CONTOURLEVELS, ModelsTree::ForeignDispatcher)
EVT_MENU(ID_MAP_FFT, ModelsTree::ForeignDispatcher)
EVT_MENU(ID_MAP_SFCALC, ModelsTree::ForeignDispatcher)
EVT_MENU(ID_MAP_REINDEX, ModelsTree::ForeignDispatcher)
EVT_MENU(ID_MAP_ADDFREE, ModelsTree::ForeignDispatcher)
EVT_MENU(ID_FIND_DENSITY, ModelsTree::ForeignDispatcher)

END_EVENT_TABLE()
}

ModelsTree::~ModelsTree() {
  delete _menu;
}

QTreeWidgetItem* ModelsTree::FindLastModelItem() {
  QTreeWidgetItem* lastModel = 0;
  QTreeWidgetItemIterator it(this);
  for (; *it; ++it) {
    TreeData* data = (TreeData*) GetItemData(*it);
    if (data == NULL || !validTreeData(data)) {
      continue;
    }
    if (data->model != NULL) {
      lastModel=*it;
    }
  }

  /* Not found */
  return lastModel;
}


void ModelsTree::contextMenuEvent(QContextMenuEvent* event) {
  ResetMenu(false);

  QPoint pos = event->pos();
  QTreeWidgetItem* item = itemAt(pos);
  if (item && !selectedItems().contains(item)) {
    clearSelection();
    item->setSelected(true);
  }
  TreeData* data = (TreeData*) GetItemData(item);
  if (data && data->map) {
    displaylist()->SetCurrentMap(data->map);
  }
  if (selectedItems().size()) {
    _menu->doExec(QCursor::pos());
  }

  ResetMenu(true); // must re-enable all items to make history happy
}

void ModelsTree::OnItemClicked(QTreeWidgetItem *item, int) {
  if (_working) {
    return;
  }
  TreeData* data = (TreeData*) GetItemData(item);
  if (data == NULL || !validTreeData(data)) {
    return;
  }

  if (data->model != NULL) {
    Molecule* model = data->model;
    displaylist()->SetCurrent(model);
  } else if (data->map != NULL) {
    EMap* map = data->map;
    displaylist()->SetCurrentMap(map);
  } else if (data->chain != NULL) {
    setCurrentChain(data->chain);
    if (syncView) {
      view->select(residuesTree->getModel(), residuesTree->getCurrentResidue(), atomsTree->getCurrentAtom());
    }
  }
}

void ModelsTree::ForeignDispatcher(const MIActionEvent &evt) {

  switch (evt.GetId()) {
    case ID_MAP_CONTOUR:
      view->OnMapContour();
      break;
    case ID_MAP_CONTOURLEVELS:
      view->OnMapContourLevels();
      break;
    case ID_MAP_FFT:
      view->OnMapFFT();
      break;
    case ID_MAP_SFCALC:
      view->OnMapSFCalc();
      break;
    case ID_MAP_REINDEX:
      view->OnMapReindex();
      break;
    case ID_MAP_ADDFREE:
      view->OnMapAddFree();
      break;
    case ID_FIND_DENSITY:
      view->OnFindLigandDensity();
      break;
  }
}

void ModelsTree::updateSolidSurfMenu() {
  QList<QTreeWidgetItem *> selected;
  GetSelections(selected);
  if (selected.size() == 0) {
    return;
  }
  if (view) {
    view->updateSolidsurfMenu();

    foreach (QAction *action, solidSurfActionGroup->actions()) {
      switch (action->data().toInt()) {

        case ID_SOLIDSURFACE_COLOR:
        case ID_SOLIDSURFACE_COLOR_BY_ATOM:
          action->setEnabled(MISurfaceCount() != 0);
          break;

        case ID_SOLIDSURFACE_BUILD:
          action->setEnabled(true);
          break;

        default:
          break;
      }
    }
  }
}

void ModelsTree::solidSurfaceActionTriggered(QAction* action) {
  QList<QTreeWidgetItem *> selected;
  GetSelections(selected);
  if (selected.size() == 0) {
    return;
  }

  std::set<Molecule*> mols;
  std::set<RESIDUE *> chains;

  for (int i = 0; i < selected.size(); ++i) {
    QTreeWidgetItem* item = selected[i];
    if (!item) {
      continue;
    }
    TreeData* data = (TreeData*) GetItemData(item);
    if (data == NULL || !validTreeData(data)) {
      continue;
    }

    if (data->model != NULL) {
      mols.insert(data->model);
    } else if (data->chain != NULL) {
      chains.insert(data->chain);
    }
  }

  if (mols.size() == 0 && chains.size()==0) {
    return;
  }

  //if all chains in model selected, select model
  bool invalid_iterator=false;
  do {
    invalid_iterator=false;

    for (std::set<RESIDUE*>::iterator i=chains.begin(); i != chains.end(); ++i) {
      RESIDUE *chain=*i;
      Molecule *mol=chainToModel[chain];

      std::vector<RESIDUE*> chains_for_mol;
      for (std::map<RESIDUE*, Molecule*>::iterator j=chainToModel.begin(); j!=chainToModel.end(); ++j) {
        if (j->second == mol) {
          chains_for_mol.push_back(j->first);
        }
      }

      bool found_unselected_chain=false;
      for (unsigned int j=0; j < chains_for_mol.size(); ++j) {
        if (chains.find(chains_for_mol[j]) == chains.end()) {
          found_unselected_chain=true;
          break;
        }
      }

      if (!found_unselected_chain) {
        // if we're here, we have all chains in the model, so...
        // add to models; remove chains from set
        mols.insert(mol);
        for (unsigned int j=0; j < chains_for_mol.size(); ++j) {
          chains.erase(chains_for_mol[j]);
        }
        invalid_iterator=true;
        break; // break out of now-invalid iteration loop for "i"
      }
    }
  } while (invalid_iterator);

  // at this point we have mols==>complete molecules and chains==>partial molecules
  std::set<Molecule*> partial_mols;
  for (std::set<RESIDUE*>::iterator i=chains.begin(); i != chains.end(); ++i) {
    RESIDUE *chain=*i;
    Molecule *mol=chainToModel[chain];
    partial_mols.insert(mol);
  }


  // count total atoms
  unsigned int count=0;
  std::map<Molecule*, unsigned int> mol_atom_count;
  std::vector<Molecule *> surface_mols;
  for (std::set<Molecule*>::iterator i=mols.begin(); i != mols.end(); ++i) {
    mol_atom_count[*i]=CountAtoms(*i);
    count += mol_atom_count[*i];
    surface_mols.push_back(*i);
  }
  for (std::set<Molecule*>::iterator i=partial_mols.begin(); i != partial_mols.end(); ++i) {
    mol_atom_count[*i]=CountAtoms(*i);
    count += mol_atom_count[*i];
    surface_mols.push_back(*i);

  }

  if (count == 0) {
    return;
  }
  std::vector<unsigned int> sel;
  sel.resize(count);
  unsigned int offset=0;
  for (std::set<Molecule*>::iterator i=mols.begin(); i != mols.end(); ++i) {
    unsigned int start=offset;
    unsigned int end=mol_atom_count[*i] + offset;

    unsigned int *sp=&sel[start]; // avoid deref in loop
    for (unsigned int j=start; j < end; ++j) {
      *sp++=1;
    }
    offset+=mol_atom_count[*i];
  }

  for (std::set<Molecule*>::iterator i=partial_mols.begin(); i != partial_mols.end(); ++i) {
    Molecule *mol=*i;

    unsigned short last_chain_id;
    bool in_sel_chain=false;

    for (MIIter<RESIDUE> currRes = mol->GetResidues(); currRes; ++currRes) {

      if (chains.find(currRes) != chains.end()) {
        // begin of chain
        in_sel_chain=true;
        last_chain_id=currRes->chain_id();
      } else if (in_sel_chain && currRes->chain_id() != last_chain_id) {
        // end of chain
        in_sel_chain=false;
      }

      if (in_sel_chain) {
        for (int j=0; j < currRes->atomCount(); ++j) {
          sel[offset++]=1;
        }
      } else {
        offset += currRes->atomCount();
      }
    }
  }

  view->solidSurfaceCommand(action->data().toInt(), surface_mols, sel);
}

void ModelsTree::setCurrentChain(RESIDUE* residue) {
  if (residue == currentChain) {
    return;
  }
  RESIDUE* oldChain = currentChain;
  currentChain = residue;
  if (currentChain != NULL) {
    if (chainToModel.find(currentChain) == chainToModel.end()) {
      chainToModel[currentChain] = findModel(residue);
    }
    currentModel = chainToModel[currentChain];
  } else {
    currentModel = NULL;
  }
  stylizeChain(oldChain);
  stylizeChain(currentChain);
  residuesTree->setModel(currentModel);
  residuesTree->setChain(currentChain);
}

void ModelsTree::goToResidue(const std::string& nameRef) {
  std::string name = nameRef;
  MIStringTrim(name,true);
  MIStringTrim(name,false);

  std::string residueName = name;
  char chain_id;

  RESIDUE* chain = NULL;
  RESIDUE* residue = NULL;
  if (currentChain != NULL) {
    // Search current chain using full string
    chain = currentChain;
    chain_id = currentChain->chain_id()&255;
    residue = residue_from_name(chain, residueName.c_str(), chain_id);
  }
  if (residue == NULL && currentModel != NULL) {
    // Search model using full string
    chain = currentModel->getResidues();
    chain_id = '*';
    residue = residue_from_name(chain, residueName.c_str(), chain_id);
  }
  if (residue == NULL) {
    // Search model using first character as chain_id and rest as name
    char chain_id = name.size()>0?name[0]:' ';
    if (chain_id == '_') {
      chain_id = ' ';
    }
    residueName = &name[1];
    MIStringTrim(residueName,false);
    MIStringTrim(residueName,true);
    if (residueName.empty() && residuesTree->getCurrentResidue() != NULL) {
      residueName = residuesTree->getCurrentResidue()->name().c_str();
    }
    if (currentModel != NULL) {
      chain = currentModel->getResidues();
      residue = residue_from_name(chain, residueName.c_str(), chain_id);
    }
  }

  if (residue != NULL) {
    setCurrentChain(findChain(residue));
    residuesTree->setCurrentResidue(residue);
    if (syncView) {
      view->select(residuesTree->getModel(), residuesTree->getCurrentResidue(), atomsTree->getCurrentAtom());
      view->setFocusResidue(residuesTree->getCurrentResidue(), true);
    }
  } else if (name != MIToUpper(name)) {
    goToResidue(MIToUpper(name));
  }
}

RESIDUE* ModelsTree::getCurrentChain() {
  return currentChain;
}

void ModelsTree::OnItemActivated(QTreeWidgetItem* item, int) {
  if (view == NULL) {
    return;
  }
  TreeData* data = (TreeData*) GetItemData(item);
  if (data == NULL || !validTreeData(data)) {
    return;
  }

  if (data->chain != NULL) {
    RESIDUE* residue = data->chain;
    view->setFocusResidue(residue);
  }
}

void ModelsTree::ResetMenu(bool show_all) {
  //remove all items
  for (size_t i = 0; i < _menu_item_ids.size(); ++i) {
    QAction *act=_menu->findAction(_menu_item_ids[i]);
    if (act) {
      act->setVisible(false);
    }
  }

  if (show_all) {
    // add all items back in, in order
    for (size_t i = 0; i < _menu_item_ids.size(); ++i) {
      _menu_map[_menu_item_ids[i]]->setVisible(true);
      _menu->Enable(_menu_item_ids[i], true);
    }
    return;
  }

  QList<QTreeWidgetItem *> selected;
  GetSelections(selected);
  bool modelSelected = false;
  bool chainSelected = false;
  bool mapSelected = false;
  bool mapHasPhases = false;
  bool crystalSelected = false;
  for (int i = 0; i < selected.size(); ++i) {
    QTreeWidgetItem* item = selected[i];
    if (!item) {
      continue;
    }
    TreeData* data = (TreeData*) GetItemData(item);
    if (data == NULL || !validTreeData(data)) {
      continue;
    }

    if (data->model != NULL) {
      modelSelected = true;
    } else if (data->chain != NULL) {
      chainSelected = true;
    } else if (data->map != NULL) {
      EMap* map = data->map;
      mapHasPhases = map->HasPhases();
      mapSelected = true;
    } else if (data->mapHeader != NULL) {
      crystalSelected = true;
    }
  }

  if (modelSelected || chainSelected || mapSelected) {
    _menu_map[ID_MODELSVIEW_MODELSTREE_SHOW]->setVisible(true);
    _menu_map[ID_MODELSVIEW_MODELSTREE_DELETE]->setVisible(true);
    if (!mapSelected) {
      _menu_map[ID_MODELSVIEW_MODELSTREE_COLOR]->setVisible(true);
    }
  }
  _menu_map[ID_MODELSVIEW_MODELSTREE_EDIT]->setVisible(true);
  if (modelSelected || mapSelected) {
    _menu->Enable(ID_MODELSVIEW_MODELSTREE_EDIT, false);
  }
  if (chainSelected) {
    _menu_map[ID_MODELSVIEW_MODELSTREE_SORTCHAINS]->setVisible(true);
  }
  if (modelSelected || chainSelected) {
    _menu_map[ID_MODELSVIEW_MODELSTREE_COPY]->setVisible(true);
    _menu_map[ID_MODELSVIEW_MODELSTREE_INSERT]->setVisible(true);
    if (Application::instance()->GetResidueBuffer() == NULL) {
      _menu->Enable(ID_MODELSVIEW_MODELSTREE_INSERT, false);
    }
    //_menu_map[ID_SOLIDSURFACE_MENU]->setVisible(true);
  }
  if (mapSelected) {
    _menu_map[ID_MODELSVIEW_MODELSTREE_MAPPROPERTIES]->setVisible(true);
    _menu_map[ID_MAP_CONTOUR]->setVisible(true);
    _menu_map[ID_MAP_CONTOURLEVELS]->setVisible(true);
    _menu_map[ID_MAP_FFT]->setVisible(true);
    _menu_map[ID_MAP_SFCALC]->setVisible(true);
    _menu_map[ID_MAP_REINDEX]->setVisible(true);
    _menu_map[ID_MAP_ADDFREE]->setVisible(true);
    _menu_map[ID_FIND_DENSITY]->setVisible(true);
    _menu_map[ID_MODELSVIEW_MODELSTREE_MAPEXPORT]->setVisible(true);
    _menu->Enable(ID_MAP_FFT, mapHasPhases);
    _menu->Enable(ID_MODELSVIEW_MODELSTREE_MAPEXPORT, mapHasPhases);
  }
  if (modelSelected) {
    _menu_map[ID_MODELSVIEW_MODELSTREE_MODELEXPORT]->setVisible(true);
  }
  if (crystalSelected) {
    _menu_map[ID_MODELSVIEW_MODELSTREE_SAVETOCRYSTALS]->setVisible(true);
  }
}



void ModelsTree::OnSaveToCrystals() {

  QList<QTreeWidgetItem *> selected;
  GetSelections(selected);
  if (selected.size() == 0) {
    return;
  }
  std::set<Molecule*> chainsModified;
  for (int i = 0; i < selected.size(); ++i) {
    QTreeWidgetItem* item = selected[i];
    if (!item) {
      continue;
    }
    TreeData* data = (TreeData*) GetItemData(item);
    if (data == NULL || !validTreeData(data)) {
      continue;
    }

    if (data->mapHeader != NULL) {
		MIGetStringDialog dlg(NULL, "Crystal name","Crystal name:");
      std::string crystal;
      if (dlg.GetValue("",crystal) && crystal.size()) {
        data->mapHeader->SaveCrystal(std::string(crystal.c_str()));
      }
    }
  }
}

void ModelsTree::ShowItem() {
  QList<QTreeWidgetItem *> selected;
  GetSelections(selected);
  if (selected.size() == 0) {
    return;
  }
  for (int i = 0; i < selected.size(); ++i) {
    QTreeWidgetItem* item = selected[i];
    if (!item) {
      continue;
    }
    TreeData* data = (TreeData*) GetItemData(item);
    if (data == NULL || !validTreeData(data)) {
      continue;
    }

    if (data->model != NULL) {
      Molecule* model = data->model;
      if (model->Visible()) {
        model->Hide();
      } else {
        model->Show();
      }
    } else if (data->map != NULL) {
      EMap* map = data->map;
      if (map->Visible()) {
        map->Hide();
      } else {
        map->Show();
      }
    } else if (data->chain != NULL) {
      RESIDUE* chain = data->chain;
      Molecule* model = chainToModel[chain];
      model->toggleChainHidden(chain);
    }
  }

}

void ModelsTree::ColorItem() {
  QList<QTreeWidgetItem *> selected;
  GetSelections(selected);
  if (selected.size() == 0) {
    return;
  }
  std::vector<Molecule*> models;
  std::vector<RESIDUE*> chains;
  std::vector<EMap*> maps;
  for (int i = 0; i < selected.size(); ++i) {
    QTreeWidgetItem* item = selected[i];
    if (!item) {
      continue;
    }
    TreeData* data = (TreeData*) GetItemData(item);
    if (data == NULL || !validTreeData(data)) {
      continue;
    }

    if (data->model != NULL) {
      models.push_back(data->model);
    } else if (data->chain != NULL) {
      chains.push_back(data->chain);
    } else if (data->map != NULL) {
      maps.push_back(data->map);
    }
  }
  int color;
  int colorMethod;
  if (models.size() > 0 || chains.size() > 0) {
    if (!colorDialog(color, colorMethod)) {
      return;
    }
  }
  if (models.size() > 0) {
    std::vector<Molecule*>::iterator iter = models.begin();
    for (; iter != models.end(); ++iter) {
      Molecule* model = *iter;
      model->setColor(color, colorMethod);
    }
  }
  if (chains.size() > 0) {
    std::vector<RESIDUE*>::iterator iter = chains.begin();
    for (; iter != chains.end(); ++iter) {
      RESIDUE* chain = *iter;
      Molecule* model = chainToModel[chain];
      model->setChainColor(chain, color, colorMethod);
    }
  }
  if (maps.size() > 0) {
    std::vector<EMap*>::iterator iter = maps.begin();
    for (; iter != maps.end(); ++iter) {
      EMap* map = *iter;
      map->ContourLevels();
    }
  }
}

void ModelsTree::DeleteItem() {
  QList<QTreeWidgetItem *> selected;
  GetSelections(selected);
  if (selected.size() == 0) {
    return;
  }
  std::vector<Molecule*> models;
  std::vector<RESIDUE*> chains;
  std::vector<EMap*> maps;
  for (int i = 0; i < selected.size(); ++i) {
    QTreeWidgetItem* item = selected[i];
    if (!item) {
      continue;
    }
    TreeData* data = (TreeData*) GetItemData(item);
    if (data == NULL || !validTreeData(data)) {
      continue;
    }

    if (data->model != NULL) {
      models.push_back(data->model);
    } else if (data->chain != NULL) {
      chains.push_back(data->chain);
    } else if (data->map != NULL) {
      maps.push_back(data->map);
    }
  }

  if (models.size() > 0) {
    std::string mess;
    if (models.size() == 1) {
      Molecule* model = models[0];
      mess=::format("Are you sure you want to delete model\n%s?", model->pathname.c_str());
    } else {
      mess=::format("Are you sure you want to delete the %d models selected?", models.size());
    }
    if (MIMessageBox(mess.c_str(), "Confirm Delete Model", MIDIALOG_YES_NO) == MI_YES) {
      std::vector<Molecule*>::iterator iter = models.begin();
      for (; iter != models.end(); ++iter) {
        Molecule* model = *iter;
        mess=::format("Deleted model %s", model->pathname.c_str());
        view->SaveModelFile(model, mess.c_str());
        //view->GetDisplaylist()->DeleteItem(model);
        delete model;
      }
    }
  } else if (chains.size() > 0) {
    std::string mess;
    if (chains.size() == 1) {
      RESIDUE* chain = chains[0];
      mess=::format("Are you sure you want to delete chain\n%s?", chainstring(chain).c_str());
    } else {
      mess=::format("Are you sure you want to delete the %d chains selected?", chains.size());
    }
    if (MIMessageBox(mess.c_str(), "Confirm Delete Chains", MIDIALOG_YES_NO) == MI_YES) {
      std::vector<RESIDUE*>::iterator iter = chains.begin();
      for (; iter != chains.end(); ++iter) {
        RESIDUE* chain = *iter;
        Molecule* model = chainToModel[chain];
        mess=::format("Deleted chain %s", chainstring(chain).c_str());
        view->SaveModelFile(model, mess.c_str());

        if (currentChain == chain) {
          currentChain = NULL;
        }
        //REDUNDANT: view->PurgeChain(chain);
        model->DeleteChain(chain);
      }
    }
  } else if (maps.size() > 0) {
    std::string mess;
    if (maps.size() == 1) {
      EMap* map = maps[0];
      mess=::format("Are you sure you want to delete map\n%s?", map->mapName.c_str());
    } else {
      mess=::format("Are you sure you want to delete the %d maps selected?", maps.size());
    }
    if (MIMessageBox(mess.c_str(), "Confirm Delete Maps", MIDIALOG_YES_NO) == MI_YES) {
      std::vector<EMap*>::iterator iter = maps.begin();
      for (; iter != maps.end(); ++iter) {
        EMap* map = *iter;
        view->Purge(map);
      }
    }
  }
}

void ModelsTree::EditItem() {
  QList<QTreeWidgetItem *> selected;
  GetSelections(selected);
  if (selected.size() == 0) {
    return;
  }
  std::vector<RESIDUE*> chains;
  std::vector<EMap*> maps;
  std::vector<CMapHeaderBase*> crystals;
  for (int i = 0; i < selected.size(); ++i) {
    QTreeWidgetItem* item = selected[i];
    if (!item) {
      continue;
    }
    TreeData* data = (TreeData*) GetItemData(item);
    if (data == NULL || !validTreeData(data)) {
      continue;
    }

    if (data->chain != NULL) {
      chains.push_back(data->chain);
    } else if (data->map != NULL) {
      maps.push_back(data->map);
    } else if (data->mapHeader != NULL) {
      crystals.push_back(data->mapHeader);
    }
  }
  std::set<Molecule*> chainsModified;
  if (chains.size() > 0) {
    std::vector<RESIDUE*>::iterator iter = chains.begin();
    for (; iter != chains.end(); ++iter) {
      RESIDUE* chain = *iter;
      Molecule* model = chainToModel[chain];
      int chain_id = chain->chain_id();
      char chainid = (char)(chain_id&255);

      MIGetStringDialog dlg(0, "Chain ID", "Enter Chain ID (single char)");
      std::string s;
      if (dlg.GetValue(std::string(1, chainid), s) && s.size()) {
        unsigned char c = s.c_str()[0];
        model->setChainId(chain, c);
        chainsModified.insert(model);
      }
      int n;
      MIGetIntegerDialog dlg2(0, "Residue number", "Enter number for first residue (rest of chain will also be renumbered)");
      if (dlg2.GetValue(1, n)) {
        model->renumberChain(chain, n);
        chainsModified.insert(model);
      }
    }
  }
  if (crystals.size() > 0) {
    MIData data;
    MISelectCrystalDialog dlg(0, "Select Crystal");
    data["info"].str=crystals[0]->Label();
    dlg.GetResults(data);
    std::vector<CMapHeaderBase*>::iterator crystalIter = crystals.begin();
    for (; crystalIter != crystals.end(); ++crystalIter) {
      CMapHeaderBase* mapHeader = *crystalIter;
      CMapHeaderBase mh(data["info"].str);
      mapHeader->updateSymmetryAndCell(mh);

      for (int i=0; i < displaylist()->MapCount(); ++i) {
        EMapBase *map=displaylist()->GetMap(i);
        if (map->mapheader == mapHeader) {
          map->RecalcResolution();
          if (map->HasPhases() &&  map->FFTMap())
            view->doMapContour(map);
        }
      }

    }
  }
  std::set<Molecule*>::iterator iter = chainsModified.begin();
  while (iter != chainsModified.end()) {
    Molecule* model = *iter;
    ++iter;
    refreshChains(model);
  }
}

void ModelsTree::SortChains() {
  QList<QTreeWidgetItem *> selected;
  GetSelections(selected);
  if (selected.size() == 0) {
    return;
  }
  for (int i = 0; i < selected.size(); ++i) {
    QTreeWidgetItem* item = selected[i];
    if (!item) {
      continue;
    }
    TreeData* data = (TreeData*) GetItemData(item);
    if (data == NULL || !validTreeData(data)) {
      continue;
    }

    if (data->chain != NULL) {
      RESIDUE* chain = data->chain;
      Molecule* model = chainToModel[chain];
      view->SaveModelFile(model, "Sorted chains");
      model->SortChains();
      refreshChains(model);
    }
  }
}

void ModelsTree::CopyItem() {
  QList<QTreeWidgetItem *> selected;
  GetSelections(selected);
  if (selected.size() == 0) {
    return;
  }
  Application::instance()->ClearResidueBuffer();
  for (int i = 0; i < selected.size(); ++i) {
    QTreeWidgetItem* item = selected[i];
    if (!item) {
      continue;
    }
    TreeData* data = (TreeData*) GetItemData(item);
    if (data == NULL || !validTreeData(data)) {
      continue;
    }
    if (data->chain != NULL) {
      RESIDUE* chain = data->chain;
      RESIDUE* res = chain;
      while (res != NULL && res->chain_id() == chain->chain_id()) {
        Application::instance()->CopyResidueBuffer(res);
        res = res->next();
      }
    } else if (data->residue != NULL) {
      RESIDUE* residue = data->residue;
      Application::instance()->CopyResidueBuffer(residue);
    }
  }
}

void ModelsTree::InsertItem() {
  QList<QTreeWidgetItem *> selected;
  GetSelections(selected);
  if (selected.size() == 0) {
    return;
  }
  for (int i = 0; i < selected.size(); ++i) {
    QTreeWidgetItem* item = selected[i];
    if (!item) {
      continue;
    }
    TreeData* data = (TreeData*) GetItemData(item);
    if (data == NULL || !validTreeData(data)) {
      continue;
    }
    if (data->model != NULL) {
      Molecule* model = data->model;
      std::string chainId("N");
      MIGetStringDialog dlg(0, "New chain ID", "Chain ID for new chain");
      if (!dlg.GetValue(chainId, chainId) || !chainId.size()) {
        return;
      }
      unsigned char c = chainId.c_str()[0];
      RESIDUE* buffer = Application::instance()->GetResidueBuffer();
      model->InsertResidues(model->getResidues(), buffer, 3, (c&255) + 1*256);
      model->Build();
      model->SortChains();
      addChains(model);
      setFocus(Qt::MouseFocusReason);
      break;
    } else if (data->chain != NULL) {
      RESIDUE* chain = data->chain;
      Molecule* model = chainToModel[chain];
      unsigned short chain_id = chain->chain_id();
      RESIDUE* buffer = Application::instance()->GetResidueBuffer();
      model->InsertResidues(chain, buffer, 3, chain_id);
      model->Build();
      model->SortChains();
      addChains(model);
      setFocus(Qt::MouseFocusReason);
      break;
    }
  }
}

void ModelsTree::MapProperties() {
  QList<QTreeWidgetItem *> selected;
  GetSelections(selected);
  if (selected.size() == 0) {
    return;
  }
  for (int i = 0; i < selected.size(); ++i) {
    QTreeWidgetItem* item = selected[i];
    if (!item) {
      continue;
    }
    TreeData* data = (TreeData*) GetItemData(item);
    if (data == NULL || !validTreeData(data)) {
      continue;
    }
    if (data->map != NULL) {
      EMap* map = data->map;
      MIMessageBox(map->Info().c_str(), "Map Properties", MIDIALOG_ICON_INFORMATION);
    }
  }
}

void ModelsTree::MapExport() {
  QList<QTreeWidgetItem *> selected;
  GetSelections(selected);
  if (selected.size() == 0) {
    return;
  }
  for (int i = 0; i < selected.size(); ++i) {
    QTreeWidgetItem* item = selected[i];
    if (!item) {
      continue;
    }
    TreeData* data = (TreeData*) GetItemData(item);
    if (data == NULL || !validTreeData(data)) {
      continue;
    }
    if (data->map != NULL) {
      EMap* map = data->map;
      map->Export();
    }
  }
}

void ModelsTree::ModelExport() {
  QList<QTreeWidgetItem *> selected;
  GetSelections(selected);
  if (selected.size() == 0) {
    return;
  }
  for (int i = 0; i < selected.size(); ++i) {
    QTreeWidgetItem* item = selected[i];
    if (!item) {
      continue;
    }
    TreeData* data = (TreeData*) GetItemData(item);
    if (data == NULL || !validTreeData(data)) {
      continue;
    }
    if (data->model != NULL) {
      Molecule* model = data->model;
      const std::string& s = MIFileSelector("Choose a name for the model file", "", "", "pdb",
                               "PDB file (*.pdb)|*.pdb|All files (*.*)|*.*", MI_SAVE_MODE);
      if (s.size()) {
        model->SavePDBFile(s.c_str());
      }
    }
  }
}

void ModelsTree::setView(MIGLWidget* view) {
  if (view == NULL) {
    return;
  }
  this->view = view;
  connect(view, SIGNAL(focusResidueChanged(chemlib::RESIDUE*)),
          this, SLOT(focusResidueChanged(chemlib::RESIDUE*)));
  addModels(displaylist());
  addMaps(displaylist());
  solidSurfMenu->clear();
  if (view) {
    solidSurfActionGroup = new QActionGroup(this);
    connect(solidSurfActionGroup, SIGNAL(triggered(QAction*)), this, SLOT(solidSurfaceActionTriggered(QAction*)));

    solidsurf_menu_action(solidSurfMenu, solidSurfActionGroup, ID_SOLIDSURFACE_BUILD, "Build surface");
    solidsurf_menu_action(solidSurfMenu, solidSurfActionGroup, ID_SOLIDSURFACE_COLOR, "Color surface");
    solidsurf_menu_action(solidSurfMenu, solidSurfActionGroup, ID_SOLIDSURFACE_COLOR_BY_ATOM, "Color surface by atom type");
  }
  fillSolidSurfMenu(this, view, solidSurfMenu);
}

void ModelsTree::setResiduesTree(ResiduesTree* tree) {
  residuesTree = tree;
}

void ModelsTree::setAtomsTree(AtomsTree* tree) {
  atomsTree = tree;
}

RESIDUE* ModelsTree::findChain(RESIDUE* residue) {
  if (!residue)
    return NULL;

  RESIDUE* chain = NULL;
  bool found = false;
  Displaylist::ModelList::iterator modelIter = displaylist()->begin();
  while (!found && modelIter != displaylist()->end()) {
    Molecule* model = *modelIter;
    ++modelIter;
    unsigned short chain_id = 0;
    for (RESIDUE* res = model->getResidues(); res != NULL; res = res->next()) {
      if ((res->chain_id()) != chain_id) {
        chain = res;
        chain_id = chain->chain_id();
      }
      if (res == residue) {
        found = true;
        break;
      }
    }
  }
  if (!found) {
    return NULL;
  }
  return chain;
}

Molecule* ModelsTree::findModel(RESIDUE* residue) {
  Molecule* model = NULL;
  bool found = false;
  Displaylist::ModelList::iterator modelIter = displaylist()->begin();
  while (!found && modelIter != displaylist()->end()) {
    model = *modelIter;
    ++modelIter;
    for (RESIDUE* res = model->getResidues(); res != NULL; res = res->next()) {
      if (res == residue) {
        found = true;
        break;
      }
    }
  }
  if (!found) {
    return NULL;
  }
  return model;
}

void ModelsTree::focusResidueChanged(RESIDUE* residue) {
  if (syncView) {
    RESIDUE* chain = findChain(residue);
    if (chain != NULL) {
      setCurrentChain(chain);
      residuesTree->setCurrentResidue(residue);
    }
  }
}

void ModelsTree::selectionChanged(Molecule*, RESIDUE* residue, MIAtom* atom) {
  if (syncView) {
    RESIDUE* chain = findChain(residue);
    if (chain != NULL) {
      setCurrentChain(chain);
      residuesTree->setCurrentResidue(residue);
      atomsTree->setCurrentAtom(atom);
    }
  }
}

Displaylist* ModelsTree::displaylist() {
  return view->GetDisplaylist();
}

void ModelsTree::addModels(Displaylist* displaylist) {
  Displaylist::ModelList::iterator modelIter = displaylist->begin();
  while (modelIter != displaylist->end()) {
    Molecule* model = *modelIter;
    ++modelIter;
    addModel(model);
  }
  connect(displaylist, SIGNAL(modelAdded(Molecule*)),
          this, SLOT(modelAdded(Molecule*)));
  connect(displaylist, SIGNAL(currentMoleculeChanged(Molecule*,Molecule*)),
          this, SLOT(currentModelChanged(Molecule*,Molecule*)));
  connect(displaylist, SIGNAL(selectionChanged(Molecule*,chemlib::RESIDUE*,chemlib::MIAtom*)),
          this, SLOT(selectionChanged(Molecule*,chemlib::RESIDUE*,chemlib::MIAtom*)));
}

void ModelsTree::addMaps(Displaylist* displaylist) {
  Displaylist::MapList::iterator mapIter = displaylist->getMaps().begin();
  while (mapIter != displaylist->getMaps().end()) {
    EMap* map = *mapIter;
    ++mapIter;
    addMap(map);
  }
  connect(displaylist, SIGNAL(mapAdded(EMap*)),
          this, SLOT(mapAdded(EMap*)));
  connect(displaylist, SIGNAL(mapToBeDeleted(EMap*)),
          this, SLOT(mapToBeDeleted(EMap*)));
  connect(displaylist, SIGNAL(currentMapChanged(EMap*,EMap*)),
          this, SLOT(currentMapChanged(EMap*,EMap*)));
}

void ModelsTree::addModel(Molecule* model) {
  QTreeWidgetItem* previousItem = FindLastModelItem();
  TreeData* data = new TreeData;
  data->model = model;
  QTreeWidgetItem* item;
  if (previousItem) {
    item = insertItem(rootId, previousItem, model->pathname.c_str(), 0, 0, data);
  } else {
    item = prependItem(rootId, model->pathname.c_str(), 0, 0, data);
  }
  if (!item) {
    return;
  }
  modelToData[model] = data;
  stylizeModel(model);

  model->SortChains();
  addModelCrystal(model);
  addChains(model);

  connect(model, SIGNAL(moleculeToBeDeleted(chemlib::MIMoleculeBase*)),
          this, SLOT(modelToBeDeleted(chemlib::MIMoleculeBase*)));
  connect(model, SIGNAL(residuesToBeDeleted(chemlib::MIMoleculeBase*,std::vector<chemlib::RESIDUE*>&)),
          this, SLOT(residuesToBeDeleted(chemlib::MIMoleculeBase*,std::vector<chemlib::RESIDUE*>&)));
  connect(model, SIGNAL(residuesDeleted(chemlib::MIMoleculeBase*)),
          this, SLOT(residuesDeleted(chemlib::MIMoleculeBase*)));
  connect(model, SIGNAL(moleculeChanged(chemlib::MIMoleculeBase*)),
          this, SLOT(moleculeChanged(chemlib::MIMoleculeBase*)));
  connect(model, SIGNAL(atomChanged(chemlib::MIMoleculeBase*,chemlib::MIAtomList&)),
          this, SLOT(atomChanged(chemlib::MIMoleculeBase*,chemlib::MIAtomList&)));

  expandItem(item);
}

void ModelsTree::moleculeChanged(MIMoleculeBase* model) {
  Molecule* m = dynamic_cast<Molecule*>(model);
  if (m != NULL) {
    refreshChains(m);
  }
}

void ModelsTree::residuesToBeDeleted(MIMoleculeBase*, std::vector<RESIDUE*>& residues) {
  std::vector<RESIDUE*>::iterator iter;
  for (iter = residues.begin(); iter != residues.end(); ++iter) {
    RESIDUE* residue = *iter;
    if (residueToChainData.find(residue) != residueToChainData.end()) {
      if (currentChain == residue) {
        setCurrentChain(NULL);
      }
      TreeData* data = residueToChainData[residue];
      if (data == NULL || !validTreeData(data)) {
        residueToChainData.erase(residue);
        continue;
      }
      QTreeWidgetItem* item = data->GetId();
      if (item) {
        residueToChainData.erase(residue);
        chainToModel.erase(residue);
        Delete(item);
      }
    }
  }
}

void ModelsTree::residuesDeleted(MIMoleculeBase* model) {
  Molecule* m = dynamic_cast<Molecule*>(model);
  if (m != NULL) {
    refreshChains(m);
  }
}

void ModelsTree::addModelCrystal(Molecule* model) {
  TreeData* modelData = modelToData[model];
  if (modelData == NULL || !validTreeData(modelData)) {
    modelToData.erase(model);
    return;
  }
  QTreeWidgetItem* modelItem = modelData->GetId();
  if (!modelItem) {
    return;
  }
  CMapHeaderBase& mapHeader = model->GetMapHeader();
  connect(&mapHeader, SIGNAL(mapHeaderChanged(CMapHeaderBase*)),
          this, SLOT(mapHeaderChanged(CMapHeaderBase*)));
  TreeData* data = new TreeData;
  data->mapHeader = &mapHeader;
  QTreeWidgetItem* item = insertItem(modelItem, 0, mapHeader.Label().c_str(), 4, 4, data);
  if (!item) {
    return;
  }
  modelToCrystalData[model] = data;
  mapHeaderToCrystalData[&mapHeader] = data;
  stylizeCrystal(&mapHeader);
}

Molecule* ModelsTree::getCurrentModel() {
  return currentModel;
}

void ModelsTree::refreshChains(Molecule* model) {
  model->SortChains();
  // Delete all chain items for this model
  std::vector<TreeData*> chainDataForModel;
  ResidueToDataMap::iterator iter = residueToChainData.begin();
  while (iter != residueToChainData.end()) {
    RESIDUE* residue = iter->first;
    TreeData* data = iter->second;
    ++iter;
    if (validTreeData(data)) {
      if (chainToModel.find(residue) != chainToModel.end()
          && chainToModel[residue] == model) {
        chainDataForModel.push_back(data);
      }
    }
  }

  std::vector<TreeData*>::iterator iter2;
  for (iter2 = chainDataForModel.begin(); iter2 != chainDataForModel.end(); ++iter2) {
    TreeData* data = *iter2;
    QTreeWidgetItem* item = data->GetId();
    if (item && data->chain != NULL) {
      RESIDUE* chain = data->chain;
      Delete(item);
      residueToChainData.erase(chain);
      chainToModel.erase(chain);
      if (currentChain == chain) {
        setCurrentChain(NULL);
      }
    }
  }
  addChains(model);
}

void ModelsTree::addChains(Molecule* model) {
  TreeData* modelData = modelToData[model];
  if (modelData == NULL || !validTreeData(modelData)) {
    modelToData.erase(model);
    return;
  }
  QTreeWidgetItem* modelItem = modelData->GetId();
  if (!modelItem) {
    return;
  }
  std::set<TreeData*> currentChainDataSet;
  int chain_id = -9999;
  QTreeWidgetItem* previousChainItem;
  for (RESIDUE* res = model->getResidues(); res != NULL; res = res->next()) {
    if (res->chain_id() != chain_id) {
      chain_id = res->chain_id();
      if (residueToChainData.find(res) == residueToChainData.end()) {
        std::string text = chainstring(res).c_str();
        QTreeWidgetItem* chainItem;
        TreeData* data = new TreeData;
        data->chain = res;
        if (previousChainItem) {
          chainItem = insertItem(modelItem, previousChainItem, text, 5, 5, data);
        } else {
          chainItem = appendItem(modelItem, text, 5, 5, data);
        }
        if (chainItem) {
          previousChainItem = chainItem;
          residueToChainData[res] = data;
          chainToModel[res] = model;
          stylizeChain(res);
          currentChainDataSet.insert(data);
        }
      } else {
        TreeData* data = residueToChainData[res];
        if (data == NULL || !validTreeData(data)) {
          residueToChainData.erase(res);
          continue;
        }
        previousChainItem = data->GetId();
        if (previousChainItem) {
          stylizeChain(res);
          currentChainDataSet.insert(data);
        }
      }
    }
  }
  // Delete any items for chains not found
  ResidueToDataMap::iterator iter = residueToChainData.begin();
  while (iter != residueToChainData.end()) {
    TreeData* data = iter->second;
    ++iter;
    if (data == NULL || !validTreeData(data)) {
      continue;
    }
    QTreeWidgetItem* chainItem = data->GetId();
    if (chainItem && data->chain != NULL) {
      RESIDUE* chain = data->chain;
      Molecule* chainModel = chainToModel[chain];
      if (model == chainModel
          && currentChainDataSet.find(data) == currentChainDataSet.end()) {
        Delete(chainItem);
      }
    }
  }
  if (currentChain == NULL) {
    setCurrentChain(model->getResidues());
  }
}

void ModelsTree::restylize() {
  ModelToDataMap::iterator iter = modelToData.begin();
  while (iter != modelToData.end()) {
    Molecule* model = iter->first;
    TreeData* data = iter->second;
    ++iter;
    if (data == NULL || !validTreeData(data)) {
      continue;
    }
    stylizeModel(model);
  }
  ResidueToDataMap::iterator iter2 = residueToChainData.begin();
  while (iter2 != residueToChainData.end()) {
    RESIDUE* chain = iter2->first;
    TreeData* data = iter2->second;
    ++iter2;
    if (data == NULL || !validTreeData(data)) {
      continue;
    }
    stylizeChain(chain);
  }
  MapToDataMap::iterator iter3 = mapToData.begin();
  while (iter3 != mapToData.end()) {
    EMap* map = iter3->first;
    TreeData* data = iter3->second;
    ++iter3;
    if (data == NULL || !validTreeData(data)) {
      continue;
    }
    stylizeMap(map);
  }
  MapHeaderToDataMap::iterator iter4 = mapHeaderToCrystalData.begin();
  while (iter4 != mapHeaderToCrystalData.end()) {
    CMapHeaderBase* mapHeader = iter4->first;
    TreeData* data = iter4->second;
    ++iter4;
    if (data == NULL || !validTreeData(data)) {
      continue;
    }
    stylizeCrystal(mapHeader);
  }
}

void ModelsTree::setTruncateWidth(int width) {
  truncateWidth = width - 64;
}

std::string ModelsTree::truncateLeft(const std::string& name, int shorter) {
  if (!name.size())
    return name;

  QFontMetrics fm = fontMetrics();
  size_t width = truncateWidth - shorter;
  std::string prefix("...");
  std::string result = name;
  size_t stringWidth = fm.width(result.c_str());
  for (size_t i = 1; stringWidth > width && i < name.size(); ++i) {
    result = prefix + &name[i];
    stringWidth = fm.width(result.c_str());
  }

  return result;
}

std::string ModelsTree::truncateRight(const std::string& name, int shorter) {
  if (!name.size())
    return name;

  QFontMetrics fm = fontMetrics();
  size_t width = truncateWidth - shorter;
  std::string suffix("...");
  std::string shorterName(name);
  std::string result(name);
  size_t stringWidth = fm.width(result.c_str());
  for (size_t i = 1; stringWidth > width && i < name.size(); ++i) {
    shorterName.resize(name.size() - i);
    result = shorterName + suffix;
    stringWidth = fm.width(result.c_str());
  }

  return result;
}

void ModelsTree::stylizeModel(Molecule* model) {
  std::string text = truncateLeft(model->pathname.c_str());
  TreeData* data = modelToData[model];
  if (data == NULL || !validTreeData(data)) {
    modelToData.erase(model);
    return;
  }
  QTreeWidgetItem* item = data->GetId();
  if (!item) {
    return;
  }
  item->setText(0, text.c_str());

}

void ModelsTree::stylizeMap(EMap* map) {
  TreeData* data = mapToData[map];
  if (data == NULL || !validTreeData(data)) {
    mapToData.erase(map);
    return;
  }
  QTreeWidgetItem* item = data->GetId();
  if (!item) {
    return;
  }

  std::string mapName = map->mapName.c_str();
  int mapType = map->mapheader->maptype;
  if (mapType == (int)MIMapType::DirectFFT) {
    if (map->fColumnName.size() > 0) {
      mapName += ": ";
      mapName += map->fColumnName.c_str();
    }
  } else {
    mapName += ": ";
    mapName += StringForMapType(mapType);
  }

  std::string text = truncateLeft(mapName.c_str());
  item->setText(0, text.c_str());
}

void ModelsTree::atomChanged(MIMoleculeBase* model, MIAtomList& atoms) {
  std::set<RESIDUE*> chains;
  MIAtom_iter iter = atoms.begin();
  while (iter != atoms.end()) {
    MIAtom* atom = *iter;
    ++iter;
    RESIDUE* res = residue_from_atom(model->getResidues(), atom);
    if (res != NULL) {
      RESIDUE* chain = findChain(res);
      if (chain != NULL) {
        chains.insert(chain);
      }
    }
  }
  std::set<RESIDUE*>::iterator iter2 = chains.begin();
  while (iter2 != chains.end()) {
    RESIDUE* chain = *iter2;
    ++iter2;
    stylizeChain(chain);
  }
}

void ModelsTree::stylizeChain(RESIDUE* residue) {
  if (residue == NULL || residueToChainData.find(residue) == residueToChainData.end()) {
    return;
  }
  TreeData* data = residueToChainData[residue];
  if (data == NULL || !validTreeData(data)) {
    residueToChainData.erase(residue);
    return;
  }
  QTreeWidgetItem* item = data->GetId();
  if (!item) {
    return;
  }

  bool partial = false;
  bool hidden = true;
  int shownAtoms = 0;
  int hiddenAtoms = 0;
  RESIDUE* res = residue;
  int chainId = res->chain_id();
  while (res != NULL && chainId == res->chain_id()) {
    for (int ia = 0; ia < res->atomCount(); ++ia) {
      MIAtom* atom = res->atom(ia);
      if (atom->isHidden()) {
        ++hiddenAtoms;
        if (shownAtoms > 0) {
          hidden = false;
          partial = true;
          break;
        }
      } else {
        ++shownAtoms;
        hidden = false;
        if (hiddenAtoms > 0) {
          partial = true;
          break;
        }
      }
    }
    res = res->next();
  }

  int image = 5;
  if (residue == currentChain) {
    if (hidden) {
      image = 10;
    } else if (partial) {
      image = 8;
    } else {
      image = 6;
    }
  } else {
    if (hidden) {
      image = 9;
    } else if (partial) {
      image = 7;
    } else {
      image = 5;
    }
  }

  std::string name = chainstring(residue).c_str();
  std::string text = truncateRight(name, 24);
  item->setText(0, text.c_str());
  item->setIcon(0,_imageList[ image]);
  item->setIcon(0,_imageList[ image ]);
}

void ModelsTree::addMap(EMap* map) {
  TreeData* data = new TreeData;
  data->map = map;
  QTreeWidgetItem* item = appendItem(rootId, map->mapName.c_str(), 2, 2, data);
  if (!item) {
    return;
  }
  mapToData[map] = data;
  connect(map, SIGNAL(mapFftRecalculated(EMap*)),
          this, SLOT(mapFftRecalculated(EMap*)));
  stylizeMap(map);

  addMapCrystal(map);

  expandItem(item);
}

void ModelsTree::addMapCrystal(EMap* map) {
  CMapHeaderBase* mapHeader = map->GetMapHeader();
  if (mapHeader != NULL) {
    connect(mapHeader, SIGNAL(mapHeaderChanged(CMapHeaderBase*)),
            this, SLOT(mapHeaderChanged(CMapHeaderBase*)));
    TreeData* mapData = mapToData[map];
    if (mapData == NULL || !validTreeData(mapData)) {
      mapToData.erase(map);
      return;
    }
    QTreeWidgetItem* mapItem = mapData->GetId();
    if (!mapItem) {
      return;
    }
    TreeData* data = new TreeData;
    data->mapHeader = mapHeader;
    QTreeWidgetItem* item = appendItem(mapItem, mapHeader->Label().c_str(), 4, 4, data);
    if (!item) {
      return;
    }
    mapToCrystalData[map] = data;
    mapHeaderToCrystalData[mapHeader] = data;
    stylizeCrystal(mapHeader);
  }
}

void ModelsTree::stylizeCrystal(CMapHeaderBase* mapHeader) {
  std::string text = mapHeader->Label().c_str();
  text = truncateRight(text, 16);
  TreeData* data = mapHeaderToCrystalData[mapHeader];
  if (data == NULL || !validTreeData(data)) {
    mapHeaderToCrystalData.erase(mapHeader);
    return;
  }
  QTreeWidgetItem* item = data->GetId();
  if (!item) {
    return;
  }
  item->setText(0, text.c_str());
}

void ModelsTree::mapHeaderChanged(CMapHeaderBase* mapHeader) {
  if (mapHeaderToCrystalData.find(mapHeader) != mapHeaderToCrystalData.end()) {
    stylizeCrystal(mapHeader);
  }
}

void ModelsTree::mapFftRecalculated(EMap* map) {
  stylizeMap(map);
}

void ModelsTree::modelAdded(Molecule* model) {
  addModel(model);
  update();
}

void ModelsTree::modelToBeDeleted(MIMoleculeBase* mol) {
  Molecule *model=(Molecule*)mol;
  TreeData* data = modelToData[model];
  if (!validTreeData(data)) {
    return;
  }
  QTreeWidgetItem* item = data->GetId();
  if (!item) {
    return;
  }
  modelToData.erase(model);
  if (model == currentModel) {
    setCurrentChain(NULL);
  }
  for (RESIDUE* res = model->getResidues(); RESIDUE::isValid(res); res = res->next()) {
    residueToChainData.erase(res);
  }
  mapHeaderToCrystalData.erase(&(model->GetMapHeader()));
  Delete(item);
  update();
}

void ModelsTree::currentModelChanged(Molecule* oldModel, Molecule* newModel) {
  TreeData* data;
  QTreeWidgetItem* item;
  if (oldModel != NULL) {
    data = modelToData[oldModel];
    if (data != NULL) {
      if (!validTreeData(data)) {
        modelToData.erase(oldModel);
      } else {
        item = data->GetId();
        if (item) {
          item->setIcon(0,_imageList[ 0]);
          item->setIcon(0,_imageList[ 0 ]);
        }
      }
    }
  }
  if (newModel != NULL) {
    data = modelToData[newModel];
    if (data != NULL) {
      if (!validTreeData(data)) {
        modelToData.erase(oldModel);
      } else {
        item = data->GetId();
        if (item) {
          item->setIcon(0,_imageList[ 1]);
          item->setIcon(0,_imageList[ 1 ]);
        }
      }
    }
  }
}

void ModelsTree::mapAdded(EMap* map) {
  addMap(map);
  update();
}

void ModelsTree::mapToBeDeleted(EMap* map) {
  TreeData* data = mapToData[map];
  if (data == NULL) {
    return;
  }
  if (!validTreeData(data)) {
    mapToData.erase(map);
    mapToCrystalData.erase(map);
    mapHeaderToCrystalData.erase(map->GetMapHeader());
    return;
  }
  QTreeWidgetItem* item = data->GetId();
  if (!item) {
    return;
  }
  Delete(item);
  mapToData.erase(map);
  mapToCrystalData.erase(map);
  mapHeaderToCrystalData.erase(map->GetMapHeader());
  update();
}

void ModelsTree::currentMapChanged(EMap* oldMap, EMap* newMap) {
  TreeData* data;
  QTreeWidgetItem* item;
  if (oldMap != NULL) {
    data = mapToData[oldMap];
    if (data == NULL || !validTreeData(data)) {
      mapToData.erase(oldMap);
    } else {
      item = data->GetId();
      if (item) {
        item->setIcon(0,_imageList[ 2]);
        item->setIcon(0,_imageList[ 2 ]);
      }
    }
    stylizeMap(oldMap);
  }
  if (newMap != NULL) {
    data = mapToData[newMap];
    if (data == NULL || !validTreeData(data)) {
      mapToData.erase(newMap);
    } else {
      item = data->GetId();
      if (item) {
        item->setIcon(0,_imageList[ 3]);
        item->setIcon(0,_imageList[ 3 ]);
      }
    }
    stylizeMap(newMap);
  }
}

bool ModelsTree::colorDialog(int& color, int& colorMethod) {
  MIGenericDialog dlg(this, "Choose Color");
  MIData data;
  data["Color:"].color[0] = (unsigned char)view->WhenShownColor;
  data["Color:"].isColorIndex = true;

  data["Method:"].radio = view->WhenShownColorMethod;
  data["Method:"].radio_count = 7;
  data["Method:"].radio_labels.push_back("Carbon Only");
  data["Method:"].radio_labels.push_back("All Atoms");
  data["Method:"].radio_labels.push_back("Secondary Structure");
  data["Method:"].radio_labels.push_back("B-Value");
  data["Method:"].radio_labels.push_back("Atom Type");
  data["Method:"].radio_labels.push_back("Hydrophobicity");
  data["Method:"].radio_labels.push_back("Shapley");

  if (!dlg.GetResults(data)) {
    return false;
  }
  color = (int)data["Color:"].color[0];
  colorMethod = (int)data["Method:"].radio;
  return true;
}


ModelsView::LineEditToModelsTreeMap ModelsView::lineEditToModelsTree;
ModelsView::PanelToLineEditMap ModelsView::panelToLineEdit;
ModelsView::PanelToToolButtonMap ModelsView::panelToToolButton;
ModelsView::ButtonCtrlToModelsTreeMap ModelsView::buttonCtrlToModelsTree;
ModelsView::ButtonCtrlToLineEditMap ModelsView::buttonCtrlToLineEdit;
ModelsView::PanelToButtonCtrlMap ModelsView::panelToButtonCtrl;

class ModelsViewPanel : public QWidget {
public:

  ModelsTree* modelsTree;
  ResiduesTree* residuesTree;
  AtomsTree* atomsTree;

  ModelsViewPanel(QWidget* parent);
  virtual ~ModelsViewPanel();
};

ModelsViewPanel::ModelsViewPanel(QWidget* parent)
  : QWidget(parent), modelsTree(NULL), residuesTree(NULL), atomsTree(NULL) {
}

ModelsViewPanel::~ModelsViewPanel() {
  //printf("In ModelsViewPanel dtor\n");
}

ModelsView::ModelsView(QWidget* parent) : ViewSyncedPanel(parent) {
}

ModelsView::~ModelsView() {
  //printf("In ModelsView dtor\n");
}


class MyLineEdit : public QLineEdit {

Q_OBJECT
public:
    MyLineEdit(const QString &str, QWidget *parent) : QLineEdit(str, parent) {}

    void focusInEvent(QFocusEvent *evt) {
      if (text() == QString("Go to residue"))
        clear();
      selectAll();
      QLineEdit::focusInEvent(evt);
    }

    void enterEvent(QEvent *evt) {
      setFocus(Qt::MouseFocusReason);
      QLineEdit::enterEvent(evt);
    }

};


QWidget* ModelsView::createPanelForView(MIGLWidget* view, QWidget* parent) {
  parent->layout()->setContentsMargins(0, 0, 0, 0);

  ModelsViewPanel* panel = new ModelsViewPanel(parent);

  // create vertical layout
  QHBoxLayout *hlayout=new QHBoxLayout();
  hlayout->setContentsMargins(0, 0, 0, 0);
  hlayout->setSpacing(2);

  QLineEdit* goToResidueLineEdit = new MyLineEdit("Go to residue", panel);
  goToResidueLineEdit->setToolTip("Go to residue (A 123)");

  panelToLineEdit[panel] = goToResidueLineEdit;

  QPushButton* gotoButton = new QPushButton("GoTo",panel);
  panelToButtonCtrl[panel] = gotoButton;

  QToolButton* toolButton = new QToolButton(panel);
  toolButton->setToolTip("Sync with selection");
  toolButton->setIcon(QIcon(QPixmap(synced_xpm)));
  toolButton->setCheckable(true);
  toolButton->setChecked(MIConfig::Instance()->GetProfileInt("ModelsView", "syncView", 1) != 0);
  connect(toolButton, SIGNAL(toggled(bool)), this, SLOT(updateSyncView(bool)));

  panelToToolButton[panel] = toolButton;
  updateSyncView(toolButton);

  QSpacerItem *spacer=new QSpacerItem(20,5,QSizePolicy::Expanding);
  hlayout->addItem(spacer);
  hlayout->addWidget(goToResidueLineEdit);
  hlayout->addWidget(gotoButton);
  hlayout->addWidget(toolButton);

  QVBoxLayout *vbox=new QVBoxLayout(panel);
  vbox->addLayout(hlayout);
  vbox->setContentsMargins(0, 0, 0, 0);
  vbox->setSpacing(2);

  QSplitter *splitter=new QSplitter(Qt::Vertical, panel);
  splitter->setContentsMargins(0, 0, 0, 0);
  vbox->addWidget(splitter);

  panel->setLayout(vbox);
  ModelsTree* modelsTree = new ModelsTree(splitter);
  ResiduesTree* residuesTree = new ResiduesTree(splitter);
  AtomsTree* atomsTree = new AtomsTree(splitter);

  modelsTree->setView(view);
  modelsTree->setTruncateWidth(width());
  residuesTree->setView(view);
  atomsTree->setView(view);

  splitter->addWidget(modelsTree);
  splitter->addWidget(residuesTree);
  splitter->addWidget(atomsTree);

  connect(splitter, SIGNAL(splitterMoved(int, int)),
          this, SLOT(OnSplitterChanged(int, int)));

  QSettings *settings=MIGetQSettings(); // could use MIConfig (it's the same file), but this api is easier here
  if (settings->contains("ModelsViewSplitter"))
    splitter->restoreState(settings->value("ModelsViewSplitter",splitter->saveState()).toByteArray());

  connect(gotoButton, SIGNAL(clicked()), this, SLOT(OnGoToResidueReturnPressed()));
  connect(goToResidueLineEdit, SIGNAL(returnPressed()), this, SLOT(OnGoToResidueReturnPressed()));

  lineEditToModelsTree[goToResidueLineEdit] = modelsTree;
  buttonCtrlToModelsTree[gotoButton] = modelsTree;
  buttonCtrlToLineEdit[gotoButton] = goToResidueLineEdit;

  modelsTree->setResiduesTree(residuesTree);
  modelsTree->setAtomsTree(atomsTree);
  residuesTree->setAtomsTree(atomsTree);

  panel->modelsTree = modelsTree;
  panel->residuesTree = residuesTree;
  panel->atomsTree = atomsTree;

  return panel;
}

void ModelsView::destroyContentsForView(MIGLWidget*, QWidget* panel) {
  if (panelToToolButton.find(panel) != panelToToolButton.end()) {
    panelToToolButton.erase(panel);
  }
  if (panelToLineEdit.find(panel) != panelToLineEdit.end()) {
    QLineEdit* lineEdit = panelToLineEdit[panel];
    panelToLineEdit.erase(panel);
    lineEditToModelsTree.erase(lineEdit);
  }
  if (panelToButtonCtrl.find(panel) != panelToButtonCtrl.end()) {
    QPushButton* buttonCtrl = panelToButtonCtrl[panel];
    panelToButtonCtrl.erase(panel);
    buttonCtrlToModelsTree.erase(buttonCtrl);
    buttonCtrlToLineEdit.erase(buttonCtrl);
  }
}

void ModelsView::updateSyncView(bool state) {
  syncView = state;
  PanelToToolButtonMap::iterator iter = panelToToolButton.begin();
  while (iter != panelToToolButton.end()) {
    QToolButton* tb = iter->second;
    ++iter;
    tb->setChecked(syncView);
  }
  MIConfig::Instance()->WriteProfileInt("ModelsView", "syncView", syncView);
}

void ModelsView::OnSplitterChanged(int, int) {
  QSplitter *splitter=dynamic_cast<QSplitter*>(sender());
  if (!splitter) {
    //printf("Splitter changed signal not from splitter\n");
    return;
  }
  QSettings *settings=MIGetQSettings(); // could use MIConfig (it's the same file), but this api is easier here
  settings->setValue("ModelsViewSplitter",splitter->saveState());
}

void ModelsView::resizeEvent(QResizeEvent *event) {
  LineEditToModelsTreeMap::iterator iter;
  for (iter = lineEditToModelsTree.begin(); iter != lineEditToModelsTree.end(); ++iter) {
    ModelsTree* modelsTree = iter->second;
    modelsTree->setTruncateWidth(event->size().width());
    modelsTree->restylize();
  }
  ViewSyncedPanel::resizeEvent(event);
}

void ModelsView::OnGoToResidueReturnPressed() {
  QLineEdit* lineEdit = dynamic_cast<QLineEdit*>(sender());
  if (lineEdit) {
    lineEditToModelsTree[lineEdit]->goToResidue(lineEdit->text().toStdString());
    lineEdit->selectAll();
    return;
  }
  QPushButton* buttonCtrl = dynamic_cast<QPushButton*>(sender());
  if (buttonCtrl) {
    QLineEdit *lineEdit=buttonCtrlToLineEdit[buttonCtrl];
    buttonCtrlToModelsTree[buttonCtrl]->goToResidue(lineEdit->text().toStdString());
    lineEdit->selectAll();
    return;
  }
}

ModelsTree* ModelsView::GetCurrentModelsTree() const {
  ModelsViewPanel* panel = dynamic_cast<ModelsViewPanel*>(stackedLayout->currentWidget());
  if (panel == NULL) {
    return NULL;
  }
  return panel->modelsTree;
}

ResiduesTree* ModelsView::GetCurrentResiduesTree() const {
  ModelsViewPanel* panel = dynamic_cast<ModelsViewPanel*>(stackedLayout->currentWidget());
  if (panel == NULL) {
    return NULL;
  }
  return panel->residuesTree;
}

AtomsTree* ModelsView::GetCurrentAtomsTree() const {
  ModelsViewPanel* panel = dynamic_cast<ModelsViewPanel*>(stackedLayout->currentWidget());
  if (panel == NULL) {
    return NULL;
  }
  return panel->atomsTree;
}

#include "ModelsView.moc"
