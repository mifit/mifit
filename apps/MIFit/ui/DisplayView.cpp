#include <QApplication>
#include <QVBoxLayout>

#include <set>

#include "core/corelib.h"
#include <chemlib/chemlib.h>
#include <chemlib/RESIDUE_.h>

#include "ui/MIDialog.h"

#include "Displaylist.h"
#include "DisplayView.h"
#include "MIQTreeWidget.h"
#include "TreeData.h"
#include "MIEventHandler.h"
#include "MIEventHandlerMacros.h"
#include "MIGLWidget.h"
#include "MIMenu.h"
#include "id.h"

#include <images/model.xpm>
#include <images/modelSelected.xpm>
#include <images/annot.xpm>
#include <images/annotSelected.xpm>
#include <images/annotHidden.xpm>
#include <images/annotHiddenSelected.xpm>
#include <images/label.xpm>
#include <images/labelSelected.xpm>
#include <images/labelHidden.xpm>
#include <images/labelHiddenSelected.xpm>
#include <images/surf.xpm>
#include <images/surfHidden.xpm>
#include <images/surfGreyed.xpm>

const int ID_DISPLAYVIEW_IMPORTERRORS = 10000;
const int ID_DISPLAYVIEW_AUTOSHOW_IMPORTED_ERRORS = 10001;

using namespace chemlib;

class DisplayTree : public MIQTreeWidget, public MIEventHandler {
  Q_OBJECT

  MIGLWidget* view;

  Displaylist* displaylist();

public:
  DisplayTree(QWidget* parent);
  virtual ~DisplayTree();
  void setView(MIGLWidget* view);

private Q_SLOTS:
  void OnItemClicked(QTreeWidgetItem *item, int column); // single click
  void OnItemActivated(QTreeWidgetItem *item, int column); // double click
  void OnItemPressed(QTreeWidgetItem *item, int column); // possible right click

  void DeleteItem();
  void ShowItem();
  void EditItem();
  void ImportErrors();
  void AutoShowErrors();

private:
  void addModels(Displaylist* displaylist);
  void addModel(Molecule* model);
  void addAnnotations(Molecule* model);
  void addAnnotation(Molecule* model, Annotation* annotation);
  void addAtomLabels(Molecule* model);
  void addAtomLabel(Molecule* model, ATOMLABEL* label);
  void addSurface(Molecule* model);

private Q_SLOTS:
  void modelAdded(Molecule* model);
  void modelToBeDeleted(chemlib::MIMoleculeBase* model);
  void currentModelChanged(Molecule* oldModel, Molecule* newModel);
  void annotationAdded(Molecule* model, Annotation* annotation);
  void annotationToBeDeleted(Molecule* model, Annotation* annotation);
  void annotationChanged(Annotation* annotation);
  void atomLabelAdded(Molecule* model, ATOMLABEL* label);
  void atomLabelToBeDeleted(Molecule* model, Molecule::AtomLabelList labels);
  void atomLabelChanged(Molecule* model, ATOMLABEL* label);
  void surfaceChanged(Molecule* model);

private:
  void stylizeItem(QTreeWidgetItem* item, Annotation* annotation);
  void stylizeItem(QTreeWidgetItem* item, ATOMLABEL* label);
  void stylizeSurfaceItem(QTreeWidgetItem* item, Molecule* model);

  QTreeWidgetItem* rootId;
  QTreeWidgetItem* annotationListItem;
  QTreeWidgetItem* atomLabelListItem;

  typedef std::map<Molecule*, TreeData*> ModelToDataMap;
  ModelToDataMap modelToData;
  ModelToDataMap modelToAnnotationListData;
  ModelToDataMap modelToAtomLabelListData;
  ModelToDataMap modelToSurfaceData;

  typedef std::map<Annotation*, TreeData*> AnnotationToDataMap;
  AnnotationToDataMap annotationToData;

  typedef std::map<ATOMLABEL*, TreeData*> AtomLabelToDataMap;
  AtomLabelToDataMap atomLabelToData;

};



DisplayTree::DisplayTree(QWidget* parent) : MIQTreeWidget(parent), MIEventHandler(this), view(NULL) {

  setSelectionMode(QAbstractItemView::ExtendedSelection);


  std::vector<QIcon> imageList;
  QIcon modelImage=QIcon(QPixmap(model_xpm));
  imageList.push_back(modelImage);
  QIcon modelSelectedImage=QIcon(QPixmap(modelSelected_xpm));
  imageList.push_back(modelSelectedImage);
  QIcon annotationImage=QIcon(QPixmap(annot_xpm));
  imageList.push_back(annotationImage);
  QIcon annotationSelectedImage=QIcon(QPixmap(annotSelected_xpm));
  imageList.push_back(annotationSelectedImage);
  QIcon annotationHiddenImage=QIcon(QPixmap(annotHidden_xpm));
  imageList.push_back(annotationHiddenImage);
  QIcon annotationHiddenSelectedImage=QIcon(QPixmap(annotHiddenSelected_xpm));
  imageList.push_back(annotationHiddenSelectedImage);
  QIcon labelImage=QIcon(QPixmap(label_xpm));
  imageList.push_back(labelImage);
  QIcon labelSelectedImage=QIcon(QPixmap(labelSelected_xpm));
  imageList.push_back(labelSelectedImage);
  QIcon labelHiddenImage=QIcon(QPixmap(labelHidden_xpm));
  imageList.push_back(labelHiddenImage);
  QIcon labelHiddenSelectedImage=QIcon(QPixmap(labelHiddenSelected_xpm));
  imageList.push_back(labelHiddenSelectedImage);
  QIcon surfaceImage=QIcon(QPixmap(surf_xpm));
  imageList.push_back(surfaceImage);
  QIcon surfaceHiddenImage=QIcon(QPixmap(surfHidden_xpm));
  imageList.push_back(surfaceHiddenImage);
  QIcon surfaceGreyedImage=QIcon(QPixmap(surfGreyed_xpm));
  imageList.push_back(surfaceGreyedImage);

  AssignImageList(imageList);

  std::string rootText = std::string("Display List");
  setHeaderLabel(rootText.c_str());
  rootId = invisibleRootItem();

  connect(this, SIGNAL(itemPressed(QTreeWidgetItem *, int)),
          this, SLOT(OnItemPressed(QTreeWidgetItem *, int)));

  connect(this, SIGNAL(itemClicked(QTreeWidgetItem *, int)),
          this, SLOT(OnItemClicked(QTreeWidgetItem *, int)));

  connect(this, SIGNAL(itemActivated(QTreeWidgetItem *, int)),
          this, SLOT(OnItemActivated(QTreeWidgetItem *, int)));

BEGIN_EVENT_TABLE(this, none)
EVT_MENU(ID_DISPLAYVIEW_DELETE, DisplayTree::DeleteItem)
EVT_MENU(ID_DISPLAYVIEW_SHOW, DisplayTree::ShowItem)
EVT_MENU(ID_DISPLAYVIEW_EDIT, DisplayTree::EditItem)
EVT_MENU(ID_DISPLAYVIEW_IMPORTERRORS, DisplayTree::ImportErrors)
EVT_MENU(ID_DISPLAYVIEW_AUTOSHOW_IMPORTED_ERRORS, DisplayTree::AutoShowErrors)
END_EVENT_TABLE()
}

DisplayTree::~DisplayTree() {
}

void DisplayTree::OnItemClicked(QTreeWidgetItem *, int) {
}

void DisplayTree::OnItemActivated(QTreeWidgetItem *item, int) {
  TreeData* data = (TreeData*) GetItemData(item);
  if (data == NULL || !validTreeData(data)) {
    return;
  }

  if (data->annotation != NULL) {
    Annotation* annotation = data->annotation;
    view->moveTo(annotation->GetX(), annotation->GetY(), annotation->GetZ());
  } else if (data->atomLabel != NULL) {
    ATOMLABEL* label = data->atomLabel;
    view->moveTo(label->atom()->x(), label->atom()->y(), label->atom()->z());
  }
}

void DisplayTree::OnItemPressed(QTreeWidgetItem *id, int) {
  if (QApplication::mouseButtons() != Qt::RightButton)
    return;

  QList<QTreeWidgetItem*> selected;
  GetSelections(selected);
  if (selected.size() <= 1) {
    //set selected item to right-clicked item
    clearSelection();
    id->setSelected(true);
  }

  bool hasAnnotationList = false;
  bool hasAtomLabelList = false;
  bool hasSurface = false;
  bool hasAnnotation = false;
  bool hasAtomLabel = false;
  Molecule* surfaceModel = NULL;
  for (int i = 0; i < selected.size(); ++i) {
    QTreeWidgetItem* item = selected[i];
    TreeData* data = (TreeData*) GetItemData(item);
    if (data == NULL || !validTreeData(data)) {
      continue;
    }
    if (data->annotationList != NULL) {
      hasAnnotationList = true;
    }
    if (data->atomLabelList != NULL) {
      hasAtomLabelList = true;
    }
    if (data->surface != NULL) {
      hasSurface = true;
      surfaceModel = data->surface;
    }
    if (data->annotation != NULL) {
      hasAnnotation = true;
    }
    if (data->atomLabel != NULL) {
      hasAtomLabel = true;
    }
  }

  MIMenu* menu = new MIMenu(*this);
  menu->Append(ID_DISPLAYVIEW_SHOW, "Show/Hide", "Show or hide this display item", false);
  menu->Append(ID_DISPLAYVIEW_EDIT, "Edit", "Edit this display item", false);
  menu->Append(ID_DISPLAYVIEW_DELETE, "Delete", "Delete this display item", false);
  menu->Append(ID_DISPLAYVIEW_IMPORTERRORS, "Import from error list", "Import from error list", false);
  if (!(hasAnnotationList || hasAtomLabelList || hasAnnotation || hasAtomLabel)) {
    menu->Enable(ID_DISPLAYVIEW_DELETE,  false);
  }
  if (selected.size() > 1) {
    menu->Enable(ID_DISPLAYVIEW_EDIT, false);
  } else if (hasAnnotationList || hasAtomLabelList || hasSurface) {
    menu->Enable(ID_DISPLAYVIEW_EDIT, false);
  }
  if (hasSurface && Molecule::isValid(surfaceModel) && surfaceModel->getDots().size() == 0) {
    menu->Enable(ID_DISPLAYVIEW_SHOW, false);
    menu->Enable(ID_DISPLAYVIEW_EDIT, false);
    menu->Enable(ID_DISPLAYVIEW_DELETE,  false);
  }

  menu->AppendCheckItem(ID_DISPLAYVIEW_AUTOSHOW_IMPORTED_ERRORS, "Auto-show imported errors", "Automatically show errors imported from PDB file");
  menu->Check(ID_DISPLAYVIEW_AUTOSHOW_IMPORTED_ERRORS,
              MIConfig::Instance()->GetProfileInt("DisplayView", "autoShowError", 1) != 0);

  QPoint pos=QCursor::pos();
  menu->doExec(pos);
  delete menu;
}

void DisplayTree::DeleteItem() {
  QList<QTreeWidgetItem*> selected;
  GetSelections(selected);
  if (selected.size() == 0) {
    return;
  }

  for (int i = 0; i < selected.size(); ++i) {
    QTreeWidgetItem* item = selected[i];
    TreeData* data = (TreeData*) GetItemData(item);
    if (data == NULL || !validTreeData(data)) {
      continue;
    }
    if (Molecule::isValid(data->annotationList)) {
      Molecule* model = data->annotationList;
      model->clearAnnotations();
    } else if (Molecule::isValid(data->atomLabelList)) {
      Molecule* model = data->atomLabelList;
      model->clearAtomLabels();
    } else if (Molecule::isValid(data->surface)) {
      data->surface->FreeDots();
    } else if (data->annotation != NULL && Molecule::isValid(data->annotationModel)) {
      Annotation* annotation = data->annotation;
      Molecule* model = data->annotationModel;
      model->deleteAnnotation(annotation);
    } else if (data->atomLabel != NULL && Molecule::isValid(data->atomLabelModel)) {
      ATOMLABEL* label = data->atomLabel;
      Molecule* model = data->atomLabelModel;
      model->deleteAtomLabel(label);
    }
  }
}

void DisplayTree::ShowItem() {
  QList<QTreeWidgetItem*> selected;
  GetSelections(selected);
  if (selected.size() == 0) {
    return;
  }
  std::set<Annotation*> annotations;
  std::set<std::pair<ATOMLABEL*, Molecule*> > labels;
  for (int i = 0; i < selected.size(); ++i) {
    QTreeWidgetItem* item = selected[i];
    TreeData* data = (TreeData*) GetItemData(item);
    if (data == NULL || !validTreeData(data)) {
      continue;
    }
    if (Molecule::isValid(data->annotationList)) {
      Molecule* model = data->annotationList;
      Molecule::AnnotationList::iterator iter;
      for (iter = model->getAnnotations().begin(); iter != model->getAnnotations().end(); ++iter) {
        Annotation* annotation = *iter;
        annotations.insert(annotation);
      }
    } else if (Molecule::isValid(data->atomLabelList)) {
      Molecule* model = data->atomLabelList;
      Molecule::AtomLabelList::iterator iter;
      for (iter = model->getAtomLabels().begin(); iter != model->getAtomLabels().end(); ++iter) {
        ATOMLABEL* label = *iter;
        labels.insert(std::make_pair(label, model));
      }
    } else if (Molecule::isValid(data->surface)) {
      Molecule* model = data->surface;
      if (model->DotsVisible()) {
        model->HideDots();
      } else {
        model->ShowDots();
      }
    } else if (data->annotation != NULL && Molecule::isValid(data->annotationModel)) {
      Annotation* annotation = data->annotation;
      annotations.insert(annotation);
    } else if (data->atomLabel != NULL && Molecule::isValid(data->atomLabelModel)) {
      ATOMLABEL* label = data->atomLabel;
      Molecule* model = data->atomLabelModel;
      labels.insert(std::make_pair(label, model));
    }
  }

  std::set<Annotation*>::iterator iter;
  for (iter = annotations.begin(); iter != annotations.end(); ++iter) {
    Annotation* annotation = *iter;
    annotation->setHidden(!annotation->isHidden());
  }

  std::set<std::pair<ATOMLABEL*, Molecule*> >::iterator iter2;
  for (iter2 = labels.begin(); iter2 != labels.end(); ++iter2) {
    ATOMLABEL* label = iter2->first;
    Molecule* model = iter2->second;
    model->setAtomLabelVisible(label, !label->isVisible());
  }
}

void DisplayTree::EditItem() {
  QList<QTreeWidgetItem*> selected;
  GetSelections(selected);
  if (selected.size() == 0) {
    return;
  }

  std::set<Annotation*> annotations;
  std::set<std::pair<ATOMLABEL*, Molecule*> > labels;
  std::set<Molecule*> surfaces;
  for (int i = 0; i < selected.size(); ++i) {
    QTreeWidgetItem* item = selected[i];
    TreeData* data = (TreeData*) GetItemData(item);
    if (data == NULL || !validTreeData(data)) {
      continue;
    }
    if (Molecule::isValid(data->annotationList)) {
      Molecule* model = data->annotationList;
      Molecule::AnnotationList::iterator iter;
      for (iter = model->getAnnotations().begin(); iter != model->getAnnotations().end(); ++iter) {
        Annotation* annotation = *iter;
        annotations.insert(annotation);
      }
    } else if (Molecule::isValid(data->atomLabelList)) {
      Molecule* model = data->atomLabelList;
      Molecule::AtomLabelList::iterator iter;
      for (iter = model->getAtomLabels().begin(); iter != model->getAtomLabels().end(); ++iter) {
        ATOMLABEL* label = *iter;
        labels.insert(std::make_pair(label, model));
      }
    } else if (Molecule::isValid(data->surface) && data->surface->getDots().size() > 0) {
      Molecule* model = data->surface;
      surfaces.insert(model);
    } else if (data->annotation != NULL && Molecule::isValid(data->annotationModel)) {
      Annotation* annotation = data->annotation;
      annotations.insert(annotation);
    } else if (data->atomLabel != NULL && Molecule::isValid(data->atomLabelModel)) {
      ATOMLABEL* label = data->atomLabel;
      Molecule* model = data->atomLabelModel;
      labels.insert(std::make_pair(label, model));
    }
  }

  if (annotations.size() > 0) {
    Annotation* annotation = *annotations.begin();
    MIGenericDialog dlg(NULL, "Edit Annotation");
    MIData data;
    data["Text:"].str = std::string(annotation->GetText());
    data["Color:"].color[0] = annotation->m_color.red;
    data["Color:"].color[1] = annotation->m_color.green;
    data["Color:"].color[2] = annotation->m_color.blue;
    data["Color:"].isColor = true;
    dlg.order("Text:");
    dlg.order("Color:");
    if (dlg.GetResults(data)) {
      std::set<Annotation*>::iterator iter;
      for (iter = annotations.begin(); iter != annotations.end(); ++iter) {
        annotation = *iter;
        annotation->setText(data["Text:"].str.c_str());
        annotation->setColor(PaletteColor(
                               (unsigned char)data["Color:"].color[0],
                               (unsigned char)data["Color:"].color[1],
                               (unsigned char)data["Color:"].color[2]));
      }
    }
  }

  if (labels.size() > 0) {
    ATOMLABEL* label = labels.begin()->first;
    MIGenericDialog dlg(NULL, "Edit atom label");
    MIData data;
    data["Text:"].str = std::string(label->label().c_str());
    data["Color:"].color[0] = label->red();
    data["Color:"].color[1] = label->green();
    data["Color:"].color[2] = label->blue();
    data["Color:"].isColor = true;
    dlg.order("Text:");
    dlg.order("Color:");
    if (dlg.GetResults(data)) {
      std::set<std::pair<ATOMLABEL*, Molecule*> >::iterator iter;
      for (iter = labels.begin(); iter != labels.end(); ++iter) {
        label = iter->first;
        Molecule* model = iter->second;
        model->setAtomLabelText(label, data["Text:"].str.c_str());
        model->setAtomLabelColor(label,
                                 (unsigned char)data["Color:"].color[0],
                                 (unsigned char)data["Color:"].color[1],
                                 (unsigned char)data["Color:"].color[2]);
      }
    }
  }

  if (surfaces.size() > 0) {
    Molecule* model = *surfaces.begin();
    int color = MIColorChooser(model->getDots()[0].color);
    std::set<Molecule*>::iterator iter;
    for (iter = surfaces.begin(); iter != surfaces.end(); ++iter) {
      Molecule* model = *iter;
      model->setDotsColor(color);
    }
  }
}

//# (G)eometry, (V)an der Waals, (O)mega, (P)hi-psi, (C)is peptide,
//# (R)otamer chi-1, (D)ensity

std::string toErrorDescription(const std::string& s) {
  if (s.size() && s[0] == 'G') {
    return "Geometry";
  } else if (s.size() && s[0] == 'V') {
    return "Van der Waals";
  } else if (s.size() && s[0] == 'O') {
    return "Omega";
  } else if (s.size() && s[0] == 'P') {
    return "Phi-psi";
  } else if (s.size() && s[0] == 'C') {
    return "Cis peptide";
  } else if (s.size() && s[0] == 'R') {
    return "Rotamer chi-1";
  } else if (s.size() && s[0] == 'D') {
    return "Density";
  }
  return "";
}

void DisplayTree::AutoShowErrors() {
  bool state=MIConfig::Instance()->GetProfileInt("DisplayView", "autoShowError", 1) != 0;
  state = !state;
  MIConfig::Instance()->WriteProfileInt("DisplayView", "autoShowError", (int)state);
}

void DisplayTree::ImportErrors() {
  std::string file = MIFileSelector("Select error list file", "", "", "", "*.txt", 0, 0);
  if (file.size() == 0) {
    return;
  }
  FILE* errorFile = fopen(file.c_str(), "r");
  Displaylist* models = displaylist();
  if (models == NULL) {
    return;
  }
  char buffer[1024];
  while (fgets(buffer, sizeof(buffer), errorFile) != NULL) {
    if (buffer[0] == '#') {
      continue;
    }
    std::string line = buffer;
    MIStringTrim(line, true);
    MIStringTrim(line, false);

    std::vector<std::string> results;
    MIStringSplit(line," \t\r\n",results);

    if (results.size() <= 3) {
      continue;
    }
    std::string chain = results[0];
    std::string residueName = results[1];
    std::string residueType = results[2];
    std::string errors;
    for (unsigned int i=3; i< results.size(); ++i) {
      std::string error = toErrorDescription(results[i]);
      if (!error.empty()) {
        if (!errors.empty()) {
          errors += ", ";
        }
        errors += error;
      }
    }
    std::list<Molecule*>::iterator modelIter = models->begin();
    while (modelIter != models->end()) {
      Molecule* model = *modelIter;
      ++modelIter;
      if (model == NULL) {
        continue;
      }
      RESIDUE* res = residue_from_name(model->getResidues(), residueName.c_str(), chain.size() ? chain[0] : '\0');
      if (res != NULL) {
        MIAtom* atom = atom_from_name("CA", *res);
        if (atom == NULL && res->atomCount() > 0) {
          atom = res->atom(0);
        }
        if (atom != NULL) {
          std::string text;
          text=::format("Error in %s %s %s: %s", chain.c_str(), residueName.c_str(), residueType.c_str(), errors.c_str());
          Annotation* annotation = new Annotation(text.c_str(), atom->x(), atom->y(), atom->z());
          model->addAnnotation(annotation);
        }
      }
    }
  }
  fclose(errorFile);
}

Displaylist* DisplayTree::displaylist() {
  if (view == NULL) {
    return NULL;
  }
  return view->GetDisplaylist();
}

void DisplayTree::setView(MIGLWidget* view) {
  if (view == NULL) {
    return;
  }
  this->view = view;
  addModels(displaylist());
}

void DisplayTree::addModels(Displaylist* displaylist) {
  std::list<Molecule*>::iterator modelIter;
  for (modelIter = displaylist->begin(); modelIter != displaylist->end(); ++modelIter) {
    Molecule* model = *modelIter;
    addModel(model);
  }
  connect(displaylist, SIGNAL(modelAdded(Molecule*)),
          this, SLOT(modelAdded(Molecule*)));
  connect(displaylist, SIGNAL(currentMoleculeChanged(Molecule*,Molecule*)),
          this, SLOT(currentModelChanged(Molecule*,Molecule*)));
}

void DisplayTree::addModel(Molecule* model) {
  TreeData* data = new TreeData;
  data->model = model;
  QTreeWidgetItem* item = appendItem(rootId, model->pathname.c_str(), 0, 0, data);
  modelToData[model] = data;
  addAnnotations(model);
  addAtomLabels(model);
  addSurface(model);

  connect(model, SIGNAL(moleculeToBeDeleted(MIMoleculeBase*)),
          this, SLOT(modelToBeDeleted(MIMoleculeBase*)));
  connect(model, SIGNAL(annotationAdded(Molecule*,Annotation*)),
          this, SLOT(annotationAdded(Molecule*,Annotation*)));
  connect(model, SIGNAL(annotationToBeDeleted(Molecule*,Annotation*)),
          this, SLOT(annotationToBeDeleted(Molecule*,Annotation*)));
  connect(model, SIGNAL(atomLabelAdded(Molecule*,ATOMLABEL*)),
          this, SLOT(atomLabelAdded(Molecule*,ATOMLABEL*)));
  connect(model, SIGNAL(atomLabelToBeDeleted(Molecule*,AtomLabelList)),
          this, SLOT(atomLabelToBeDeleted(Molecule*,Molecule::AtomLabelList)));
  connect(model, SIGNAL(atomLabelChanged(Molecule*,ATOMLABEL*)),
          this, SLOT(atomLabelChanged(Molecule*,ATOMLABEL*)));
  connect(model, SIGNAL(surfaceChanged(Molecule*)),
          this, SLOT(surfaceChanged(Molecule*)));

  item->setExpanded(true);
}

void DisplayTree::addAnnotations(Molecule* model) {
  TreeData* modelData = modelToData[model];
  if (modelData == NULL || !validTreeData(modelData)) {
    modelToData.erase(model);
    return;
  }
  QTreeWidgetItem* modelItem = modelData->GetId();
  QTreeWidgetItem* item;
  if (modelToAnnotationListData.find(model) == modelToAnnotationListData.end()) {
    TreeData* data = new TreeData;
    data->annotationList = model;
    item = appendItem(modelItem, "Annotations", 2, 3, data);
    modelToAnnotationListData[model] = data;
  } else {
    TreeData* data = modelToAnnotationListData[model];
    if (data == NULL || !validTreeData(data)) {
      modelToAnnotationListData.erase(model);
      return;
    }
    item = data->GetId();
  }
  Molecule::AnnotationList::iterator iter;
  for (iter = model->getAnnotations().begin(); iter != model->getAnnotations().end(); ++iter) {
    Annotation* annotation = *iter;
    addAnnotation(model, annotation);
  }
}

void DisplayTree::addAnnotation(Molecule* model, Annotation* annotation) {
  TreeData* listData = modelToAnnotationListData[model];
  if (listData == NULL || !validTreeData(listData)) {
    modelToAnnotationListData.erase(model);
    return;
  }
  QTreeWidgetItem* listId = listData->GetId();
  TreeData* data = new TreeData;
  data->annotation = annotation;
  data->annotationModel = model;
  QTreeWidgetItem* item = appendItem(listId, annotation->GetText(), 2, 3, data);
  annotationToData[annotation] = data;
  connect(annotation, SIGNAL(annotationChanged(Annotation*)),
          this, SLOT(annotationChanged(Annotation*)));
  stylizeItem(item, annotation);
}

void DisplayTree::addAtomLabels(Molecule* model) {
  TreeData* modelData = modelToData[model];
  if (modelData == NULL || !validTreeData(modelData)) {
    modelToData.erase(model);
    return;
  }
  QTreeWidgetItem* modelItem = modelData->GetId();
  QTreeWidgetItem* item;
  if (modelToAtomLabelListData.find(model) == modelToAtomLabelListData.end()) {
    TreeData* data = new TreeData;
    data->atomLabelList = model;
    item = appendItem(modelItem, "Atom Labels", 6, 7, data);
    modelToAtomLabelListData[model] = data;
  } else {
    TreeData* data = modelToAtomLabelListData[model];
    if (data == NULL || !validTreeData(data)) {
      modelToAtomLabelListData.erase(model);
      return;
    }
    item = data->GetId();
  }
  Molecule::AtomLabelList::iterator iter;
  for (iter = model->getAtomLabels().begin(); iter != model->getAtomLabels().end(); ++iter) {
    ATOMLABEL* label = *iter;
    addAtomLabel(model, label);
  }
}

void DisplayTree::addAtomLabel(Molecule* model, ATOMLABEL* label) {
  TreeData* listData = modelToAtomLabelListData[model];
  if (listData == NULL || !validTreeData(listData)) {
    modelToAtomLabelListData.erase(model);
    return;
  }
  QTreeWidgetItem* listId = listData->GetId();
  TreeData* data = new TreeData;
  data->atomLabel = label;
  data->atomLabelModel = model;
  QTreeWidgetItem* item = appendItem(listId, label->label().c_str(), 6, 7, data);
  atomLabelToData[label] = data;
  stylizeItem(item, label);
}

void DisplayTree::addSurface(Molecule* model) {
  TreeData* modelData = modelToData[model];
  if (modelData == NULL || !validTreeData(modelData)) {
    modelToData.erase(model);
    return;
  }
  QTreeWidgetItem* modelItem = modelData->GetId();
  QTreeWidgetItem* item;
  if (modelToSurfaceData.find(model) == modelToSurfaceData.end()) {
    TreeData* data = new TreeData;
    data->surface = model;
    item = appendItem(modelItem, "Surface", 10, 10, data);
    modelToSurfaceData[model] = data;
  } else {
    TreeData* data = modelToSurfaceData[model];
    if (data == NULL || !validTreeData(data)) {
      modelToSurfaceData.erase(model);
      return;
    }
    item = data->GetId();
  }
  stylizeSurfaceItem(item, model);
}

void DisplayTree::stylizeSurfaceItem(QTreeWidgetItem* item, Molecule* model) {
  QBrush brush=item->foreground(0);

  if (model->getDots().size() > 0) {
    std::string text;
    text=::format("Surface (%d dots)", model->getDots().size());
    item->setText(0, text.c_str());
    brush.setColor(Qt::black);
    item->setForeground(0,brush);
    if (model->DotsVisible()) {
      item->setIcon(0, _imageList[10]);
      item->setIcon(0, _imageList[10]);
    } else {
      item->setIcon(0, _imageList[11]);
      item->setIcon(0, _imageList[11]);
    }
  } else {
    item->setText(0, "Surface");
    brush.setColor(Qt::lightGray);
    item->setForeground(0,brush);
    item->setIcon(0, _imageList[12]);
    item->setIcon(0, _imageList[12]);
  }
}

void DisplayTree::stylizeItem(QTreeWidgetItem* item, Annotation* annotation) {
  int image = 2;
  int selectedImage = 3;
  if (annotation->isHidden()) {
    image = 4;
    selectedImage = 5;
  }
  item->setIcon(0, _imageList[image]);
  item->setIcon(0, _imageList[selectedImage]);
}

void DisplayTree::stylizeItem(QTreeWidgetItem* item, ATOMLABEL* label) {
  int image = 6;
  int selectedImage = 7;
  if (!label->isVisible()) {
    image = 8;
    selectedImage = 9;
  }
  item->setIcon(0, _imageList[image]);
  item->setIcon(0, _imageList[selectedImage]);
}

void DisplayTree::annotationChanged(Annotation* annotation) {
  TreeData* data = annotationToData[annotation];
  if (data == NULL || !validTreeData(data)) {
    annotationToData.erase(annotation);
    return;
  }
  QTreeWidgetItem* item = data->GetId();
  if (item->text(0).toStdString() != annotation->m_text) {
    item->setText(0, annotation->m_text.c_str());
  }
  stylizeItem(item, annotation);
  update();
}

void DisplayTree::atomLabelChanged(Molecule*, ATOMLABEL* label) {
  TreeData* data = atomLabelToData[label];
  if (data == NULL || !validTreeData(data)) {
    atomLabelToData.erase(label);
    return;
  }
  QTreeWidgetItem* item = data->GetId();
  item->setText(0, label->label().c_str());
  stylizeItem(item, label);
  update();
}

void DisplayTree::modelAdded(Molecule* model) {
  addModel(model);
  update();
}

void DisplayTree::modelToBeDeleted(chemlib::MIMoleculeBase* mol) {
  Molecule* model = dynamic_cast<Molecule*>(mol);
  if (model == NULL) {
    return;
  }

  TreeData* data = modelToData[model];
  if (!validTreeData(data)) {
    return;
  }
  QTreeWidgetItem* item = data->GetId();
  Delete(item);
  modelToData.erase(model);
  modelToAnnotationListData.erase(model);
  modelToAtomLabelListData.erase(model);
  modelToSurfaceData.erase(model);
  update();
}

void DisplayTree::currentModelChanged(Molecule* oldModel, Molecule* newModel) {
  QTreeWidgetItem* item;
  if (oldModel != NULL) {
    TreeData* data = modelToData[oldModel];
    if (!validTreeData(data)) {
      return;
    }
    item = data->GetId();
    item->setIcon(0, _imageList[0]);
    item->setIcon(0, _imageList[0]);
  }
  if (newModel != NULL) {
    TreeData* data = modelToData[newModel];
    if (!validTreeData(data)) {
      return;
    }
    item = data->GetId();
    item->setIcon(0, _imageList[1]);
    item->setIcon(0, _imageList[1]);
  }
}

void DisplayTree::annotationAdded(Molecule* model, Annotation* annotation) {
  addAnnotation(model, annotation);
  update();
}

void DisplayTree::annotationToBeDeleted(Molecule*, Annotation* annotation) {
  if (annotationToData.find(annotation) == annotationToData.end()) {
    return;
  }
  TreeData* data = annotationToData[annotation];
  if (!validTreeData(data)) {
    return;
  }
  QTreeWidgetItem* item = data->GetId();
  Delete(item);
  annotationToData.erase(annotation);
  update();
}

void DisplayTree::atomLabelAdded(Molecule* model, ATOMLABEL* label) {
  addAtomLabel(model, label);
  update();
}

void DisplayTree::atomLabelToBeDeleted(Molecule*, Molecule::AtomLabelList labels) {
  Molecule::AtomLabelList::iterator iter = labels.begin();
  for (; iter != labels.end(); ++iter) {
    ATOMLABEL* label = *iter;
    if (atomLabelToData.find(label) == atomLabelToData.end()) {
      return;
    }
    TreeData* data = atomLabelToData[label];
    if (!validTreeData(data)) {
      return;
    }
    QTreeWidgetItem* item = data->GetId();
    Delete(item);
    atomLabelToData.erase(label);
  }
  update();
}

void DisplayTree::surfaceChanged(Molecule* model) {
  addSurface(model);
  update();
}

DisplayView::DisplayView(QWidget* parent) : ViewSyncedPanel(parent) {
}

DisplayView::~DisplayView() {
}

QWidget* DisplayView::createPanelForView(MIGLWidget* view, QWidget* parent) {
  DisplayTree* displayTree = new DisplayTree(parent);
  displayTree->setView(view);
  return displayTree;
}

void DisplayView::destroyContentsForView(MIGLWidget*, QWidget*) {
}

#include "DisplayView.moc"
