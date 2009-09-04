#ifndef mifit_ui_TreeData_h
#define mifit_ui_TreeData_h

#include <set>
#include <chemlib/chemlib.h>
#include <map/maplib.h>
#include "core/corelib.h"
#include "EMap.h"
#include "core/Annotation.h"
#include "core/ATOMLABEL.h"

class TreeData;

extern std::set<TreeData*> treeDataRegistry;

extern bool validTreeData(TreeData* data);
extern void invalidateTreeData(TreeData* data);

class QTreeWidgetItem;

class TreeData {
    QTreeWidgetItem *item;

public:
  chemlib::MIAtom* atom;
  chemlib::RESIDUE* residue;
  chemlib::RESIDUE* chain;
  CMapHeaderBase* mapHeader;
  EMap* map;
  Molecule* model;
  Molecule* annotationList;
  Molecule* atomLabelList;
  Molecule* surface;
  Annotation* annotation;
  Molecule* annotationModel;
  ATOMLABEL* atomLabel;
  Molecule* atomLabelModel;
  
    TreeData() : item(NULL), atom(NULL), residue(NULL), chain(NULL), mapHeader(NULL), map(NULL), model(NULL),
      annotationList(NULL), atomLabelList(NULL), surface(NULL),
      annotation(NULL), annotationModel(NULL),
      atomLabel(NULL), atomLabelModel(NULL) {
    treeDataRegistry.insert(this);
  }

  ~TreeData() {
    invalidateTreeData(this);
  }

    QTreeWidgetItem *GetId() const { return item; }
    void SetId(QTreeWidgetItem* id) { item=id; }

};


#endif
