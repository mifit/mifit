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
#include <QMessageBox>
#include <QInputDialog>
#include <QFileDialog>

#include <set>

#include <map/maplib.h>
#include "core/corelib.h"
#include <nongui/nonguilib.h>
#include <chemlib/Monomer.h>

#include "ModelsView.h"
#include "EMap.h"
#include "Displaylist.h"
#include "TreeData.h"
#include "Application.h"

#include "MIMenu.h"
#include "MIGLWidget.h"
#include "MIMainWindow.h"
#include "MIQTreeWidget.h"
#include "GenericDataDialog.h"
#include "ui/SelectCrystal.h"
#include "ui/MIColorPickerDlg.h"

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

using namespace chemlib;

bool syncView;
typedef std::map<const char*, int> ImageIndexMap;

class AtomsTree
    : public MIQTreeWidget
{
    Q_OBJECT

    QTreeWidgetItem *rootId;
    MIGLWidget *view;
    Molecule *model;
    Residue *residue;
    MIAtom *currentAtom;

    ImageIndexMap imageIndex;
    typedef std::map<MIAtom*, TreeData*> AtomToDataMap;
    AtomToDataMap atomToData;

    QMenu *_menu;
    bool _working;

    virtual void contextMenuEvent(QContextMenuEvent *event);

public:

    AtomsTree(QWidget *parent);
    virtual ~AtomsTree();

    void setView(MIGLWidget *view);
    void setModel(Molecule *model);
    void setResidue(Residue *residue);
    void setCurrentAtom(MIAtom *atom);
    MIAtom *getCurrentAtom();
    void stylizeAtom(MIAtom *atom);

private slots:
    void atomChanged(chemlib::MIMoleculeBase *model, chemlib::MIAtomList &atom);
    void atomsToBeDeleted(chemlib::MIMoleculeBase *model, const chemlib::MIAtomList &atoms);

    void OnItemClicked(QTreeWidgetItem *item, int column); // single click
    void OnItemActivated(QTreeWidgetItem *item, int column); // double click

    void ShowItem();
    void ColorItem();
    void DeleteItem();
    void EditItem();

};

AtomsTree::AtomsTree(QWidget *parent)
    : MIQTreeWidget(parent),
      view(NULL),
      model(NULL),
      residue(NULL)
{

    setSelectionMode(QAbstractItemView::ExtendedSelection);

    std::vector<QIcon> imageList;
    imageIndex["atom"] = imageList.size();
    QIcon atomImage = QIcon(QPixmap(atom_xpm));
    imageList.push_back(atomImage);
    imageIndex["atomSelected"] = imageList.size();
    QIcon atomSelectedImage = QIcon(QPixmap(atomSelected_xpm));
    imageList.push_back(atomSelectedImage);
    imageIndex["atomHidden"] = imageList.size();
    QIcon atomHiddenImage = QIcon(QPixmap(atomHidden_xpm));
    imageList.push_back(atomHiddenImage);
    imageIndex["atomHiddenSelected"] = imageList.size();
    QIcon atomHiddenSelectedImage = QIcon(QPixmap(atomHiddenSelected_xpm));
    imageList.push_back(atomHiddenSelectedImage);

    AssignImageList(imageList);

    std::string rootText = std::string("Atoms List");
    setHeaderLabel(rootText.c_str());
    rootId = invisibleRootItem();
    _menu = new QMenu(this);
    _menu->addAction("Show/Hide", this, SLOT(ShowItem()));
    _menu->addAction("Delete", this, SLOT(DeleteItem()));
    _menu->addAction("Color", this, SLOT(ColorItem()));
    _menu->addAction("Edit", this, SLOT(EditItem()));

    _working = false;


    connect(this, SIGNAL(itemClicked(QTreeWidgetItem *, int)),
            this, SLOT(OnItemClicked(QTreeWidgetItem *, int)));

    connect(this, SIGNAL(itemActivated(QTreeWidgetItem *, int)),
            this, SLOT(OnItemActivated(QTreeWidgetItem *, int)));

}

AtomsTree::~AtomsTree()
{
    setVisible(false);
    delete _menu;
}

void AtomsTree::setView(MIGLWidget *view)
{
    this->view = view;
}

void AtomsTree::atomChanged(MIMoleculeBase*, MIAtomList &atoms)
{
    MIAtom_iter iter = atoms.begin();
    while (iter != atoms.end())
    {
        MIAtom *atom = *iter;
        ++iter;
        stylizeAtom(atom);
    }
}

void AtomsTree::atomsToBeDeleted(MIMoleculeBase*, const MIAtomList &atoms)
{
    bool setCurrent = false;
    for (unsigned int i = 0; i<atoms.size(); ++i)
    {
        MIAtom *atom = atoms[i];
        if (atomToData.find(atom) != atomToData.end())
        {
            TreeData *data = atomToData[atom];
            if (data == NULL || !validTreeData(data))
            {
                atomToData.erase(atom);
                continue;
            }
            QTreeWidgetItem *item = data->GetId();
            if (item)
            {
                if (currentAtom == atom)
                {
                    setCurrent = true;
                }
                atomToData.erase(atom);
                Delete(item);
            }
        }
    }

    if (setCurrent)
        setCurrentAtom(NULL);
}

void AtomsTree::setModel(Molecule *model)
{
    if (this->model != NULL)
    {
        disconnect(this->model, SIGNAL(atomChanged(chemlib::MIMoleculeBase*, chemlib::MIAtomList&)));
        disconnect(this->model, SIGNAL(atomsToBeDeleted(chemlib::MIMoleculeBase*, chemlib::MIAtomList)));
    }
    this->model = model;
    if (model != NULL)
    {
        connect(model, SIGNAL(atomChanged(chemlib::MIMoleculeBase*, chemlib::MIAtomList&)),
                this, SLOT(atomChanged(chemlib::MIMoleculeBase*, chemlib::MIAtomList&)));
        connect(model, SIGNAL(atomsToBeDeleted(chemlib::MIMoleculeBase*, chemlib::MIAtomList)),
                this, SLOT(atomsToBeDeleted(chemlib::MIMoleculeBase*, chemlib::MIAtomList)));
    }
}

void AtomsTree::setResidue(Residue *residue)
{
    currentAtom = NULL;
    atomToData.clear();
    this->residue = residue;
    setVisible(false);
    DeleteChildren(rootId);
    if (residue != NULL)
    {
        for (int ia = 0; ia < residue->atomCount(); ++ia)
        {
            MIAtom *atom = residue->atom(ia);
            TreeData *data = new TreeData;
            data->atom = atom;
            QTreeWidgetItem *item = appendItem(rootId, MIAtom::liststring(atom).c_str(), imageIndex["atom"], imageIndex["atomSelected"], data);
            if (!item)
            {
                Logger::log("Error adding atom to atoms list");
                continue;
            }
            atomToData[atom] = data;
            stylizeAtom(atom);
        }
        if (currentAtom == NULL)
        {
            setCurrentAtom(atom_from_name("CA", *residue));
        }
        if (currentAtom == NULL && residue->atomCount() > 0)
        {
            setCurrentAtom(residue->atom(0));
        }

    }
    setVisible(true);
    update();
}

void AtomsTree::setCurrentAtom(MIAtom *atom)
{
    if (atom == currentAtom)
    {
        return;
    }
    MIAtom *oldAtom = currentAtom;
    currentAtom = atom;
    stylizeAtom(oldAtom);
    stylizeAtom(atom);
    if (atomToData.find(currentAtom) != atomToData.end())
    {
        TreeData *data = atomToData[currentAtom];
        if (!validTreeData(data))
        {
            atomToData.erase(currentAtom);
            return;
        }
        QTreeWidgetItem *item = data->GetId();
        if (item)
        {
            scrollToItem(item);
        }
    }
}

MIAtom*AtomsTree::getCurrentAtom()
{
    return currentAtom;
}

void AtomsTree::stylizeAtom(MIAtom *atom)
{
    if (atom == NULL)
    {
        return;
    }
    TreeData *data = atomToData[atom];
    if (!validTreeData(data))
    {
        atomToData.erase(atom);
        return;
    }
    QTreeWidgetItem *item = data->GetId();
    if (!item)
    {
        return;
    }
    int image = imageIndex["atom"];
    if (atom == currentAtom)
    {
        if (atom->isHidden())
        {
            image = imageIndex["atomHiddenSelected"];
        }
        else
        {
            image = imageIndex["atomSelected"];
        }
    }
    else if (atom->isHidden())
    {
        image = imageIndex["atomHidden"];
    }
    item->setText(0, MIAtom::liststring(atom).c_str());
    item->setIcon(0, _imageList[ image]);
    item->setIcon(0, _imageList[ image ]);
}

void AtomsTree::OnItemClicked(QTreeWidgetItem *item, int)
{
    if (_working)
    {
        return;
    }

    TreeData *data = (TreeData*) GetItemData(item);
    if (data == NULL || !validTreeData(data))
    {
        return;
    }

    if (data->atom != NULL)
    {
        MIAtom *atom = data->atom;
        setCurrentAtom(atom);
        if (syncView && currentAtom != NULL)
        {
            view->select(model, residue, currentAtom);
        }
    }
}

void AtomsTree::OnItemActivated(QTreeWidgetItem *item, int)
{
    TreeData *data = (TreeData*) GetItemData(item);
    if (data == NULL || !validTreeData(data))
    {
        return;
    }
    if (data->atom != NULL)
    {
        MIAtom *atom = data->atom;
        view->moveTo(atom->x(), atom->y(), atom->z());
    }
}

void AtomsTree::contextMenuEvent(QContextMenuEvent *event)
{
    QPoint pos = event->pos();
    QTreeWidgetItem *item = itemAt(pos);
    if (item && !selectedItems().contains(item))
    {
        clearSelection();
        item->setSelected(true);
    }
    if (selectedItems().size())
    {
        _menu->exec(QCursor::pos());
    }
}

void AtomsTree::ShowItem()
{
    QList<QTreeWidgetItem*> selected;
    GetSelections(selected);
    if (selected.size() == 0)
    {
        return;
    }
    MIAtomList atoms;
    for (int i = 0; i < selected.size(); ++i)
    {
        QTreeWidgetItem *item = selected[i];
        if (!item)
        {
            continue;
        }
        TreeData *data = (TreeData*) GetItemData(item);
        if (data == NULL || !validTreeData(data))
        {
            continue;
        }
        if (data->atom != NULL)
        {
            MIAtom *atom = data->atom;
            atoms.push_back(atom);
        }
    }
    if (atoms.size())
        model->toggleAtomsHidden(atoms);
}

void AtomsTree::ColorItem()
{
    QList<QTreeWidgetItem*> selected;
    GetSelections(selected);
    if (selected.size() == 0)
    {
        return;
    }
    MIAtomList atoms;
    for (int i = 0; i < selected.size(); ++i)
    {
        QTreeWidgetItem *item = selected[i];
        if (!item)
        {
            continue;
        }
        TreeData *data = (TreeData*) GetItemData(item);
        if (data == NULL || !validTreeData(data))
        {
            continue;
        }
        if (MIAtom::isValid(data->atom))
        {
            atoms.push_back(data->atom);
        }
    }
    if (atoms.size() == 0)
    {
        return;
    }
    int color = MIColorPickerDlg::getColor(this, atoms[0]->color());
    model->setAtomsColor(atoms, color);
}

void AtomsTree::DeleteItem()
{
    QList<QTreeWidgetItem*> selected;
    GetSelections(selected);
    if (selected.size() == 0)
    {
        return;
    }
    MIAtomList atoms;
    for (int i = 0; i < selected.size(); ++i)
    {
        QTreeWidgetItem *item = selected[i];
        if (!item)
        {
            continue;
        }
        TreeData *data = (TreeData*) GetItemData(item);
        if (data == NULL || !validTreeData(data))
        {
            continue;
        }
        if (MIAtom::isValid(data->atom))
        {
            atoms.push_back(data->atom);
        }
    }
    if (atoms.size() == 0)
    {
        return;
    }
    model->DeleteAtoms(atoms);
}

unsigned int CountAtoms(Molecule *model)
{
    unsigned int count = 0;
    ResidueListIterator currRes = model->residuesBegin();
    for (; currRes != model->residuesEnd(); ++currRes)
    {
        count += currRes->atomCount();
    }
    return count;
}

void AtomsTree::EditItem()
{
    QList<QTreeWidgetItem*> selected;
    GetSelections(selected);
    if (selected.size() == 0)
    {
        return;
    }

    MIAtomList atoms;
    for (int i = 0; i < selected.size(); ++i)
    {
        QTreeWidgetItem *item = selected[i];
        if (!item)
        {
            continue;
        }
        TreeData *data = (TreeData*) GetItemData(item);
        if (data == NULL || !validTreeData(data))
        {
            continue;
        }
        if (MIAtom::isValid(data->atom))
        {
            atoms.push_back(data->atom);
        }
    }
    if (atoms.size() == 0)
    {
        return;
    }
    std::string s;
    s = ::format("%0.4f %0.4f", atoms[0]->BValue(), atoms[0]->occ());
    float bvalue;
    float occ;
    QString str = QInputDialog::getText(this, "Edit atom", "Edit B-Value and occupancy", QLineEdit::Normal, s.c_str());
    if (!str.isEmpty())
    {
        sscanf(str.toAscii().constData(), "%f%f", &bvalue, &occ);
    }
    else
    {
        return;
    }
    model->setAtomsBValueAndOccupancy(atoms, bvalue, occ);
}

class ResiduesTree
    : public MIQTreeWidget
{
    Q_OBJECT

    QTreeWidgetItem *rootId;
    AtomsTree *atomsTree;
    MIGLWidget *view;
    Molecule *model;
    Residue *currentResidue;

    ImageIndexMap imageIndex;
    typedef std::map<Residue*, TreeData*> ResidueToDataMap;
    ResidueToDataMap residueToData;

    QMenu *_menu;
    QAction *_insertAction;
    QAction *_pasteAction;
    bool _working;

    virtual void contextMenuEvent(QContextMenuEvent *event);

public:

    ResiduesTree(QWidget *parent);
    virtual ~ResiduesTree();

    void setView(MIGLWidget *view);
    void setAtomsTree(AtomsTree *tree);
    void setModel(Molecule *model);
    Molecule *getModel();
    void setChain(Residue *chain);
    void setCurrentResidue(Residue *residue);
    Residue *getCurrentResidue();
    void stylizeResidue(Residue *residue);

private slots:
    void residuesToBeDeleted(chemlib::MIMoleculeBase *model, std::vector<chemlib::Residue*> &residues);
    void atomChanged(chemlib::MIMoleculeBase *model, chemlib::MIAtomList &atoms);

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
};


ResiduesTree::ResiduesTree(QWidget *parent)
    : MIQTreeWidget(parent),
      atomsTree(NULL),
      view(NULL),
      model(NULL),
      currentResidue(NULL)
{

    setSelectionMode(QAbstractItemView::ExtendedSelection);

    std::vector<QIcon> imageList;
    imageIndex["residue"] = imageList.size();
    QIcon residueImage = QIcon(QPixmap(residue_xpm));
    imageList.push_back(residueImage);
    imageIndex["residueSelected"] = imageList.size();
    QIcon residueSelectedImage = QIcon(QPixmap(residueSelected_xpm));
    imageList.push_back(residueSelectedImage);
    imageIndex["residuePartial"] = imageList.size();
    QIcon residuePartialImage = QIcon(QPixmap(residuePartial_xpm));
    imageList.push_back(residuePartialImage);
    imageIndex["residuePartialSelected"] = imageList.size();
    QIcon residuePartialSelectedImage = QIcon(QPixmap(residuePartialSelected_xpm));
    imageList.push_back(residuePartialSelectedImage);
    imageIndex["residueHidden"] = imageList.size();
    QIcon residueHiddenImage = QIcon(QPixmap(residueHidden_xpm));
    imageList.push_back(residueHiddenImage);
    imageIndex["residueHiddenSelected"] = imageList.size();
    QIcon residueHiddenSelectedImage = QIcon(QPixmap(residueHiddenSelected_xpm));
    imageList.push_back(residueHiddenSelectedImage);

    AssignImageList(imageList);

    std::string rootText = std::string("Residues List");
    setHeaderLabel(rootText.c_str());
    rootId = invisibleRootItem();

    _menu = new QMenu(this);
    _menu->addAction("Show/Hide", this, SLOT(ShowItem()));
    _menu->addAction("Delete", this, SLOT(DeleteItem()));
    _menu->addAction("Color", this, SLOT(ColorItem()));
    _menu->addAction("Edit", this, SLOT(EditItem()));
    _menu->addAction("Edit atoms", this, SLOT(EditItemAtoms()));
    _menu->addAction("Copy", this, SLOT(CopyItem()));
    _insertAction = _menu->addAction("Insert", this, SLOT(InsertItem()));
    _pasteAction = _menu->addAction("Paste", this, SLOT(PasteItem()));

    _working = false;


    connect(this, SIGNAL(itemClicked(QTreeWidgetItem *, int)),
            this, SLOT(OnItemClicked(QTreeWidgetItem *, int)));

    connect(this, SIGNAL(itemActivated(QTreeWidgetItem *, int)),
            this, SLOT(OnItemActivated(QTreeWidgetItem *, int)));

}

ResiduesTree::~ResiduesTree()
{
    setVisible(false);
    delete _menu;
}

void ResiduesTree::setView(MIGLWidget *view)
{
    this->view = view;
}

void ResiduesTree::setAtomsTree(AtomsTree *tree)
{
    atomsTree = tree;
}

void ResiduesTree::setModel(Molecule *model)
{
    if (this->model != NULL)
    {
        disconnect(this->model, SIGNAL(residuesToBeDeleted(chemlib::MIMoleculeBase*, std::vector<chemlib::Residue*>&)));
        disconnect(this->model, SIGNAL(atomChanged(chemlib::MIMoleculeBase*, chemlib::MIAtomList&)));
    }
    this->model = model;
    if (model != NULL)
    {
        connect(model, SIGNAL(residuesToBeDeleted(chemlib::MIMoleculeBase*, std::vector<chemlib::Residue*>&)),
                this, SLOT(residuesToBeDeleted(chemlib::MIMoleculeBase*, std::vector<chemlib::Residue*>&)));
        connect(model, SIGNAL(atomChanged(chemlib::MIMoleculeBase*, chemlib::MIAtomList&)),
                this, SLOT(atomChanged(chemlib::MIMoleculeBase*, chemlib::MIAtomList&)));
    }
}

Molecule*ResiduesTree::getModel()
{
    return model;
}

void ResiduesTree::residuesToBeDeleted(MIMoleculeBase*, std::vector<Residue*> &residues)
{
    bool setCurrent = false;
    std::vector<Residue*>::iterator iter;
    for (iter = residues.begin(); iter != residues.end(); ++iter)
    {
        Residue *residue = *iter;
        if (residueToData.find(residue) != residueToData.end())
        {
            TreeData *data = residueToData[residue];
            if (data == NULL || !validTreeData(data))
            {
                residueToData.erase(residue);
                continue;
            }
            QTreeWidgetItem *item = data->GetId();
            if (item)
            {
                if (currentResidue == residue)
                {
                    setCurrent = true;
                }
                residueToData.erase(residue);
                _working = true;
                Delete(item);
                _working = false;
            }
        }
    }
    if (setCurrent)
        setCurrentResidue(NULL);
}


void ResiduesTree::setChain(Residue *chain)
{
    setVisible(false);
    setCurrentResidue(NULL);
    ResidueToDataMap::iterator iter = residueToData.begin();
    for (iter = residueToData.begin(); iter != residueToData.end(); ++iter)
    {
        invalidateTreeData(iter->second);
    }
    residueToData.clear();
    DeleteChildren(rootId);
    if (chain != NULL)
    {
        int chain_id = chain->chain_id();
        Residue *res = chain;
        while ((res != NULL) && res->chain_id() == chain_id)
        {
            TreeData *data = new TreeData;
            data->residue = res;
            QTreeWidgetItem *item = appendItem(rootId, Residue::liststring(res).c_str(), imageIndex["residue"], imageIndex["residueSelected"], data);
            if (!item)
            {
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

void ResiduesTree::setCurrentResidue(Residue *residue)
{
    if (residue == currentResidue)
    {
        return;
    }
    Residue *oldResidue = currentResidue;
    currentResidue = residue;
    stylizeResidue(oldResidue);
    stylizeResidue(residue);
    atomsTree->setModel(model);
    atomsTree->setResidue(residue);
    if (residueToData.find(currentResidue) != residueToData.end())
    {
        TreeData *data = residueToData[currentResidue];
        if (!validTreeData(data))
        {
            residueToData.erase(currentResidue);
            return;
        }
        QTreeWidgetItem *item = data->GetId();
        if (item)
        {
            scrollToItem(item);
        }
    }
}

Residue*ResiduesTree::getCurrentResidue()
{
    return currentResidue;
}

void ResiduesTree::atomChanged(MIMoleculeBase *model, MIAtomList &atoms)
{
    std::set<Residue*> residues;
    MIAtom_iter iter = atoms.begin();
    while (iter != atoms.end())
    {
        MIAtom *atom = *iter;
        ++iter;
        Residue *res = residue_from_atom(model->residuesBegin(), atom);
        residues.insert(res);
    }
    std::set<Residue*>::iterator iter2 = residues.begin();
    while (iter2 != residues.end())
    {
        Residue *res = *iter2;
        ++iter2;
        stylizeResidue(res);
    }
}

void ResiduesTree::stylizeResidue(Residue *residue)
{
    if (residue == NULL)
    {
        return;
    }
    if (residueToData.find(residue) == residueToData.end())
    {
        return;
    }
    TreeData *data = residueToData[residue];
    if (!validTreeData(data))
    {
        residueToData.erase(residue);
        return;
    }
    QTreeWidgetItem *item = data->GetId();
    if (!item)
    {
        return;
    }

    bool partial = false;
    int shownAtoms = 0;
    int hiddenAtoms = 0;
    for (int ia = 0; ia < residue->atomCount(); ++ia)
    {
        MIAtom *atom = residue->atom(ia);
        if (atom->isHidden())
        {
            ++hiddenAtoms;
            if (shownAtoms > 0)
            {
                partial = true;
                break;
            }
        }
        else
        {
            ++shownAtoms;
            if (hiddenAtoms > 0)
            {
                partial = true;
                break;
            }
        }
    }
    bool hidden = hiddenAtoms == residue->atomCount();

    int image = 0;
    if (residue == currentResidue)
    {
        if (hidden)
        {
            image = imageIndex["residueHiddenSelected"];
        }
        else if (partial)
        {
            image = imageIndex["residuePartialSelected"];
        }
        else
        {
            image = imageIndex["residueSelected"];
        }
    }
    else
    {
        if (hidden)
        {
            image = imageIndex["residueHidden"];
        }
        else if (partial)
        {
            image = imageIndex["residuePartial"];
        }
        else
        {
            image = imageIndex["residue"];
        }
    }
    item->setIcon(0, _imageList[ image]);
    item->setIcon(0, _imageList[ image ]);
    item->setText(0, Residue::liststring(residue).c_str());
}

void ResiduesTree::OnItemClicked(QTreeWidgetItem *item, int)
{
    if (_working)
    {
        return;
    }

    TreeData *data = (TreeData*) GetItemData(item);
    if (data == NULL || !validTreeData(data))
    {
        return;
    }

    if (data->residue != NULL)
    {
        Residue *residue = data->residue;
        setCurrentResidue(residue);
        if (syncView && currentResidue != NULL)
        {
            view->select(model, currentResidue, atomsTree->getCurrentAtom());
        }
    }
}

void ResiduesTree::OnItemActivated(QTreeWidgetItem *item, int)
{
    if (view != NULL)
    {
        TreeData *data = (TreeData*) GetItemData(item);
        if (data == NULL || !validTreeData(data))
        {
            return;
        }
        if (data->residue != NULL)
        {
            Residue *residue = data->residue;
            view->setFocusResidue(residue);
        }
    }
}

void ResiduesTree::contextMenuEvent(QContextMenuEvent *event)
{
    if (Application::instance()->GetResidueBuffer() == NULL)
    {
        _insertAction->setEnabled(false);
        _pasteAction->setEnabled(false);
    }
    else
    {
        _insertAction->setEnabled(true);
        _pasteAction->setEnabled(true);
    }

    QPoint pos = event->pos();
    QTreeWidgetItem *item = itemAt(pos);
    if (item && !selectedItems().contains(item))
    {
        clearSelection();
        item->setSelected(true);
    }
    if (selectedItems().size())
    {
        _menu->exec(QCursor::pos());
    }
}

void ResiduesTree::ShowItem()
{
    QList<QTreeWidgetItem*> selected;
    GetSelections(selected);
    if (selected.size() == 0)
    {
        return;
    }
    std::vector<Residue*> residues;
    for (int i = 0; i < selected.size(); ++i)
    {
        QTreeWidgetItem *item = selected[i];
        if (!item)
        {
            continue;
        }
        TreeData *data = (TreeData*) GetItemData(item);
        if (data == NULL || !validTreeData(data))
        {
            continue;
        }
        if (data->residue != NULL)
        {
            Residue *residue = data->residue;
            residues.push_back(residue);
        }
    }
    if (residues.size())
        model->toggleResiduesHidden(residues);
}

void ResiduesTree::ColorItem()
{
    QList<QTreeWidgetItem*> selected;
    GetSelections(selected);
    if (selected.size() == 0)
    {
        return;
    }
    std::vector<Residue*> residues;
    for (int i = 0; i < selected.size(); ++i)
    {
        QTreeWidgetItem *item = selected[i];
        if (!item)
        {
            continue;
        }
        TreeData *data = (TreeData*) GetItemData(item);
        if (data == NULL || !validTreeData(data))
        {
            continue;
        }
        if (data->residue != NULL)
        {
            Residue *residue = data->residue;
            residues.push_back(residue);
        }
    }


    GenericDataDialog dlg(this);
    dlg.setWindowTitle("Choose Color");

    dlg.addColorIndexField("Color:", view->WhenShownColor);
    QStringList methods;
    methods << "Carbon Only" << "All Atoms" << "Secondary Structure" << "B-Value"
            << "Atom Type" << "Hydrophobicity" << "Shapley";
    dlg.addComboField("Method:", methods, view->WhenShownColorMethod);

    if (dlg.exec() != QDialog::Accepted)
    {
        return;
    }
    int color = dlg.value(0).toInt();
    int colorMethod = dlg.value(1).toInt();

    model->setResiduesColor(residues, color, colorMethod);
}

void ResiduesTree::DeleteItem()
{
    QList<QTreeWidgetItem*> selected;
    GetSelections(selected);
    if (selected.size() == 0)
    {
        return;
    }
    std::vector<Residue*> residues;
    for (int i = 0; i < selected.size(); ++i)
    {
        QTreeWidgetItem *item = selected[i];
        if (!item)
        {
            continue;
        }
        TreeData *data = (TreeData*) GetItemData(item);
        if (data == NULL || !validTreeData(data))
        {
            continue;
        }
        if (data->residue != NULL)
        {
            Residue *residue = data->residue;
            residues.push_back(residue);
        }
    }
    if (residues.size() == 1)
    {
        std::string mess;
        Residue *residue = residues[0];
        mess = ::format("Are you sure you want to delete residue %s?", residue->name().c_str());
        if (QMessageBox::question(this, "Confirm Delete Residue", mess.c_str(), QMessageBox::Yes | QMessageBox::No) == QMessageBox::Yes)
        {
            //REDUNDANT: view->Purge(residue);
            model->DeleteRes(residue);
        }
    }
    else if (residues.size() > 1)
    {
        std::string mess;
        mess = ::format("Are you sure you want to delete the %d residues selected?", residues.size());
        if (QMessageBox::question(this, "Confirm Delete Residues", mess.c_str(), QMessageBox::Yes | QMessageBox::No) == QMessageBox::Yes)
        {
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

void ResiduesTree::EditItem()
{
    QList<QTreeWidgetItem*> selected;
    GetSelections(selected);
    if (selected.size() == 0)
    {
        return;
    }

    std::string newname("");
    bool first = true;
    std::vector<Residue*> residues;

    for (int i = 0; i < selected.size(); ++i)
    {
        QTreeWidgetItem *item = selected[i];
        if (!item)
        {
            continue;
        }
        TreeData *data = (TreeData*) GetItemData(item);
        if (data == NULL || !validTreeData(data))
        {
            continue;
        }
        if (data->residue != NULL)
        {
            residues.push_back(data->residue);
            if (first)
            {
                first = false;
                Residue *residue = data->residue;
                std::string s;
                s = ::format("Enter new starting residue number (%d digits max)", (int)(MAXNAME -1));
                QString str = QInputDialog::getText(this, "Edit residue number", s.c_str(), QLineEdit::Normal, residue->name().c_str());
                if (str.isEmpty())
                {
                    return;
                }
                newname = str.toStdString();
            }
        }
    }

    if (residues.size())
    {
        model->setResidueNames(residues, newname.c_str());
    }
}

void ResiduesTree::EditItemAtoms()
{
    QList<QTreeWidgetItem*> selected;
    GetSelections(selected);
    if (selected.size() == 0)
    {
        return;
    }
    // Collect atoms from all residues
    MIAtomList atoms;
    for (int i = 0; i < selected.size(); ++i)
    {
        QTreeWidgetItem *item = selected[i];
        if (!item)
        {
            continue;
        }
        TreeData *data = (TreeData*) GetItemData(item);
        if (data == NULL || !validTreeData(data))
        {
            continue;
        }
        Residue *res = data->residue;
        if (Monomer::isValid(res) && res->atomCount() > 0)
        {
            for (int j = 0; j < res->atomCount(); ++j)
            {
                atoms.push_back(res->atom(j));
            }
        }
    }
    if (atoms.size() == 0)
    {
        return;
    }
    // Prompt for new b-value and occupancy
    MIAtom *atom = atoms[0];
    std::string s;
    s = ::format("%0.4f %0.4f", atom->BValue(), atom->occ());
    float bvalue;
    float occ;
    QString str = QInputDialog::getText(this, "Edit atom", "Edit B-Value and occupancy", QLineEdit::Normal, s.c_str());
    if (!str.isEmpty())
    {
        sscanf(str.toAscii().constData(), "%f%f", &bvalue, &occ);
    }
    else
    {
        return;
    }
    // Set new values on all atoms
    model->setAtomsBValueAndOccupancy(atoms, bvalue, occ);
}

void ResiduesTree::CopyItem()
{
    QList<QTreeWidgetItem*> selected;
    GetSelections(selected);
    if (selected.size() == 0)
    {
        return;
    }
    Application::instance()->ClearResidueBuffer();
    for (int i = 0; i < selected.size(); ++i)
    {
        QTreeWidgetItem *item = selected[i];
        if (!item)
        {
            continue;
        }
        TreeData *data = (TreeData*) GetItemData(item);
        if (data == NULL || !validTreeData(data))
        {
            continue;
        }
        if (data->residue != NULL)
        {
            Residue *residue = data->residue;
            Application::instance()->CopyResidueBuffer(residue);
        }
    }
}

void ResiduesTree::InsertItem()
{
    QList<QTreeWidgetItem*> selected;
    GetSelections(selected);
    if (selected.size() == 0)
    {
        return;
    }
    for (int i = 0; i < selected.size(); ++i)
    {
        QTreeWidgetItem *item = selected[i];
        if (!item)
        {
            continue;
        }
        TreeData *data = (TreeData*) GetItemData(item);
        if (data == NULL || !validTreeData(data))
        {
            continue;
        }
        if (data->residue != NULL)
        {
            Residue *residue = data->residue;
            unsigned short chain_id = residue->chain_id();
            Residue *buffer = Application::instance()->GetResidueBuffer();
            model->InsertResidues(residue, buffer, 0, chain_id);
            model->SortChains();
            break;
        }
    }

    model->Build();

}

void ResiduesTree::PasteItem()
{
    QList<QTreeWidgetItem*> selected;
    GetSelections(selected);
    if (selected.size() == 0)
    {
        return;
    }
    std::vector<Residue*> toDelete;
    for (int i = 0; i < selected.size(); ++i)
    {
        QTreeWidgetItem *item = selected[i];
        if (!item)
        {
            continue;
        }
        TreeData *data = (TreeData*) GetItemData(item);
        if (data == NULL || !validTreeData(data))
        {
            continue;
        }
        if (data->residue != NULL)
        {
            Residue *residue = data->residue;
            toDelete.push_back(residue);
        }
    }
    std::vector<Residue*> insertedResidues;
    for (int i = 0; i < selected.size(); ++i)
    {
        QTreeWidgetItem *item = selected[i];
        if (!item)
        {
            continue;
        }
        TreeData *data = (TreeData*) GetItemData(item);
        if (data == NULL || !validTreeData(data))
        {
            continue;
        }
        if (data->residue != NULL)
        {
            Residue *residue = data->residue;
            unsigned short chain_id = residue->chain_id();
            Residue *buffer = Application::instance()->GetResidueBuffer();
            model->InsertResidues(residue, buffer, 0, chain_id);
            break;
        }
    }

    model->DeleteResidues(toDelete);

    model->Build();
    std::vector<Residue*>::iterator iter = insertedResidues.begin();
    while (iter != insertedResidues.end())
    {
        Residue *res = *iter;
        if (iter == insertedResidues.begin())
        {
            setCurrentResidue(res);
        }
        ++iter;
        TreeData *data = residueToData[res];
        if (data == NULL || !validTreeData(data))
        {
            residueToData.erase(res);
            continue;
        }
        QTreeWidgetItem *item = data->GetId();
        if (!item)
        {
            continue;
        }
        item->setSelected(true);
    }
    setFocus(Qt::MouseFocusReason);
}


class ModelsTree
    : public MIQTreeWidget
{
    Q_OBJECT

    MIGLWidget *view;
    ResiduesTree *residuesTree;
    AtomsTree *atomsTree;
    Molecule *currentModel;
    Residue *currentChain;

    bool colorDialog(int &color, int &colorMethod);
    void refreshChains(Molecule *model);
    Residue *findChain(Residue *residue);
    Molecule *findModel(Residue *residue);

    int truncateWidth;

    QMenu *_menu;

    void ResetMenu(bool show_all = true);

    bool _working;

    virtual void contextMenuEvent(QContextMenuEvent *event);

public:

    ModelsTree(QWidget *parent);
    virtual ~ModelsTree();

    void setView(MIGLWidget *view);
    void setResiduesTree(ResiduesTree *tree);
    void setAtomsTree(AtomsTree *tree);
    Displaylist *displaylist();
    void addModels(Displaylist *displaylist);
    void addModel(Molecule *model);
    Molecule *getCurrentModel();
    void addModelCrystal(Molecule *model);
    void addMaps(Displaylist *displaylist);
    void addMap(EMap *map);
    void addMapCrystal(EMap *map);
    void addChains(Molecule *model);
    void setCurrentChain(Residue *residue);
    Residue *getCurrentChain();
    void stylizeModel(Molecule *model);
    void stylizeMap(EMap *map);
    void stylizeChain(Residue *residue);
    void stylizeCrystal(CMapHeaderBase *mapHeader);
    void restylize();
    void setTruncateWidth(int width);
    std::string truncateLeft(const std::string &name, int shorter = 0);
    std::string truncateRight(const std::string &name, int shorter = 0);
    void goToResidue(const std::string &name);

private slots:
    void modelAdded(Molecule *model);
    void moleculeChanged(chemlib::MIMoleculeBase *model);
    void atomChanged(chemlib::MIMoleculeBase *model, chemlib::MIAtomList &atoms);
    void modelToBeDeleted(chemlib::MIMoleculeBase *model);
    void currentModelChanged(Molecule *oldModel, Molecule *newModel);
    void mapAdded(EMap *model);
    void mapToBeDeleted(EMap *model);
    void currentMapChanged(EMap *oldModel, EMap *newModel);
    void mapFftRecalculated(EMap *map);
    void residuesToBeDeleted(chemlib::MIMoleculeBase *model, std::vector<chemlib::Residue*> &residues);
    void residuesDeleted(chemlib::MIMoleculeBase *model);
    void mapHeaderChanged(CMapHeaderBase *mapHeader);
    void focusResidueChanged(chemlib::Residue *residue);
    void selectionChanged(Molecule *model, chemlib::Residue *residue, chemlib::MIAtom *atom);

private:
    QTreeWidgetItem *rootId;
    QTreeWidgetItem *FindLastModelItem();

    ImageIndexMap imageIndex;

    typedef std::map<Molecule*, TreeData*> ModelToDataMap;
    ModelToDataMap modelToData;
    ModelToDataMap modelToCrystalData;

    typedef std::map<EMap*, TreeData*> MapToDataMap;
    MapToDataMap mapToData;
    MapToDataMap mapToCrystalData;

    typedef std::map<Residue*, Molecule*> ResidueToModelMap;
    ResidueToModelMap chainToModel;
    typedef std::map<Residue*, TreeData*> ResidueToDataMap;
    ResidueToDataMap residueToChainData;

    typedef std::map<CMapHeaderBase*, TreeData*> MapHeaderToDataMap;
    MapHeaderToDataMap mapHeaderToCrystalData;

    QAction *showAction_;
    QAction *deleteAction_;
    QAction *colorAction_;
    QAction *editAction_;
    QAction *sortChainsAction_;
    QAction *copyAction_;
    QAction *insertAction_;
    QAction *propertiesAction_;
    QAction *contourAction_;
    QAction *contourOptionsAction_;
    QAction *fftPhasesAction_;
    QAction *sfcalcAction_;
    QAction *reindexAction_;
    QAction *addFreeRAction_;
    QAction *findLigandDensityAction_;
    QAction *mapExportAction_;
    QAction *modelExportAction_;
    QAction *saveCrystalAction_;


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

    void contextActionTriggered(QAction *action);

};


ModelsTree::ModelsTree(QWidget *parent)
    : MIQTreeWidget(parent),
      view(NULL),
      currentModel(NULL),
      currentChain(NULL)
{

    setSelectionMode(QAbstractItemView::ExtendedSelection);


    std::vector<QIcon> imageList;
    QIcon modelImage = QIcon(QPixmap(model_xpm));
    imageList.push_back(modelImage);
    QIcon modelSelectedImage = QIcon(QPixmap(modelSelected_xpm));
    imageList.push_back(modelSelectedImage);
    QIcon mapImage = QIcon(QPixmap(map_xpm));
    imageList.push_back(mapImage);
    QIcon mapSelectedImage = QIcon(QPixmap(mapSelected_xpm));
    imageList.push_back(mapSelectedImage);
    QIcon crystalImage = QIcon(QPixmap(crystal_xpm));
    imageList.push_back(crystalImage);
    QIcon chainImage = QIcon(QPixmap(chain_xpm));
    imageList.push_back(chainImage);
    QIcon chainSelectedImage = QIcon(QPixmap(chainSelected_xpm));
    imageList.push_back(chainSelectedImage);
    QIcon chainPartialImage = QIcon(QPixmap(chainPartial_xpm));
    imageList.push_back(chainPartialImage);
    QIcon chainPartialSelectedImage = QIcon(QPixmap(chainPartialSelected_xpm));
    imageList.push_back(chainPartialSelectedImage);
    QIcon chainHiddenImage = QIcon(QPixmap(chainHidden_xpm));
    imageList.push_back(chainHiddenImage);
    QIcon chainHiddenSelectedImage = QIcon(QPixmap(chainHiddenSelected_xpm));
    imageList.push_back(chainHiddenSelectedImage);

    AssignImageList(imageList);

    std::string rootText = std::string("Models List");
    setHeaderLabel(rootText.c_str());
    rootId = invisibleRootItem();

    setTruncateWidth(width());

    _menu = new QMenu(this);
    showAction_ = _menu->addAction("Show/Hide");
    connect(showAction_, SIGNAL(triggered()), SLOT(ShowItem()));
    deleteAction_ = _menu->addAction("Delete");
    connect(deleteAction_, SIGNAL(triggered()), SLOT(DeleteItem()));
    colorAction_ = _menu->addAction("Color");
    connect(colorAction_, SIGNAL(triggered()), SLOT(ColorItem()));
    editAction_ = _menu->addAction("Edit");
    connect(editAction_, SIGNAL(triggered()), SLOT(EditItem()));
    sortChainsAction_ = _menu->addAction("Sort Chains");
    connect(sortChainsAction_, SIGNAL(triggered()), SLOT(SortChains()));
    copyAction_ = _menu->addAction("Copy");
    connect(copyAction_, SIGNAL(triggered()), SLOT(CopyItem()));
    insertAction_ = _menu->addAction("Insert");
    connect(insertAction_, SIGNAL(triggered()), SLOT(InsertItem()));
    propertiesAction_ = _menu->addAction("Properties");
    connect(propertiesAction_, SIGNAL(triggered()), SLOT(MapProperties()));
    contourAction_ = _menu->addAction("Contour");
    contourOptionsAction_ = _menu->addAction("Contour options...");
    fftPhasesAction_ = _menu->addAction("FFT Phases...");
    sfcalcAction_ = _menu->addAction("Calculate Structure Factors...");
    reindexAction_ = _menu->addAction("Reindex Reflections");
    addFreeRAction_ = _menu->addAction("Add Free R-flag");
    findLigandDensityAction_ = _menu->addAction("Find Ligand Density");

    mapExportAction_ = _menu->addAction("Export");
    connect(mapExportAction_, SIGNAL(triggered()), SLOT(MapExport()));
    modelExportAction_ = _menu->addAction("Export PDB...");
    connect(modelExportAction_, SIGNAL(triggered()), SLOT(ModelExport()));
    saveCrystalAction_ = _menu->addAction("Save to crystals");
    connect(saveCrystalAction_, SIGNAL(triggered()), SLOT(OnSaveToCrystals()));

    connect(_menu, SIGNAL(triggered(QAction*)), SLOT(contextActionTriggered(QAction*)));

    _working = false;

    connect(this, SIGNAL(itemClicked(QTreeWidgetItem *, int)),
            this, SLOT(OnItemClicked(QTreeWidgetItem *, int)));

    connect(this, SIGNAL(itemActivated(QTreeWidgetItem *, int)),
            this, SLOT(OnItemActivated(QTreeWidgetItem *, int)));

}

ModelsTree::~ModelsTree()
{
    delete _menu;
}

QTreeWidgetItem*ModelsTree::FindLastModelItem()
{
    QTreeWidgetItem *lastModel = 0;
    QTreeWidgetItemIterator it(this);
    for (; *it; ++it)
    {
        TreeData *data = (TreeData*) GetItemData(*it);
        if (data == NULL || !validTreeData(data))
        {
            continue;
        }
        if (data->model != NULL)
        {
            lastModel = *it;
        }
    }

    /* Not found */
    return lastModel;
}


void ModelsTree::contextMenuEvent(QContextMenuEvent *event)
{
    ResetMenu(false);

    QPoint pos = event->pos();
    QTreeWidgetItem *item = itemAt(pos);
    if (item && !selectedItems().contains(item))
    {
        clearSelection();
        item->setSelected(true);
    }
    TreeData *data = (TreeData*) GetItemData(item);
    if (data && data->map)
    {
        displaylist()->SetCurrentMap(data->map);
    }
    if (selectedItems().size())
    {
        _menu->exec(QCursor::pos());
    }

    ResetMenu(true); // must re-enable all items to make history happy
}

void ModelsTree::OnItemClicked(QTreeWidgetItem *item, int)
{
    if (_working)
    {
        return;
    }
    TreeData *data = (TreeData*) GetItemData(item);
    if (data == NULL || !validTreeData(data))
    {
        return;
    }

    if (data->model != NULL)
    {
        Molecule *model = data->model;
        displaylist()->SetCurrent(model);
    }
    else if (data->map != NULL)
    {
        EMap *map = data->map;
        displaylist()->SetCurrentMap(map);
    }
    else if (data->chain != NULL)
    {
        setCurrentChain(data->chain);
        if (syncView)
        {
            view->select(residuesTree->getModel(), residuesTree->getCurrentResidue(), atomsTree->getCurrentAtom());
        }
    }
}

void ModelsTree::contextActionTriggered(QAction *action)
{
    if (action == contourAction_)
        view->OnMapContour();
    else if (action == contourOptionsAction_)
        view->OnMapContourLevels();
    else if (action == fftPhasesAction_)
        view->OnMapFFT();
    else if (action == sfcalcAction_)
        view->OnMapSFCalc();
    else if (action == reindexAction_)
        view->OnMapReindex();
    else if (action == addFreeRAction_)
        view->OnMapAddFree();
    else if (action == findLigandDensityAction_)
        view->OnFindLigandDensity();
}

void ModelsTree::setCurrentChain(Residue *residue)
{
    if (residue == currentChain)
    {
        return;
    }
    Residue *oldChain = currentChain;
    currentChain = residue;
    if (currentChain != NULL)
    {
        if (chainToModel.find(currentChain) == chainToModel.end())
        {
            chainToModel[currentChain] = findModel(residue);
        }
        currentModel = chainToModel[currentChain];
    }
    else
    {
        currentModel = NULL;
    }
    stylizeChain(oldChain);
    stylizeChain(currentChain);
    residuesTree->setModel(currentModel);
    residuesTree->setChain(currentChain);
}

void ModelsTree::goToResidue(const std::string &nameRef)
{
    std::string name = nameRef;
    MIStringTrim(name, true);
    MIStringTrim(name, false);

    std::string residueName = name;
    char chain_id;

    Residue *chain = NULL;
    Residue *residue = NULL;
    if (currentChain != NULL)
    {
        // Search current chain using full string
        chain = currentChain;
        chain_id = currentChain->chain_id()&255;
        residue = residue_from_name(chain, residueName.c_str(), chain_id);
    }
    if (residue == NULL && currentModel != NULL)
    {
        // Search model using full string
        chain = currentModel->residuesBegin();
        chain_id = '*';
        residue = residue_from_name(chain, residueName.c_str(), chain_id);
    }
    if (residue == NULL)
    {
        // Search model using first character as chain_id and rest as name
        char chain_id = name.size()>0 ? name[0] : ' ';
        if (chain_id == '_')
        {
            chain_id = ' ';
        }
        residueName = &name[1];
        MIStringTrim(residueName, false);
        MIStringTrim(residueName, true);
        if (residueName.empty() && residuesTree->getCurrentResidue() != NULL)
        {
            residueName = residuesTree->getCurrentResidue()->name().c_str();
        }
        if (currentModel != NULL)
        {
            chain = currentModel->residuesBegin();
            residue = residue_from_name(chain, residueName.c_str(), chain_id);
        }
    }

    if (residue != NULL)
    {
        setCurrentChain(findChain(residue));
        residuesTree->setCurrentResidue(residue);
        if (syncView)
        {
            view->select(residuesTree->getModel(), residuesTree->getCurrentResidue(), atomsTree->getCurrentAtom());
            view->setFocusResidue(residuesTree->getCurrentResidue(), true);
        }
    }
    else if (name != MIToUpper(name))
    {
        goToResidue(MIToUpper(name));
    }
}

Residue*ModelsTree::getCurrentChain()
{
    return currentChain;
}

void ModelsTree::OnItemActivated(QTreeWidgetItem *item, int)
{
    if (view == NULL)
    {
        return;
    }
    TreeData *data = (TreeData*) GetItemData(item);
    if (data == NULL || !validTreeData(data))
    {
        return;
    }

    if (data->chain != NULL)
    {
        Residue *residue = data->chain;
        view->setFocusResidue(residue);
    }
}

void ModelsTree::ResetMenu(bool show_all)
{
    foreach (QAction *a, _menu->actions())
    {
        a->setVisible(show_all);
        a->setEnabled(true);
    }

    if (show_all) return;

    QList<QTreeWidgetItem*> selected;
    GetSelections(selected);
    bool modelSelected = false;
    bool chainSelected = false;
    bool mapSelected = false;
    bool mapHasPhases = false;
    bool crystalSelected = false;
    for (int i = 0; i < selected.size(); ++i)
    {
        QTreeWidgetItem *item = selected[i];
        if (!item)
            continue;

        TreeData *data = (TreeData*) GetItemData(item);
        if (data == NULL || !validTreeData(data))
            continue;

        if (data->model != NULL)
            modelSelected = true;
        else if (data->chain != NULL)
            chainSelected = true;
        else if (data->map != NULL)
        {
            EMap *map = data->map;
            mapHasPhases = map->HasPhases();
            mapSelected = true;
        }
        else if (data->mapHeader != NULL)
            crystalSelected = true;
    }

    if (modelSelected || chainSelected || mapSelected)
    {
        showAction_->setVisible(true);
        deleteAction_->setVisible(true);
        if (!mapSelected)
            colorAction_->setVisible(true);
    }
    editAction_->setVisible(true);
    if (modelSelected || mapSelected)
        editAction_->setEnabled(false);

    if (chainSelected)
        sortChainsAction_->setVisible(true);

    if (modelSelected || chainSelected)
    {
        copyAction_->setVisible(true);
        insertAction_->setVisible(true);
        if (Application::instance()->GetResidueBuffer() == NULL)
            insertAction_->setEnabled(false);
    }
    if (mapSelected)
    {
        propertiesAction_->setVisible(true);
        contourAction_->setVisible(true);
        contourOptionsAction_->setVisible(true);
        fftPhasesAction_->setVisible(true);
        sfcalcAction_->setVisible(true);
        reindexAction_->setVisible(true);
        addFreeRAction_->setVisible(true);
        findLigandDensityAction_->setVisible(true);
        mapExportAction_->setVisible(true);
        fftPhasesAction_->setEnabled(mapHasPhases);
        mapExportAction_->setEnabled(mapHasPhases);
    }

    if (modelSelected)
        modelExportAction_->setVisible(true);

    if (crystalSelected)
        saveCrystalAction_->setVisible(true);
}



void ModelsTree::OnSaveToCrystals()
{

    QList<QTreeWidgetItem*> selected;
    GetSelections(selected);
    if (selected.size() == 0)
    {
        return;
    }
    std::set<Molecule*> chainsModified;
    for (int i = 0; i < selected.size(); ++i)
    {
        QTreeWidgetItem *item = selected[i];
        if (!item)
        {
            continue;
        }
        TreeData *data = (TreeData*) GetItemData(item);
        if (data == NULL || !validTreeData(data))
        {
            continue;
        }

        if (data->mapHeader != NULL)
        {
            QString crystal = QInputDialog::getText(this, "Crystal name", "Crystal name:");
            if (!crystal.isEmpty())
            {
                data->mapHeader->SaveCrystal(crystal.toStdString());
            }
        }
    }
}

void ModelsTree::ShowItem()
{
    QList<QTreeWidgetItem*> selected;
    GetSelections(selected);
    if (selected.size() == 0)
    {
        return;
    }
    for (int i = 0; i < selected.size(); ++i)
    {
        QTreeWidgetItem *item = selected[i];
        if (!item)
        {
            continue;
        }
        TreeData *data = (TreeData*) GetItemData(item);
        if (data == NULL || !validTreeData(data))
        {
            continue;
        }

        if (data->model != NULL)
        {
            Molecule *model = data->model;
            if (model->Visible())
            {
                model->Hide();
            }
            else
            {
                model->Show();
            }
        }
        else if (data->map != NULL)
        {
            EMap *map = data->map;
            if (map->Visible())
            {
                map->Hide();
            }
            else
            {
                map->Show();
            }
        }
        else if (data->chain != NULL)
        {
            Residue *chain = data->chain;
            Molecule *model = chainToModel[chain];
            model->toggleChainHidden(chain);
        }
    }

}

void ModelsTree::ColorItem()
{
    QList<QTreeWidgetItem*> selected;
    GetSelections(selected);
    if (selected.size() == 0)
    {
        return;
    }
    std::vector<Molecule*> models;
    std::vector<Residue*> chains;
    std::vector<EMap*> maps;
    for (int i = 0; i < selected.size(); ++i)
    {
        QTreeWidgetItem *item = selected[i];
        if (!item)
        {
            continue;
        }
        TreeData *data = (TreeData*) GetItemData(item);
        if (data == NULL || !validTreeData(data))
        {
            continue;
        }

        if (data->model != NULL)
        {
            models.push_back(data->model);
        }
        else if (data->chain != NULL)
        {
            chains.push_back(data->chain);
        }
        else if (data->map != NULL)
        {
            maps.push_back(data->map);
        }
    }
    int color;
    int colorMethod;
    if (models.size() > 0 || chains.size() > 0)
    {
        if (!colorDialog(color, colorMethod))
        {
            return;
        }
    }
    if (models.size() > 0)
    {
        std::vector<Molecule*>::iterator iter = models.begin();
        for (; iter != models.end(); ++iter)
        {
            Molecule *model = *iter;
            model->setColor(color, colorMethod);
        }
    }
    if (chains.size() > 0)
    {
        std::vector<Residue*>::iterator iter = chains.begin();
        for (; iter != chains.end(); ++iter)
        {
            Residue *chain = *iter;
            Molecule *model = chainToModel[chain];
            model->setChainColor(chain, color, colorMethod);
        }
    }
    if (maps.size() > 0)
    {
        std::vector<EMap*>::iterator iter = maps.begin();
        for (; iter != maps.end(); ++iter)
        {
            EMap *map = *iter;
            map->ContourLevels();
        }
    }
}

void ModelsTree::DeleteItem()
{
    QList<QTreeWidgetItem*> selected;
    GetSelections(selected);
    if (selected.size() == 0)
    {
        return;
    }
    std::vector<Molecule*> models;
    std::vector<Residue*> chains;
    std::vector<EMap*> maps;
    for (int i = 0; i < selected.size(); ++i)
    {
        QTreeWidgetItem *item = selected[i];
        if (!item)
        {
            continue;
        }
        TreeData *data = (TreeData*) GetItemData(item);
        if (data == NULL || !validTreeData(data))
        {
            continue;
        }

        if (data->model != NULL)
        {
            models.push_back(data->model);
        }
        else if (data->chain != NULL)
        {
            chains.push_back(data->chain);
        }
        else if (data->map != NULL)
        {
            maps.push_back(data->map);
        }
    }

    if (models.size() > 0)
    {
        std::string mess;
        if (models.size() == 1)
        {
            Molecule *model = models[0];
            mess = ::format("Are you sure you want to delete model\n%s?", model->pathname.c_str());
        }
        else
        {
            mess = ::format("Are you sure you want to delete the %d models selected?", models.size());
        }
        if (QMessageBox::question(this, "Confirm Delete Model", mess.c_str(), QMessageBox::Yes | QMessageBox::No) == QMessageBox::Yes)
        {
            std::vector<Molecule*>::iterator iter = models.begin();
            for (; iter != models.end(); ++iter)
            {
                Molecule *model = *iter;
                mess = ::format("Deleted model %s", model->pathname.c_str());
                view->SaveModelFile(model, mess.c_str());
                //view->GetDisplaylist()->DeleteItem(model);
                delete model;
            }
        }
    }
    else if (chains.size() > 0)
    {
        std::string mess;
        if (chains.size() == 1)
        {
            Residue *chain = chains[0];
            mess = ::format("Are you sure you want to delete chain\n%s?", chainstring(chain).c_str());
        }
        else
        {
            mess = ::format("Are you sure you want to delete the %d chains selected?", chains.size());
        }
        if (QMessageBox::question(this, "Confirm Delete Chains", mess.c_str(), QMessageBox::Yes | QMessageBox::No) == QMessageBox::Yes)
        {
            std::vector<Residue*>::iterator iter = chains.begin();
            for (; iter != chains.end(); ++iter)
            {
                Residue *chain = *iter;
                Molecule *model = chainToModel[chain];
                mess = ::format("Deleted chain %s", chainstring(chain).c_str());
                view->SaveModelFile(model, mess.c_str());

                if (currentChain == chain)
                {
                    currentChain = NULL;
                }
                //REDUNDANT: view->PurgeChain(chain);
                model->DeleteChain(chain);
            }
        }
    }
    else if (maps.size() > 0)
    {
        std::string mess;
        if (maps.size() == 1)
        {
            EMap *map = maps[0];
            mess = ::format("Are you sure you want to delete map\n%s?", map->mapName.c_str());
        }
        else
        {
            mess = ::format("Are you sure you want to delete the %d maps selected?", maps.size());
        }
        if (QMessageBox::question(this, "Confirm Delete Maps", mess.c_str(), QMessageBox::Yes | QMessageBox::No) == QMessageBox::Yes)
        {
            std::vector<EMap*>::iterator iter = maps.begin();
            for (; iter != maps.end(); ++iter)
            {
                EMap *map = *iter;
                view->Purge(map);
            }
        }
    }
}

void ModelsTree::EditItem()
{
    QList<QTreeWidgetItem*> selected;
    GetSelections(selected);
    if (selected.size() == 0)
    {
        return;
    }
    std::vector<Residue*> chains;
    std::vector<EMap*> maps;
    std::vector<CMapHeaderBase*> crystals;
    for (int i = 0; i < selected.size(); ++i)
    {
        QTreeWidgetItem *item = selected[i];
        if (!item)
        {
            continue;
        }
        TreeData *data = (TreeData*) GetItemData(item);
        if (data == NULL || !validTreeData(data))
        {
            continue;
        }

        if (data->chain != NULL)
        {
            chains.push_back(data->chain);
        }
        else if (data->map != NULL)
        {
            maps.push_back(data->map);
        }
        else if (data->mapHeader != NULL)
        {
            crystals.push_back(data->mapHeader);
        }
    }
    std::set<Molecule*> chainsModified;
    if (chains.size() > 0)
    {
        std::vector<Residue*>::iterator iter = chains.begin();
        for (; iter != chains.end(); ++iter)
        {
            Residue *chain = *iter;
            Molecule *model = chainToModel[chain];
            int chain_id = chain->chain_id();
            char chainid = (char)(chain_id&255);

            QString s = QInputDialog::getText(this, "Chain ID", "Enter Chain ID (single char)", QLineEdit::Normal, QString(chainid));
            if (!s.isEmpty())
            {
                unsigned char c = s.toAscii().at(0);
                model->setChainId(chain, c);
                chainsModified.insert(model);
            }
            bool ok;
            int n = QInputDialog::getInt(0, "Residue number", "Enter number for first residue (rest of chain will also be renumbered)", 1, 1, 999, 1, &ok);
            if (ok)
            {
                model->renumberChain(chain, n);
                chainsModified.insert(model);
            }
        }
    }
    if (crystals.size() > 0)
    {
        MIData data;
        data["info"].str = crystals[0]->Label();
        SelectCrystal::doSelectCrystal(data);
        std::vector<CMapHeaderBase*>::iterator crystalIter = crystals.begin();
        for (; crystalIter != crystals.end(); ++crystalIter)
        {
            CMapHeaderBase *mapHeader = *crystalIter;
            CMapHeaderBase mh(data["info"].str);
            mapHeader->updateSymmetryAndCell(mh);

            for (int i = 0; i < displaylist()->MapCount(); ++i)
            {
                EMapBase *map = displaylist()->GetMap(i);
                if (map->mapheader == mapHeader)
                {
                    map->RecalcResolution();
                    if (map->HasPhases() &&  map->FFTMap())
                        view->doMapContour(map);
                }
            }

        }
    }
    std::set<Molecule*>::iterator iter = chainsModified.begin();
    while (iter != chainsModified.end())
    {
        Molecule *model = *iter;
        ++iter;
        refreshChains(model);
    }
}

void ModelsTree::SortChains()
{
    QList<QTreeWidgetItem*> selected;
    GetSelections(selected);
    if (selected.size() == 0)
    {
        return;
    }
    for (int i = 0; i < selected.size(); ++i)
    {
        QTreeWidgetItem *item = selected[i];
        if (!item)
        {
            continue;
        }
        TreeData *data = (TreeData*) GetItemData(item);
        if (data == NULL || !validTreeData(data))
        {
            continue;
        }

        if (data->chain != NULL)
        {
            Residue *chain = data->chain;
            Molecule *model = chainToModel[chain];
            view->SaveModelFile(model, "Sorted chains");
            model->SortChains();
            refreshChains(model);
        }
    }
}

void ModelsTree::CopyItem()
{
    QList<QTreeWidgetItem*> selected;
    GetSelections(selected);
    if (selected.size() == 0)
    {
        return;
    }
    Application::instance()->ClearResidueBuffer();
    for (int i = 0; i < selected.size(); ++i)
    {
        QTreeWidgetItem *item = selected[i];
        if (!item)
        {
            continue;
        }
        TreeData *data = (TreeData*) GetItemData(item);
        if (data == NULL || !validTreeData(data))
        {
            continue;
        }
        if (data->chain != NULL)
        {
            Residue *chain = data->chain;
            Residue *res = chain;
            while (res != NULL && res->chain_id() == chain->chain_id())
            {
                Application::instance()->CopyResidueBuffer(res);
                res = res->next();
            }
        }
        else if (data->residue != NULL)
        {
            Residue *residue = data->residue;
            Application::instance()->CopyResidueBuffer(residue);
        }
    }
}

void ModelsTree::InsertItem()
{
    QList<QTreeWidgetItem*> selected;
    GetSelections(selected);
    if (selected.size() == 0)
    {
        return;
    }
    for (int i = 0; i < selected.size(); ++i)
    {
        QTreeWidgetItem *item = selected[i];
        if (!item)
        {
            continue;
        }
        TreeData *data = (TreeData*) GetItemData(item);
        if (data == NULL || !validTreeData(data))
        {
            continue;
        }
        if (data->model != NULL)
        {
            Molecule *model = data->model;
            QString chainId = QInputDialog::getText(this, "New chain ID", "Chain ID for new chain", QLineEdit::Normal, "N");
            if (chainId.isEmpty())
            {
                return;
            }
            unsigned char c = chainId.toAscii().at(0);
            Residue *buffer = Application::instance()->GetResidueBuffer();
            model->InsertResidues(model->residuesBegin(), buffer, 3, (c&255) + 1*256);
            model->Build();
            model->SortChains();
            addChains(model);
            setFocus(Qt::MouseFocusReason);
            break;
        }
        else if (data->chain != NULL)
        {
            Residue *chain = data->chain;
            Molecule *model = chainToModel[chain];
            unsigned short chain_id = chain->chain_id();
            Residue *buffer = Application::instance()->GetResidueBuffer();
            model->InsertResidues(chain, buffer, 3, chain_id);
            model->Build();
            model->SortChains();
            addChains(model);
            setFocus(Qt::MouseFocusReason);
            break;
        }
    }
}

void ModelsTree::MapProperties()
{
    QList<QTreeWidgetItem*> selected;
    GetSelections(selected);
    if (selected.size() == 0)
    {
        return;
    }
    for (int i = 0; i < selected.size(); ++i)
    {
        QTreeWidgetItem *item = selected[i];
        if (!item)
        {
            continue;
        }
        TreeData *data = (TreeData*) GetItemData(item);
        if (data == NULL || !validTreeData(data))
        {
            continue;
        }
        if (data->map != NULL)
        {
            EMap *map = data->map;
            QMessageBox::information(this, "Map Properties", map->Info().c_str());
        }
    }
}

void ModelsTree::MapExport()
{
    QList<QTreeWidgetItem*> selected;
    GetSelections(selected);
    if (selected.size() == 0)
    {
        return;
    }
    for (int i = 0; i < selected.size(); ++i)
    {
        QTreeWidgetItem *item = selected[i];
        if (!item)
        {
            continue;
        }
        TreeData *data = (TreeData*) GetItemData(item);
        if (data == NULL || !validTreeData(data))
        {
            continue;
        }
        if (data->map != NULL)
        {
            EMap *map = data->map;
            map->Export();
        }
    }
}

void ModelsTree::ModelExport()
{
    QList<QTreeWidgetItem*> selected;
    GetSelections(selected);
    if (selected.size() == 0)
    {
        return;
    }
    for (int i = 0; i < selected.size(); ++i)
    {
        QTreeWidgetItem *item = selected[i];
        if (!item)
        {
            continue;
        }
        TreeData *data = (TreeData*) GetItemData(item);
        if (data == NULL || !validTreeData(data))
        {
            continue;
        }
        if (data->model != NULL)
        {
            Molecule *model = data->model;
            QString s = QFileDialog::getSaveFileName(this, "Choose a name for the model file", "",
                                                  "PDB file (*.pdb);;All files (*.*)");
            if (!s.isEmpty())
            {
                model->SavePDBFile(s.toAscii().constData());
            }
        }
    }
}

void ModelsTree::setView(MIGLWidget *view)
{
    if (view == NULL)
    {
        return;
    }
    this->view = view;
    connect(view, SIGNAL(focusResidueChanged(chemlib::Residue*)),
            this, SLOT(focusResidueChanged(chemlib::Residue*)));
    addModels(displaylist());
    addMaps(displaylist());
}

void ModelsTree::setResiduesTree(ResiduesTree *tree)
{
    residuesTree = tree;
}

void ModelsTree::setAtomsTree(AtomsTree *tree)
{
    atomsTree = tree;
}

Residue*ModelsTree::findChain(Residue *residue)
{
    if (!residue)
        return NULL;

    Residue *chain = NULL;
    bool found = false;
    Displaylist::ModelList::iterator modelIter = displaylist()->begin();
    while (!found && modelIter != displaylist()->end())
    {
        Molecule *model = *modelIter;
        ++modelIter;
        unsigned short chain_id = 0;
        for (Residue *res = model->residuesBegin(); res != NULL; res = res->next())
        {
            if ((res->chain_id()) != chain_id)
            {
                chain = res;
                chain_id = chain->chain_id();
            }
            if (res == residue)
            {
                found = true;
                break;
            }
        }
    }
    if (!found)
    {
        return NULL;
    }
    return chain;
}

Molecule*ModelsTree::findModel(Residue *residue)
{
    Molecule *model = NULL;
    bool found = false;
    Displaylist::ModelList::iterator modelIter = displaylist()->begin();
    while (!found && modelIter != displaylist()->end())
    {
        model = *modelIter;
        ++modelIter;
        for (Residue *res = model->residuesBegin(); res != NULL; res = res->next())
        {
            if (res == residue)
            {
                found = true;
                break;
            }
        }
    }
    if (!found)
    {
        return NULL;
    }
    return model;
}

void ModelsTree::focusResidueChanged(Residue *residue)
{
    if (syncView)
    {
        Residue *chain = findChain(residue);
        if (chain != NULL)
        {
            setCurrentChain(chain);
            residuesTree->setCurrentResidue(residue);
        }
    }
}

void ModelsTree::selectionChanged(Molecule*, Residue *residue, MIAtom *atom)
{
    if (syncView)
    {
        Residue *chain = findChain(residue);
        if (chain != NULL)
        {
            setCurrentChain(chain);
            residuesTree->setCurrentResidue(residue);
            atomsTree->setCurrentAtom(atom);
        }
    }
}

Displaylist*ModelsTree::displaylist()
{
    return view->GetDisplaylist();
}

void ModelsTree::addModels(Displaylist *displaylist)
{
    Displaylist::ModelList::iterator modelIter = displaylist->begin();
    while (modelIter != displaylist->end())
    {
        Molecule *model = *modelIter;
        ++modelIter;
        addModel(model);
    }
    connect(displaylist, SIGNAL(modelAdded(Molecule*)),
            this, SLOT(modelAdded(Molecule*)));
    connect(displaylist, SIGNAL(currentMoleculeChanged(Molecule*, Molecule*)),
            this, SLOT(currentModelChanged(Molecule*, Molecule*)));
    connect(displaylist, SIGNAL(selectionChanged(Molecule*, chemlib::Residue*, chemlib::MIAtom*)),
            this, SLOT(selectionChanged(Molecule*, chemlib::Residue*, chemlib::MIAtom*)));
}

void ModelsTree::addMaps(Displaylist *displaylist)
{
    Displaylist::MapList::iterator mapIter = displaylist->getMaps().begin();
    while (mapIter != displaylist->getMaps().end())
    {
        EMap *map = *mapIter;
        ++mapIter;
        addMap(map);
    }
    connect(displaylist, SIGNAL(mapAdded(EMap*)),
            this, SLOT(mapAdded(EMap*)));
    connect(displaylist, SIGNAL(mapToBeDeleted(EMap*)),
            this, SLOT(mapToBeDeleted(EMap*)));
    connect(displaylist, SIGNAL(currentMapChanged(EMap*, EMap*)),
            this, SLOT(currentMapChanged(EMap*, EMap*)));
}

void ModelsTree::addModel(Molecule *model)
{
    QTreeWidgetItem *previousItem = FindLastModelItem();
    TreeData *data = new TreeData;
    data->model = model;
    QTreeWidgetItem *item;
    if (previousItem)
    {
        item = insertItem(rootId, previousItem, model->pathname.c_str(), 0, 0, data);
    }
    else
    {
        item = prependItem(rootId, model->pathname.c_str(), 0, 0, data);
    }
    if (!item)
    {
        return;
    }
    modelToData[model] = data;
    stylizeModel(model);

    model->SortChains();
    addModelCrystal(model);
    addChains(model);

    connect(model, SIGNAL(moleculeToBeDeleted(chemlib::MIMoleculeBase*)),
            this, SLOT(modelToBeDeleted(chemlib::MIMoleculeBase*)));
    connect(model, SIGNAL(residuesToBeDeleted(chemlib::MIMoleculeBase*, std::vector<chemlib::Residue*>&)),
            this, SLOT(residuesToBeDeleted(chemlib::MIMoleculeBase*, std::vector<chemlib::Residue*>&)));
    connect(model, SIGNAL(residuesDeleted(chemlib::MIMoleculeBase*)),
            this, SLOT(residuesDeleted(chemlib::MIMoleculeBase*)));
    connect(model, SIGNAL(moleculeChanged(chemlib::MIMoleculeBase*)),
            this, SLOT(moleculeChanged(chemlib::MIMoleculeBase*)));
    connect(model, SIGNAL(atomChanged(chemlib::MIMoleculeBase*, chemlib::MIAtomList&)),
            this, SLOT(atomChanged(chemlib::MIMoleculeBase*, chemlib::MIAtomList&)));

    expandItem(item);
}

void ModelsTree::moleculeChanged(MIMoleculeBase *model)
{
    Molecule *m = dynamic_cast<Molecule*>(model);
    if (m != NULL)
    {
        refreshChains(m);
    }
}

void ModelsTree::residuesToBeDeleted(MIMoleculeBase*, std::vector<Residue*> &residues)
{
    std::vector<Residue*>::iterator iter;
    for (iter = residues.begin(); iter != residues.end(); ++iter)
    {
        Residue *residue = *iter;
        if (residueToChainData.find(residue) != residueToChainData.end())
        {
            if (currentChain == residue)
            {
                setCurrentChain(NULL);
            }
            TreeData *data = residueToChainData[residue];
            if (data == NULL || !validTreeData(data))
            {
                residueToChainData.erase(residue);
                continue;
            }
            QTreeWidgetItem *item = data->GetId();
            if (item)
            {
                residueToChainData.erase(residue);
                chainToModel.erase(residue);
                Delete(item);
            }
        }
    }
}

void ModelsTree::residuesDeleted(MIMoleculeBase *model)
{
    Molecule *m = dynamic_cast<Molecule*>(model);
    if (m != NULL)
    {
        refreshChains(m);
    }
}

void ModelsTree::addModelCrystal(Molecule *model)
{
    TreeData *modelData = modelToData[model];
    if (modelData == NULL || !validTreeData(modelData))
    {
        modelToData.erase(model);
        return;
    }
    QTreeWidgetItem *modelItem = modelData->GetId();
    if (!modelItem)
    {
        return;
    }
    CMapHeaderBase &mapHeader = model->GetMapHeader();
    connect(&mapHeader, SIGNAL(mapHeaderChanged(CMapHeaderBase*)),
            this, SLOT(mapHeaderChanged(CMapHeaderBase*)));
    TreeData *data = new TreeData;
    data->mapHeader = &mapHeader;
    QTreeWidgetItem *item = insertItem(modelItem, 0, mapHeader.Label().c_str(), 4, 4, data);
    if (!item)
    {
        return;
    }
    modelToCrystalData[model] = data;
    mapHeaderToCrystalData[&mapHeader] = data;
    stylizeCrystal(&mapHeader);
}

Molecule*ModelsTree::getCurrentModel()
{
    return currentModel;
}

void ModelsTree::refreshChains(Molecule *model)
{
    model->SortChains();
    // Delete all chain items for this model
    std::vector<TreeData*> chainDataForModel;
    ResidueToDataMap::iterator iter = residueToChainData.begin();
    while (iter != residueToChainData.end())
    {
        Residue *residue = iter->first;
        TreeData *data = iter->second;
        ++iter;
        if (validTreeData(data))
        {
            if (chainToModel.find(residue) != chainToModel.end()
                && chainToModel[residue] == model)
            {
                chainDataForModel.push_back(data);
            }
        }
    }

    std::vector<TreeData*>::iterator iter2;
    for (iter2 = chainDataForModel.begin(); iter2 != chainDataForModel.end(); ++iter2)
    {
        TreeData *data = *iter2;
        QTreeWidgetItem *item = data->GetId();
        if (item && data->chain != NULL)
        {
            Residue *chain = data->chain;
            Delete(item);
            residueToChainData.erase(chain);
            chainToModel.erase(chain);
            if (currentChain == chain)
            {
                setCurrentChain(NULL);
            }
        }
    }
    addChains(model);
}

void ModelsTree::addChains(Molecule *model)
{
    TreeData *modelData = modelToData[model];
    if (modelData == NULL || !validTreeData(modelData))
    {
        modelToData.erase(model);
        return;
    }
    QTreeWidgetItem *modelItem = modelData->GetId();
    if (!modelItem)
    {
        return;
    }
    std::set<TreeData*> currentChainDataSet;
    int chain_id = -9999;
    QTreeWidgetItem *previousChainItem;
    for (Residue *res = model->residuesBegin(); res != NULL; res = res->next())
    {
        if (res->chain_id() != chain_id)
        {
            chain_id = res->chain_id();
            if (residueToChainData.find(res) == residueToChainData.end())
            {
                std::string text = chainstring(res).c_str();
                QTreeWidgetItem *chainItem;
                TreeData *data = new TreeData;
                data->chain = res;
                if (previousChainItem)
                {
                    chainItem = insertItem(modelItem, previousChainItem, text, 5, 5, data);
                }
                else
                {
                    chainItem = appendItem(modelItem, text, 5, 5, data);
                }
                if (chainItem)
                {
                    previousChainItem = chainItem;
                    residueToChainData[res] = data;
                    chainToModel[res] = model;
                    stylizeChain(res);
                    currentChainDataSet.insert(data);
                }
            }
            else
            {
                TreeData *data = residueToChainData[res];
                if (data == NULL || !validTreeData(data))
                {
                    residueToChainData.erase(res);
                    continue;
                }
                previousChainItem = data->GetId();
                if (previousChainItem)
                {
                    stylizeChain(res);
                    currentChainDataSet.insert(data);
                }
            }
        }
    }
    // Delete any items for chains not found
    ResidueToDataMap::iterator iter = residueToChainData.begin();
    while (iter != residueToChainData.end())
    {
        TreeData *data = iter->second;
        ++iter;
        if (data == NULL || !validTreeData(data))
        {
            continue;
        }
        QTreeWidgetItem *chainItem = data->GetId();
        if (chainItem && data->chain != NULL)
        {
            Residue *chain = data->chain;
            Molecule *chainModel = chainToModel[chain];
            if (model == chainModel
                && currentChainDataSet.find(data) == currentChainDataSet.end())
            {
                Delete(chainItem);
            }
        }
    }
    if (currentChain == NULL)
    {
        setCurrentChain(model->residuesBegin());
    }
}

void ModelsTree::restylize()
{
    ModelToDataMap::iterator iter = modelToData.begin();
    while (iter != modelToData.end())
    {
        Molecule *model = iter->first;
        TreeData *data = iter->second;
        ++iter;
        if (data == NULL || !validTreeData(data))
        {
            continue;
        }
        stylizeModel(model);
    }
    ResidueToDataMap::iterator iter2 = residueToChainData.begin();
    while (iter2 != residueToChainData.end())
    {
        Residue *chain = iter2->first;
        TreeData *data = iter2->second;
        ++iter2;
        if (data == NULL || !validTreeData(data))
        {
            continue;
        }
        stylizeChain(chain);
    }
    MapToDataMap::iterator iter3 = mapToData.begin();
    while (iter3 != mapToData.end())
    {
        EMap *map = iter3->first;
        TreeData *data = iter3->second;
        ++iter3;
        if (data == NULL || !validTreeData(data))
        {
            continue;
        }
        stylizeMap(map);
    }
    MapHeaderToDataMap::iterator iter4 = mapHeaderToCrystalData.begin();
    while (iter4 != mapHeaderToCrystalData.end())
    {
        CMapHeaderBase *mapHeader = iter4->first;
        TreeData *data = iter4->second;
        ++iter4;
        if (data == NULL || !validTreeData(data))
        {
            continue;
        }
        stylizeCrystal(mapHeader);
    }
}

void ModelsTree::setTruncateWidth(int width)
{
    truncateWidth = width - 64;
}

std::string ModelsTree::truncateLeft(const std::string &name, int shorter)
{
    if (!name.size())
        return name;

    QFontMetrics fm = fontMetrics();
    size_t width = truncateWidth - shorter;
    std::string prefix("...");
    std::string result = name;
    size_t stringWidth = fm.width(result.c_str());
    for (size_t i = 1; stringWidth > width && i < name.size(); ++i)
    {
        result = prefix + &name[i];
        stringWidth = fm.width(result.c_str());
    }

    return result;
}

std::string ModelsTree::truncateRight(const std::string &name, int shorter)
{
    if (!name.size())
        return name;

    QFontMetrics fm = fontMetrics();
    size_t width = truncateWidth - shorter;
    std::string suffix("...");
    std::string shorterName(name);
    std::string result(name);
    size_t stringWidth = fm.width(result.c_str());
    for (size_t i = 1; stringWidth > width && i < name.size(); ++i)
    {
        shorterName.resize(name.size() - i);
        result = shorterName + suffix;
        stringWidth = fm.width(result.c_str());
    }

    return result;
}

void ModelsTree::stylizeModel(Molecule *model)
{
    std::string text = truncateLeft(model->pathname.c_str());
    TreeData *data = modelToData[model];
    if (data == NULL || !validTreeData(data))
    {
        modelToData.erase(model);
        return;
    }
    QTreeWidgetItem *item = data->GetId();
    if (!item)
    {
        return;
    }
    item->setText(0, text.c_str());

}

void ModelsTree::stylizeMap(EMap *map)
{
    TreeData *data = mapToData[map];
    if (data == NULL || !validTreeData(data))
    {
        mapToData.erase(map);
        return;
    }
    QTreeWidgetItem *item = data->GetId();
    if (!item)
    {
        return;
    }

    std::string mapName = map->mapName.c_str();
    int mapType = map->mapheader->maptype;
    if (mapType == (int)MIMapType::DirectFFT)
    {
        if (map->fColumnName.size() > 0)
        {
            mapName += ": ";
            mapName += map->fColumnName.c_str();
        }
    }
    else
    {
        mapName += ": ";
        mapName += StringForMapType(mapType);
    }

    std::string text = truncateLeft(mapName.c_str());
    item->setText(0, text.c_str());
}

void ModelsTree::atomChanged(MIMoleculeBase *model, MIAtomList &atoms)
{
    std::set<Residue*> chains;
    MIAtom_iter iter = atoms.begin();
    while (iter != atoms.end())
    {
        MIAtom *atom = *iter;
        ++iter;
        Residue *res = residue_from_atom(model->residuesBegin(), atom);
        if (res != NULL)
        {
            Residue *chain = findChain(res);
            if (chain != NULL)
            {
                chains.insert(chain);
            }
        }
    }
    std::set<Residue*>::iterator iter2 = chains.begin();
    while (iter2 != chains.end())
    {
        Residue *chain = *iter2;
        ++iter2;
        stylizeChain(chain);
    }
}

void ModelsTree::stylizeChain(Residue *residue)
{
    if (residue == NULL || residueToChainData.find(residue) == residueToChainData.end())
    {
        return;
    }
    TreeData *data = residueToChainData[residue];
    if (data == NULL || !validTreeData(data))
    {
        residueToChainData.erase(residue);
        return;
    }
    QTreeWidgetItem *item = data->GetId();
    if (!item)
    {
        return;
    }

    bool partial = false;
    bool hidden = true;
    int shownAtoms = 0;
    int hiddenAtoms = 0;
    Residue *res = residue;
    int chainId = res->chain_id();
    while (res != NULL && chainId == res->chain_id())
    {
        for (int ia = 0; ia < res->atomCount(); ++ia)
        {
            MIAtom *atom = res->atom(ia);
            if (atom->isHidden())
            {
                ++hiddenAtoms;
                if (shownAtoms > 0)
                {
                    hidden = false;
                    partial = true;
                    break;
                }
            }
            else
            {
                ++shownAtoms;
                hidden = false;
                if (hiddenAtoms > 0)
                {
                    partial = true;
                    break;
                }
            }
        }
        res = res->next();
    }

    int image = 5;
    if (residue == currentChain)
    {
        if (hidden)
        {
            image = 10;
        }
        else if (partial)
        {
            image = 8;
        }
        else
        {
            image = 6;
        }
    }
    else
    {
        if (hidden)
        {
            image = 9;
        }
        else if (partial)
        {
            image = 7;
        }
        else
        {
            image = 5;
        }
    }

    std::string name = chainstring(residue).c_str();
    std::string text = truncateRight(name, 24);
    item->setText(0, text.c_str());
    item->setIcon(0, _imageList[ image]);
    item->setIcon(0, _imageList[ image ]);
}

void ModelsTree::addMap(EMap *map)
{
    TreeData *data = new TreeData;
    data->map = map;
    QTreeWidgetItem *item = appendItem(rootId, map->mapName.c_str(), 2, 2, data);
    if (!item)
    {
        return;
    }
    mapToData[map] = data;
    connect(map, SIGNAL(mapFftRecalculated(EMap*)),
            this, SLOT(mapFftRecalculated(EMap*)));
    stylizeMap(map);

    addMapCrystal(map);

    expandItem(item);
}

void ModelsTree::addMapCrystal(EMap *map)
{
    CMapHeaderBase *mapHeader = map->GetMapHeader();
    if (mapHeader != NULL)
    {
        connect(mapHeader, SIGNAL(mapHeaderChanged(CMapHeaderBase*)),
                this, SLOT(mapHeaderChanged(CMapHeaderBase*)));
        TreeData *mapData = mapToData[map];
        if (mapData == NULL || !validTreeData(mapData))
        {
            mapToData.erase(map);
            return;
        }
        QTreeWidgetItem *mapItem = mapData->GetId();
        if (!mapItem)
        {
            return;
        }
        TreeData *data = new TreeData;
        data->mapHeader = mapHeader;
        QTreeWidgetItem *item = appendItem(mapItem, mapHeader->Label().c_str(), 4, 4, data);
        if (!item)
        {
            return;
        }
        mapToCrystalData[map] = data;
        mapHeaderToCrystalData[mapHeader] = data;
        stylizeCrystal(mapHeader);
    }
}

void ModelsTree::stylizeCrystal(CMapHeaderBase *mapHeader)
{
    std::string text = mapHeader->Label().c_str();
    text = truncateRight(text, 16);
    TreeData *data = mapHeaderToCrystalData[mapHeader];
    if (data == NULL || !validTreeData(data))
    {
        mapHeaderToCrystalData.erase(mapHeader);
        return;
    }
    QTreeWidgetItem *item = data->GetId();
    if (!item)
    {
        return;
    }
    item->setText(0, text.c_str());
}

void ModelsTree::mapHeaderChanged(CMapHeaderBase *mapHeader)
{
    if (mapHeaderToCrystalData.find(mapHeader) != mapHeaderToCrystalData.end())
    {
        stylizeCrystal(mapHeader);
    }
}

void ModelsTree::mapFftRecalculated(EMap *map)
{
    stylizeMap(map);
}

void ModelsTree::modelAdded(Molecule *model)
{
    addModel(model);
    update();
}

void ModelsTree::modelToBeDeleted(MIMoleculeBase *mol)
{
    Molecule *model = (Molecule*)mol;
    TreeData *data = modelToData[model];
    if (!validTreeData(data))
    {
        return;
    }
    QTreeWidgetItem *item = data->GetId();
    if (!item)
    {
        return;
    }
    modelToData.erase(model);
    if (model == currentModel)
    {
        setCurrentChain(NULL);
    }
    for (Residue *res = model->residuesBegin(); Residue::isValid(res); res = res->next())
    {
        residueToChainData.erase(res);
    }
    mapHeaderToCrystalData.erase(&(model->GetMapHeader()));
    Delete(item);
    update();
}

void ModelsTree::currentModelChanged(Molecule *oldModel, Molecule *newModel)
{
    TreeData *data;
    QTreeWidgetItem *item;
    if (oldModel != NULL)
    {
        data = modelToData[oldModel];
        if (data != NULL)
        {
            if (!validTreeData(data))
            {
                modelToData.erase(oldModel);
            }
            else
            {
                item = data->GetId();
                if (item)
                {
                    item->setIcon(0, _imageList[ 0]);
                    item->setIcon(0, _imageList[ 0 ]);
                }
            }
        }
    }
    if (newModel != NULL)
    {
        data = modelToData[newModel];
        if (data != NULL)
        {
            if (!validTreeData(data))
            {
                modelToData.erase(oldModel);
            }
            else
            {
                item = data->GetId();
                if (item)
                {
                    item->setIcon(0, _imageList[ 1]);
                    item->setIcon(0, _imageList[ 1 ]);
                }
            }
        }
    }
}

void ModelsTree::mapAdded(EMap *map)
{
    addMap(map);
    update();
}

void ModelsTree::mapToBeDeleted(EMap *map)
{
    TreeData *data = mapToData[map];
    if (data == NULL)
    {
        return;
    }
    if (!validTreeData(data))
    {
        mapToData.erase(map);
        mapToCrystalData.erase(map);
        mapHeaderToCrystalData.erase(map->GetMapHeader());
        return;
    }
    QTreeWidgetItem *item = data->GetId();
    if (!item)
    {
        return;
    }
    Delete(item);
    mapToData.erase(map);
    mapToCrystalData.erase(map);
    mapHeaderToCrystalData.erase(map->GetMapHeader());
    update();
}

void ModelsTree::currentMapChanged(EMap *oldMap, EMap *newMap)
{
    TreeData *data;
    QTreeWidgetItem *item;
    if (oldMap != NULL)
    {
        data = mapToData[oldMap];
        if (data == NULL || !validTreeData(data))
        {
            mapToData.erase(oldMap);
        }
        else
        {
            item = data->GetId();
            if (item)
            {
                item->setIcon(0, _imageList[ 2]);
                item->setIcon(0, _imageList[ 2 ]);
            }
        }
        stylizeMap(oldMap);
    }
    if (newMap != NULL)
    {
        data = mapToData[newMap];
        if (data == NULL || !validTreeData(data))
        {
            mapToData.erase(newMap);
        }
        else
        {
            item = data->GetId();
            if (item)
            {
                item->setIcon(0, _imageList[ 3]);
                item->setIcon(0, _imageList[ 3 ]);
            }
        }
        stylizeMap(newMap);
    }
}

bool ModelsTree::colorDialog(int &color, int &colorMethod)
{
    GenericDataDialog dlg(this);
    dlg.setWindowTitle("Choose Color");
    dlg.addColorIndexField("Color:", view->WhenShownColor);

    QStringList methods;
    methods << "Carbon Only" << "All Atoms" << "Secondary Structure" << "B-Value" << "Atom Type"
            << "Hydrophobicity" << "Shapley";
    dlg.addComboField("Method:", methods, view->WhenShownColorMethod);

    if (dlg.exec() != QDialog::Accepted)
    {
        return false;
    }
    color = dlg.value(0).toInt();
    colorMethod = dlg.value(1).toInt();
    return true;
}


ModelsView::LineEditToModelsTreeMap ModelsView::lineEditToModelsTree;
ModelsView::PanelToLineEditMap ModelsView::panelToLineEdit;
ModelsView::PanelToToolButtonMap ModelsView::panelToToolButton;
ModelsView::ButtonCtrlToModelsTreeMap ModelsView::buttonCtrlToModelsTree;
ModelsView::ButtonCtrlToLineEditMap ModelsView::buttonCtrlToLineEdit;
ModelsView::PanelToButtonCtrlMap ModelsView::panelToButtonCtrl;

class ModelsViewPanel
    : public QWidget
{
public:

    ModelsTree *modelsTree;
    ResiduesTree *residuesTree;
    AtomsTree *atomsTree;

    ModelsViewPanel(QWidget *parent);
    virtual ~ModelsViewPanel();
};

ModelsViewPanel::ModelsViewPanel(QWidget *parent)
    : QWidget(parent),
      modelsTree(NULL),
      residuesTree(NULL),
      atomsTree(NULL)
{
}

ModelsViewPanel::~ModelsViewPanel()
{
    //printf("In ModelsViewPanel dtor\n");
}

ModelsView::ModelsView(QWidget *parent)
    : ViewSyncedPanel(parent)
{
}

ModelsView::~ModelsView()
{
    //printf("In ModelsView dtor\n");
}


class MyLineEdit
    : public QLineEdit
{

    Q_OBJECT
public:
    MyLineEdit(const QString &str, QWidget *parent)
        : QLineEdit(str, parent)
    {
    }

    void focusInEvent(QFocusEvent *evt)
    {
        if (text() == QString("Go to residue"))
            clear();
        selectAll();
        QLineEdit::focusInEvent(evt);
    }

    void enterEvent(QEvent *evt)
    {
        setFocus(Qt::MouseFocusReason);
        QLineEdit::enterEvent(evt);
    }

};


QWidget*ModelsView::createPanelForView(MIGLWidget *view, QWidget *parent)
{
    parent->layout()->setContentsMargins(0, 0, 0, 0);

    ModelsViewPanel *panel = new ModelsViewPanel(parent);

    // create vertical layout
    QHBoxLayout *hlayout = new QHBoxLayout();
    hlayout->setContentsMargins(0, 0, 0, 0);
    hlayout->setSpacing(2);

    QLineEdit *goToResidueLineEdit = new MyLineEdit("Go to residue", panel);
    goToResidueLineEdit->setToolTip("Go to residue (A 123)");

    panelToLineEdit[panel] = goToResidueLineEdit;

    QPushButton *gotoButton = new QPushButton("GoTo", panel);
    panelToButtonCtrl[panel] = gotoButton;

    QToolButton *toolButton = new QToolButton(panel);
    toolButton->setToolTip("Sync with selection");
    toolButton->setIcon(QIcon(QPixmap(synced_xpm)));
    toolButton->setCheckable(true);
    toolButton->setChecked(MIConfig::Instance()->GetProfileInt("ModelsView", "syncView", 1) != 0);
    connect(toolButton, SIGNAL(toggled(bool)), this, SLOT(updateSyncView(bool)));

    panelToToolButton[panel] = toolButton;
    updateSyncView(toolButton);

    QSpacerItem *spacer = new QSpacerItem(20, 5, QSizePolicy::Expanding);
    hlayout->addItem(spacer);
    hlayout->addWidget(goToResidueLineEdit);
    hlayout->addWidget(gotoButton);
    hlayout->addWidget(toolButton);

    QVBoxLayout *vbox = new QVBoxLayout(panel);
    vbox->addLayout(hlayout);
    vbox->setContentsMargins(0, 0, 0, 0);
    vbox->setSpacing(2);

    QSplitter *splitter = new QSplitter(Qt::Vertical, panel);
    splitter->setContentsMargins(0, 0, 0, 0);
    vbox->addWidget(splitter);

    panel->setLayout(vbox);
    ModelsTree *modelsTree = new ModelsTree(splitter);
    ResiduesTree *residuesTree = new ResiduesTree(splitter);
    AtomsTree *atomsTree = new AtomsTree(splitter);

    modelsTree->setView(view);
    modelsTree->setTruncateWidth(width());
    residuesTree->setView(view);
    atomsTree->setView(view);

    splitter->addWidget(modelsTree);
    splitter->addWidget(residuesTree);
    splitter->addWidget(atomsTree);

    connect(splitter, SIGNAL(splitterMoved(int, int)),
            this, SLOT(OnSplitterChanged(int, int)));

    QSettings *settings = MIGetQSettings(); // could use MIConfig (it's the same file), but this api is easier here
    if (settings->contains("ModelsViewSplitter"))
        splitter->restoreState(settings->value("ModelsViewSplitter", splitter->saveState()).toByteArray());

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

void ModelsView::destroyContentsForView(MIGLWidget*, QWidget *panel)
{
    if (panelToToolButton.find(panel) != panelToToolButton.end())
    {
        panelToToolButton.erase(panel);
    }
    if (panelToLineEdit.find(panel) != panelToLineEdit.end())
    {
        QLineEdit *lineEdit = panelToLineEdit[panel];
        panelToLineEdit.erase(panel);
        lineEditToModelsTree.erase(lineEdit);
    }
    if (panelToButtonCtrl.find(panel) != panelToButtonCtrl.end())
    {
        QPushButton *buttonCtrl = panelToButtonCtrl[panel];
        panelToButtonCtrl.erase(panel);
        buttonCtrlToModelsTree.erase(buttonCtrl);
        buttonCtrlToLineEdit.erase(buttonCtrl);
    }
}

void ModelsView::updateSyncView(bool state)
{
    syncView = state;
    PanelToToolButtonMap::iterator iter = panelToToolButton.begin();
    while (iter != panelToToolButton.end())
    {
        QToolButton *tb = iter->second;
        ++iter;
        tb->setChecked(syncView);
    }
    MIConfig::Instance()->WriteProfileInt("ModelsView", "syncView", syncView);
}

void ModelsView::OnSplitterChanged(int, int)
{
    QSplitter *splitter = dynamic_cast<QSplitter*>(sender());
    if (!splitter)
    {
        //printf("Splitter changed signal not from splitter\n");
        return;
    }
    QSettings *settings = MIGetQSettings(); // could use MIConfig (it's the same file), but this api is easier here
    settings->setValue("ModelsViewSplitter", splitter->saveState());
}

void ModelsView::resizeEvent(QResizeEvent *event)
{
    LineEditToModelsTreeMap::iterator iter;
    for (iter = lineEditToModelsTree.begin(); iter != lineEditToModelsTree.end(); ++iter)
    {
        ModelsTree *modelsTree = iter->second;
        modelsTree->setTruncateWidth(event->size().width());
        modelsTree->restylize();
    }
    ViewSyncedPanel::resizeEvent(event);
}

void ModelsView::OnGoToResidueReturnPressed()
{
    QLineEdit *lineEdit = dynamic_cast<QLineEdit*>(sender());
    if (lineEdit)
    {
        lineEditToModelsTree[lineEdit]->goToResidue(lineEdit->text().toStdString());
        lineEdit->selectAll();
        return;
    }
    QPushButton *buttonCtrl = dynamic_cast<QPushButton*>(sender());
    if (buttonCtrl)
    {
        QLineEdit *lineEdit = buttonCtrlToLineEdit[buttonCtrl];
        buttonCtrlToModelsTree[buttonCtrl]->goToResidue(lineEdit->text().toStdString());
        lineEdit->selectAll();
        return;
    }
}

ModelsTree*ModelsView::GetCurrentModelsTree() const
{
    ModelsViewPanel *panel = dynamic_cast<ModelsViewPanel*>(stackedLayout->currentWidget());
    if (panel == NULL)
    {
        return NULL;
    }
    return panel->modelsTree;
}

ResiduesTree*ModelsView::GetCurrentResiduesTree() const
{
    ModelsViewPanel *panel = dynamic_cast<ModelsViewPanel*>(stackedLayout->currentWidget());
    if (panel == NULL)
    {
        return NULL;
    }
    return panel->residuesTree;
}

AtomsTree*ModelsView::GetCurrentAtomsTree() const
{
    ModelsViewPanel *panel = dynamic_cast<ModelsViewPanel*>(stackedLayout->currentWidget());
    if (panel == NULL)
    {
        return NULL;
    }
    return panel->atomsTree;
}

#include "ModelsView.moc"
