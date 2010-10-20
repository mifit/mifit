#ifndef mifit_model_Stack_h
#define mifit_model_Stack_h

class Stack;

#include <stack>
#include <deque>
#include <QtCore/QObject>
#include <QtCore/QStringList>
#include "CRect.h"
#include "Molecule.h"

typedef struct
{
    chemlib::MIAtom *atom;
    chemlib::Residue *residue;
    Molecule *molecule;
} StackItem;


/**
 * Implements an atom stack for storing atom picks.
 */
class Stack : public QObject
{
    Q_OBJECT
    Q_PROPERTY(int size READ size NOTIFY sizeChanged)
    Q_PROPERTY(QStringList stringList READ toStringList NOTIFY sizeChanged);

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
    void Purge(chemlib::MIMoleculeBase *model);
    void Purge(chemlib::Residue *res);
    void Purge(chemlib::MIAtom *atom);

public:

    Stack();
    ~Stack();

    int size();
    bool empty();
    StackItem top();
    const DataContainer& getData();
    QStringList toStringList() const;

    void Pop(chemlib::MIAtom* &atom, chemlib::Residue* &res);
    void Pop(chemlib::MIAtom* &atom, chemlib::Residue* &res, Molecule* &m);
    chemlib::MIAtom *Pop();
    void Peek(chemlib::MIAtom* &atom, chemlib::Residue* &res, Molecule* &m);
    void Push(chemlib::MIAtom *atom, chemlib::Residue *res, Molecule *m);

    void Clear();
    void ExpandTopAllAtoms();
    void ExpandTop2AllAtoms();
    void ExpandTop2Range();
    bool InStack(chemlib::Residue *res);

    bool StackChanged();
    void ClearChanged();

    CRect&getClearBox();
    CRect&getHideBox();
    CRect&getPopBox();
    bool PickClearBox(int sx, int sy);
    bool PickHideBox(int sx, int sy);
    bool PickPopBox(int sx, int sy);
    void Minimize();
    void Maximize();
    bool isMinimized();
    void ToggleMinMax();

    // note: in practice these slots are called by CMolwView, not directly by the signal
    void residuesToBeDeleted(chemlib::MIMoleculeBase *model, std::vector<chemlib::Residue*> &residues);
    void atomsToBeDeleted(chemlib::MIMoleculeBase *model, const chemlib::MIAtomList &atoms);
    void moleculeToBeDeleted(chemlib::MIMoleculeBase *model);

public slots:
    void pop();
    void clear();

signals:
    void emptyChanged(bool);
    void sizeChanged();
};

#endif // ifndef mifit_model_Stack_h

