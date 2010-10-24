#ifndef mifit_model_Stack_h
#define mifit_model_Stack_h

class Stack;

#include <deque>
#include <QtCore/QObject>
#include <QtCore/QStringList>
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
    Q_PROPERTY(bool visible READ isVisible WRITE setVisible NOTIFY visibilityChanged)

public:
    typedef std::deque<StackItem> DataContainer;

    Stack();
    ~Stack();

    bool empty();
    int size();
    bool contains(chemlib::Residue *res);

    void Push(chemlib::MIAtom *atom, chemlib::Residue *res, Molecule *m);
    void Peek(chemlib::MIAtom *&atom, chemlib::Residue *&res, Molecule *&m);

    void Pop(chemlib::MIAtom *&atom, chemlib::Residue *&res);
    void Pop(chemlib::MIAtom *&atom, chemlib::Residue *&res, Molecule *&m);
    chemlib::MIAtom *Pop();

    const DataContainer &getData();

    QStringList toStringList() const;

    bool isVisible() const;

    // note: in practice these slots are called by CMolwView, not directly by the signal
    void residuesToBeDeleted(chemlib::MIMoleculeBase *model, std::vector<chemlib::Residue*> &residues);
    void atomsToBeDeleted(chemlib::MIMoleculeBase *model, const chemlib::MIAtomList &atoms);
    void moleculeToBeDeleted(chemlib::MIMoleculeBase *model);

public slots:
    void pop();
    void clear();
    void setVisible(bool visible);

signals:
    void sizeChanged();
    void visibilityChanged(bool);

private:
    DataContainer data;
    bool visible_;

    void dataPop();
    void removeAll(chemlib::MIMoleculeBase *model);
    void removeAll(chemlib::Residue *res);
    void removeAll(chemlib::MIAtom *atom);

};

#endif // ifndef mifit_model_Stack_h

