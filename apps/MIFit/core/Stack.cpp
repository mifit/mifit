#include <chemlib/Monomer.h>
#include "Stack.h"
#include "Molecule.h"


using namespace chemlib;

Stack::Stack()
    : visible_(true)
{
}

Stack::~Stack()
{
}

void Stack::moleculeToBeDeleted(MIMoleculeBase *mol)
{
    removeAll(mol);

    // removeAll may not delete everything, if there was an atom pushed on the
    // stack where mol and/or res was null, the stack might still have a ref
    // to it, so we have to try harder
    std::vector<Residue*> residues;
    ResidueListIterator res = mol->residuesBegin();
    for (; res != mol->residuesEnd(); ++res)
        residues.push_back(res);

    for (res = mol->symmResiduesBegin(); res != mol->symmResiduesEnd(); ++res)
        residues.push_back(res);

    residuesToBeDeleted(mol, residues);
}

void Stack::residuesToBeDeleted(MIMoleculeBase *mol, std::vector<Residue*> &residues)
{
    // removeAll may not delete everything, if there was an atom pushed on the
    // stack where mol and/or res was null, the stack might still have a ref
    // to it, so we have to try harder
    MIAtomList atoms;
    foreach (Residue *residue, residues)
    {
        removeAll(residue);
        atoms.insert(atoms.end(), residue->atoms().begin(), residue->atoms().end());
    }
    atomsToBeDeleted(mol, atoms);
}

void Stack::atomsToBeDeleted(MIMoleculeBase*, const MIAtomList &atoms)
{
    foreach (MIAtom *atom, atoms)
        removeAll(atom);
}


void Stack::Push(MIAtom *atom, Residue *res, Molecule *m)
{
    if (!MIAtom::isValid(atom) || !Monomer::isValid(res) || !MIMoleculeBase::isValid(m))
        return;

    StackItem item;
    item.atom = atom;
    item.residue = res;
    item.molecule = m;
    data.push_back(item);
    sizeChanged();
}

void Stack::Pop(MIAtom* &atom, Residue* &res)
{
    if (data.size() == 0)
    {
        atom = NULL;
        res = NULL;
        return;
    }
    StackItem item = data.back();
    dataPop();
    atom = item.atom;
    res = item.residue;
}

void Stack::Pop(MIAtom* &atom, Residue* &res, Molecule* &m)
{
    if (data.size() == 0)
    {
        atom = NULL;
        res = NULL;
        m = NULL;
        return;
    }
    StackItem item = data.back();
    dataPop();
    atom = item.atom;
    res = item.residue;
    m = item.molecule;
}

MIAtom *Stack::Pop()
{
    if (data.size() == 0)
        return NULL;

    StackItem item = data.back();
    dataPop();
    return item.atom;
}

void Stack::dataPop()
{
    bool wasEmpty = empty();
    data.pop_back();
    if (!wasEmpty)
        emit sizeChanged();
}

void Stack::Peek(MIAtom* &atom, Residue* &res, Molecule* &m)
{
    if (data.size() == 0)
    {
        atom = NULL;
        res = NULL;
        m = NULL;
        return;
    }
    StackItem item = data.back();
    atom = item.atom;
    res = item.residue;
    m = item.molecule;
}

bool Stack::contains(Residue *res)
{
    foreach (StackItem item, data)
        if (item.residue == res)
            return true;

    return false;
}

void Stack::removeAll(MIMoleculeBase *model)
{
    bool changed = false;
    DataContainer::iterator iter = data.begin();
    while (iter != data.end())
    {
        if (iter->molecule == model)
        {
            changed = true;
            iter = data.erase(iter);
        }
        else
            ++iter;
    }
    if (changed)
        emit sizeChanged();
}

void Stack::removeAll(Residue *res)
{
    bool changed = false;
    DataContainer::iterator iter = data.begin();
    while (iter != data.end())
    {
        if (iter->residue == res)
        {
            changed = true;
            iter = data.erase(iter);
        }
        else
            ++iter;
    }
    if (changed)
        emit sizeChanged();
}

void Stack::removeAll(MIAtom *atom)
{
    bool changed = false;
    DataContainer::iterator iter = data.begin();
    while (iter != data.end())
    {
        if (iter->atom == atom)
        {
            changed = true;
            iter = data.erase(iter);
        }
        else
            ++iter;
    }
    if (changed)
        emit sizeChanged();
}

int Stack::size()
{
    return data.size();
}

bool Stack::empty()
{
    return data.empty();
}

const Stack::DataContainer &Stack::getData()
{
    return data;
}

QStringList Stack::toStringList() const
{
    QStringList list;
    DataContainer::const_reverse_iterator iter = data.rbegin();
    int n = 0;
    while (iter != data.rend() && n < 4)
    {
        StackItem item = *iter;
        ++iter;
        if (MIAtom::isValid(item.atom) && Monomer::isValid(item.residue) && MIMoleculeBase::isValid(item.molecule) && item.molecule->Visible())
        {
            ++n;
            list += tr("%1: %2 %3 %4")
                    .arg(n)
                    .arg(item.residue->type().c_str())
                    .arg(item.residue->name().c_str())
                    .arg(item.atom->name());
        }
    }
    if (data.size() > 4)
        list += tr(" + %1 more...").arg(data.size() - 4);

    return list;
}

void Stack::pop()
{
    dataPop();
}

void Stack::clear()
{
    bool wasEmpty = empty();
    DataContainer().swap(data);
    if (!wasEmpty)
        emit sizeChanged();
}

bool Stack::isVisible() const
{
    return visible_;
}

void Stack::setVisible(bool visible)
{
    if (visible_ == visible)
        return;

    visible_ = visible;
    emit visibilityChanged(visible_);
}
