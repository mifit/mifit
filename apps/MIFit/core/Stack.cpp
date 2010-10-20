#include <stdio.h>

#include <nongui/nonguilib.h>
#include <chemlib/chemlib.h>
#include <chemlib/Monomer.h>

#include "Stack.h"
#include "RESIDUE.h"
#include "Molecule.h"


using namespace chemlib;

Stack::Stack()
{
    changed = false;
    minimized = false;
}

Stack::~Stack()
{
}

void Stack::moleculeToBeDeleted(MIMoleculeBase *mol)
{
    Purge(mol);

    // Purge may not delete everything, if there was an atom pushed on the
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
    std::vector<Residue*>::iterator iter;
    for (iter = residues.begin(); iter != residues.end(); ++iter)
    {
        Residue *residue = *iter;
        Purge(residue);
    }

    // Purge may not delete everything, if there was an atom pushed on the
    // stack where mol and/or res was null, the stack might still have a ref
    // to it, so we have to try harder
    MIAtomList atoms;
    for (size_t i = 0; i < residues.size(); ++i)
    {
        atoms.insert(atoms.end(), residues[i]->atoms().begin(), residues[i]->atoms().end());
    }
    atomsToBeDeleted(mol, atoms);

}

void Stack::atomsToBeDeleted(MIMoleculeBase*, const MIAtomList &atoms)
{
    for (unsigned int i = 0; i<atoms.size(); ++i)
    {
        Purge(atoms[i]);
    }
}


void Stack::Push(MIAtom *atom, Residue *res, Molecule *m)
{
    if (!MIAtom::isValid(atom) || !Monomer::isValid(res) || !MIMoleculeBase::isValid(m))
    {
        return;
    }
    StackItem item;
    item.atom = atom;
    item.residue = res;
    item.molecule = m;
    bool wasEmpty = empty();
    data.push_back(item);
    sizeChanged();
    changed = true;
    if (wasEmpty != empty())
    {
        emptyChanged(empty());
    }
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
    {
        return NULL;
    }
    StackItem item = data.back();
    dataPop();
    return item.atom;
}

void Stack::dataPop()
{
    bool wasEmpty = empty();
    data.pop_back();
    changed = true;
    emit sizeChanged();
    if (wasEmpty != empty())
    {
        emptyChanged(empty());
    }
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

bool Stack::StackChanged()
{
    return changed;
}

void Stack::ClearChanged()
{
    changed = false;
}

void Stack::ExpandTopAllAtoms()
{
    MIAtom *a;
    Molecule *m;
    Residue *r;
    Pop(a, r, m);

    if (a && a->type() & AtomType::SYMMATOM)
    {
        Push(a, r, m);
        Logger::message("Can't expand symmetry atom(s) on stack");
        return;
    }

    for (int i = 0; i < r->atomCount(); i++)
    {
        Push(r->atom(i), r, m);
    }
}

void Stack::ExpandTop2AllAtoms()
{
    MIAtom *a1;
    MIAtom *a2;
    Molecule *m1, *m2;
    Residue *r1, *r2;
    Pop(a1, r1, m1);
    Pop(a2, r2, m2);
    if (m1 != m2)
    {
        Logger::message("Both atoms must be in same molecule");
        Push(a2, r2, m2);
        Push(a1, r1, m1);
        return;
    }

    if ((a1 && a1->type() & AtomType::SYMMATOM)
        || (a2 && a2->type() & AtomType::SYMMATOM))
    {
        Logger::message("Can't expand symmetry atom(s) on stack");
        Push(a2, r2, m2);
        Push(a1, r1, m1);
        return;
    }

    if (r1 == r2)
    {
        for (int i = 0; i < r2->atomCount(); i++)
        {
            Push(r2->atom(i), r2, m2);
        }
        return;
    }
    std::vector<ResidueListIterator> marks;
    ResidueListIterator res;
    for (res = m1->residuesBegin(); res != m1->residuesBegin(); ++res)
    {
        if (res == ResidueListIterator(r1)
            || res == ResidueListIterator(r2))
        {
            marks.push_back(res);
        }
    }
    if (marks.size() != 2)
    {
        Logger::message("Error expanding stack: residue not found in model");
        Push(a2, r2, m2);
        Push(a1, r1, m1);
    }

    ++marks[1];
    for (res = marks[0]; res != marks[1]; ++res)
    {
        for (int i = 0; i < res->atomCount(); i++)
        {
            Push(res->atom(i), res, m2);
        }
    }
}

void Stack::ExpandTop2Range()
{
    MIAtom *a1, *a2, *a;
    Molecule *m1, *m2;
    Residue *r1, *r2;
    Pop(a1, r1, m1);
    Pop(a2, r2, m2);
    if (m1 != m2)
    {
        Logger::message("Both atoms must be in same molecule");
        Push(a2, r2, m2);
        Push(a1, r1, m1);
        return;
    }

    if ((a1 && a1->type() & AtomType::SYMMATOM)
        || (a2 && a2->type() & AtomType::SYMMATOM))
    {
        Logger::message("Can't expand symmetry atom(s) on stack");
        Push(a2, r2, m2);
        Push(a1, r1, m1);
        return;
    }

    if (r1 == r2)
    {
        if ((a = atom_from_name("CA", *r2)) == NULL)
        {
            Push(a, r2, m2);
        }
        else
        {
            Push(r2->atom(0), r2, m2);
        }
        return;
    }

    ResidueListIterator res;
    std::vector<ResidueListIterator> marks;
    for (res = m1->residuesBegin(); res != m1->residuesEnd(); ++res)
    {
        if (res == ResidueListIterator(r1) || res == ResidueListIterator(r2))
            marks.push_back(res);
        if (marks.size() == 2)
            break;
    }
    if (marks.size() != 2)
    {
        Logger::message("Error expanding stack: residue not found in model");
        Push(a2, r2, m2);
        Push(a1, r1, m1);
        return;
    }

    ++marks[1];
    for (res = marks[0]; res != marks[1]; ++res)
    {
        if ((a = atom_from_name("CA", *res)) != NULL)
        {
            Push(a, res, m2);
        }
        else
        {
            Push(res->atom(0), res, m2);
        }
    }
}

bool Stack::InStack(Residue *res)
{
    DataContainer::iterator iter = data.begin();
    while (iter != data.end())
    {
        StackItem item = *iter;
        if (item.residue == res)
        {
            return true;
        }
        ++iter;
    }
    return false;
}

void Stack::Purge(MIMoleculeBase *model)
{
    DataContainer::iterator iter = data.begin();
    while (iter != data.end())
    {
        StackItem item = *iter;
        if (item.molecule == model)
        {
            changed = true;
            data.erase(iter);
            iter = data.begin();
            continue;
        }
        ++iter;
    }
    if (changed)
        emit sizeChanged();
}

void Stack::Purge(Residue *res)
{
    DataContainer::iterator iter = data.begin();
    while (iter != data.end())
    {
        StackItem item = *iter;
        if (item.residue == res)
        {
            changed = true;
            data.erase(iter);
            iter = data.begin();
            continue;
        }
        ++iter;
    }
    if (changed)
        emit sizeChanged();
}

void Stack::Purge(MIAtom *atom)
{
    DataContainer::iterator iter = data.begin();
    while (iter != data.end())
    {
        StackItem item = *iter;
        if (item.atom == atom)
        {
            changed = true;
            data.erase(iter);
            iter = data.begin();
            continue;
        }
        ++iter;
    }
    if (changed)
        emit sizeChanged();
}


bool Stack::isMinimized()
{
    return minimized;
}

int Stack::size()
{
    return data.size();
}

bool Stack::empty()
{
    return data.empty();
}

StackItem Stack::top()
{
    return data.back();
}

const Stack::DataContainer &Stack::getData()
{
    return data;
}

CRect &Stack::getClearBox()
{
    return ClearBox;
}

CRect &Stack::getHideBox()
{
    return HideBox;
}

CRect &Stack::getPopBox()
{
    return PopBox;
}

void Stack::Clear()
{
    bool wasEmpty = empty();
    DataContainer().swap(data);
    if (!wasEmpty)
    {
        emit sizeChanged();
        emit emptyChanged(true);
    }
}

bool Stack::PickClearBox(int sx, int sy)
{
    return ClearBox.Within(sx, sy);
}

bool Stack::PickHideBox(int sx, int sy)
{
    return HideBox.Within(sx, sy);
}

bool Stack::PickPopBox(int sx, int sy)
{
    return PopBox.Within(sx, sy);
}

void Stack::Minimize()
{
    minimized = true;
}

void Stack::Maximize()
{
    minimized = false;
}

void Stack::ToggleMinMax()
{
    minimized = !minimized;
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
    {
        list += tr(" + %1 more...").arg(data.size() - 4);
    }
    return list;
}

void Stack::pop()
{
    dataPop();
}

void Stack::clear()
{
    Clear();
}
