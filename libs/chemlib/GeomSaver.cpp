#include "model.h"
#include "Bond.h"
#include "MIMoleculeBase.h"
#include "GeomSaver.h"
#include <chemlib/RESIDUE_.h>


using namespace std;

namespace chemlib
{

GeomSaver::GeomSaver()
{
    // the first position is empty an a value of 0 for a save token
    // indicates an error condition
    SaveSets.push_back(SaveItem());
}

GeomSaver::~GeomSaver()
{
}

unsigned int GeomSaver::Save(Residue *res, int nres, MIMoleculeBase *model)
{
    SaveSets.push_back(SaveItem(model, std::string("")));
    unsigned int token = SaveSets.size()-1;
    int n = 0;
    Residue *lastres = res;
    Residue *firstres = res;
    vector<SaveAtom> *Save = &(SaveSets[token].SaveSet);
    MIAtom_const_iter atom, endAtom;
    for (MIIter<Residue> resIter = Residue::getIter(res); resIter && n < nres; ++resIter)
    {
        const MIAtomList &atoms = resIter->atoms();
        endAtom = atoms.end();
        for (atom = atoms.begin(); atom != endAtom; ++atom)
        {
            Save->push_back(SaveAtom(*atom));
        }
        n++;
        lastres = resIter;
    }
    char buf[1024];
    sprintf(buf, "%s_%s%c to %s_%s%c, %d atoms, Model: %s",
            firstres->type().c_str(), firstres->name().c_str(), (char)(firstres->chain_id()&255),
            lastres->type().c_str(), lastres->name().c_str(), (char)(lastres->chain_id()&255),
            Save->size(), model->compound.c_str());
    SaveSets[token].Title = buf;
    return unique();
}

unsigned int GeomSaver::Save(const MIAtomList &atoms, MIMoleculeBase *model)
{
    SaveSets.push_back(SaveItem(model, std::string("")));
    unsigned int token = SaveSets.size()-1;
    unsigned int i;
    vector<SaveAtom> *Save = &(SaveSets[token].SaveSet);
    for (i = 0; i < atoms.size(); i++)
    {
        Save->push_back(atoms[i]);
    }
    char buf[1024];
    sprintf(buf, "%d atoms, Model: %s",
            Save->size(), model->compound.c_str());
    SaveSets[token].Title = buf;
    return unique();
}

int GeomSaver::Restore(unsigned int token) const
{
    if (token < 1 || token > SaveSets.size())
    {
        return 0;
    }
    const std::vector<SaveAtom> *Save = &(SaveSets[token].SaveSet);
    if (!Save)
    {
        return 0;
    }
    for (unsigned int i = 0; i < Save->size(); i++)
    {
        (*Save)[i].Restore();
    }
    return 1;
}

bool GeomSaver::RestoreColor(unsigned int token, unsigned int mask)
{
    if (token < 1 || token > SaveSets.size())
    {
        return false;
    }
    std::vector<SaveAtom> *Save = &(SaveSets[token].SaveSet);
    if (!Save)
    {
        return false;
    }
    for (unsigned int i = 0; i < Save->size(); i++)
    {
        (*Save)[i].RestoreColor(mask);
    }
    return true;
}

void GeomSaver::Purge(MIAtom *atom)
{
    std::vector<SaveAtom> *Save;
    for (unsigned int i = 0; i < SaveSets.size(); i++)
    {
        Save = &(SaveSets[i].SaveSet);
        unsigned int j = 0;
        while (j < Save->size())
        {
            if ((*Save)[j].from == atom)
            {
                (*Save).erase((*Save).begin()+j);
            }
            else
            {
                j++;
            }
        }
    }
}

void GeomSaver::Purge(MIMoleculeBase *model)
{
    for (unsigned int i = 0; i < SaveSets.size(); i++)
    {
        MIMoleculeBase *m_test = SaveSets[i].SaveMolecule;
        if (m_test == model)
        {
            SaveSets[i].Title = "Empty";
            SaveSets[i].SaveMolecule = NULL;
            SaveSets[i].SaveSet.clear();
        }
    }
}

unsigned int GeomSaver::RestoreLast(MIMoleculeBase *model)
{
    for (unsigned int i = SaveSets.size()-1; i > 0; i--)
    {
        MIMoleculeBase *m_test = SaveSets[i].SaveMolecule;
        if (m_test == model)
        {
            Restore(i);
            return i;
        }
    }
    return 0;
}

unsigned int GeomSaver::unique()
{
    // check that the list is unique and remove adjacent duplicates
    unsigned int i = 1;
    while (i < SaveSets.size()-1)
    {
        if (SaveSets[i].SaveMolecule == SaveSets[i+1].SaveMolecule
            && SaveSets[i].SaveSet == SaveSets[i+1].SaveSet)
        {
            SaveSets.erase(SaveSets.begin()+i);
            continue;
        }
        i++;
    }
    return SaveSets.size()-1;
}

MIMoleculeBase*GeomSaver::Model(unsigned int token)
{
    if (token > 0 && token < SaveSets.size())
    {
        return SaveSets[token].SaveMolecule;
    }
    return NULL;
}

void GeomSaver::Clear()
{
    SaveSets.clear();
    SaveSets.push_back(SaveItem());
}

}


