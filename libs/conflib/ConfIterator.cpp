#include <chemlib/chemlib.h>
#include <chemlib/Monomer.h>

#include "ConfIterator.h"
#include "bumps.h"

using namespace chemlib;
using namespace std;

namespace conflib
{

//ConfIterator::ConfIterator(RESIDUE *res, vector<Bond> &bonds) {}
//ConfIterator::~ConfItereator() {}

ConfEnumerator::ConfEnumerator(Residue *res, vector<Bond> &bonds, vector<TORSION> &torsions)
    : _res(res)
{

    //Create flexible torsions
    unsigned int i;
    for (i = 0; i < torsions.size(); ++i)
    {
        if (torsions[i].nideal == 0)
        {
            _flexors.push_back(conflib::FlexTorsion(torsions[i], bonds, false));
        }
    }

    //Create interatomic bumps
    GetBumps(res, _bumps, bonds);

    //Compute the number of theoretical conformations
    _nTheory = 1;
    for (unsigned int i = 0; i < _flexors.size(); ++i)
    {
        _nTheory *= _flexors[i].NumAngles();
    }


    //Mark all of the atoms as modified (i.e. none of the atom pairs are guaranteed bump-free)
    for (int i = 0; i < _res->atomCount(); ++i)
    {
        _res->atom(i)->set_search_flag(1);
    }

    _nTry = 0;
    _init = false;
}

ConfEnumerator::~ConfEnumerator()
{
}

bool ConfEnumerator::Next()
{
    if (_flexors.size() == 0)
    {
        return false;
    }

    ConfSaver orig(_res);
    orig.Save();

    if (_init == false)
    {
        _init = true;

        //Set the initial conformation
        for_each(_flexors.begin(), _flexors.end(), void_mem_fun_ref(&conflib::FlexTorsion::Set));
    }

    vector<FlexTorsion>::iterator ft;
    for ( ; _nTry < _nTheory; ++_nTry)
    {

        //Create the next conformation by advancing torsions to
        //their next value...usually we will only increment the first one, but
        //if Advance() returns false, that torsion has "rolled over" to its
        //initial value, and we need to increment the next one, and so on
        for (ft = _flexors.begin(); ft != _flexors.end() && !ft->Advance(); ft++)
        {
        }

        //Check for clashes or redundancy. Our work is done if there are none.
        if (CheckConf())
        {
            for (int i = 0; i < _res->atomCount(); ++i)
            {
                _res->atom(i)->set_search_flag(0);
            }
            return true;
        }
    }

    //if we reach here, we went through all the confs without finding a bump-free one
    orig.Restore(1);
    return false;
}

int ConfEnumerator::GenerateConfs(ConfSaver &confs, int max)
{
    if (_flexors.size() == 0)
    {
        return 0;
    }

    //Save the original conformation, since we will be mucking with it.
    ConfSaver orig(_res);
    orig.Save();

    //Set all the torsion angles in the molecule to those currently held in
    //the flexors
    _init = true;
    for_each(_flexors.begin(), _flexors.end(), void_mem_fun_ref(&conflib::FlexTorsion::Set));

    int nFound = 0;
    vector<FlexTorsion>::iterator ft;
    for ( ; _nTry < _nTheory && nFound < max; ++_nTry)
    {

        //Create the next conformation by advancing torsions to
        //their next value...usually we will only increment the first one, but
        //if Advance() returns false, that torsion has "rolled over" to its
        //initial value, and we need to increment the next one, and so on
        for (ft = _flexors.begin(); ft != _flexors.end() && !ft->Advance(); ft++)
        {
        }

        if (CheckConf())
        {
            confs.Save();
            nFound++;

            for (int i = 0; i < _res->atomCount(); ++i)
            {
                _res->atom(i)->set_search_flag(0);
            }
        }
    }

    //Restore the initial conformation
    orig.Restore(1);
    return nFound;
}

bool ConfEnumerator::CheckConf()
{
    std::vector<Bond>::iterator i = _bumps.begin();
    std::vector<Bond>::iterator e = _bumps.end();
    while (i < e)
    {
        //Only check atom pair for bumps if at least one has been moved
        if ((i->getAtom1()->search_flag() || i->getAtom2()->search_flag())
            && SquaredAtomDist(*i->getAtom1(), *i->getAtom2()) < i->ideal_length)
        {
            return false;
        }
        i++;
    }

    ConfFingerprint fp(_bumps);
    if (find(_prints.begin(), _prints.end(), fp) != _prints.end())
    {
        return false;
    }

    _prints.push_back(fp);
    return true;
}

ConfSampler::ConfSampler(Residue *res, vector<Bond> &bonds, vector<TORSION> &torsions)
    : _res(res)
{

    //Create flexible torsions
    unsigned int i;
    for (i = 0; i < torsions.size(); ++i)
    {
        if (torsions[i].nideal == 0)
        {
            _flexors.push_back(conflib::FlexTorsion(torsions[i], bonds, false));
        }
    }

    //Create interatomic bumps
    GetBumps(res, _bumps, bonds);

    //Mark all of the atoms as modified (i.e. none of the atom pairs are guaranteed bump-free)
    for (int i = 0; i < _res->atomCount(); ++i)
    {
        _res->atom(i)->set_search_flag(1);
    }
}

ConfSampler::~ConfSampler()
{
}

bool ConfSampler::Next()
{
    if (_flexors.size() == 0)
    {
        return false;
    }

    ConfSaver orig(_res);
    orig.Save();

    vector<FlexTorsion>::iterator ft;
    int i, max = 1000;
    for (i = 0; i < max; ++i)
    {

        //Create the next conformation
        for (ft = _flexors.begin(); ft != _flexors.end(); ft++)
        {
            ft->Pick(irand_approx(ft->NumAngles()));
        }

        //Check for clashes or redundancy. Our work is done if there are none.
        if (CheckConf())
        {
            for (int i = 0; i < _res->atomCount(); ++i)
            {
                _res->atom(i)->set_search_flag(0);
            }
            return true;
        }
    }

    //if we reach here, we went through max attempts without finding a bump-free one
    orig.Restore(1);
    return false;
}

int ConfSampler::GenerateConfs(ConfSaver &confs, int max)
{
    if (_flexors.size() == 0)
    {
        return 0;
    }

    //Save the original conformation, since we will be mucking with it.
    ConfSaver orig(_res);
    orig.Save();

    int n, nFound = 0;
    vector<FlexTorsion>::iterator ft;
    for (n = 0 ; n < max; ++n)
    {

        //Create the next conformation
        for (ft = _flexors.begin(); ft != _flexors.end(); ft++)
        {
            ft->Pick(irand_approx(ft->NumAngles()));
        }

        if (CheckConf())
        {
            confs.Save();
            nFound++;

            for (int i = 0; i < _res->atomCount(); ++i)
            {
                _res->atom(i)->set_search_flag(0);
            }
        }
    }

    //Restore the initial conformation
    orig.Restore(1);
    return nFound;
}

bool ConfSampler::CheckConf()
{
    std::vector<Bond>::iterator i = _bumps.begin();
    std::vector<Bond>::iterator e = _bumps.end();
    while (i < e)
    {
        //Only check atom pair for bumps if at least one has been moved
        if ((i->getAtom1()->search_flag() || i->getAtom2()->search_flag())
            && SquaredAtomDist(*i->getAtom1(), *i->getAtom2()) < i->ideal_length)
        {
            return false;
        }
        i++;
    }

    ConfFingerprint fp(_bumps);
    if (find(_prints.begin(), _prints.end(), fp) != _prints.end())
    {
        return false;
    }

    _prints.push_back(fp);
    return true;
}

} // namespace conflib
