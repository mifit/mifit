#include <chemlib/chemlib.h>
#include "sdg.h"
#include "sdg_parms.h"
#include <algorithm>

using namespace chemlib;
using namespace std;

namespace conflib
{

sdgEngine::sdgEngine(vector<sdgDistance> &dists, vector<sdgVolume> &vols, const vector<MIAtom*> &atoms)
    : _distances(dists),
      _volumes(vols),
      _atoms(atoms)
{
    _nDist = _distances.size();
    _nVol = _volumes.size();
    _nSteps = sdg_params::SDG_NSTEP_FACTOR * atoms.size();

    if (_nVol > 0)
    {
        _vol_odds = std::max((double) sdg_params::SDG_MIN_VOL_ODDS, 1 - 8.0 * _nVol / (_nDist + _nVol));
    }
    else
    {
        _vol_odds = 0;
    }
}

double sdgEngine::DoOptimize()
{
    float dist_pace = sdg_params::SDG_DIST_TEMP;
    float vol_pace = sdg_params::SDG_VOL_TEMP;

    if (_nDist == 0)
    {
        return 0;
    }
    unsigned int i, j;

    double score;
    double *d = new double[_distances.size()];
    double *v = new double[_volumes.size()];
    for (i = 0; i < _distances.size(); ++i)
    {
        d[i] = _distances[i].Measure();
    }

    for (i = 0; i < _volumes.size(); ++i)
    {
        v[i] = _volumes[i].Measure();
    }

    for (int i = 0; i < sdg_params::SDG_NCYCLES; ++i)
    {
        DoCycle(dist_pace, vol_pace);
        score = 0;
        for (j = 0; j < _distances.size(); ++j)
        {
            d[j] = _distances[j].Measure();
            score += _distances[j].Score();
        }
        for (j = 0; j < _volumes.size(); ++j)
        {
            v[j] = _volumes[j].Measure();
        }
        dist_pace -= sdg_params::SDG_ANNEAL_DIST;
        vol_pace -= sdg_params::SDG_ANNEAL_VOL;
    }
    delete[] d;
    delete[] v;
    return score;
}

void sdgEngine::DoCycle(float dist_pace, float vol_pace)
{
    for (int i = 0; i < _nSteps; ++i)
    {
        DoStep(dist_pace, vol_pace);
    }
}

void sdgEngine::DoStep(float dist_pace, float vol_pace)
{
    int x;

    //Decide whether to tweak a volume or a distance
    float r = (float)rand();
    float rm =  RAND_MAX;
    if (_nVol > 0 && r / rm < _vol_odds)
    {
        x = irand_approx(_nVol);
        _volumes[x].Tweak(vol_pace);
    }
    else
    {
        x = irand_approx(_nDist);
        _distances[x].Tweak(dist_pace);
    }
}

void sdgEngine::Explode(float factor)
{
    double com[3];                      //residue center-of-mass
    CenterOfMass(_atoms, com);

    MIAtom_const_iter i, e = _atoms.end();

    for (i = _atoms.begin(); i != e; ++i)
    {
        (*i)->setPosition( (float)(com[0] + factor * ((*i)->x() - com[0])),
                           (float)(com[1] + factor * ((*i)->y() - com[1])),
                           (float)(com[2] + factor * ((*i)->z() - com[2])));
    }
}

sdgDistance::sdgDistance(MIAtom &a1, MIAtom &a2, double lower, double upper, bool isRange)
    : _a1(a1),
      _a2(a2)
{
    _isOpenEnd = false;
    _isRange = isRange;
    _lower = lower;
    _upper = upper;
}

sdgDistance::sdgDistance(MIAtom &a1, MIAtom &a2, double ideal_dist, bool isRange)
    : _a1(a1),
      _a2(a2)
{
    _isOpenEnd = isRange;
    _isRange = isRange;
    _lower = ideal_dist;
}

sdgDistance::sdgDistance(const sdgDistance &orig)
    : _a1(orig._a1),
      _a2(orig._a2)
{
    _isOpenEnd = orig._isOpenEnd;
    _isRange = orig._isRange;
    _lower = orig._lower;
    _upper = orig._upper;
}

sdgDistance&sdgDistance::operator=(const sdgDistance &orig)
{
    _a1.copyShallow(orig._a1);
    _a2.copyShallow(orig._a2);
    _isOpenEnd = orig._isOpenEnd;
    _isRange = orig._isRange;
    _lower = orig._lower;
    _upper = orig._upper;
    return *this;
}

bool sdgDistance::operator==(const sdgDistance &d2) const
{

    return ((&d2._a1 == &_a1 && &d2._a2 == &_a2)
            || (&d2._a1 == &_a2 && &d2._a2 == &_a1));
}

double sdgDistance::Measure()
{
    return AtomDist(_a1, _a2);
}

double sdgDistance::Score()
{
    double sqdist = (_a1.x() - _a2.x()) * (_a1.x() - _a2.x())
                    +(_a1.y() - _a2.y()) * (_a1.y() - _a2.y())
                    +(_a1.z() - _a2.z()) * (_a1.z() - _a2.z());

    if (!_isRange)
    {
        return (sqdist - _lower * _lower) * (sqdist - _lower * _lower) / (_lower * _lower);
    }
    else if (_isOpenEnd && sqdist >= _lower * _lower)
    {
        return 0;
    }
    else if (_isOpenEnd)
    {
        return (sqdist - _lower * _lower) * (sqdist - _lower * _lower) / (_lower * _lower);
    }
    //Finite range: adjust to min or max if we're outside the range
    else if (sqdist < _lower * _lower)
    {
        return (sqdist - _lower * _lower) * (sqdist - _lower * _lower) / (_lower * _lower);
    }
    else if (sqdist > _upper * _upper)
    {
        return (sqdist - _upper * _upper) * (sqdist - _upper * _upper) / (_upper * _upper);
    }
    else
    {
        return 0;
    }
}

void sdgDistance::Tweak(float pace)
{
    double cd = AtomDist(_a1, _a2);              //cd = Current Distance

    //1. Calculate the size of the step for each atom, or return if the
    //distance is in the allowed interval
    double step;
    //Exact distance: just calculate error
    if (!_isRange)
    {
        step = pace * 0.5 * (_lower - cd) / (cd + sdg_params::SDG_NONZERO);
    }
    //Open-ended range: if the distance is too short adjust to min,
    //otherwise we're fine
    else if (_isOpenEnd && cd >= _lower)
    {
        return;
    }
    else if (_isOpenEnd)
    {
        step = pace * 0.5 * (_lower - cd) / (cd + sdg_params::SDG_NONZERO);
    }
    //Finite range: adjust to min or max if we're outside the range
    else if (cd < _lower)
    {
        step = pace * 0.5 * (_lower - cd) / (cd + sdg_params::SDG_NONZERO);
    }
    else if (cd > _upper)
    {
        step = pace * 0.5 * (_upper - cd) / (cd + sdg_params::SDG_NONZERO);
    }
    else
    {
        return;
    }

    //2. Calcluate the vector  direction
    double v[3];
    BondVector(&_a1, &_a2, v);
    if (cd < sdg_params::SDG_NONZERO)                               //Step in an arbitrary direction if atoms
    {
        v[0] = rand() * sdg_params::SDG_NONZERO / RAND_MAX;         //are on top of each other
        v[1] = rand() * sdg_params::SDG_NONZERO / RAND_MAX;
        v[2] = rand() * sdg_params::SDG_NONZERO / RAND_MAX;
    }

    //3. Move the atoms
    AtomStep(&_a1, v, -step);
    AtomStep(&_a2, v, step);

}

double sdgDistance::GetIdeal()
{
    if (!_isRange)
    {
        return _lower;
    }
    else if (_isOpenEnd)
    {
        return _lower;
    }
    else
    {
        return 0.5 * (_lower + _upper);
    }
}

sdgVolume::sdgVolume(MIAtom &center, MIAtom &a1, MIAtom &a2, MIAtom &a3, double lower, double upper, int openEnd)
    : _center(center),
      _a1(a1),
      _a2(a2),
      _a3(a3)
{
    _isRange = true;
    _openEnd = openEnd;
    _lower = lower;
    _upper = upper;
}

sdgVolume::sdgVolume(MIAtom &center, MIAtom &a1, MIAtom &a2, MIAtom &a3, double ideal_volume, int openEnd)
    : _center(center),
      _a1(a1),
      _a2(a2),
      _a3(a3)
{
    if (openEnd == sdg_params::SDG_OPENEND_LOW)
    {
        _isRange = true;
        _openEnd = sdg_params::SDG_OPENEND_LOW;
        _upper = ideal_volume;
    }
    else if (openEnd == sdg_params::SDG_OPENEND_HIGH)
    {
        _isRange = true;
        _openEnd = sdg_params::SDG_OPENEND_HIGH;
        _lower = ideal_volume;
    }
    else if (openEnd == sdg_params::SDG_OPENEND_NONE)
    {
        _isRange = false;
        _openEnd = sdg_params::SDG_OPENEND_NONE;
        _lower = ideal_volume;
    }
}

sdgVolume::sdgVolume(MIAtom &center, MIAtom &a1, MIAtom &a2, MIAtom &a3, double ideal_volume)
    : _center(center),
      _a1(a1),
      _a2(a2),
      _a3(a3)
{
    _isRange = false;
    _openEnd = sdg_params::SDG_OPENEND_NONE;
    _lower = ideal_volume;
}

sdgVolume::sdgVolume(const sdgVolume &orig)
    : _center(orig._center),
      _a1(orig._a1),
      _a2(orig._a2),
      _a3(orig._a3)
{
    _isRange = orig._isRange;
    _openEnd = orig._openEnd;
    _lower = orig._lower;
    _upper = orig._upper;
}

sdgVolume&sdgVolume::operator=(const sdgVolume &orig)
{
    _center.copyShallow(orig._center);
    _a1.copyShallow(orig._a1);
    _a2.copyShallow(orig._a2);
    _a3.copyShallow(orig._a3);
    _isRange = orig._isRange;
    _openEnd = orig._openEnd;
    _lower = orig._lower;
    _upper = orig._upper;
    return *this;
}

bool sdgVolume::operator==(const sdgVolume &vol2) const
{

    if (&vol2._center != &_center
        && &vol2._center != &_a1
        && &vol2._center != &_a2
        && &vol2._center != &_a3)
    {
        return false;
    }
    else if (&vol2._a1 != &_a1
             && &vol2._a1 != &_a2
             && &vol2._a1 != &_a3
             && &vol2._a1 != &_center)
    {
        return false;
    }
    else if (&vol2._a2 != &_a1
             && &vol2._a2 != &_a2
             && &vol2._a2 != &_a3
             && &vol2._a2 != &_center)
    {
        return false;
    }
    else if (&vol2._a3 != &_a1
             && &vol2._a3 != &_a2
             && &vol2._a3 != &_a3
             && &vol2._a3 != &_center)
    {
        return false;
    }
    return true;
}

double sdgVolume::Measure()
{
    return SignedAtomVolume(_center, _a1, _a2, _a3);

    //1. Calculate the bond vectors from the center to each atom,
    //	double b1[3], b2[3], b3[3];				//Bond vectors from center to atoms 1, 2, and 3
    //	BondVector(&_center, &_a1, b1);
    //	BondVector(&_center, &_a2, b2);
    //	BondVector(&_center, &_a3, b3);

    //2. Compute the signed volume
    //	return Cross_and_dot_3D(b1, b2, b3) / 6.0;			//AKA the determinant of the 3x3 matrix[b1,b2,b3]
}

void sdgVolume::Tweak(float pace)
{

    //1. Compute the signed volume
    double cv = SignedAtomVolume(_center, _a1, _a2, _a3);
    double step, vagn;

    //2. Check if the volume is in the allowed interval
    if (!_isRange)
    {
        step = pace * (_lower - cv) / (VolAtomGradNorm(_center, _a1, _a2, _a3) + sdg_params::SDG_NONZERO);
    }
    //Open-ended range: if the distance doesn't meet the cutoff, use
    //the cutoff as the target distance, otherwise we're fine
    else if (_openEnd == sdg_params::SDG_OPENEND_LOW && (cv <= _upper))
    {
        return;
    }
    else if (_openEnd == sdg_params::SDG_OPENEND_LOW)
    {
        step = pace * (_upper - cv) / (VolAtomGradNorm(_center, _a1, _a2, _a3) + sdg_params::SDG_NONZERO);
        //		step = pace * (_upper - cv);
    }
    else if (_openEnd == sdg_params::SDG_OPENEND_HIGH && (cv >= _lower))
    {
        return;
    }
    else if (_openEnd == sdg_params::SDG_OPENEND_HIGH)
    {
        step = pace * (_lower - cv) / (VolAtomGradNorm(_center, _a1, _a2, _a3) + sdg_params::SDG_NONZERO);
        //		step = pace * (_lower - cv);
    }
    //Finite range: adjust to min or max (whichever is closer) if we're outside the range
    else if (cv < _lower)
    {
        vagn = VolAtomGradNorm(_center, _a1, _a2, _a3);
        step = pace * (_lower - cv) / (VolAtomGradNorm(_center, _a1, _a2, _a3) + sdg_params::SDG_NONZERO);
    }
    else if (cv > _upper)
    {
        vagn = VolAtomGradNorm(_center, _a1, _a2, _a3);
        step = pace * (_upper - cv) / (VolAtomGradNorm(_center, _a1, _a2, _a3) + sdg_params::SDG_NONZERO);
    }
    else
    {
        return;
    }

    //3. Calculate the bond vectors from the center to each atom,
    double b1[3], b2[3], b3[3];             //Bond vectors from center to atoms 1, 2, and 3
    BondVector(&_center, &_a1, b1);
    BondVector(&_center, &_a2, b2);
    BondVector(&_center, &_a3, b3);

    //4. Calculate the gradients (of the signed volume) for each of the four atoms
    double grad[4][3];
    //	grad[1][0] = b2[1]*b3[2] - b2[2]*b3[1];
    //	grad[1][1] = b2[2]*b3[0] - b2[0]*b3[2];
    //	grad[1][2] = b2[0]*b3[1] - b2[1]*b3[0];
    //	grad[2][0] = b1[2]*b3[1] - b1[1]*b3[2];
    //	grad[2][1] = b1[0]*b3[2] - b1[2]*b3[0];
    //	grad[2][2] = b1[1]*b3[0] - b1[0]*b3[1];
    //	grad[3][0] = b1[1]*b2[2] - b1[2]*b2[1];
    //	grad[3][1] = b1[2]*b2[0] - b1[0]*b2[2];
    //	grad[3][2] = b1[0]*b2[1] - b1[1]*b2[0];
    //	grad[0][0] = -(grad[1][0] + grad[2][0] + grad[3][0]);
    //	grad[0][1] = -(grad[1][1] + grad[2][1] + grad[3][1]);
    //	grad[0][2] = -(grad[1][2] + grad[2][2] + grad[3][2]);

    MIAtom &a1 = _center;
    MIAtom &a2 = _a1;
    MIAtom &a3 = _a2;
    MIAtom &a4 = _a3;

    grad[0][0] = a2.y()*a3.z() - a2.z()*a3.y() - (a2.y()*a4.z() - a2.z()*a4.y()) + a3.y()*a4.z() - a3.z()*a4.y();
    grad[0][1] = -(a2.x()*a3.z() - a2.z()*a3.x() - (a2.x()*a4.z() - a2.z()*a4.x()) + a3.x()*a4.z() - a3.z()*a4.x());
    grad[0][2] = a2.x()*a3.y() - a2.y()*a3.x() - (a2.x()*a4.y() - a2.y()*a4.x()) + a3.x()*a4.y() - a3.y()*a4.x();
    grad[1][0] = -(a1.y()*a3.z() - a1.z()*a3.y() - (a1.y()*a4.z() - a1.z()*a4.y()) + a3.y()*a4.z() - a3.z()*a4.y());
    grad[1][1] = a1.x()*a3.z() - a1.z()*a3.x() - (a1.x()*a4.z() - a1.z()*a4.x()) + a3.x()*a4.z() - a3.z()*a4.x();
    grad[1][2] = -(a1.x()*a3.y() - a1.y()*a3.x() - (a1.x()*a4.y() - a1.y()*a4.x()) + a3.x()*a4.y() - a3.y()*a4.x());
    grad[2][0] = a1.y()*a2.z() - a1.z()*a2.y() - (a1.y()*a4.z() - a1.z()*a4.y()) + a2.y()*a4.z() - a2.z()*a4.y();
    grad[2][1] = -(a1.x()*a2.z() - a1.z()*a2.x() - (a1.x()*a4.z() - a1.z()*a4.x()) + a2.x()*a4.z() - a2.z()*a4.x());
    grad[2][2] = a1.x()*a2.y() - a1.y()*a2.x() - (a1.x()*a4.y() - a1.y()*a4.x()) + a2.x()*a4.y() - a2.y()*a4.x();
    grad[3][0] = -(a1.y()*a2.z() - a1.z()*a2.y() - (a1.y()*a3.z() - a1.z()*a3.y()) + a2.y()*a3.z() - a2.z()*a3.y());
    grad[3][1] = a1.x()*a2.z() - a1.z()*a2.x() - (a1.x()*a3.z() - a1.z()*a3.x()) + a2.x()*a3.z() - a2.z()*a3.x();
    grad[3][2] = -(a1.x()*a2.y() - a1.y()*a2.x() - (a1.x()*a3.y() - a1.y()*a3.x()) + a2.x()*a3.y() - a2.y()*a3.x());

    //5. Move the atoms
    AtomStep(&_center, grad[0], step / 6.0);
    AtomStep(&_a1, grad[1], step / 6.0);
    AtomStep(&_a2, grad[2], step / 6.0);
    AtomStep(&_a3, grad[3], step / 6.0);

    //  JAC this was unused, so I commented it out
    //  double gn = (grad[1][0] * grad[1][0] +
    //                grad[1][1] * grad[1][1] +
    //                grad[1][2] * grad[1][2] +
    //                grad[2][0] * grad[2][0] +
    //                grad[2][1] * grad[2][1] +
    //                grad[2][2] * grad[2][2] +
    //                grad[3][0] * grad[3][0] +
    //                grad[3][1] * grad[3][1] +
    //                grad[3][2] * grad[3][2] +
    //                grad[0][0] * grad[0][0] +
    //                grad[0][1] * grad[0][1] +
    //                grad[0][2] * grad[0][2]) / 36.0;


    // JAC this is pointless, since the value is not used and there are no side-effects
    // cv = Measure();
}

}
