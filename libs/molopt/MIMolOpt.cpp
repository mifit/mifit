#include <vector>
#include <algorithm>
#include <cfloat>

#include <nongui/nonguilib.h>
#include <math/mathlib.h>
#include <chemlib/chemlib.h>
#include <chemlib/RESIDUE_.h>
#include <chemlib/MIMoleculeBase.h>
#include <conflib/conflib.h>
#include <map/maplib.h>

#include "MIMolOpt.h"


//#include "mifit_algorithm.h"

using namespace chemlib;
using namespace std;

#define X 0
#define Y 1
#define Z 2

#define DEFAULTDIST 1.9F /* connect points closer than this */

// constants for refinement DE algorithm
#define IM1 2147483563
#define IM2 2147483399
#define AM (1.0f/IM1)
#define IA1 40014
#define IA2 40692
#define IQ1 53668
#define IQ2 52774
#define IR1 12211
#define IR2 3791
#define IMM1 (IM1-1)
#define NTAB 32
#define NDIV (1+IMM1/NTAB)
#define EPS 1.2e-7
#define RNMX (1.0f-EPS)


//////////////////////////////////////////////////////////////////////
// Construction/Destruction
//////////////////////////////////////////////////////////////////////

MIMolOpt::MIMolOpt()
    : refineTargetLocked(false)
{

    RefiRes = NULL;
    nucleic = 0;
    AutoFit = false;
    CurrentMap = NULL;
    // Position 0 is empty and is used to mean no save set
    SaveToken = 0;
    BondWeight = AngleWeight = PlaneWeight = MapWeight = BumpWeight = TorsionWeight = 1.0F;
    nCycles = 20;
    fit_while_refine = true;
    RefiVerbose = false;
}

MIMolOpt::~MIMolOpt()
{
}

void MIMolOpt::Refine()
{
    if (!IsRefining())
    {
        return;
    }

    float dt, dp, da, db;

    dt = StdevTorsions();
    dp = StdevPlanes();
    da = StdevAngles();
    db = StdevBonds();
    Logger::log("Before: stdev bonds=%0.3f angles=%0.3f planes=%0.3f torsions=%0.3f",
                db, da, dp, dt);

    for (int i = 0; i < nCycles; i++)
    {
        resetderivatives();
        if (dict.RefiBonds.size() > 0)
        {
            minimize_bonds(dict.RefiBonds, dict.RefiBonds.size());
        }
        if (dict.RefiAngles.size() > 0)
        {
            minimize_angles(dict.RefiAngles, dict.RefiAngles.size());
        }
        if (dict.RefiPlanes.size() > 0)
        {
            minimize_planes(dict.RefiPlanes, dict.RefiPlanes.size());
        }
        if (dict.RefiConstraints.size() > 0)
        {
            minimize_constraints(dict.RefiConstraints, dict.RefiConstraints.size());
        }
        minimize_phipsi();
        minimize_torsions();
        if (dict.RefiBumps.size() > 0)
        {
            minimize_bumps(dict.RefiBumps, dict.RefiBumps.size());
        }
        minimize_map();
        applyderivatives();
    }

    dt = StdevTorsions();
    dp = StdevPlanes();
    da = StdevAngles();
    db = StdevBonds();
    Logger::log("After: stdev bonds=%0.3f angles=%0.3f planes=%0.3f torsions=%0.3f",
                db, da, dp, dt);
}

int MIMolOpt::minimize_bonds(std::vector<Bond> &bonds, unsigned int nbonds)
{
    unsigned int i;
    MIAtom *a1, *a2;
    float r, d, dx, dy, dz;
    float bweight = 1.0, weight;
    float dsum = 0.0, ssum = 0.0;
    float dd, ss;
    int nsum = 0;
    int nmoved = 0;

    if (BondWeight == 0 || nbonds == 0)
    {
        return (0);
    }
    bweight = BondWeight;

    /* find deviations */
    for (i = 0; i < nbonds; i++)
    {
        a1 = bonds[i].getAtom1();
        a2 = bonds[i].getAtom2();
        /* build derivatives */
        r = (float)AtomDist(*a1, *a2);
        if (r <= 0.01)
        {
            continue;           // to avoid divide by 0
        }
        d = r - bonds[i].ideal_length;
        dd = d*d;
        ss = bonds[i].tolerance*bonds[i].tolerance;
        dsum += dd;
        ssum += ss;
        nsum++;
        weight = dd/ss;
        if (weight < 0.10)
        {
            continue;               /* don't waste time on very small shifts */
        }
        if (weight > 9.0)   /* 5.0 in unsquared units */
        {
            weight = 9.0; /* clamp weight */
        }
        weight *= bweight;
        dx = d*(a2->x() - a1->x())/r/2.0f; /*move each atom half the distance*/
        dy = d*(a2->y() - a1->y())/r/2.0f; /*so that the total movement is d */
        dz = d*(a2->z() - a1->z())/r/2.0f;
        dx *= weight;
        dy *= weight;
        dz *= weight;
        a1->addDelta(dx, dy, dz);
        a2->addDelta(-dx, -dy, -dz);
        a1->addWeight(weight);
        a2->addWeight(weight);
        //sumw += weight;
        //sumd += dx*dx + dy*dy + dz*dz;
        nmoved++;
    }
    if (RefiVerbose)
    {
        Logger::log("RMS error in bond lengths: %0.3f (sigma =%0.3f)", sqrt(dsum/(float)nsum), sqrt(ssum/(float)nsum));
    }
    return (nmoved);
}

int MIMolOpt::minimize_angles(std::vector<ANGLE> &angles, unsigned int nangles)
{
    unsigned int i;
    MIAtom *a1, *a2;
    float r, d, dx, dy, dz;
    float aweight = 1.0;
    float weight;
    float dsum = 0.0, ssum = 0.0;
    float dd, ss;
    int nsum = 0;
    int nmoved = 0;

    if (AngleWeight == 0 || nangles == 0)
    {
        return (0);
    }
    aweight = AngleWeight;

    /* find deviations */
    for (i = 0; i < dict.RefiAngles.size(); i++)
    {

        if (angles[i].ideal_angle <= 0.01F)
        {
            // don't minimize unset angles
            continue;
        }
        /* angles are constrained by 1-3 distance rather
         * than by angle value
         */
        a1 = angles[i].getAtom1();
        a2 = angles[i].atom3;
        /* build derivatives */
        r = (float)AtomDist(*a1, *a2);
        if (r <= 0.01F)
        {
            continue;            // to avoid divide by 0
        }
        d = r - angles[i].ideal_angle;
        dd = d*d;
        ss = angles[i].tolerance*angles[i].tolerance;
        dsum += dd;
        ssum += ss;
        nsum++;
        weight = dd/ss;
        if (weight < 0.10F)
        {
            continue;                 /* don't waste time on very small shifts */
        }
        if (weight > 9.0F)
        {
            weight = 9.0F; /* clamp weight */
            if (RefiVerbose)
            {
                Logger::log("Very bad angle: %s %s: %s -> %s: ideal: %6.2f actual: %6.2f",
                            angles[i].res->type().c_str(), angles[i].res->name().c_str(), angles[i].getAtom1()->name(),
                            angles[i].atom3->name(), angles[i].ideal_angle, r);
            }
        }
        weight *= aweight; /* times overall angle weight */
        dx = d*(a2->x() - a1->x())/r/2.0f;
        dy = d*(a2->y() - a1->y())/r/2.0f;
        dz = d*(a2->z() - a1->z())/r/2.0f;
        dx *= weight;
        dy *= weight;
        dz *= weight;
        a1->addDelta(dx, dy, dz);
        a2->addDelta(-dx, -dy, -dz);
        a1->addWeight(weight);
        a2->addWeight(weight);
        //sumw += weight;
        //sumd += dx*dx + dy*dy + dz*dz;
        nmoved++;
    }
    if (RefiVerbose)
    {
        Logger::log("RMS error in angle lengths: %0.3f (sigma =%0.3f)",
                    sqrt(dsum/nsum), sqrt(ssum/nsum));
    }
    //Logger::log("angles moved %f weight %f", (float)sqrt(sumd/(double)dict.RefiAngles.size()), sumw);
    return (nmoved);
}

/* fit a least squares plane to a set of points
 * method from schomaker et al, acta cryst, 12, p. 600, 1959.
 */
void lsqplane(PLANE &plane)
{
    float xs[3], xxs[3][3], b[3][3];
    float a[3][3], vmi[3], bv[3], vm[3], d;
    float zip = 1.0E-5F;
    float orm, vm0, ratio0, ratio1, ratio2, rat01, rat02;
    int i, j, kk, nnn, k, n;

    n = plane.natoms;
    for (i = 0; i < 3; i++)
    {
        xs[i] = 0.0;
    }
    for (k = 0; k < n; k++)
    {
        /* zip is added to prevent numerical instability if
         * atoms in a plane = 0
         */
        xs[0] += plane.atoms[k]->x()+zip;
        xs[1] += plane.atoms[k]->y()+zip;
        xs[2] += plane.atoms[k]->z()+zip;
    }
    for (i = 0; i < 3; i++)
    {
        for (j = 0; j < 3; j++)
        {
            xxs[i][j] = 0.0;
        }
    }

    for (k = 0; k < n; k++)
    {
        xxs[0][0] += plane.atoms[k]->x() * plane.atoms[k]->x();
        xxs[0][1] += plane.atoms[k]->x() * plane.atoms[k]->y();
        xxs[0][2] += plane.atoms[k]->x() * plane.atoms[k]->z();
        xxs[1][0] += plane.atoms[k]->y() * plane.atoms[k]->x();
        xxs[1][1] += plane.atoms[k]->y() * plane.atoms[k]->y();
        xxs[1][2] += plane.atoms[k]->y() * plane.atoms[k]->z();
        xxs[2][0] += plane.atoms[k]->z() * plane.atoms[k]->x();
        xxs[2][1] += plane.atoms[k]->z() * plane.atoms[k]->y();
        xxs[2][2] += plane.atoms[k]->z() * plane.atoms[k]->z();
    }
    for (i = 0; i < 3; i++)
    {
        for (j = 0; j < 3; j++)
        {
            a[i][j] = xxs[i][j] - xs[i]*xs[j]/(float)n;
        }
    }

    /* evaluate matrix */
    b[0][0] = a[1][1] * a[2][2] - a[1][2] * a[2][1];
    b[1][0] = a[2][0] * a[1][2] - a[1][0] * a[2][2];
    b[2][0] = a[1][0] * a[2][1] - a[2][0] * a[1][1];
    b[0][1] = a[2][1] * a[0][2] - a[0][1] * a[2][2];
    b[1][1] = a[0][0] * a[2][2] - a[2][0] * a[0][2];
    b[2][1] = a[2][0] * a[0][1] - a[0][0] * a[2][1];
    b[0][2] = a[0][1] * a[1][2] - a[0][2] * a[1][1];
    b[1][2] = a[1][0] * a[0][2] - a[0][0] * a[1][2];
    b[2][2] = a[0][0] * a[1][1] - a[1][0] * a[0][1];

    /* choose the largest column vector of b as initial solution */
    bv[0] = b[0][0]*b[0][0] + b[1][0]*b[1][0] + b[2][0]*b[2][0];
    bv[1] = b[0][1]*b[0][1] + b[1][1]*b[1][1] + b[2][1]*b[2][1];
    bv[2] = b[0][0]*b[0][2] + b[1][2]*b[1][2] + b[2][2]*b[2][2];
    kk = 0;
    if (bv[1] > bv[0])
    {
        kk = 1;
    }
    if (bv[2] > bv[kk])
    {
        kk = 2;
    }
    vm0 = b[0][kk];
    for (i = 0; i < 3; i++)
    {
        vmi[i] = b[i][kk]/vm0;
    }
    /* solve to convergence by iteration of M(I)=B*M(I-1) */
    for (nnn = 0; nnn < 10; nnn++)
    {
        vm[0] = b[0][0]*vmi[0]+b[0][1]*vmi[1]+b[0][2]*vmi[2];
        vm[1] = b[1][0]*vmi[0]+b[1][1]*vmi[1]+b[1][2]*vmi[2];
        vm[2] = b[2][0]*vmi[0]+b[2][1]*vmi[1]+b[2][2]*vmi[2];
        ratio0 = vm[0]/vmi[0];
        ratio1 = vm[1]/vmi[1];
        ratio2 = vm[2]/vmi[2];
        rat01 = (float)fabs(ratio1/ratio0-1.0f);
        rat02 = (float)fabs(ratio2/ratio0-1.0f);
        if (rat01 < zip && rat02 < zip)
        {
            break;
        }
        else
        {
            for (i = 0; i < 3; i++)
            {
                vmi[i] = vm[i]/vm[0];
            }
        }
    } /* 100*/
    orm = 0.0;
    /* normalize the solution vectors */
    for (i = 0; i < 3; i++)
    {
        orm += vm[i]*vm[i];
    }
    orm = sqrt(orm);
    for (i = 0; i < 3; i++)
    {
        vm[i] /= orm;
    }
    d = (vm[0]*xs[0]+vm[1]*xs[1]+vm[2]*xs[2])/(float)n;
    plane.vm[0] = vm[0];
    plane.vm[1] = vm[1];
    plane.vm[2] = vm[2];
    plane.d = d;
}

int MIMolOpt::minimize_planes(std::vector<PLANE> &planes, unsigned int nplanes)
{
    unsigned int i, j;
    MIAtom *a1;
    float d, dx, dy, dz;
    float bweight = 1.0, weight;
    float dsum = 0.0, ssum = 0.0;
    float dd, ss;
    int nsum = 0;

    if (PlaneWeight == 0 || nplanes == 0)
    {
        return (0);
    }
    bweight = PlaneWeight;

    /* find deviations */
    for (i = 0; i < nplanes; i++)
    {
        /* build derivatives */
        lsqplane(planes[i]);
        for (j = 0; j < (unsigned int)planes[i].natoms; j++)
        {
            a1 = planes[i].atoms[j];
            d =  a1->x() * planes[i].vm[0];
            d += a1->y() * planes[i].vm[1];
            d += a1->z() * planes[i].vm[2];
            d -= planes[i].d;
            dd = d*d;
            ss = planes[i].tolerance*planes[i].tolerance;
            dsum += dd;
            ssum += ss;
            nsum++;
            weight = dd/ss;
            if (weight < 0.05F)
            {
                /* don't waste time on very small shifts */
                continue;
            }
            if (weight > 9.0F) /* 3.0 in unsquared units */
            {
                if (RefiVerbose)
                {
                    if (planes[i].res)
                        Logger::log("Very bad out of plane: %s %s: %s: %6.2f",
                                    planes[i].res->type().c_str(), planes[i].res->name().c_str(), planes[i].atoms[j]->name(), d);
                    else
                        Logger::log("Very bad out of plane in plane %d: %s: %6.2f", i,
                                    planes[i].atoms[j]->name(), d);

                }
                weight = 9.0F; /* clamp weight */
            }
            weight *= bweight;
            dx = d * planes[i].vm[0];
            dy = d * planes[i].vm[1];
            dz = d * planes[i].vm[2];
            dx *= weight;
            dy *= weight;
            dz *= weight;
            a1->addDelta(-dx, -dy, -dz);
            a1->addWeight(weight);
        }
    }
    if (RefiVerbose)
    {
        Logger::log("RMS error in planes: %0.3f (sigma =%0.3f)",
                    sqrt(dsum/nsum), sqrt(ssum/nsum));
    }
    return (nplanes);
}

int MIMolOpt::minimize_constraints(std::vector<Bond> &bonds, unsigned int nbonds)
{
    /* like minimize bonds except for handling of frozen atoms */
    unsigned int i;
    MIAtom *a1, *a2;
    extern float distance();
    float r, d, dx, dy, dz;
    float bweight = 1.0, weight;
    float dsum = 0.0, ssum = 0.0;
    float dd, ss;
    int nsum = 0;
    if (BondWeight == 0 || nbonds == 0)
    {
        return (0);
    }
    bweight = BondWeight;

    /* find deviations */
    for (i = 0; i < nbonds; i++)
    {
        a1 = bonds[i].getAtom1();
        a2 = bonds[i].getAtom2();
        /* build derivatives */
        r = (float)AtomDist(*a1, *a2);
        if (r < 0.01)
        {
            continue;          // to avoid divide by zero
        }
        d = r - bonds[i].ideal_length;
        dd = d*d;
        ss = bonds[i].tolerance*bonds[i].tolerance;
        dsum += dd;
        ssum += ss;
        nsum++;
        weight = dd/ss;
        if (weight < 0.10)
        {
            continue;               /* don't waste time on very small shifts */
        }
        if ((weight = 9.0) > 9.0) /* 5.0 in unsquared units */
        {
            weight = 9.0; /* clamp weight */
        }
        weight *= bweight;
        if ( (a2->type()) & AtomType::FREEZE)
        {
            dx = d*(a2->x() - a1->x())/r;
            dy = d*(a2->y() - a1->y())/r;
            dz = d*(a2->z() - a1->z())/r;
            dx *= weight;
            dy *= weight;
            dz *= weight;
            a1->addDelta(dx, dy, dz);
            a1->addWeight(weight);
        }
        else
        {
            dx = d*(a2->x() - a1->x())/r/2.0f; /*move each atom half the distance*/
            dy = d*(a2->y() - a1->y())/r/2.0f; /*so that the total movement is d */
            dz = d*(a2->z() - a1->z())/r/2.0f;
            dx *= weight;
            dy *= weight;
            dz *= weight;
            a1->addDelta(dx, dy, dz);
            a2->addDelta(-dx, -dy, -dz);
            a1->addWeight(weight);
            a2->addWeight(weight);
        }
    }
    if (RefiVerbose)
    {
        Logger::log("RMS error in constraints: %0.3f (sigma =%0.3f)",
                    sqrt(dsum/nsum), sqrt(ssum/nsum));
    }
    return (nbonds);
}

int MIMolOpt::minimize_bumps(std::vector<Bond> &bumps, unsigned int nbumps)
{
    unsigned int i;
    MIAtom *a1, *a2;
    float r, d, dx, dy, dz;
    float bweight = 1.0, weight;
    float dsum = 0.0, ssum = 0.0;
    float dd, ss;
    int nsum = 0;
    //int color = Colors::WHITE;
    //float mx,my,mz;
    if (BumpWeight <= 0 || nbumps == 0)
    {
        return (0);
    }
    bweight = BumpWeight;

    //Vu.clear();
    /* find deviations */
    for (i = 0; i < nbumps; i++)
    {
        a1 = bumps[i].getAtom1();
        a2 = bumps[i].getAtom2();
        /* build derivatives */
        r = (float)AtomDist(*a1, *a2);
        if (r < 0.01F)
        {
            r = 0.01F;           // to avoid the unlikely event of a divide by 0
        }
        d = r - bumps[i].ideal_length;
        if (d >= (-0.05))
        {
            continue;
        }
        dd = d*d;
        ss = bumps[i].tolerance*bumps[i].tolerance;
        weight = dd/ss;
        //if(d < -2.0F*bumps[i].tolerance){
        //color = RED;
        //}else{
        //color = WHITE;
        //}
        //mx = (a2->x+a1->x)/2.0F;
        //my = (a2->y+a1->y)/2.0F;
        //mz = (a2->z+a1->z)/2.0F;
        dsum += dd;
        ssum += ss;
        nsum++;
        /* don't waste time on very small shifts */
        if (weight < 0.10F)
        {
            continue;
        }
        // clamp large shifts
        if (weight > 9.0F)
        {
            weight = 9.0F;
        }
        weight *= bweight;
        dx = d*(a2->x() - a1->x())/r/2.0F; /*move each atom half the distance*/
        dy = d*(a2->y() - a1->y())/r/2.0F; /*so that the total movement is d */
        dz = d*(a2->z() - a1->z())/r/2.0F;
        dx *= weight;
        dy *= weight;
        dz *= weight;
        a1->addDelta(dx, dy, dz);
        a2->addDelta(-dx, -dy, -dz);
        a1->addWeight(weight);
        a2->addWeight(weight);
    }
    return (nbumps);
}

//FIXME: this was called just dTorsion, but it conflicted with the version in chemlib.
//       possibly they do the same thing, now that ANATOM and MIAtom have been merged
/* return dx dy dz needed to move a3 about a1 a2 bond by dangle */
static int my_dTorsion(MIAtom *a1, MIAtom *a2, MIAtom *a3, float dangle, float *dx, float *dy, float *dz)
{
    float tmat[4][3];
    initrotate(a1->x(), a1->y(), a1->z(), a2->x() -a1->x(),
               a2->y() - a1->y(), a2->z() - a1->z(),
               dangle, tmat);
    xl_rotate(a3->x(), a3->y(), a3->z(), dx, dy, dz, tmat);
    *dx = *dx - a3->x();
    *dy = *dy - a3->y();
    *dz = *dz - a3->z();
    return 1;
}

int MIMolOpt::minimize_torsions()
{
    //TORSION pp[];
    //int n;
    float w;
    float dx, dy, dz;
    float chi, dchi, d;
    //int sweight = 10;
    float bweight;
    int i, j;
    int n = dict.RefiTorsions.size();
    float ideal;
    char type[10];
    bweight = TorsionWeight;
    for (i = 0; i < n; i++)
    {
        chi = (float)CalcAtomTorsion(dict.RefiTorsions[i].getAtom1(), dict.RefiTorsions[i].getAtom2(), dict.RefiTorsions[i].atom3, dict.RefiTorsions[i].atom4);
        if (chi < 0.0)
        {
            chi += 360.0;
        }
        dchi = 0.0;
        for (j = 0; j < dict.RefiTorsions[i].nideal; j++)
        {
            d = dict.RefiTorsions[i].ideal[j] - chi;
            if (d < -180.0)
            {
                d += 360.0;
            }
            if (d >  180.0)
            {
                d -= 360.0;
            }
            if (fabs(d) < fabs(dchi) || j == 0)
            {
                dchi = d;
                ideal = dict.RefiTorsions[i].ideal[j];
            }
        }

        if (fabs(dchi) > 0.1F)
        {
            /* clamp to avoid excessive movements in a single cycle
             * and to preserve small angle approximation of dx,dy,dz */
            w  = dchi*dchi/(dict.sigmatorsion*dict.sigmatorsion)*bweight;
            if (dchi > 5.0)
            {
                dchi = 5.0;
            }
            if (dchi < -5.0)
            {
                dchi = -5.0;
            }
            strcpy(type, dict.RefiTorsions[i].type);
            my_dTorsion(dict.RefiTorsions[i].getAtom2(), dict.RefiTorsions[i].atom3, dict.RefiTorsions[i].atom4, dchi, &dx, &dy, &dz);
            dx *= w;
            dy *= w;
            dz *= w;
            dict.RefiTorsions[i].atom4->addDelta(dx, dy, dz);
            dict.RefiTorsions[i].atom4->addWeight(w);
        }
    }
    return n;
}

float MIMolOpt::phipsi_energy(float phi, float psi)
{
    /* allowed phi-psi region for all residues except Gly and Pro */
    /* table is every 10 degrees from -180 to 180 with phi fastest */
    static float phipsidat[] /*1296*/  =
    {
        3.48F, 3.74F, 4.03F, 4.10F, 4.18F, 4.24F, 4.26F, 4.24F, 4.18F, 4.10F, 3.91F, 3.62F, 3.24F, 2.88F, 2.40F, 1.92F, 1.44F, 1.20F, 1.08F, 1.08F, 1.32F, 1.56F, 1.80F, 1.80F, 1.80F, 1.68F, 1.44F, 1.20F, 1.08F, 1.08F, 1.32F, 1.56F, 1.92F, 2.16F, 2.64F, 3.00F,
        3.60F, 3.98F, 4.39F, 4.46F, 4.54F, 4.60F, 4.62F, 4.60F, 4.54F, 4.46F, 4.27F, 3.86F, 3.48F, 3.12F, 2.76F, 2.28F, 1.80F, 1.56F, 1.32F, 1.20F, 1.20F, 1.32F, 1.44F, 1.44F, 1.44F, 1.32F, 1.20F, 1.08F, 1.08F, 1.08F, 1.44F, 1.80F, 2.28F, 2.52F, 2.88F, 3.12F,
        3.48F, 3.86F, 4.27F, 4.46F, 4.54F, 4.60F, 4.62F, 4.60F, 4.54F, 4.46F, 4.39F, 4.10F, 3.72F, 3.36F, 3.00F, 2.76F, 2.28F, 2.04F, 1.56F, 1.32F, 1.08F, 1.08F, 1.08F, 1.08F, 1.08F, 1.08F, 1.08F, 1.08F, 1.08F, 1.08F, 1.44F, 1.80F, 2.28F, 2.64F, 3.00F, 3.24F,
        3.36F, 3.72F, 4.08F, 4.32F, 4.32F, 4.32F, 4.32F, 4.32F, 4.32F, 4.32F, 4.32F, 4.20F, 3.96F, 3.60F, 3.24F, 2.88F, 2.52F, 2.28F, 1.80F, 1.44F, 1.08F, 1.08F, 1.08F, 1.08F, 1.08F, 1.08F, 1.08F, 1.08F, 1.08F, 1.08F, 1.44F, 1.80F, 2.28F, 2.52F, 2.88F, 3.12F,
        3.24F, 3.60F, 3.96F, 4.32F, 4.32F, 4.32F, 4.32F, 4.32F, 4.32F, 4.32F, 4.32F, 4.32F, 4.08F, 3.72F, 3.36F, 3.00F, 2.64F, 2.28F, 1.80F, 1.44F, 1.08F, 1.08F, 1.08F, 1.08F, 1.08F, 1.08F, 1.08F, 1.08F, 1.08F, 1.08F, 1.32F, 1.68F, 2.04F, 2.28F, 2.64F, 3.00F,
        3.24F, 3.48F, 3.84F, 4.20F, 4.32F, 4.32F, 4.32F, 4.32F, 4.32F, 4.32F, 4.32F, 4.20F, 3.96F, 3.60F, 3.36F, 2.88F, 2.52F, 2.16F, 1.80F, 1.44F, 1.08F, 1.20F, 1.32F, 1.44F, 1.44F, 1.44F, 1.44F, 1.32F, 1.20F, 1.08F, 1.20F, 1.44F, 1.80F, 2.04F, 2.52F, 2.88F,
        3.24F, 3.36F, 3.60F, 3.96F, 4.20F, 4.32F, 4.32F, 4.32F, 4.32F, 4.32F, 4.20F, 3.96F, 3.60F, 3.36F, 3.24F, 3.00F, 2.64F, 2.28F, 1.92F, 1.68F, 1.44F, 1.56F, 1.68F, 1.80F, 1.80F, 1.80F, 1.80F, 1.56F, 1.32F, 1.08F, 1.08F, 1.20F, 1.56F, 1.92F, 2.52F, 2.88F,
        3.24F, 3.24F, 3.36F, 3.60F, 3.96F, 4.08F, 4.08F, 4.08F, 4.08F, 4.08F, 3.84F, 3.60F, 3.36F, 3.24F, 3.24F, 3.12F, 2.76F, 2.40F, 2.04F, 1.92F, 1.80F, 1.92F, 2.04F, 2.28F, 2.40F, 2.40F, 2.28F, 1.80F, 1.44F, 1.08F, 1.08F, 1.08F, 1.44F, 1.80F, 2.52F, 2.88F,
        3.12F, 3.24F, 3.24F, 3.36F, 3.60F, 3.72F, 3.72F, 3.72F, 3.72F, 3.72F, 3.48F, 3.36F, 3.24F, 3.12F, 3.00F, 2.88F, 2.64F, 2.40F, 2.16F, 2.16F, 2.28F, 2.40F, 2.52F, 2.64F, 2.64F, 2.52F, 2.28F, 1.80F, 1.44F, 1.08F, 1.08F, 1.08F, 1.44F, 1.80F, 2.40F, 2.76F,
        3.00F, 3.24F, 3.24F, 3.24F, 3.36F, 3.36F, 3.36F, 3.36F, 3.36F, 3.36F, 3.24F, 3.12F, 3.00F, 2.76F, 2.64F, 2.52F, 2.40F, 2.28F, 2.16F, 2.16F, 2.40F, 2.64F, 2.88F, 3.00F, 2.88F, 2.64F, 2.28F, 1.80F, 1.44F, 1.08F, 1.08F, 1.08F, 1.44F, 1.80F, 2.28F, 2.64F,
        2.76F, 3.12F, 3.24F, 3.24F, 3.24F, 3.24F, 3.24F, 3.24F, 3.24F, 3.24F, 3.24F, 3.12F, 2.88F, 2.52F, 2.16F, 1.92F, 1.80F, 1.80F, 1.92F, 2.04F, 2.40F, 2.76F, 3.12F, 3.24F, 2.88F, 2.52F, 2.16F, 1.92F, 1.56F, 1.20F, 1.08F, 1.08F, 1.44F, 1.80F, 2.16F, 2.40F,
        2.52F, 2.88F, 3.12F, 3.24F, 3.24F, 3.24F, 3.24F, 3.24F, 3.24F, 3.24F, 3.24F, 3.12F, 2.76F, 2.40F, 1.92F, 1.68F, 1.44F, 1.44F, 1.68F, 1.92F, 2.40F, 2.76F, 3.12F, 3.24F, 3.00F, 2.76F, 2.40F, 2.16F, 1.68F, 1.32F, 1.08F, 1.08F, 1.44F, 1.80F, 2.16F, 2.28F,
        2.28F, 2.64F, 3.00F, 3.24F, 3.24F, 3.24F, 3.24F, 3.24F, 3.24F, 3.24F, 3.24F, 3.12F, 2.76F, 2.40F, 1.80F, 1.44F, 1.08F, 1.08F, 1.44F, 1.80F, 2.40F, 2.76F, 3.24F, 3.36F, 3.24F, 2.88F, 2.52F, 2.28F, 1.80F, 1.44F, 1.08F, 1.08F, 1.44F, 1.80F, 2.16F, 2.16F,
        2.28F, 2.64F, 3.00F, 3.24F, 3.24F, 3.24F, 3.24F, 3.24F, 3.24F, 3.24F, 3.24F, 3.00F, 2.64F, 2.28F, 1.80F, 1.44F, 1.08F, 1.08F, 1.44F, 1.80F, 2.40F, 2.76F, 3.36F, 3.60F, 3.60F, 3.12F, 2.64F, 2.28F, 1.80F, 1.44F, 1.08F, 1.08F, 1.44F, 1.80F, 2.16F, 2.16F,
        2.40F, 2.76F, 3.12F, 3.24F, 3.24F, 3.24F, 3.36F, 3.36F, 3.36F, 3.24F, 3.24F, 3.00F, 2.64F, 2.28F, 1.92F, 1.56F, 1.20F, 1.08F, 1.44F, 1.80F, 2.28F, 2.64F, 3.24F, 3.60F, 3.60F, 3.12F, 2.64F, 2.28F, 1.80F, 1.44F, 1.08F, 1.08F, 1.44F, 1.80F, 2.16F, 2.16F,
        2.52F, 2.88F, 3.24F, 3.24F, 3.24F, 3.36F, 3.60F, 3.72F, 3.72F, 3.48F, 3.36F, 3.12F, 2.76F, 2.40F, 2.04F, 1.80F, 1.44F, 1.20F, 1.32F, 1.68F, 2.04F, 2.40F, 2.88F, 3.36F, 3.48F, 3.12F, 2.64F, 2.28F, 1.80F, 1.44F, 1.08F, 1.08F, 1.44F, 1.80F, 2.16F, 2.16F,
        2.52F, 2.88F, 3.24F, 3.24F, 3.24F, 3.48F, 3.84F, 4.08F, 4.08F, 3.84F, 3.60F, 3.36F, 3.00F, 2.64F, 2.28F, 2.04F, 1.68F, 1.32F, 1.20F, 1.56F, 1.92F, 2.28F, 2.52F, 2.88F, 3.12F, 3.12F, 2.76F, 2.40F, 1.80F, 1.44F, 1.08F, 1.08F, 1.44F, 1.80F, 2.16F, 2.16F,
        2.52F, 2.88F, 3.24F, 3.24F, 3.36F, 3.72F, 4.08F, 4.32F, 4.32F, 4.20F, 3.96F, 3.60F, 3.24F, 2.88F, 2.52F, 2.28F, 1.92F, 1.56F, 1.20F, 1.44F, 1.80F, 2.16F, 2.40F, 2.64F, 3.00F, 3.00F, 2.76F, 2.40F, 1.80F, 1.44F, 1.08F, 1.08F, 1.44F, 1.80F, 2.16F, 2.16F,
        2.52F, 2.88F, 3.24F, 3.24F, 3.36F, 3.72F, 4.08F, 4.32F, 4.32F, 4.32F, 4.20F, 3.96F, 3.60F, 3.24F, 2.76F, 2.40F, 2.04F, 1.80F, 1.44F, 1.56F, 1.80F, 2.16F, 2.28F, 2.40F, 2.64F, 2.76F, 2.64F, 2.40F, 1.80F, 1.44F, 1.08F, 1.08F, 1.44F, 1.80F, 2.16F, 2.16F,
        2.40F, 2.76F, 3.12F, 3.24F, 3.36F, 3.72F, 4.08F, 4.32F, 4.32F, 4.32F, 4.32F, 4.20F, 3.96F, 3.60F, 3.12F, 2.64F, 2.28F, 2.04F, 1.68F, 1.68F, 1.80F, 2.16F, 2.28F, 2.28F, 2.40F, 2.40F, 2.40F, 2.28F, 1.80F, 1.44F, 1.08F, 1.08F, 1.44F, 1.80F, 2.16F, 2.16F,
        2.28F, 2.64F, 3.00F, 3.24F, 3.24F, 3.48F, 3.84F, 4.20F, 4.32F, 4.32F, 4.32F, 4.32F, 4.20F, 3.84F, 3.36F, 2.88F, 2.52F, 2.28F, 1.92F, 1.80F, 1.68F, 1.80F, 1.80F, 1.80F, 1.80F, 1.80F, 1.80F, 1.80F, 1.56F, 1.32F, 1.08F, 1.08F, 1.44F, 1.80F, 2.16F, 2.16F,
        2.16F, 2.40F, 2.76F, 3.12F, 3.24F, 3.36F, 3.60F, 3.96F, 4.20F, 4.32F, 4.32F, 4.32F, 4.32F, 4.08F, 3.72F, 3.24F, 2.76F, 2.40F, 2.04F, 1.80F, 1.56F, 1.44F, 1.44F, 1.44F, 1.44F, 1.44F, 1.44F, 1.44F, 1.32F, 1.20F, 1.08F, 1.08F, 1.44F, 1.80F, 2.16F, 2.16F,
        2.28F, 2.52F, 2.88F, 3.12F, 3.24F, 3.24F, 3.36F, 3.60F, 3.96F, 4.20F, 4.32F, 4.32F, 4.32F, 4.20F, 3.84F, 3.48F, 3.00F, 2.64F, 2.28F, 1.92F, 1.56F, 1.20F, 1.08F, 1.08F, 1.08F, 1.08F, 1.08F, 1.08F, 1.08F, 1.08F, 1.08F, 1.08F, 1.44F, 1.80F, 2.16F, 2.16F,
        2.52F, 2.64F, 2.88F, 3.00F, 3.24F, 3.24F, 3.24F, 3.36F, 3.60F, 3.84F, 4.08F, 4.20F, 4.32F, 4.32F, 3.96F, 3.60F, 3.12F, 2.76F, 2.40F, 2.04F, 1.68F, 1.32F, 1.08F, 1.08F, 1.08F, 1.08F, 1.08F, 1.08F, 1.08F, 1.08F, 1.08F, 1.08F, 1.44F, 1.80F, 2.16F, 2.28F,
        2.52F, 2.64F, 2.88F, 3.00F, 3.24F, 3.24F, 3.24F, 3.24F, 3.36F, 3.48F, 3.72F, 3.84F, 3.96F, 3.96F, 3.72F, 3.48F, 3.24F, 3.00F, 2.64F, 2.28F, 1.80F, 1.44F, 1.08F, 1.08F, 1.08F, 1.08F, 1.08F, 1.08F, 1.08F, 1.08F, 1.08F, 1.08F, 1.44F, 1.80F, 2.16F, 2.28F,
        2.40F, 2.40F, 2.52F, 2.64F, 3.00F, 3.12F, 3.24F, 3.24F, 3.24F, 3.24F, 3.36F, 3.48F, 3.60F, 3.60F, 3.48F, 3.36F, 3.24F, 3.00F, 2.64F, 2.28F, 1.80F, 1.44F, 1.08F, 1.08F, 1.08F, 1.08F, 1.08F, 1.08F, 1.08F, 1.08F, 1.08F, 1.08F, 1.44F, 1.80F, 2.16F, 2.28F,
        1.92F, 2.04F, 2.28F, 2.40F, 2.76F, 3.00F, 3.24F, 3.24F, 3.24F, 3.24F, 3.24F, 3.24F, 3.24F, 3.24F, 3.12F, 3.00F, 2.88F, 2.76F, 2.52F, 2.28F, 1.80F, 1.44F, 1.08F, 1.08F, 1.08F, 1.08F, 1.08F, 1.08F, 1.08F, 1.08F, 1.08F, 1.08F, 1.32F, 1.56F, 1.80F, 1.80F,
        1.56F, 1.68F, 1.92F, 2.04F, 2.40F, 2.64F, 2.88F, 2.88F, 2.88F, 2.88F, 3.00F, 3.00F, 3.00F, 2.88F, 2.76F, 2.64F, 2.52F, 2.40F, 2.28F, 2.04F, 1.68F, 1.32F, 1.08F, 1.08F, 1.08F, 1.08F, 1.08F, 1.08F, 1.08F, 1.08F, 1.08F, 1.08F, 1.20F, 1.32F, 1.44F, 1.44F,
        1.20F, 1.32F, 1.68F, 1.92F, 2.28F, 2.40F, 2.64F, 2.76F, 2.76F, 2.64F, 2.64F, 2.64F, 2.64F, 2.52F, 2.40F, 2.28F, 2.04F, 1.92F, 1.80F, 1.68F, 1.44F, 1.20F, 1.08F, 1.08F, 1.08F, 1.08F, 1.08F, 1.08F, 1.08F, 1.08F, 1.08F, 1.08F, 1.08F, 1.08F, 1.08F, 1.08F,
        1.32F, 1.44F, 1.68F, 1.92F, 2.28F, 2.40F, 2.52F, 2.52F, 2.40F, 2.28F, 2.28F, 2.28F, 2.28F, 2.04F, 1.92F, 1.80F, 1.68F, 1.56F, 1.44F, 1.32F, 1.32F, 1.32F, 1.44F, 1.44F, 1.44F, 1.32F, 1.20F, 1.08F, 1.08F, 1.08F, 1.08F, 1.08F, 1.08F, 1.08F, 1.08F, 1.20F,
        1.56F, 1.80F, 1.92F, 2.04F, 2.28F, 2.52F, 2.76F, 2.88F, 2.76F, 2.52F, 2.28F, 2.16F, 2.16F, 1.92F, 1.68F, 1.44F, 1.32F, 1.20F, 1.08F, 1.08F, 1.32F, 1.56F, 1.80F, 1.80F, 1.80F, 1.56F, 1.32F, 1.08F, 1.08F, 1.08F, 1.08F, 1.08F, 1.08F, 1.08F, 1.08F, 1.32F,
        1.92F, 2.40F, 2.52F, 2.40F, 2.52F, 2.76F, 3.00F, 3.00F, 2.88F, 2.64F, 2.40F, 2.16F, 2.16F, 1.92F, 1.56F, 1.20F, 1.08F, 1.08F, 1.08F, 1.08F, 1.44F, 1.80F, 2.28F, 2.28F, 2.28F, 1.92F, 1.56F, 1.20F, 1.08F, 1.08F, 1.08F, 1.08F, 1.08F, 1.08F, 1.08F, 1.44F,
        2.16F, 2.64F, 2.88F, 2.76F, 2.76F, 2.88F, 3.12F, 3.24F, 3.24F, 3.00F, 2.64F, 2.28F, 2.16F, 2.04F, 1.68F, 1.32F, 1.08F, 1.08F, 1.08F, 1.08F, 1.44F, 1.80F, 2.40F, 2.40F, 2.40F, 2.04F, 1.68F, 1.32F, 1.08F, 1.08F, 1.08F, 1.08F, 1.20F, 1.32F, 1.44F, 1.68F,
        2.40F, 2.88F, 3.24F, 3.12F, 3.12F, 3.12F, 3.24F, 3.24F, 3.24F, 3.12F, 2.88F, 2.64F, 2.40F, 2.28F, 1.80F, 1.44F, 1.08F, 1.08F, 1.08F, 1.08F, 1.44F, 1.80F, 2.52F, 2.64F, 2.64F, 2.28F, 1.80F, 1.44F, 1.08F, 1.08F, 1.08F, 1.08F, 1.32F, 1.56F, 1.80F, 1.92F,
        2.76F, 3.00F, 3.24F, 3.24F, 3.24F, 3.24F, 3.24F, 3.24F, 3.24F, 3.24F, 3.12F, 3.00F, 2.64F, 2.40F, 1.80F, 1.44F, 1.08F, 1.08F, 1.08F, 1.08F, 1.44F, 1.80F, 2.40F, 2.52F, 2.52F, 2.28F, 1.80F, 1.44F, 1.08F, 1.08F, 1.08F, 1.08F, 1.44F, 1.80F, 2.28F, 2.40F,
        3.12F, 3.36F, 3.60F, 3.60F, 3.60F, 3.60F, 3.60F, 3.60F, 3.60F, 3.60F, 3.48F, 3.36F, 3.00F, 2.64F, 2.04F, 1.56F, 1.20F, 1.08F, 1.08F, 1.08F, 1.44F, 1.80F, 2.28F, 2.40F, 2.40F, 2.16F, 1.68F, 1.32F, 1.08F, 1.08F, 1.20F, 1.32F, 1.68F, 1.92F, 2.40F, 2.64F
    };
    int index;

    /* table start at 180 */
    phi += 180.0F;
    /* table is flipped for psi */
    psi = -psi;
    psi += 180.0F;
    if (psi >= 360.0F)
    {
        psi -= 360.0F;
    }
    if (phi >= 360.0F)
    {
        phi -= 360.0F;
    }
    if (psi < 0.0)
    {
        psi += 360.0F;
    }
    if (phi < 0.0)
    {
        phi += 360.0F;
    }
    if (phi > 360.0F)
    {
        return 0.0;
    }
    if (psi > 360.0F)
    {
        return 0.0;
    }
    index = (int)((psi+5.0F)/10.0F)*36 + (int)((phi+5.0F)/10.0F);
    if (index < 0)
    {
        index = 0;
    }
    if (index >= 1296)
    {
        index = 1295;
    }
    return phipsidat[index];
}

int MIMolOpt::minimize_phipsi()
{
    /* get energy at current position and +/- 10 degrees */
    //TORSION pp[];
    int n = dict.RefiPhiPsis.size();
    int i, ii, jj;
    int bestii = 1;
    int bestjj = 1;
    float w;
    float beste = -1.0, dele;
    float e[3][3], phi, psi, dphi, dpsi, omega, domega;
    float dx, dy, dz;
    float sumi, sumj, sume;
    for (i = 0; i < n; i += 4)
    {
        beste = -1.0;
        bestii = 1;
        bestjj = 1;
        sumi = sumj = sume = 0.0;
        if (strcmp("GLY", dict.RefiPhiPsis[i].res->type().c_str()))
        {
            phi = (float)CalcAtomTorsion(dict.RefiPhiPsis[i].getAtom1(), dict.RefiPhiPsis[i].getAtom2(), dict.RefiPhiPsis[i].atom3, dict.RefiPhiPsis[i].atom4);
            psi = (float)CalcAtomTorsion(dict.RefiPhiPsis[i+1].getAtom1(), dict.RefiPhiPsis[i+1].getAtom2(), dict.RefiPhiPsis[i+1].atom3, dict.RefiPhiPsis[i+1].atom4);
            for (ii = 0; ii < 3; ii++)
            {
                for (jj = 0; jj < 3; jj++)
                {
                    e[ii][jj] = phipsi_energy(
                        phi+(ii-1)*30.0f, psi+(jj-1)*30.0f);
                    sumi += (float)ii*e[ii][jj];
                    sumj += (float)jj*e[ii][jj];
                    sume += e[ii][jj];
                    if (e[ii][jj] > beste)
                    {
                        beste = e[ii][jj];
                        bestii = ii;
                        bestjj = jj;
                    }
                }
            }
            /* if at peak skip */
            if (beste > e[1][1]+0.1f)
            {
                sumi = sumi/sume;
                sumj = sumj/sume;
                dele = beste-e[1][1];
                /* 90 = 9 points * 30.0 between points */
                dphi = (sumi-1.0f)*270.0f;
                if (dphi > 5.0f)
                {
                    dphi = 5.0f;
                }
                if (dphi < -5.0f)
                {
                    dphi = -5.0f;
                }
                dpsi = (sumj-1.0f)*270.0f;
                if (dpsi > 5.0f)
                {
                    dpsi = 5.0f;
                }
                if (dpsi < -5.0f)
                {
                    dpsi = -5.0f;
                }
                /*
                   printf("bestii = %d, bestjj = %d sumi=%f sumj=%f\n",bestii, bestjj, sumi, sumj);
                   printf("for %s phi-psi=%f %f, dphi-psi=%f %f\n", dict.RefiPhiPsis[i].res->name(),phi,psi,dphi,dpsi);
                 */
                /* figure out change to atoms */
                w = dphi*dphi/(dict.sigmatorsion*dict.sigmatorsion);
                if (w > 0.10F)
                {
                    my_dTorsion(dict.RefiPhiPsis[i].getAtom2(), dict.RefiPhiPsis[i].atom3, dict.RefiPhiPsis[i].atom4, dphi, &dx, &dy, &dz);
                    dx *= w;
                    dy *= w;
                    dz *= w;
                    dict.RefiPhiPsis[i].atom4->addDelta(dx, dy, dz);
                    dict.RefiPhiPsis[i].atom4->addWeight(w);
                }
                /*
                   my_dTorsion(dict.RefiPhiPsis[i].getAtom2(), dict.RefiPhiPsis[i].atom3, dict.RefiPhiPsis[i].getAtom1(), -dphi/2.0, &dx, &dy, &dz);
                   dx *= w;
                   dy *= w;
                   dz *= w;
                   dict.RefiPhiPsis[i].getAtom1()->dx += dx;
                   dict.RefiPhiPsis[i].getAtom1()->dy += dy;
                   dict.RefiPhiPsis[i].getAtom1()->dz += dz;
                   dict.RefiPhiPsis[i].getAtom1()->weight += w;
                 */
                w = dpsi*dpsi/(dict.sigmatorsion*dict.sigmatorsion);
                if (w > 0.10F)
                {
                    my_dTorsion(dict.RefiPhiPsis[i+1].getAtom2(), dict.RefiPhiPsis[i+1].atom3, dict.RefiPhiPsis[i+1].atom4, dpsi, &dx, &dy, &dz);
                    dx *= w;
                    dy *= w;
                    dz *= w;
                    dict.RefiPhiPsis[i+1].atom4->addDelta(dx, dy, dz);
                    dict.RefiPhiPsis[i+1].atom4->addWeight(w);
                }
                /*
                   my_dTorsion(dict.RefiPhiPsis[i+1].getAtom2(), dict.RefiPhiPsis[i+1].atom3, dict.RefiPhiPsis[i+1].getAtom1(), -dpsi/2.0, &dx, &dy, &dz);
                   dx *= w;
                   dy *= w;
                   dz *= w;
                   dict.RefiPhiPsis[i+1].getAtom1()->dx += dx;
                   dict.RefiPhiPsis[i+1].getAtom1()->dy += dy;
                   dict.RefiPhiPsis[i+1].getAtom1()->dz += dz;
                   dict.RefiPhiPsis[i+1].getAtom1()->weight += w;
                 */
            }
        }
        /* put omega at 180.0 and omega' at 0.0*/
        for (ii = 2; ii <= 3; ii++)
        {
            omega = (float)CalcAtomTorsion(dict.RefiPhiPsis[i+ii].getAtom1(), dict.RefiPhiPsis[i+ii].getAtom2(), dict.RefiPhiPsis[i+ii].atom3, dict.RefiPhiPsis[i+ii].atom4);
            domega = dict.RefiPhiPsis[i+ii].ideal[0] - omega;
            if (domega < -180.0)
            {
                domega += 360.0;
            }
            if (domega >  180.0)
            {
                domega -= 360.0;
            }
            if (fabs(domega) > 1.0)
            {
                /* clamp to avoid excessive movements in a single cycle
                 * and to preserve small angle aRefiPhiPsisroximation of dx,dy,dz */
                if (domega > 5.0)
                {
                    domega = 5.0;
                }
                if (domega < -5.0)
                {
                    domega = -5.0;
                }
                w  = domega*domega/(dict.sigmatorsion*dict.sigmatorsion);
                my_dTorsion(dict.RefiPhiPsis[i+ii].getAtom2(), dict.RefiPhiPsis[i+ii].atom3, dict.RefiPhiPsis[i+ii].atom4, domega/2.0f, &dx, &dy, &dz);
                dx *= w;
                dy *= w;
                dz *= w;
                dict.RefiPhiPsis[i+ii].atom4->addDelta(dx, dy, dz);
                dict.RefiPhiPsis[i+ii].atom4->addWeight(w);
                my_dTorsion(dict.RefiPhiPsis[i+ii].getAtom2(), dict.RefiPhiPsis[i+ii].atom3, dict.RefiPhiPsis[i+ii].getAtom1(), -domega/2.0f, &dx, &dy, &dz);
                dx *= w;
                dy *= w;
                dz *= w;
                dict.RefiPhiPsis[i+ii].getAtom1()->addDelta(dx, dy, dz);
                dict.RefiPhiPsis[i+ii].getAtom1()->addWeight(w);
                /*
                   printf("%s omega %f, domega %f (%f %f %f to %s)\n",dict.RefiPhiPsis[i+ii].res->name(),omega, domega, dx, dy, dz, dict.RefiPhiPsis[i+ii].getAtom1()->name);
                 */
            }
        }
    }
    return n;
}

int MIMolOpt::minimize_map()
{
    if (!CurrentMap)
    {
        return 0;
    }
    if (!CurrentMap->HasDensity())
    {
        return 0;
    }
    CMapHeaderBase *mh = CurrentMap->GetMapHeader();
    Residue *res;
    Residue *reslist = RefiRes;
    int nres = nRefiRes;
    float mweight = 1.0;
    int sliderweight = GetMapWeightI();
    float weight;
    float dsum = 0;
    MIAtom *a;
    float dx, dy, dz;
    float dxyz;
    float x, y, z, fx1, fy1, fz1, fx2, fy2, fz2;
    float zweight;
    float r1, r2, r;
    int n = 0, i;

    if (sliderweight == 0)
    {
        return (0);
    }
    mweight = (float)sliderweight/1000.0f;
    dxyz = mh->resmin/3.0f;

    res = reslist;
    while (Monomer::isValid(res) && n < nres)
    {
        for (i = 0; i < res->atomCount(); i++)
        {
            a = res->atom(i);
            x = a->x();
            y = a->y();
            z = a->z();
            fx1 = x -dxyz;
            fx2 = x +dxyz;
            fy1 = y -dxyz;
            fy2 = y +dxyz;
            fz1 = z -dxyz;
            fz2 = z +dxyz;
            transform(mh->ctof, &fx1, &fy1, &fz1);
            transform(mh->ctof, &fx2, &fy2, &fz2);
            transform(mh->ctof, &x, &y, &z);
            r1 = CurrentMap->avgrho(fx1, y, z);
            r2 = CurrentMap->avgrho(fx2, y, z);
            dx = (r2-r1)/10.0f * dxyz/(float)mh->nx;
            r1 = CurrentMap->avgrho(x, fy1, z);
            r2 = CurrentMap->avgrho(x, fy2, z);
            dy = (r2-r1)/10.0f * dxyz/(float)mh->ny;
            r1 = CurrentMap->avgrho(x, y, fz1);
            r2 = CurrentMap->avgrho(x, y, fz2);
            dz = (r2-r1)/10.0f * dxyz/(float)mh->nz;
            transform(mh->ftoc, &dx, &dy, &dz);
            if (mh->resmin > 2.5)
            {
                r = CurrentMap->avgrho(x, y, z);
                if (r > 75.0)
                {
                    dx /= 2.0;
                    dy /= 2.0;
                    dz /= 2.0;
                }
            }
            /* weight by slider and size relative to average
             * protein atom */
            zweight = ZByName(a->name())/6.7f;
            weight = mweight * zweight * zweight;
            dx *= weight;
            dy *= weight;
            dz *= weight;
            a->addDelta(dx, dy, dz);
            dsum += dx*dx + dy*dy + dz*dz;
            a->addWeight(weight);
        }
        n++;
        res = res->next();
    }
    if (RefiVerbose)
    {
        Logger::log("RMS map movement: %0.3f", sqrt(dsum/n));
    }
    return n;
}

//Begin ligand refinement via DE

void copy_pos(int dim, float *a, float *b)
{
    int j;
    for (j = 0; j < dim; j++)
    {
        a[j] = b[j];
    }
}

//*C*F****************************************************************
//                                                                  **
// SRC-FUNCTION   :rnd_uni()                                        **
// LONG_NAME      :random_uniform                                   **
// DESCRIPTION    :rnd_uni() generates an equally distributed ran-  **
//                 dom number in the interval [0,1]. For further    **
//                 reference see Press, W.H. et alii, Numerical     **
//                 Recipes in C, Cambridge University Press, 1992.  **
//                                                                  **
// PARAMETERS     :*idum    serves as a seed value                  **
//                                                                  **
// PRECONDITIONS  :*idum must be negative on the first call.        **
//                                                                  **
// POSTCONDITIONS :*idum will be changed                            **
//                                                                  **
//*C*F*E**************************************************************
float rnd_uni(long *idum)
{
    long j;
    long k;
    static long idum2 = 123456789;
    static long iy = 0;
    static long iv[NTAB];
    float temp;

    if (*idum <= 0)
    {
        if (-(*idum) < 1)
        {
            *idum = 1;
        }
        else
        {
            *idum = -(*idum);
        }
        idum2 = (*idum);
        for (j = NTAB+7; j >= 0; j--)
        {
            k = (*idum)/IQ1;
            *idum = IA1*(*idum-k*IQ1)-k*IR1;
            if (*idum < 0)
            {
                *idum += IM1;
            }
            if (j < NTAB)
            {
                iv[j] = *idum;
            }
        }
        iy = iv[0];
    }
    k = (*idum)/IQ1;
    *idum = IA1*(*idum-k*IQ1)-k*IR1;
    if (*idum < 0)
    {
        *idum += IM1;
    }
    k = idum2/IQ2;
    idum2 = IA2*(idum2-k*IQ2)-k*IR2;
    if (idum2 < 0)
    {
        idum2 += IM2;
    }
    j = iy/NDIV;
    iy = iv[j]-idum2;
    iv[j] = *idum;
    if (iy < 1)
    {
        iy += IMM1;
    }
    if ((temp = AM*iy) > RNMX)
    {
        return ((float)RNMX);
    }
    else
    {
        return (temp);
    }
}

float MIMolOpt::scorestep(Bond *bonds, ANGLE *angles, PLANE *planes, Bond *cons, Bond *bumps, MIAtomList atoms, int /* count */)
{
    MIAtom *a, *a1, *a2;
    int i, j;
    //int ii,jj;
    int ndim = 0;
    float r = 0;
    float d_squared = 0;
    float s_squared = 0;
    float bondscore = 0;
    float anglescore = 0;
    float planescore = 0;
    float constraintscore = 0;
    float torsionscore = 0;
    float bumpscore = 0;
    float mapscore = 0;
    float totalscore = 0;
    float chi, dchi, d;
    float *angle_s, *plane_s, *cons_s, *tors_s, *bump_s;
    float mean = 0;
    float stdev = 0;
    //float e[3][3], phi, psi, dphi, dpsi, dele;
    //float sumi, sumj, sume;

    //  FILE  *fpout = fopen("scores.txt","a");
    //  if (fpout == NULL)
    //  {
    //          printf("\nCannot open output file\n");
    //          exit(1);
    //  }
    //  fprintf(fpout,"new score ------------------\n");

    //compute bondscore
    ndim = dict.RefiBonds.size();
    float dsum = 0.0;
    if (ndim != 0)
    {
        //float bs;
        //bond_s = new float[ndim];
        for (i = 0; i < ndim; i++)
        {
            a1 = bonds[i].getAtom1();
            a2 = bonds[i].getAtom2();
            r = (float)AtomDist(*a1, *a2);
            d_squared = fabs(r - bonds[i].ideal_length);
            d_squared *= d_squared;
            dsum += d_squared;
            s_squared = bonds[i].tolerance;
            //s_squared *= s_squared;
            //if(d_squared > 0 )
            //	bs = s_squared/d_squared;
            //else
            //	bs = 100.0* s_squared;
            //fprintf(fpout, "d_sq = %f s_sq = %0f score = %f\n", d_squared, s_squared, bs);
            //bond_s[i] = d_squared/s_squared;
            //	mean += bond_s[i];
            //bondscore += bond_s[i];
        }
        bondscore = dsum / (float)ndim;
        /*
           mean = mean/ndim;
           for(i=0; i < ndim; i++)
            stdev += (bond_s[i]-mean) * (bond_s[i] - mean);
           if (stdev != 0){
            stdev = stdev / (ndim-1);
            stdev = sqrt(stdev);
            bondscore = bondscore/stdev;
           }
           //delete[] bond_s;
           //mean = 0;
           stdev = 0;
         */
    }

    //compute anglescore
    ndim = dict.RefiAngles.size();
    if (ndim != 0)
    {
        angle_s = new float[ndim];
        for (i = 0; i < ndim; i++)
        {
            a1 = angles[i].getAtom1();
            a2 = angles[i].getAtom2();
            r = (float)AtomDist(*a1, *a2);
            d_squared = r - angles[i].ideal_angle;
            d_squared *= d_squared;
            s_squared = angles[i].tolerance;
            s_squared *= s_squared;
            angle_s[i] = d_squared/s_squared;
            mean += angle_s[i];
            anglescore += angle_s[i];
        }
        mean = mean/ndim;
        for (i = 0; i < ndim; i++)
        {
            stdev += (angle_s[i]-mean) * (angle_s[i] - mean);
        }
        if (stdev != 0)
        {
            stdev = stdev / (ndim-1);
            stdev = sqrt(stdev);
            anglescore = anglescore/stdev;
        }
        delete[] angle_s;
        mean = 0;
        stdev = 0;
    }

    //compute planescore
    ndim = dict.RefiPlanes.size();
    if (ndim != 0)
    {
        plane_s = new float[ndim];
        for (i = 0; i < ndim; i++)
        {
            lsqplane(planes[i]);
            for (j = 0; j < planes[i].natoms; j++)
            {
                a = planes[i].atoms[j];
                d_squared = a->x() * planes[i].vm[0];
                d_squared += a->y() * planes[i].vm[1];
                d_squared += a->z() * planes[i].vm[2];
                d_squared -= planes[i].d;
                d_squared *= d_squared;
            }
            s_squared = planes[i].tolerance;
            s_squared *= s_squared;
            plane_s[i] = d_squared/s_squared;
            mean += plane_s[i];
            planescore += plane_s[i];
        }
        mean = mean/ndim;
        for (i = 0; i < ndim; i++)
        {
            stdev += (plane_s[i]-mean) * (plane_s[i] - mean);
        }
        if (stdev != 0)
        {
            stdev = stdev / (ndim-1);
            stdev = sqrt(stdev);
            planescore = planescore/stdev;
        }
        delete[] plane_s;
        mean = 0;
        stdev = 0;
    }

    //compute constraintscore
    ndim = dict.RefiConstraints.size();
    if (ndim != 0)
    {
        cons_s = new float[ndim];
        for (i = 0; i < ndim; i++)
        {
            a1 = cons[i].getAtom1();
            a2 = cons[i].getAtom2();
            r = (float)AtomDist(*a1, *a2);
            d_squared = r - cons[i].ideal_length;
            d_squared *= d_squared;
            s_squared = cons[i].tolerance;
            s_squared *= s_squared;
            cons_s[i] = d_squared/s_squared;
            mean += cons_s[i];
            constraintscore += cons_s[i];
        }
        mean = mean/ndim;
        for (i = 0; i < ndim; i++)
        {
            stdev += (cons_s[i]-mean) * (cons_s[i] - mean);
        }
        if (stdev != 0)
        {
            stdev = stdev / (ndim-1);
            stdev = sqrt(stdev);
            constraintscore = constraintscore/stdev;
        }
        delete[] cons_s;
        mean = 0;
        stdev = 0;
    }


    //compute phipsiscore
    /*
       for(i=0; i< RefiPhiPsis.size(); i+=4){
        beste = -1.0;
        bestii = 1;
        bestjj = 1;
        sumi = sumj = sume = 0.0;
        if(strcmp("GLY",RefiPhiPsis[i].res->type())){
            phi = CalcAtomTorsion(RefiPhiPsis[i].getAtom1(),RefiPhiPsis[i].getAtom2(),RefiPhiPsis[i].atom3,RefiPhiPsis[i].atom4);
            psi = CalcAtomTorsion(RefiPhiPsis[i+1].getAtom1(),RefiPhiPsis[i+1].getAtom2(),RefiPhiPsis[i+1].atom3,RefiPhiPsis[i+1].atom4);
            for(ii=0; ii<3; ii++){
                for(jj=0; jj<3; jj++){
                    e[ii][jj] = phipsi_energy(phi+(ii-1)*30.0, psi+(jj-1)*30.0);
                    sumi += (float)ii*e[ii][jj];
                    sumj += (float)jj*e[ii][jj];
                    sume += e[ii][jj];
                    if(e[ii][jj]>beste){
                        beste = e[ii][jj];
                        bestii = ii;
                        bestjj = jj;
                    }
                }
            }
            phipsiscore -= sume;
        }
       }
     */

    //compute torsionscore
    ndim = dict.RefiTorsions.size();
    if (ndim != 0)
    {
        tors_s = new float[ndim];
        for (i = 0; i < ndim; i++)
        {
            chi = (float)CalcAtomTorsion(dict.RefiTorsions[i].getAtom1(), dict.RefiTorsions[i].getAtom2(), dict.RefiTorsions[i].atom3, dict.RefiTorsions[i].atom4);
            if (chi < 0.0)
            {
                chi += 360.0;
            }
            dchi = 0.0;
            for (j = 0; j < dict.RefiTorsions[i].nideal; j++)
            {
                d = dict.RefiTorsions[i].ideal[j] - chi;
                if (d < -180.0)
                {
                    d += 360.0;
                }
                if (d >  180.0)
                {
                    d -= 360.0;
                }
                if (fabs(d) < fabs(dchi) || j == 0)
                {
                    dchi = d;
                }
            }
            if (fabs(dchi) > 0.1F)
            {
                d_squared = dchi * dchi;
                s_squared = dict.sigmatorsion * dict.sigmatorsion * TorsionWeight;
                tors_s[i] = d_squared/s_squared;
                mean += tors_s[i];
                torsionscore += tors_s[i];
            }
        }
        mean = mean/ndim;
        for (i = 0; i < ndim; i++)
        {
            stdev += (tors_s[i]-mean) * (tors_s[i] - mean);
        }
        if (stdev != 0)
        {
            stdev = stdev / (ndim-1);
            stdev = sqrt(stdev);
            torsionscore = torsionscore/stdev;
        }
        delete[] tors_s;
        mean = 0;
        stdev = 0;
    }

    //compute bumpscore
    ndim = dict.RefiBumps.size();
    if (ndim != 0)
    {
        bump_s = new float[ndim];
        for (i = 0; i < ndim ; i++)
        {
            a1 = bumps[i].getAtom1();
            a2 = bumps[i].getAtom2();
            r = (float)AtomDist(*a1, *a2);
            d_squared = r - bumps[i].ideal_length;
            d_squared *= d_squared;
            s_squared = bumps[i].tolerance;
            s_squared *= s_squared;
            bump_s[i] = d_squared/s_squared;
            mean += bump_s[i];
            bumpscore += bump_s[i];
        }
        mean = mean/ndim;
        for (i = 0; i < ndim; i++)
        {
            stdev += (bump_s[i]-mean) * (bump_s[i] - mean);
        }
        if (stdev != 0)
        {
            stdev = stdev / (ndim-1);
            stdev = sqrt(stdev);
            bumpscore = bumpscore/stdev;
        }
        delete[] bump_s;
        mean = 0;
        stdev = 0;

    }
    //if (bumpscore > 9.0) bumpscore = 9.0;

    //compute mapscore
    if (CurrentMap && CurrentMap->HasDensity())
    {
        mapscore = -(CurrentMap->RDensity(atoms));
    }

    totalscore = -bondscore  //+ anglescore //+ planescore + constraintscore //+ phipsiscore
                 /*+ torsionscore + bumpscore + mapscore*/;

    // fprintf(fpout,"%f,%f,%f,%f,",bondscore,anglescore,planescore,constraintscore);
    // fprintf(fpout,"%f,%f,%f,%f,%d\n",torsionscore,bumpscore,mapscore,totalscore,count);
    // fclose(fpout);

    return (totalscore);
}

int MIMolOpt::takestep(int seed)
{
    const float inibound_h = 0.05F;
    const float inibound_l = -0.05F;
    int i, j, k, r1, r2, r3;
    int n = 0;
    int imin;
    int gen = 0;
    int maxgen = 10;
    int numparams;        //D in DE jargon
    int popsize = 20;
    float F = 0.8F;
    float CR = 0.9F;
    float cmin = 0;
    float trial_cost, cmean, cvar;
    float *best, *bestit, *cost, *tmp;
    float *pold;
    float *pnew;
    float *pswap;
    long init;
    vector<MIAtom*> refiAtoms;
    MIAtomList atoms;
    Residue *res = RefiRes;
    MIAtom *a;
    //FILE	*fpin_ptr, *fpout;

    /*-----Read input data------------------------------------------------*/

    //fpin_ptr   = fopen("parameters.txt","r");

    //if (fpin_ptr == NULL)
    //{
    //	printf("\nCannot open input file\n");
    //	exit(1);                                 /* input file is necessary */
    //}

    //fpout   = fopen("output.txt","w");

    //if (fpout == NULL)
    //{
    //	printf("\nCannot open output file\n");
    //	exit(1);                                 /* output file is necessary */
    //}

    //fscanf(fpin_ptr,"%d",&strategy);       /*---choice of strategy-----------------*/
    //fscanf(fpin_ptr,"%d,",&maxgen);        /*---maximum number of generations------*/
    //fscanf(fpin_ptr,"%d,",&refresh);       /*---output refresh cycle---------------*/
    //fscanf(fpin_ptr,"%d,",&popsize);       /*---population size.-------------------*/
    //fscanf(fpin_ptr,"%f,",&F);             /*---weight factor----------------------*/
    //fscanf(fpin_ptr,"%f",&CR);             /*---crossing over factor---------------*/

    //fclose(fpin_ptr);

    while (Monomer::isValid(res) && n < nRefiRes)
    {
        for (i = 0; i < res->atomCount(); i++)
        {
            refiAtoms.push_back(res->atom(i));
        }
        n++;
        res = res->next();
    }

    numparams = refiAtoms.size()*3;
    pold = new float[popsize*numparams];
    pnew = new float[popsize*numparams];
    pswap = new float[popsize*numparams];
    best = new float[numparams];
    bestit = new float[numparams];
    cost = new float[popsize];
    tmp = new float[numparams];

    /*------Initialization------------------------------------------------*/

    init = -(long) seed;
    int y;
    for (i = 0; i < popsize; i++)
    {
        for (j = 0; j < numparams; j++) /* spread initial population members */
        {
            y = (i*numparams) + j;
            // the first daughter stays at the start otherwise add a random number
            if (i != 0)
            {
                pold[y] = inibound_l + rnd_uni(&init)*(inibound_h - inibound_l);
            }
            else
            {
                pold[y] = 0;
            }
        }
    }
    for (j = 0; j < popsize; j++)
    {
        k = j * numparams;
        for (unsigned int i = 0; i < refiAtoms.size(); i++)
        {
            a = refiAtoms[i];
            a->translate(pold[k], pold[k+1], pold[k+2]);
            k += 3;
        }
        cost[j] = scorestep(&dict.RefiBonds[0], &dict.RefiAngles[0], &dict.RefiPlanes[0], &dict.RefiConstraints[0], &dict.RefiBumps[0], refiAtoms, gen); /* obj. funct. value */
        k = j * numparams;
        for (unsigned int i = 0; i < refiAtoms.size(); i++)
        {
            a = refiAtoms[i];
            a->translate( -pold[k], -pold[k+1], -pold[k+2]);
            k += 3;
        }
    }
    cmin = cost[0];
    imin = 0;
    for (j = 1; j < popsize; j++)
    {
        if (cost[j] < cmin)
        {
            cmin = cost[j];
            imin = j;
        }
    }
    copy_pos(numparams, best, &pold[imin]);          /* save best member ever          */
    copy_pos(numparams, bestit, &pold[imin]);        /* save best member of generation */
    //fprintf(fpout,"Generation,Lowest Score,Score Variance\n");
    //fprintf(fpout,"%d,%f\n",gen,cmin);

    /*=======================================================================*/
    /*=========Iteration loop================================================*/
    /*=======================================================================*/

    for (gen = 1; gen < maxgen; gen++)
    {
        imin = 0;

        for (j = 0; j < popsize; j++)       // Start of loop through ensemble
        {
            copy_pos(numparams, tmp, &pold[j*numparams]); //Copy old member of population into new trial

            do                              //Pick a random population member
            {
                r1 = (int)(rnd_uni(&init)*popsize);
            } while (r1 == j);

            do                              // Pick another random population member
            {
                r2 = (int)(rnd_uni(&init)*popsize);
            } while ((r2 == j) || (r2 == r1));

            do                              // Pick another random population member
            {
                r3 = (int)(rnd_uni(&init)*popsize);
            } while ((r3 == j) || (r3 == r1) || (r3 == r2));
            r1 *= numparams;
            r2 *= numparams;
            r3 *= numparams;
            n = (int)(rnd_uni(&init)*numparams);

            for (k = 0; (rnd_uni(&init) < CR) && (k < numparams); k++)
            {
                tmp[n] = pold[r1+n] + F * (pold[r2+n] - pold[r3+n]);
                n = (n+1)%numparams;
            }

            //=======Trial mutation now in tmp[]. Test how good this choice really was.==================

            for (unsigned int i = 0; i < refiAtoms.size(); i++)
            {
                a = refiAtoms[i];
                k = i * 3;
                a->translate(tmp[k], tmp[k+1], tmp[k+2]);
            }

            trial_cost = scorestep(&dict.RefiBonds[0], &dict.RefiAngles[0], &dict.RefiPlanes[0], &dict.RefiConstraints[0], &dict.RefiBumps[0], refiAtoms, gen); // Evaluate new vector in tmp[]

            for (unsigned int i = 0; i < refiAtoms.size(); i++)
            {
                a = refiAtoms[i];
                k = i * 3;
                a->translate(-tmp[k], -tmp[k+1], -tmp[k+2]);
            }

            if (trial_cost < cost[j])   // improved objective function value ?
            {
                cost[j] = trial_cost;
                copy_pos(numparams, &pnew[j*numparams], tmp);
                if (trial_cost < cmin)         // Was this a new minimum?
                {
                    cmin = trial_cost;         // if so, reset cmin to new low...
                    imin = j;
                    copy_pos(numparams, best, tmp);
                }
            }
            else
            {
                copy_pos(numparams, &pnew[j*numparams], &pold[j*numparams]); // replace target with old value
            }
        }   // End mutation loop through pop.

        copy_pos(numparams, bestit, best); // Save best population member of current iteration

        // swap population arrays. New generation becomes old one

        copy_pos(popsize*numparams, pswap, pold);
        copy_pos(popsize*numparams, pold, pnew);
        copy_pos(popsize*numparams, pnew, pswap);

        //----Compute the energy variance (just for monitoring purposes)-----------

        cmean = 0;          // compute the mean value first
        for (j = 0; j < popsize; j++)
        {
            cmean += cost[j];
        }
        cmean /= popsize;

        cvar = 0;           // now the variance
        for (j = 0; j < popsize; j++)
        {
            cvar += (cost[j] - cmean)*(cost[j] - cmean);
        }
        cvar /= popsize-1;

        //if (gen%refresh == 0) fprintf(fpout,"\n%d,%f,%f\n",gen,cmin,cvar);
    }
    //fclose(fpout);

    for (unsigned int i = 0; i < refiAtoms.size(); i++)
    {
        a = refiAtoms[i];
        k = i * 3;
        a->translate(bestit[k], bestit[k+1], bestit[k+2]);
    }

    delete[] pold;
    delete[] pnew;
    delete[] pswap;
    delete[] best;
    delete[] bestit;
    delete[] cost;
    delete[] tmp;
    return (0);
}

int MIMolOpt::resetatomderivatives(MIAtom *atoms[], int natoms)
{
    int i;
    for (i = 0; i < natoms; i++)
    {
        atoms[i]->resetDelta();
        atoms[i]->resetWeight();
    }
    return 1;
}

int MIMolOpt::resetderivatives()
{
    if (nRefiRes == 0 || RefiRes == NULL)
    {
        return 0;
    }
    Residue *res = RefiRes;
    int i;
    int n = 0;
    while (Monomer::isValid(res) && n < nRefiRes)
    {
        for (i = 0; i < res->atomCount(); i++)
        {
            res->atom(i)->resetDelta();
            res->atom(i)->resetWeight();
        }
        n++;
        res = res->next();
    }
    return 1;
}

int MIMolOpt::applyderivatives()
{
    if (nRefiRes == 0 || RefiRes == NULL)
    {
        return 0;
    }
    Residue *res = RefiRes;
    int i;
    int n = 0;
    MIAtom *a;
    float dx, dy, dz;
    float maxd = 1.0F;
    while (Monomer::isValid(res) && n < nRefiRes)
    {
        for (i = 0; i < res->atomCount(); i++)
        {
            a = res->atom(i);
            if (a->weight() < 0.001F)
            {
                a->resetDelta();
                a->resetWeight();
                a->addWeight(1.0f);
                continue;
            }
            dx = a->dx()/a->weight();
            if (dx > maxd)
            {
                dx = maxd;
            }
            if (dx < -maxd)
            {
                dx = -maxd;
            }
            dy = a->dy()/a->weight();
            if (dy > maxd)
            {
                dy = maxd;
            }
            if (dy < -maxd)
            {
                dy = -maxd;
            }
            dz = a->dz()/a->weight();
            if (dz > maxd)
            {
                dz = maxd;
            }
            if (dz < -maxd)
            {
                dz = -maxd;
            }
            a->translate(dx, dy, dz);
            a->resetDelta();
            a->addDelta(dx, dy, dz);
            a->resetWeight();
            a->addWeight(1.0f);
        }
        n++;
        res = res->next();
    }
    return 1;
}

void MIMolOpt::ConnectTo(MIMoleculeBase *mol)
{
    connect(mol, SIGNAL(moleculeDeleted(chemlib::MIMoleculeBase*)),
            this, SLOT(moleculeDeleted(chemlib::MIMoleculeBase*)));
    connect(mol, SIGNAL(moleculeToBeDeleted(chemlib::MIMoleculeBase*)),
            this, SLOT(moleculeToBeDeleted(chemlib::MIMoleculeBase*)));
    connect(mol, SIGNAL(residuesToBeDeleted(chemlib::MIMoleculeBase*, std::vector<chemlib::Residue*>&)),
            this, SLOT(residuesToBeDeleted(chemlib::MIMoleculeBase*, std::vector<chemlib::Residue*>&)));
    connect(mol, SIGNAL(atomsToBeDeleted(chemlib::MIMoleculeBase*, chemlib::MIAtomList)),
            this, SLOT(atomsToBeDeleted(chemlib::MIMoleculeBase*, chemlib::MIAtomList)));
}


long MIMolOpt::SetRefiRes(Residue *res1, Residue *res2, MIMoleculeBase *node, EMapBase *emap)
{
    Residue *res;
    Residue *start = NULL;
    int nres;
    if (dict.EmptyDictCheck() == false)
    {
        return 0;
    }

    if (!res1 || !res2)
    {
        return (0);
    }
    /* find correct order */
    if (res1 == res2)
    {
        nres = 1;
        start = res1;
    }
    else
    {
        res = res1;
        nres = 0;
        while (res != NULL)
        {
            nres++;
            if (res == res2)
            {
                start = res1;
                break;
            }
            res = res->next();
        }
        if (!start)
        {
            res = res2;
            nres = 0;
            while (res != NULL)
            {
                nres++;
                if (res == res1)
                {
                    start = res2;
                    break;
                }
                res = res->next();
            }
        }
    }
    if (!start)
    {
        Logger::message("Program can't make sense of start and end - not from same model?");
        return (0);
    }
    res = start;
    int i = 0, j;
    // don't let symm atoms be fit - they are not permanent nor saveable!
    while (res != NULL && i < nres)
    {
        for (j = 0; j < res->atomCount(); j++)
        {
            if (res->atom(j)->type() & AtomType::SYMMATOM)
            {
                return 0;
            }
        }
        i++;
        res = res->next();
    }
    internalSetRefiRes(start, nres);
    while (start != NULL)
    {
        if (dict.GetDictResidue(start == NULL))
        {
            return 0;
        }
        nres--;
        if (nres <= 0)
        {
            start = NULL;
        }
        else
        {
            start = start->next();
        }
    }
    CurrentModel = node;

    CurrentMap = emap;
    ResActiveModel = node->getResidues();
    dict.FindGeom(RefiRes, nRefiRes, ResActiveModel);
    if (dict.constrain_CA)
    {
        dict.ConstrainCalpha(RefiRes, nRefiRes);
    }
    if (dict.constrain_Ends)
    {
        dict.RestrainEnds(RefiRes, nRefiRes);
    }
    dict.BuildBumps(RefiRes, nRefiRes);
    Do();
    return nRefiRes;
}

void MIMolOpt::moleculeDeleted(MIMoleculeBase *molecule)
{
}

void MIMolOpt::moleculeToBeDeleted(MIMoleculeBase *molecule)
{
    Purge(molecule);

    std::vector<Residue*> residues;
    for (MIIter<Residue> res = molecule->GetResidues(); res; ++res)
    {
        residues.push_back(res);
    }
    residuesToBeDeleted(molecule, residues);
}

void MIMolOpt::residuesToBeDeleted(MIMoleculeBase *mol, std::vector<Residue*> &residues)
{

    MIAtomList atoms;
    for (size_t i = 0; i < residues.size(); ++i)
    {
        Purge(residues[i]);
        atoms.insert(atoms.end(), residues[i]->atoms().begin(), residues[i]->atoms().end());
    }
    atomsToBeDeleted(mol, atoms);
}

void MIMolOpt::atomsToBeDeleted(MIMoleculeBase* /* mol */, const MIAtomList &atoms)
{
    for (size_t i = 0; i < atoms.size(); ++i)
    {
        Purge(atoms[i]);
    }
}

void MIMolOpt::RefiAllTorsions(Residue *reslist)
{
    int i;
    TORSION *t;
    while (Monomer::isValid(reslist))
    {
        if (dict.DictMap.find(reslist->type()) == dict.DictMap.end())
        {
            Logger::log("residue type %s not found in dictionary",
                        reslist->type().c_str());
        }

        for (i = 0; (unsigned int) i < dict.TorsDict.size(); i++)
        {
            if (strcmp(dict.TorsDict[i].restype, reslist->type().c_str()))
            {
                continue;
            }
            if (find_if(dict.RefiTorsions.begin(), dict.RefiTorsions.end(),
                        bind2nd(TorsionMatch(), dict.TorsDict[i])) !=
                dict.RefiTorsions.end())
            {
                continue;
            }

            t = dict.getTORSION(reslist, dict.TorsDict[i].type);
            if (t != NULL)
            {
                strcpy(t->type, dict.TorsDict[i].type);
                dict.RefiTorsions.push_back(*t);
            }

            if (t)
            {
                free(t);
                t = NULL;
            }
        } //Loop over dictionary torsions
        reslist = reslist->next();
    } //Loop over residues
}

float bond_dist(MIAtom *atom1, MIAtom *atom2)
{
    float ideal = 0.0;

    ideal = bond_radius(atom1)+ bond_radius(atom2);
    return (ideal);
}

float bond_radius(MIAtom *atom)
{
    /* radius within which an atom can be considered covalently bound */
    float r = 1.3F; /* long default betting that most light atoms are covered*/
    if (strncmp(atom->name(), "C", 1) == 0)
    {
        r = 0.95F;
    }
    else if (strncmp(atom->name(), "c", 1) == 0)
    {
        r = 0.95F;
    }
    else if (strncmp(atom->name(), "N", 1) == 0)
    {
        r = 0.95F;
    }
    else if (strncmp(atom->name(), "n", 1) == 0)
    {
        r = 0.95F;
    }
    else if (strncmp(atom->name(), "O", 1) == 0)
    {
        r = 0.95F;
    }
    else if (strncmp(atom->name(), "o", 1) == 0)
    {
        r = 0.95F;
    }
    else if (strncmp(atom->name(), "F", 1) == 0)
    {
        r = 1.28F;
    }
    else if (strncmp(atom->name(), "f", 1) == 0)
    {
        r = 1.28F;
    }
    else if (strncmp(atom->name(), "S", 1) == 0)
    {
        r = 1.20F;
    }
    else if (strncmp(atom->name(), "s", 1) == 0)
    {
        r = 1.20F;
    }
    else if (strncmp(atom->name(), "H", 1) == 0)
    {
        r = 0.50F;
    }
    else if (strncmp(atom->name(), "h", 1) == 0)
    {
        r = 0.50F;
    }
    else if (strncmp(atom->name(), "P", 1) == 0)
    {
        r = 1.20F;
    }
    else if (strncmp(atom->name(), "p", 1) == 0)
    {
        r = 1.20F;
    }
    else if (strncmp(atom->name(), "CU", 2) == 0)
    {
        r = 1.25F;
    }
    else if (strncmp(atom->name(), "cu", 2) == 0)
    {
        r = 1.25F;
    }
    else if (strncmp(atom->name(), "FE", 2) == 0)
    {
        r = 1.32F;
    }
    else if (strncmp(atom->name(), "fe", 2) == 0)
    {
        r = 1.32F;
    }
    else if (strncmp(atom->name(), "ZN", 2) == 0)
    {
        r = 1.32F;
    }
    else if (strncmp(atom->name(), "zn", 2) == 0)
    {
        r = 1.32F;
    }
    else if (strncmp(atom->name(), "MN", 2) == 0)
    {
        r = 1.10F;
    }
    else if (strncmp(atom->name(), "mn", 2) == 0)
    {
        r = 1.10F;
    }
    else if (strncmp(atom->name(), "1H", 2) == 0)
    {
        r = 0.5F;
    }
    else if (strncmp(atom->name(), "2H", 2) == 0)
    {
        r = 0.5F;
    }
    else if (strncmp(atom->name(), "3H", 2) == 0)
    {
        r = 0.5F;
    }
    else if (strncmp(atom->name(), "4H", 2) == 0)
    {
        r = 0.5F;
    }
    else if (strncmp(atom->name(), "5H", 2) == 0)
    {
        r = 0.5F;
    }
    return (r);
}

int MIMolOpt::getbonddist(Residue *res, Bond *bond)
{
    MIAtom *a1 = NULL, *a2 = NULL;
    int i;
    //char buf[100];
    /* special case N - C peptide bond */
    if ((!strcmp(bond->getAtom1()->name(), "C") && !strcmp(bond->getAtom2()->name(), "N"))
        || (!strcmp(bond->getAtom1()->name(), "N") && !strcmp(bond->getAtom2()->name(), "C")))
    {
        bond->ideal_length = 1.32F;
        bond->tolerance = 1.32F*0.20F;
        return (1);
    }
    for (i = 0; i < res->atomCount(); i++)
    {
        if (!strcmp(res->atom(i)->name(), bond->getAtom1()->name()))
        {
            a1 = res->atom(i);
        }
        if (!strcmp(res->atom(i)->name(), bond->getAtom2()->name()))
        {
            a2 = res->atom(i);
        }
    }
    if (a1 == NULL || a2 == NULL)
    {
        /* in the absence of any information just use
         * the incoming distance */
        bond->ideal_length = 0.0;
        bond->tolerance = 0.2F;
        Logger::log("Can't find bond %s %s: %s to %s in dictionary\n", res->type().c_str(), res->name().c_str(), bond->getAtom1()->name(), bond->getAtom2()->name());
        return (0);
    }
    bond->ideal_length = (float)AtomDist(*a1, *a2);
    bond->tolerance = 0.2f*bond->ideal_length;
    return (1);
}

void MIMolOpt::Do()
{
    ConnectTo(CurrentModel); //Connect up signals
    SaveToken = geomsaver.Save(RefiRes, nRefiRes, CurrentModel);

    Residue *res = RefiRes;
    int i = 0, j;
    while (res != NULL && i < nRefiRes)
    {
        for (j = 0; j < res->atomCount(); j++)
        {
            res->atom(j)->addType(AtomType::REFIATOM);
        }
        i++;
        res = res->next();
    }
}

bool MIMolOpt::Undo()
{
    if (SaveToken != 0)
    {
        if (geomsaver.Restore(SaveToken))
            SaveToken--;
        Logger::log("SaveToken is now %d", (int)SaveToken);
        return true;
    }
    return false;
}

bool MIMolOpt::Redo()
{
    if (SaveToken && (int)(SaveToken+1) < geomsaver.NumberSets())
    {
        if (geomsaver.Restore(SaveToken) && MIMoleculeBase::isValid(geomsaver.Model(SaveToken)) && geomsaver.Model(SaveToken)->Build())
            SaveToken++;
        Logger::log("SaveToken is now %d", (int)SaveToken);
        return true;
    }
    return false;
}

void MIMolOpt::Accept()
{
    if (!IsRefining())
    {
        return;
    }

    ConnectTo(CurrentModel); //Connect up signals
    unsigned long newSaveToken = geomsaver.Save(RefiRes, nRefiRes, CurrentModel);
    geomsaver.RestoreColor(SaveToken, AtomType::REFIATOM);
    SaveToken = newSaveToken;
    if (CurrentModel)
    {
        CurrentModel->SetCoordsChanged(true);
        CurrentModel->SetModified(true);
    }
    clearRefineTarget();
}

void MIMolOpt::Cancel()
{
    if (!IsRefining())
    {
        return;
    }
    geomsaver.RestoreColor(SaveToken, AtomType::REFIATOM);
    if (geomsaver.Restore(SaveToken))
        SaveToken--;
    clearRefineTarget();
}

void MIMolOpt::lockRefineTarget()
{
    refineTargetLocked = true;
}

void MIMolOpt::unlockRefineTarget()
{
    refineTargetLocked = false;
}

void MIMolOpt::clearRefineTarget()
{
    if (!refineTargetLocked)
    {
        internalSetRefiRes(NULL, 0);
        CurrentModel = NULL;
        dict.Clear();
    }
}

void MIMolOpt::Reset()
{
    if (!IsRefining())
    {
        return;
    }
    geomsaver.Restore(SaveToken);
}

void MIMolOpt::Purge(MIMoleculeBase *node)
{
    if (CurrentModel == node)
    {
        clearRefineTarget();
    }
    geomsaver.Purge(node);
}

void MIMolOpt::Purge(Residue *res)
{
    if (IsRefining())
    {
        if (CurrentModel->Contains(res))
        {
            clearRefineTarget();
        }
    }
    for (int i = 0; i < res->atomCount(); i++)
    {
        geomsaver.Purge(res->atom(i));
    }
}

void MIMolOpt::Purge(MIAtom *atom)
{
    if (IsRefining())
    {
        //this will only get hit the for the first atom in the mol, after which
        // IsRefining will be false, so it's not too expensive to do this
        Residue *res = residue_from_atom(CurrentModel->getResidues(), atom);
        if (res)
        {
            clearRefineTarget();
        }
    }

    // and this isn't too expensive either; if the residue is purged first,
    // there's nothing for this to do
    geomsaver.Purge(atom);
}






void MIMolOpt::Purge(EMapBase *emap)
{
    if (CurrentMap == emap)
    {
        CurrentMap = NULL;
    }
}

static void GetCenter(MIAtomList &atoms, float &cx, float &cy, float &cz)
{
    cx = cy = cz = 0;
    for (unsigned int i = 0; i < atoms.size(); i++)
    {
        cx += atoms[i]->x();
        cy += atoms[i]->y();
        cz += atoms[i]->z();
    }
    cx /= (float)atoms.size();
    cy /= (float)atoms.size();
    cz /= (float)atoms.size();
}

struct trial
{
    double *p;
    double score;
};


static bool trial_compare(const trial &lt, const trial &rt)
{
    return lt.score > rt.score;
}

#define X 0
#define Y 1
#define Z 2
static void RotateAtomVec(float rx, float ry, float rz, float cx, float cy, float cz, MIAtomList *atoms)
{
    float mat[3][3];
    buildmat(rx, ry, rz, mat);
    orthomatrix(mat, mat);
    float xrot, yrot, zrot, xdir, ydir, zdir;
    MIAtom *a;
    if (atoms)
    {
        for (unsigned int i = 0; i < atoms->size(); i++)
        {
            a = (*atoms)[i];
            xdir = a->x() - cx;
            ydir = a->y() - cy;
            zdir = a->z() - cz;
            xrot =  (xdir*mat[X][X]+ydir*mat[X][Y]+zdir*mat[X][Z]);
            yrot =  (xdir*mat[Y][X]+ydir*mat[Y][Y]+zdir*mat[Y][Z]);
            zrot =  (xdir*mat[Z][X]+ydir*mat[Z][Y]+zdir*mat[Z][Z]);
            a->setPosition(xrot + cx, yrot + cy, zrot + cz);
        }
    }
}

#undef Z
#undef Y
#undef X

//these overloads are to avoid type truncation warnings
static void RotateAtomVec(double rx, double ry, double rz, float cx, float cy, float cz, MIAtomList *atoms)
{
    RotateAtomVec((float)rx, (float)ry, (float)rz, (float)cx, (float)cy, (float)cz, atoms);
}

//static void RotateAtomVec(double rx, double ry, double rz, double cx, double cy, double cz, MIAtomList * atoms) {
//  RotateAtomVec((float)rx, (float)ry, (float)rz, (float)cx, (float)cy, (float)cz, atoms);
//}

static void TranslateAtomVec(float x, float y, float z, MIAtomList *atoms)
{
    MIAtom *atom;
    for (unsigned int i = 0; i < atoms->size(); i++)
    {
        atom = (*atoms)[i];
        atom->translate(x, y, z);
    }
}

static void TranslateAtomVec(double x, double y, double z, MIAtomList *atoms)
{
    TranslateAtomVec((float)x, (float)y, (float)z, atoms);
}

void score(trial &t, MIAtomList &atoms, float cx, float cy, float cz, EMapBase *emap)
{
    RotateAtomVec(t.p[0], t.p[1], t.p[2], cx, cy, cz, &atoms);
    TranslateAtomVec(t.p[3], t.p[4], t.p[5], &atoms);
    t.score = emap->RDensity(atoms);
}

void copy_trial(trial *t1, trial *t2, unsigned int td)
{
    t1->score = t2->score;
    for (unsigned int i = 0; i < td; i++)
    {
        t1->p[i] = t2->p[i];
    }
}

void MIMolOpt::RigidOptimize(MIAtomList &CurrentAtoms, MIMoleculeBase *fitmol, EMapBase *emap)
{
    CurrentMap = emap;
    float maxt = 1.0F;
    float maxr = 1.5F;
    float cx = 0, cy = 0, cz = 0;
    unsigned int i, j, k;
    trial t;
    unsigned int p1, p2, p3;
    trial *ti, *ta, *tb, *tc, *t2;
    float F = 0.7F;
    float CR = 0.8F;
    int D = 6;
    unsigned int npop = 20*D;

    if (!fitmol || !emap)
    {
        return;
    }
    if (!emap->HasDensity() || CurrentAtoms.size() == 0)
    {
        return;
    }

    float rdens_start = CurrentMap->RDensity(CurrentAtoms);

    ConnectTo(fitmol); //Connect up signals
    SaveToken = geomsaver.Save(CurrentAtoms, fitmol);

    // find the center of mass of the atoms
    GetCenter(CurrentAtoms, cx, cy, cz);

    float damp = 1.0, param;
    vector<float> var_start;
    var_start.reserve(D);
    for (j = 0; j < 3; j++)
    {
        param = maxr/damp;
        var_start.push_back(param);
    }
    for (j = 3; j < 6; j++)
    {
        param = maxt/damp;
        var_start.push_back(param);
    }

    vector<trial> population;
    vector<trial> population2;
    vector<trial>::iterator p;
    double best_r, worst_r;

    for (i = 0; i < npop; i++)
    {
        t.p = new double[D];
        for (j = 0; j < (unsigned int)D; j++)
        {
            t.p[j] = grand(var_start[j]);
        }
        population.push_back(t);
        t.p = new double[D];
        population2.push_back(t);
    }

    for (p = population.begin(); p != population.end(); p++)
    {
        score(*p, CurrentAtoms, cx, cy, cz, emap);
        geomsaver.Restore(SaveToken);
    }
    sort(population.begin(), population.end(), trial_compare);
    best_r = population[0].score;
    worst_r = population[population.size()-1].score;

    for (int igen = 0; igen < 40; igen++)
    {
        for (i = 0; i < npop; i++)
        {
            ti = &population[i];
            do
            {
                p1 = irand(npop);
            } while (p1 == i);
            ta = &population[p1];
            do
            {
                p2 = irand(npop);
            } while (p2 == i || p2 == p1);
            tb = &population[p2];
            do
            {
                p3 = irand(npop);
            } while (p3 == i || p3 == p1 || p3 == p2);
            tc = &population[p3];
            j = irand(D);
            for (k = 1; k <= (unsigned int)D; k++)
            {
                if (frand() < CR || k == (unsigned int)D)
                {
                    t.p[j] = tc->p[j] + F*(ta->p[j]-tb->p[j]);
                }
                else
                {
                    t.p[j] = ti->p[j];
                }
                j = j+1;
                if (j >= (unsigned int)D)
                {
                    j = j-(unsigned int)D;
                }
            }
            score(t, CurrentAtoms, cx, cy, cz, emap);
            geomsaver.Restore(SaveToken);
            t2 = &population2[i];
            if (t.score >= ti->score)
            {
                copy_trial(t2, &t, D);
            }
            else
            {
                copy_trial(t2, ti, D);
            }
        }
        for (i = 0; i < npop; i++)
        {
            copy_trial(&population[i], &population2[i], D);
        }
        sort(population.begin(), population.end(), trial_compare);
        best_r = population[0].score;
        worst_r = population[population.size()-1].score;
        //Logger::log("Generation %d: best=%0.2f worst=%0.2f", igen, best_r, worst_r);
    }

    if (best_r > rdens_start)
    {
        RotateAtomVec(population[0].p[0], population[0].p[1], population[0].p[2], cx, cy, cz, &CurrentAtoms);
        TranslateAtomVec(population[0].p[3], population[0].p[4], population[0].p[5], &CurrentAtoms);
    }
    else
    {
        Logger::log("Final model no better or worse than start - model not moved");
    }

    Logger::log("Start Rdens=%0.2f Final=%0.2f\nRot: %0.2f %0.2f %0.2f Trans: %0.2f %0.2f %0.2f", rdens_start, best_r,
                population[0].p[0], population[0].p[1], population[0].p[2], population[0].p[3], population[0].p[4], population[0].p[5]);

    for (i = 0; i < population.size(); i++)
    {
        delete[] population[i].p;
        delete[] population2[i].p;
    }
}

float collision_score(MIAtomList &atoms, MIAtomList nabors)
{
    float sum = 0;
    unsigned int i, j;
    float dx, dy, dz, d;
    float toocloseCC = 3.5F*3.5F;
    float toocloseNO = 2.7F*2.7F;
    float toocloseself = 2.3F*2.3F;
    float tooclose;
    MIAtom *a1, *a2;
    for (i = 0; i < atoms.size(); i++)
    {
        a1 = atoms[i];
        if (a1->type() & AtomType::SIDECHAIN)
        {
            if (a1->name()[0] != 'C')
            {
                tooclose = toocloseNO;
            }
            else
            {
                tooclose = toocloseCC;
            }
            for (j = 0; j < nabors.size(); j++)
            {
                dx = a1->x() - nabors[j]->x();
                dy = a1->y() - nabors[j]->y();
                dz = a1->z() - nabors[j]->z();
                if ((d = dx*dx+dy*dy+dz*dz) < tooclose)
                {
                    if (d != 0)
                    {
                        sum += tooclose/d;
                    }
                    else
                    {
                        sum += 20.0;
                    }
                }
            }
        }
        // self collisions between main and side
        if (i < atoms.size()-1)
        {
            for (j = i+1; j < atoms.size(); j++)
            {
                a2 = atoms[j];
                if ((a2->type() & AtomType::SIDECHAIN && a1->type() & AtomType::MAINCHAIN)
                    || (a1->type() & AtomType::SIDECHAIN && a2->type() & AtomType::MAINCHAIN))
                {
                    dx = a1->x() - a2->x();
                    dy = a1->y() - a2->y();
                    dz = a1->z() - a2->z();
                    if ((d = dx*dx+dy*dy+dz*dz) < toocloseself)
                    {
                        if (d != 0)
                        {
                            sum += toocloseself/d;
                        }
                        else
                        {
                            sum += 20.0;
                        }
                    }
                }
            }
        }
    }
    return sum;
}

float collision_score_simple(MIAtomList &atoms, MIAtomList nabors)
{
    float sum = 0;
    unsigned int i, j;
    float dx, dy, dz, d;
    float toocloseCC = 3.5F*3.5F;
    float toocloseNO = 2.7F*2.7F;
    float tooclose;
    MIAtom *a1;
    for (i = 0; i < atoms.size(); i++)
    {
        a1 = atoms[i];
        if (a1->type() & AtomType::SIDECHAIN)
        {
            if (a1->name()[0] != 'C')
            {
                tooclose = toocloseNO;
            }
            else
            {
                tooclose = toocloseCC;
            }
            for (j = 0; j < nabors.size(); j++)
            {
                dx = a1->x() - nabors[j]->x();
                dy = a1->y() - nabors[j]->y();
                dz = a1->z() - nabors[j]->z();
                if ((d = dx*dx+dy*dy+dz*dz) < tooclose)
                {
                    if (d != 0)
                    {
                        sum += tooclose/d;
                    }
                    else
                    {
                        sum += 20.0;
                    }
                }
            }
        }
    }
    return sum;
}

static void score_torsion(trial &t, MIAtomList atoms, MIMoleculeBase *model, EMapBase *emap,
                          vector<TORSION> &torsions, unsigned int td,
                          unsigned short *flags, MIAtomList &nabors)
{
    unsigned int i, j, k;
    for (i = 0; i < td; i++)
    {
        for (j = 0; j < atoms.size(); j++)
        {
            k = atoms.size()*i + j;
            atoms[j]->setType(flags[k]);
        }
        model->SetTatom(torsions[i].getAtom2(), torsions[i].atom3);
        while (t.p[i] < 0.0)
        {
            t.p[i] += 360.0F;
        }
        while (t.p[i] >= 360.0F)
        {
            t.p[i] -= 360.0F;
        }
        model->RotateTorsion((float)t.p[i]);
    }
    float collisions = collision_score(atoms, nabors);
    t.score = emap->RDensity(atoms) - collisions;
}

static inline float diff_angle(float angle1, float angle2)
{
    float diff = angle1 - angle2;
    if (diff < -180.0F)
    {
        diff += 360.0F;
    }
    if (diff > 180.0F)
    {
        diff -= 360.0F;
    }
    return diff;
}

void MIMolOpt::TorsionOptimize(MIAtomList &CurrentAtoms, MIMoleculeBase *fitmol,
                               EMapBase *emap, vector<TORSION> &torsions, bool do_setup)
{
    CurrentMap = emap;
    float maxr = 360.0F;
    unsigned int i, j, k;
    trial t;
    unsigned int p1, p2, p3;
    trial *ti, *ta, *tb, *tc, *t2;
    double rdens_start;
    float F = 0.7F;
    float CR = 0.8F;
    unsigned int td = torsions.size();
    unsigned int npop = 20 * td;
    int natoms = CurrentAtoms.size();
    MIAtomList Neighbours;
    MIAtom *a;
    vector<trial> population;
    vector<trial> population2;
    vector<trial>::iterator p;
    std::string s;
    double best_r, worst_r;

    if (!fitmol || !emap)
    {
        return;
    }
    if (!emap->HasDensity() || CurrentAtoms.size() == 0 || torsions.size() == 0)
    {
        return;
    }

    FindNeighbours(CurrentAtoms, Neighbours, fitmol, 10.0F);

    for (i = 0; i < CurrentAtoms.size(); i++)
    {
        a = CurrentAtoms[i];
        a->removeType(AtomType::SIDECHAIN);
        a->removeType(AtomType::MAINCHAIN);
        if (a->name()[0] != 'H' && strcmp(a->name(), "CB") != 0)
        {
            if (MIAtom::MIIsSideChainAtom(a))
            {
                a->addType(AtomType::SIDECHAIN);
            }
            else
            {
                a->addType(AtomType::MAINCHAIN);
            }
        }
    }
    Logger::log("nNeighbours = %d", (int)Neighbours.size());

    unsigned short *tatomflags = new unsigned short[td * natoms];
    for (i = 0; i < (unsigned int)td; i++)
    {
        if (do_setup)
        {
            fitmol->SetupTorsion(torsions[i].getAtom2(), torsions[i].atom3, &CurrentAtoms);
        }
        for (j = 0; j < CurrentAtoms.size(); j++)
        {
            k = CurrentAtoms.size()*i + j;
            tatomflags[k] = CurrentAtoms[j]->type();
        }
    }

    ConnectTo(fitmol); //Connect up signals
    SaveToken = geomsaver.Save(CurrentAtoms, fitmol);

    for (i = 0; i < npop; i++)
    {
        t.p = new double[td];
        if (i == 0)
        {
            for (j = 0; j < (unsigned int)td; j++)
            {
                t.p[j] = 0.0;
            }
        }
        else
        {
            for (j = 0; j < (unsigned int)td; j++)
            {
                if (i != 0)
                {
                    t.p[j] = frand(maxr);
                }
                else
                {
                    t.p[j] = 0.0;
                }
            }
        }
        population.push_back(t);
        t.p = new double[td];
        population2.push_back(t);
    }

    for (p = population.begin(); p != population.end(); p++)
    {
        score_torsion(*p, CurrentAtoms, fitmol, emap, torsions, td, tatomflags, Neighbours);
        geomsaver.Restore(SaveToken);
    }
    rdens_start = population[0].score;

    sort(population.begin(), population.end(), trial_compare);
    best_r = population[0].score;
    worst_r = population[population.size()-1].score;

    for (int igen = 0; igen < 40; igen++)
    {
        for (i = 0; i < npop; i++)
        {
            ti = &population[i];
            do
            {
                p1 = irand(npop);
            } while (p1 == i);
            ta = &population[p1];
            do
            {
                p2 = irand(npop);
            } while (p2 == i || p2 == p1);
            tb = &population[p2];
            do
            {
                p3 = irand(npop);
            } while (p3 == i || p3 == p1 || p3 == p2);
            tc = &population[p3];
            j = irand(td);
            for (k = 1; k <= (unsigned int)td; k++)
            {
                if (frand() < CR || k == (unsigned int)td)
                {
                    t.p[j] = tc->p[j] + F*(diff_angle((float)ta->p[j], (float)tb->p[j]));
                }
                else
                {
                    t.p[j] = ti->p[j];
                }
                j = j+1;
                if (j >= (unsigned int)td)
                {
                    j = j-td;
                }
            }
            score_torsion(t, CurrentAtoms, fitmol, emap, torsions, td, tatomflags, Neighbours);
            geomsaver.Restore(SaveToken);
            t2 = &population2[i];
            if (t.score >= ti->score)
            {
                copy_trial(t2, &t, td);
            }
            else
            {
                copy_trial(t2, ti, td);
            }
        }
        for (i = 0; i < npop; i++)
        {
            copy_trial(&population[i], &population2[i], td);
        }
        sort(population.begin(), population.end(), trial_compare);
        best_r = population[0].score;
        worst_r = population[population.size()-1].score;
        //Logger::log("Generation %d: best=%0.2f worst=%0.2f", igen, best_r, worst_r);
    }

    if (best_r > rdens_start)
    {
        score_torsion(population[0], CurrentAtoms, fitmol, emap, torsions, td, tatomflags, Neighbours);
    }
    else
    {
        Logger::log("Final model no better or worse than start - model not moved");
    }

    Logger::log("Start Score=%0.2f Final=%0.2f\n", rdens_start, best_r);

    for (i = 0; i < population.size(); i++)
    {
        delete[] population[i].p;
        delete[] population2[i].p;
    }
    delete[] tatomflags;
    fitmol->ClearTorsion();

}

static void score_full(trial &t, MIAtomList &atoms, float cx, float cy, float cz, /* center of atoms */
                       MIMoleculeBase *model,
                       vector<TORSION> &torsions, unsigned int td, unsigned short *flags,
                       float sx, float sy, float sz, /* screen center */
                       InterpBox &box)
{
    while (t.p[0] < 0.0)
    {
        t.p[0] += 360.0F;
    }
    while (t.p[0] >= 360.0F)
    {
        t.p[0] -= 360.0F;
    }
    while (t.p[1] < 0.0)
    {
        t.p[1] += 360.0F;
    }
    while (t.p[1] >= 360.0F)
    {
        t.p[1] -= 360.0F;
    }
    while (t.p[2] < 0.0)
    {
        t.p[2] += 360.0F;
    }
    while (t.p[2] >= 360.0F)
    {
        t.p[2] -= 360.0F;
    }
    RotateAtomVec(t.p[0], t.p[1], t.p[2], cx, cy, cz, &atoms);

    TranslateAtomVec(t.p[3], t.p[4], t.p[5], &atoms);

    // Apply the torsions
    if (td > 6)
    {
        unsigned int i, j, k, it = 0;
        for (i = 6; i < td; i++)
        {
            for (j = 0; j < atoms.size(); j++)
            {
                k = atoms.size()*it + j;
                atoms[j]->setType(flags[k]);
            }
            model->SetTatom(torsions[it].getAtom2(), torsions[it].atom3);
            while (t.p[i] < 0.0)
            {
                t.p[i] += 360.0F;
            }
            while (t.p[i] >= 360.0F)
            {
                t.p[i] -= 360.0F;
            }
            model->RotateTorsion((float)t.p[i]);
            it++;
        }
    }
    // keep the model from drifting off into neverland
    float tx, ty, tz;
    // penalize moving too far from center
    GetCenter(atoms, tx, ty, tz);
    tx = tx - sx;
    ty = ty - sy;
    tz = tz - sz;
    float dmoved = tx*tx + ty*ty + tz*tz;
    // allows things to move 5 A or so without penalty;
    //	std::string mov;
    //	mov.Printf("%6.2f", dmoved);
    //	Logger::log(mov);

    dmoved -= 125.F;
    if (dmoved < 0.0)
    {
        dmoved = 0.0;
    }

    t.score = box.RDensity(atoms) - dmoved/3.0F;
}

float BumpScore(vector<Bond> &bumps)
{
    float score = 0.0, d;
    for (size_t i = 0; i < bumps.size(); i++)
    {
        d = (float)AtomDist(*bumps[i].getAtom1(), *bumps[i].getAtom2());
        if (d < bumps[i].ideal_length)
        {
            d = bumps[i].ideal_length - d;
            score += d;
            if (d < bumps[i].ideal_length/2.0F)
            {
                score += 10.0F;
            }
        }
    }
    return score;
}

void MIMolOpt::FullOptimize(MIAtomList &CurrentAtoms, MIMoleculeBase *fitmol, EMapBase *emap, const float *center,
                            InterpBox &box, unsigned int refine_level,
                            MIMolOptCheckPoint *checkpoint)
{

    // Full optimize - optimize all angles, translations and torsions
    // i.e place a ligand in density.  Assumes ligand overlaps wanted density.
    CurrentMap = emap;
    float maxt = 3.0F;
    float maxr = 180.0F;
    float cx = 0, cy = 0, cz = 0;
    float dx, dy, dz;
    unsigned int i, j, k;
    trial t;

    unsigned int p1, p2, p3;
    trial *ti, *ta, *tb, *tc, *t2;


    // These parameters are the key to optimizing the DE Solver
    float F = 0.8F;
    float CR = 0.9F;
    const int D = 6;
    float good_enough = 3.0F;
    unsigned int maxgen = 400;
    unsigned int n_per_param = 300, n_per_torsion = 360;

    unsigned int npop;
    unsigned int itrial, maxtrials = 10, nconftries = 0;
    float conf_angle = 180.0;
    vector<TORSION>  torsions;
    Residue *res;
    unsigned short *tatomflags;
    GeomSaver best_solution;
    float screen_center_x = center[0];
    float screen_center_y = center[1];
    float screen_center_z = center[2];
    unsigned int SaveToken, ConformerToken;
    vector<unsigned int>SaveTokens;
    vector<trial> population;
    vector<trial> population2;
    vector<trial>::iterator p;
    vector<float> var_start;
    double best_r, worst_r, best_r_yet, rdens_start;
    vector<Bond> bumps;
    bool flag = 0;
    char buf[2000];

    //FILE	*fpout_ptr;
    //FILE	*fpin_ptr;

    if (!fitmol || !emap)
    {
        return;
    }
    if (!emap->HasDensity() || CurrentAtoms.size() == 0)
    {
        return;
    }

    // I/O for testing purposes
    /*
       fpin_ptr = fopen("C:\\MIFit_Src\\sti_example\\parameters.txt","r");
       if (fpin_ptr == NULL)
       {
       printf("\nCannot open input file\n");
       exit(1);
       }
       fscanf(fpin_ptr,"%f,%f,%d,%d",&F,&CR,&n_per_param,&maxgen);
       fclose(fpin_ptr);
     */
    //fpout_ptr = fopen("DEtest.txt","a");


    //find the residue associated with CurrentAtoms
    // ( might want to add ability to search for more than one residue in future)
    res = residue_from_atom(fitmol->getResidues(), CurrentAtoms[0]);
    if (res != NULL && refine_level != Refine_Level::None)
    {
        // get the torsions associated with these atoms
        dict.GetResidueTorsions(res, torsions);
    }

    // if there are torsions we need to worry about internal bumps
    if (torsions.size() > 0)
    {
        dict.BuildInternalBumpBonds(CurrentAtoms, bumps);
    }

    //the number of trials per torsion depends on the refine level;
    std::string what("Full Search:");
    if (refine_level == Refine_Level::Quick)
    {
        n_per_torsion = 10;
        n_per_param = 20;
        maxgen = 400;
        good_enough = 4.7F;
        what = "Quick Search:";
    }
    else if (refine_level == Refine_Level::Optimize)
    {
        n_per_torsion = 10;
        n_per_param = 20;
        maxgen = 500;
        conf_angle = 45.0F;
        maxr = 5.0F;
        maxt = 0.2F;
        good_enough = 6.0F;
        //		good_enough = 60.0F;
        what = "Optimizing:";
    }
    else
    {
        n_per_torsion = 30;
        //		n_per_torsion = 100;
        good_enough = 8.0F;
        //		good_enough = 6.0F;   JB increased from 7.3 to 8.0
        n_per_param = 40;
        maxgen = 150;
        //		maxgen = 200;  JB decreased back from 700 to 150
        what = "Full Search:";
    }
    Logger::log("Found %d atoms and %d torsions to be optimized",
                (int)CurrentAtoms.size(), (int)torsions.size());

    ConformerToken = SaveToken = best_solution.Save(CurrentAtoms, fitmol);
    SaveTokens.push_back(SaveToken);


    Logger::footer("Setting up torsions...");
    // find the atoms to be torsioned and save their flags for each torsion
    if (torsions.size() > 0)
    {
        tatomflags = new unsigned short[(torsions.size())*CurrentAtoms.size()];
        for (i = 0; i < (unsigned int)(torsions.size()); i++)
        {
            fitmol->SetupTorsion(torsions[i].getAtom2(), torsions[i].atom3, &CurrentAtoms);
            for (j = 0; j < CurrentAtoms.size(); j++)
            {
                k = CurrentAtoms.size()*i + j;
                tatomflags[k] = CurrentAtoms[j]->type();
            }
        }
    }

    for (j = 0; j < 3; j++)
    {
        var_start.push_back(maxr);
    }
    // translations +/- 2 A
    for (j = 3; j < 6; j++)
    {
        var_start.push_back(maxt);
    }

    npop = n_per_param*D;
    maxtrials = (torsions.size())*n_per_torsion;
    sprintf(buf, "Number of Generations = %d, Population size = %d, Maximum Trials = %d\n\n", maxgen, npop, maxtrials);
    Logger::log(buf);
    var_start.reserve(D);
    population.reserve(npop);
    population2.reserve(npop);
    for (itrial = 0; itrial < maxtrials; itrial++)
    {
        //if not the first trial, generate a new conformation
        if (itrial > 0)
        {
            unsigned int j, k, it;
            nconftries = 0;
            // look for a random conformation w/o bad bumps
            do
            {
                for (it = 0; it < torsions.size(); it++)
                {
                    for (j = 0; j < CurrentAtoms.size(); j++)
                    {
                        k = CurrentAtoms.size()*it + j;
                        CurrentAtoms[j]->setType(tatomflags[k]);
                    }
                    fitmol->SetTatom(torsions[it].getAtom2(), torsions[it].atom3);
                    fitmol->RotateTorsion(frand2(conf_angle));
                }
                // this counter is used to prevent an infinite loop
                nconftries++;
            } while (nconftries < 100 && BumpScore(bumps) > 10.0F);
        }
        // ideally we will rarely see this message
        if (nconftries >= 100)
        {
            Logger::log("Warning: max conformer tries reached");
        }
        // find the center of mass of the atoms
        GetCenter(CurrentAtoms, cx, cy, cz);
        // move to center
        dx = cx - screen_center_x;
        dy = cy - screen_center_y;
        dz = cz - screen_center_z;
        for (j = 0; j < CurrentAtoms.size(); j++)
        {
            CurrentAtoms[j]->translate(-dx, -dy, -dz);
        }
        ConformerToken = best_solution.Save(CurrentAtoms, fitmol);

        // find the best rotation-translation for the conformation
        for (i = 0; i < npop; i++)
        {
            t.p = new double[D];
            if (i == 0 && itrial == 0)
            {
                for (j = 0; j < (unsigned int)D; j++)
                {
                    t.p[j] = 0.0;
                }
            }
            else
            {
                for (j = 0; j < (unsigned int)D; j++)
                {
                    t.p[j] = frand2(var_start[j]);
                }
            }
            population.push_back(t);
            t.p = new double[D];
            population2.push_back(t);
        }

        for (p = population.begin(); p != population.end(); p++)
        {
            score_full(*p, CurrentAtoms, cx, cy, cz, fitmol,
                       torsions, D, tatomflags, screen_center_x, screen_center_y,
                       screen_center_z, box);

            best_solution.Restore(ConformerToken);
        }
        if (itrial == 0)
        {
            rdens_start = best_r_yet = population[0].score;
        }

        sort(population.begin(), population.end(), trial_compare);
        best_r = population[0].score;
        worst_r = population[population.size()-1].score;

        for (unsigned int igen = 0; igen < maxgen; igen++)
        {
            for (i = 0; i < npop; i++)
            {
                ti = &population[i];
                do
                {
                    p1 = irand(npop);
                } while (p1 == i);
                ta = &population[p1];
                do
                {
                    p2 = irand(npop);
                } while (p2 == i || p2 == p1);
                tb = &population[p2];
                do
                {
                    p3 = irand(npop);
                } while (p3 == i || p3 == p1 || p3 == p2);
                tc = &population[p3];
                //The following code introduces different DE strategies
                //Added by Daniel Pick, Jan. 2005
                /*
                   //strategy Rand1Bin
                   j = irand(D);
                   for(k=1; k<=(unsigned int)D; k++){
                    if(frand() < CR || k==(unsigned int)D){
                        t.p[j]= tc->p[j] + F*(ta->p[j]-tb->p[j]);
                    } else {
                        t.p[j]= ti->p[j];
                    }
                    j=j+1;
                    if(j>=(unsigned int)D)j=j-(unsigned int)D;
                   }
                 */
                //strategy Rand1Exp
                j = irand(D);
                flag = 0;
                for (k = 1; k <= (unsigned int)D; k++)
                {

                    if (frand() < CR || k == (unsigned int)D)
                    {
                        flag = 1;
                    }
                    if (flag == 1)
                    {
                        t.p[j] = tc->p[j] + F*(ta->p[j]-tb->p[j]);
                    }
                    else
                    {
                        t.p[j] = ti->p[j];
                    }
                    j = (j+1)%D;
                    //if(j>=(unsigned int)D)j=j-(unsigned int)D;
                }
                /*
                   //strategy RandtoBest1Exp
                   j = irand(D);
                   flag = 0;
                   for(k=1; k<=(unsigned int)D; k++){
                    if(frand() < CR || k==(unsigned int)D) flag = 1;
                    if (flag == 1) {
                        tc = &population[0];
                        t.p[j] += F * (tc->p[j] - t.p[j]) + F*(ta->p[j]-tb->p[j]);
                    } else {
                        t.p[j]= ti->p[j];
                    }
                    j= (j+1)%D;
                    //if(j>=(unsigned int)D)j=j-(unsigned int)D;
                   }
                 */
                score_full(t, CurrentAtoms, cx, cy, cz, fitmol,
                           torsions, D, tatomflags, screen_center_x, screen_center_y,
                           screen_center_z, box);
                best_solution.Restore(ConformerToken);
                t2 = &population2[i];
                if (t.score >= ti->score)
                {
                    copy_trial(t2, &t, D);
                }
                else
                {
                    copy_trial(t2, ti, D);
                }
            }
            for (i = 0; i < npop; i++)
            {
                copy_trial(&population[i], &population2[i], D);
            }
            sort(population.begin(), population.end(), trial_compare);
            best_r = population[0].score;
            worst_r = population[population.size()-1].score;

            if (igen%20 == 0)
            {
                sprintf(buf, "Generation %d: best=%0.2f worst=%0.2f", igen, best_r, worst_r);
                Logger::log(buf);
            }
            if (best_r >= good_enough)
            {
                Logger::log("Good enough after Generation %d: best=%0.2f", igen, best_r);
                break;
            }
            //float percent = 100.0F*((float)igen/(float)maxgen);
            // replace the lowest 10% with random every 10th generation

            // This if statement was commented out in original code
            //if(igen%10==0 && igen > 0){
            //	for(i=ROUND(0.9*npop); i<npop; i++){
            //		for(j=0;j<(unsigned int)D;j++) {
            //			population[i].p[j]=frand2(var_start[j]);
            //		}
            //	}
            //}
        }

        Logger::footer("%s Trial %d of %d: %d percent done; best score = %0.1f; best yet %0.1f", what.c_str(),
                       itrial+1, maxtrials, ROUND(100.0F*((float)(itrial+1)/(float)maxtrials)), best_r, best_r_yet);
        //fprintf(fpout_ptr,"%s Trial %d of %d: %d percent done; best score = %0.1f; best yet %0.1f \n\n", what.c_str(),
        //	itrial+1, maxtrials, ROUND(100.0F*((float)(itrial+1)/(float)maxtrials)), best_r, best_r_yet);
        Logger::log(buf);

        if (best_r > best_r_yet)
        {
            score_full(population[0], CurrentAtoms, cx, cy, cz, fitmol,
                       torsions, D, tatomflags, screen_center_x, screen_center_y,
                       screen_center_z, box);
            if (checkpoint)
            {
                (*checkpoint)(fitmol);
            }
            best_r_yet = best_r;
            SaveToken = best_solution.Save(CurrentAtoms, fitmol);
            SaveTokens.push_back(SaveToken);
        }
        best_solution.Restore(SaveTokens[0]);
        //s.Printf("Trial %d: Best score=%0.2f\nRot: %0.2f %0.2f %0.2f Trans: %0.2f %0.2f %0.2f", itrial+1, best_r,
        //	population[0].p[0], population[0].p[1], population[0].p[2],population[0].p[3], population[0].p[4], population[0].p[5]);
        //Logger::log(s);
        for (i = 0; i < population.size(); i++)
        {
            delete[] population[i].p;
            delete[] population2[i].p;
        }
        population.clear();
        population2.clear();
        // Original code
        if (best_r >= good_enough)
        {
            break;
        }
    }


    if (best_r_yet > rdens_start)
    {
        SaveToken = SaveTokens[SaveTokens.size()-1];
        best_solution.Restore(SaveToken);
        Logger::log("Final score = %f", (float)best_r_yet);
    }
    else
    {
        Logger::log("Final model no better or worse than start - model not moved");
    }


    if (torsions.size() > 0)
    {
        delete[] tatomflags;
    }
    Logger::footer("");
}

void MIMolOpt::LigandOptimize(MIAtomList &CurrentAtoms, MIMoleculeBase *fitmol, EMapBase *emap, const float *center,
                              InterpBox &box, unsigned int refine_level, GeomSaver &confs,
                              MIMolOptCheckPoint *checkpoint)
{

    // Full optimize - optimize all angles, translations and torsions
    // i.e place a ligand in density.  Assumes ligand overlaps wanted density.
    CurrentMap = emap;
    float maxt = 3.0F;
    float maxr = 180.0F;
    float cx = 0, cy = 0, cz = 0;
    float dx, dy, dz;
    unsigned int i, j, k;
    trial t;

    unsigned int p1, p2, p3;
    trial *ti, *ta, *tb, *tc, *t2;


    // These parameters are the key to optimizing the DE Solver
    float F = 0.1F;
    float CR = 0.1F;
    //float F=0.8F; original mutation factor
    //float CR=0.9F; original cross-over factor
    const int D = 6;
    float good_enough = 3.0F;
    unsigned int maxgen = 400;
    unsigned int n_per_param = 300, n_per_torsion = 360;

    unsigned int npop;
    unsigned int itrial, maxtrials = 10;
    float conf_angle = 180.0;
    vector<TORSION>  torsions;                  //Not used any more, but score_full() needs a placeholder
    Residue *res;
    unsigned short *tatomflags = NULL;
    GeomSaver best_solution;
    float screen_center_x = center[0];
    float screen_center_y = center[1];
    float screen_center_z = center[2];
    unsigned int SaveToken, ConformerToken;
    vector<unsigned int>SaveTokens;
    vector<trial> population;
    vector<trial> population2;
    vector<trial>::iterator p;
    vector<float> var_start;
    double best_r, worst_r, rdens_start;
    double best_r_yet = -DBL_MAX;
    //	vector<Bond> bumps;
    bool flag = 0;
    char buf[2000];

    //FILE	*fpout_ptr;
    //FILE	*fpin_ptr;

    if (!fitmol || !emap)
    {
        return;
    }
    if (!emap->HasDensity() || CurrentAtoms.size() == 0)
    {
        return;
    }

    // I/O for testing purposes
    /*
       fpin_ptr = fopen("C:\\MIFit_Src\\sti_example\\parameters.txt","r");
       if (fpin_ptr == NULL)
       {
       printf("\nCannot open input file\n");
       exit(1);
       }
       fscanf(fpin_ptr,"%f,%f,%d,%d",&F,&CR,&n_per_param,&maxgen);
       fclose(fpin_ptr);
     */
    //fpout_ptr = fopen("DEtest.txt","a");


    //find the residue associated with CurrentAtoms
    // ( might want to add ability to search for more than one residue in future)
    res = residue_from_atom(fitmol->getResidues(), CurrentAtoms[0]);
    /*
        if(res !=NULL && refine_level != Refine_Level::None){
            // get the torsions associated with these atoms
            GetResidueTorsions(res, torsions);
        }
     */

    //Let the geomrefiner know that we're going to refine
    dict.Clear();
    //Store state of constraint prefs to restore and
    bool tmp_ca = dict.GetConstrainCA();
    bool tmp_ends = dict.GetConstrainEnds();

    //Don't add constraints to refinement in the editor
    dict.SetConstrainCA(false);
    dict.SetConstrainEnds(false);
    SetRefiRes(res, res, fitmol, emap);

    //Restore state
    dict.SetConstrainCA(tmp_ca);
    dict.SetConstrainEnds(tmp_ends);
    /*
        // if there are torsions we need to worry about internal bumps
        int nt = torsions.size();
        if(torsions.size() > 0){
            BuildInternalBumpBonds(CurrentAtoms, bumps);
        }
     */
    dict.GetFlexibleTorsions(torsions, res);
    if (confs.NumberSets() <= 2 && !torsions.empty())
    {
        conflib::GenerateEnsemble(res,
                                  fitmol,
                                  dict.RefiBonds,
                                  torsions,
                                  confs);
    }

    //the number of trials per torsion depends on the refine level;
    std::string what("Full Search:");
    if (refine_level == Refine_Level::Quick)
    {
        n_per_torsion = 10;
        n_per_param = 20;
        maxgen = 400;
        good_enough = 4.7F;
        what = "Quick Search:";
    }
    else if (refine_level == Refine_Level::Optimize)
    {
        n_per_torsion = 10;
        n_per_param = 20;
        maxgen = 500;
        conf_angle = 45.0F;
        maxr = 5.0F;
        maxt = 0.2F;
        good_enough = 6.0F;
        //		good_enough = 60.0F;
        what = "Optimizing:";
    }
    else
    {
        n_per_torsion = 40;
        //		n_per_torsion = 100; JBtest increase from 30 to 100
        good_enough = 1000.0F;
        //		good_enough = 6.0F; JB increased to 8.0 JBtest disable, 1000
        n_per_param = 40;
        maxgen = 100;
        //		maxgen = 200; JB decreased back to 150 JBtest 100
        what = "Full Search:";
    }
    Logger::log("Found %d atoms and %d torsions to be optimized",
                //		(int)CurrentAtoms.size(), (int)torsions.size());
                (int)CurrentAtoms.size(), (int)torsions.size());

    ConformerToken = SaveToken = best_solution.Save(CurrentAtoms, fitmol);
    SaveTokens.push_back(SaveToken);

    /*
        Logger::footer("Setting up torsions...");
        // find the atoms to be torsioned and save their flags for each torsion
        if(torsions.size() > 0){
            tatomflags = new unsigned short[(torsions.size())*CurrentAtoms.size()];
            for(i=0; i<(unsigned int)(torsions.size()); i++){
                fitmol->SetupTorsion(torsions[i].getAtom2(), torsions[i].atom3, &CurrentAtoms);
                for(j=0;j<CurrentAtoms.size();j++){
                    k = CurrentAtoms.size()*i + j;
                    tatomflags[k] = CurrentAtoms[j]->type;
                }
            }
        }
     */
    for (j = 0; j < 3; j++)
    {
        var_start.push_back(maxr);
    }
    // translations +/- 2 A
    for (j = 3; j < 6; j++)
    {
        var_start.push_back(maxt);
    }


    //Hack a "trial" with no rotation or translation, to simply
    //score an existing structure

    trial noMove;
    noMove.p = new double[D];
    for (int xDim = 0; xDim < D; xDim++)
    {
        noMove.p[xDim] = 0.0;
    }

    int c;
    npop = n_per_param*D;
    maxtrials = (torsions.size())*n_per_torsion;
    sprintf(buf, "Number of Generations = %d, Population size = %d, Maximum Trials = %d\n\n", maxgen, npop, maxtrials);
    Logger::log(buf);
    var_start.reserve(D);
    population.reserve(npop);
    population2.reserve(npop);
    for (itrial = 0; itrial < maxtrials; itrial++)
    {

        // Choose a random conformation from the set in "confs"
        c = chemlib::irand_approx(confs.NumberSets() - 1) + 1;
        confs.Restore(c);


        // find the center of mass of the atoms
        GetCenter(CurrentAtoms, cx, cy, cz);
        // move to center
        dx = cx - screen_center_x;
        dy = cy - screen_center_y;
        dz = cz - screen_center_z;
        for (j = 0; j < CurrentAtoms.size(); j++)
        {
            CurrentAtoms[j]->translate(-dx, -dy, -dz);
        }
        ConformerToken = best_solution.Save(CurrentAtoms, fitmol);

        // find the best rotation-translation for the conformation
        for (i = 0; i < npop; i++)
        {
            t.p = new double[D];
            if (i == 0 && itrial == 0)
            {
                for (j = 0; j < (unsigned int)D; j++)
                {
                    t.p[j] = 0.0;
                }
            }
            else
            {
                for (j = 0; j < (unsigned int)D; j++)
                {
                    t.p[j] = frand2(var_start[j]);
                }
            }
            population.push_back(t);
            t.p = new double[D];
            population2.push_back(t);
        }

        for (p = population.begin(); p != population.end(); p++)
        {
            score_full(*p, CurrentAtoms, screen_center_x, screen_center_y, screen_center_z, fitmol,
                       torsions, D, tatomflags, screen_center_x, screen_center_y,
                       screen_center_z, box);

            //			Refine();
            best_solution.Restore(ConformerToken);
            //			score_full( *p, CurrentAtoms, screen_center_x, screen_center_y, screen_center_z, fitmol,
            //				torsions, D, tatomflags, screen_center_x, screen_center_y,
            //				screen_center_z, box);
        }
        if (itrial == 0)
        {
            rdens_start = best_r_yet = population[0].score;
        }

        sort(population.begin(), population.end(), trial_compare);
        best_r = population[0].score;
        worst_r = population[population.size()-1].score;

        for (unsigned int igen = 0; igen < maxgen; igen++)
        {
            for (i = 0; i < npop; i++)
            {
                ti = &population[i];
                do
                {
                    p1 = irand(npop);
                } while (p1 == i);
                ta = &population[p1];
                do
                {
                    p2 = irand(npop);
                } while (p2 == i || p2 == p1);
                tb = &population[p2];
                do
                {
                    p3 = irand(npop);
                } while (p3 == i || p3 == p1 || p3 == p2);
                tc = &population[p3];
                //The following code introduces different DE strategies
                //Added by Daniel Pick, Jan. 2005
                /*
                   //strategy Rand1Bin
                   j = irand(D);
                   for(k=1; k<=(unsigned int)D; k++){
                    if(frand() < CR || k==(unsigned int)D){
                        t.p[j]= tc->p[j] + F*(ta->p[j]-tb->p[j]);
                    } else {
                        t.p[j]= ti->p[j];
                    }
                    j=j+1;
                    if(j>=(unsigned int)D)j=j-(unsigned int)D;
                   }
                 */
                //strategy Rand1Exp
                j = irand(D);
                flag = 0;
                for (k = 1; k <= (unsigned int)D; k++)
                {

                    if (frand() < CR || k == (unsigned int)D)
                    {
                        flag = 1;
                    }
                    if (flag == 1)
                    {
                        t.p[j] = tc->p[j] + F*(ta->p[j]-tb->p[j]);
                    }
                    else
                    {
                        t.p[j] = ti->p[j];
                    }
                    j = (j+1)%D;
                    //if(j>=(unsigned int)D)j=j-(unsigned int)D;
                }
                /*
                   //strategy RandtoBest1Exp
                   j = irand(D);
                   flag = 0;
                   for(k=1; k<=(unsigned int)D; k++){
                    if(frand() < CR || k==(unsigned int)D) flag = 1;
                    if (flag == 1) {
                        tc = &population[0];
                        t.p[j] += F * (tc->p[j] - t.p[j]) + F*(ta->p[j]-tb->p[j]);
                    } else {
                        t.p[j]= ti->p[j];
                    }
                    j= (j+1)%D;
                    //if(j>=(unsigned int)D)j=j-(unsigned int)D;
                   }
                 */
                score_full(t, CurrentAtoms, screen_center_x, screen_center_y, screen_center_z, fitmol,
                           torsions, D, tatomflags, screen_center_x, screen_center_y,
                           screen_center_z, box);
                //				Refine();
                //				score_full( t, CurrentAtoms, screen_center_x, screen_center_y, screen_center_z, fitmol,
                //					torsions, D, tatomflags, screen_center_x, screen_center_y,
                //					screen_center_z, box);
                best_solution.Restore(ConformerToken);
                t2 = &population2[i];
                if (t.score >= ti->score)
                {
                    copy_trial(t2, &t, D);
                }
                else
                {
                    copy_trial(t2, ti, D);
                }
            }
            for (i = 0; i < npop; i++)
            {
                copy_trial(&population[i], &population2[i], D);
            }
            sort(population.begin(), population.end(), trial_compare);
            best_r = population[0].score;
            worst_r = population[population.size()-1].score;

            //			if(igen%20==0){
            //				sprintf(buf,"Generation %d: best=%0.2f worst=%0.2f", igen, best_r, worst_r);
            //				Logger::log(buf);
            //			}
            if (best_r >= good_enough)
            {
                Logger::log("Good enough after Generation %d: best=%0.2f", igen, best_r);
                break;
            }
            //float percent = 100.0F*((float)igen/(float)maxgen);
            // replace the lowest 10% with random every 10th generation

            // This if statement was commented out in original code
            //if(igen%10==0 && igen > 0){
            //	for(i=ROUND(0.9*npop); i<npop; i++){
            //		for(j=0;j<(unsigned int)D;j++) {
            //			population[i].p[j]=frand2(var_start[j]);
            //		}
            //	}
            //}
        } //End of loop over generations (igen)
        score_full(population[0], CurrentAtoms, screen_center_x, screen_center_y, screen_center_z, fitmol,
                   torsions, D, tatomflags, screen_center_x, screen_center_y,
                   screen_center_z, box);
        Refine();
        score_full(noMove, CurrentAtoms, screen_center_x, screen_center_y, screen_center_z, fitmol,
                   torsions, D, tatomflags, screen_center_x, screen_center_y,
                   screen_center_z, box);

        //We don't assume that the refinement improved the score
        if (best_r > noMove.score)
        {
            best_solution.Restore(ConformerToken);
            score_full(population[0], CurrentAtoms, screen_center_x, screen_center_y, screen_center_z, fitmol,
                       torsions, D, tatomflags, screen_center_x, screen_center_y,
                       screen_center_z, box);
        }
        else
        {
            best_r = noMove.score;
        }



        Logger::footer("%s Trial %d of %d: %d percent done; best score = %0.1f; best yet %0.1f", what.c_str(),
                       itrial+1, maxtrials, ROUND(100.0F*((float)(itrial+1)/(float)maxtrials)), best_r, best_r_yet);
        Logger::log("%s Trial %d of %d: conformation %d best score = %0.2f", what.c_str(),
                    itrial+1, maxtrials, c, best_r);
        //fprintf(fpout_ptr,"%s Trial %d of %d: %d percent done; best score = %0.1f; best yet %0.1f \n\n", what.c_str(),
        //	itrial+1, maxtrials, ROUND(100.0F*((float)(itrial+1)/(float)maxtrials)), best_r, best_r_yet);
        Logger::log(buf);

        if (best_r > best_r_yet)
        {
            //			score_full( population[0], CurrentAtoms, screen_center_x, screen_center_y, screen_center_z, fitmol,
            //					torsions, D, tatomflags, screen_center_x, screen_center_y,
            //					screen_center_z, box);
            if (checkpoint)
            {
                (*checkpoint)(fitmol);
            }
            best_r_yet = best_r;
            SaveToken = best_solution.Save(CurrentAtoms, fitmol);
            SaveTokens.push_back(SaveToken);
        }
        //		best_solution.Restore(SaveTokens[0]);
        best_solution.Restore(ConformerToken);

        //s.Printf("Trial %d: Best score=%0.2f\nRot: %0.2f %0.2f %0.2f Trans: %0.2f %0.2f %0.2f", itrial+1, best_r,
        //	population[0].p[0], population[0].p[1], population[0].p[2],population[0].p[3], population[0].p[4], population[0].p[5]);
        //Logger::log(s);
        for (i = 0; i < population.size(); i++)
        {
            delete[] population[i].p;
            delete[] population2[i].p;
        }
        population.clear();
        population2.clear();
        // Original code
        if (best_r >= good_enough)
        {
            break;
        }
    } //End of loop over trials


    if (best_r_yet > rdens_start)
    {
        SaveToken = SaveTokens[SaveTokens.size()-1];
        best_solution.Restore(SaveToken);
        Logger::log("Final score = %f", (float)best_r_yet);
    }
    else
    {
        Logger::log("Final model no better or worse than start - model not moved");
    }

    //When I added refinement to the ligand fitting, this became necessary to
    //end the refine.
    if (IsRefining())
    {
        geomsaver.RestoreColor(this->SaveToken, AtomType::REFIATOM);
        internalSetRefiRes(NULL, 0);
    }

    Logger::footer("");
}

bool MIMolOpt::BuildMainchain(Residue *res, MIMoleculeBase *model, EMapBase *emap, bool addAtomsToNextResidue)
{
    // build a peptide plane between two residues, the input and the next one.
    // and then find the best rotation along the CA-CA axis according to the map
    if (!res)
    {
        return false;
    }
    Residue *next = res->next();
    if (!next)
    {
        return false;
    }
    MIAtom *CA = MIAtomFromNameIncludingSynonyms("CA", res);
    MIAtom *CAnext = MIAtomFromNameIncludingSynonyms("CA", next);
    if (!CA || !CAnext)
    {
        return false;
    }
    int i;

    float N1[]  = {18.749F,  14.921F,  11.009F};
    float CA1[]  = {18.476F,  13.539F,  11.387F};
    float C1[]  = {18.653F,  13.366F,  12.905F};
    float O1[]  = {19.643F,  13.901F,  13.431F};
    float N2[]  = {17.785F, 12.699F, 13.672F};
    float CA2[]  = {18.053F, 12.491F, 15.093F};
    float C2[]  = {18.776F, 11.149F, 15.201F};
    float O2[]  = {18.347F,  10.169F,  14.571F};
    float dv[3];
    dv[0] = CA1[0] - CA->x();
    dv[1] = CA1[1] - CA->y();
    dv[2] = CA1[2] - CA->z();

    // move model peptide on to CA1
    for (i = 0; i < 3; i++)
    {
        N1[i] -= dv[i];
        CA1[i] -= dv[i];
        C1[i] -= dv[i];
        O1[i] -= dv[i];
        N2[i] -= dv[i];
        CA2[i] -= dv[i];
        C2[i] -= dv[i];
        O2[i] -= dv[i];
    }

    // v1 = CA->CAnext and v2 is CA1->CA2
    float v1[3], v2[3];
    v1[0] = CAnext->x() - CA->x();
    v1[1] = CAnext->y() - CA->y();
    v1[2] = CAnext->z() - CA->z();
    for (i = 0; i < 3; i++)
    {
        v2[i] = CA2[i] - CA1[i];
    }
    // get the cross product = we will rotate about this to bring vectors coincident
    float c[3];
    cross(v1, v2, c);

    // angle between v1 and v2 (in radians)
    float angle = vectorangle(v1, v2);
    float degtor = (float)acos(-1.0)/180.0F;
    angle /= degtor;

    float mat[4][3];

    initrotate(CA->x(), CA->y(), CA->z(), c[0], c[1], c[2], -angle, mat);

    rotate(N1, mat);
    rotate(C1, mat);
    rotate(O1, mat);
    rotate(CA2, mat);
    rotate(N2, mat);
    rotate(C2, mat);
    rotate(O2, mat);

    MIAtomList atoms;

    MIAtom *aC1 = MIAtomFromNameIncludingSynonyms("C", res);
    MIAtom *aO1 = MIAtomFromNameIncludingSynonyms("O", res);
    MIAtom *aN2 = MIAtomFromNameIncludingSynonyms("N", next);
    MIAtom *aN1 = MIAtomFromNameIncludingSynonyms("N", res);
    MIAtom *aO2 = MIAtomFromNameIncludingSynonyms("O", next);
    MIAtom *aC2 = MIAtomFromNameIncludingSynonyms("C", next);
    if (!aC1)
    {
        aC1 = new MIAtom;
        aC1->copyShallow(*CA);
        aC1->setName("C");
        aC1->setType(MIAtom::MIGetAtomTypeFromName(aC1->name()));
        if (MIGetColorSetter())
        {
            (*MIGetColorSetter())(aC1);
        }
        res->addAtom(aC1);
    }
    aC1->setPosition(C1[0], C1[1], C1[2]);
    if (!aN1)
    {
        aN1 = new MIAtom;
        aN1->copyShallow(*CA);
        aN1->setName("N");
        aN1->setType(MIAtom::MIGetAtomTypeFromName(aN1->name()));
        if (MIGetColorSetter())
        {
            (*MIGetColorSetter())(aN1);
        }
        res->addAtom(aN1);
        aN1->setPosition(N1[0], N1[1], N1[2]);
    }

    if (!aO1)
    {
        aO1 = new MIAtom;
        aO1->copyShallow(*CA);
        aO1->setName("O");
        aO1->setType(MIAtom::MIGetAtomTypeFromName(aO1->name()));
        if (MIGetColorSetter())
        {
            (*MIGetColorSetter())(aO1);
        }
        res->addAtom(aO1);
    }
    aO1->setPosition(O1[0], O1[1], O1[2]);
    if (!aN2)
    {
        aN2 = new MIAtom;
        aN2->copyShallow(*CAnext);
        aN2->setName("N");
        aN2->setType(MIAtom::MIGetAtomTypeFromName(aN2->name()));
        if (MIGetColorSetter())
        {
            (*MIGetColorSetter())(aN2);
        }
        if (addAtomsToNextResidue)
        {
            next->addAtom(aN2);
        }
    }
    aN2->setPosition(N2[0], N2[1], N2[2]);
    if (!aO2)
    {
        aO2 = new MIAtom;
        aO2->copyShallow(*CAnext);
        aO2->setName("O");
        aO2->setType(MIAtom::MIGetAtomTypeFromName(aO2->name()));
        if (MIGetColorSetter())
        {
            (*MIGetColorSetter())(aO2);
        }
        if (addAtomsToNextResidue)
        {
            next->addAtom(aO2);
        }
        aO2->setPosition(O2[0], O2[1], O2[2]);
    }
    if (!aC2)
    {
        aC2 = new MIAtom;
        aC2->copyShallow(*CAnext);
        aC2->setName("C");
        aC2->setType(MIAtom::MIGetAtomTypeFromName(aC2->name()));
        if (MIGetColorSetter())
        {
            (*MIGetColorSetter())(aC2);
        }
        if (addAtomsToNextResidue)
        {
            next->addAtom(aC2);
        }
        aC2->setPosition(C2[0], C2[1], C2[2]);
    }

    aC1->addType(AtomType::TORSIONATOM);
    aO1->addType(AtomType::TORSIONATOM);
    aN2->addType(AtomType::TORSIONATOM);


    atoms.push_back(CA);
    atoms.push_back(aC1);
    atoms.push_back(aO1);
    atoms.push_back(aN2);
    atoms.push_back(CAnext);
    vector<TORSION> torsions;
    TORSION torsion;
    torsion.setAtom2(CA);
    torsion.atom3 = CAnext;
    torsion.nideal = 0;
    torsion.tolerance = 5.0F;
    torsions.push_back(torsion);

    if (emap)
    {
        TorsionOptimize(atoms, model, emap, torsions, false);
    }

    aC1->removeType(AtomType::TORSIONATOM);
    aO1->removeType(AtomType::TORSIONATOM);
    aN2->removeType(AtomType::TORSIONATOM);

    if (strcmp(res->type().c_str(), "MRK") == 0 || strcmp(res->type().c_str(), "VEC") == 0)
    {
        res->setType("ALA");
    }

    return true;

}

int MIMolOpt::FindNeighbours(MIAtomList &CurrentAtoms, MIAtomList &Neighbours,
                             MIMoleculeBase *fitmol, float distance)
{
    Residue *res;
    MIAtom *a;
    unsigned int i, j;
    for (res = fitmol->getResidues(); Monomer::isValid(res); res = res->next())
    {
        for (i = 0; i < (unsigned int)res->atomCount(); i++)
        {
            a = res->atom(i);
            if (a->name()[0] != 'H')
            {
                for (j = 0; j < CurrentAtoms.size(); j++)
                {
                    if (CurrentAtoms[j] == a)
                    {
                        //skip rest of this residue
                        goto nextresidue;
                    }
                }
                for (j = 0; j < CurrentAtoms.size(); j++)
                {
                    if (AtomDist(*CurrentAtoms[j], *a) < distance)
                    {
                        Neighbours.push_back(a);
                        break;
                    }
                }
            }
        }
nextresidue:;
    }
    res = fitmol->getSymmResidues();
    while (Monomer::isValid(res))
    {
        for (i = 0; i < (unsigned int)res->atomCount(); i++)
        {
            a = res->atom(i);
            if (a->name()[0] != 'H')
            {
                for (j = 0; j < CurrentAtoms.size(); j++)
                {
                    if (CurrentAtoms[j] == a)
                    {
                        //skip rest of this residue
                        goto nextsymmresidue;
                    }
                }
                for (j = 0; j < CurrentAtoms.size(); j++)
                {
                    if (AtomDist(*CurrentAtoms[j], *a) < distance)
                    {
                        Neighbours.push_back(a);
                        break;
                    }
                }
            }
        }
nextsymmresidue:;
        res = res->next();
    }
    return (int)CurrentAtoms.size();
}

float MIMolOpt::StdevBonds()
{
    float sumd = 0, d;
    if (dict.RefiBonds.size() == 0)
    {
        return 0.0;
    }
    for (unsigned int i = 0; i < dict.RefiBonds.size(); i++)
    {
        d = (float) AtomDist(*dict.RefiBonds[i].getAtom1(), *dict.RefiBonds[i].getAtom2());
        d = d - dict.RefiBonds[i].ideal_length;
        sumd += d*d;
    }
    return (float)sqrt(sumd/(double)dict.RefiBonds.size());
}

float MIMolOpt::StdevAngles()
{
    float sumd = 0, d;
    if (dict.RefiAngles.size() == 0)
    {
        return 0;
    }
    for (unsigned int i = 0; i < dict.RefiAngles.size(); i++)
    {
        d = (float)AtomDist(*dict.RefiAngles[i].getAtom1(), *dict.RefiAngles[i].atom3);
        d = d - dict.RefiAngles[i].ideal_angle;
        sumd += d*d;
    }
    return (float)sqrt(sumd/(double)dict.RefiAngles.size());
}

float MIMolOpt::StdevPlanes()
{
    unsigned int i;
    float sumplane = 0;
    float diff, d;
    MIAtom *a1;
    if (dict.RefiPlanes.size() == 0)
    {
        return 0;
    }
    for (i = 0; i < dict.RefiPlanes.size(); i++)
    {
        /* build derivatives */
        lsqplane(dict.RefiPlanes[i]);
        int nsum = 0;
        float dsum = 0;
        for (int j = 0; j < (int)dict.RefiPlanes[i].natoms; j++)
        {
            a1 = dict.RefiPlanes[i].atoms[j];
            d =  a1->x() * dict.RefiPlanes[i].vm[0];
            d += a1->y() * dict.RefiPlanes[i].vm[1];
            d += a1->z() * dict.RefiPlanes[i].vm[2];
            d -= dict.RefiPlanes[i].d;
            d = d*d;
            dsum += d;
            nsum++;
        }
        diff = (float)sqrt(dsum/(double)nsum);
        sumplane += diff*diff;
    }
    return (float)sqrt(sumplane/(double)dict.RefiPlanes.size());
}

float MIMolOpt::StdevTorsions()
{
    if (dict.RefiTorsions.size() == 0)
    {
        return 0;
    }
    float chi, dchi, d, ideal;
    float sumd = 0;
    unsigned int i, j;
    for (i = 0; i < dict.RefiTorsions.size(); i++)
    {
        chi = (float)CalcAtomTorsion(dict.RefiTorsions[i].getAtom1(), dict.RefiTorsions[i].getAtom2(), dict.RefiTorsions[i].atom3, dict.RefiTorsions[i].atom4);
        if (chi < 0.0)
        {
            chi += 360.0;
        }
        dchi = 0.0;
        for (j = 0; j < (unsigned int)dict.RefiTorsions[i].nideal; j++)
        {
            d = dict.RefiTorsions[i].ideal[j] - chi;
            if (d < -180.0)
            {
                d += 360.0;
            }
            if (d >  180.0)
            {
                d -= 360.0;
            }
            if (fabs(d) < fabs(dchi) || j == 0)
            {
                dchi = d;
                ideal = dict.RefiTorsions[i].ideal[j];
            }
        }
        sumd += dchi*dchi;
    }
    return (float)sqrt(sumd/(double)dict.RefiTorsions.size());
}

static void molrep_score(trial &t, MIAtomList &atoms, EMapBase *emap)
{
    RotateAtomVec(t.p[0], t.p[1], t.p[2], 0.0f, 0.0f, 0.0f, &atoms);
    float fx, fy, fz;
    fx = (float)t.p[3];
    fy = (float)t.p[4];
    fz = (float)t.p[5];
    emap->mapheader->FtoC(&fx, &fy, &fz);
    TranslateAtomVec(fx, fy, fz, &atoms);
    //emap->SFCalc(model->getResidues());
    t.score = emap->CorrScore(atoms);
}

void MIMolOpt::MolecularReplace(MIMoleculeBase *model, EMapBase *emap)
{
    MIAtomList CurrentAtoms;
    unsigned int SaveToken, StartToken;
    CurrentMap = emap;
    float maxt = 0.5F;
    float maxr = 180.0F;
    float cx = 0, cy = 0, cz = 0;
    unsigned int i, j, k;
    trial t;
    unsigned int p1, p2, p3;
    trial *ti, *ta, *tb, *tc, *t2;
    float F = 0.8F;
    float CR = 0.9F;
    int D = 6;
    unsigned int npop = 20*D;
    float rms;

    if (!model || !emap)
    {
        return;
    }

    Residue *res = model->getResidues();
    while (Monomer::isValid(res))
    {
        for (int i = 0; i < res->atomCount(); i++)
        {
            CurrentAtoms.push_back(res->atom(i));
        }
        res = res->next();
    }

    if (!emap->HasPhases() || CurrentAtoms.size() == 0)
    {
        return;
    }

    double resmin = 5.0F;
    double resmax = 1000.0F;
    /* swap if sizes reversed */
    if (resmin > resmax)
    {
        double t = resmax;
        resmax = resmin;
        resmin = t;
    }
    emap->mapheader->resmin = (float)std::max(resmin, (0.5/emap->refls_stholmax));
    emap->mapheader->resmax = (float)std::min(resmax, (0.5/emap->refls_stholmin));
    emap->mapheader->maptype = MIMapType::Fo;
    float grid = 2.5F;
    emap->mapheader->nx = MIMapFactor((int)(grid*emap->mapheader->a/emap->mapheader->resmin), FFT_PRIME, EVEN, 2);
    emap->mapheader->ny = MIMapFactor((int)(grid*emap->mapheader->b/emap->mapheader->resmin), FFT_PRIME, EVEN, 2);
    emap->mapheader->nz = MIMapFactor((int)(grid*emap->mapheader->c/emap->mapheader->resmin), FFT_PRIME, EVEN, 2);
    emap->mapheader->hmax = (int)(emap->mapheader->a/emap->mapheader->resmin);
    emap->mapheader->kmax = (int)(emap->mapheader->b/emap->mapheader->resmin);
    emap->mapheader->lmax = (int)(emap->mapheader->c/emap->mapheader->resmin);
    while (emap->mapheader->nx <= 2*emap->mapheader->hmax)
    {
        emap->mapheader->nx = MIMapFactor(emap->mapheader->nx+1, FFT_PRIME, EVEN, 2);
    }
    while (emap->mapheader->ny <= 2*emap->mapheader->kmax)
    {
        emap->mapheader->ny = MIMapFactor(emap->mapheader->ny+1, FFT_PRIME, EVEN, 2);
    }
    while (emap->mapheader->nz <= 2*emap->mapheader->lmax)
    {
        emap->mapheader->nz = MIMapFactor(emap->mapheader->nz+1, FFT_PRIME, EVEN, 2);
    }

    ConnectTo(model); //Connect up signals
    StartToken = geomsaver.Save(CurrentAtoms, model);

    // find the center of mass of the atoms
    GetCenter(CurrentAtoms, cx, cy, cz);

    // translate to origin
    model->Translate(-cx, -cy, -cz, &CurrentAtoms);
    GetCenter(CurrentAtoms, cx, cy, cz);

    SaveToken = geomsaver.Save(CurrentAtoms, model);

    vector<float> var_start;
    var_start.reserve(D);
    for (j = 0; j < 3; j++)
    {
        var_start.push_back(maxr);
    }
    for (j = 3; j < 6; j++)
    {
        var_start.push_back(maxt);
    }

    vector<trial> population;
    vector<trial> population2;
    vector<trial>::iterator p;
    std::string s;
    double best_r, worst_r;
    int igen;

    for (i = 0; i < npop; i++)
    {
        t.p = new double[D];
        for (j = 0; j < (unsigned int)D; j++)
        {
            t.p[j] = frand2(var_start[j]);
        }
        population.push_back(t);
        t.p = new double[D];
        population2.push_back(t);
    }

    float score_sum = 0;

    for (p = population.begin(); p != population.end(); p++)
    {
        molrep_score(*p, CurrentAtoms, emap);
        score_sum += (float) (p->score*p->score);
        geomsaver.Restore(SaveToken);
    }
    sort(population.begin(), population.end(), trial_compare);
    best_r = population[0].score;
    worst_r = population[population.size()-1].score;
    rms = (float)sqrt(score_sum/(double)population.size());

    for (igen = 0; igen < 500; igen++)
    {
        for (i = 0; i < npop; i++)
        {
            ti = &population[i];
            do
            {
                p1 = irand(npop);
            } while (p1 == i);
            ta = &population[p1];
            do
            {
                p2 = irand(npop);
            } while (p2 == i || p2 == p1);
            tb = &population[p2];
            do
            {
                p3 = irand(npop);
            } while (p3 == i || p3 == p1 || p3 == p2);
            tc = &population[p3];
            j = irand(D);
            for (k = 1; k <= (unsigned int)D; k++)
            {
                if (frand() < CR || k == (unsigned int)D)
                {
                    t.p[j] = tc->p[j] + F*(ta->p[j]-tb->p[j]);
                }
                else
                {
                    t.p[j] = ti->p[j];
                }
                j = j+1;
                if (j >= (unsigned int)D)
                {
                    j = j-(unsigned int)D;
                }
            }
            molrep_score(t, CurrentAtoms, emap);
            geomsaver.Restore(SaveToken);
            t2 = &population2[i];
            if (t.score >= ti->score)
            {
                copy_trial(t2, &t, D);
            }
            else
            {
                copy_trial(t2, ti, D);
            }
        }
        for (i = 0; i < npop; i++)
        {
            copy_trial(&population[i], &population2[i], D);
        }
        sort(population.begin(), population.end(), trial_compare);
        best_r = population[0].score;
        worst_r = population[population.size()-1].score;
        Logger::footer("Generation %d: best=%0.2f worst=%0.2f rms=%0.2f", igen, best_r, worst_r, rms);
        if (best_r > 5.0 * rms)
        {
            Logger::log("Ended after generation %d with score 5*rms", igen);
            break;
        }
    }

    if (best_r > 0.0)
    {
        RotateAtomVec(population[0].p[0], population[0].p[1], population[0].p[2], cx, cy, cz, &CurrentAtoms);
        TranslateAtomVec(population[0].p[3], population[0].p[4], population[0].p[5], &CurrentAtoms);
    }
    else
    {
        Logger::log("Final model no better or worse than start - model not moved");
    }

    Logger::log("RFinal=%0.2f\nRot: %0.2f %0.2f %0.2f Trans: %0.2f %0.2f %0.2f", best_r,
                population[0].p[0], population[0].p[1], population[0].p[2], population[0].p[3], population[0].p[4], population[0].p[5]);
    for (i = 0; i < population.size(); i++)
    {
        delete[] population[i].p;
        delete[] population2[i].p;
    }
}

// checks to see is an atom is in those being currently refined
bool MIMolOpt::IsBeingRefined(MIAtom *atom)
{
    if (nRefiRes == 0 || RefiRes == NULL)
    {
        return false;
    }
    Residue *res = RefiRes;
    int i;
    int n = 0;
    while (Monomer::isValid(res) && n < nRefiRes)
    {
        for (i = 0; i < res->atomCount(); i++)
        {
            if (atom == res->atom(i))
            {
                return true;
            }
        }
        n++;
        res = res->next();
    }
    return false;
}

bool TorsionMatch::operator ()(const chemlib::TORSION &t, const chemlib::TORSDICT &td) const
{
    return strncmp(t.res->type().c_str(), td.restype, MAXNAME) == 0
           && strncmp(t.type, td.type, 11) == 0;
}

void MIMolOpt::internalSetRefiRes(chemlib::Residue *residue, int nResidues)
{
    bool wasRefining = IsRefining();
    RefiRes = residue;
    nRefiRes = nResidues;
    if (wasRefining != IsRefining())
    {
        isRefiningChanged(IsRefining());
    }
}
