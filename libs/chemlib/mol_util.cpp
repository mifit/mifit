#include <cstring>
#include <cmath>
#include <algorithm>

#include <math/mathlib.h>
#include <chemlib/RESIDUE_.h>

#include "mol_util.h"
#include "atom_util.h"
#include "CovalentGeom.h"
#include "valence.h"
#include "mol_util_private.h"

#ifdef _WIN32
#define strncasecmp strnicmp
#endif

using namespace std;

namespace chemlib
{

/////////////////////////////////////////////////////////////////////////////
// Function:    Atomic_Name
// Purpose:		Look up the symbol for an element from the periodic table
// Input:       The atomic number of the element
// Output:      A two-character string, right-justified, all-caps, containing the symbol
//				Returns carbon if the atomic number is zero or is greater than NELEMENTS
// Requires:
/////////////////////////////////////////////////////////////////////////////
const char *Atomic_Name(int atomic_number)
{
    static char per_tab[NELEMENTS+1][3] =
    {
        " C",
        " H", "HE",
        "LI", "BE", " B", " C", " N", " O", " F", "NE",
        "NA", "MG", "AL", "SI", " P", " S", "CL", "AR",
        " K", "CA", "SC", "TI", " V", "CR", "MN", "FE", "CO",
        "NI", "CU", "ZN", "GA", "GE", "AS", "SE", "BR", "KR",
        "RB", "SR", " Y", "ZR", "NB", "MO", "TC", "RU", "RH",
        "PD", "AG", "CD", "IN", "SN", "SB", "TE", " I", "XE",
        "CS", "BA", "LA",
        "CE", "PR", "ND", "PM", "SM", "EU", "GD", "TB", "DY",
        "HO", "ER", "TM", "YB", "LU",
        "HF", "TA", " W", "RE", "OS", "IR",
        "PT", "AU", "HG", "TL", "PB", "BI", "PO", "AT", "RN",
        "FR", "RA", "AC",
        "TH", "PA", " U", "NP", "PU", "AM", "CM", "BK", "CF",
        "ES", "FM", "MD", "NO", "LR",
        "RF"
    };
    if (atomic_number > NELEMENTS)
    {
        return per_tab[6];
    }
    return per_tab[atomic_number];
}

/////////////////////////////////////////////////////////////////////////////
// Function:    Left_Atomic_Name
// Purpose:		Look up the symbol for an element from the periodic table
// Input:       The atomic number of the element
// Output:      A left-justified string (one or two chars), all-caps, containing the symbol
//				Returns carbon if the atomic number is zero or is greater than NELEMENTS
// Requires:
/////////////////////////////////////////////////////////////////////////////
const char *Left_Atomic_Name(int atomic_number)
{
    static char per_tab[NELEMENTS+1][3] =
    {
        "C",
        "H",  "HE",
        "LI", "BE", "B",  "C",  "N",  "O",  "F",  "NE",
        "NA", "MG", "AL", "SI", "P",  "S",  "CL", "AR",
        "K",  "CA", "SC", "TI", "V",  "CR", "MN", "FE", "CO",
        "NI", "CU", "ZN", "GA", "GE", "AS", "SE", "BR", "KR",
        "RB", "SR", "Y",  "ZR", "NB", "MO", "TC", "RU", "RH",
        "PD", "AG", "CD", "IN", "SN", "SB", "TE", "I",  "XE",
        "CS", "BA", "LA",
        "CE", "PR", "ND", "PM", "SM", "EU", "GD", "TB", "DY",
        "HO", "ER", "TM", "YB", "LU",
        "HF", "TA", "W",  "RE", "OS", "IR",
        "PT", "AU", "HG", "TL", "PB", "BI", "PO", "AT", "RN",
        "FR", "RA", "AC",
        "TH", "PA", "U",  "NP", "PU", "AM", "CM", "BK", "CF",
        "ES", "FM", "MD", "NO", "LR",
        "RF"
    };
    if (atomic_number > NELEMENTS)
    {
        return per_tab[6];
    }
    return per_tab[atomic_number];
}

/////////////////////////////////////////////////////////////////////////////
// Function:    Atomic_Number
// Purpose:		Look up an element using its two-letter symbol
// Input:       A string with the right-justified, all caps symbol
// Output:      The atomic number of the element, or -1 if not found
// Requires:
/////////////////////////////////////////////////////////////////////////////
int Atomic_Number(const char *symbol)
{
    static char per_tab[NELEMENTS+1][3] =
    {
        " *",
        " H", "HE",
        "LI", "BE", " B", " C", " N", " O", " F", "NE",
        "NA", "MG", "AL", "SI", " P", " S", "CL", "AR",
        " K", "CA", "SC", "TI", " V", "CR", "MN", "FE", "CO",
        "NI", "CU", "ZN", "GA", "GE", "AS", "SE", "BR", "KR",
        "RB", "SR", " Y", "ZR", "NB", "MO", "TC", "RU", "RH",
        "PD", "AG", "CD", "IN", "SN", "SB", "TE", " I", "XE",
        "CS", "BA", "LA",
        "CE", "PR", "ND", "PM", "SM", "EU", "GD", "TB", "DY",
        "HO", "ER", "TM", "YB", "LU",
        "HF", "TA", " W", "RE", "OS", "IR",
        "PT", "AU", "HG", "TL", "PB", "BI", "PO", "AT", "RN",
        "FR", "RA", "AC",
        "TH", "PA", " U", "NP", "PU", "AM", "CM", "BK", "CF",
        "ES", "FM", "MD", "NO", "LR",
        "RF"
    };
    for (int i = 0; i <= NELEMENTS; ++i)
    {
        if (symbol[0] == per_tab[i][0]          //Compare the first two chars
            && symbol[1] == per_tab[i][1])
        {
            return i;
        }
    }

    if (symbol[0] == ' ' && symbol[1] == 'D')
    {
        return 1;
    }
    return -1;
}

int Atomic_Number_Nformat(const std::string &atname)
{
    switch (atname.length())
    {
    case 0:
        return 0;
    case 1:
        char symbol[3];
        symbol[0] = ' ';
        symbol[1] = atname[0];
        symbol[2] = '\0';
        return Atomic_Number(symbol);
    default:
        return Atomic_Number(atname.c_str());
    }
}

//Covalent radii from CSD website
//http://www.ccdc.cam.ac.uk/products/csd/radii/
float CovalentRadius(int atomic_number)
{
    static float rad_tab[NELEMENTS+1] =
    {
        1.50F,
        0.23F, 1.50F,
        0.68F, 0.35F, 0.83F, 0.68F, 0.68F, 0.68F, 0.64F, 1.50F,
        0.97F, 1.10F, 1.35F, 1.20F, 1.05F, 1.02F, 0.99F, 1.51F,
        1.33F, 0.99F, 1.44F, 1.47F, 1.33F, 1.35F, 1.35F, 1.34F, 1.33F,
        1.50F, 1.52F, 1.45F, 1.22F, 1.17F, 1.21F, 1.22F, 1.21F, 1.50F,
        1.47F, 1.12F, 1.78F, 1.56F, 1.48F, 1.47F, 1.35F, 1.40F, 1.45F,
        1.50F, 1.59F, 1.69F, 1.63F, 1.46F, 1.46F, 1.47F, 1.40F, 1.50F,
        1.67F, 1.34F,
        1.87F, 1.83F, 1.82F, 1.81F, 1.80F, 1.80F, 1.99F, 1.79F, 1.76F,
        1.75F, 1.74F, 1.73F, 1.72F, 1.94F,
        1.72F, 1.57F, 1.43F, 1.37F, 1.35F, 1.37F, 1.32F,
        1.50F, 1.50F, 1.70F, 1.55F, 1.54F, 1.54F, 1.68F, 1.21F, 1.50F,
        1.50F, 1.90F,
        1.88F, 1.79F, 1.61F, 1.58F, 1.55F, 1.53F, 1.51F, 0.99F, 1.54F,
        1.83F, 1.50F, 1.50F, 1.50F, 1.50F, 1.50F, 1.50F
    };

    if (atomic_number > NELEMENTS || atomic_number < 0)
    {
        return 1.50F;
    }
    else
    {
        return rad_tab[atomic_number];
    }
}

/////////////////////////////////////////////////////////////////////////////
// Function:    Electroneg
// Purpose:		Look up the electronegativity of an atom
// Input:       The atomic number of the atom
// Output:      100x the electronegativity, on the Pauling scale
// Requires:	Valid atomic numbers range from 1 to 104
/////////////////////////////////////////////////////////////////////////////
int Electroneg(int atomic_number)
{
    static int elneg_tab[NELEMENTS+1] =
    {
        130,
        220,   0,                                       //H -> He
        98, 157, 204, 255, 304, 344, 398, 426,
        93, 131, 161, 190, 219, 258, 316, 311,          //Na -> Ar
        82, 100, 136, 154, 163, 166, 155, 183, 188,
        191, 190, 165, 181, 201, 218, 255, 296, 300,
        82,  95, 122, 133, 160, 216, 190, 220, 228, //Rb -> Rh
        220, 193, 169, 178, 196, 205, 210, 266, 260,
        79,  89,
        110, 112, 113, 114, 113, 117, 120, 120, 110,    //Lanthanides
        122, 123, 124, 125, 110,                        //
        127, 130, 150, 236, 190, 220, 220,
        228, 254, 200, 162, 233, 202, 200, 220, 212,
        70,  90,
        110, 130, 150, 138, 136, 128, 130, 130, 130,
        130, 130, 130, 130, 130, 130, 130
    };

    if (atomic_number > NELEMENTS || atomic_number < 0)
    {
        return 130;
    }
    else
    {
        return elneg_tab[atomic_number];
    }
}

/////////////////////////////////////////////////////////////////////////////
// Function:    MaxValence
// Purpose:		Determine the maximum sum of bond orders that can be attached
//				to a give atom.
// Input:       The atomic number of the atom
// Output:      A value from 1 to 12 for valid input, or 0 for invalid input
// Requires:	Valid atomic numbers range from 1 to 104
/////////////////////////////////////////////////////////////////////////////
int MaxValence(int atomic_number)
{
    static int valnc_tab[NELEMENTS+1] =
    {
        0,
        1, 0,                   //H and He
        1, 4, 4, 4, 5, 2, 1, 0,     //Li to Ne
        1, 4, 6, 6, 5, 6, 7, 0,     //Na to Ar
        1, 2, 6, 6, 6, 6, 8, 6, 6,  //K  to Co
        6, 6, 6, 3, 4, 5, 6, 5, 6,  //Ni to Kr
        1, 2, 6, 6, 6, 6, 8, 6, 6,  //Rb to Rh
        6, 6, 6, 3, 4, 5, 6, 7, 8,  //Pd to Xe
        1, 2, 12,               //Cs to La
        6, 6, 6, 6, 6, 6, 6, 6, 6,  //Lanthanides
        6, 6, 6, 6, 6,
        6, 8, 6, 6, 6, 6,       //W to Au
        6, 6, 6, 3, 4, 5, 4, 1, 0,  //Hg to Rn
        1, 2, 6,                //Fr to Ac
        6, 6, 6, 6, 6, 6, 6, 6, 6,  //Actinides
        6, 6, 6, 6,   6,            //
        6                       //Transactinides
    };

    if (atomic_number > NELEMENTS || atomic_number < 0)
    {
        return 0;
    }
    else
    {
        return valnc_tab[atomic_number];
    }
}

double BondLimit(const char *s)
{
    double d = 0.95;
    if (!strncmp(s, "C", 1))
    {
        d = 0.90;
    }
    if (!strncmp(s, "N", 1))
    {
        d = 0.95;
    }
    if (!strncmp(s, "O", 1))
    {
        d = 0.95;
    }
    if (!strncmp(s, "S", 1))
    {
        d = 1.20;
    }
    if (!strncmp(s, "P", 1))
    {
        d = 1.10;
    }
    if (!strncmp(s, "H", 1))
    {
        d = 0.40;
    }
    if (isdigit(s[0]))
    {
        d = 0.40;
    }
    if (!strncmp(s, "FE", 2))
    {
        d = 1.30;
    }
    if (!strncmp(s, "BR", 2))
    {
        d = 1.20;
    }
    if (!strncmp(s, "I", 2))
    {
        d = 1.30;
    }
    if (!strncmp(s, "CL", 2))
    {
        d = 1.20;
    }
    if (!strncmp(s, "MG", 2))
    {
        d = 1.20;
    }
    if (!strncmp(s, "MN", 2))
    {
        d = 1.20;
    }
    if (!strncmp(s, "PT", 2))
    {
        d = 1.50;
    }
    if (!strncmp(s, "HG+", 3))
    {
        d = 1.50;
    }
    if (!strncmp(s, "U", 1))
    {
        d = 1.40;
    }
    if (!strncmp(s, "CA+", 3))
    {
        d = 1.20;
    }
    return (d);
}

bool IsBondRotatable::operator()(const Bond &bond) const
{
    if (bond.iscyclic)                          //If this bond is in a ring, it is
    {
        return false;                       //not rotatable
    }
    else if (bond.getAtom1()->nabors().size() == 1)     //If this bond is terminal, it is not
    {
        return false;                       //rotatable
    }
    else if (bond.getAtom2()->nabors().size() == 1)
    {
        return false;
    }

    return true;                            //If we reach this line, the bond is rotatable

}

/////////////////////////////////////////////////////////////////////////////
// Function:    GatherAtmPtrs
// Purpose:		Compile pointers to all the atoms in a group of residues
// Input:       A vector of pointers to residues
// Output:      Writes the atom pointers to a vector
// Requires:
/////////////////////////////////////////////////////////////////////////////
void GatherAtmPtrs(MIAtomList &target, const std::vector<Residue*> &source)
{

    target.reserve((source.front())->atoms().size() * source.size());

    std::vector<Residue*>::const_iterator res;
    for (res = source.begin(); res != source.end(); ++res)
    {
        for (unsigned int i = 0; i < (*res)->atoms().size(); ++i)
        {
            target.push_back((*res)->atom(i));
        }
    }
}

/////////////////////////////////////////////////////////////////////////////
// Function:    EnumerateTorsions
// Purpose:     Generates all connected four-atom sequences centered on a given bond
//              in a (molecular) graph
// Input:       A pointer to a bond, and a vector of atom-vectors to be filled on output
// Output:      Fills the vector with four-atom vectors
//              containing the bond
// Requires:
/////////////////////////////////////////////////////////////////////////////
void EnumerateTorsions(const Bond *bond, std::vector< MIAtomList > &torsions)
{
    MIAtom *atom2, *atom3;
    MIAtomList torsion;                   //Vector of four atoms

    atom2 = bond->getAtom1();
    atom3 = bond->getAtom2();

    MIAtom_const_iter nbr1;
    MIAtom_const_iter nbr2;
    for (nbr1 = atom2->nabors().begin(); nbr1 != atom2->nabors().end(); ++nbr1)
    {
        if (*nbr1 == atom3)
        {
            continue;
        }

        for (nbr2 = atom3->nabors().begin(); nbr2 != atom3->nabors().end(); ++nbr2)
        {
            if (*nbr2 == atom2 || *nbr2 == *nbr1)
            {
                continue;
            }

            torsion.push_back(*nbr1);   //End atoms are variable
            torsion.push_back(atom2);
            torsion.push_back(atom3);   //Two center atoms are fixed
            torsion.push_back(*nbr2);

            torsions.push_back(torsion);
            torsion.clear();
        }
    }
}

double CalcAtomTorsion(const MIAtom *na, const MIAtom *nb, const MIAtom *nc, const MIAtom *nd)
{
    double bax, bay, baz;
    double bcx, bcy, bcz;
    double cbx, cby, cbz;
    double cdx, cdy, cdz;
    double babcx, babcy, babcz;
    double cbcdx, cbcdy, cbcdz;
    double xx, xy, xz;
    double bcdx, cost, sint;

    bax =  (nb->x() - na->x());
    bay =  (nb->y() - na->y());
    baz =  (nb->z() - na->z());

    bcx =  (nb->x() - nc->x());
    bcy =  (nb->y() - nc->y());
    bcz =  (nb->z() - nc->z());

    cbx =  (nc->x() - nb->x());
    cby =  (nc->y() - nb->y());
    cbz =  (nc->z() - nb->z());

    cdx =  (nc->x() - nd->x());
    cdy =  (nc->y() - nd->y());
    cdz =  (nc->z() - nd->z());

    babcx = bay*bcz - bcy*baz;
    babcy = baz*bcx - bcz*bax;
    babcz = bax*bcy - bcx*bay;

    cbcdx = cby*cdz - cdy*cbz;
    cbcdy = cbz*cdx - cdz*cbx;
    cbcdz = cbx*cdy - cdx*cby;

    xx = cbcdy*babcz - babcy*cbcdz;
    xy = cbcdz*babcx - babcz*cbcdx;
    xz = cbcdx*babcy - babcx*cbcdy;

    bcdx = bcx*xx + bcy*xy + bcz*xz;
    cost = babcx*cbcdx + babcy*cbcdy + babcz*cbcdz;
    sint = (double)sqrt(xx*xx + xy*xy + xz*xz);
    if (bcdx < 0)
    {
        sint *= -1;
    }
    return ((double)(RAD2DEG * atan2(sint, cost)));
}

double CalcAtomTorsion(MIAtom *na, MIAtom *nb, MIAtom *nc, MIAtom *nd)
{
    double bax, bay, baz;
    double bcx, bcy, bcz;
    double cbx, cby, cbz;
    double cdx, cdy, cdz;
    double babcx, babcy, babcz;
    double cbcdx, cbcdy, cbcdz;
    double xx, xy, xz;
    double bcdx, cost, sint;
    bax =  (nb->x() - na->x());
    bay =  (nb->y() - na->y());
    baz =  (nb->z() - na->z());

    bcx =  (nb->x() - nc->x());
    bcy =  (nb->y() - nc->y());
    bcz =  (nb->z() - nc->z());

    cbx =  (nc->x() - nb->x());
    cby =  (nc->y() - nb->y());
    cbz =  (nc->z() - nb->z());

    cdx =  (nc->x() - nd->x());
    cdy =  (nc->y() - nd->y());
    cdz =  (nc->z() - nd->z());

    babcx = bay*bcz - bcy*baz;
    babcy = baz*bcx - bcz*bax;
    babcz = bax*bcy - bcx*bay;

    cbcdx = cby*cdz - cdy*cbz;
    cbcdy = cbz*cdx - cdz*cbx;
    cbcdz = cbx*cdy - cdx*cby;

    xx = cbcdy*babcz - babcy*cbcdz;
    xy = cbcdz*babcx - babcz*cbcdx;
    xz = cbcdx*babcy - babcx*cbcdy;

    bcdx = bcx*xx + bcy*xy + bcz*xz;
    cost = babcx*cbcdx + babcy*cbcdy + babcz*cbcdz;
    sint = (double)sqrt(xx*xx + xy*xy + xz*xz);
    if (bcdx < 0)
    {
        sint *= -1;
    }
    return ((double)((180.0/M_PI)*atan2(sint, cost)));
}

/* fit a least squares plane to a set of points
 * method from schomaker et al, acta cryst, 12, p. 600, 1959.
 */
void lsqplane(PLANE *plane)
{
    float xs[3], xxs[3][3], b[3][3];
    float a[3][3], vmi[3], bv[3], vm[3], d;
    float zip = 1.0E-5F;
    float orm, vm0, ratio0, ratio1, ratio2, rat01, rat02;
    int i, j, kk, nnn, k, n;

    n = plane->natoms;
    for (i = 0; i < 3; i++)
    {
        xs[i] = 0.0;
    }
    for (k = 0; k < n; k++)
    {
        /* zip is added to prevent numerical instability if
         * atoms in a plane = 0
         */
        xs[0] += plane->atoms[k]->x()+zip;
        xs[1] += plane->atoms[k]->y()+zip;
        xs[2] += plane->atoms[k]->z()+zip;
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
        xxs[0][0] += plane->atoms[k]->x() * plane->atoms[k]->x();
        xxs[0][1] += plane->atoms[k]->x() * plane->atoms[k]->y();
        xxs[0][2] += plane->atoms[k]->x() * plane->atoms[k]->z();
        xxs[1][0] += plane->atoms[k]->y() * plane->atoms[k]->x();
        xxs[1][1] += plane->atoms[k]->y() * plane->atoms[k]->y();
        xxs[1][2] += plane->atoms[k]->y() * plane->atoms[k]->z();
        xxs[2][0] += plane->atoms[k]->z() * plane->atoms[k]->x();
        xxs[2][1] += plane->atoms[k]->z() * plane->atoms[k]->y();
        xxs[2][2] += plane->atoms[k]->z() * plane->atoms[k]->z();
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
        rat01 = (float)fabs(ratio1/ratio0-1.0);
        rat02 = (float)fabs(ratio2/ratio0-1.0);
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
    plane->vm[0] = vm[0];
    plane->vm[1] = vm[1];
    plane->vm[2] = vm[2];
    plane->d = d;
}

double CalcAtomAngle(const MIAtom &a1, const MIAtom &a2, const MIAtom &a3)
{
    //angle between vectorsfrom a2->a1 to a2->a3
    double v1[3], v2[3], a;
    v1[0] = a1.x() - a2.x();
    v1[1] = a1.y() - a2.y();
    v1[2] = a1.z() - a2.z();
    v2[0] = a3.x() - a2.x();
    v2[1] = a3.y() - a2.y();
    v2[2] = a3.z() - a2.z();
    double d12 = AtomDist(a1, a2);
    double d23 = AtomDist(a2, a3);
    if (d12 != 0.0 && d23 != 0.0)
    {
        a = RAD2DEG * acos((v1[0]*v2[0] + v1[1]*v2[1] + v1[2]*v2[2])
                           /(d12*d23));
    }
    else
    {
        a = 0.0;
    }
    if (a < 0.0)
    {
        a = a+360.0;
    }
    if (a >= 360.0)
    {
        a = a -360.0;
    }
    if (a >= 180.0)
    {
        a = 360.0-a;
    }
    return (a);

    return a;
}

/* fit a least squares plane to a set of points
 * method from schomaker et al, acta cryst, 12, p. 600, 1959.
 */
void lsqplane(Plane &plane, float *normal, float *displace)
{
    float xs[3], xxs[3][3], b[3][3];
    float a[3][3], vmi[3], bv[3], vm[3], d;
    float zip = 1.0E-5F;
    float orm, vm0, ratio0, ratio1, ratio2, rat01, rat02;
    int i, j, kk, nnn, k, n;

    n = plane.NumAtoms();
    for (i = 0; i < 3; i++)
    {
        xs[i] = 0.0;
    }
    for (k = 0; k < n; k++)
    {
        /* zip is added to prevent numerical instability if
         * atoms in a plane = 0
         */
        xs[0] += plane.GetAtom(k)->x() + zip;
        xs[1] += plane.GetAtom(k)->y() + zip;
        xs[2] += plane.GetAtom(k)->z() + zip;
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
        xxs[0][0] += plane.GetAtom(k)->x() * plane.GetAtom(k)->x();
        xxs[0][1] += plane.GetAtom(k)->x() * plane.GetAtom(k)->y();
        xxs[0][2] += plane.GetAtom(k)->x() * plane.GetAtom(k)->z();
        xxs[1][0] += plane.GetAtom(k)->y() * plane.GetAtom(k)->x();
        xxs[1][1] += plane.GetAtom(k)->y() * plane.GetAtom(k)->y();
        xxs[1][2] += plane.GetAtom(k)->y() * plane.GetAtom(k)->z();
        xxs[2][0] += plane.GetAtom(k)->z() * plane.GetAtom(k)->x();
        xxs[2][1] += plane.GetAtom(k)->z() * plane.GetAtom(k)->y();
        xxs[2][2] += plane.GetAtom(k)->z() * plane.GetAtom(k)->z();
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
        rat01 = (float)fabs(ratio1/ratio0-1.0);
        rat02 = (float)fabs(ratio2/ratio0-1.0);
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
    normal[0] = vm[0];
    normal[1] = vm[1];
    normal[2] = vm[2];
    *displace = d;
}

void DepthFirstSearch(MIAtom *root,
                      const MIAtom *prev,
                      const MIAtom *block,
                      std::vector <MIAtom*> &aggregate)
{

    if (std::find(aggregate.begin(), aggregate.end(), root) == aggregate.end())
    {
        aggregate.push_back(root);              //Add this atom to the aggregation
    }
    else
    {
        return;
    }

    std::vector <MIAtom*>::const_iterator nabor;
    for (nabor = root->nabors().begin(); nabor != root->nabors().end(); ++nabor)
    {
        if (*nabor != prev && *nabor != block)
        {
            DepthFirstSearch(*nabor, root, block, aggregate);
        }
    }
}

void DepthFirstSearch(MIAtom *root,
                      const MIAtom *prev,
                      const MIAtom *block,
                      const vector<Bond> &bonds,
                      vector<MIAtom*> &aggregate)
{

    if (std::find(aggregate.begin(), aggregate.end(), root) == aggregate.end())
    {
        aggregate.push_back(root);              //Add this atom to the aggregation
    }
    else
    {
        return;                                 //Don't visit the same atom twice
    }

    std::vector <MIAtom*> nabors;
    GetNabors(root, bonds, nabors);

    std::vector <MIAtom*>::const_iterator nabor;
    for (nabor = nabors.begin(); nabor != nabors.end(); ++nabor)
    {
        if (*nabor != prev && *nabor != block)
        {
            DepthFirstSearch(*nabor, root, block, bonds, aggregate);
        }
    }
}

//Front-end for a depth-first search whose results will exclude the
//initial atom
void ExclDepthFirstSearch(MIAtom *root,
                          const MIAtom *block,
                          const vector<Bond> &bonds,
                          vector<MIAtom*> &aggregate)
{

    std::vector <MIAtom*> nabors;
    GetNabors(root, bonds, nabors);

    std::vector <MIAtom*>::iterator nabor;
    for (nabor = nabors.begin(); nabor != nabors.end(); ++nabor)
    {
        if (*nabor != block)
        {
            DepthFirstSearch(*nabor, root, root, bonds, aggregate);
        }
    }
}

void TrimBonds(std::vector<Bond> &trimmed_bonds,
               const std::vector<Bond> &orig_bonds,
               const MIAtomList &atoms)
{
    std::vector<Bond>::const_iterator bnd;
    for (bnd = orig_bonds.begin(); bnd != orig_bonds.end(); ++bnd)
    {
        if (std::find(atoms.begin(), atoms.end(), bnd->getAtom1()) == atoms.end())
        {
            continue;
        }
        if (std::find(atoms.begin(), atoms.end(), bnd->getAtom2()) == atoms.end())
        {
            continue;
        }
        trimmed_bonds.push_back(*bnd);
    }
}

int CountResBonds(const RESIDUE &res, const std::vector<Bond> &bonds)
{
    return count_if(bonds.begin(), bonds.end(), bind2nd(ResContainsBond(), res));
}

//	int CountResBonds(const RESIDUE &res, const std::vector<Bond> &bonds) {
//		int c = 0;
//
//		for (int i=0; i<bonds.size(); ++i) {
//
//			if (find(res.atoms(), res.atoms() + res.natoms, bonds[i].getAtom1()) != res.atoms() + res.natoms &&
//				find(res.atoms(), res.atoms() + res.natoms, bonds[i].getAtom2()) != res.atoms() + res.natoms) {
//				c++;
//			}
//		}
//		return c;
//	}


int Build(Ligand &lig)
{
    int i;
    MIAtom *atom;
    MIAtom *atom2;
    float x, y, z;
    //		int ix,iy,iz;
    double dx, dy, dz, dist = 1.90, d1;
    std::vector<Residue*>::iterator ri;
    std::vector<Residue*>::iterator nter = lig.residues.begin();
    short nextid;
    //		nresidues = 0;
    Bond bond;

    long nth = 0;
    char link_here[MAXNAME];
    char link_next[MAXNAME]; //MAXNAME is in Xguicryst.h

    //		nlinks = 0;
    lig.bonds.clear();
    //		InitSeqPos();
    lig.FixAtomicNumbers();


    for (ri = lig.residues.begin(); ri != lig.residues.end(); ++ri)
    {
        Residue *res = *ri;
        for (unsigned int i = 0; i < res->atoms().size(); i++)
        {
            res->atom(i)->removeType(AtomType::BONDED);
        }
    }


    //		xmax = ymax = zmax = -9999;
    //		xmin = ymin = zmin = 9999;

    for (ri = lig.residues.begin(); ri != lig.residues.end(); ++ri)
    {
        Residue *res = *ri;
        /*   gather atoms in residue and "pointtovu" them */
        for (unsigned int i = 1; i < res->atoms().size(); i++)
        {
            x = res->atom(i)->x();
            y = res->atom(i)->y();
            z = res->atom(i)->z();
            d1 = BondLimit(res->atom(i)->name());
            for (unsigned int j = 0; j < i; j++)
            {
                dx = res->atom(j)->x() -x;
                dy = res->atom(j)->y() -y;
                dz = res->atom(j)->z() -z;
                if (res->atom(j)->altloc() != ' ')
                {
                    /* check to see if OK to bond */
                    if (res->atom(i)->altloc() != ' ')
                    {
                        if (res->atom(j)->altloc() != res->atom(i)->altloc())
                        {
                            continue;
                        }
                    }
                }
                dist = d1 + BondLimit(res->atom(j)->name());
                if ( (dx < dist || dx > -dist) && (dy < dist || dy > -dist)
                     && hypot(dz, hypot(dx, dy)) < dist)
                {
                    //if(!alreadybonded(res->atom(i),res->atom(j))){
                    res->atom(i)->addNabor(res->atom(j));
                    res->atom(j)->addNabor(res->atom(i));
                    res->atom(i)->addBondnumber(lig.bonds.size());
                    res->atom(j)->addBondnumber(lig.bonds.size());

                    bond.setAtom1(res->atom(i));
                    bond.setAtom2(res->atom(j));
                    bond.getAtom1()->addType(AtomType::BONDED);
                    bond.getAtom2()->addType(AtomType::BONDED);
                    bond.type = B_NORMAL;
                    lig.bonds.push_back(bond);
                }
            }
        }
        // find the N atom of next res and look for bond-
        int is_dna = IsDna(*res);
        if (is_dna == 1)
        {
            strcpy(link_here, "O3*"); // dna link thru phosphate bond
            strcpy(link_next, "P");
        }
        else if (is_dna == 2)
        {
            strcpy(link_here, "O3\'"); // dna link thru phosphate bond
            strcpy(link_next, "P");
        }
        else
        {
            strcpy(link_here, "C"); // protein link thru peptide bond
            strcpy(link_next, "N"); //
        }
        if (ri+1 != lig.residues.end()
            && (atom = atom_from_name(link_next, **(ri+1))) != NULL
            && (atom2 = atom_from_name(link_here, *res)) != NULL)
        {
            x = atom->x();
            y = atom->y();
            z = atom->z();
            dx = atom2->x() - x;
            dy = atom2->y() - y;
            dz = atom2->z() - z;
            dist = BondLimit(atom->name())+BondLimit(atom2->name());
            if (atom2->altloc() != ' ')
            {
                /* check to see if OK to bond */
                if (atom->altloc() != ' ')
                {
                    if (atom->altloc() != atom2->altloc())
                    {
                        goto skip;
                    }
                }
            }
            if ( (dx < dist || dx > -dist) && (dy < dist || dy > -dist)
                 && hypot(dz, hypot(dx, dy)) < dist)
            {
                if (!lig.AlreadyBonded(atom, atom2))
                {
                    atom->addNabor(atom2);
                    atom2->addNabor(atom);
                    atom->addBondnumber(lig.bonds.size());
                    atom2->addBondnumber(lig.bonds.size());

                    bond.setAtom1(atom);
                    bond.setAtom2(atom2);
                    bond.getAtom1()->addType(AtomType::BONDED);
                    bond.getAtom2()->addType(AtomType::BONDED);
                    bond.type = B_NORMAL;
                    lig.bonds.push_back(bond);
                }
            }
        }
        // if the residue is last in chain look for cyclic peptide
skip:       if (ri+1 == lig.residues.end())
        {
            nextid = -1;
        }
        else
        {
            nextid = (res+1)->chain_id();
        }
        if ((res->chain_id()) != nextid)
        {
            if ((atom = atom_from_name(link_here, *res)) != NULL
                && (atom2 = atom_from_name(link_next, **nter)) != NULL)
            {
                x = atom->x();
                y = atom->y();
                z = atom->z();
                dx = atom2->x() -x;
                dy = atom2->y() -y;
                dz = atom2->z() -z;
                if (atom2->altloc() != ' ')
                {
                    /* check to see if OK to bond */
                    if (atom->altloc() != ' ')
                    {
                        if (atom->altloc() != atom2->altloc())
                        {
                            goto skip2;
                        }
                    }
                }
                dist = BondLimit(atom->name())+BondLimit(atom2->name());
                if ( (dx < dist || dx > -dist) && (dy < dist || dy > -dist)
                     && hypot(dz, hypot(dx, dy)) < dist)
                {
                    if (!lig.AlreadyBonded(atom, atom2))
                    {
                        atom->addNabor(atom2);
                        atom2->addNabor(atom);
                        atom->addBondnumber(lig.bonds.size());
                        atom2->addBondnumber(lig.bonds.size());

                        bond.setAtom1(atom);
                        bond.setAtom2(atom2);
                        bond.getAtom1()->addType(AtomType::BONDED);
                        bond.getAtom2()->addType(AtomType::BONDED);
                        bond.type = B_NORMAL;
                        lig.bonds.push_back(bond);
                    }
                }
            }


skip2:            if (ri+1 != lig.residues.end())
            {
                nter = ri+1;
            }
        }
        // finally if only 1 atom and its a CA link it to next CA
        if (ri+1 != lig.residues.end())
        {
            MIAtom *atom1;
            if (res->atoms().size() <= 2 && (atom1 = atom_from_name("CA", *res)) != NULL
                && (atom2 = atom_from_name("CA", *(res+1))) != NULL
                && (res->chain_id()) == ((res+1)->chain_id()))
            {
                if (!lig.AlreadyBonded(res->atom(0), atom2))
                {
                    res->atom(0)->addNabor(atom2);
                    atom2->addNabor(res->atom(0));
                    res->atom(0)->addBondnumber(lig.bonds.size());
                    atom2->addBondnumber(lig.bonds.size());

                    bond.setAtom1(res->atom(0));
                    bond.setAtom2(atom2);
                    bond.getAtom1()->addType(AtomType::BONDED);
                    bond.getAtom2()->addType(AtomType::BONDED);
                    bond.type = B_NORMAL;
                    lig.bonds.push_back(bond);
                }
            }
            nth++;
        }
    } //End loop over residues

    // Add the CONECT records to the list
    for (i = 0; (unsigned int) i < lig.connects.size(); i++)
    {
        lig.Connect(lig.connects[i]);
    }
    return (lig.bonds.size());
}

bool IsPeptide(const Residue &res)
{
    // does residue have minimal requirements for peptide?
    int found = 0;
    for (int i = 0; i < res.atomCount(); i++)
    {
        if (!strcmp(res.atom(i)->name(), "C"))
        {
            found++;
        }
        if (!strcmp(res.atom(i)->name(), "O"))
        {
            found++;
        }
        if (!strcmp(res.atom(i)->name(), "OT1"))
        {
            found++;
        }
        if (!strcmp(res.atom(i)->name(), "OT2"))
        {
            found++;
        }
        if (!strcmp(res.atom(i)->name(), "OTX"))
        {
            found++;
        }
        if (!strcmp(res.atom(i)->name(), "CA"))
        {
            found++;
        }
        if (!strcmp(res.atom(i)->name(), "N"))
        {
            found++;
        }
        if (found == 4)
        {
            break;
        }
    }
    if (found == 4)
    {
        return true;
    }
    else
    {
        return false;
    }
}

bool IsNucleic(Residue *res)
{

    static char type3[][4] =
    {
        "A", "T", "U", "G", "C",
        "ADE", "THY", "URA", "GUA", "CYT"
    };
    for (unsigned int i = 0; i < sizeof(type3)/4; i++)
    {
        if (!strncmp(res->type().c_str(), type3[i], strlen(type3[i])))
        {
            /* require that it have at least P o3* o4* c4* */
            if (atom_from_name("P", *res)
                && (atom_from_name("O3'", *res) || atom_from_name("O3*", *res))
                && (atom_from_name("O4'", *res) || atom_from_name("O4*", *res))
                && (atom_from_name("C4'", *res) || atom_from_name("C4*", *res)))
            {
                res->add_linkage_type(NUCLEIC);
                return (1);
            }
            else
            {
                return (0);
            }
        }
    }
    return (0);
}

bool IsNucleic(RESIDUE *res)
{

    static char type3[][4] =
    {
        "A", "T", "U", "G", "C",
        "ADE", "THY", "URA", "GUA", "CYT"
    };
    for (unsigned int i = 0; i < sizeof(type3)/4; i++)
    {
        if (!strncasecmp(res->type().c_str(), type3[i], strlen(type3[i])))
        {
            /* require that it have at least P o3* o4* c4* */
            if (atom_from_name("P", *res)
                && (atom_from_name("O3'", *res) || atom_from_name("O3*", *res))
                && (atom_from_name("O4'", *res) || atom_from_name("O4*", *res))
                && (atom_from_name("C4'", *res) || atom_from_name("C4*", *res)))
            {
                res->add_linkage_type(NUCLEIC);
                return true;
            }
            else
            {
                return false;
            }
        }
    }
    return false;
}

bool IsWater(const RESIDUE *res)
{
    if (!strcmp(res->type().c_str(), "HOH"))
    {
        return 1;
    }
    if (!strcmp(res->type().c_str(), "WAT"))
    {
        return 1;
    }
    if (!strcmp(res->type().c_str(), "hoh"))
    {
        return 1;
    }
    if (!strcmp(res->type().c_str(), "wat"))
    {
        return 1;
    }
    return 0;
}

int IsDna(const Residue &res)     // does residue have minimal requirements for dna?
{
    int found = 0, primefound = 1;
    for (unsigned int i = 0; i < res.atoms().size(); i++)
    {
        if (!strcmp(res.atom(i)->name(), "C2*"))
        {
            found++;
        }
        if (!strcmp(res.atom(i)->name(), "C1*"))
        {
            found++;
        }
        if (!strcmp(res.atom(i)->name(), "O4*"))
        {
            found++;
        }
        if (!strcmp(res.atom(i)->name(), "O3*"))
        {
            found++;
        }
        if (!strcmp(res.atom(i)->name(), "C4*"))
        {
            found++;
        }
        if (!strcmp(res.atom(i)->name(), "C3*"))
        {
            found++;
        }
        if (found >= 5)
        {
            break;
        }
        if (!strcmp(res.atom(i)->name(), "C2\'"))
        {
            found++;
        }
        if (!strcmp(res.atom(i)->name(), "C1\'"))
        {
            found++;
        }
        if (!strcmp(res.atom(i)->name(), "O1\'"))
        {
            found++;
            primefound = 2;
        }
        if (!strcmp(res.atom(i)->name(), "C4\'"))
        {
            found++;
        }
        if (!strcmp(res.atom(i)->name(), "C3\'"))
        {
            found++;
        }
        if (found == 5)
        {
            break;
        }
    }
    if (found >= 5)
    {
        return (primefound);
    }
    else
    {
        return (0);
    }
}

int IsDna(const RESIDUE &res)     // does residue have minimal requirements for dna?
{
    int found = 0, primefound = 1;
    for (int i = 0; i < res.atomCount(); i++)
    {
        if (!strcmp(res.atom(i)->name(), "C2*"))
        {
            found++;
        }
        if (!strcmp(res.atom(i)->name(), "C1*"))
        {
            found++;
        }
        if (!strcmp(res.atom(i)->name(), "O4*"))
        {
            found++;
        }
        if (!strcmp(res.atom(i)->name(), "O3*"))
        {
            found++;
        }
        if (!strcmp(res.atom(i)->name(), "C4*"))
        {
            found++;
        }
        if (!strcmp(res.atom(i)->name(), "C3*"))
        {
            found++;
        }
        if (found >= 5)
        {
            break;
        }
        if (!strcmp(res.atom(i)->name(), "C2\'"))
        {
            found++;
        }
        if (!strcmp(res.atom(i)->name(), "C1\'"))
        {
            found++;
        }
        if (!strcmp(res.atom(i)->name(), "O1\'"))
        {
            found++;
            primefound = 2;
        }
        if (!strcmp(res.atom(i)->name(), "C4\'"))
        {
            found++;
        }
        if (!strcmp(res.atom(i)->name(), "C3\'"))
        {
            found++;
        }
        if (found == 5)
        {
            break;
        }
    }
    if (found >= 5)
    {
        return (primefound);
    }
    else
    {
        return (0);
    }
}

bool IsAminoAcid(RESIDUE *res)
{

    static char type3[][4] =
    {
        /*
         # load up the translation table before we start...
         # The upper and lower case forms are common and
         # mixed-case form is as given in Dickerson and Geis. mp
         */
        "ala", "Ala", "ALA", "A",
        "dal", "Dal", "DAL", /* d-alinine for synthetic peptides */
        "arg", "Arg", "ARG", "R",
        "asn", "Asn", "ASN", "N",
        "asp", "Asp", "ASP", "D",
        "cys", "Cys", "CYS", "C",
        "gln", "Gln", "GLN", "Q",
        "glu", "Glu", "GLU", "E",
        "gly", "Gly", "GLY", "G",
        "his", "His", "HIS", "H",
        "ile", "Ile", "ILE", "I",
        "leu", "Leu", "LEU", "L",
        "lys", "Lys", "LYS", "K",
        "met", "Met", "MET", "M",
        "phe", "Phe", "PHE", "F",
        "pro", "Pro", "PRO", "P",
        "ser", "Ser", "SER", "S",
        "thr", "Thr", "THR", "T",
        "trp", "Trp", "TRP", "W",
        "tyr", "Tyr", "TYR", "Y",
        "val", "Val", "VAL", "V",
        /*
           "sex", "Sex", "SEX", "X",
           "thx", "Thx", "THX", "X"
         */
        "mse", "Mse", "MSE", "M",
        "ptr", "Ptr", "PTR", "Y",
        "sep", "Sep", "SEP", "S",
        "tpo", "Tpo", "TPO", "T",
        "cse", "Cse", "CSE", "C",
        "soc", "Soc", "SOC", "C",
    };
    MIAtom *CA, *N, *C;
    for (unsigned int i = 0; i < sizeof(type3)/4; i++)
    {
        if (!strcmp(res->type().c_str(), type3[i]))
        {
            /* require that it have CA, N, C (O could be OT1 if C-terminal)*/
            CA = atom_from_name("CA", *res);
            N  = atom_from_name("N", *res);
            C  = atom_from_name("C", *res);
            if (CA != NULL && N != NULL && C != NULL)
            {
                res->add_linkage_type(PEPTIDE);
                return true;
            }
            else
            {
                return false;
            }
        }
    }
    return false;
}

void AminoOrNucleic(RESIDUE *res1, int reset)
{
    RESIDUE *res;
    for (res = res1; Residue::isValid(res); res = res->next())
    {
        if (reset)
        {
            res->set_linkage_type(0);
        }
        IsAminoAcid(res);
        IsNucleic(res);
    }
}

bool hbondable(const MIAtom &a1, const MIAtom &a2, const Residue &r1, const Residue &r2)
{
    /* compares a bunch of rules for h-bonding and returns
       true if an h-bond is plausible and false if its not
       note: this routine works without H-s which means its only
       approximate at best */
    // no intra-residue h-bonds
    if (&r1 == &r2)
    {
        return false;
    }
    if (&a1 == &a2)
    {
        return false;
    }
    // must be Nitrogen or Oxygen
    if (a2.name()[0] != 'N' && a2.name()[0] != 'O')
    {
        return false;
    }
    if (a1.name()[0] != 'N' && a1.name()[0] != 'O')
    {
        return false;
    }
    // first distance cut - faster than doing full distance
    float dcut = 3.30F;
    if (fabs(a1.x() - a2.x()) > dcut)
    {
        return false;
    }
    if (fabs(a1.y() - a2.y()) > dcut)
    {
        return false;
    }
    if (fabs(a1.z() - a2.z()) > dcut)
    {
        return false;
    }
    // Nitrogen and Nitrogen never form h-bonds
    if (a2.name()[0] == 'N' && a1.name()[0] == 'N')
    {
        return false;
    }
    // N never bonds N and O never bonds O
    if (!strcmp(a1.name(), "O") && !strcmp(a2.name(), "O"))
    {
        // check that both are peptides
        // i.e. water is also called "O"
        if (IsPeptide(r1) && IsPeptide(r2))
        {
            return false;
        }
    }
    if (!strcmp(a1.name(), "N") && !strcmp(a2.name(), "N"))
    {
        return false;
    }
    float d = (float)AtomDist(a1, a2);
    // second distance cut for likelies...
    if (d < 2.3F || d > 3.4F)
    {
        return false;
    }

    MIAtom *a3;
    float angle;
    // check for C-O -- X of acute angles
    if (strcmp(a1.name(), "O") == 0 && (a3 = atom_from_name("C", r1)) != NULL)
    {
        angle = (float)CalcAtomAngle(a2, a1, *a3);
        if (angle < 100.0F)
        {
            return false;
        }
        if (angle < 130.0F && d > 3.0F)
        {
            return false;
        }
    }
    if (!strcmp(a2.name(), "O") == 0 && (a3 = atom_from_name("C", r2)) != NULL)
    {
        angle = (float)CalcAtomAngle(a1, a2, *a3);
        if (angle < 100.0F)
        {
            return false;
        }
        if (angle < 130.0F && d > 3.0F)
        {
            return false;
        }
    }
    return true;
}

MIAtom *atom_from_name(const char *name, const Residue &residue)
{
    for (size_t i = 0; i < residue.atoms().size(); i++)
    {
        if (strcmp(name, residue.atom(i)->name()) == 0)
        {
            return residue.atom(i);
        }
    }
    return NULL;
}

void BisectAtom(const MIAtom *source, double *v)
{
    if (source->nabors().size() == 0)
    {
        v[0] = 1;
        v[1] = 0;
    }
    else
    {
        v[0] = source->x() - source->nabors().front()->x();
        v[1] = source->y() - source->nabors().front()->y();
        NormVect2D(v);
    }
}

RESIDUE *residue_from_name(RESIDUE *res, const char *resid, const char c)
{
    while (Residue::isValid(res))
    {
        if (!strcmp(resid, res->name().c_str()))
        {
            if (c == (char)(res->chain_id()&255) || c == '*')
            {
                return (res);
            }
        }
        res = res->next();
    }
    return (NULL);
}

static RESIDUE *LAST_RES = NULL;
static RESIDUE *LAST_RES_LIST = NULL;
void clear_residue_from_atom_cache(const Residue *res)
{

    if (res == LAST_RES)
    {
        LAST_RES = NULL;
        LAST_RES_LIST = NULL;
    }
}

RESIDUE *residue_from_atom(RESIDUE *res, MIAtom *a)
{

    // simple cache: assume most likely place for atom is last residue found,
    // or the residue after or before that, but only if using the same
    // residue list

    // if LAST_RES is not valid, clear the cache
    if (!Residue::isValid(LAST_RES))
    {
        LAST_RES = NULL;
        LAST_RES_LIST = NULL;
    }
    if (LAST_RES_LIST && LAST_RES && LAST_RES_LIST == res)
    {

        RESIDUE *r = LAST_RES;
        if (Residue::isValid(r))
        {
            for (int i = 0; i< r->atomCount(); ++i)
            {
                if (r->atom(i) == a)
                {
                    LAST_RES = r;
                    LAST_RES_LIST = res;
                    return r;
                }
            }
            r = r->next(); // nope, try next
        }

        if (Residue::isValid(r))
        {
            for (int i = 0; i< r->atomCount(); ++i)
            {
                if (r->atom(i) == a)
                {
                    LAST_RES = r;
                    LAST_RES_LIST = res;
                    return r;
                }
            }
            r = r->prev(); // nope, try prev
        }

        if (Residue::isValid(r))
        {
            for (int i = 0; i< r->atomCount(); ++i)
            {
                if (r->atom(i) == a)
                {
                    LAST_RES = r;
                    LAST_RES_LIST = res;
                    return r;
                }
            }
        }
    }

    //not in cache; do full search
    RESIDUE *r = res;
    while (Residue::isValid(r))
    {
        for (int i = 0; i < r->atomCount(); i++)
        {
            if (r->atom(i) == a)
            {
                LAST_RES = r;
                LAST_RES_LIST = res;
                return r;
            }
        }
        r = r->next();
    }

    return NULL;
}

std::vector<Residue*>::iterator
ResSearch(const MIAtom *query,
          std::vector<Residue*>::iterator res_begin,
          std::vector<Residue*>::iterator res_end)
{

    std::vector<Residue*>::iterator ri;
    unsigned int i;
    for (ri = res_begin; ri != res_end; ++ri)
    {
        Residue *res = *ri;
        for (i = 0; i < res->atoms().size(); ++i)
        {
            if (query == res->atom(i))
            {
                return ri;
            }
        }
    }
    return res_end;
}

int MaxNumBonds(const MIAtom *atom)
{
    return MaxValence(atom->atomicnumber());
}

int CurrentValence(const MIAtom &atom,
                   const std::vector<Bond> &bonds)
{

    unsigned char order;
    int nPart = 0;
    int valence = 0;
    std::vector<int>::const_iterator xBnd;
    for (xBnd = atom.bondnumbers().begin(); xBnd != atom.bondnumbers().end(); ++xBnd)
    {
        order = bonds[*xBnd].getOrder();

        if (order == PARTIALDOUBLEBOND)
        {
            nPart++;
        }
        else
        {
            valence += order;
        }
    }

    switch (nPart)
    {
    case 1: valence += 2;
        break;
    case 2: valence += 3;
        break;
    case 3: valence += 4;
        break;
    case 4: valence += 5;
        break;
    case 5: valence += 6;
        break;
    case 6: valence += 7;
        break;
    }

    return valence;

}

int UnusedValences(const MIAtom &atom,
                   const std::vector<Bond> &bonds)
{

    int v = CurrentValence(atom, bonds);

    std::vector<int> pssbl_valences = GetValenceStates(atom.atomicnumber());

    std::vector<int>::iterator i, e;
    i = pssbl_valences.begin();
    e = pssbl_valences.end();

    while (i != e)
    {
        if (v <= *i)
        {
            return *i - v;
        }
        ++i;
    }

    return MaxValence(atom.atomicnumber());
}


int CountAromaticBonds(const MIAtom &atom, const std::vector<Bond> &bonds)
{
    std::vector<int>::const_iterator i, e;
    i = atom.bondnumbers().begin();
    e = atom.bondnumbers().end();

    int n = 0;
    while (i != e)
    {
        if (bonds[*i].isaromatic == 1)
        {
            ++n;
        }
        ++i;
    }
    return n;
}

bool AtSixSixFusion(const MIAtom &atom, const std::vector<Bond> &bonds)
{
    std::vector<int>::const_iterator i, e;
    i = atom.bondnumbers().begin();
    e = atom.bondnumbers().end();

    int nSix = 0;
    while (i != e)
    {
        if (bonds[*i].smallest_aromatic_ring == 6)
        {
            ++nSix;
        }
        ++i;
    }

    return nSix == 3;
}

bool AtFiveSixFusion(const MIAtom &atom, const std::vector<Bond> &bonds)
{
    std::vector<int>::const_iterator i, e;
    i = atom.bondnumbers().begin();
    e = atom.bondnumbers().end();

    int nSix = 0;
    while (i != e)
    {
        if (bonds[*i].smallest_aromatic_ring == 6)
        {
            ++nSix;
        }
        ++i;
    }

    i = atom.bondnumbers().begin();
    int nFive = 0;
    while (i != e)
    {
        if (bonds[*i].smallest_aromatic_ring == 5)
        {
            ++nFive;
        }
        ++i;
    }

    return nFive == 2 && nSix == 1;
}

bool AtFiveFiveFusion(const MIAtom &atom, const std::vector<Bond> &bonds)
{
    std::vector<int>::const_iterator i, e;
    i = atom.bondnumbers().begin();
    e = atom.bondnumbers().end();

    int nFive = 0;
    while (i != e)
    {
        if (bonds[*i].smallest_aromatic_ring == 5)
        {
            ++nFive;
        }
        ++i;
    }

    return nFive == 3;
}

void GuessBondOrders(RESIDUE *res, vector<Bond> &bonds)
{

    Ligand mol(*res, bonds);
    mol.ClearBondOrders();
    mol.FindRingSystems();
    mol.GuessBondOrders();

    MIAtom *a1_old, *a2_old;
    MIAtom *a1_new, *a2_new;
    Bond *bond_old;
    Bond *bond_new;
    int i, j;
    unsigned int k;

    for (i = 0; i < res->atomCount(); ++i)
    {
        a1_old = res->atom(i);
        a1_new = mol.residues.front()->atom(i);
        a1_old->setHybrid(a1_new->hybrid());
        a1_old->setIsaromatic(a1_new->isaromatic());
        a1_old->setIscyclic(a1_new->iscyclic());
        a1_old->setHcount(a1_new->hcount());

        for (j = 0; j < i; ++j)
        {
            a2_old = res->atom(j);
            a2_new = mol.residues.front()->atom(j);

            bond_old = 0;
            for (k = 0; k < bonds.size(); k++)
            {
                if ((bonds[k].getAtom1() == a1_old && bonds[k].getAtom2() == a2_old)
                    || (bonds[k].getAtom1() == a2_old && bonds[k].getAtom2() == a1_old) )
                {
                    bond_old = &bonds[k];
                    break;
                }
            }

            bond_new = mol.GetBond(a1_new, a2_new);

            if (bond_old != 0 && bond_new != 0)
            {
                bond_old->setOrder(bond_new->getOrder());
                bond_old->isaromatic = bond_new->isaromatic;
                bond_old->iscyclic = bond_new->iscyclic;
            }
        }
    }
}

void HybridizeTerminalAtom(MIAtom &atom, std::vector<Bond> &bonds)
{
    Bond &bond = bonds[atom.bondnumbers().front()];

    unsigned int order;
    order = PredictBondOrder(atom.atomicnumber(),
                             atom.nabors().front()->atomicnumber(),
                             (float)AtomDist(atom, *atom.nabors()[0]));
    if (order == PARTIALDOUBLEBOND)
    {
        bond.setOrder(DOUBLEBOND);
    }
    else
    {
        bond.setOrder(order);
    }

    switch (order)
    {
    case TRIPLEBOND: atom.setHybrid(1);
        break;
    case DOUBLEBOND: atom.setHybrid(2);
        break;
    case PARTIALDOUBLEBOND: atom.setHybrid(2);
        break;
    case SINGLEBOND: atom.setHybrid(3);
        break;
    }
}

float AverageBondAngle(MIAtom &atom)
{
    MIAtom_const_iter n1, n2, ne;
    ne = atom.nabors().end();
    n1 = atom.nabors().begin();

    int n_angles = 0;
    float sum_angles = 0.0F;

    for (; n1 != ne; ++n1)
    {
        for (n2 = n1; n2 != ne; ++n2)
        {
            if (n1 == n2)
            {
                continue;
            }
            sum_angles += (float)CalcAtomAngle(**n1, atom, **n2);
            n_angles++;
        }
    }

    return sum_angles / n_angles;
}

void HybridFromGeom(MIAtom &atom)
{
    MIAtom_const_iter n1, n2, ne;
    ne = atom.nabors().end();
    n1 = atom.nabors().begin();

    int n_angles = 0;
    float sum_angles = 0.0F;

    for (; n1 != ne; ++n1)
    {
        for (n2 = n1; n2 != ne; ++n2)
        {
            if (n1 == n2)
            {
                continue;
            }
            sum_angles += (float)CalcAtomAngle(**n1, atom, **n2);
            n_angles++;
        }
    }

    float average_angle = sum_angles / n_angles;
    if (average_angle > 155)
    {
        atom.setHybrid(1);
    }
    else if (average_angle > 116)
    {
        atom.setHybrid(2);
    }
    else
    {
        atom.setHybrid(3);
    }
}

int PredictValenceFrom3D(MIAtom &atom)
{
    MIAtom_const_iter a, ae;

    unsigned char order;
    int valence = 0;
    int nPart = 0;

    ae = atom.nabors().end();
    a = atom.nabors().begin();

    while (a != ae)
    {
        order = PredictBondOrder(atom.atomicnumber(),
                                 (*a)->atomicnumber(),
                                 (float)AtomDist(atom, *atom.nabors()[0]));

        if (order == PARTIALDOUBLEBOND)
        {
            nPart++;
        }
        else
        {
            valence += order;
        }

        ++a;
    }

    if (nPart > 0)
    {
        valence += nPart + 1;
    }

    return valence;
}

static char three[][4] = {"ALA", "CYS", "ASP", "GLU", "PHE",
                          "GLY", "HIS", "ILE", "LYS", "LEU",
                          "MET", "ASN", "PRO", "GLN", "ARG",
                          "SER", "THR", "VAL", "TRP", "TYR",
                          "WAT", "HOH", "HEM", "CU ", "ZN ",
                          "MSE", "PTR", "SEP", "TPO", "CSE", "SOC", };
static char one[] = {'A', 'C', 'D', 'E', 'F',
                     'G', 'H', 'I', 'K', 'L',
                     'M', 'N', 'P', 'Q', 'R',
                     'S', 'T', 'V', 'W', 'Y',
                     'O', 'O', 'X', 'U', 'Z',
                     'M', 'Y', 'S', 'T', 'C', 'C'};

char singleletter(const char *t)
{
    // convert three letter residue names to single letter names
    for (unsigned int i = 0; i < sizeof(three)/4; i++)
    {
        if (!strcmp(t, three[i]))
        {
            return (one[i]);
        }
    }
    return ('?');
}

const char *tripleletter(const char t)
{
    // convert three letter residue names to single letter names
    for (unsigned int i = 0; i < sizeof(one); i++)
    {
        if (t == one[i])
        {
            return (three[i]);
        }
    }
    return (NULL);
}

int BuildCALinks(vector<Bond> &bonds, const RESIDUE *res)
{
    MIAtom *a1, *a2;
    Bond bond;
    const RESIDUE *firstres = res, *nextres;
    int nlinks = 0;
    while (Residue::isValid(res) && (nextres = res->next()) != NULL)
    {
        if (res->chain_id() != nextres->chain_id())
        {
            nextres = firstres;
            firstres = res->next();
        }

        if ((a1 = atom_from_name("CA", *res)) != NULL
            && (a2 = atom_from_name("CA", *nextres)) != NULL && AtomDist(*a1, *a2) < 4.5)
        {
            bond.setAtom1(a1);
            bond.setAtom2(a2);
            bond.type = B_LINK;
            bonds.push_back(bond);
            nlinks++;
        }

        res = res->next();
    }
    return nlinks;
}

bool IsPolarH(const MIAtom *atom, const RESIDUE *res)
{
    static char hnames[][5] =
    {
        "HE", "ARG",
        "HH11", "ARG",
        "HH12", "ARG",
        "HH21", "ARG",
        "HH22", "ARG",
        "HD21", "ASN",
        "HD22", "ASN",
        "HE21", "GLN",
        "HE22", "GLN",
        "HD1",  "HIS",
        "HE2",  "HIS",
        "HZ1",  "LYS",
        "HZ2",  "LYS",
        "HZ3",  "LYS",
        "HG", "SER",
        "HG1",  "THR",
        "HE1",  "TRP",
        "HH", "TYR"
    };

    if (strcmp(atom->name(), "H") == 0)
    {
        return 1;
    }
    // a shelx name for the same hydrogen
    if (strcmp(atom->name(), "H0") == 0)
    {
        return 1;
    }
    else
    {
        for (unsigned int i = 0; i < sizeof(hnames)/5; i += 2)
        {
            if (strcmp(atom->name(), hnames[i]) == 0 && strcmp(res->type().c_str(), hnames[i+1]) == 0)
            {
                return 1;
            }
        }
    }
    return 0;
}

/////////////////////////////////////////////////////////////////////////////
// Function:  CopyCoords
// Purpose:   Copies the (real-space, not screen) coordinates from the atoms
//        in one residue to the corresponding atoms in another residue,
//        overwriting the x, y, and z fields of the atoms in the target
//        RESIDUE
// Input:       Two residue ptrs
// Output:    True if all the atoms in the target have corresponding atoms
//        in the source residue
// Requires:  That the atom names correspond...this comparison is case-sensitive
//        to support the (unlikely) possibility that the names of two
//        atoms differ only in case.
//        Also note, if this returns false, the coordinates of the
//        target residue are now garbage.
/////////////////////////////////////////////////////////////////////////////
bool MICopyCoords(RESIDUE *target, const RESIDUE *source)
{
    int i;
    if (!Residue::isValid(source))
    {
        return false;
    }

    std::map<std::string, const MIAtom*> atm_map;
    for (i = 0; i < source->atomCount(); ++i)
    {
        atm_map.insert(pair<std::string, const MIAtom*> (source->atom(i)->name(), source->atom(i)));
    }

    std::map<std::string, const MIAtom*>::iterator a;
    for (i = 0; i < target->atomCount(); ++i)
    {
        a = atm_map.find(std::string(target->atom(i)->name()));
        if (a != atm_map.end())
        {
            target->atom(i)->copyPosition(*a->second);
        }
        else
        {
            return false;
        }
    }
    return true;
}

}   //namespace chemlib


