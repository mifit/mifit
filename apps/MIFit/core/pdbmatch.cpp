#include <chemlib/chemlib.h>
#include <chemlib/Monomer.h>
#include "corelib.h"

#include "rotlsq.h"
#include "RESIDUE.h"

#define X 0
#define Y 1
#define Z 2

using namespace chemlib;

//@{
// A hit for the LSQ function.
//@}
class hit
{
public:
    int sum;
    int astart, aend;
    char rchain;
    char file[1024];
};

int atomvector(const MIAtom *atom1, const MIAtom *atom2)
{
    float d = 0;
    float dx, dy, dz;
    if (!atom1 || !atom2)
    {
        return (0);
    }
    dx = atom1->x() - atom2->x();
    dy = atom1->y() - atom2->y();
    dz = atom1->z() - atom2->z();
    d = (float)sqrt(dx*dx+dy*dy+dz*dz);
    d *= 100;
    return ((int)d);
}

int hcompare(const void *i, const void *j)
{
    hit *hj = (hit*)j;
    hit *hi = (hit*)i;
    return (hj->sum - hi->sum);
}


static void copya(double a[3], MIAtom *atom)
{
    a[0] = atom->x();
    a[1] = atom->y();
    a[2] = atom->z();
}

Residue *pdbvec(std::string &pentdir, const MIAtomList &CA, const MIAtomList &CB)
{
    FILE *fvec;
    std::string vectorfile = pentdir + "/pdbvec.list";
    if ((fvec = fopen(vectorfile.c_str(), "r")) == NULL)
    {
        Logger::message("Error: Unable to find pdbvec.list in data directory - incorrect installation?");
        return NULL;
    }

    Residue *hit = matchvec(CA, CB, pentdir, fvec);
    fclose(fvec);
    if (!hit)
    {
        return NULL;
    }

    // least-squares fit hit onto atom vectors
    double (*a)[3], (*b)[3];
    double *w;
    double r[3][3];
    double v[3];
    MIAtom *atom;
    Residue *sres;
    std::vector<float> x_target;
    std::vector<float> y_target;
    std::vector<float> z_target;
    int n = CA.size(), i, j;

    for (i = 0; (unsigned int)i < CB.size(); i++)
    {
        if (CB[i])
        {
            n++;
        }
    }

    if (n < 3)
    {
        Logger::message("Must have at least 3 matching pairs of atoms to calculate least-quares fit");
        return NULL;
    }
    a = new double[n][3];
    b = new double[n][3];
    w = new double[n];
    if (!w)
    {
        Logger::message("Can't allocate memory - try fewer matches");
        if (a)
        {
            delete[] a;
        }
        if (b)
        {
            delete[] b;
        }
        return NULL;
    }
    sres = hit;
    n = 0;
    for (j = 0; (unsigned int)j < CA.size(); j++)
    {
        if ((atom = atom_from_name("CA", *sres)) != NULL)
        {
            copya(b[n], CA[j]);
            x_target.push_back(CA[j]->x());
            y_target.push_back(CA[j]->y());
            z_target.push_back(CA[j]->z());
            copya(a[n], atom);
            w[n] = 1.0;
            n++;
            sres = sres->next();
        }
    }
    sres = hit;
    for (j = 0; (unsigned int)j < CB.size(); j++)
    {
        if (CB[j] && (atom = atom_from_name("CB", *sres)) != NULL)
        {
            copya(b[n], CA[j]);
            x_target.push_back(CA[j]->x());
            y_target.push_back(CA[j]->y());
            z_target.push_back(CA[j]->z());
            copya(a[n], atom);
            w[n] = 0.5;
            n++;
            sres = sres->next();
        }
    }
    if (n < 3)
    {
        Logger::message("Must have at least 3 matching pairs of atoms to calculate least-quares fit");
        if (a)
        {
            delete[] a;
        }
        if (b)
        {
            delete[] b;
        }
        if (w)
        {
            delete[] w;
        }
        return 0;
    }
    rotlsqfit(a, b, w, n, r, v);

    double x, y, z, tx, ty, tz, sum = 0.0, d;

    sres = hit;
    n = 0;
    for (j = 0; (unsigned int)j < CA.size(); j++)
    {
        if ((atom = atom_from_name("CA", *sres)) != NULL)
        {
            tx = atom->x();
            ty = atom->y();
            tz = atom->z();
            x = r[0][0]*tx + r[0][1]*ty + r[0][2]*tz + v[0];
            y = r[1][0]*tx + r[1][1]*ty + r[1][2]*tz + v[1];
            z = r[2][0]*tx + r[2][1]*ty + r[2][2]*tz + v[2];
            tx = x - x_target[n];
            ty = y - y_target[n];
            tz = z - z_target[n];
            d = tx*tx + ty*ty + tz*tz;
            sum += d;
            n++;
            sres = sres->next();
        }
    }
    sres = hit;
    for (j = 0; (unsigned int)j < CB.size(); j++)
    {
        if (CB[j] && (atom = atom_from_name("CB", *sres)) != NULL)
        {
            tx = atom->x();
            ty = atom->y();
            tz = atom->z();
            x = r[0][0]*tx + r[0][1]*ty + r[0][2]*tz + v[0];
            y = r[1][0]*tx + r[1][1]*ty + r[1][2]*tz + v[1];
            z = r[2][0]*tx + r[2][1]*ty + r[2][2]*tz + v[2];
            tx = x - x_target[n];
            ty = y - y_target[n];
            tz = z - z_target[n];
            d = tx*tx + ty*ty + tz*tz;
            sum += d;
            n++;
            sres = sres->next();
        }
    }

    sum = sqrt(sum/(double)n);

    Logger::log("RMS Distances of matches = %0.2f\n", (float)sum);

    sres = hit;
    while (Monomer::isValid(sres))
    {
        for (i = 0; i < sres->atomCount(); i++)
        {
            tx = sres->atom(i)->x();
            ty = sres->atom(i)->y();
            tz = sres->atom(i)->z();
            x = r[0][0]*tx + r[0][1]*ty + r[0][2]*tz + v[0];
            y = r[1][0]*tx + r[1][1]*ty + r[1][2]*tz + v[1];
            z = r[2][0]*tx + r[2][1]*ty + r[2][2]*tz + v[2];
            sres->atom(i)->setPosition(x, y, z);
            if (sres->atom(i)->name()[0] == 'C')
            {
                sres->atom(i)->setColor(Colors::MAGENTA);
            }
            else
            {
                sres->atom(i)->setColor(color_by_name(sres->atom(i)->name()));
            }
        }
        sres = sres->next();
    }
    delete[] a;
    delete[] b;
    delete[] w;

    return hit;
}

static MIAtom *scan_atom(const char *linebuf)
{
#define column(n) (linebuf+(n)-1)
    MIAtom *atom = new MIAtom;
    double x, y, z;
    double Bvalue;
    double occupancy;
    char abuf[chemlib::MAXATOMNAME];
    int lb = strlen(linebuf);
    if (lb < 53)
    {
        delete atom;
        return NULL;
    }
    strncpy(abuf, column(13), 5);
    abuf[5] = '\0';
    int i = 4;
    while (abuf[i] == ' ' && i > 0)
    {
        abuf[i] = '\0';
        i--;
    }
    int ii = 0;
    while (i != 0 && abuf[ii] == ' ' && ii <= 4)
    {
        ii++;
    }
    atom->setName(abuf+ii);
    atom->setColor(color_by_name(atom->name()));
    if (atom->color() == Colors::BLACK)
    {
        atom->setColor(Colors::WHITE);
    }
    x  = atof(column(30));
    y  = atof(column(39));
    z  = atof(column(47));
    // catch nonsense atoms sometimes added by XPLOR
    if (x >= 9999.0 || y >= 9999.0 || z >= 9999.0)
    {
        delete atom;
        return NULL;
    }
    atom->setPosition((float)x, (float)y, (float)z);
    if (lb > 55)
    {
        occupancy = atof(column(55));
    }
    else
    {
        occupancy = 1.0;
    }
    atom->setOcc((float)occupancy);
    if (lb > 61)
    {
        Bvalue = atof(column(61));
    }
    else
    {
        Bvalue = 15.0;
    }
    if (Bvalue < 0.0)
    {
        Bvalue = 15.0;
    }
    atom->setBValue((float)Bvalue);
    atom->set_radius_type(0);
    atom->setType(MIAtom::MIGetAtomTypeFromName(atom->name()));
    atom->setMass(atoi(column(7)));   // temporary storage for atomnumber
    // used for figuring CONECT info later on
    char chainid = *column(22);
    atom->setAltloc(linebuf[16]);
    strncpy(abuf, column(13), 2);
    atom->setAtomicnumber(Atomic_Number(abuf));
    atom->setAtomnumber(chainid);
    atom->setSymmop(0);

    return atom;
#undef column
}

#define AMINOACID 1
#define LINKED 2
#define LINKEDAA 3
#define MAXHITS 1
#define column(n) (buf+(n)-1)

Residue *matchvec(const MIAtomList &CA, const MIAtomList &CB, std::string &pentdir, FILE *fvec)
{
    FILE *fpdb;
    int nvec = 0, i, n, d, best = 9999999;
    //int noCB = false;
    MIAtom *a1, *a2, *a3, *a4, *a5;
    MIAtom *b1, *b2, *b3, *b4, *b5;
    char file[1024];
    hit top[MAXHITS];
    char buf[1024];
    char seq[10];
    int found;
    int d1 = 0, d2 = 0, d3 = 0, d4 = 0, d5 = 0, d6 = 0;
    int d7 = 0, d8 = 0, d9 = 0, d10 = 0, d11 = 0, d12 = 0;
    int v1 = 0, v2 = 0, v3 = 0, v4 = 0, v5 = 0, v6 = 0;
    int v7 = 0, v8 = 0, v9 = 0, v10 = 0, v11 = 0, v12 = 0;
    int is;
    char rchain, chain;
    int atomnumber, astart, aend;

    if (CA.size() < 4 || CB.size() < 4)
    {
        Logger::log("matchvec: error: need at least 5 atom positions");
        return NULL;
    }

    for (i = 0; i < MAXHITS; i++)
    {
        top[i].sum = 9999999;
    }

    a1 = CA[0];
    a2 = CA[1];
    a3 = CA[2];
    a4 = CA[3];
    a5 = CA[4];
    if (!a1 || !a2 || !a3 || !a4 || !a5)
    {
        Logger::log("Can't find 5 Calpha's");
        return (0);
    }
    b1 = CB[0];
    b2 = CB[1];
    b3 = CB[2];
    b4 = CB[3];
    b5 = CB[4];
    /* if no CB's at all then we switch how we do match */
    //if( !b1 && !b2 && !b3 && !b4 && !b5)noCB = true;
    /* if no CB (i.e. GLY) use CA instead */
    //if(b1 == NULL)b1=a1;
    //if(b2 == NULL)b2=a2;
    //if(b3 == NULL)b3=a3;
    //if(b4 == NULL)b4=a4;
    //if(b5 == NULL)b5=a5;
    d1 = atomvector(a1, a3);
    d2 = atomvector(a1, a4);
    d3 = atomvector(a1, a5);
    d4 = atomvector(a2, a4);
    d5 = atomvector(a2, a5);
    d6 = atomvector(a3, a5);
    //if( !noCB){
    if (b1 && b2)
    {
        d7 = atomvector(b1, b2);
    }
    if (b5 && b3)
    {
        d8 = atomvector(b5, b3);
    }
    if (b2 && a4)
    {
        d9 = atomvector(b2, a4);
    }
    if (b3 && a1)
    {
        d10 = atomvector(b3, a1);
    }
    if (b3 && a5)
    {
        d11 = atomvector(b3, a5);
    }
    if (b4 && a2)
    {
        d12 = atomvector(b4, a2);
    }
    //}
    while (fgets(buf, 1024, fvec) != NULL)
    {
        if (!strncmp(buf, "START", 5))
        {
            sscanf(buf, "%*s%s", file);
            continue;
        }
        sscanf(buf, "%d%d%d%d%d%d%d%d%d%d%d%d%d", &n, &v1, &v2, &v3, &v4, &v5, &v6, &v7, &v8, &v9, &v10, &v11, &v12);
        d =  (d1-v1)*(d1-v1);
        d += (d2-v2)*(d2-v2);
        d += (d3-v3)*(d3-v3);
        d += (d4-v4)*(d4-v4);
        d += (d5-v5)*(d5-v5);
        int n = 5;
        //if( !noCB){
        if (d6 != 0)
        {
            d += (d6-v6)*(d6-v6);
            n++;
        }
        /* non-use of d7 is deliberate */
        if (d8 != 0)
        {
            d += (d8-v8)*(d8-v8);
            n++;
        }
        if (d9 != 0)
        {
            d += (d9-v9)*(d9-v9);
            n++;
        }
        if (d10 != 0)
        {
            d += (d10-v10)*(d10-v10);
            n++;
        }
        if (d11 != 0)
        {
            d += (d11-v11)*(d11-v11);
            n++;
        }
        if (d12 != 0)
        {
            d += (d12-v12)*(d12-v12);
            n++;
        }
        d /= n;
        //}
        /* if no CB's then
         * calculate gly,pro penalty see if still best */
        sscanf(buf, "%*d%*d%*d%*d%*d%*d%*d%*d%*d%*d%*d%*d%*d%*d%*d %*c %s",
               seq);
        for (is = 0; is < 5; is++)
        {
            if (seq[is] == 'G' || seq[is] == 'P')
            {
                d +=  40;
            }
        }
        if (d < best)
        {
            sscanf(buf, "%*d%*d%*d%*d%*d%*d%*d%*d%*d%*d%*d%*d%*d%d%d %c %s", &astart, &aend, &rchain, seq);
            top[0].sum = d;
            top[0].astart = astart;
            top[0].aend = aend;
            top[0].rchain = rchain;
            strcpy(top[0].file, file);
            qsort(top, MAXHITS, sizeof(hit), hcompare);
            best = top[0].sum;
        }
        nvec++;
    }
    sprintf(buf, "Checked %d pentamer vectors", nvec);
    Logger::log(buf);
    found = 0;
    /* write out in reverse order so that best hit is first in xfit
     * for(i=0;i<MAXHITS;i++){
     */
    n = 0;
    Residue *res = (Residue*)NULL, *reslist = 0;
    MIAtom *atom;
    MIAtomList atoms;
    for (i = MAXHITS-1; i >= 0; i--)
    {
        /* now find this sequence and write it out to stdout */
        char type[20], rname[20], aname[20];
        strcpy(buf, pentdir.c_str());
        strcat(buf, "/");
        strcat(buf, top[i].file);
        fpdb = fopen(buf, "r");
        if (!fpdb)
        {
            Logger::log("cannot open pdb source file");
            return (0);
        }
        if (top[i].rchain == '*')
        {
            top[i].rchain = ' ';
        }
        atoms.clear();
        char oldname[20], oldtype[20];
        oldname[0] = '\0';
        while (fgets(buf, 1024, fpdb) != NULL)
        {
            if (strncmp("ATOM", buf, 4))
            {
                continue;
            }
            sscanf(column(7), "%d", &atomnumber);
            chain = *column(22);
            sscanf(column(18), "%3s", type);
            sscanf(column(23), "%4s", rname);
            sscanf(column(13), "%4s", aname);
            if (atoms.size() > 0)
            {
                // first residue
                if (oldname[0] == '\0')
                {
                    strcpy(oldname, rname);
                }
                if (strcmp(oldname, rname) != 0)
                {
                    if (reslist)
                    {
                        res = res->insertResidue(make_res(atoms));
                        res->set_linkage_type(MIDDLE);
                    }
                    else
                    {
                        res = reslist = make_res(atoms);
                        res->set_linkage_type(NTERMINUS);
                    }
                    res->setName(std::string(oldname, 4));
                    res->setType(std::string(oldtype, 4));
                    strcpy(oldname, rname);
                    res->setSecstr('U');
                    atoms.clear();
                }
            }
            if (atomnumber <= top[i].aend && atomnumber >= top[i].astart && chain == top[i].rchain)
            {
                if (!(strcmp(type, "PRO") == 0 || strcmp(type, "GLY") == 0))
                {
                    if (strcmp(aname, "CA") == 0 || strcmp(aname, "CB") == 0
                        || strcmp(aname, "O") == 0 || strcmp(aname, "C") == 0
                        || strcmp(aname, "N") == 0)
                    {
                        buf[17] = 'A';
                        buf[18] = 'L';
                        buf[19] = 'A';
                        buf[21] = '1'+n;
                        if (buf[strlen(buf)-1] == '\n')
                        {
                            buf[strlen(buf)-1] = '\0';
                        }
                        atom = scan_atom(buf);
                        sscanf(column(18), "%3s", oldtype);
                        if (atom)
                        {
                            atoms.push_back(atom);
                        }
                        found++;
                    }
                }
                else
                {
                    buf[21] = '1'+n;
                    if (buf[strlen(buf)-1] == '\n')
                    {
                        buf[strlen(buf)-1] = '\0';
                    }
                    atom = scan_atom(buf);
                    sscanf(column(18), "%3s", oldtype);
                    if (atom)
                    {
                        atoms.push_back(atom);
                    }
                    found++;
                }
            }
        }
        // finish last residue
        if (atoms.size() > 0)
        {
            res = res->insertResidue(make_res(atoms));
            res->setName(std::string(oldname, 4));
            res->setType(std::string(oldtype, 4));
            strcpy(oldname, rname);
            res->setSecstr('U');
            res->set_linkage_type(CTERMINUS);
        }
        fclose(fpdb);
        n++;
    }
    printf("END\n");
    /* for(i=0;i<MAXHITS;i++){ */
    for (i = MAXHITS-1; i >= 0; i--)
    {
        sprintf(buf, "%3d: hit=%5d start=%5d end=%5d %c %s", i, top[i].sum, top[i].astart, top[i].aend, top[i].rchain, top[i].file);
        Logger::log(buf);
    }
    sprintf(buf, "Found %d atoms", found);
    Logger::log(buf);
    return (reslist);
}

