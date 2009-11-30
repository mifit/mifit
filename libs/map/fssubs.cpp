#include <cstdio>
#include <cstring>
#include <cmath>

#include <nongui/nonguilib.h>
#include <math/mathlib.h>

#include "fssubs.h"
#include "CMapHeaderBase.h"

/* readin a binary integer with conditional byte-swabbing */
int fssubs_in_int(FILE *fp, unsigned int swab)
{
    unsigned char in[4];
    fread(&(in[0]), 1, 4, fp);
    if (!swab)
    {
        return (int)((unsigned int)in[0]+((unsigned int)(in[1])<<8)+((unsigned int)(in[2])<<16)+((unsigned int)(in[3])<<24));
    }
    else
    {
        return (int)(((unsigned int)(in[0])>>24)+((unsigned int)(in[1])>>16)+((unsigned int)(in[2])>>8)+((unsigned int)(in[3])));
    }
}

/* readin a binary float with conditional byte-swabbing */
float fssubs_in_float(FILE *fp, unsigned int swab)
{
    unsigned char in[4];
    union
    {
        unsigned char buf[4];
        float f;
    } t;
    fread(&(in[0]), 1, 4, fp);
    if (!swab)
    {
        t.buf[0] = in[0];
        t.buf[1] = in[1];
        t.buf[2] = in[2];
        t.buf[3] = in[3];
    }
    else
    {
        t.buf[0] = in[3];
        t.buf[1] = in[2];
        t.buf[2] = in[1];
        t.buf[3] = in[0];
    }
    return t.f;
}

int
fsread(int *rho, int *nx, int *ny, int *nz, float *rhoscale, float *rmsrho, FILE *fp, char *err)
{
    int irow, icol, iplane, ii, ix, iy, iz, i;
    int nrow, ncol, nplane, n;
    int nsym, norn, index;
    int header, trailer;
    int sint = sizeof(int), irho;
    float rhomin = 9999., rhomax = (-9999.);
    int buf[200];

    /* read first record itle(20),junk,para(6),nsym*/
    fread(&header, sizeof(int), 1, fp);
    n = 20 + 1 + 6 ;
    fread(buf, sizeof(int), n, fp);
    header = header - sizeof(int)*n;
    fread(&nsym, sizeof(int), 1, fp);
    header = header -sizeof(int);
    /* go to  next record */
    while (header)
    {
        if (((int)fread(buf, 1, 1, fp)) <= 0)
        {
            Logger::log("fsread: premature end-of-file");
            return (0);
        }
        header--;
    }
    fread(&trailer, sizeof(int), 1, fp);
    /* there are nsym records with 9 numbers next */
    for (i = 0; i < nsym; i++)
    {
        fread(&header, sizeof(header), 1, fp);
        while (header)
        {
            if (((int)fread(buf, 1, 1, fp)) <= 0)
            {
                Logger::log("fsread: premature end-of-file");
                return (0);
            }
            header--;
        }
        fread(&trailer, sizeof(int), 1, fp);
    }
    /* read dimension info */
    fread(&header, sizeof(header), 1, fp);
    fread(&nrow, sizeof(int), 1, fp);
    header = header - sizeof(int);
    fread(&ncol, sizeof(int), 1, fp);
    header = header - sizeof(int);
    fread(&nplane, sizeof(int), 1, fp);
    header = header - sizeof(int);
    fread(buf, sizeof(int), 1, fp);
    header = header - sizeof(int);
    fread(&norn, sizeof(int), 1, fp);
    header = header - sizeof(int);
    /* go to  next record */
    while (header)
    {
        if (((int)fread(buf, 1, 1, fp)) <= 0)
        {
            Logger::log("fsread: premature end-of-file");
            return (0);
        }
        header--;
    }
    if (((int)fread(buf, sizeof(int), 1, fp)) <= 0)
    {
        Logger::log("fsread: premature end-of-file");
        return (0);
    }

    if (norn < 0 || norn > 2)
    {
        sprintf(err, "Error in reading fsfor format map file\n");
        return (0);
    }
    if (norn == 0)
    {
        Logger::log("Map in rows of x, columns in z, planes of y");
        *nx = nrow;
        *nz = ncol;
        *ny = nplane;
    }
    else if (norn == 1)
    {
        Logger::log("Map in rows of y, columns in z, planes of x");
        *ny = nrow;
        *nz = ncol;
        *nx = nplane;
    }
    else if (norn == 2)
    {
        Logger::log("Map in rows of x, columns in y, planes of z");
        *nx = nrow;
        *ny = ncol;
        *nz = nplane;
    }
    else
    {
        Logger::log("Illegal sort order - error in map file");
        return (0);
    }

    ii = 0;
    *rmsrho = 0.0;
    for (iplane = 0; iplane < nplane; iplane++)
    {
        for (icol = 0; icol < ncol; icol++)
        {
            fread(buf, sint, 1, fp);
            for (irow = 0; irow < nrow; irow++)
            {
                ii++;
                fread(&irho, sint, 1, fp);
                if (norn == 0)
                {
                    ix = irow;
                    iy = iplane;
                    iz = icol;
                }
                if (norn == 1)
                {
                    ix = iplane;
                    iy = irow;
                    iz = icol;
                }
                if (norn == 2)
                {
                    ix = irow;
                    iy = icol;
                    iz = iplane;
                }
                index = *nx*(*ny*iz+iy) + ix;
                rho[index] = irho;
                *rmsrho = *rmsrho + irho*irho;
                if (irho < rhomin)
                {
                    rhomin = (float)irho;
                }
                if (irho > rhomax)
                {
                    rhomax = (float)irho;
                }
            }
            fread(buf, sint, 1, fp);
        }
    }
    *rmsrho = (float)sqrt((double)*rmsrho/(double)ii);
    *rhoscale = 1.0;
    sprintf(err, "Total points = %d  Rhomin = %f Rhomax = %f Rmsrho = %f\n", ii, rhomin, rhomax, *rmsrho);
    return (ii);
}

int
fsread_uni(std::vector<int> &rho, int *nx, int *ny, int *nz, float *rhoscale, float *rmsrho, FILE *fp, char *err, int swab)
{
    int irow, icol, iplane, ii, ix, iy, iz, i;
    int nrow, ncol, nplane, n;
    int nsym, norn, index;
    int header, trailer;
    int sint = sizeof(int), irho;
    float rhomin = 9999., rhomax = (-9999.);
    int buf[200];

    /* read first record itle(20),junk,para(6),nsym*/
    header = fssubs_in_int(fp, swab);
    n = 20 + 1 + 6 ;
    fread(buf, sizeof(int), n, fp);
    header = header - sizeof(int)*n;
    nsym = fssubs_in_int(fp, swab);
    header = header -sizeof(int);
    /* go to  next record */
    while (header)
    {
        if (((int)fread(buf, 1, 1, fp)) <= 0)
        {
            Logger::log("fsread: premature end-of-file");
            return (0);
        }
        header--;
    }
    fread(&trailer, sizeof(int), 1, fp);
    /* there are nsym records with 9 numbers next */
    for (i = 0; i < nsym; i++)
    {
        header = fssubs_in_int(fp, swab);
        while (header)
        {
            if (((int)fread(buf, 1, 1, fp)) <= 0)
            {
                Logger::log("fsread: premature end-of-file");
                return (0);
            }
            header--;
        }
        fread(&trailer, sizeof(int), 1, fp);
    }
    /* read dimension info */
    header = fssubs_in_int(fp, swab);
    nrow = fssubs_in_int(fp, swab);
    header = header - sizeof(int);
    ncol = fssubs_in_int(fp, swab);
    header = header - sizeof(int);
    nplane = fssubs_in_int(fp, swab);
    header = header - sizeof(int);
    fread(buf, sizeof(int), 1, fp);
    header = header - sizeof(int);
    norn = fssubs_in_int(fp, swab);
    header = header - sizeof(int);
    /* go to  next record */
    while (header)
    {
        if (((int)fread(buf, 1, 1, fp)) <= 0)
        {
            Logger::log("fsread: premature end-of-file");
            return (0);
        }
        header--;
    }
    if (((int)fread(buf, sizeof(int), 1, fp)) <= 0)
    {
        Logger::log("fsread: premature end-of-file");
        return (0);
    }

    if (norn < 0 || norn > 2)
    {
        sprintf(err, "Error in reading fsfor format map file\n");
        return (0);
    }
    if (norn == 0)
    {
        Logger::log("Map in rows of x, columns in z, planes of y");
        *nx = nrow;
        *nz = ncol;
        *ny = nplane;
    }
    else if (norn == 1)
    {
        Logger::log("Map in rows of y, columns in z, planes of x");
        *ny = nrow;
        *nz = ncol;
        *nx = nplane;
    }
    else if (norn == 2)
    {
        Logger::log("Map in rows of x, columns in y, planes of z");
        *nx = nrow;
        *ny = ncol;
        *nz = nplane;
    }
    else
    {
        Logger::log("Illegal sort order - error in map file");
        return (0);
    }

    ii = 0;
    *rmsrho = 0.0;
    for (iplane = 0; iplane < nplane; iplane++)
    {
        for (icol = 0; icol < ncol; icol++)
        {
            fread(buf, sint, 1, fp);
            for (irow = 0; irow < nrow; irow++)
            {
                ii++;
                irho = fssubs_in_int(fp, swab);
                if (norn == 0)
                {
                    ix = irow;
                    iy = iplane;
                    iz = icol;
                }
                if (norn == 1)
                {
                    ix = iplane;
                    iy = irow;
                    iz = icol;
                }
                if (norn == 2)
                {
                    ix = irow;
                    iy = icol;
                    iz = iplane;
                }
                index = *nx*(*ny*iz+iy) + ix;
                rho[index] = irho;
                *rmsrho = *rmsrho + irho*irho;
                if (irho < rhomin)
                {
                    rhomin = (float)irho;
                }
                if (irho > rhomax)
                {
                    rhomax = (float)irho;
                }
            }
            fread(buf, sint, 1, fp);
        }
    }
    *rmsrho = (float)sqrt((double)*rmsrho/(double)ii);
    *rhoscale = 1.0;
    sprintf(err, "Total points = %d  Rhomin = %f Rhomax = %f Rmsrho = %f\n", ii, rhomin, rhomax, *rmsrho);
    return (ii);
}

int
fssize(FILE *fp)
{
    /*  reads in the size info from a fsfour format
        fortran unformatted map file see fsfor.f */
    int nrow, ncol, nplane, n;
    int header, trailer;
    int nsym, i;
    int buf[200];

    /* read first record itle(20),junk,para(6),nsym*/
    fread(&header, sizeof(int), 1, fp);
    n = 20 + 1 + 6 ;
    fread(buf, sizeof(int), n, fp);
    header = header - sizeof(int)*n;
    fread(&nsym, sizeof(int), 1, fp);
    header = header -sizeof(int);
    if (nsym < 0 || nsym > (int)MISymmop::MAXSYMMOPS)
    {
        fprintf(stderr, "fssize: impossible nsym\n");
        return (0);
    }
    /* go to  next record */
    while (header)
    {
        fread(buf, 1, 1, fp);
        header--;
    }
    fread(&trailer, sizeof(int), 1, fp);
    /* there are nsym records next */
    for (i = 0; i < nsym; i++)
    {
        fread(&header, sizeof(int), 1, fp);
        while (header)
        {
            fread(buf, 1, 1, fp);
            header--;
        }
        fread(&trailer, sizeof(int), 1, fp);
    }
    /* read dimension info */
    fread(&header, sizeof(int), 1, fp);
    fread(&nrow, sizeof(int), 1, fp);
    fread(&ncol, sizeof(int), 1, fp);
    if (fread(&nplane, sizeof(int), 1, fp) == 1)
    {
        return (ncol*nplane*nrow);
    }
    else
    {
        return (0);
    }
}

int
fssize_uni(FILE *fp, int *swab)
{
    /*  reads in the size info from a fsfour format
        fortran unformatted map file see fsfor.f */
    /* uni version of fssize - attempts to deal with byte-swabbed maps
     *  that is a map written on big-endian machine  but read on a little endian
     * or vice-versa */
    int nrow, ncol, nplane, n;
    int header, trailer, save_header;
    int nsym, i;
    int buf[200];

    *swab = false;
    rewind(fp);
    /* read first record itle(20),junk,para(6),nsym*/
    header = save_header = fssubs_in_int(fp, *swab);
    if (header != 0)
    {
        n = 20 + 1 + 6 ;
        fread(buf, sizeof(int), n, fp);
        header = header - sizeof(int)*n;
        nsym = fssubs_in_int(fp, *swab);
        header = header -sizeof(int);
    }
    if (nsym < 0 || nsym > (int)MISymmop::MAXSYMMOPS || save_header == 0)
    {
        *swab = true;
        rewind(fp);
        header = save_header = fssubs_in_int(fp, *swab);
        n = 20 + 1 + 6 ;
        fread(buf, sizeof(int), n, fp);
        header = header - sizeof(int)*n;
        nsym = fssubs_in_int(fp, *swab);
        header = header -sizeof(int);
        if (nsym < 0 || nsym > (int)MISymmop::MAXSYMMOPS)
        {
            fprintf(stderr, "fssize_uni: impossible nsym %d\n", nsym);
            return (0);
        }
        else
        {
            fprintf(stderr, "map will swabbed nsym = %d\n", nsym);
        }
    }
    /* go to  next record */
    while (header)
    {
        fread(buf, 1, 1, fp);
        header--;
    }
    trailer = fssubs_in_int(fp, *swab);
    if (trailer != save_header)
    {
        fprintf(stderr, "fssize_uni: error reading FSFOUR map - record sizes don't match\n");
        return 0;
    }
    /* there are nsym records next */
    for (i = 0; i < nsym; i++)
    {
        header = fssubs_in_int(fp, *swab);
        while (header)
        {
            fread(buf, 1, 1, fp);
            header--;
        }
        fread(&trailer, sizeof(int), 1, fp);
    }
    /* read dimension info */
    header = fssubs_in_int(fp, *swab);
    nrow = fssubs_in_int(fp, *swab);
    ncol = fssubs_in_int(fp, *swab);
    nplane = fssubs_in_int(fp, *swab);
    return (ncol*nplane*nrow);
}

int
fswrite(float *rho, FILE *fp, CMapHeaderBase *mh, int norn, int ncent)
{
    /* float * rho; rho[nx*ny*nz] is sorted x fast, y medium, z slow in C rho[nz][ny][nx] */
    /* FILE *fp; the file to output to must be opened in calling routine */
    /* CMapHeaderBase *mh;  see Xguicryst.h for definition */
    /* int norn;  0 = y planes; 1 = x planes; 2 = zplanes */
    /* int ncent;  0 for acentric; 1 for centric (i.e. Patterson) */
    int irow, icol, iplane, ii, ix, iy, iz, i;
    int j;
    int nrow, ncol, nplane;
    int nsym, index;
    int header;
    int sint = sizeof(int), irho;
    //int buf[256];
    int nx, ny, nz;
    float para[6];
    int mnx[6];
    float scale;
    float rspmx = 0.0;
    int nbit = 0;
    int latt = 1; /* fsfour code for p lattice */
    int npic = 10;
    int noset = 0;
    float zero = 0.0;
    char itle[80];

    for (i = 0; i < 80; i++)
    {
        itle[i] = ' ';
    }
    strcpy(itle, "Fsfour format map from xfft (XtalView)");
    scale = (float)mh->scale;
    nsym = mh->nsym;
    nx = mh->nx;
    ny = mh->ny;
    nz = mh->nz;
    para[0] = mh->a;
    para[1] = mh->b;
    para[2] = mh->c;
    para[3] = mh->alpha;
    para[4] = mh->beta;
    para[5] = mh->gamma;
    for (i = 0; i < 6; i++)
    {
        mnx[i] = 0;
    }
    if (norn < 0 || norn > 2)
    {
        norn = 0;
    }
    if (ncent < 0 || ncent > 1)
    {
        ncent = 0;
    }
    if (norn == 0)
    {
        nrow = nx;
        ncol = nz;
        nplane = ny;
    }
    if (norn == 1)
    {
        nrow = ny;
        ncol = nz;
        nplane = nx;
    }
    if (norn == 2)
    {
        nrow = nx;
        ncol = ny;
        nplane = nz;
    }
    header = (20 + 1 + 6 + 4*1)*4;
    fwrite(&header, sint, 1, fp);
    fwrite(itle, 1, 80, fp);
    fwrite(&noset, sint, 1, fp);
    fwrite(para, sint, 6, fp);
    fwrite(&nsym, sint, 1, fp);
    fwrite(&ncent, sint, 1, fp);
    fwrite(&latt, sint, 1, fp);
    fwrite(&npic, sint, 1, fp);
    fwrite(&header, sint, 1, fp);

    /* there are nsym records with 9 numbers next */
    header = 9 * 4;
    for (i = 0; i < nsym; i++)
    {
        fwrite(&header, sint, 1, fp);
        for (j = 0; j < 9; j++)
        {
            fwrite(&zero, sint, 1, fp);
        }
        fwrite(&header, sint, 1, fp);
    }
    header = 4 * ( 7 + 6 );
    fwrite(&header, sint, 1, fp);
    fwrite(&nrow, sizeof(int), 1, fp);
    fwrite(&ncol, sizeof(int), 1, fp);
    fwrite(&nplane, sizeof(int), 1, fp);
    fwrite(&scale, sizeof(int), 1, fp);
    fwrite(&norn, sizeof(int), 1, fp);
    fwrite(mnx, sint, 6, fp);
    fwrite(&rspmx, sint, 1, fp);
    fwrite(&nbit, sint, 1, fp);
    fwrite(&header, sint, 1, fp);

    ii = 0;
    header = 4 * nrow;
    for (iplane = 0; iplane < nplane; iplane++)
    {
        for (icol = 0; icol < ncol; icol++)
        {
            fwrite(&header, sint, 1, fp);
            for (irow = 0; irow < nrow; irow++)
            {
                ii++;
                if (norn == 0)
                {
                    ix = irow;
                    iy = iplane;
                    iz = icol;
                }
                if (norn == 1)
                {
                    ix = iplane;
                    iy = irow;
                    iz = icol;
                }
                if (norn == 2)
                {
                    ix = irow;
                    iy = icol;
                    iz = iplane;
                }
                index = nx*(ny*iz+iy) + ix;
                irho = ROUND(rho[index]);
                fwrite(&irho, sint, 1, fp);
            }
            fwrite(&header, sint, 1, fp);
        }
    }
    return (ii);
}

