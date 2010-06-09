#include <nongui/nonguilib.h>
#include <math/mathlib.h>
#include <chemlib/chemlib.h>
#include <chemlib/Monomer.h>


#include "maplib.h"
#include "InterpBox.h"
#include "EMapBase.h"
#include "fft.h"



using namespace chemlib;

InterpBox::InterpBox(MIAtomList &vatoms, EMapBase *from_map)
{
    atoms = &(vatoms[0]);
    natoms = vatoms.size();
    emap = from_map;
    Init();
}

//Increased size factor from 1.5 to 2.0 to accomodate refinements while ligand-fitting
void InterpBox::Init()
{
    spacing = 0.2F;
    MIAtom *a;
    a = atoms[0];
    min_x = a->x();
    min_y = a->y();
    min_z = a->z();
    max_x = a->x();
    max_y = a->y();
    max_z = a->z();
    for (int i = 0; i < natoms; i++)
    {
        a = atoms[i];
        min_x = std::min(a->x(), min_x);
        min_y = std::min(a->y(), min_y);
        min_z = std::min(a->z(), min_z);
        max_x = std::max(a->x(), max_x);
        max_y = std::max(a->y(), max_y);
        max_z = std::max(a->z(), max_z);
    }
    float x_size = max_x - min_x;
    float y_size = max_y - min_y;
    float z_size = max_z - min_z;
    float centerx = (min_x+max_x)/2.0F;
    float centery = (min_y+max_y)/2.0F;
    float centerz = (min_z+max_z)/2.0F;
    min_x = centerx - 2.00F*x_size;
    max_x = centerx + 2.00F*x_size;
    min_y = centery - 2.00F*y_size;
    max_y = centery + 2.00F*y_size;
    min_z = centerz - 2.00F*z_size;
    max_z = centerz + 2.00F*z_size;
    nx = ROUND((max_x-min_x)/spacing);
    ny = ROUND((max_y-min_y)/spacing);
    nz = ROUND((max_z-min_z)/spacing);
    grid_points.resize(nx*ny*nz);
    int index;
    float fx, fy, fz;
    float rho;
    CMapHeaderBase *mh = emap->mapheader;
    bool diff_map = false;
    if (emap->mapheader->maptype == (int)MIMapType::FoFc)
    {
        diff_map = true;
    }
    if (emap->mapheader->maptype == (int)MIMapType::FoFcSigmaA)
    {
        diff_map = true;
    }
    max_x = min_x + spacing*nx;
    max_y = min_y + spacing*ny;
    max_z = min_z + spacing*nz;
    for (int ix = 0; ix < nx; ix++)
    {
        for (int iy = 0; iy < ny; iy++)
        {
            for (int iz = 0; iz < nz; iz++)
            {
                fx = min_x + (float)ix *spacing;
                fy = min_y + (float)iy *spacing;
                fz = min_z + (float)iz *spacing;
                index = (nx*(ny*iz+iy)+ix);
                transform(mh->ctof, &fx, &fy, &fz);
                rho = emap->avgrho(fx, fy, fz);
                if (!diff_map)
                {
                    rho *= 2.5F;
                }
                if (rho > 450.0F)
                {
                    rho = 450.0F + 0.25F*(rho-450.0F);
                }
                // penalize heavily breaks
                if (rho < 0)
                {
                    rho *= 2.0F;
                }
                grid_points[index] = rho/50.0F;
            }
        }
    }
}

void InterpBox::ZeroModel(Residue *res)
{
    // zeroes out the density in a radius around each atom
    MIAtom *a;
    float cx, cy, cz;
    int radius = ROUND(0.75F/spacing); // one-half of a C-C bond
    //int radius = ROUND(1.81F/spacing);
    int radius2 = radius*radius;
    int d, dx, dy, dz;
    int ix, iy, iz;
    int iix, iiy, iiz;
    int index;
    float r;
    while (Monomer::isValid(res))
    {
        for (int i = 0; i < res->atomCount(); i++)
        {
            a = res->atom(i);
            if (!(a->type() & AtomType::FITATOM))
            {
                cx = a->x() - min_x;
                cy = a->y() - min_y;
                cz = a->z() - min_z;
                cx /= spacing;
                cy /= spacing;
                cz /= spacing;
                ix = ROUND(cx);
                iy = ROUND(cy);
                iz = ROUND(cz);
                for (iix = ix-radius; iix <= ix+radius; iix++)
                {
                    if (iix >= 0 && iix < nx)
                    {
                        dx = ix-iix;
                        for (iiy = iy-radius; iiy <= iy+radius; iiy++)
                        {
                            if (iiy >= 0 && iiy < ny)
                            {
                                dy = iy-iiy;
                                for (iiz = iz-radius; iiz <= iz+radius; iiz++)
                                {
                                    if (iiz >= 0 && iiz < nz)
                                    {
                                        dz = iz-iiz;
                                        d = dx*dx + dy*dy + dz*dz;
                                        if (d <= radius2)
                                        {
                                            index = (nx*(ny*iiz+iiy)+iix);
                                            r = (-5.0F) + 5.0F*d/radius2;
                                            grid_points[index] = r;
                                        }
                                    }
                                }
                            }
                        }
                    }
                }
            }
        }
        res = res->next();
    }
}

float InterpBox::RDensity(MIAtomList atoms)
{
    float rho = 0;
    int ix, iy, iz, index;
    float cx, cy, cz;
    MIAtom *a;
    for (size_t i = 0; i < atoms.size(); i++)
    {
        a = atoms[i];
        cx = a->x() - min_x;
        if (cx < 0)
        {
            return 0.0;
        }
        cy = a->y() - min_y;
        if (cy < 0)
        {
            return 0.0;
        }
        cz = a->z() - min_z;
        if (cz < 0)
        {
            return 0.0;
        }
        cx /= spacing;
        cy /= spacing;
        cz /= spacing;
        ix = ROUND(cx);
        if (ix >= nx)
        {
            return 0;
        }
        iy = ROUND(cy);
        if (iy >= ny)
        {
            return 0;
        }
        iz = ROUND(cz);
        if (iz >= nz)
        {
            return 0;
        }
        index = (nx*(ny*iz+iy)+ix);
        rho += grid_points[index];
    }
    return rho/(float)atoms.size();
}

