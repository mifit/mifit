#include <vector>
#include <set>
#include <algorithm>
#include <cfloat>

#include "core/corelib.h"
#include "uilib.h"
#include <chemlib/Residue.h>

//#include <GL/gl.h>
//#include <GL/glu.h>

#include <jacgrid/jacgrid.h>
#include "surf.h"

#include <math/Matrix4.h>
#include <math/Tuple4.h>

using namespace chemlib;
using namespace mi::math;
using namespace std;

class MISurface;

struct sort_element
{
    unsigned int v[3];
    float z;
};

static int compare(const void *e1, const void *e2)
{
    const sort_element *E1 = (const sort_element*)e1;
    const sort_element *E2 = (const sort_element*)e2;
    if (E1->z < E2->z)
        return 1;
    else if (E1->z == E2->z)
        return 0;
    else
        return -1;
}

static float distsq(const float *a, const float *b)
/*- the square of the distance -*/
{
    return (a[0]-b[0])*(a[0]-b[0])
           +(a[1]-b[1])*(a[1]-b[1])
           +(a[2]-b[2])*(a[2]-b[2]);
}


class MISurface
    : public surfaceT
{
public:
    MISurface()
        : _atom_count(0),
          _hash(0.0f),
          _valid(0)
    {
    }

    MISurface(unsigned int surface_type,
              const std::vector<Molecule*> &mols,
              unsigned int count,
              const unsigned int *selected,
              unsigned int mask = 1)
        : _atom_count(0),
          _hash(0.0f),
          _valid(0)
    {

        // get and verify atom count
        for (unsigned int i = 0; i < mols.size(); ++i)
        {
            Molecule *mol = mols[i];
            for (ResidueListIterator currRes = mol->residuesBegin();
                 currRes != mol->residuesEnd(); ++currRes)
            {
                _atom_count += currRes->atomCount();
            }
        }
        if (count != _atom_count)
        {
            return;
        }

        std::vector<float> xyzr;
        xyzr.resize(_atom_count*4);


        // build xyzr array for JAC surfacing routines
        unsigned int idx = 0;
        for (unsigned int i = 0; i < mols.size(); ++i)
        {
            Molecule *mol = mols[i];
            for (ResidueListIterator currRes = mol->residuesBegin();
                 currRes != mol->residuesEnd(); ++currRes)
            {
                for (int i = 0; i < currRes->atomCount(); ++i)
                {
                    MIAtom *atom = currRes->atom(i);
                    xyzr[idx++] = atom->x();
                    xyzr[idx++] = atom->y();
                    xyzr[idx++] = atom->z();
                    xyzr[idx++] = atom->getRadius();

                    //and calculate hash
                    _hash += atom->x();
                    _hash += atom->y();
                    _hash += atom->z();
                    _hash += atom->getRadius();
                }
            }

        }
        JACAtomsXYZR jacAtoms(&xyzr[0], _atom_count, selected, mask);

        if (surface_type == MISurfaceType::Accessible)
        {
            JACMakeAccessibleSurface(*this, jacAtoms);
        }
        else
        {
            JACMakeMolecularSurface(*this, jacAtoms);
        }

        for (unsigned int ss = 0; ss < MIGetSurfaceSmoothingSteps(); ++ss)
        {
            JACSmoothSurface(*this);
            Reduce(gridUnit[0]*0.5f);
        }
        JACSurfaceNormalize(*this);

        // make owner-to-vertex map; each owner maps to a list of vertices
        size_t max_own = owner.size();
        if (max_own)
        {
            // should be true unless catastrophic failure of surface generation
            unsigned int *own = &owner[0];
            for (size_t i = 0; i < max_own; ++i)
            {
                _o2v[own[i]].insert(i);
            }
        }

        // init colors
        SetColor(Colors::WHITE);
        _valid = true;
    }


    // color surface to a single color
    bool SetColor(short color)
    {
        unsigned int size = vertices.size()/3;
        _color = std::vector<short>(size, color);
        return true;
    }

    // color selected bits of surface with given color
    bool SetColor(short color,
                  const std::vector<Molecule*> &mols,
                  unsigned int count = 0,
                  const unsigned int *selected = 0,
                  unsigned int mask = 1)
    {

        if (!CheckHash(mols, count))
        {
            return false;
        }

        unsigned int aidx = 0;
        for (unsigned int i = 0; i < mols.size(); ++i)
        {
            Molecule *mol = mols[i];
            for (ResidueListIterator currRes = mol->residuesBegin();
                 currRes != mol->residuesEnd(); ++currRes)
            {
                for (int i = 0; i < currRes->atomCount(); ++i, ++aidx)
                {

                    // if (no selected array or atom selected) and atom is
                    // represented in surface, then color all corresponding
                    // vertices
                    if ((!selected || (selected[aidx] & mask))
                        && _o2v.find(aidx) != _o2v.end())
                    {
                        for (std::set<unsigned int>::iterator j = _o2v[aidx].begin();
                             j != _o2v[aidx].end(); ++j)
                        {
                            _color[*j] = color;
                        }
                    }
                }
            }
        }
        return true;
    }


    // color based on atom type;
    bool ColorByAtomType(const std::vector<Molecule*> &mols,
                         unsigned int count = 0,
                         const unsigned int *selected = 0,
                         unsigned int mask = 1)
    {
        if (!CheckHash(mols, count))
        {
            return false;
        }

        unsigned int aidx = 0;
        for (unsigned int i = 0; i < mols.size(); ++i)
        {
            Molecule *mol = mols[i];
            for (ResidueListIterator currRes = mol->residuesBegin();
                 currRes != mol->residuesEnd(); ++currRes)
            {
                for (int i = 0; i < currRes->atomCount(); ++i, ++aidx)
                {

                    // if (no selected array or atom selected) and atom is
                    // represented in surface, then color all corresponding
                    // vertices
                    if ((!selected || (selected[aidx] & mask))
                        && _o2v.find(aidx) != _o2v.end())
                    {
                        for (std::set<unsigned int>::iterator j = _o2v[aidx].begin();
                             j != _o2v[aidx].end(); ++j)
                        {
                            _color[*j] = currRes->atom(i)->color();
                        }
                    }
                }
            }
        }
        return true;
    }

    void Resort()
    {
        float third = 1.0f/3.0f;

        float view[16];
        glGetFloatv(GL_MODELVIEW_MATRIX, view);
        Tuple3<float> rot(view[8], view[9], view[10]); // get z-vector of rotation mat
        if (rot.epsilonEquals(last_rot, 0.01f))
        {
            return;
        }
        last_rot = rot;

        if (!triangles.size() || !vertices.size())
            return;

        Matrix4<float> viewmat(view);
        viewmat.transpose();
        unsigned int ntriangles = triangles.size()/3;

        std::vector<sort_element> elements;
        elements.resize(ntriangles);
        unsigned int *tp = &triangles[0];
        float *verts = &vertices[0];

        for (unsigned int i = 0; i < ntriangles; ++i)
        {
            sort_element &se = elements[i];
            se.v[0] = *tp++;
            se.v[1] = *tp++;
            se.v[2] = *tp++;

            // transform z by viewing matrix
            float v[4] = {0.0f, 0.0f, 0.0f, 1.0f};
            for (int j = 0; j<3; ++j)
            {
                unsigned int vn = se.v[j]*3;
                v[0] += verts[vn+0];
                v[1] += verts[vn+1];
                v[2] += verts[vn+2];
            }
            v[0] *= third;
            v[1] *= third;
            v[2] *= third;

            Tuple4<float> vec(v[0], v[1], v[2], v[3]);
            viewmat.transform(vec);
            //printf("Transforming triange %d, center: %.4f  %.4f  %.4f, z=%.4f\n", i,v[0],v[1],v[2],vec.z);
            se.z = vec.z;
        }
        qsort(&elements[0], ntriangles, sizeof(sort_element), compare);

        //copy back results
        triangles.clear();
        triangles.resize(ntriangles*3);
        tp = &triangles[0];
        for (unsigned int i = 0; i < ntriangles; ++i)
        {
            sort_element &se = elements[i];
            *tp++ = se.v[0];
            *tp++ = se.v[1];
            *tp++ = se.v[2];
        }

    }


    void Draw()
    {
        register unsigned int c1, c2, c3;
        register unsigned int i;

        if (!triangles.size() || !vertices.size())
            return;

        bool do_transparent = (fabs(MIGetSurfaceAlpha() - 1.0f) > 0.001f);
        if (do_transparent)
            Resort();

        glPushAttrib(GL_ALL_ATTRIB_BITS);
        glDisable(GL_CULL_FACE);  // we want to see the inside of clipped surfaces

        // set up color material for surfaces
        float ambient[] = { 0.2f, 0.2f, 0.2f, 1.0f };
        float diffuse[] = { 0.0f, 0.0f, 0.0f, 1.0f };
        float emission[] = { 0.0f, 0.0f, 0.0f, 1.0f };
        float specular[] = { 0.2f, 0.2f, 0.2f, 1.0f };
        int shininess = 128;
        glDisable(GL_COLOR_MATERIAL);
        glMaterialfv(GL_FRONT, GL_AMBIENT, ambient);
        glMaterialfv(GL_FRONT, GL_DIFFUSE, diffuse);
        glMaterialfv(GL_FRONT, GL_EMISSION, emission);
        glMaterialfv(GL_FRONT, GL_SPECULAR, specular);
        glMateriali(GL_FRONT, GL_SHININESS, shininess);

        glColorMaterial(GL_BACK, GL_AMBIENT_AND_DIFFUSE);
        glColor4f(0.75, 0.75, 0.75, 1.0);
        glColorMaterial(GL_FRONT, GL_AMBIENT_AND_DIFFUSE);
        glEnable(GL_COLOR_MATERIAL);

        if (do_transparent)
        {
            glPolygonMode(GL_FRONT_AND_BACK, GL_FILL);
            glDepthMask(GL_FALSE);
            glBlendFunc(GL_SRC_ALPHA, GL_ONE_MINUS_SRC_ALPHA);
            glEnable(GL_BLEND);
            glEnable(GL_STENCIL_TEST);
            glClearStencil(0);
            glClear(GL_STENCIL_BUFFER_BIT);
            glStencilOp(GL_KEEP, GL_KEEP, GL_INCR);
            glStencilFunc(GL_EQUAL, 0, 1);
            glDepthFunc(GL_LEQUAL);
        }

        glBegin(GL_TRIANGLES);
        // avoid vector operator[] fn call overhead in loop
        unsigned int *tri = &triangles[0];
        float *verts = &vertices[0];
        float *norms = &normals[0];
        short *colr = &_color[0];

        if (MIGetSurfaceUseGradColor() && _prop_based_color.size())
        {
            float *colors = &_prop_based_color[0];
            size_t ntri = triangles.size();
            for (i = 0; i < ntri; i += 3)
            {
                c1 = tri[i];
                c2 = tri[i+1];
                c3 = tri[i+2];

                glColor4fv(&colors[c1*4]);
                glNormal3fv(&norms[c1*3]);
                glVertex3fv(&verts[c1*3]);

                glColor4fv(&colors[c2*4]);
                glNormal3fv(&norms[c2*3]);
                glVertex3fv(&verts[c2*3]);

                glColor4fv(&colors[c3*4]);
                glNormal3fv(&norms[c3*3]);
                glVertex3fv(&verts[c3*3]);
            }
        }
        else
        {
            size_t ntri = triangles.size();
            for (i = 0; i < ntri; i += 3)
            {
                c1 = tri[i];
                c2 = tri[i+1];
                c3 = tri[i+2];

                glColor4fv(getColor(colr[c1]));
                glNormal3fv(&norms[c1*3]);
                glVertex3fv(&verts[c1*3]);

                glColor4fv(getColor(colr[c2]));
                glNormal3fv(&norms[c2*3]);
                glVertex3fv(&verts[c2*3]);

                glColor4fv(getColor(colr[c3]));
                glNormal3fv(&norms[c3*3]);
                glVertex3fv(&verts[c3*3]);
            }
        }
        glEnd();
        glPopAttrib();
    }

    bool IsValid() const
    {
        return _valid;
    }



    bool CalcDistance(const std::vector<Molecule*> &mols,
                      unsigned int count,
                      const unsigned int *selected,
                      unsigned int mask = 1)
    {

        //FIXME: wxProgress dialog?

        // for each vertex, find the closest atom
        size_t nv = vertices.size()/3;
        if (!nv)
            return false;
        _color_property = std::vector<float>(nv, FLT_MAX);

        unsigned int aidx = 0;
        for (unsigned int i = 0; i < mols.size(); ++i)
        {
            Molecule *mol = mols[i];
            for (ResidueListIterator currRes = mol->residuesBegin();
                 currRes != mol->residuesEnd(); ++currRes)
            {
                for (int i = 0; i < currRes->atomCount(); ++i, ++aidx)
                {
                    MIAtom *atom = currRes->atom(i);

                    if (!selected || (aidx < count && selected[aidx] & mask))
                    {
                        float *vrtx = &vertices[0];
                        float xyz[3] = { atom->x(), atom->y(), atom->z() };
                        float *prop = &_color_property[0];
                        for (size_t v = 0; v < nv; ++v, vrtx += 3, ++prop)
                        {
                            float d = distsq(xyz, vrtx);
                            if (d < *prop)
                                *prop = d;
                        }
                    }
                }
            }
        }


        float *prop = &_color_property[0];
        for (size_t v = 0; v < nv; ++v, ++prop)
        {
            *prop = (float)sqrt(*prop);
        }
        return true;
    }

protected:

    bool CheckHash(const std::vector<Molecule*> &mols, unsigned int count)
    {
        if (count && count != _atom_count)
        {
            return false;
        }

        float hash = 0.0f;
        for (unsigned int i = 0; i < mols.size(); ++i)
        {
            Molecule *mol = mols[i];
            for (ResidueListIterator currRes = mol->residuesBegin();
                 currRes != mol->residuesEnd(); ++currRes)
            {
                for (int i = 0; i < currRes->atomCount(); ++i)
                {
                    MIAtom *atom = currRes->atom(i);
                    hash += atom->x();
                    hash += atom->y();
                    hash += atom->z();
                    hash += atom->getRadius();
                }
            }
        }
        return hash == _hash;
    }

    const float *getColor(short index)
    {
        static float color[4];
        int c = PaletteIndex(index);
        color[0] = Colors::RPallette[c] / 255.0f;
        color[1] = Colors::GPallette[c] / 255.0f;
        color[2] = Colors::BPallette[c] / 255.0f;
        color[3] = MIGetSurfaceAlpha();
        return color;
    }


    void BuildPropColors()
    {
        size_t prop_siz = _color_property.size();

        if (MIGetSurfaceUseGradColor() && prop_siz)
        {
            float minval = FLT_MAX;
            float maxval = FLT_MIN;
            float *prop = &_color_property[0];
            for (size_t i = 0; i < prop_siz; ++i)
            {
                if (prop[i] > maxval)
                    maxval = prop[i];
                if (prop[i] < minval)
                    minval = prop[i];
            }

            if (fabs(MIGetSurfaceMinValue() - MIGetSurfaceMaxValue()) > 0.01)
            {
                minval = MIGetSurfaceMinValue();
                maxval = MIGetSurfaceMaxValue();
            }

            Logger::log("Surface gradient values range from %f to %f\n", minval, maxval);

            if (maxval == minval)
            {
                maxval += 1.0f; // allow for dividing by maxval-minval below, with out test in inner loop
            }

            _prop_based_color.resize(prop_siz*4);
            float *pc = &_prop_based_color[0];

            const float *c = getColor(MIGetSurfaceMinColor());
            float minc[4] = { c[0], c[1], c[2], c[3] };
            c = getColor(MIGetSurfaceMaxColor());
            float maxc[4] = { c[0], c[1], c[2], c[3] };
            float alpha = MIGetSurfaceAlpha();
            for (size_t i = 0; i < prop_siz; ++i, pc += 4, ++prop)
            {
                float ctmp[4];
                float *color;
                if (*prop < minval)
                {
                    color = minc;
                }
                else if (*prop > maxval)
                {
                    color = maxc;
                }
                else
                {
                    color = ctmp;
                    float f = (*prop - minval)/(maxval-minval);
                    for (int j = 0; j<3; j++)
                    {
                        color[j] = minc[j] + f*(maxc[j] - minc[j]);
                    }
                }

                pc[0] = color[0];
                pc[1] = color[1];
                pc[2] = color[2];
                pc[3] = alpha;
            }
        }
    }





    unsigned int _atom_count;
    float _hash;
    bool _valid;
    Tuple3<float> last_rot;

public:
    std::vector<short> _color;
    std::vector<float> _color_property;
    std::vector<float> _prop_based_color; // 4 * size of _color
    std::map<unsigned int, std::set<unsigned int> >_o2v; // owner to vertex map
};


// since transparency requires all surfaces for a view to be sorted
// together to get the proper triangle ordering, we create this
// "SuperSurface" class to combine the data from all surfaces
class SuperSurface
    : public MISurface
{

public:
    void Rebuild()
    {
        //clear current data, if any
        std::vector<unsigned int>().swap(triangles);
        std::vector<float>().swap(vertices);
        std::vector<float>().swap(normals);
        std::vector<unsigned int>().swap(owner);
        std::vector<short>().swap(_color);
        std::vector<float>().swap(_color_property);

        // reserve room in vectors for efficiency
        size_t tri_size = 0, vert_size = 0, norm_size = 0, owner_size = 0;
        for (size_t i = 0; i < surfaces.size(); ++i)
        {
            MISurface &s = *surfaces[i];
            tri_size += s.triangles.size();
            vert_size += s.vertices.size();
            norm_size += s.normals.size();
            owner_size += s.owner.size();
        }

        triangles.reserve(tri_size);
        vertices.reserve(vert_size);
        normals.reserve(norm_size);
        owner.reserve(owner_size);
        _color.reserve(owner_size);
        _color_property.reserve(owner_size);

        for (size_t i = 0; i < surfaces.size(); ++i )
        {
            MISurface &s = *surfaces[i];
            size_t offset = vertices.size()/3; // store curent number of verts for offseting new vertex indices in triangle list
            size_t jmax = s.triangles.size();
            for (size_t j = 0; j < jmax; ++j)
            {
                triangles.push_back(s.triangles[j] + offset);
            }

            vertices.insert(vertices.end(), s.vertices.begin(), s.vertices.end());
            normals.insert(normals.end(), s.normals.begin(), s.normals.end());
            owner.insert(owner.end(), s.owner.begin(), s.owner.end());
            _color.insert(_color.end(), s._color.begin(), s._color.end());
            _color_property.insert(_color_property.end(), s._color_property.begin(), s._color_property.end());
        }
        BuildPropColors();
    }

    unsigned int Count()
    {
        return surfaces.size();
    }

    void Delete(unsigned int surfnum)
    {
        if (surfnum >= surfaces.size())
            return;
        delete surfaces[surfnum];
        surfaces.erase(surfaces.begin() + surfnum);
        Rebuild();
    }

    unsigned int Append(MISurface *surf)
    {
        surfaces.push_back(surf);
        Rebuild();
        return Count()-1; // index number of last element
    }


    std::vector<MISurface*> surfaces;
};


// each view maintains its own set of surfaces, so...
static void *CURRENT_VIEW = 0;
static std::map<void*, SuperSurface> SuperSurfaces;


//////////////////////////////////////////////////////////////
//
// public api implementation
//
//////////////////////////////////////////////////////////////

unsigned int MIMakeSurface(const std::vector<Molecule*> &mols,
                           unsigned int surface_type,
                           unsigned int atom_count,
                           const unsigned int *selected,
                           unsigned int mask)
{
    MISurface *surf = new MISurface(surface_type, mols, atom_count, selected, mask);
    if (!surf->IsValid())
    {
        delete surf;
        return 0;
    }

    return SuperSurfaces[CURRENT_VIEW].Append(surf);
}

bool MIColorSurface(unsigned int surfnum, short color)
{
    if (surfnum >= SuperSurfaces[CURRENT_VIEW].Count())
        return false;


    bool retval = SuperSurfaces[CURRENT_VIEW].surfaces[surfnum]->SetColor(color);
    SuperSurfaces[CURRENT_VIEW].Rebuild();
    return retval;
}

bool MIColorSurface(unsigned int surfnum,
                    short color,
                    const std::vector<Molecule*> &mols,
                    unsigned int atom_count,
                    const unsigned int *selected,
                    unsigned int mask)
{
    if (surfnum >= SuperSurfaces[CURRENT_VIEW].Count())
        return false;

    bool retval = SuperSurfaces[CURRENT_VIEW].surfaces[surfnum]->SetColor(color, mols, atom_count, selected, mask);
    SuperSurfaces[CURRENT_VIEW].Rebuild();
    return retval;
}


bool MISurfaceCalcDistance(unsigned int surfnum,
                           const std::vector<Molecule*> &mols,
                           unsigned int atom_count,
                           const unsigned int *selected,
                           unsigned int selected_mask)
{
    if (surfnum >= SuperSurfaces[CURRENT_VIEW].Count())
        return false;

    bool retval = SuperSurfaces[CURRENT_VIEW].surfaces[surfnum]->CalcDistance(mols, atom_count, selected, selected_mask);
    SuperSurfaces[CURRENT_VIEW].Rebuild();
    return retval;
}




bool MIAtomColorSurface(unsigned int surfnum,
                        const std::vector<Molecule*> &mols,
                        unsigned int atom_count,
                        const unsigned int *selected,
                        unsigned int mask)
{
    if (surfnum >= SuperSurfaces[CURRENT_VIEW].Count())
        return false;

    bool retval = SuperSurfaces[CURRENT_VIEW].surfaces[surfnum]->ColorByAtomType(mols, atom_count, selected, mask);
    SuperSurfaces[CURRENT_VIEW].Rebuild();
    return retval;
}


void MIDrawSurfaces()
{
    if (SuperSurfaces.find(CURRENT_VIEW) != SuperSurfaces.end())
        SuperSurfaces[CURRENT_VIEW].Draw();
}


void MIDeleteSurface(unsigned int surfnum)
{
    SuperSurfaces[CURRENT_VIEW].Delete(surfnum);
}

unsigned int MISurfaceCount()
{
    return SuperSurfaces[CURRENT_VIEW].Count();
}

void MISurfaceVectorsFromStack(std::vector<Molecule*> &mols,
                               std::vector<unsigned int> &selected,
                               unsigned int selection_type,
                               Stack *AtomStack)
{

    Residue *res;
    Molecule *mol;
    MIAtom *a;
    std::set<MIAtom*> atoms;
    std::set<Residue*> residues;
    std::set<Molecule*> mol_set;

    std::vector<Molecule*>().swap(mols); //was mols.clear();
    std::vector<unsigned int>().swap(selected); // was selected.clear();

    AtomStack->Pop(a, res, mol);
    while (a && res && mol)
    {
        mol_set.insert(mol);
        atoms.insert(a);
        residues.insert(res);

        if (selection_type == MISurfaceSelectionMode::SingleResidue
            || selection_type == MISurfaceSelectionMode::Peptide)
        {
            break;
        }

        AtomStack->Pop(a, res, mol);
    }

    // now, build mols vector from mol_set
    for (std::set<Molecule*>::iterator i = mol_set.begin(); i!=mol_set.end(); ++i)
        mols.push_back(*i);

    // get atom count in all molecules
    unsigned int count = 0;
    for (unsigned int i = 0; i < mols.size(); ++i)
    {
        Molecule *mol = mols[i];
        for (ResidueListIterator currRes = mol->residuesBegin();
             currRes != mol->residuesEnd(); ++currRes)
        {
            count += currRes->atomCount();
        }
    }
    selected.resize(count);

    // build selection vector
    unsigned int atom_index = 0;
    for (unsigned int i = 0; i < mols.size(); ++i)
    {
        Molecule *mol = mols[i];

        for (ResidueListIterator currRes = mol->residuesBegin();
             currRes != mol->residuesEnd(); ++currRes)
        {
            for (int i = 0; i < currRes->atomCount(); ++i)
            {
                MIAtom *atom = currRes->atom(i);
                switch (selection_type)
                {
                case MISurfaceSelectionMode::AtomsOnly:
                    selected[atom_index] = (atoms.find(atom) != atoms.end());
                    break;
                case MISurfaceSelectionMode::SingleResidue:
                case MISurfaceSelectionMode::Residues:
                    selected[atom_index] = (residues.find(currRes) != residues.end());
                    break;
                case MISurfaceSelectionMode::Peptide:
                    selected[atom_index] = IsPeptide(*currRes);
                    break;
                case MISurfaceSelectionMode::EntireMolecule:
                    selected[atom_index] = 1;
                    break;
                }
                atom_index++;
            }
        }
    }
}

static unsigned int CURRENT_SURF_TYPE = MISurfaceType::Molecular;
unsigned int MIGetSurfaceType()
{
    return CURRENT_SURF_TYPE;
}

void MISetSurfaceType(unsigned int s)
{
    CURRENT_SURF_TYPE = s;
}

static unsigned int CURRENT_SURF_SELECTION_MODE = MISurfaceSelectionMode::Peptide;
unsigned int MIGetSurfaceSelectionMode()
{
    return CURRENT_SURF_SELECTION_MODE;
}

void MISetSurfaceSelectionMode(unsigned int s)
{
    CURRENT_SURF_SELECTION_MODE = s;
}

unsigned int SURFACE_QUALITY = 0;
bool MISetSurfaceQuality(unsigned int q)
{
    switch (q)
    {
    default:
        return false;
        break;
    case 0:
        SURFACE_QUALITY = 0;
        JACSetGridDimension(65);
        break;
    case 1:
        SURFACE_QUALITY = 1;
        JACSetGridDimension(97);
        break;
    case 2:
        SURFACE_QUALITY = 2;
        JACSetGridDimension(129);
        break;
    }
    return true;
}

unsigned int MIGetSurfaceQuality()
{
    return SURFACE_QUALITY;
}

unsigned int SMOOTHING_STEPS = 1;
bool MISetSurfaceSmoothingSteps(unsigned int s)
{
    if (SMOOTHING_STEPS > 5)
    {
        return false;
    }
    SMOOTHING_STEPS = s;
    return true;
}

unsigned int MIGetSurfaceSmoothingSteps()
{
    return SMOOTHING_STEPS;
}

void MISurfaceSetCurrentView(void *v)
{
    CURRENT_VIEW = v;
}

float SURFACE_ALPHA = -1.0f;
float MIGetSurfaceAlpha()
{
    if (SURFACE_ALPHA == -1.0f)
        SURFACE_ALPHA = MIConfig::Instance()->GetProfileInt("View Parameters", "surfaceAlpha", 100)/100.0f;
    return SURFACE_ALPHA;
}

void MISetSurfaceAlpha(float alpha)
{
    if (alpha >= 0.0f && alpha <= 1.0f)
    {
        SURFACE_ALPHA = alpha;
        MIConfig::Instance()->WriteProfileInt("View Parameters", "surfaceAlpha", (int)(SURFACE_ALPHA*100.0f));
    }
    if (MIGetSurfaceUseGradColor())
        SuperSurfaces[CURRENT_VIEW].Rebuild();
}


float SURFACE_MINVAL = 0.0f;
float MIGetSurfaceMinValue()
{
    return SURFACE_MINVAL;
}
void MISetSurfaceMinValue(float minval)
{
    SURFACE_MINVAL = minval;
    if (MIGetSurfaceUseGradColor())
        SuperSurfaces[CURRENT_VIEW].Rebuild();
}

float SURFACE_MAXVAL = 0.0f;
float MIGetSurfaceMaxValue()
{
    return SURFACE_MAXVAL;
}
void MISetSurfaceMaxValue(float maxval)
{
    SURFACE_MAXVAL = maxval;
    if (MIGetSurfaceUseGradColor())
        SuperSurfaces[CURRENT_VIEW].Rebuild();
}

int SURFACE_MINCOLOR = -1;
int MIGetSurfaceMinColor()
{
    if (SURFACE_MINCOLOR == -1)
        SURFACE_MINCOLOR = MIConfig::Instance()->GetProfileInt("View Parameters", "surfaceMinColor", Colors::GREEN);
    return SURFACE_MINCOLOR;
}
void MISetSurfaceMinColor(int mincolor)
{
    SURFACE_MINCOLOR = mincolor;
    MIConfig::Instance()->WriteProfileInt("View Parameters", "surfaceMincolor", SURFACE_MINCOLOR);
    if (MIGetSurfaceUseGradColor())
        SuperSurfaces[CURRENT_VIEW].Rebuild();
}


int SURFACE_MAXCOLOR = -1;
int MIGetSurfaceMaxColor()
{
    if (SURFACE_MAXCOLOR == -1)
        SURFACE_MAXCOLOR = MIConfig::Instance()->GetProfileInt("View Parameters", "surfaceMaxColor", Colors::WHITE);
    return SURFACE_MAXCOLOR;
}
void MISetSurfaceMaxColor(int maxcolor)
{
    SURFACE_MAXCOLOR = maxcolor;
    MIConfig::Instance()->WriteProfileInt("View Parameters", "surfaceMaxColor", SURFACE_MAXCOLOR);
    if (MIGetSurfaceUseGradColor())
        SuperSurfaces[CURRENT_VIEW].Rebuild();
}

bool SURFACE_USEGRADCOLOR = false;
bool MIGetSurfaceUseGradColor()
{
    return SURFACE_USEGRADCOLOR;
}
void MISetSurfaceUseGradColor(bool UseGradColor)
{
    SURFACE_USEGRADCOLOR = UseGradColor;
    if (UseGradColor)
        SuperSurfaces[CURRENT_VIEW].Rebuild();
}
