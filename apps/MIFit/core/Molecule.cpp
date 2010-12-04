#include <nongui/nonguilib.h>
#include <math/mathlib.h>
#include <chemlib/chemlib.h>
#include <chemlib/Residue.h>
#include <map/maplib.h>


#include "Molecule.h"
#include "SecondaryStructure.h"
#include "Stack.h"
#include "RESIDUE.h"
#include "ViewPoint.h"
#include "ui/CMapHeader.h"

#include <climits>
#include <vector>
#include <utility>

using namespace chemlib;
using namespace std;

#define X 0
#define Y 1
#define Z 2

Molecule::Molecule(int type)
    : MIMoleculeBase(),
      mapheader(new CMapHeader),
      symmatoms_visible(false),
      m_pSecondaryStructure(NULL),
      secStruc_ribbon(false),
      secStruc_schematic(false)

{
    ModelType = type;
    visible = dots_visible = labels_visible = annots_visible = 1;
    HVisible = 1;
    modelnumber = 0; // 0 means uninitialized valid numbers are above 0
    nribbons = 0;
    ribbon_coloring = 0;
    undoable = false;
    s_main = s_sides = s_waters = s_nonprotein = 1;
    s_link = 0;                                               // the default startup settings
    srfdotsper = 0.5;
    srfboxsize = 3.0;
    symm_center[0] = 0;
    symm_center[1] = 0;
    symm_center[2] = 0;
    symm_radius = 0.0;
}

Molecule::Molecule(Residue *reslist, std::string cmpd, FILE *fp, Bond *conns, int nconns, int type)
    : MIMoleculeBase(reslist, cmpd.c_str(), conns, nconns),
      mapheader(new CMapHeader),
      symmatoms_visible(false),
      m_pSecondaryStructure(NULL),
      secStruc_ribbon(false),
      secStruc_schematic(false)
{
    ModelType = type;
    visible = dots_visible = labels_visible = annots_visible = 1;
    pathname = cmpd;
    compound = cmpd;
    HVisible = 1;
    modelnumber = 0; // 0 means uninitialized valid numbers are above 0
    nribbons = 0;
    ribbon_coloring = 0;
    undoable = false;

    if (fp && ModelType == MoleculeType::PDB)
    {
        rewind(fp);
        GetPDBInfo(fp);
    }

    s_main = s_sides = s_waters = s_nonprotein = 1;
    s_link = 0;                                               // the default startup settings
    srfdotsper = 0.5;
    srfboxsize = 3.0;
    symm_center[0] = 0;
    symm_center[1] = 0;
    symm_center[2] = 0;
    symm_radius = 0.0;
}

void Molecule::SetMapHeader(const CMapHeaderBase &mh)
{
    *mapheader = mh;
}

void Molecule::GetPDBInfo(FILE *fp)
{
    char buf[1000];
    char startname[20], endname[20];
    char startchain, endchain;
    int i, nsecstr = 0;
    int inhead = true;
    compound = "";
    source = "";
    author = "";
    CMapHeaderBase mapHeader;
    while (fgets(buf, sizeof buf, fp) != NULL)
    {
        // truncate PDB code at end of record
        // get rid of trailing whitespace
        buf[80] = '\0';
        for (i = strlen(buf)-1; i > 0; i--)
        {
            if (buf[i] == ' ' || buf[i] == '\t' || buf[i] == '\r' || buf[i] == '\n')
            {
                buf[i] = '\0';
            }
            else
            {
                break;
            }
        }
        if (strncmp(buf, "ATOM ", 5) == 0
            || strncmp(buf, "HETATM", 6) == 0
            || strncmp(buf, "ANISOU", 6) == 0
            || strncmp(buf, "TER ", 4) == 0
            || strncmp(buf, "CONECT", 6) == 0
            || strncmp(buf, "CRYST1", 6) == 0
            || strncmp(buf, "SCALE", 5) == 0)
        {

            inhead = false;
        }
        else
        {
            std::string b = std::string(buf);
            MIStringTrim(b);
            if (inhead)
            {
                FileHead.push_back(b.c_str());
            }
            else
            {
                FileTail.push_back(b.c_str());
            }
        }
        if (!strncmp(buf, "COMPND", 6))
        {
            compound += buf + 10;
        }
        else if (!strncmp(buf, "SOURCE", 6))
        {
            source += buf + 10;
        }
        else if (!strncmp(buf, "AUTHOR", 6))
        {
            author += buf + 10;
        }
        else if (strncmp(buf, "CRYST1", 6) == 0)
        {
            std::string data = buf;
            if (data.length() > 55)
            {
                MIStringToNumber(std::string(data, 6, 9), mapHeader.a);
                MIStringToNumber(std::string(data, 15, 9), mapHeader.b);
                MIStringToNumber(std::string(data, 24, 9), mapHeader.c);
                MIStringToNumber(std::string(data, 33, 7), mapHeader.alpha);
                MIStringToNumber(std::string(data, 40, 7), mapHeader.beta);
                MIStringToNumber(std::string(data, 47, 7), mapHeader.gamma);
                std::string foo = std::string(data, 55, 11);
                MIStringTrim(foo);
                int spaceGroupNumber = mapHeader.FindHMSymbol(foo.c_str());
                if (spaceGroupNumber > 0)
                {
                    mapHeader.spgpno = spaceGroupNumber;
                    mapHeader.SetSymmOps();
                    SetMapHeader(mapHeader);
                }
            }
        }
        else if (strncmp(buf, "SCALE1", 6) == 0)
        {
            std::string data = buf;
            if (data.length() > 30)
            {
                MIStringToNumber(std::string(data, 10, 10), mapHeader.ctof[0][0]);
                MIStringToNumber(std::string(data, 20, 10), mapHeader.ctof[0][1]);
                MIStringToNumber(std::string(data, 30, 10), mapHeader.ctof[0][2]);
                SetMapHeader(mapHeader);
            }
        }
        else if (strncmp(buf, "SCALE2", 6) == 0)
        {
            std::string data = buf;
            if (data.length() > 30)
            {
                MIStringToNumber(std::string(data, 10, 10), mapHeader.ctof[1][0]);
                MIStringToNumber(std::string(data, 20, 10), mapHeader.ctof[1][1]);
                MIStringToNumber(std::string(data, 30, 10), mapHeader.ctof[1][2]);
                SetMapHeader(mapHeader);
            }
        }
        else if (strncmp(buf, "SCALE3", 6) == 0)
        {
            std::string data = buf;
            if (data.length() > 30)
            {
                MIStringToNumber(std::string(data, 10, 10), mapHeader.ctof[2][0]);
                MIStringToNumber(std::string(data, 20, 10), mapHeader.ctof[2][1]);
                MIStringToNumber(std::string(data, 30, 10), mapHeader.ctof[2][2]);
                SetMapHeader(mapHeader);
            }
            //}else if (!strncmp(buf,"REMARK",6)){
            //	if(buf[9]=='2' && strstr(buf,"RESOLUTION")) resolution += buf+11;
            //	if(buf[9]=='3' && (strstr(buf,"R-FACTOR")||strstr(buf,"R-VALUE")||
            //		strstr(buf,"R FACTOR")||strstr(buf,"R VALUE"))) rfactor += buf+11;
        }
        else if (!strncmp(buf, "HELIX", 5))
        {
            nsecstr++;
            sscanf(buf+21, "%s", startname);
            startchain = buf[19];
            sscanf(buf+33, "%s", endname);
            endchain = buf[31];
            ResidueListIterator res1 = residue_from_name(residues, startname, startchain);
            ResidueListIterator res2 = residue_from_name(residues, endname, endchain);

            ResidueListIterator res = residuesBegin();
            for (; res != residuesEnd(); ++res)
            {
                if (res == res1)
                {
                    break;
                }
            }
            if (res != residuesEnd())
            {
                do
                {
                    res->setSecstr('H');
                    ++res;
                } while (res != residuesEnd() && res != res2);
            }
        }
        else if (!strncmp(buf, "SHEET", 5))
        {
            nsecstr++;
            sscanf(buf+22, "%s", startname);
            startchain = buf[21];
            sscanf(buf+33, "%s", endname);
            endchain = buf[32];
            ResidueListIterator res1 = residue_from_name(residues, startname, startchain);
            ResidueListIterator res2 = residue_from_name(residues, endname, endchain);

            ResidueListIterator res = residuesBegin();
            for (; res != residuesEnd(); ++res)
            {
                if (res == res1)
                {
                    break;
                }
            }
            if (res != residuesEnd())
            {
                do
                {
                    res->setSecstr('S');
                    ++res;
                } while (res != residuesEnd() && res != res2);
            }
        }
    }
    if (compound.size() < 4)
    {
        compound = pathname;
    }
    if (nsecstr == 0)
    {
        SecStrFromAngles();
    }
}

Molecule::~Molecule()
{
    moleculeToBeDeleted(this);
    // follow residue list freeing as we go;
    if (residues)
    {
        FreeResidueList(residues);
        residues = 0;
    }
    if (SymmResidues)
    {
        FreeResidueList(SymmResidues);
        SymmResidues = 0;
    }
    if (mapheader != NULL)
    {
        delete mapheader;
    }
    AnnotationList::iterator iter = annotations.begin();
    while (iter != annotations.end())
    {
        Annotation *annotation = *iter;
        ++iter;
        delete annotation;
    }
    std::deque<MIAtom*>::iterator ribbonatomIter;
    for (ribbonatomIter = ribbonatoms.begin(); ribbonatomIter != ribbonatoms.end(); ++ribbonatomIter)
    {
        delete *ribbonatomIter;
    }
}

bool Molecule::VisibleBounds(mi::opengl::Viewpoint *camera, float &axmin, float &axmax,
                             float &aymin, float &aymax, float &azmin, float &azmax)
{
    axmin = std::numeric_limits<float>::max();
    aymin = axmin;
    azmin = axmin;
    axmax = -std::numeric_limits<float>::max();
    aymax = axmax;
    azmax = axmax;

    bool result = false;

    ResidueListIterator res = residuesBegin();
    for (; res != residuesEnd(); ++res)
    {
        for (int i = 0; i < res->atomCount(); i++)
        {
            MIAtom *a = res->atom(i);
            if (a->color() >= 0)
            {
                // changed so that it returns the world coordinates if visible
                // instead of the screen coordinates
                mi::math::Vector3<float> av(a->pos().x(), a->pos().y(), a->pos().z());
                camera->untransformVector(av);
                QVector3D pos(av.x, -av.y, -av.z);
                axmin = std::min(axmin, static_cast<float>(pos.x()));
                aymin = std::min(aymin, static_cast<float>(pos.y()));
                azmin = std::min(azmin, static_cast<float>(pos.z()));
                axmax = std::max(axmax, static_cast<float>(pos.x()));
                aymax = std::max(aymax, static_cast<float>(pos.y()));
                azmax = std::max(azmax, static_cast<float>(pos.z()));
            }
            result = true;
        }
    }
    for (size_t i = 0; i < ribbonatoms.size(); i++)
    {
        MIAtom *a = ribbonatoms[i];
        mi::math::Vector3<float> av(a->pos().x(), a->pos().y(), a->pos().z());
        camera->untransformVector(av);
        QVector3D pos(av.x, -av.y, -av.z);
        axmin = std::min(axmin, static_cast<float>(pos.x()));
        aymin = std::min(aymin, static_cast<float>(pos.y()));
        azmin = std::min(azmin, static_cast<float>(pos.z()));
        axmax = std::max(axmax, static_cast<float>(pos.x()));
        aymax = std::max(aymax, static_cast<float>(pos.y()));
        azmax = std::max(azmax, static_cast<float>(pos.z()));
        result = true;
    }

    return result;
}

void Molecule::FreeDots()
{
    std::vector<SURFDOT>().swap(dots); // was dots.resize(0);
    surfaceChanged(this);
}

void Molecule::Translate(float x, float y, float z, MIAtomList *atoms)
{
    MIMoleculeBase::Translate(x, y, z, atoms);
}

void Molecule::Translate(float x, float y, float z, MIAtomList *atoms, SurfaceDots *dots)
{
    MIAtom *atom;
    SURFDOT *dot;
    unsigned int i;
    for (i = 0; i < atoms->size(); i++)
    {
        atom = (*atoms)[i];
        atom->translate(x, y, z);
    }
    // now move the surface if any
    if (dots)
    {
        for (i = 0; i < dots->size(); i++)
        {
            dot = (SURFDOT*)&dots[i];
            dot->x += x;
            dot->y += y;
            dot->z += z;
        }
    }
    SetCoordsChanged(true);
}

void Molecule::Rotate(float rx, float ry, float rz, float cx, float cy, float cz, ViewPoint *viewpoint,
                      MIAtomList *atoms, SurfaceDots *dots)
{
    float rotmat[3][3], viewmat[3][3];
    float umat[3][3], mat[3][3];
    buildmat(rx, ry, rz, rotmat);
    orthomatrix(rotmat, rotmat);
    viewpoint->copymatrix(viewmat);
    uinv(viewmat, umat);
    matmul(umat, rotmat, mat);
    matmul(mat, viewmat, mat);
    orthomatrix(mat, mat);
    long i;
    float xrot, yrot, zrot, xdir, ydir, zdir;
    SURFDOT *dot;
    MIAtom *a;
    if (atoms)
    {
        for (i = 0; (unsigned int) i < atoms->size(); i++)
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
    // now move the surface if any
    if (dots)
    {
        for (i = 0; (unsigned int) i < dots->size(); i++)
        {
            dot = &(*dots)[i];
            xdir = dot->x - cx;
            ydir = dot->y - cy;
            zdir = dot->z - cz;
            xrot =  (xdir*mat[X][X]+ydir*mat[X][Y]+zdir*mat[X][Z]);
            yrot =  (xdir*mat[Y][X]+ydir*mat[Y][Y]+zdir*mat[Y][Z]);
            zrot =  (xdir*mat[Z][X]+ydir*mat[Z][Y]+zdir*mat[Z][Z]);
            dot->x = xrot + cx;
            dot->y = yrot + cy;
            dot->z = zrot + cz;
        }
    }
    SetCoordsChanged(true);
}

static int color_by_charge(const char *name)
{
    if (!strcmp(name, "ALA"))
    {
        return (Colors::GREEN);
    }
    if (!strcmp(name, "CYS"))
    {
        return (Colors::GREEN);
    }
    if (!strcmp(name, "ASP"))
    {
        return (Colors::RED);
    }
    if (!strcmp(name, "GLU"))
    {
        return (Colors::RED);
    }
    if (!strcmp(name, "PHE"))
    {
        return (Colors::GREEN);
    }
    if (!strcmp(name, "GLY"))
    {
        return (Colors::GREEN);
    }
    if (!strcmp(name, "HIS"))
    {
        return (Colors::CYAN);
    }
    if (!strcmp(name, "ILE"))
    {
        return (Colors::GREEN);
    }
    if (!strcmp(name, "LYS"))
    {
        return (Colors::BLUE);
    }
    if (!strcmp(name, "LEU"))
    {
        return (Colors::GREEN);
    }
    if (!strcmp(name, "MET"))
    {
        return (Colors::GREEN);
    }
    if (!strcmp(name, "ASN"))
    {
        return (Colors::WHITE);
    }
    if (!strcmp(name, "PRO"))
    {
        return (Colors::GREEN);
    }
    if (!strcmp(name, "GLN"))
    {
        return (Colors::WHITE);
    }
    if (!strcmp(name, "ARG"))
    {
        return (Colors::BLUE);
    }
    if (!strcmp(name, "SER"))
    {
        return (Colors::WHITE);
    }
    if (!strcmp(name, "THR"))
    {
        return (Colors::GREEN);
    }
    if (!strcmp(name, "VAL"))
    {
        return (Colors::GREEN);
    }
    if (!strcmp(name, "TRP"))
    {
        return (Colors::WHITE);
    }
    if (!strcmp(name, "TYR"))
    {
        return (Colors::WHITE);
    }
    return Colors::WHITE;
}

static int color_by_phobic(const char *name)
{
    if (!strcmp(name, "ALA"))
    {
        return (Colors::YELLOW);
    }
    if (!strcmp(name, "CYS"))
    {
        return (Colors::YELLOW);
    }
    if (!strcmp(name, "ASP"))
    {
        return (Colors::CYAN);
    }
    if (!strcmp(name, "GLU"))
    {
        return (Colors::CYAN);
    }
    if (!strcmp(name, "PHE"))
    {
        return (Colors::YELLOW);
    }
    if (!strcmp(name, "GLY"))
    {
        return (Colors::BLUE);
    }
    if (!strcmp(name, "HIS"))
    {
        return (Colors::CYAN);
    }
    if (!strcmp(name, "ILE"))
    {
        return (Colors::YELLOW);
    }
    if (!strcmp(name, "LYS"))
    {
        return (Colors::CYAN);
    }
    if (!strcmp(name, "LEU"))
    {
        return (Colors::YELLOW);
    }
    if (!strcmp(name, "MET"))
    {
        return (Colors::YELLOW);
    }
    if (!strcmp(name, "ASN"))
    {
        return (Colors::BLUE);
    }
    if (!strcmp(name, "PRO"))
    {
        return (Colors::YELLOW);
    }
    if (!strcmp(name, "GLN"))
    {
        return (Colors::BLUE);
    }
    if (!strcmp(name, "ARG"))
    {
        return (Colors::CYAN);
    }
    if (!strcmp(name, "SER"))
    {
        return (Colors::BLUE);
    }
    if (!strcmp(name, "THR"))
    {
        return (Colors::GREEN);
    }
    if (!strcmp(name, "VAL"))
    {
        return (Colors::YELLOW);
    }
    if (!strcmp(name, "TRP"))
    {
        return (Colors::YELLOW);
    }
    if (!strcmp(name, "TYR"))
    {
        return (Colors::BLUE);
    }
    return Colors::WHITE;
}

static int color_by_shapely(const char *name)
{
    if (!strcmp(name, "ALA"))
    {
        return (Colors::WHITE);
    }
    if (!strcmp(name, "CYS"))
    {
        return (Colors::YELLOW);
    }
    if (!strcmp(name, "ASP"))
    {
        return (Colors::RED);
    }
    if (!strcmp(name, "GLU"))
    {
        return (Colors::RED);
    }
    if (!strcmp(name, "PHE"))
    {
        return (Colors::MAGENTA);
    }
    if (!strcmp(name, "GLY"))
    {
        return (Colors::WHITE);
    }
    if (!strcmp(name, "HIS"))
    {
        return (Colors::CUSTOM3);
    }
    if (!strcmp(name, "ILE"))
    {
        return (Colors::GREEN);
    }
    if (!strcmp(name, "LYS"))
    {
        return (Colors::BLUE);
    }
    if (!strcmp(name, "LEU"))
    {
        return (Colors::GREEN);
    }
    if (!strcmp(name, "MET"))
    {
        return (Colors::YELLOW);
    }
    if (!strcmp(name, "ASN"))
    {
        return (Colors::CYAN);
    }
    if (!strcmp(name, "PRO"))
    {
        return (Colors::CUSTOM4);
    }
    if (!strcmp(name, "GLN"))
    {
        return (Colors::CYAN);
    }
    if (!strcmp(name, "ARG"))
    {
        return (Colors::BLUE);
    }
    if (!strcmp(name, "SER"))
    {
        return (Colors::CUSTOM5);
    }
    if (!strcmp(name, "THR"))
    {
        return (Colors::CUSTOM5);
    }
    if (!strcmp(name, "VAL"))
    {
        return (Colors::GREEN);
    }
    if (!strcmp(name, "TRP"))
    {
        return (Colors::CUSTOM1);
    }
    if (!strcmp(name, "TYR"))
    {
        return (Colors::MAGENTA);
    }
    return Colors::WHITE;
}

static void setradiustype(MIAtom *a, int type, std::string atoms)
{
    if (type == 0)
    {
        return;
    }
    if (!atoms.empty() && atoms.c_str()[0] != '*')
    {

        std::string aname = a->name();
        aname += ";";
        if (atoms.find(aname) == std::string::npos)
        {
            return;
        }
        else
        {
            a->set_radius_type(type);
        }
    }
    a->set_radius_type(type);
}

int Molecule::getcolor(ResidueListIterator res, MIAtom *a, bool docolor, int color, int colormethod, std::string atoms)
{
    if (!atoms.empty() && atoms[0] != '*')
    {
        std::string aname = a->name();
        aname += ";";
        if (atoms.find(aname) ==std::string::npos)
        {
            return (a->color());
        }
    }
    if (!docolor)
    {
        return (abs(a->color()));
    }
    if (colormethod == Colors::COLORALL)
    {
        return (color);
    }
    else if (colormethod == Colors::COLORC)
    {
        if (a->name()[0] == 'C' && a->atomicnumber() != 17)
        {
            return color;
        }
        else if (a->atomicnumber() == 17)
        {
            return Colors::MAGENTA;
        }
        else
        {
            return (color_by_name(a->name()));
        }
    }
    else if (colormethod == Colors::COLORBVALUE)
    {
        return bvaluecolor(a->BValue());
    }
    else if (colormethod == Colors::COLORSECOND)
    {
        return (secstrcolor(res->secstr()));
    }
    else if (colormethod == Colors::COLORTYPE)
    {
        return (color_by_name(a->name()));
    }
    else if (colormethod == Colors::COLORCHARGE)
    {
        return (color_by_charge(res->type().c_str()));
    }
    else if (colormethod == Colors::COLORHYDRO)
    {
        return (color_by_phobic(res->type().c_str()));
    }
    else if (colormethod == Colors::COLORSHAPELY)
    {
        return (color_by_shapely(res->type().c_str()));
    }
    else if (colormethod == Colors::COLOROFF)
    {
        return (-abs(a->color()));
    }
    else if (colormethod == Colors::COLORON)
    {
        return (abs(a->color()));
    }
    return (a->color());
}

void Molecule::Select(int clear, int main, int sides, int nonprotein,
                      std::string resshow, std::string atoms, Residue *start, Residue *end,
                      bool docolor, int color, int colormethod, int link, int waters, int listtype,
                      int radiustype)
{
    // traverse the residue list selecting residues and atoms
    // as appropiate

    Residue *last = NULL;
    int i, inrange = 0;
    bool dna;
    Bond bond;
    char buf[100];
    s_main = main;
    s_sides = sides;
    s_waters = waters;
    s_nonprotein = nonprotein;
    s_link = link;
    s_radiustype = radiustype;
    for (int j = 0; j < (int)atoms.length(); j++)
    {
        if (atoms[j] == ' ' || atoms[j] == ',')
        {
            atoms[j] = ';';
        }
    }
    atoms += ";";
    for (ResidueListIterator res = residuesBegin(); res != residuesEnd(); ++res)
    {
        if (!start)
        {
            start = res;
        }
        if (clear)
        {
            for (i = 0; i < res->atomCount(); i++)
            {
                res->atom(i)->setColor(-abs(res->atom(i)->color()));
            }
            if (nlinks > 0)
            {
                // since the bonds get sorted the only way to
                // get rid of the old links is to rebuild
                Build();

            }
            nlinks = 0;
        }
        last = res;
    }
    if (end == NULL)
    {
        end = last;
    }
    for (ResidueListIterator res = residuesBegin(); res != residuesEnd(); ++res)
    {
        if (!inrange)
        {
            if (start == res)
            {
                inrange = 1;
            }
            else
            {
                continue;
            }
        }
        if (main && inrange)
        {
            if (IsPeptide(*res))
            {
                for (i = 0; i < res->atomCount(); i++)
                {
                    if (MIAtom::MIIsMainChainAtom(res->atom(i)))
                    {
                        res->atom(i)->setColor(
                            getcolor(res, res->atom(i), docolor, color, colormethod, atoms));
                        setradiustype(res->atom(i), radiustype, atoms);
                    }
                }
            }
            else if (IsDna(*res))
            {
                for (i = 0; i < res->atomCount(); i++)
                {
                    if (MIAtom::MIIsMainChainDNAAtom(res->atom(i)))
                    {
                        res->atom(i)->setColor(
                            getcolor(res, res->atom(i), docolor, color, colormethod, atoms));
                        setradiustype(res->atom(i), radiustype, atoms);
                    }
                }
            }
        }
        if (sides && inrange)
        {
            if (IsPeptide(*res))
            {
                for (i = 0; i < res->atomCount(); i++)
                {
                    if (MIAtom::MIIsSideChainAtom(res->atom(i)))
                    {
                        res->atom(i)->setColor(getcolor(res, res->atom(i), docolor, color, colormethod, atoms));
                        setradiustype(res->atom(i), radiustype, atoms);
                    }
                }
            }
            else if (IsDna(*res))
            {
                for (i = 0; i < res->atomCount(); i++)
                {
                    if (MIAtom::MIIsSideChainDNAAtom(res->atom(i)))
                    {
                        res->atom(i)->setColor(getcolor(res, res->atom(i), docolor, color, colormethod, atoms));
                        setradiustype(res->atom(i), radiustype, atoms);
                    }
                }
            }
        }
        if (nonprotein && inrange)
        {
            if (!IsPeptide(*res) && !IsWater(res) && !IsDna(*res))
            {
                //if(!IsDna(res))
                for (i = 0; i < res->atomCount(); i++)
                {
                    res->atom(i)->setColor(
                        getcolor(res, res->atom(i), docolor, color, colormethod, atoms));
                    setradiustype(res->atom(i), radiustype, atoms);
                }
            }
            //else if(!link)
            //for(i=0;i<res->atomCount();i++)
            //res->atom(i)->color =
            //getcolor(res,res->atom(i),docolor,color,colormethod,atoms);
        }
        if (waters && inrange)
        {
            if (IsWater(res))
            {
                for (i = 0; i < res->atomCount(); i++)
                {
                    res->atom(i)->setColor(
                        getcolor(res, res->atom(i), docolor, color, colormethod, atoms));
                    setradiustype(res->atom(i), radiustype, atoms);
                }
            }
        }
        if (!resshow.empty() && inrange)
        {
            strcpy(buf, ",");
            if (listtype != 1)
            {
                strcat(buf, res->type().c_str());
            }
            else
            {
                strcat(buf, res->name().c_str());
            }
            strcat(buf, ",");
            if (resshow.find(buf) != std::string::npos)
            {
                for (i = 0; i < res->atomCount(); i++)
                {
                    if (MIAtom::MIIsSideChainAtom(res->atom(i)))
                    {
                        res->atom(i)->setColor(
                            getcolor(res, res->atom(i), docolor, color, colormethod, atoms));
                        setradiustype(res->atom(i), radiustype, atoms);
                    }
                }
            }
        }
        if (inrange)
        {
            if (end == res)
            {
                inrange = 0;
            }
        }
    }

    // add links if specified
    if (link)
    {
        MIAtom *a1, *a2, *a3;
        Residue *firstres = residues, *nextres;
        char linkname[10];
        double maxdist;
        Residue *res = residues;
        //if(nlinks > 0 && nlinks <= nbonds){
        // since the bonds get sorted the only way to
        // get rid of the old links is to rebuild
        //Build();
        //}

        //nlinks =0;
        inrange = 0;
        while ((res != NULL) && (nextres = res->next()) != NULL)
        {
            if (!inrange)
            {
                if (start == res)
                {
                    if (res->next() != end)
                    {
                        inrange = 1;
                        firstres = res;
                    }
                }
            }
            if (res->chain_id() != nextres->chain_id())
            {
                nextres = firstres;
                firstres = res->next();
            }
            dna = IsDna(*res) != 0;
            if (dna)
            {
                strcpy(linkname, "P");
                maxdist = 8.0;
            }
            else
            {
                strcpy(linkname, "CA");
                maxdist = 4.5;
            }
            if (inrange && (a1 = atom_from_name(linkname, *res)) != NULL
                && (a2 = atom_from_name(linkname, *nextres)) != NULL)
            {
                if (!linked(a1, a2) && AtomDist(*a1, *a2) < maxdist)
                {
                    // if already linked don't add it again
                    bond.setAtom1(a1);
                    bond.setAtom2(a2);
                    bond.type = B_LINK;
                    bonds.push_back(bond);
                    nlinks++;
                    //if(!checkbonds()) break;
                    if (dna)
                    {
                        if (!strcmp(res->type().c_str(), "A") || !strcmp(res->type().c_str(), "G"))
                        {
                            strcpy(linkname, "N1");
                        }
                        else
                        {
                            strcpy(linkname, "N3");
                        }
                        a3 = atom_from_name(linkname, *res);
                        if (a1 && a3 && !linked(a1, a3))
                        {
                            bond.setAtom1(a1);
                            bond.setAtom2(a3);
                            bond.type = B_LINK;
                            bonds.push_back(bond);
                            a3->setColor(getcolor(res->next(), a3, docolor, color, colormethod, atoms));
                            setradiustype(a3, radiustype, atoms);
                            nlinks++;
                            //if(!checkbonds())break;
                        }
                    }
                }
                a1->setColor(getcolor(res, a1, docolor, color, colormethod, atoms));
                a2->setColor(getcolor(res->next(), a2, docolor, color, colormethod, atoms));
                setradiustype(a1, radiustype, atoms);
                setradiustype(a2, radiustype, atoms);
            }
            if (inrange)
            {
                if (end == res || end == res->next())
                {
                    inrange = 0;
                }
            }
            res = res->next();
        }
    }
}

void Molecule::GenSymmAtoms(ViewPoint *viewpoint)
{
    symm_center[0] = viewpoint->getcenter(0);
    symm_center[1] = viewpoint->getcenter(1);
    symm_center[2] = viewpoint->getcenter(2);
    ClearSymmList();
    symm_radius = (float)viewpoint->width()/viewpoint->scale()*0.7F;
    SymmResidues = SymmResidue(residuesBegin(), mapheader, symm_center, symm_radius);
    if (SymmResidues != NULL)
    {
        Build(true);
    }
}

bool Molecule::CheckCenter(float x, float y, float z)
{
    // if we have symm atoms and we have moved substantially, make a new list

    if (!symmatoms_visible)
        return false;
    float center[3] = {x, y, z};
    double d;
    double dx = center[0] - symm_center[0];
    double dy = center[1] - symm_center[1];
    double dz = center[2] - symm_center[2];
    d = sqrt(dx*dx + dy*dy + dz*dz);
    if (d > 0.60*symm_radius)
    {
        symm_center[0] = x;
        symm_center[1] = y;
        symm_center[2] = z;
        return true;
    }
    return false;
}

MIAtom*Molecule::GetAtom(int natom)
{
    MIAtom *a = MIMoleculeBase::GetAtom(natom);
    if (a)
    {
        return a;
    }


    unsigned int i;
    for (i = 0; i < ribbonatoms.size(); i++)
    {
        if (ribbonatoms[i]->atomnumber() == natom)
        {
            return ribbonatoms[i];
        }
    }
    return NULL;
}

Residue*Molecule::MatchPentamer(std::string &pentdir, Residue *start)
{
    MIAtomList vCA, vCB;
    if (start == NULL)
    {
        return NULL;
    }

    // check that start is in residue list
    ResidueListIterator res = residuesBegin();
    for (; res != residuesEnd(); ++res)
    {
        if (res == ResidueListIterator(start))
        {
            break;
        }
    }
    if (!res)
    {
        return NULL;
    }

    MIAtom *CA, *CB;
    for (; (bool)res && vCA.size() < 5; ++res)
    {
        CA = atom_from_name("CA", *res);
        CB = atom_from_name("CB", *res);
        if (CA)
        {
            vCA.push_back(CA);
            vCB.push_back(CB);
        }
        else
        {
            return NULL;
        }
    }
    return pdbvec(pentdir, vCA, vCB);
}

void Molecule::TranslateResidueToCenter(Residue *res, ViewPoint *viewpoint)
{
    float cx = viewpoint->getcenter(0);
    float cy = viewpoint->getcenter(1);
    float cz = viewpoint->getcenter(2);
    float dx = 0, dy = 0, dz = 0;
    int i;
    for (i = 0; i < res->atomCount(); i++)
    {
        dx += res->atom(i)->x()-cx;
        dy += res->atom(i)->y()-cy;
        dz += res->atom(i)->z()-cz;
    }
    dx /= (float)res->atomCount();
    dy /= (float)res->atomCount();
    dz /= (float)res->atomCount();
    for (i = 0; i < res->atomCount(); i++)
    {
        res->atom(i)->translate(-dx, -dy, -dz);
    }
    SetCoordsChanged(true);
}

void Molecule::TranslateAtomsToCenter(MIAtomList &atoms, ViewPoint *viewpoint)
{
    float cx = viewpoint->getcenter(0);
    float cy = viewpoint->getcenter(1);
    float cz = viewpoint->getcenter(2);
    float dx = 0, dy = 0, dz = 0;
    unsigned int i;
    for (i = 0; i < atoms.size(); i++)
    {
        dx += atoms[i]->x()-cx;
        dy += atoms[i]->y()-cy;
        dz += atoms[i]->z()-cz;
    }
    dx /= (float)atoms.size();
    dy /= (float)atoms.size();
    dz /= (float)atoms.size();
    for (i = 0; i < atoms.size(); i++)
    {
        atoms[i]->translate(-dx, -dy, -dz);
    }
    SetCoordsChanged(true);
}

void Molecule::ShowAll()
{
    MIAtomList changed;
    ResidueListIterator res = residuesBegin();
    for (; res != residuesEnd(); ++res)
    {
        for (int i = 0; i < res->atomCount(); i++)
        {
            if (res->atom(i)->isHidden())
            {
                res->atom(i)->setColor(-res->atom(i)->color());
                changed.push_back(res->atom(i));
            }
        }
    }
    atomChanged(this, changed);
}

void Molecule::HideAll()
{
    // printf("In HideAll\n");

    MIAtomList changed;
    ResidueListIterator res = residuesBegin();
    for (; res != residuesEnd(); ++res)
    {
        for (int i = 0; i < res->atomCount(); i++)
        {
            if (res->atom(i)->color() > 0)
            {
                res->atom(i)->setColor(-res->atom(i)->color());
                changed.push_back(res->atom(i));
            }
        }
    }
    atomChanged(this, changed);
}

void Molecule::ClearRibbons()
{
    // go through bonds deleting  any that are of the the type B_RIBBON
    size_t i = 0;
    while (i < bonds.size())
    {
        if (bonds[i].type == B_RIBBON)
        {
            bonds.erase(bonds.begin()+i);
        }
        else
        {
            i++;
        }
    }
    std::deque<MIAtom*>::iterator iter;
    for (iter = ribbonatoms.begin(); iter != ribbonatoms.end(); ++iter)
    {
        delete *iter;
    }
    std::deque<MIAtom*>().swap(ribbonatoms); // was ribbonatoms.resize(0);
}

MIAtom*Molecule::AllocRibbonAtom()
{
    MIAtom *atom = new MIAtom;
    atom->setColor(Colors::BLUE);
    ribbonatoms.push_back(atom);
    return ribbonatoms.back();
}

void Molecule::BuildRibbons()
{
    Residue *res = residues;
    Residue *nextres = NULL;
    Residue *lastres = NULL;
    MIAtom *O1;
    MIAtom *CA1;
    MIAtom *lastCA1 = NULL;
    MIAtom *lastlastCA1 = NULL;
    MIAtom *nextCA1;
    MIAtom *C1;
    MIAtom *p1;
    MIAtom *p2;
    MIAtom *p3;
    MIAtom *p4;
    MIAtom *pm1;
    MIAtom *pm2;
    MIAtom *pm3;
    MIAtom *pm4;
    MIAtom *t;
    MIAtom *l1;
    MIAtom *l2;
    MIAtom *l3;
    MIAtom *l4;
    float dx, dy, dz;
    float mx, my, mz;
    float dx1 = 0, dy1 = 0, dz1 = 0;
    int i, color;
    float angle;
    int smooth = true;
    Bond bond;
    //bond.color=BLUE;

    ClearRibbons();

    l1 = l2 = l3 = l4 = NULL;

    while (res != NULL)
    {
        if (IsPeptide(*res))
        {
            CA1 = O1 = C1 = nextCA1 = NULL;
            nextres = res->next();
            for (i = 0; i < res->atomCount(); i++)
            {
                if (!strcmp(res->atom(i)->name(), "CA"))
                {
                    CA1 = res->atom(i);
                }
                if (!strcmp(res->atom(i)->name(), "O"))
                {
                    O1 = res->atom(i);
                }
                if (!strcmp(res->atom(i)->name(), "C"))
                {
                    C1 = res->atom(i);
                }
            }
            if (nextres)
            {
                for (i = 0; i < nextres->atomCount(); i++)
                {
                    if (!strcmp(nextres->atom(i)->name(), "CA"))
                    {
                        nextCA1 = nextres->atom(i);
                    }
                }
            }
            if (CA1 && O1 && C1)
            {
                dx = (O1->x() - C1->x())*3/4;
                dy = (O1->y() - C1->y())*3/4;
                dz = (O1->z() - C1->z())*3/4;
                // create four atoms to be used as bonding points
                if (smooth && res->secstr() == 'H' && lastres && lastres->secstr() == 'H'
                    && lastCA1 && lastlastCA1 && nextCA1)
                {
                    // create four atoms to be used as bonding points
                    pm1 = AllocRibbonAtom();
                    pm2 = AllocRibbonAtom();
                    pm3 = AllocRibbonAtom();
                    pm4 = AllocRibbonAtom();
                    if (!pm1 || !pm2 || !pm3 || !pm4)
                    {
                        return;
                    }
                }
                else
                {
                    pm1 = pm2 = pm3 = pm4 = NULL;
                }
                p1 = AllocRibbonAtom();
                p2 = AllocRibbonAtom();
                p3 = AllocRibbonAtom();
                p4 = AllocRibbonAtom();
                if (!p1 || !p2 || !p3 || !p4)
                {
                    return;                  // bail!! error already given
                }
                p1->setColor(secstrcolor(res->secstr()));
                p2->setColor(secstrcolor(res->secstr()));
                p3->setColor(secstrcolor(res->secstr()));
                p4->setColor(secstrcolor(res->secstr()));
                if (ribbon_coloring == 1)
                {
                    color = abs(CA1->color());
                    if (color == 0)
                    {
                        color = color_by_name(CA1->name());
                    }
                    p1->setColor(color);
                    p2->setColor(color);
                    p3->setColor(color);
                    p4->setColor(color);
                }
                p1->setName(res->name().c_str());
                p2->setName(res->name().c_str());
                p3->setName(res->name().c_str());
                p4->setName(res->name().c_str());
                p1->setType(AtomType::RIBBONATOM);
                p2->setType(AtomType::RIBBONATOM);
                p3->setType(AtomType::RIBBONATOM);
                p4->setType(AtomType::RIBBONATOM);
                p1->setOcc(res->secstr());
                p2->setOcc(res->secstr());
                p3->setOcc(res->secstr());
                p4->setOcc(res->secstr());
                if (res->secstr() == 'S')
                {
                    p1->setPosition(CA1->x() + dx, CA1->y() + dy, CA1->z() + dz);
                    p2->setPosition(CA1->x() - dx, CA1->y() - dy, CA1->z() - dz);
                    dx1 = dy1 = dz1 = 0;
                    if (smooth && lastCA1 && nextCA1)
                    {
                        mx = (lastCA1->x() + nextCA1->x())/2;
                        my = (lastCA1->y() + nextCA1->y())/2;
                        mz = (lastCA1->z() + nextCA1->z())/2;
                        dx1 = (mx - CA1->x())/2;
                        dy1 = (my - CA1->y())/2;
                        dz1 = (mz - CA1->z())/2;
                        p1->translate(dx1, dy1, dz1);
                        p2->translate(dx1, dy1, dz1);
                    }
                }
                else if (res->secstr() == 'H')
                {
                    p1->setPosition(CA1->x() + dx, CA1->y() + dy, CA1->z() + dz);
                    p2->setPosition(CA1->x() - dx, CA1->y() - dy, CA1->z() - dz);
                    if (pm1 && pm2 && pm3 && pm4)
                    {
                        pm1->setColor(p1->color());
                        pm2->setColor(p3->color());
                        pm3->setColor(p3->color());
                        pm4->setColor(p4->color());
                        pm1->setName(res->name().c_str());
                        pm2->setName(res->name().c_str());
                        pm3->setName(res->name().c_str());
                        pm4->setName(res->name().c_str());
                        pm1->setType(AtomType::RIBBONATOM);
                        pm2->setType(AtomType::RIBBONATOM);
                        pm3->setType(AtomType::RIBBONATOM);
                        pm4->setType(AtomType::RIBBONATOM);
                        pm1->setOcc(res->secstr());
                        pm2->setOcc(res->secstr());
                        pm3->setOcc(res->secstr());
                        pm4->setOcc(res->secstr());
                        mx = (CA1->x() + lastCA1->x())/2;
                        my = (CA1->y() + lastCA1->y())/2;
                        mz = (CA1->z() + lastCA1->z())/2;
                        dx1 = (CA1->x() - nextCA1->x())+(lastCA1->x() - lastlastCA1->x());
                        dy1 = (CA1->y() - nextCA1->y())+(lastCA1->y() - lastlastCA1->y());
                        dz1 = (CA1->z() - nextCA1->z())+(lastCA1->z() - lastlastCA1->z());
                        dx1 /= 5;
                        dy1 /= 5;
                        dz1 /= 5;
                        mx += dx1;
                        my += dy1;
                        mz += dz1;
                        pm1->setPosition(mx + dx, my + dy, mz + dz);
                        pm2->setPosition(mx - dx, my - dy, mz - dz);
                        pm3->setPosition(pm1->x() + (pm2->x() - pm1->x())/3,
                                         pm1->y() + (pm2->y() - pm1->y())/3,
                                         pm1->z() + (pm2->z() - pm1->z())/3);
                        pm4->setPosition(pm1->x() + (pm2->x() - pm1->x())*2/3,
                                         pm1->y() + (pm2->y() - pm1->y())*2/3,
                                         pm1->z() + (pm2->z() - pm1->z())*2/3);
                    }
                }
                else
                {
                    p1->setPosition(CA1->x() + dx*3/10,
                                    CA1->y() + dy*3/10,
                                    CA1->z() + dz*3/10);
                    p2->setPosition(CA1->x() - dx*3/10,
                                    CA1->y() - dy*3/10,
                                    CA1->z() - dz*3/10);
                    if (smooth && lastCA1 && nextCA1 && lastres && nextres
                        && nextres->secstr() != 'H' && lastres->secstr() != 'H'
                        && nextres->secstr() != 'S' && lastres->secstr() != 'S')
                    {
                        mx = (lastCA1->x() + nextCA1->x())/2;
                        my = (lastCA1->y() + nextCA1->y())/2;
                        mz = (lastCA1->z() + nextCA1->z())/2;
                        dx1 = (mx - CA1->x())/2;
                        dy1 = (my - CA1->y())/2;
                        dz1 = (mz - CA1->z())/2;
                        p1->translate(dx1, dy1, dz1);
                        p2->translate(dx1, dy1, dz1);
                    }
                }
                p3->setPosition(p1->x() + (p2->x() - p1->x())/3,
                                p1->y() + (p2->y() - p1->y())/3,
                                p1->z() + (p2->z() - p1->z())/3);
                p4->setPosition(p1->x() + (p2->x() - p1->x())*2/3,
                                p1->y() + (p2->y() - p1->y())*2/3,
                                p1->z() + (p2->z() - p1->z())*2/3);
                if (l1 && l2 && l3 && l4 && AtomDist(*p3, *l3) < 4.50)
                {
                    // find the order that causes least twisting
                    angle = CalcAtomTorsion(l1, l2, p2, p1);
                    if (fabs(angle) > 90.0)
                    {
                        t = p1;
                        p1 = p2;
                        p2 = t;
                        t = p3;
                        p3 = p4;
                        p4 = t;
                        if (pm1 && pm2 && pm3 && pm4)
                        {
                            t = pm1;
                            pm1 = pm2;
                            pm2 = t;
                            t = pm3;
                            pm3 = pm4;
                            pm4 = t;
                        }
                    }
                    p1->setBValue(1);
                    p2->setBValue(2);
                    p3->setBValue(3);
                    p4->setBValue(4);
                    if (pm1 && pm2 && pm3 && pm4)
                    {
                        pm1->setBValue(1);
                        pm2->setBValue(2);
                        pm3->setBValue(3);
                        pm4->setBValue(4);
                        nribbons++;
                        bond.type = B_RIBBON;
                        bond.setAtom1(pm1);
                        bond.setAtom2(l1);
                        bonds.push_back(bond);
                        bond.type = B_RIBBON;
                        bond.setAtom1(pm2);
                        bond.setAtom2(l2);
                        bonds.push_back(bond);
                        bond.type = B_RIBBON;
                        bond.setAtom1(pm3);
                        bond.setAtom2(l3);
                        bonds.push_back(bond);
                        bond.type = B_RIBBON;
                        bond.setAtom1(pm4);
                        bond.setAtom2(l4);
                        bonds.push_back(bond);
                        nribbons++;
                        bond.type = B_RIBBON;
                        bond.setAtom1(pm1);
                        bond.setAtom2(p1);
                        bonds.push_back(bond);
                        bond.type = B_RIBBON;
                        bond.setAtom1(pm2);
                        bond.setAtom2(p2);
                        bonds.push_back(bond);
                        bond.type = B_RIBBON;
                        bond.setAtom1(pm3);
                        bond.setAtom2(p3);
                        bonds.push_back(bond);
                        bond.type = B_RIBBON;
                        bond.setAtom1(pm4);
                        bond.setAtom2(p4);
                        bonds.push_back(bond);
                    }
                    else
                    {
                        nribbons++;
                        bond.type = B_RIBBON;
                        bond.setAtom1(p1);
                        bond.setAtom2(l1);
                        bonds.push_back(bond);
                        bond.type = B_RIBBON;
                        bond.setAtom1(p2);
                        bond.setAtom2(l2);
                        bonds.push_back(bond);
                        bond.type = B_RIBBON;
                        bond.setAtom1(p3);
                        bond.setAtom2(l3);
                        bonds.push_back(bond);
                        bond.type = B_RIBBON;
                        bond.setAtom1(p4);
                        bond.setAtom2(l4);
                        bonds.push_back(bond);
                    }
                }
                l1 = p1;
                l2 = p2;
                l3 = p3;
                l4 = p4;
                lastCA1 = CA1;
                if (res->next())
                {
                    if ((p1->occ() == 'S' || p1->occ() == 'H')
                        && (res->next()->secstr() !=  p1->occ()))
                    {
                        // add an arrowhead
                        p1 = AllocRibbonAtom();
                        p2 = AllocRibbonAtom();
                        p3 = AllocRibbonAtom();
                        p4 = AllocRibbonAtom();
                        if (!p1 || !p2 || !p3 || !p4)
                        {
                            return;            // bail!! error already given
                        }
                        if (ribbon_coloring == 1)
                        {
                            color = abs(CA1->color());
                            if (color == 0)
                            {
                                color = color_by_name(CA1->name());
                            }
                        }
                        else
                        {
                            color = secstrcolor(res->secstr());
                        }
                        p1->setColor(color);
                        p2->setColor(color);
                        p3->setColor(color);
                        p4->setColor(color);
                        p1->setName(res->name().c_str());
                        p2->setName(res->name().c_str());
                        p3->setName(res->name().c_str());
                        p4->setName(res->name().c_str());
                        p1->setType(AtomType::RIBBONATOM);
                        p2->setType(AtomType::RIBBONATOM);
                        p3->setType(AtomType::RIBBONATOM);
                        p4->setType(AtomType::RIBBONATOM);
                        p1->setOcc(res->secstr());
                        p2->setOcc(res->secstr());
                        p3->setOcc(res->secstr());
                        p4->setOcc(res->secstr());
                        p1->setPosition(CA1->x() + dx*2,
                                        CA1->y() + dy*2,
                                        CA1->z() + dz*2);
                        p2->setPosition(CA1->x() - dx*2,
                                        CA1->y() - dy*2,
                                        CA1->z() - dz*2);
                        if (res->secstr() == 'S')
                        {
                            p1->translate(dx1, dy1, dz1);
                            p2->translate(dx1, dy1, dz1);
                        }
                        p3->setPosition(p1->x() + (p2->x() - p1->x())/3,
                                        p1->y() + (p2->y() - p1->y())/3,
                                        p1->z() + (p2->z() - p1->z())/3);
                        p4->setPosition(p1->x() + (p2->x() - p1->x())*2/3,
                                        p1->y() + (p2->y() - p1->y())*2/3,
                                        p1->z() + (p2->z() - p1->z())*2/3);
                        p1->setBValue(1);
                        p2->setBValue(2);
                        p3->setBValue(3);
                        p4->setBValue(4);
                        nribbons++;
                        bond.type = B_RIBBON;
                        bond.setAtom1(p1);
                        bond.setAtom2(l1);
                        bonds.push_back(bond);
                        bond.type = B_RIBBON;
                        bond.setAtom1(p2);
                        bond.setAtom2(l2);
                        bonds.push_back(bond);
                        bond.type = B_RIBBON;
                        bond.setAtom1(p3);
                        bond.setAtom2(l3);
                        bonds.push_back(bond);
                        bond.type = B_RIBBON;
                        bond.setAtom1(p4);
                        bond.setAtom2(l4);
                        bonds.push_back(bond);
                        l1 = p1;
                        l2 = p2;
                        l3 = p3;
                        l4 = p4;
                    }
                }
            }
            else
            {
                l1 = l2 = l3 = l4 = NULL;
                lastCA1 = CA1;
            }
        }
        else
        {
            l1 = l2 = l3 = l4 = NULL;
            lastCA1 = NULL;
        }
        lastres = res;
        lastlastCA1 = lastCA1;
        res = res->next();
    }
}

void Molecule::Save(XMLArchive &ar)
{
    // give all the atoms a unique number
    size_t i;
    int n = 0;
    for (ResidueListIterator res = residuesBegin(); res != residuesEnd(); ++res)
    {
        for (int i = 0; i < res->atomCount(); i++)
        {
            res->atom(i)->setAtomnumber(n++);
        }
    }
    for (ResidueListIterator res = symmResiduesBegin(); res != symmResiduesEnd(); ++res)
    {
        for (int i = 0; i < res->atomCount(); i++)
        {
            res->atom(i)->setAtomnumber(n++);
        }
    }
    for (i = 0; i < ribbonatoms.size(); i++)
    {
        ribbonatoms[i]->setAtomnumber(n++);
    }
    ar.BeginTag("Molecule");
    ar.WriteField("Compound", (const char*)compound.c_str());
    ar.WriteField("Pathname", (const char*)pathname.c_str());
    ar.WriteField("Author", (const char*)author.c_str());
    ar.WriteField("Source", (const char*)source.c_str());
    ar.WriteField("dots_visible", dots_visible);
    ar.WriteField("visible", visible);
    ar.WriteField("HVisible", HVisible);
    ar.WriteField("labels_visible", labels_visible);
    ar.WriteField("link_here", link_here);
    ar.WriteField("link_next", link_next);
    ar.WriteField("modelnumber", modelnumber);
    ar.WriteField("nresidues", nresidues);
    ar.WriteField("nlinks", nlinks);
    ar.WriteField("nribbons", nribbons);
    ar.WriteField("ribbon_coloring", ribbon_coloring);
    ar.WriteField("s_link", s_link);
    ar.WriteField("s_main", s_main);
    ar.WriteField("s_nonprotein", s_nonprotein);
    ar.WriteField("s_radiustype", s_radiustype);
    ar.WriteField("s_sides", s_sides);
    ar.WriteField("s_waters", s_waters);
    ar.WriteField("srfboxsize", srfboxsize);
    ar.WriteField("srfdotsper", srfdotsper);
    ar.WriteField("symm_radius", symm_radius);
    ar.WriteField("symmatoms_visible", symmatoms_visible);
    ar.WriteSymmCenter(symm_center);
    ar.BeginTag("ResidueList");
    for (ResidueListIterator res = residuesBegin(); res != residuesEnd(); ++res)
    {
        ar.WriteField(*res);
    }
    ar.EndTag();
    // write connects
    for (i = 0; i < connects.size(); i++)
    {
        ar.WriteConnect(connects[i]);
    }
    // write mapheader
    ar.WriteMapHeader(*mapheader);
    // write SymmResidues
    ar.BeginTag("SymmResidueList");
    for (ResidueListIterator res = symmResiduesBegin(); res != symmResiduesEnd(); ++res)
    {
        ar.WriteField(*res);
    }
    ar.EndTag();
    // write ribbonatoms
    if (ribbonatoms.size() > 0)
    {
        ar.BeginTag("RibbonAtoms");
        for (i = 0; i < ribbonatoms.size(); i++)
        {
            ar.WriteField(ribbonatoms[i]);
        }
        ar.EndTag();
    }
    if (hbonds.size() > 0)
    {
        for (i = 0; i < hbonds.size(); i++)
        {
            ar.WriteHBond(hbonds[i]);
        }
    }
    // write labels
    if (atomLabels.size() > 0)
    {
        AtomLabelList::iterator iter = atomLabels.begin();
        while (iter != atomLabels.end())
        {
            ATOMLABEL *label = *iter;
            ++iter;
            ar.WriteLabel(*label);
        }
    }
    // write Annotation
    if (annotations.size() > 0)
    {
        AnnotationList::iterator iter = annotations.begin();
        while (iter != annotations.end())
        {
            Annotation *annotation = *iter;
            ++iter;
            ar.WriteAnnotation(annotation);
        }
    }
    // write bonds
    for (i = 0; i < bonds.size(); i++)
    {
        ar.WriteBond(bonds[i]);
    }
    // write dots
    for (i = 0; i < dots.size(); i++)
    {
        ar.WriteDot(dots[i]);
    }
    ar.EndTag();
    SetCoordsChanged(false);
}

void Molecule::Do()
{
    short savec;
    for (ResidueListIterator res = residuesBegin(); res != residuesEnd(); ++res)
    {
        for (int i = 0; i < res->atomCount(); i++)
        {
            savec = (short)res->atom(i)->color()*16;
            savec |= (short)res->atom(i)->radius_type();
            savecolors.push_back(savec);
        }
    }
    undoable = true;
}

void Molecule::UnDo()
{
    unsigned int j = 0;
    if (savecolors.size() <= 0)
    {
        return;
    }
    for (ResidueListIterator res = residuesBegin(); res != residuesEnd(); ++res)
    {
        for (int i = 0; i < res->atomCount(); i++)
        {
            if (j < savecolors.size())
            {
                res->atom(i)->setColor(savecolors[j]/16);
                res->atom(i)->set_radius_type(savecolors[j]& (short)(!16));
                j++;
            }
        }
    }
    undoable = false;
}

long Molecule::SolventSurface(ViewPoint*, float dotsper)
{
    int j, k;
    MIAtom *b[100];
    MIAtom *atom;
    extern int SurfResult;
    SurfResult = 1;
    srfdotsper = dotsper;
    extern float spacing;
    float r1;
    SURFDOT *adots;
    long nadots = 0;
    long maxadots = 0;
    void *hdots = NULL;

    std::vector<SURFDOT>().swap(dots); // was dots.clear();
    spacing = dotsper;

    float radius_offset = 1.4f;
    float radius[NTYPES];
    for (unsigned int i = 0; i < NTYPES; ++i)
    {
        radius[i] = MIAtom::MIAtomRadiusForType(i)+radius_offset;
    }
    initspheres(radius);

    for (ResidueListIterator res = residuesBegin(); res != residuesEnd() && SurfResult != 0; ++res)
    {
        if (!IsWater(res))
        {
            for (k = 0; k < res->atomCount(); k++)
            {
                int nb = 0;
                atom = res->atom(k);
                r1 = atom->getRadius() + radius_offset;
                nb = 0;
                for (ResidueListIterator res2 = residuesBegin(); res2 != residuesEnd(); ++res2)
                {
                    if (!IsWater(res2))
                    {
                        for (j = 0; j < res2->atomCount(); j++)
                        {
                            if (res2->atom(j) != atom)
                            {
                                if (AtomDist(*atom, *res2->atom(j)) < (r1 + res2->atom(j)->getRadius()+radius_offset))
                                {
                                    if (nb < 100)
                                    {
                                        b[nb] = res2->atom(j);
                                        nb++;
                                    }
                                }
                            }
                        }
                    }
                }
                nadots = 0;
                nadots = atomsurf(atom, 1.0F, b, nb, &adots, nadots, &maxadots, hdots, srfdotsper, radius_offset);
                if (SurfResult == 0)
                {
                    break;
                }
                for (int l = 0; l < nadots; l++)
                {
                    dots.push_back(adots[l]);
                }
            }
        }
    }
    if (hdots != NULL)
    {
        free(hdots);
        hdots = NULL;
    }
    surfaceChanged(this);
    return (dots.size());
}

long Molecule::SurfaceResidues(float dotsper, bool ignore_hidden)
{
    extern int SurfResult;
    SurfResult = 1;
    srfdotsper = dotsper;
    for (ResidueListIterator res = residuesBegin(); res != residuesEnd() && SurfResult != 0; ++res)
    {
        if (res->flags()&128)
        {
            for (int i = 0; i < res->atomCount(); ++i)
            {
                Surface(res->atom(i), ignore_hidden, false);
                if (SurfResult == 0)
                    break;
            }
            res->setFlags(res->flags()&(!128));
        }
    }
    surfaceChanged(this);
    return dots.size();
}

long Molecule::SurfaceResidue(Residue *res, float dotsper, bool ignore_hidden)
{
    int i;
    long nadots = 0;
    long maxadots = 0;
    void *hdots = NULL;
    SURFDOT *adots;

    float r1;
    srfdotsper = dotsper;
    extern int SurfResult;
    SurfResult = 1;
    if (res != NULL && SurfResult != 0)
    {
        MIAtom *b[100];
        MIAtom *atom;
        int nb = 0;
        for (int j = 0; j < res->atomCount(); j++)
        {
            atom = res->atom(j);
            r1 = atom->getRadius();
            nb = 0;
            for (i = 0; i < res->atomCount(); i++)
            {
                if (res->atom(i) != atom && (res->atom(i)->color() >= 0 || (!ignore_hidden)))
                {
                    if (AtomDist(*atom, *res->atom(i)) < (r1 + res->atom(i)->getRadius()))
                    {
                        if (nb < 100)
                        {
                            b[nb] = res->atom(i);
                            nb++;
                        }
                    }
                }
            }
            nadots = 0;
            nadots = atomsurf(atom, 1.0F, b, nb, &adots, nadots, &maxadots, hdots, srfdotsper);
            if (SurfResult == 0)
            {
                break;
            }
            for (int l = 0; l < nadots; l++)
            {
                dots.push_back(adots[l]);
            }
        }
    }
    if (hdots != NULL)
    {
        free(hdots);
        hdots = NULL;
    }
    surfaceChanged(this);
    return (dots.size());
}

long Molecule::SurfaceAtom(MIAtom *atom, float dotsper, bool ignore_hidden)
{
    srfdotsper = dotsper;
    long result = Surface(atom, ignore_hidden);
    return (result);
}

long Molecule::SurfaceAroundAtom(MIAtom *atom, float dotsper, float radius)
{
    SURFDOT *adots;
    long nadots = 0;
    long maxadots = 0;
    void *hdots = NULL;

    srfdotsper = dotsper;
    nadots = atomsurfradius(atom, radius, &adots, nadots, &maxadots, hdots, dotsper);
    for (int l = 0; l < nadots; l++)
    {
        dots.push_back(adots[l]);
    }
    if (hdots != NULL)
    {
        free(hdots);
        hdots = NULL;
    }
    surfaceChanged(this);
    return (dots.size());
}

long Molecule::Surface(MIAtom *atom, bool ignore_hidden, bool send_signal)
{
    int i;
    MIAtom *b[100];
    int nb = 0;
    float r1 = atom->getRadius();
    SURFDOT *adots;
    long nadots = 0;
    long maxadots = 0;
    void *hdots = NULL;

    for (ResidueListIterator res = residuesBegin(); res != residuesEnd(); ++res)
    {
        for (i = 0; i < res->atomCount(); i++)
        {
            if (res->atom(i) != atom && (res->atom(i)->color() >= 0 || (!ignore_hidden)))
            {
                if (AtomDist(*atom, *res->atom(i)) < (r1 + res->atom(i)->getRadius()))
                {
                    if (nb < 100)
                    {
                        b[nb] = res->atom(i);
                        nb++;
                    }
                }
            }
        }
    }
    nadots = 0;
    nadots = atomsurf(atom, 1.0F, b, nb, &adots, nadots, &maxadots, hdots, srfdotsper);
    for (int l = 0; l < nadots; l++)
    {
        dots.push_back(adots[l]);
    }
    if (hdots != NULL)
    {
        free(hdots);
        hdots = NULL;
    }
    if (send_signal)
    {
        surfaceChanged(this);
    }
    return dots.size();
}

void Molecule::ShowHydrogens(bool Show)
{
    for (ResidueListIterator res = residuesBegin(); res != residuesEnd(); ++res)
    {
        for (int i = 0; i < res->atomCount(); i++)
        {
            if (MIAtom::MIIsHydrogen(res->atom(i)) )
            {
                if (Show)
                {
                    res->atom(i)->setColor(abs(res->atom(i)->color()));
                }
                else
                {
                    res->atom(i)->setColor(-abs(res->atom(i)->color()));
                }
            }
        }
    }
}

long Molecule::SurfaceCenter(ViewPoint *viewpoint, float dotsper, float boxsize, bool ignore_hidden)
{
    float cx = viewpoint->getcenter(0);
    float cy = viewpoint->getcenter(1);
    float cz = viewpoint->getcenter(2);
    extern int SurfResult;
    SurfResult = 1;
    srfdotsper = dotsper;
    srfboxsize = boxsize;
    boxsize *= boxsize;
    std::vector<SURFDOT>().swap(dots);
    for (ResidueListIterator res = residuesBegin(); res != residuesEnd() && SurfResult != 0; ++res)
    {
        for (int i = 0; i < res->atomCount(); ++i)
        {
            float x = res->atom(i)->x() - cx;
            float y = res->atom(i)->y() - cy;
            float z = res->atom(i)->z() - cz;
            if (x*x + y*y + z*z < boxsize)
                Surface(res->atom(i), ignore_hidden, false);

            if (SurfResult == 0)
                break;
        }
    }
    surfaceChanged(this);
    return (dots.size());
}

SecondaryStructure*Molecule::getSecondaryStructure()
{
    return m_pSecondaryStructure;
}

void Molecule::MakeSecondaryStructure(bool bRibbon, bool bSchematic)
{

    secStruc_ribbon = bRibbon;
    secStruc_schematic = bSchematic;
    if (m_pSecondaryStructure)
    {
        delete m_pSecondaryStructure;
        m_pSecondaryStructure = NULL;
    }
    if (!bRibbon && !bSchematic)
    {
        return;
    }
    m_pSecondaryStructure = new SecondaryStructure();

    try
    {
        if (bRibbon)
        {
            Residue *res = residues;
            while (res != NULL)
            {
                while ((res != NULL) && atom_from_name("CA", *res) == NULL)
                {
                    res = res->next();
                }
                if (res == NULL)
                {
                    break;
                }
                Residue *startResidue = res;
                int residueCount = 0;
                unsigned short currentChain = res->chain_id();
                while ((res != NULL) && atom_from_name("CA", *res) != NULL && res->chain_id() == currentChain)
                {
                    residueCount++;
                    res = res->next();
                }
                m_pSecondaryStructure->AddRibbonSegment(startResidue, residueCount-1);
            }
        }

        if (bSchematic)
        {
            std::vector<std::pair<Residue*, Residue*> > pHelix;
            std::vector<std::pair<Residue*, Residue*> > pSheet;
            std::vector<std::pair<Residue*, Residue*> > pTurn;
            std::vector<std::pair<Residue*, Residue*> > pRandom;

            Residue *res = residues;
            int residueCount = 0;
            while (res != NULL)
            {
                while ((res != NULL) && atom_from_name("CA", *res) == NULL)
                {
                    res = res->next();
                }
                if (res == NULL)
                {
                    break;
                }
                residueCount = 0;
                unsigned short currentChain = res->chain_id();
                int currentSS;
                char startSS = res->secstr();
                Residue *startRes = res;
                Residue *endRes = NULL;
                Residue *prevRes = res;
                while (res != NULL)
                {
                    Residue *nextRes = res->next();
                    if ((nextRes != NULL) && ((atom_from_name("CA", *nextRes) == NULL) || nextRes->chain_id() != currentChain))
                    {
                        res = res->next();
                        residueCount++;
                        break;
                    }
                    currentSS = res->secstr();
                    if (currentSS != 'H' && currentSS != 'S' && currentSS != 'T')
                    {
                        currentSS = 'U';
                    }
                    if (startSS != currentSS)
                    {
                        endRes = prevRes;
                        Logger::debug("end %s %s %c %d", endRes->type().c_str(), endRes->name().c_str(), endRes->chain_id(), residueCount);
                        switch (startSS)
                        {
                        default:
                        case 'U':
                            pRandom.push_back(std::make_pair(startRes, endRes));
                            break;
                        case 'H':
                            pHelix.push_back(std::make_pair(startRes, endRes));
                            break;
                        case 'S':
                            pSheet.push_back(std::make_pair(startRes, endRes));
                            break;
                        case 'T':
                            pTurn.push_back(std::make_pair(startRes, endRes));
                            break;
                        }
                        startSS = currentSS;
                        startRes = prevRes;
                        Logger::debug("start %s %s %c", startRes->type().c_str(), startRes->name().c_str(), startRes->chain_id());
                        residueCount = 0;
                    }
                    prevRes = res;
                    res = res->next();
                    residueCount++;
                }
                endRes = prevRes;
                Logger::debug("end %s %s %c %d", endRes->type().c_str(), endRes->name().c_str(), endRes->chain_id(), residueCount);
                switch (startSS)
                {
                default:
                case 'U':
                    pRandom.push_back(std::make_pair(startRes, endRes));
                    break;
                case 'H':
                    pHelix.push_back(std::make_pair(startRes, endRes));
                    break;
                case 'S':
                    pSheet.push_back(std::make_pair(startRes, endRes));
                    break;
                case 'T':
                    pTurn.push_back(std::make_pair(startRes, endRes));
                    break;
                }
            }
            m_pSecondaryStructure->AddSchematic(this, pHelix, pSheet, pTurn, pRandom);
        }
    }
    catch (std::string e)
    {
        std::string s;
        s = format("Unable to create secondary structure:\n%s", e.c_str());
        Logger::message(s);
    }

}

void Molecule::DeleteSecondaryStructure()
{
    if (m_pSecondaryStructure)
    {
        delete m_pSecondaryStructure;
        m_pSecondaryStructure = NULL;
        secStruc_ribbon = false;
        secStruc_schematic = false;
    }
}

Molecule::AnnotationList&Molecule::getAnnotations()
{
    return annotations;
}

void Molecule::addAnnotation(Annotation *annotation)
{
    if (annotation->m_id == 0)
    {
        annotation->m_id = annotations.size()+1;
    }
    annotations.push_back(annotation);
    annotationAdded(this, annotation);
}

void Molecule::deleteAnnotation(Annotation *annotation)
{
    AnnotationList::iterator pos = find(annotations.begin(), annotations.end(), annotation);
    if (pos != annotations.end())
    {
        Annotation *annotation = *pos;
        annotations.erase(pos);
        doAnnotationDelete(annotation);
    }
}

void Molecule::clearGeomAnnotations()
{
    AnnotationList::iterator iter = annotations.begin();
    do
    {
        while (iter != annotations.end())
        {
            Annotation *annotation = *iter;
            if (annotation->m_type == Annotation::Geom_error)
            {
                annotations.erase(iter);
                doAnnotationDelete(annotation);
                iter = annotations.begin();
                break;
            }
            else
            {
                ++iter;
            }
        }
    } while (iter != annotations.end());
}

void Molecule::clearAnnotations()
{
    AnnotationList::iterator iter = annotations.begin();
    while (iter != annotations.end())
    {
        Annotation *annotation = *iter;
        ++iter;
        annotationToBeDeleted(this, annotation);
    }
    iter = annotations.begin();
    while (iter != annotations.end())
    {
        Annotation *annotation = *iter;
        ++iter;
        delete annotation;
    }
    std::vector<Annotation*>().swap(annotations); // was annotations.clear();
    annotationDeleted(this);
}

void Molecule::addAnnotation(const std::string &s, float x, float y, float z, int id)
{
    // later we may want this to generate a unique GUID
    if (id == 0)
    {
        id = annotations.size()+1;
    }
    addAnnotation(new Annotation(s.c_str(), x, y, z, id));

}

void Molecule::doAnnotationDelete(Annotation *annotation)
{
    annotationToBeDeleted(this, annotation);
    delete annotation;
    annotationDeleted(this);
}

Molecule::AtomLabelList&Molecule::getAtomLabels()
{
    return atomLabels;
}

void Molecule::unlabelAtom(MIAtom *atom)
{
    ATOMLABEL *label = findLabelForAtom(atom);
    if (label != NULL)
    {
        deleteAtomLabel(label);
    }
}

ATOMLABEL*Molecule::findLabelForAtom(MIAtom *atom)
{
    AtomLabelList::iterator iter = atomLabels.begin();
    while (iter != atomLabels.end())
    {
        ATOMLABEL *label = *iter;
        ++iter;
        if (label->atom() == atom)
        {
            return label;
        }
    }
    return NULL;
}

bool Molecule::isAtomLabeled(MIAtom *atom)
{
    return findLabelForAtom(atom) != NULL;
}

void Molecule::labelAtom(MIAtom *atom, Residue *res)
{
    if (!MIAtom::isValid(atom) || !Monomer::isValid(res))
    {
        return;
    }

    ATOMLABEL *label = findLabelForAtom(atom);
    if (label == NULL)
    {
        label = new ATOMLABEL(res, atom);
        addAtomLabel(label);
    }
    else
    {
        label->visible(true);
        atomLabelChanged(this, label);
    }
}

void Molecule::labelAtomStyle(MIAtom *atom, int style)
{
    if (!MIAtom::isValid(atom))
    {
        return;
    }

    ATOMLABEL *label = findLabelForAtom(atom);
    if (label != NULL)
    {
        label->style(style);
    }
}

void Molecule::labelEveryNthResidue(int n)
{
    int ir = 0;
    for (ResidueListIterator res = residuesBegin(); res != residuesEnd(); ++res)
    {
        ir++;
        if (ir%n == 0)
        {
            labelAtom(atom_default(res), res);
        }
    }
}

void Molecule::addAtomLabel(ATOMLABEL *label)
{
    atomLabels.push_back(label);
    atomLabelAdded(this, label);
}

void Molecule::updateAtomLabels()
{
    AtomLabelList::iterator iter = atomLabels.begin();
    while (iter != atomLabels.end())
    {
        ATOMLABEL *label = *iter;
        ++iter;
        label->label();
        atomLabelChanged(this, label);
    }
}

void Molecule::deleteAtomLabel(ATOMLABEL *label)
{
    AtomLabelList::iterator pos = find(atomLabels.begin(), atomLabels.end(), label);
    if (pos != atomLabels.end())
    {
        ATOMLABEL *lab = *pos;
        atomLabels.erase(pos);
        doAtomLabelDelete(lab);
    }
}

void Molecule::clearAtomLabels()
{
    atomLabelToBeDeleted(this, atomLabels);
    AtomLabelList::iterator iter = atomLabels.begin();
    while (iter != atomLabels.end())
    {
        ATOMLABEL *label = *iter;
        ++iter;
        delete label;
    }
    std::vector<ATOMLABEL*>().swap(atomLabels); // was atomLabels.clear();
    atomLabelDeleted(this);
}

void Molecule::doAtomLabelDelete(ATOMLABEL *label)
{
    AtomLabelList labels;
    labels.push_back(label);
    atomLabelToBeDeleted(this, labels);
    delete label;
    atomLabelDeleted(this);
}

void Molecule::setAtomLabelVisible(ATOMLABEL *label, bool visible)
{
    label->visible(visible);
    atomLabelChanged(this, label);
}

void Molecule::setAtomLabelColor(ATOMLABEL *label, unsigned char red, unsigned char green, unsigned char blue)
{
    label->red(red);
    label->green(green);
    label->blue(blue);
    atomLabelChanged(this, label);
}

void Molecule::setAtomLabelText(ATOMLABEL *label, const char *text)
{
    label->label(text);
    atomLabelChanged(this, label);
}

void Molecule::setDotsColor(int color)
{
    if (dots.size() > 0)
    {
        SurfaceDots::iterator iter = dots.begin();
        while (iter != dots.end())
        {
            (*iter).color = color;
            ++iter;
        }
        surfaceChanged(this);
    }
}

void Molecule::setAtomBValueAndOccupancy(MIAtom *atom, float bvalue, float occ)
{
    doAtomBValueAndOccupancy(atom, bvalue, occ);
    MIAtomList changed;
    changed.push_back(atom);
    atomChanged(this, changed);
}

void Molecule::setAtomsBValueAndOccupancy(MIAtomList atoms, float bvalue, float occ)
{
    MIAtom_iter iter;
    for (iter = atoms.begin(); iter != atoms.end(); ++iter)
    {
        doAtomBValueAndOccupancy(*iter, bvalue, occ);
    }
    atomChanged(this, atoms);
}

void Molecule::doAtomBValueAndOccupancy(MIAtom *atom, float bvalue, float occ)
{
    atom->setBValue(bvalue);
    atom->setOcc(occ);
}

void Molecule::doAtomColor(MIAtom *atom, int color)
{
    atom->setColor(color);
}

void Molecule::setAtomColor(MIAtom *atom, int color)
{
    doAtomColor(atom, color);
    MIAtomList changed;
    changed.push_back(atom);
    atomChanged(this, changed);
}

void Molecule::setAtomsColor(MIAtomList atoms, int color)
{
    MIAtom_iter iter;
    for (iter = atoms.begin(); iter != atoms.end(); ++iter)
    {
        doAtomColor(*iter, color);
    }
    atomChanged(this, atoms);
}

void Molecule::toggleChainHidden(Residue *chain)
{
    Residue *residue = chain;
    int chain_id = chain->chain_id();
    MIAtomList changed;
    while ((residue != NULL) && chain_id == residue->chain_id())
    {
        for (int i = 0; i < residue->atomCount(); i++)
        {
            residue->atom(i)->setColor(-residue->atom(i)->color());
            changed.push_back(residue->atom(i));
        }
        residue = residue->next();
    }
    atomChanged(this, changed);
}

void Molecule::toggleResidueHidden(Residue *residue)
{
    MIAtomList changed;
    for (int i = 0; i < residue->atomCount(); i++)
    {
        residue->atom(i)->setColor(-residue->atom(i)->color());
        changed.push_back(residue->atom(i));
    }
    atomChanged(this, changed);
}

void Molecule::toggleResiduesHidden(std::vector<Residue*> &residues)
{
    MIAtomList changed;
    for (size_t j = 0; j< residues.size(); ++j)
    {
        Residue *residue = residues[j];
        for (int i = 0; i < residue->atomCount(); i++)
        {
            residue->atom(i)->setColor(-residue->atom(i)->color());
            changed.push_back(residue->atom(i));
        }
    }
    atomChanged(this, changed);
}

void Molecule::toggleAtomHidden(MIAtom *atom)
{
    atom->setColor(-atom->color());
    MIAtomList changed;
    changed.push_back(atom);
    atomChanged(this, changed);
}

void Molecule::toggleAtomsHidden(MIAtomList &atoms)
{
    MIAtomList changed;
    for (size_t i = 0; i < atoms.size(); ++i)
    {
        MIAtom *atom = atoms[i];
        atom->setColor(-atom->color());
        changed.push_back(atom);
    }
    atomChanged(this, changed);
}


void Molecule::setResidueColor(Residue *residue, int color, int colorMethod)
{
    MIAtomList changed;
    for (int i = 0; i < residue->atomCount(); i++)
    {
        MIAtom *atom = residue->atom(i);
        doAtomColor(atom, getcolor(residue, atom, true, color, colorMethod, "*"));
        changed.push_back(atom);
    }
    atomChanged(this, changed);
}

void Molecule::setResiduesColor(std::vector<Residue*> &residues, int color, int colorMethod)
{
    MIAtomList changed;

    for (size_t j = 0; j< residues.size(); ++j)
    {
        Residue *residue = residues[j];
        for (int i = 0; i < residue->atomCount(); i++)
        {
            MIAtom *atom = residue->atom(i);
            doAtomColor(atom, getcolor(residue, atom, true, color, colorMethod, "*"));
            changed.push_back(atom);
        }
    }
    atomChanged(this, changed);
}

void Molecule::setChainColor(Residue *chain, int color, int colorMethod)
{
    MIAtomList changed;
    Residue *residue = chain;
    short chain_id = chain->chain_id();
    while ((residue != NULL) && residue->chain_id() == chain_id)
    {
        for (int i = 0; i < residue->atomCount(); i++)
        {
            MIAtom *atom = residue->atom(i);
            doAtomColor(atom, getcolor(residue, atom, true, color, colorMethod, "*"));
            changed.push_back(atom);
        }
        residue = residue->next();
    }
    atomChanged(this, changed);
}

void Molecule::setColor(int color, int colorMethod)
{
    MIAtomList changed;
    for (ResidueListIterator res = residuesBegin(); res != residuesEnd(); ++res)
    {
        for (int i = 0; i < res->atomCount(); i++)
        {
            MIAtom *atom = res->atom(i);
            doAtomColor(atom, getcolor(res, atom, true, color, colorMethod, "*"));
            changed.push_back(atom);
        }
    }
    atomChanged(this, changed);
}

void Molecule::Show()
{
    visible = 1;
    moleculeChanged(this);
}

void Molecule::Hide()
{
    visible = 0;
    moleculeChanged(this);
}

void Molecule::FixHeaders(std::vector<std::string> &headers)
{
    std::string s;
    int i = 0;
    bool foundcryst = false;
    bool foundscale = false;
    bool foundorigin = false;
    bool foundcompound = false;
    float zero = 0;
    float one = 1.0F;

    i = 0;
    std::vector<std::string>::iterator iter = headers.begin();
    while (iter != headers.end())
    {
        std::string s = std::string(iter->c_str());
        ++iter;
        if (strncmp(s.c_str(), "SCALE", 5)==0)
        {
            foundscale = true;
        }
        if (strncmp(s.c_str(), "COMPND", 6)==0)
        {
            foundcompound = true;
        }
        if (strncmp(s.c_str(), "ORIGX", 5)==0)
        {
            foundorigin = true;
        }
        if (strncmp(s.c_str(), "CRYST1", 6)==0)
        {
            foundcryst = true;
        }
    }
    if (!foundcompound)
    {
        std::string cmpnd, buf;
        char cont = ' ';
        int rec_no = 1;
        buf = std::string(compound.c_str());
        cmpnd = MIBeforeFirst(buf, ';');
        buf = MIAfterFirst(buf, ';');
        while (!cmpnd.empty())
        {
            if (rec_no == 1)
            {
                cont = ' ';
            }
            else
            {
                cont = (char)rec_no+'0';
            }
            s = format("COMPND   %c %s;", cont, cmpnd.c_str());
            headers.push_back(s.c_str());
            cmpnd = MIBeforeFirst(buf, ';');
            buf = MIAfterFirst(buf, ';');
        }
    }
    if (!foundcryst)
    {
        s = format("CRYST1%9.3f%9.3f%9.3f%7.2f%7.2f%7.2f %-11s%4i          ",
                   mapheader->a, mapheader->b, mapheader->c, mapheader->alpha, mapheader->beta,
                   mapheader->gamma, mapheader->GetHMSymbol(), mapheader->nsym);
        headers.push_back(s.c_str());
    }
    if (!foundorigin)
    {
        s = format("ORIGX1    %10.6f%10.6f%10.6f     %10.5f", one, zero, zero, zero);
        headers.push_back(s.c_str());
        s = format("ORIGX2    %10.6f%10.6f%10.6f     %10.5f", zero, one, zero, zero);
        headers.push_back(s.c_str());
        s = format("ORIGX3    %10.6f%10.6f%10.6f     %10.5f", zero, zero, one, zero);
        headers.push_back(s.c_str());
    }
    if (!foundscale)
    {
        s = format("SCALE1    %10.6f%10.6f%10.6f     %10.5f",
                   mapheader->ctof[0][0], mapheader->ctof[0][1], mapheader->ctof[0][2], zero);
        headers.push_back(s.c_str());
        s = format("SCALE2    %10.6f%10.6f%10.6f     %10.5f",
                   mapheader->ctof[1][0], mapheader->ctof[1][1], mapheader->ctof[1][2], zero);
        headers.push_back(s.c_str());
        s = format("SCALE3    %10.6f%10.6f%10.6f     %10.5f",
                   mapheader->ctof[2][0], mapheader->ctof[2][1], mapheader->ctof[2][2], zero);
        headers.push_back(s.c_str());
    }
}

int Molecule::SaveSymmMolecule(MIAtom *symatom, FILE *fp)
{
    /* write a model at the symm position in atom to file fp */
    ResidueListIterator res = residuesBegin();
    if (res == residuesEnd())
    {
        return 0;
    }
    int nres = 0;
    float tx = symatom->dx();
    float ty = symatom->dy();
    float tz = symatom->dz();
    float x1, y1, z1;
    //float symmat[3][4];
    int isym = symatom->symmop(); // symm number saved in the lowest byte
    int j;

    if (isym > mapheader->nsym)
    {
        Logger::message("bad symm atom");
        return 0;
    }
    /* copy the symmops into symmat
       for(j=0;j<3;j++)
        for(k=0;k<4;k++)
            symmat[j][k]=mapheader->symops[j][k][isym];*/

    fprintf(fp, "REM Symmetry molecule at %3.0f %3.0f %3.0f symmop # %d\n",
            tx, ty, tz, isym);
    for (; res != residuesEnd(); ++res)
    {
        if (res->atomCount() > 0)
        {
            Residue *sres = new Residue(*res);
            //for(j=0;j<MAXNAME;j++)if(sres.name[j]=='#')sres.name[j]=0;
            for (j = 0; j < res->atomCount(); j++)
            {
                x1 = res->atom(j)->x();
                y1 = res->atom(j)->y();
                z1 = res->atom(j)->z();
                transform(mapheader->ctof, &x1, &y1, &z1);
                mapheader->symm_mh(x1, y1, z1, &x1, &y1, &z1, isym);
                x1 += tx;
                y1 += ty;
                z1 += tz;
                transform(mapheader->ftoc, &x1, &y1, &z1);
                sres->atom(j)->setPosition(x1, y1, z1);
            }
            /* now write it to disk */
            SavePDB(fp, sres, NULL, 0, false);
            delete sres;
        }
        nres++;
    }
    fwrite("TER\n", sizeof(char), 4, fp);
    fwrite("END\n", sizeof(char), 4, fp);
    return nres;
}

void Molecule::Load(FILE *fp)
{
    char buf[1024];
    int r, n;
    Residue *res;
    MIAtom *atom;
    char resid[20], atomid[20], c;
    while (fgets(buf, sizeof buf, fp))
    {
        if (!strncmp(buf, "atomlabel", 9))
        {
            r = sscanf(buf, "%*s%d%s%s %c", &n, resid, atomid, &c);
            if (r > 2 && n == modelnumber)
            {
                if (r == 4)
                {
                    // chainid found
                    if ((res = residue_from_name(residues, resid, c)) != NULL)
                    {
                        if ((atom = atom_from_name(atomid, *res)) != NULL)
                        {
                            labelAtom(atom, res);
                        }
                    }
                }
                else
                {
                    // chain id = space
                    if ((res = residue_from_name(residues, resid, ' ')) != NULL)
                    {
                        if ((atom = atom_from_name(atomid, *res)) != NULL)
                        {
                            labelAtom(atom, res);
                        }
                    }
                }
            }
        }
        else if (!strncmp(buf, "buildlinks", 9))
        {
            if (sscanf(buf, "%*s%d", &n) == 1)
            {
                if (n == modelnumber)
                {
                    BuildLinks();
                }
            }
        }
        else if (!strncmp(buf, "colorinfo", 9))
        {
            r = sscanf(buf, "%*s%d", &n);
            if (n == modelnumber)
            {
                res = residues;
                while (Monomer::isValid(res) && fgets(buf, sizeof buf, fp))
                {
                    if (!strncmp(buf, "endcolor", 8))
                    {
                        break;
                    }
                    if (!strncmp(buf, "res ", 4))
                    {
                        read_colors(res, buf, sizeof buf);
                    }
                    res = res->next();
                }
            }
        }
        else if (!strncmp(buf, "radiusinfo", 9))
        {
            r = sscanf(buf, "%*s%d", &n);
            if (n == modelnumber)
            {
                res = residues;
                while (Monomer::isValid(res) && fgets(buf, sizeof buf, fp))
                {
                    if (!strncmp(buf, "endradius", 8))
                    {
                        break;
                    }
                    if (!strncmp(buf, "res ", 4))
                    {
                        read_radii(res, buf, sizeof buf);
                    }
                    res = res->next();
                }
            }
        }
    }
}

void Molecule::PurgeAllAtoms()
{
    clearAtomLabels();
    MIMoleculeBase::PurgeAllAtoms();
}

void Molecule::PurgeAtom(MIAtom *atom)
{
    unlabelAtom(atom);
    MIMoleculeBase::PurgeAtom(atom);
}

void Molecule::PurgeSymmetryAtom(MIAtom *atom)
{
    unlabelAtom(atom);
    MIMoleculeBase::PurgeSymmetryAtom(atom);
}

void Molecule::PurgeReferences(MIAtom *atom)
{
    unlabelAtom(atom);
    MIMoleculeBase::PurgeReferences(atom);
}

static MolPrefsHandler *MPH = 0;
void MISetMolPrefsHandler(MolPrefsHandler *mph)
{
    MPH = mph;
}

void Molecule::updateFixChainOptions(bool *breakByDiscontinuity, bool *breakByNonpeptide)
{
    if (MPH)
    {
        (*MPH)(breakByDiscontinuity, breakByNonpeptide);
    }
}

bool Molecule::symmAtomsVisible() const
{
    return symmatoms_visible;
}

void Molecule::setSymmAtomsVisible(bool value)
{
    symmatoms_visible = value;
}
