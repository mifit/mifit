#include "LigPostProcessor.h"

#include <conflib/conflib.h>
#include <chemlib/chemlib.h>
#include <chemlib/RESIDUE_.h>

using namespace chemlib;
using namespace std;

void CopyChiral(CHIRALDICT *chiral, const CHIRAL &chiralSource)
{
    chiral->restype[0] = '\0';
    strncpy(chiral->center, chiralSource.center->name(), chemlib::MAXATOMNAME);
    strncpy(chiral->name[0], chiralSource.getAtom1()->name(), chemlib::MAXATOMNAME);
    strncpy(chiral->name[1], chiralSource.getAtom2()->name(), chemlib::MAXATOMNAME);
    strncpy(chiral->name[2], chiralSource.atom3->name(), chemlib::MAXATOMNAME);
    chiral->order = chiralSource.order;
}

void CopyPlane(PLANEDICT *plane, const PLANE &planeSource)
{
    memset(plane, 0, sizeof(PLANEDICT) );
    plane->natoms = planeSource.natoms;
    for (int i = 0; i < plane->natoms && i < MAXPLANE; ++i)
    {
        strncpy(plane->name[i], planeSource.atoms[i]->name(), MAXNAME);
    }
    strncpy(plane->restype, planeSource.res->type().c_str(), MAXNAME);
}

void CopyTorsion(TORSDICT *torsion, const TORSION &torSource)
{
    memset(torsion, 0, sizeof(TORSDICT) );
    strncpy(torsion->type, torSource.type, 11);

    strncpy(torsion->name[0], torSource.getAtom1()->name(), chemlib::MAXATOMNAME);
    strncpy(torsion->name[1], torSource.getAtom2()->name(), chemlib::MAXATOMNAME);
    strncpy(torsion->name[2], torSource.atom3->name(), chemlib::MAXATOMNAME);
    strncpy(torsion->name[3], torSource.atom4->name(), chemlib::MAXATOMNAME);

    memcpy(torsion->ideal, torSource.ideal, sizeof(float[3]) );
    torsion->nideal = torSource.nideal;
    strncpy(torsion->restype, torSource.res->type().c_str(), MAXNAME);
}

/////////////////////////////////////////////////////////////////////////////
// Function:    LigPostProcessor (constructor)
// Purpose:		Determine from the file type, what additional information is needed
// Input:       The LigDictEntry object, and a string with the file extension
// Output:      None
// Requires:	A code for 2D vs 3D input is passed in the "name1" field of
//				the RESIDUE struct.  As with REFMAC files, 'M' signifies 2D
//				and '.' signifies 3D
/////////////////////////////////////////////////////////////////////////////
LigPostProcessor::LigPostProcessor(LigDictEntry &mon, const char *formatString)
    : m_mon(mon),
      genCoords(false),
      guessBondOrders(false),
      genConstraints(false)
{

    char *format = new char[strlen(formatString)+1];
    format[strlen(formatString)] = '\0';
    const char *p;
    char *p2;
    for (p = formatString, p2 = format; *p != '\0'; ++p, ++p2)
    {
        *p2 = tolower(*p);
    }
    if (strcmp(format, "*.pdb") == 0)
    {
        genCoords = false;
        guessBondOrders = true;
        genConstraints = true;

    }
    else if (strcmp(format, "*.mol") == 0)
    {
        genCoords = (mon.res->name1() == 'M') ? true : false;
        guessBondOrders = false;
        genConstraints = true;

    }
    else if (strcmp(format, "*.cif") == 0)
    {
        genCoords = (mon.res->name1() == 'M') ? true : false;
        guessBondOrders = false;
        genConstraints = (mon.res->name1() == 'M') ? true : false;

    }
    else if (strcmp(format, "*.smi") == 0)
    {
        genCoords = true;
        guessBondOrders = false;
        genConstraints = true;
    }
    delete[] format;
}

/////////////////////////////////////////////////////////////////////////////
// Function:    Process
// Purpose:		Generates the additional information for the Dictionary entry
// Input:       Flags set in this LigPostProcessor object
// Output:      Info written to the LigDictEntry object
// Requires:
/////////////////////////////////////////////////////////////////////////////
void LigPostProcessor::Process()
{
    if (genCoords)
    {
        GenerateCoordinates();
    }
    if (guessBondOrders)
    {
        GuessBondOrders();
    }
    if (genConstraints)
    {
        GenerateConstraints();
    }
    else
    {
        GenerateConstraints(false);
    }
}

/////////////////////////////////////////////////////////////////////////////
// Function:    GenerateCoordinates
// Purpose:		Calls library function to generate coordinates for this monomer
// Input:       RESIDUE and bonds from LigDictEntry object
// Output:      Coordinates written to m_mon.res
// Requires:
/////////////////////////////////////////////////////////////////////////////
void LigPostProcessor::GenerateCoordinates()
{
    std::string log;
    conflib::GenerateCoordinates(m_mon.res, m_mon.bonds, log);
}

/////////////////////////////////////////////////////////////////////////////
// Function:    GuessBondOrders
// Purpose:		Calls library function to generate bond orders for this monomer
// Input:       RESIDUE and bonds from LigDictEntry object
// Output:      "order" fields updated in m_mon.bonds
// Requires:
/////////////////////////////////////////////////////////////////////////////
void LigPostProcessor::GuessBondOrders()
{
    chemlib::GuessBondOrders(m_mon.res, m_mon.bonds);
}

/////////////////////////////////////////////////////////////////////////////
// Function:    GenerateConstraints
// Purpose:		Calls library function to generate bondlengths, angles, torsions
//				planes, and chirals for this monomer
// Input:       RESIDUE and bonds from LigDictEntry object
// Output:      Constraints stored in LigDictEntry object in this
// Requires:
/////////////////////////////////////////////////////////////////////////////
void LigPostProcessor::GenerateConstraints()
{
    GenerateConstraints(true);
}

void LigPostProcessor::GenerateConstraints(bool allConstraints)
{

    vector<Bond> bondlengths;
    vector<TORSION> torsions;
    vector<PLANE> planes;
    vector<CHIRAL> chirals;
    vector<ANGLE> angles;

    conflib::GenerateDictionary(m_mon.res, m_mon.bonds, bondlengths, angles, torsions,
                                torsions, planes, chirals);

    if (allConstraints)
    {
        vector<ANGLE>::iterator angle = angles.begin();
        while (angle != angles.end())
        {
            m_mon.angles.push_back(*angle);
            ++angle;
        }

        // convert PLANEs to PLANEDICTs
        vector<PLANE>::const_iterator pln;
        PLANEDICT pd;
        for (pln = planes.begin(); pln != planes.end(); ++pln)
        {
            CopyPlane(&pd, *pln);
            m_mon.planes.push_back(pd);
        }

        // convert CHIRALs to CHIRALDICTs
        vector<CHIRAL>::const_iterator chrl;
        CHIRALDICT cd;
        for (chrl = chirals.begin(); chrl != chirals.end(); ++chrl)
        {
            CopyChiral(&cd, *chrl);
            strncpy(cd.restype, m_mon.res->type().c_str(), MAXNAME);
            m_mon.chirals.push_back(cd);
        }

        //FindChiralCenters(&m_mon.res, m_mon.bonds, m_mon.chirals);
    }

    // convert TORSIONs to TORSDICTs
    vector<TORSION>::const_iterator tor;
    TORSDICT td;
    for (tor = torsions.begin(); tor != torsions.end(); ++tor)
    {
        CopyTorsion(&td, *tor);
        m_mon.torsions.push_back(td);
    }

}

void LigPostProcessor::setDoGenerateCoordinates(bool on)
{
    genCoords = on;
}

void LigPostProcessor::setDoGuessBondOrders(bool on)
{
    guessBondOrders = on;
}

void LigPostProcessor::setDoGenerateConstraints(bool on)
{
    genConstraints = on;
}

