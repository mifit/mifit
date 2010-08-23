#include <chemlib/chemlib.h>

#include "SecondaryStructure.h"
#include "RibbonSegment.h"
#include "Helix.h"
#include "RESIDUE.h"
#include <chemlib/Monomer.h>

#include <QtCore/QSettings>

using namespace chemlib;

SecondaryStructure::SecondaryStructure()
{
    m_pRibbonSegmentList = NULL;
    m_pRibbonSegmentLast = NULL;
    m_pSheetList = NULL;
    m_pSheetLast = NULL;
    m_pTurnList = NULL;
    m_pTurnLast = NULL;
    m_pRandomList = NULL;
    m_pRandomLast = NULL;
}

SecondaryStructure::~SecondaryStructure()
{
    DeleteAllRibbonSegments();
    DeleteAllSchematic();
}

#define CLAMP(x, low, high)  (((x) > (high)) ? (high) : (((x) < (low)) ? (low) : (x)))

bool SecondaryStructure::AddRibbonSegment(Residue *pFirstResidue, int nResidues)
{
    RibbonSegment *pNewRSeg = new RibbonSegment();
    QSettings settings;
    double width = settings.value("Secondary Structure/tube/width", 0.8).toDouble();
    double thickness = settings.value("Secondary Structure/tube/thickness", 0.8).toDouble();
    bool square = settings.value("Secondary Structure/tube/square", false).toBool();
    long red = settings.value("Secondary Structure/tube/red", 0).toInt();
    long green = settings.value("Secondary Structure/tube/green", 0).toInt();
    long blue = settings.value("Secondary Structure/tube/blue", 255).toInt();

    pNewRSeg->RibbonProfile(width, thickness, square, 8, 10);
    pNewRSeg->SetColor(
        (unsigned char)CLAMP(red, 0, 255),
        (unsigned char)CLAMP(green, 0, 255),
        (unsigned char)CLAMP(blue, 0, 255));
    if (pNewRSeg->MakeRibbon(pFirstResidue, nResidues))
    {
        if (m_pRibbonSegmentLast == NULL)
        {
            m_pRibbonSegmentList = pNewRSeg;
        }
        else
        {
            m_pRibbonSegmentLast->m_pNext = pNewRSeg;
        }
        m_pRibbonSegmentLast = pNewRSeg;
        return true;
    }
    else
    {
        delete pNewRSeg;
        return false;
    }
}

bool SecondaryStructure::AddSchematic(chemlib::MIMoleculeBase *mol,
                                      std::vector<std::pair<chemlib::Residue*, chemlib::Residue*> > &pHelix,
                                      std::vector<std::pair<chemlib::Residue*, chemlib::Residue*> > &pSheet,
                                      std::vector<std::pair<chemlib::Residue*, chemlib::Residue*> > &pTurn,
                                      std::vector<std::pair<chemlib::Residue*, chemlib::Residue*> > &pRandom)
{

    QSettings settings;
    //First do the helices
    unsigned int i;
    double radius = settings.value("Secondary Structure/helix/radius", 2.5).toDouble();
    long red = settings.value("Secondary Structure/helix/red", 255).toInt();
    long green = settings.value("Secondary Structure/helix/green", 0).toInt();
    long blue = settings.value("Secondary Structure/helix/blue", 0).toInt();
    for (i = 0; i < pHelix.size(); i++)
    {
        Helix *pNewHelix = new Helix(radius);
        if (pNewHelix->MakeHelix(mol, pHelix[i].first, pHelix[i].second))
        {
            pNewHelix->SetColor(
                (unsigned char)CLAMP(red, 0, 255),
                (unsigned char)CLAMP(green, 0, 255),
                (unsigned char)CLAMP(blue, 0, 255));
            m_pHelixList.push_back(pNewHelix);
        }
    }

    //Now do the beta sheets
    double width = settings.value("Secondary Structure/sheet/width", 2.0).toDouble();
    double thickness = settings.value("Secondary Structure/sheet/thickness", 1.0).toDouble();
    bool square = settings.value("Secondary Structure/sheet/square", true).toBool();
    red = settings.value("Secondary Structure/sheet/red", 255).toInt();
    green = settings.value("Secondary Structure/sheet/green", 255).toInt();
    blue = settings.value("Secondary Structure/sheet/blue", 0).toInt();
    for (i = 0; i < pSheet.size(); i++)
    {
        //Ribbon stuff expexts start and number
        int nResidues = 1;
        ResidueListIterator res = Residue::getIterator(pSheet[i].first);
        ResidueListIterator resEnd = Residue::getIterator(pSheet[i].second);
        for (; res != resEnd; ++res)
        {
            nResidues++;
        }

        //Now use ribbons
        RibbonSegment *pNewRSeg = new RibbonSegment();
        pNewRSeg->RibbonProfile(width, thickness, square, 8, 10);
        pNewRSeg->SetColor(
            (unsigned char)CLAMP(red, 0, 255),
            (unsigned char)CLAMP(green, 0, 255),
            (unsigned char)CLAMP(blue, 0, 255));
        if (pNewRSeg->MakeRibbon(pSheet[i].first, nResidues))
        {
            if (m_pSheetLast == NULL)
            {
                m_pSheetList = pNewRSeg;
            }
            else
            {
                m_pSheetLast->m_pNext = pNewRSeg;
            }
            m_pSheetLast = pNewRSeg;
        }
        else
        {
            delete pNewRSeg;
            //return false;
        }
    }

    //Now do the turns
    width = settings.value("Secondary Structure/turn/width", 0.9).toDouble();
    thickness = settings.value("Secondary Structure/turn/thickness", 0.9).toDouble();
    square = settings.value("Secondary Structure/turn/square", false).toBool();
    red = settings.value("Secondary Structure/turn/red", 0).toInt();
    green = settings.value("Secondary Structure/turn/green", 0).toInt();
    blue = settings.value("Secondary Structure/turn/blue", 255).toInt();

    for (i = 0; i < pTurn.size(); i++)
    {
        //Ribbon stuff expexts start and number
        int nResidues = 1;
        ResidueListIterator res = Residue::getIterator(pTurn[i].first);
        ResidueListIterator resEnd = Residue::getIterator(pTurn[i].second);
        for (; res != resEnd; ++res)
        {
            nResidues++;
        }

        //Now use ribbons
        RibbonSegment *pNewRSeg = new RibbonSegment();
        pNewRSeg->RibbonProfile(width, thickness, square, 8, 10);
        pNewRSeg->SetColor(
            (unsigned char)CLAMP(red, 0, 255),
            (unsigned char)CLAMP(green, 0, 255),
            (unsigned char)CLAMP(blue, 0, 255));
        pNewRSeg->SetColor(0, 0, 255);
        if (pNewRSeg->MakeRibbon(pTurn[i].first, nResidues))
        {
            if (m_pTurnLast == NULL)
            {
                m_pTurnList = pNewRSeg;
            }
            else
            {
                m_pTurnLast->m_pNext = pNewRSeg;
            }
            m_pTurnLast = pNewRSeg;
        }
        else //Something wrong in MakeRibbon
        {
            //return false;
            delete pNewRSeg;
        }
    }

    //Now do the random coil
    width = settings.value("Secondary Structure/random/width", 1.8).toDouble();
    thickness = settings.value("Secondary Structure/random/thickness", 0.9).toDouble();
    square = settings.value("Secondary Structure/random/square", false).toBool();
    red = settings.value("Secondary Structure/random/red", 0).toInt();
    green = settings.value("Secondary Structure/random/green", 255).toInt();
    blue = settings.value("Secondary Structure/random/blue", 0).toInt();
    for (i = 0; i < pRandom.size(); i++)
    {
        //Ribbon stuff expexts start and number
        int nResidues = 1;
        ResidueListIterator res = Residue::getIterator(pRandom[i].first);
        ResidueListIterator resEnd = Residue::getIterator(pRandom[i].second);
        for (; res != resEnd; ++res)
        {
            nResidues++;
        }

        //Now use ribbons
        RibbonSegment *pNewRSeg = new RibbonSegment();
        pNewRSeg->RibbonProfile(width, thickness, square, 8, 10);
        pNewRSeg->SetColor(
            (unsigned char)CLAMP(red, 0, 255),
            (unsigned char)CLAMP(green, 0, 255),
            (unsigned char)CLAMP(blue, 0, 255));
        if (pNewRSeg->MakeRibbon(pRandom[i].first, nResidues))
        {
            if (m_pRandomLast == NULL)
            {
                m_pRandomList = pNewRSeg;
            }
            else
            {
                m_pRandomLast->m_pNext = pNewRSeg;
            }
            m_pRandomLast = pNewRSeg;
        }
        else //Something wrong in MakeRibbon
        {
            //return false;
            delete pNewRSeg;
        }
    }

    return true;
}
#undef CLAMP

bool SecondaryStructure::DeleteAllRibbonSegments()
{
    RibbonSegment *pCurrRSeg = m_pRibbonSegmentList;
    while (pCurrRSeg)
    {
        RibbonSegment *pNextRSeg = pCurrRSeg->m_pNext;
        delete pCurrRSeg;
        pCurrRSeg = pNextRSeg;
    }
    m_pRibbonSegmentList = NULL;
    m_pRibbonSegmentLast = NULL;
    return true;
}

bool SecondaryStructure::DeleteAllSchematic()
{
    std::vector<Helix*>::iterator pHelix = m_pHelixList.begin();
    for (; pHelix != m_pHelixList.end(); ++pHelix)
    {
        delete *pHelix;
    }

    RibbonSegment *pCurrRSeg = m_pSheetList;
    while (pCurrRSeg)
    {
        RibbonSegment *pNextRSeg = pCurrRSeg->m_pNext;
        delete pCurrRSeg;
        pCurrRSeg = pNextRSeg;
    }
    m_pSheetList = NULL;
    m_pSheetLast = NULL;

    pCurrRSeg = m_pTurnList;
    while (pCurrRSeg)
    {
        RibbonSegment *pNextRSeg = pCurrRSeg->m_pNext;
        delete pCurrRSeg;
        pCurrRSeg = pNextRSeg;
    }
    m_pTurnList = NULL;
    m_pTurnLast = NULL;

    pCurrRSeg = m_pRandomList;
    while (pCurrRSeg)
    {
        RibbonSegment *pNextRSeg = pCurrRSeg->m_pNext;
        delete pCurrRSeg;
        pCurrRSeg = pNextRSeg;
    }
    m_pRandomList = NULL;
    m_pRandomLast = NULL;

    return true;
}

