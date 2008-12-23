#include "chemlib.h"

#include "SecondaryStructure.h"
#include "RibbonSegment.h"
#include "Helix.h"
#include "RESIDUE.h"
#include "RESIDUE_.h"

using namespace chemlib;

SecondaryStructure::SecondaryStructure() {
  m_pRibbonSegmentList = NULL;
  m_pRibbonSegmentLast = NULL;
  m_pHelixList = NULL;
  m_pHelixLast = NULL;
  m_pSheetList = NULL;
  m_pSheetLast = NULL;
  m_pTurnList = NULL;
  m_pTurnLast = NULL;
  m_pRandomList = NULL;
  m_pRandomLast = NULL;
}

SecondaryStructure::~SecondaryStructure() {
  DeleteAllRibbonSegments();
  DeleteAllSchematic();
}

#define CLAMP(x, low, high)  (((x) > (high)) ? (high) : (((x) < (low)) ? (low) : (x)))

bool SecondaryStructure::AddRibbonSegment(RESIDUE* pFirstResidue, int nResidues) {
  RibbonSegment* pNewRSeg = new RibbonSegment();
  double width;
  double thickness;
  bool square;
  long red;
  long green;
  long blue;
  MIConfig::Instance()->Read("Secondary Structure/tube/width", &width, 0.8);
  MIConfig::Instance()->Read("Secondary Structure/tube/thickness", &thickness, 0.8);
  MIConfig::Instance()->Read("Secondary Structure/tube/square", &square, false);
  MIConfig::Instance()->Read("Secondary Structure/tube/red", &red, 0);
  MIConfig::Instance()->Read("Secondary Structure/tube/green", &green, 0);
  MIConfig::Instance()->Read("Secondary Structure/tube/blue", &blue, 255);
  
  pNewRSeg->RibbonProfile(width, thickness, square, 8, 10);
  pNewRSeg->SetColor(
      (unsigned char)CLAMP(red, 0, 255),
      (unsigned char)CLAMP(green, 0, 255),
      (unsigned char)CLAMP(blue, 0, 255));
  if (pNewRSeg->MakeRibbon(pFirstResidue, nResidues)) {
    if (m_pRibbonSegmentLast == NULL) {
      m_pRibbonSegmentList = pNewRSeg;
    } else {
      m_pRibbonSegmentLast->m_pNext = pNewRSeg;
    }
    m_pRibbonSegmentLast = pNewRSeg;
    return true;
  } else {
    delete pNewRSeg;
    return false;
  }
}

bool SecondaryStructure::AddSchematic(chemlib::MIMoleculeBase *mol,
      std::vector<std::pair<chemlib::RESIDUE*, chemlib::RESIDUE*> >& pHelix,
      std::vector<std::pair<chemlib::RESIDUE*, chemlib::RESIDUE*> >& pSheet,
      std::vector<std::pair<chemlib::RESIDUE*, chemlib::RESIDUE*> >& pTurn,
      std::vector<std::pair<chemlib::RESIDUE*, chemlib::RESIDUE*> >& pRandom) {
  //First do the helices
  unsigned int i;
  double radius;
  long red;
  long green;
  long blue;
  MIConfig::Instance()->Read("Secondary Structure/helix/radius", &radius, 2.5);
  MIConfig::Instance()->Read("Secondary Structure/helix/red", &red, 255);
  MIConfig::Instance()->Read("Secondary Structure/helix/green", &green, 0);
  MIConfig::Instance()->Read("Secondary Structure/helix/blue", &blue, 0);
  for (i = 0; i < pHelix.size(); i++) {
    Helix* pNewHelix = new Helix(radius);
    pNewHelix->MakeHelix(mol, pHelix[i].first, pHelix[i].second);
    pNewHelix->SetColor(
        (unsigned char)CLAMP(red, 0, 255),
        (unsigned char)CLAMP(green, 0, 255),
        (unsigned char)CLAMP(blue, 0, 255));
    if (m_pHelixLast == NULL) {
      m_pHelixList = pNewHelix;
    } else {
      m_pHelixLast->m_pNext = pNewHelix;
    }
    m_pHelixLast = pNewHelix;
  }

  //Now do the beta sheets
  double width;
  double thickness;
  bool square;
  MIConfig::Instance()->Read("Secondary Structure/sheet/width", &width, 2.0);
  MIConfig::Instance()->Read("Secondary Structure/sheet/thickness", &thickness, 1.0);
  MIConfig::Instance()->Read("Secondary Structure/sheet/square", &square, true);
  MIConfig::Instance()->Read("Secondary Structure/sheet/red", &red, 255);
  MIConfig::Instance()->Read("Secondary Structure/sheet/green", &green, 255);
  MIConfig::Instance()->Read("Secondary Structure/sheet/blue", &blue, 0);
  for (i = 0; i < pSheet.size(); i++) {
    //Ribbon stuff expexts start and number
    int nResidues = 1;
    MIIter<RESIDUE> res = RESIDUE::getIter(pSheet[i].first);
    for (; (bool) res && res != pSheet[i].second; ++res) {
      nResidues++;
    }

    //Now use ribbons
    RibbonSegment* pNewRSeg = new RibbonSegment();
    pNewRSeg->RibbonProfile(width, thickness, square, 8, 10);
    pNewRSeg->SetColor(
        (unsigned char)CLAMP(red, 0, 255),
        (unsigned char)CLAMP(green, 0, 255),
        (unsigned char)CLAMP(blue, 0, 255));
    if (pNewRSeg->MakeRibbon(pSheet[i].first, nResidues)) {
      if (m_pSheetLast == NULL) {
        m_pSheetList = pNewRSeg;
      } else {
        m_pSheetLast->m_pNext = pNewRSeg;
      }
      m_pSheetLast = pNewRSeg;
    }
    else {
      delete pNewRSeg;
      //return false;
    }
  }

  //Now do the turns
  MIConfig::Instance()->Read("Secondary Structure/turn/width", &width, 0.9);
  MIConfig::Instance()->Read("Secondary Structure/turn/thickness", &thickness, 0.9);
  MIConfig::Instance()->Read("Secondary Structure/turn/square", &square, false);
  MIConfig::Instance()->Read("Secondary Structure/turn/red", &red, 0);
  MIConfig::Instance()->Read("Secondary Structure/turn/green", &green, 0);
  MIConfig::Instance()->Read("Secondary Structure/turn/blue", &blue, 255);
  for (i = 0; i < pTurn.size(); i++) {
    //Ribbon stuff expexts start and number
    int nResidues = 1;
    MIIter<RESIDUE> res = RESIDUE::getIter(pTurn[i].first);
    for (; (bool) res && res != pTurn[i].second; ++res) {
      nResidues++;
    }

    //Now use ribbons
    RibbonSegment* pNewRSeg = new RibbonSegment();
    pNewRSeg->RibbonProfile(width, thickness, square, 8, 10);
    pNewRSeg->SetColor(
        (unsigned char)CLAMP(red, 0, 255),
        (unsigned char)CLAMP(green, 0, 255),
        (unsigned char)CLAMP(blue, 0, 255));
    pNewRSeg->SetColor(0, 0, 255);
    if (pNewRSeg->MakeRibbon(pTurn[i].first, nResidues)) {
      if (m_pTurnLast == NULL) {
        m_pTurnList = pNewRSeg;
      } else {
        m_pTurnLast->m_pNext = pNewRSeg;
      }
      m_pTurnLast = pNewRSeg;
    }
    else { //Something wrong in MakeRibbon
      //return false;
      delete pNewRSeg;
    }
  }

  //Now do the random coil
  MIConfig::Instance()->Read("Secondary Structure/random/width", &width, 1.8);
  MIConfig::Instance()->Read("Secondary Structure/random/thickness", &thickness, 0.9);
  MIConfig::Instance()->Read("Secondary Structure/random/square", &square, false);
  MIConfig::Instance()->Read("Secondary Structure/random/red", &red, 0);
  MIConfig::Instance()->Read("Secondary Structure/random/green", &green, 255);
  MIConfig::Instance()->Read("Secondary Structure/random/blue", &blue, 0);
  for (i = 0; i < pRandom.size(); i++) {
    //Ribbon stuff expexts start and number
    int nResidues = 1;
    MIIter<RESIDUE> res = RESIDUE::getIter(pRandom[i].first);
    for (; (bool) res && res != pRandom[i].second; ++res) {
      nResidues++;
    }

    //Now use ribbons
    RibbonSegment* pNewRSeg = new RibbonSegment();
    pNewRSeg->RibbonProfile(width, thickness, square, 8, 10);
    pNewRSeg->SetColor(
        (unsigned char)CLAMP(red, 0, 255),
        (unsigned char)CLAMP(green, 0, 255),
        (unsigned char)CLAMP(blue, 0, 255));
    if (pNewRSeg->MakeRibbon(pRandom[i].first, nResidues)) {
      if (m_pRandomLast == NULL) {
        m_pRandomList = pNewRSeg;
      } else {
        m_pRandomLast->m_pNext = pNewRSeg;
      }
      m_pRandomLast = pNewRSeg;
    }
    else { //Something wrong in MakeRibbon
      //return false;
      delete pNewRSeg;
    }
  }

  return true;
}
#undef CLAMP

bool SecondaryStructure::DeleteAllRibbonSegments() {
  RibbonSegment* pCurrRSeg = m_pRibbonSegmentList;
  while (pCurrRSeg) {
    RibbonSegment* pNextRSeg = pCurrRSeg->m_pNext;
    delete pCurrRSeg;
    pCurrRSeg = pNextRSeg;
  }
  m_pRibbonSegmentList = NULL;
  m_pRibbonSegmentLast = NULL;
  return true;
}

bool SecondaryStructure::DeleteAllSchematic() {
  Helix* pHelix = m_pHelixList;
  while (pHelix) {
    Helix*  pNextHelix = pHelix->m_pNext;
    delete pHelix;
    pHelix = pNextHelix;
  }
  m_pHelixList = NULL;
  m_pHelixLast = NULL;

  RibbonSegment* pCurrRSeg = m_pSheetList;
  while (pCurrRSeg) {
    RibbonSegment* pNextRSeg = pCurrRSeg->m_pNext;
    delete pCurrRSeg;
    pCurrRSeg = pNextRSeg;
  }
  m_pSheetList = NULL;
  m_pSheetLast = NULL;

  pCurrRSeg = m_pTurnList;
  while (pCurrRSeg) {
    RibbonSegment* pNextRSeg = pCurrRSeg->m_pNext;
    delete pCurrRSeg;
    pCurrRSeg = pNextRSeg;
  }
  m_pTurnList = NULL;
  m_pTurnLast = NULL;

  pCurrRSeg = m_pRandomList;
  while (pCurrRSeg) {
    RibbonSegment* pNextRSeg = pCurrRSeg->m_pNext;
    delete pCurrRSeg;
    pCurrRSeg = pNextRSeg;
  }
  m_pRandomList = NULL;
  m_pRandomLast = NULL;

  return true;
}

