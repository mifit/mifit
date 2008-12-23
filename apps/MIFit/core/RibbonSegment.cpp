#include "chemlib.h"
#include "utillib.h"

#include "RibbonSegment.h"
#include "RibbonSpan.h"
#include "RESIDUE.h"
#include "RESIDUE_.h"

using namespace chemlib;

RibbonSegment::RibbonSegment() {
  m_pRibbonSpanList = NULL;
  m_pRibbonSpanLast = NULL;
  m_pNext = NULL;
  for (int i = 0; i < 3; i++) {
    m_bRGB[i] = 255;
  }
  m_bRGB[3] = 0;
}

RibbonSegment::~RibbonSegment() {
  DeleteAllSpans();
}

bool RibbonSegment::DeleteAllSpans() {

  RibbonSpan* pCurrSpan = m_pRibbonSpanList;
  while (pCurrSpan) {
    RibbonSpan* pNextSpan = pCurrSpan->m_pNext;
    delete pCurrSpan;
    pCurrSpan = pNextSpan;
  }
  m_pRibbonSpanList = NULL;
  m_pRibbonSpanLast = NULL;
  return true;
}

bool RibbonSegment::RibbonProfile(double dWidth, double dThick, bool bSquare, int nPoints, int nSteps) {
  m_dWidth = dWidth;
  m_dThick = dThick;
  m_bSquare = bSquare;
  m_nPoints = nPoints;
  m_nSteps = nSteps;
  return true;
}

bool RibbonSegment::SetColor(unsigned char bRed, unsigned char bGreen, unsigned char bBlue) {
  m_bRGB[0] = bRed;
  m_bRGB[1] = bGreen;
  m_bRGB[2] = bBlue;
  return true;
}

bool RibbonSegment::MakeRibbon(RESIDUE* pFirstResidue, int nResidues) {

  //If there is already a ribbon, delete it.
  if (m_pRibbonSpanList) {
    DeleteAllSpans();
  }

  //Make sure there is a valid first residue
  if (!pFirstResidue || (( nResidues != 0) && (nResidues < 2))) {
    return false;
  }

  //Get the name
  m_csFirstResidueName = pFirstResidue->name();
  m_nResidues = 0;

  //Now generate the spans.  We need four residues.
  //Double up the first at the start of the segment and the last at the end
  RESIDUE* pRes[4];
  pRes[0] = pFirstResidue;
  pRes[1] = pFirstResidue;
  pRes[2] = pRes[1]->next();
  if (pRes[2] == NULL) {
    pRes[2] = pRes[1];
  }
  pRes[3] = pRes[2]->next();
  if (pRes[3] == NULL) {
    pRes[3] = pRes[2];
  }
  bool bFirst = true;
  int iResidues = 9999;
  if (nResidues != 0) {
    iResidues = nResidues;
  }
  for (int i = 0; i < iResidues; i++) {
    //Build the data arrays
    double ptca[4][3], ptc[4][3], pto[4][3];
    for (int j = 0; j < 4; j++) {
      //CA
      MIAtom* pAtom = atom_from_name("CA", *pRes[j]);
      if (!pAtom) {
        throw format("no CA atom in residue %s %s %c", pRes[j]->type().c_str(), pRes[j]->name().c_str(), pRes[j]->chain_id());
      }
      ptca[j][0] = pAtom->x();
      ptca[j][1] = pAtom->y();
      ptca[j][2] = pAtom->z();

      //C
      pAtom = atom_from_name("C", *pRes[j]);
      if (!pAtom) {
        throw format("no C atom in residue %s %s %c", pRes[j]->type().c_str(), pRes[j]->name().c_str(), pRes[j]->chain_id());
      }
      ptc[j][0] = pAtom->x();
      ptc[j][1] = pAtom->y();
      ptc[j][2] = pAtom->z();

      //O
      pAtom = atom_from_name("O", *pRes[j]);
      if (!pAtom) {
        throw format("no O atom in residue %s %s %c", pRes[j]->type().c_str(), pRes[j]->name().c_str(), pRes[j]->chain_id());
      }
      pto[j][0] = pAtom->x();
      pto[j][1] = pAtom->y();
      pto[j][2] = pAtom->z();
    }
    RibbonSpan* pSpan = new RibbonSpan();
    if (m_nResidues == 0) {
      pSpan->SetProfile(m_dWidth, m_dThick, m_bSquare, m_nPoints, m_nSteps);
    }
    if (!pSpan->MakeSpan(ptca, ptc, pto, bFirst)) {
      delete pSpan;
      return false;
    }
    bFirst = false;

    //Now add this Segment to the list
    if (m_pRibbonSpanLast == NULL) {
      m_pRibbonSpanList = pSpan;
    } else {
      m_pRibbonSpanLast->m_pNext = pSpan;
    }
    m_pRibbonSpanLast = pSpan;

    m_nResidues++;

    //If the last two residues match, we are done
    if (pRes[2] == pRes[3]) {
      break;
    }

    //Go on to next
    pRes[0] = pRes[1];
    pRes[1] = pRes[2];
    pRes[2] = pRes[3];

    //Don't update the last residue if this is the end of the segement
    if ((nResidues == 0) || (i < nResidues-1)) {
      pRes[3] = pRes[3]->next();
    }

    //Make sure we are not at the end
    if ((pRes[3] == 0) || (!atom_from_name("CA", *pRes[3]))) {
      pRes[3] = pRes[2];
    }

  }

  return true;
}

