#ifndef mifit_ui_SecondaryStructure_h
#define mifit_ui_SecondaryStructure_h

#include <vector>
#include <utility>

#include <core/MIData.h>

class RibbonSegment;
class Helix;
namespace chemlib {
class RESIDUE;
}


class SecondaryStructure {
public:
  SecondaryStructure();
  virtual ~SecondaryStructure();

  bool AddRibbonSegment(chemlib::RESIDUE* pFirstResidue, int nResidues);
  bool AddSchematic(chemlib::MIMoleculeBase *mol,
      std::vector<std::pair<chemlib::RESIDUE*, chemlib::RESIDUE*> >& pHelix,
      std::vector<std::pair<chemlib::RESIDUE*, chemlib::RESIDUE*> >& pSheet,
      std::vector<std::pair<chemlib::RESIDUE*, chemlib::RESIDUE*> >& pTurn,
      std::vector<std::pair<chemlib::RESIDUE*, chemlib::RESIDUE*> >& pRandom);
  bool DeleteAllRibbonSegments();
  bool DeleteAllSchematic();

private:

  friend class GLRenderer;

  MIData options;
  RibbonSegment* m_pRibbonSegmentList;
  RibbonSegment* m_pRibbonSegmentLast;
  std::vector<Helix*> m_pHelixList;
  RibbonSegment* m_pSheetList;
  RibbonSegment* m_pSheetLast;
  RibbonSegment* m_pTurnList;
  RibbonSegment* m_pTurnLast;
  RibbonSegment* m_pRandomList;
  RibbonSegment* m_pRandomLast;
};

#endif
