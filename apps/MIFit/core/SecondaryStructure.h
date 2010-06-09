#ifndef mifit_ui_SecondaryStructure_h
#define mifit_ui_SecondaryStructure_h

#include <vector>
#include <utility>

class RibbonSegment;
class Helix;
namespace chemlib
{
    class Residue;
}


class SecondaryStructure
{
public:
    SecondaryStructure();
    virtual ~SecondaryStructure();

    bool AddRibbonSegment(chemlib::Residue *pFirstResidue, int nResidues);
    bool AddSchematic(chemlib::MIMoleculeBase *mol,
                      std::vector<std::pair<chemlib::Residue*, chemlib::Residue*> > &pHelix,
                      std::vector<std::pair<chemlib::Residue*, chemlib::Residue*> > &pSheet,
                      std::vector<std::pair<chemlib::Residue*, chemlib::Residue*> > &pTurn,
                      std::vector<std::pair<chemlib::Residue*, chemlib::Residue*> > &pRandom);
    bool DeleteAllRibbonSegments();
    bool DeleteAllSchematic();

private:

    friend class GLRenderer;

    RibbonSegment *m_pRibbonSegmentList;
    RibbonSegment *m_pRibbonSegmentLast;
    std::vector<Helix*> m_pHelixList;
    RibbonSegment *m_pSheetList;
    RibbonSegment *m_pSheetLast;
    RibbonSegment *m_pTurnList;
    RibbonSegment *m_pTurnLast;
    RibbonSegment *m_pRandomList;
    RibbonSegment *m_pRandomLast;
};

#endif // ifndef mifit_ui_SecondaryStructure_h
