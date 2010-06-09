#ifndef mifit_model_RibbonSegment_h
#define mifit_model_RibbonSegment_h

#include <string>

namespace chemlib
{class Residue;
}
class RibbonSpan;

class RibbonSegment
{
public:
    RibbonSegment();
    virtual ~RibbonSegment();

    //For accessing the next ribbon segment
    RibbonSegment *m_pNext;

    //Methods
    bool DeleteAllSpans(void);
    bool RibbonProfile(double dWidth, double dThick, bool bSquare, int nPoints, int nSteps);
    bool MakeRibbon(chemlib::Residue *pFirstResidue, int nResidues);
    bool SetColor(unsigned char bRed, unsigned char bGreen, unsigned char bBlue);

private:

    friend class GLRenderer;

    std::string m_csFirstResidueName;
    int m_nResidues;
    RibbonSpan *m_pRibbonSpanList;
    RibbonSpan *m_pRibbonSpanLast;
    double m_dWidth;
    double m_dThick;
    bool m_bSquare;
    int m_nPoints; //Number of points in cross section
    int m_nSteps; //Numuber of parametric steps per span
    unsigned char m_bRGB[4];

};

#endif // ifndef mifit_model_RibbonSegment_h
