#ifndef mifit_model_RibbonSpan_h
#define mifit_model_RibbonSpan_h

class RibbonSpan
{
public:
    RibbonSpan();
    virtual ~RibbonSpan();

    //For list control
    RibbonSpan *m_pNext;

    //Methods
    bool MakeSpan(double ptca[4][3], double ptc[4][3], double pto[4][3], bool bFirst);
    int SetProfile(double dx, double dy, int square, int profile_pts, int profile_segs);

private:

    friend class GLRenderer;

    int npr;  // Number for profile points
    int nseg; // Number of segments

    //Ribbon data
#define RibbonSpan_MP 8
#define RibbonSpan_MS 10

    double ribbon_pts[RibbonSpan_MP*RibbonSpan_MS][3];
    double ribbon_norms[RibbonSpan_MP*RibbonSpan_MS][3];

};

#endif // ifndef mifit_model_RibbonSpan_h
