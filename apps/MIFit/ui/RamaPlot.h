#ifndef mifit_ui_RamaPlot_h
#define mifit_ui_RamaPlot_h

#include <chemlib/MIAtom_fwd.h>

#include "graphlib.h"

namespace chemlib
{
    class RESIDUE;
    class MIMoleculeBase;
}

class GraphWindow;
class RamaPlotDialog;
class RamaDataBase;
class MIGLWidget;

class RamaPlotMgr : public QObject, public GraphPickHandlerBase, public GraphKeyHandlerBase, public GraphMouseoverHandlerBase
{
    Q_OBJECT
private:
    RamaPlotMgr(bool show_allowed = false);
public:
    static RamaPlotMgr *instance();
    virtual ~RamaPlotMgr()
    {
    }

    GraphWindow *getGraphWin();

    bool IsShown();
    void Update(chemlib::MIMoleculeBase *mol,
                chemlib::RESIDUE *focusres,
                std::string modelname);
    void Update(chemlib::RESIDUE *fitres, unsigned int select_type);

    void SetShowAllowed(bool state);
    bool GetShowAllowed() const
    {
        return show_allowed;
    }

    void operator()(const GR_POINT &gr);
    void operator()(int keycode, bool shift);
    void Mouseover(int id);

    void SetView(MIGLWidget *v)
    {
        _view = v;
    }

public slots:
    // object changed slots
    void atomChanged(chemlib::MIMoleculeBase *mol, chemlib::MIAtomList &atom);
    void moleculeChanged(chemlib::MIMoleculeBase *mol);

    // common slot for residuesDeleted, atomsDeleted
    void modelObjectDeleted(chemlib::MIMoleculeBase *mol);
    void moleculeDeleted(chemlib::MIMoleculeBase *mol);

private:
    void CreateContours();
    void GraphResidue(chemlib::RESIDUE *prev,
                      chemlib::RESIDUE *res,
                      chemlib::RESIDUE *next,
                      int replace = -1);
    void CreateData();

    /**
     * Clear all data except contour lines
     */
    void ClearData();


    bool _inited;
    bool show_allowed;
    GraphWindow *_gw;
    chemlib::MIMoleculeBase *_mol;
    chemlib::RESIDUE *_focusres;
    RamaDataBase *ramadat[8];
    GraphColor cols[4];
    std::string _modelname;
    unsigned int contour_end;
    MIGLWidget *_view;
    bool _atom_changed;
    int _last_mouseover_id;

    static RamaPlotMgr *_instance;
};

#endif // ifndef mifit_ui_RamaPlot_h
