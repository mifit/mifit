#include <cstdio>
#include <cmath>

#include <math/mathlib.h>
#include <chemlib/Monomer.h>

#include "MIMainWindow.h"
#include "MIGLWidget.h"
#include <QDockWidget>

#include "rama.h"
#include "RamaPlot.h"
#include "molw.h"
#include "asplib.h"
#include "Xguicryst.h"



static RichardsonRamaData richardson_ramadata[8];
static OrigRamaData orig_ramadat[1];

using namespace chemlib;

OrigRamaData::OrigRamaData()
{
    _line_count = 3;
    _seg_count[0] = 15;
    _seg_count[1] = 5;
    _seg_count[2] = 4;

    /* draw contour lines of preferred regions
       note - these were taken from a figure found in JMB blown
       up and approximated - better ones are probably available! */
    static float line0_x[] = {-55.0f,  -40.0f,  -40.0f,  -84.0f,
                              -82.0f,  -60.0f,  -35.0f,  -35.0f,
                              -175.0f, -175.0f, -148.0f, -130.0f,
                              -122.0f, -175.0f, -175.0f};
    static float line0_y[] = {180.0f, 150.0f, 100.0f, 67.0f,
                              40.0f,    0.0f, -17.0f, -65.0f,
                              -67.0f,  -20.0f, -35.0f, -31.0f,
                              -19.0f, 40.0f, 180.0f};
    static float line1_x[] = { 58.0f, 35.0f, 35.0f, 58.0f,  58.0f};
    static float line1_y[] = {110.0f, 85.0f, 27.0f, 10.0f, 110.0f};
    static float line2_x[] = {-175.0f, -175.0f,  -65.0f,  -55.0f};
    static float line2_y[] = {-180.0f, -170.0f, -170.0f, -180.0f};
    static float *orig_x[] = { line0_x, line1_x, line2_x };
    static float *orig_y[] = { line0_y, line1_y, line2_y };

    _x = orig_x;
    _y = orig_y;
}

unsigned int OrigRamaData::region(float phi, float psi)
{
    if (MIIsNan(phi) || MIIsNan(psi))
        return 0;

    if (!(phi < -32.0 && phi > -180.0 && psi > -70.0 && psi < 180.0)
        && !(phi < -60.0 && phi > -180.0 && psi > -180.0 && psi < -165.0)
        && !(phi < 60.0 && phi > 30.0 && psi > 10.0 && psi < 110.0))
    {
        return 0;
    }
    return 1;
}

RichardsonRamaData::RichardsonRamaData(unsigned int type, bool preferred)
{
    _type = type;
    if (preferred)
    {
        switch (type)
        {
        case RamaPlotType::General:
            _line_count = sizeof(richardson_general_preferred_x)/sizeof(float);
            _x = richardson_general_preferred_x;
            _y = richardson_general_preferred_y;
            _seg_count = richardson_general_preferred_sizes;
            break;
        case RamaPlotType::Gly:
            _line_count = sizeof(richardson_gly_preferred_x)/sizeof(float);
            _x = richardson_gly_preferred_x;
            _y = richardson_gly_preferred_y;
            _seg_count = richardson_gly_preferred_sizes;
            break;
        case RamaPlotType::Pro:
            _line_count = sizeof(richardson_pro_preferred_x)/sizeof(float);
            _x = richardson_pro_preferred_x;
            _y = richardson_pro_preferred_y;
            _seg_count = richardson_pro_preferred_sizes;
            break;
        case RamaPlotType::PrePro:
            _line_count = sizeof(richardson_prepro_preferred_x)/sizeof(float);
            _x = richardson_prepro_preferred_x;
            _y = richardson_prepro_preferred_y;
            _seg_count = richardson_prepro_preferred_sizes;
            break;
        default:
            _line_count = 0;
            _seg_count = 0;
            break;
        }
    }
    else
    {
        switch (type)
        {
        case RamaPlotType::General:
            _line_count = sizeof(richardson_general_allowed_x)/sizeof(float);
            _x = richardson_general_allowed_x;
            _y = richardson_general_allowed_y;
            _seg_count = richardson_general_allowed_sizes;
            break;
        case RamaPlotType::Gly:
            _line_count = sizeof(richardson_gly_allowed_x)/sizeof(float);
            _x = richardson_gly_allowed_x;
            _y = richardson_gly_allowed_y;
            _seg_count = richardson_gly_allowed_sizes;
            break;
        case RamaPlotType::Pro:
            _line_count = sizeof(richardson_pro_allowed_x)/sizeof(float);
            _x = richardson_pro_allowed_x;
            _y = richardson_pro_allowed_y;
            _seg_count = richardson_pro_allowed_sizes;
            break;
        case RamaPlotType::PrePro:
            _line_count = sizeof(richardson_prepro_allowed_x)/sizeof(float);
            _x = richardson_prepro_allowed_x;
            _y = richardson_prepro_allowed_y;
            _seg_count = richardson_prepro_allowed_sizes;
            break;
        default:
            _line_count = 0;
            _seg_count = 0;
            break;

        }
    }
}

unsigned int RichardsonRamaData::region(float phi, float psi)
{
    if (MIIsNan(phi) || MIIsNan(psi))
        return 0;
    if (fabs(phi) > 180.0 || fabs(psi) > 180.0)
    {
        return 0;
    }
    int i = (int)((phi+180.0f)/2.0f)%180;
    int j = (int)((-psi+180.0f)/2.0f)%180;
    char retval;

    switch (_type)
    {
    case RamaPlotType::General:
        retval = richardson_general_data[j][i];
        break;
    case RamaPlotType::Gly:
        retval = richardson_gly_data[j][i];
        break;
    case RamaPlotType::Pro:
        retval = richardson_pro_data[j][i];
        break;
    case RamaPlotType::PrePro:
        retval = richardson_prepro_data[j][i];
        break;
    default:
        retval = ' ';
        break;
    }

    switch (retval)
    {
    case ' ': return 0;
        break;
    case 'x': return 1;
        break;
    case 'X': return 2;
        break;
    }
    return 0;
}

void RamaPlotMgr::handleKeyPress(int keycode, bool shift)
{
    if (!_view || !Monomer::isValid(_focusres) || !_view->IsFitting())
    {
        return;
    }

    float increment = 5.0f;
    if (shift)
    {
        increment += 5.0f;
    }

    std::string mode;
    switch (keycode)
    {
    case Qt::Key_Up:
        mode = "PSI";
        break;
    case Qt::Key_Down:
        mode = "PSI";
        increment = -increment;
        break;
    case Qt::Key_Left:
        mode = "PHI";
        break;
    case Qt::Key_Right:
        mode = "PHI";
        increment = -increment;
        break;
    default:
        return;
        break;
    }

    //do the torsion via the CMolView's facilities
    _view->FitTorsion(mode.c_str());
    _view->fitmol->RotateTorsion(increment);
    _view->UpdateCurrent();
    _view->ReDraw();
}

void RamaPlotMgr::handlePick(const GR_POINT &gr)
{
    Residue *res = static_cast<Residue*>(gr.data);
    if (Monomer::isValid(res))
    {
        if (_view)
        {
            _view->setFocusResidue(res, true);
        }
    }
}

void RamaPlotMgr::handleMouseOver(int id)
{
    if (id == _last_mouseover_id)
    {
        return;
    }

    bool changed = false;
    char buf[200];
    const std::vector<GR_POINT> &points = _gw->GetData();

    // clear last mouseover label
    if (_last_mouseover_id != -1)
    {
        _gw->graph_labelpoint(_last_mouseover_id, 0);
        changed = true;
    }

    //only add label if not already labeled
    if (!(points[id].flag & asplib::GR_LABEL))
    {
        Residue *res = static_cast<Residue*>(points[id].data);
        if (Monomer::isValid(res))
        {
            char c = (char)res->chain_id()&255;
            sprintf(buf, "%s%c", res->name().c_str(), c);
            _gw->graph_labelpoint(id, buf);
            _last_mouseover_id = id;
            changed = true;
        }
    }
    if (changed)
    {
        _gw->redraw();
    }
}

void RamaPlotMgr::GraphResidue(ResidueListIterator prev,
                               ResidueListIterator res,
                               ResidueListIterator next,
                               int replace)
{
    MIAtom *N, *C, *CA, *Nnext, *Cprev;
    unsigned int typ = 0;
    char buf[2000];

    if (!Monomer::isValid(res) || !Monomer::isValid(prev) || !Monomer::isValid(next))
    {
        return;
    }
    /* find atoms */
    N = NULL;
    CA = NULL;
    C = NULL;
    Nnext = NULL;
    Cprev = NULL;
    for (int i = 0; i < res->atomCount(); i++)
    {
        if (!strcmp(res->atom(i)->name(), "N"))
        {
            N = res->atom(i);
        }
        if (!strcmp(res->atom(i)->name(), "CA"))
        {
            CA = res->atom(i);
        }
        if (!strcmp(res->atom(i)->name(), "C"))
        {
            C = res->atom(i);
        }
    }
    for (int i = 0; i < prev->atomCount(); i++)
    {
        if (!strcmp(prev->atom(i)->name(), "C"))
        {
            Cprev = prev->atom(i);
        }
    }
    for (int i = 0; i < next->atomCount(); i++)
    {
        if (!strcmp(next->atom(i)->name(), "N"))
        {
            Nnext = next->atom(i);
        }
    }

    if (N && C && CA && Cprev && Nnext)
    {
        /* found all the atoms needed */
        float x = CalcAtomTorsion(Cprev, N, CA, C);
        float y = CalcAtomTorsion(N, CA, C, Nnext);

        // classify residue by type
        int region, symbol;
        bool isgly = false;
        typ = RamaPlotType::General;
        if (strcmp(res->type().c_str(), "GLY") == 0)
        {
            typ = RamaPlotType::Gly;
        }
        else if (strcmp(res->type().c_str(), "PRO") == 0)
        {
            typ = RamaPlotType::Pro;
        }
        else if (strcmp(next->type().c_str(), "PRO") == 0)
        {
            typ = RamaPlotType::PrePro;
        }
        region = ramadat[typ]->region(x, y);

        switch (region)
        {
        case 0:
            symbol = asplib::GR_SMALLBOX;
            break;
        case 1:
            symbol = asplib::GR_SMALLBOX;
            break;
        case 2:
            symbol = asplib::GR_SMALLCIRCLE;
            break;
        }

        if (isgly)
        {
            symbol = asplib::GR_SMALLCROSS;
        }

        // set or update point on graph
        int ip = 0;
        if (replace == -1)
        {
            ip = _gw->graph_point(x, y, asplib::GR_BREAK | symbol, res, cols[typ]);
        }
        else
        {
            ip = replace;
            _gw->graph_replacepoint(replace, x, y, asplib::GR_BREAK | symbol);
        }

        // print label if in bad region
        if (region != 2)
        {
            char c = (char)res->chain_id()&255;
            sprintf(buf, "%s%c", res->name().c_str(), c);
            _gw->graph_labelpoint(ip, buf);
            //Logger::log("Suspect phi-psi %0.1f %0.1f %s %s\n", x, y, res->type().c_str(), res->name().c_str());
        }

        // label focusres if in model
        if (Monomer::isValid(_focusres))
        {
            if ((strcmp(_focusres->name().c_str(), res->name().c_str()) == 0)
                && ((res->chain_id()&255) == (_focusres->chain_id()&255)))
            {
                char c = (char)res->chain_id()&255;
                if (c != ' ')
                {
                    sprintf(buf, ">%s_%s_%c<", res->type().c_str(), res->name().c_str(), c);
                }
                else
                {
                    sprintf(buf, ">%s_%s<", res->type().c_str(), res->name().c_str());
                }
                _gw->graph_labelpoint(ip, buf);
            }
        }
    }
}

void RamaPlotMgr::SetShowAllowed(bool state)
{
    if (state == show_allowed)
    {
        return;
    }
    show_allowed = state;

    // re-create plot
    CreateContours();
    CreateData();
    _gw->redraw();
}

void RamaPlotMgr::CreateContours()
{
    _gw->graph_clear();

    // set up data classes
    unsigned int PLOT_TYPE_MAX;
    for (unsigned int i = 0; i < 4; ++i)
    {
        richardson_ramadata[i] = RichardsonRamaData(i, true);
    }
    for (unsigned int i = 4; i < 8; ++i)
    {
        richardson_ramadata[i] = RichardsonRamaData(i-4, false);
    }
    orig_ramadat[0] = OrigRamaData();

    PLOT_TYPE_MAX = 4;
    for (unsigned int i = 0; i < PLOT_TYPE_MAX; ++i)
    {
        ramadat[i] = &richardson_ramadata[i];
    }
    if (show_allowed)
    {
        PLOT_TYPE_MAX = 8;
        for (unsigned int i = 4; i < PLOT_TYPE_MAX; ++i)
        {
            ramadat[i] = &richardson_ramadata[i];
        }
    }

    _gw->graph_clear();
    _gw->graph_title("Click on residue to center it. Show/hide chains from navigation tree.\n"
                     "Right click and drag to zoom. Click [here] to zoom out.");

    // draw contours
    for (unsigned int typ = 0; typ < PLOT_TYPE_MAX; ++typ)
    {
        for (unsigned int line = 0; line < ramadat[typ]->NumLines(); ++line)
        {
            for (unsigned int seg = 0; seg < ramadat[typ]->NumSegs(line)-1; ++seg)
            {
                _gw->graph_point(ramadat[typ]->X(line, seg), ramadat[typ]->Y(line, seg), 0, cols[typ%4], (typ < 4 ? 2 : 1));
            }
            contour_end = _gw->graph_point(ramadat[typ]->X(line, ramadat[typ]->NumSegs(line)-1),
                                           ramadat[typ]->Y(line, ramadat[typ]->NumSegs(line)-1),
                                           asplib::GR_BREAK, cols[typ%4], (typ < 4 ? 2 : 1));
        }
    }


    //draw labels, etc
    _gw->graph_xlabel("Phi");
    _gw->graph_ylabel("Psi");
    _gw->graph_scale((double)-180.0, (double)180.0, (double)-180.0, (double)180.0);
    _gw->graph_xdivisions(6);
    _gw->graph_ydivisions(6);
    _gw->graph_xminor(6);
    _gw->graph_yminor(6);
    _gw->graph_curvelabel(1, "General", asplib::GR_SMALLCIRCLE|asplib::GR_BREAK, cols[0]);
    _gw->graph_curvelabel(2, "Glycine", asplib::GR_SMALLCIRCLE|asplib::GR_BREAK, cols[1]);
    _gw->graph_curvelabel(3, "Proline", asplib::GR_SMALLCIRCLE|asplib::GR_BREAK, cols[2]);
    _gw->graph_curvelabel(4, "Pre-Proline", asplib::GR_SMALLCIRCLE|asplib::GR_BREAK, cols[3]);
}

void RamaPlotMgr::ClearData()
{
    _gw->graph_removepointsafter(contour_end+1);
}

void RamaPlotMgr::CreateData()
{
    if (!MIMoleculeBase::isValid(_mol))
    {
        return;
    }

    ResidueListIterator res = _mol->residuesBegin();
    if (res == _mol->residuesEnd())
        return;

    /* skip first residue */
    ResidueListIterator prev_res2 = res;
    ++res;
    ResidueListIterator prev_res = res;
    ++res;
    for (; res != _mol->residuesEnd(); ++res)
    {
        if (res->chain_id() == prev_res2->chain_id()
            && res->chain_id() == prev_res->chain_id())
        {
            MIAtom *a = atom_default(prev_res);
            if (a && !a->isHidden())
            {
                GraphResidue(prev_res2, prev_res, res);
            }
        }
        prev_res2 = prev_res;
        prev_res = res;
    }
}

void RamaPlotMgr::Update(MIMoleculeBase *mol,
                         Residue *focusres,
                         std::string modelname)
{

    if (_mol != mol)
    {
        if (_mol)
        {
            disconnect(_mol, SIGNAL(atomChanged(chemlib::MIMoleculeBase*, chemlib::MIAtomList&)));
            disconnect(_mol, SIGNAL(moleculeChanged(chemlib::MIMoleculeBase*)));
            disconnect(_mol, SIGNAL(moleculeDeleted(chemlib::MIMoleculeBase*)));
            disconnect(_mol, SIGNAL(residuesDeleted(chemlib::MIMoleculeBase*)));
            disconnect(_mol, SIGNAL(atomsDeleted(chemlib::MIMoleculeBase*)));
        }

        if (MIMoleculeBase::isValid(mol))
        {
            connect(mol, SIGNAL(atomChanged(chemlib::MIMoleculeBase*, chemlib::MIAtomList&)),
                    this, SLOT(atomChanged(chemlib::MIMoleculeBase*, chemlib::MIAtomList&)));
            connect(mol, SIGNAL(moleculeChanged(chemlib::MIMoleculeBase*)),
                    this, SLOT(moleculeChanged(chemlib::MIMoleculeBase*)));
            connect(mol, SIGNAL(moleculeDeleted(chemlib::MIMoleculeBase*)),
                    this, SLOT(moleculeDeleted(chemlib::MIMoleculeBase*)));
            connect(mol, SIGNAL(residuesDeleted(chemlib::MIMoleculeBase*)),
                    this, SLOT(modelObjectDeleted(chemlib::MIMoleculeBase*)));
            connect(mol, SIGNAL(atomsDeleted(chemlib::MIMoleculeBase*)),
                    this, SLOT(modelObjectDeleted(chemlib::MIMoleculeBase*)));
        }
    }

    _mol = mol;
    _modelname = modelname;
    _focusres = focusres;

    ClearData();
    CreateData();
    _gw->redraw();
    _atom_changed = false;

}


bool RamaPlotMgr::IsShown()
{
    return _gw->isVisible();
}

void RamaPlotMgr::Update(Residue *focusres, unsigned int select_type)
{
    const std::vector<GR_POINT> &points = _gw->GetData();

    _focusres = focusres;

    if (!Monomer::isValid(_focusres) || !MIMoleculeBase::isValid(_mol) || _atom_changed
        || (select_type != SINGLERESIDUE
            && select_type != SINGLEATOM))
    {
        //just to do a full re-creation of the plot

        // FIXME: presumably we could go through all residues/atoms looking for
        // AtomType::FITATOM in the atom's type field. only worthwhile if this
        // is heavliy used
        ClearData();
        CreateData();
        _gw->redraw();
        _atom_changed = false;
        return;
    }

    //get residue prev to focus res, and one prev to that
    Residue *prev_res = 0, *prev_res2 = 0, *thisres = 0;
    Residue *next_res = 0, *next_res2 = 0;
    for (ResidueListIterator res = _mol->residuesBegin(); res != _mol->residuesEnd(); ++res)
    {
        if (res == ResidueListIterator(_focusres))
        {
            thisres = res;
            ++res;
            next_res = res;
            ++res;
            next_res2 = res;
            break;
        }
        prev_res2 = prev_res;
        prev_res = res;
    }
    if (!thisres)
    {
        return;
    }

    //properly calc phi/psi for focusres & neighbors
    bool updated = false;
    for (unsigned int i = 0; i < points.size(); ++i)
    {
        Residue *r = static_cast<Residue*>(points[i].data);
        if (Monomer::isValid(r))
        {
            if (r == prev_res)
            {
                GraphResidue(prev_res2, prev_res, thisres, i);
                updated = true;
            }
            if (r == thisres)
            {
                GraphResidue(prev_res, thisres, next_res, i);
                updated = true;
            }
            if (r == next_res)
            {
                GraphResidue(thisres, next_res, next_res2, i);
                updated = true;
            }
        }
    }
    if (updated)
    {
        _gw->redraw();
        return;
    }
}
RamaPlotMgr*RamaPlotMgr::_instance = NULL;

RamaPlotMgr::RamaPlotMgr(bool allowed)
    : _inited(false),
      _gw(NULL),
      _mol(0),
      _focusres(0),
      _modelname("")
{
    _instance = this;

    show_allowed = allowed;
    _last_mouseover_id = -1;

    // init contour colors
    for (unsigned int typ = 0; typ < 4; ++typ)
    {
        switch (typ)
        {
        case RamaPlotType::General:
            cols[typ] = GraphColor(255, 0, 128);
            break;
        case RamaPlotType::Gly:
            cols[typ] = GraphColor(0, 255, 0);
            break;
        case RamaPlotType::Pro:
            cols[typ] = GraphColor(0, 0, 255);
            break;
        case RamaPlotType::PrePro:
            cols[typ] = GraphColor(0, 255, 255);
            break;
        }
    }

    _gw = new GraphWindow(MIMainWindow::instance());

    connect(_gw, SIGNAL(mouseOver(int)), this, SLOT(handleMouseOver(int)));
    connect(_gw, SIGNAL(pick(GR_POINT)), this, SLOT(handlePick(GR_POINT)));
    connect(_gw, SIGNAL(keyPress(int,bool)), this, SLOT(handleKeyPress(int,bool)));

    CreateContours();
}


RamaPlotMgr*RamaPlotMgr::instance()
{
    if (_instance)
        return _instance;
    new RamaPlotMgr(); // sets _instance
    return _instance;
}

GraphWindow*RamaPlotMgr::getGraphWin()
{
    return _gw;
}

// called when show/hide bvalue/occ, or color changed
void RamaPlotMgr::atomChanged(chemlib::MIMoleculeBase*, MIAtomList&)
{
    _atom_changed = true;
    Update(_focusres, 0);
}


void RamaPlotMgr::moleculeChanged(chemlib::MIMoleculeBase *mol)
{
    Update(mol, _focusres, _modelname);
}

void RamaPlotMgr::modelObjectDeleted(chemlib::MIMoleculeBase* )
{
    _focusres = 0;
    Update(_mol, 0, _modelname); // re-create for current molecule
}

void RamaPlotMgr::moleculeDeleted(chemlib::MIMoleculeBase*)
{
    _gw->show();
    _focusres = 0;
    Update(0, 0, ""); // clear
}
