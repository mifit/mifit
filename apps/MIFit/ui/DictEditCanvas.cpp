#include <math/Point3.h>
#include <opengl/Camera.h>
#include <opengl/Frustum.h>


#include <opengl/OpenGL.h>
#include <opengl/Sphere.h>
#include <opengl/StereoView.h>
#include <opengl/Viewpoint.h>
#include <opengl/Viewport.h>
#include <opengl/interact/FieldOfViewZoomCommand.h>
#include <opengl/interact/MouseArcBallOrbitor.h>
#include <opengl/interact/MousePicker.h>
#include <opengl/interact/MouseTranslator.h>
#include <opengl/interact/MouseZoomer.h>

#include <algorithm>

#include <QMenuBar>
#include <QKeyEvent>
#include <QMessageBox>
#include <QTime>
#include <QTimer>
#include <QInputDialog>

#include <nongui/nonguilib.h>
#include <chemlib/chemlib.h>
#include <chemlib/RESIDUE_.h>
#include <conflib/conflib.h>
#include <core/corelib.h>
#include <util/utillib.h>

#include "DictEditAnglePickingRenderable.h"
#include "DictEditAnnotationPickingRenderable.h"
#include "DictEditAtomPickingRenderable.h"
#include "DictEditBondPickingRenderable.h"
#include "DictEditCanvas.h"
#include "DictEditScene.h"
#include "Displaylist.h"
#include "GLRenderer.h"
#include "MIMolIO.h"
#include "Application.h"

#include "molw.h"
#include "ui/MIDialog.h"
#include "ui/DictEditDialog.h"
#include "MIMenuBar.h"
#include "MIMenu.h"
#include "MIEventHandlerMacros.h"

#include "ui_dicteditorform.h"

// for determining click vs move
static QTime mouseDownTime;

using namespace chemlib;
using namespace mi::math;
using namespace mi::opengl::interact;
using namespace mi::opengl;
using namespace std;


const int DictEditCanvas::mouseClickThreshold = 300;

const int DictEditCanvas::ID_POPUP_LINES = 0;
const int DictEditCanvas::ID_POPUP_BALL_LINES = 1;
const int DictEditCanvas::ID_POPUP_STICKS = 2;
const int DictEditCanvas::ID_POPUP_BALL_STICKS = 3;


bool ordering(const Bond &e1, const Bond &e2)
{
    string s1, s2;
    s1 = string(e1.getAtom1()->name()).append(e1.getAtom2()->name());
    s2 = string(e2.getAtom1()->name()).append(e2.getAtom2()->name());
    if (s1.compare(s2) < 0)
    {
        return true;
    }
    return false;
}


DictEditCanvas::DictEditCanvas(DictEditDialog *daddy)
    : QGLWidget(daddy),
      MIEventHandler(this),
      popupMenu(NULL),
      geomrefiner(NULL),
      isRotating(false),
      isZooming(false),
      isPanning(false)
{

    mouseTimer = new QTimer(this);
    mouseTimer->setSingleShot(true);
    connect(mouseTimer, SIGNAL(timeout()), this, SLOT(OnTimer()));


    camera = new Camera();
    frustum = new Frustum();
    stereoView = new StereoView(frustum, camera);
    stereoView->setStereo(MIConfig::Instance()->GetProfileInt("View Parameters", "stereo", 0) != 0);
    stereoView->setHardwareStereo(MIConfig::Instance()->GetProfileInt("View Parameters", "hardwareStereo", 1) != 0);

    scene = new DictEditScene();
    scene->frustum = frustum;
    scene->camera = camera;
    scene->renderer->setFrustum(frustum);
    scene->renderer->setQGLWidget(this);
    scene->renderer->setRenderStyle(RenderStyle::getDefaultBallAndStick());

    cameraMouseOrbitor = new MouseArcBallOrbitor(camera, frustum->getFocalLength());
    cameraMouseOrbitor->setInvertRotation(true);
    scene->cameraMouseOrbitor = cameraMouseOrbitor;

    cameraMouseTranslator = new MouseTranslator(camera, camera, 0.01f);
    scene->mouseTranslator = cameraMouseTranslator;

    zoomCommand = new FieldOfViewZoomCommand(frustum);
    mouseZoomer = new MouseZoomer(zoomCommand, 0.1f);

    AtomStack = new Stack;
    scene->atomStack = AtomStack;
    Models = new Displaylist;
    parent = daddy;

    SetPick((Bond*) NULL);
    SetPick((ANGLE*) NULL);
    SetPick((PLANE*) NULL);
    SetPick((TORSION*) NULL);
    SetPick((CHIRAL*) NULL);

    bondsdirty = true;
    needsrefining = false;
    ShowAngles = false;
    ShowHydrogens = true;
    ShowBonds = false;
    ShowPlanes = false;
    ShowTorsions = false;
    ShowChirals = false;
    working = false;
    initialized = false;

    geomrefiner = MIFitGeomRefiner();
    geomrefiner->lockRefineTarget();

    clearModifiedAngles();

    ((Molecule*)geomrefiner->GetCurrentModel())->DrawAnnotationBox(false);
    createPlanes();

    GetConfs(confs, geomrefiner->GetCurrentModel()->getResidues(), &geomrefiner->dict, geomrefiner->GetCurrentModel());
    conf = 1;

    mousePicker = new MousePicker();
    atomPickingRenderable = new DictEditAtomPickingRenderable(stereoView, Models, frustum);
    bondPickingRenderable = new DictEditBondPickingRenderable(stereoView, geomrefiner->dict.RefiBonds, frustum);
    anglePickingRenderable = new DictEditAnglePickingRenderable(stereoView, geomrefiner->dict.RefiAngles, frustum);
    annotationPickingRenderable = new DictEditAnnotationPickingRenderable(stereoView, Models, frustum, this);

    Bond::copyBondOrders(geomrefiner->GetCurrentModel()->getBonds(), geomrefiner->dict.RefiBonds);
    AutoFindChirals(true);

    createMenus();

    // note: form is a pointer vs a regular data member so that we can
    // prevent wxdr from having to know about it, which would require
    // revising the entire build order so that *all* uic runs were done
    // first.  This is just one undesirable side effect of the
    // inter-dependence of MIFit libraries.
    form = new Ui::DictEditorForm();

    form->setupUi(parent->getFrame());
    // auto-connection of slots doesn't work b/c the parent of the form is not "this"
    connect(form->showChirals, SIGNAL(clicked()), this, SLOT(on_showChirals_clicked()));
    connect(form->showPlanes,  SIGNAL(clicked()), this, SLOT(on_showPlanes_clicked()));
    connect(form->showAngles,  SIGNAL(clicked()), this, SLOT(on_showAngles_clicked()));
    connect(form->showAtomLabels, SIGNAL(clicked()), this, SLOT(on_showAtomLabels_clicked()));
    connect(form->showHydrogens, SIGNAL(clicked()), this, SLOT(on_showHydrogens_clicked()));
    connect(form->showBonds, SIGNAL(clicked()), this, SLOT(on_showBonds_clicked()));

    connect(form->exportButton, SIGNAL(clicked()), this, SLOT(on_exportButton_clicked()));
    connect(form->changeNameButton, SIGNAL(clicked()), this, SLOT(on_changeNameButton_clicked()));
    connect(form->optimizeButton, SIGNAL(clicked()), this, SLOT(on_optimizeButton_clicked()));
    connect(form->resetButton, SIGNAL(clicked()), this, SLOT(on_optimizeButton_clicked()));

    connect(form->conformerSpinbox, SIGNAL(valueChanged(int)), this, SLOT(on_conformerSpinbox_valueChanged(int)));

}

DictEditCanvas::~DictEditCanvas()
{

    geomrefiner->unlockRefineTarget();
    delete cameraMouseOrbitor;
    delete cameraMouseTranslator;
    delete zoomCommand;
    delete mouseZoomer;
    delete atomPickingRenderable;
    delete mousePicker;
    delete stereoView;
    delete camera;
    delete frustum;
    delete scene;
    delete AtomStack;
    AtomStack = NULL;
    modifiedAngles.clear();
    if (popupMenu == NULL)
    {
        delete popupMenu;
        popupMenu = NULL;
    }
    delete form;
}


#define ID_EXPORTBUTTON 10339
#define ID_RENAMEBUTTON 10340
#define ID_REFINE_BUTTON 10341
#define ID_REFINE_SPN 10342
#define ID_RESET_BUTTON 10343
#define ID_VATOMLABELS_CHK 10344
#define ID_VHYDROGENS_CHK 10345
#define ID_VANGLES_CHK 10346
#define ID_VBONDS_CHK 10347
#define ID_VCHIRALS_CHK 10348
#define ID_VPLANES_CHK 10349
#define ID_CONFORMER_SPIN 10350
#define ID_CONFORMER_TEXT 10351
#define ID_DICTEDITPANEL 10352

// Declare menubar functions

#define ID_MENU 10503
#define ID_SETANGLE_MENU 10504
#define ID_CHANGEATOM_MENU 10505
#define ID_RENAMEATOM_MENU 10506
#define ID_REMATOM_MENU 10507
#define ID_REMATOMS_MENU 10508
#define ID_REMHYDROGENS_MENU 10509
#define ID_BONDSETLEN_MENU 10510
#define ID_BONDCHGORDER_MENU 10511
#define ID_BONDREM_MENU 10512
#define ID_INVERTCHI_MENU 10513
#define ID_ADDCHI_MENU 10514
#define ID_REMCHI_MENU 10515
#define ID_NEXTCONFORMER_MENU 10516
#define ID_PREVCONFORMER_MENU 10517
#define ID_GENERATECONFS_MENU 10518
#define ID_PLANEADD_MENU 10519
#define ID_PLANEINCATOM_MENU 10520
#define ID_PLANEINCATOMS_MENU 10521
#define ID_PLANEREMATOM_MENU 10522
#define ID_PLANEREMATOMS_MENU 10523
#define ID_PLANEREM_MENU 10524


void DictEditCanvas::createMenus()
{
    menuBar = new MIMenuBar(parent->getMenuBar());

    BEGIN_EVENT_TABLE(this, ignored);

    MIMenu *item1 = new MIMenu(*this);
    item1->Append( ID_SETANGLE_MENU, "&Set Angle", "" );
    EVT_MENU(ID_SETANGLE_MENU, DictEditCanvas::OnSetAngle);
    menuBar->Append( item1, "&Angles" );

    MIMenu *item2 = new MIMenu(*this);
    item2->Append( ID_CHANGEATOM_MENU, "&Change type", "" );
    EVT_MENU(ID_CHANGEATOM_MENU, DictEditCanvas::OnChangeAtomType);

    item2->Append( ID_RENAMEATOM_MENU, "&Rename", "" );
    EVT_MENU(ID_RENAMEATOM_MENU, DictEditCanvas::OnRenameAtom);

    item2->AppendSeparator();
    item2->Append( ID_REMATOM_MENU, "Remove &One", "" );
    EVT_MENU(ID_REMATOM_MENU, DictEditCanvas::OnRemoveAtom);

    item2->Append( ID_REMATOMS_MENU, "Remove &Several", "" );
    EVT_MENU(ID_REMATOMS_MENU, DictEditCanvas::OnRemoveAtoms);

    item2->Append( ID_REMHYDROGENS_MENU, "Remove &Hydrogens", "" );
    EVT_MENU(ID_REMHYDROGENS_MENU, DictEditCanvas::OnRemoveHydrogens);
    menuBar->Append( item2, "A&toms" );

    MIMenu *item3 = new MIMenu(*this);
    item3->Append( ID_BONDSETLEN_MENU, "&Set length", "" );
    EVT_MENU(ID_BONDSETLEN_MENU, DictEditCanvas::OnSetBondLength);

    item3->Append( ID_BONDCHGORDER_MENU, "&Change order", "" );
    EVT_MENU(ID_BONDCHGORDER_MENU, DictEditCanvas::OnChangeBondOrder);

    item3->AppendSeparator();
    item3->Append( ID_BONDREM_MENU, "&Remove", "" );
    EVT_MENU(ID_BONDREM_MENU, DictEditCanvas::OnDeleteBond);
    menuBar->Append( item3, "&Bonds" );

    MIMenu *item4 = new MIMenu(*this);
    item4->Append( ID_INVERTCHI_MENU, "&Invert Chiral", "" );
    EVT_MENU(ID_INVERTCHI_MENU, DictEditCanvas::OnInvertCenter);

    item4->Append( ID_ADDCHI_MENU, "&Add Chiral", "" );
    EVT_MENU(ID_ADDCHI_MENU, DictEditCanvas::OnAddChiral);

    item4->AppendSeparator();
    item4->Append( ID_REMCHI_MENU, "&Remove Chiral", "" );
    EVT_MENU(ID_REMCHI_MENU, DictEditCanvas::OnRemoveChiral);

    menuBar->Append( item4, "&Chirals" );

    MIMenu *item5 = new MIMenu(*this);
    item5->Append( ID_NEXTCONFORMER_MENU, "&Next Conformer", "" );
    EVT_MENU(ID_NEXTCONFORMER_MENU, DictEditCanvas::OnNextConformer);

    item5->Append( ID_PREVCONFORMER_MENU, "&Previous Conformer", "" );
    EVT_MENU(ID_PREVCONFORMER_MENU, DictEditCanvas::OnPrevConformer);

    item5->Append( ID_GENERATECONFS_MENU, "&Generate Conformers", "" );
    EVT_MENU(ID_GENERATECONFS_MENU, DictEditCanvas::OnGenerateConformers)

    menuBar->Append( item5, "C&onformers" );

    MIMenu *item6 = new MIMenu(*this);
    item6->Append( ID_PLANEADD_MENU, "&Add", "" );
    EVT_MENU(ID_PLANEADD_MENU, DictEditCanvas::OnAddPlane);

    item6->Append( ID_PLANEINCATOM_MENU, "&Include atom", "" );
    EVT_MENU(ID_PLANEINCATOM_MENU, DictEditCanvas::OnIncludeAtomInPlane);


    item6->Append( ID_PLANEINCATOMS_MENU, "I&nclude atoms", "" );
    EVT_MENU(ID_PLANEINCATOMS_MENU, DictEditCanvas::OnIncludeAtomsInPlane);

    item6->AppendSeparator();
    item6->Append( ID_PLANEREMATOM_MENU, "&Remove atom", "" );
    EVT_MENU(ID_PLANEREM_MENU, DictEditCanvas::OnRemovePlane);

    item6->Append( ID_PLANEREMATOMS_MENU, "R&emove atoms", "" );
    EVT_MENU(ID_PLANEREMATOMS_MENU, DictEditCanvas::OnRemoveAtomsFromPlane);

    item6->Append( ID_PLANEREM_MENU, "R&emove Plane", "" );
    EVT_MENU(ID_PLANEREMATOM_MENU, DictEditCanvas::OnRemoveAtomFromPlane)

    menuBar->Append( item6, "&Planes" );


    popupMenu = new MIMenu(*this, this);

    popupMenu->Append(ID_POPUP_LINES, "Line");
    EVT_MENU(ID_POPUP_LINES, DictEditCanvas::OnPopupMenu)

    popupMenu->Append(ID_POPUP_BALL_LINES, "Ball and Line");
    EVT_MENU(ID_POPUP_BALL_LINES, DictEditCanvas::OnPopupMenu)

    popupMenu->Append(ID_POPUP_STICKS, "Stick");
    EVT_MENU(ID_POPUP_STICKS, DictEditCanvas::OnPopupMenu)

    popupMenu->Append(ID_POPUP_BALL_STICKS, "Ball and Stick");
    EVT_MENU(ID_POPUP_BALL_STICKS, DictEditCanvas::OnPopupMenu)

    END_EVENT_TABLE();
}



void DictEditCanvas::initializeForRender()
{
    model = (Molecule*)geomrefiner->GetCurrentModel();
    Models->AddItem(model);
    scene->models = Models;

    model->DrawAnnotationBox(false);
    std::string resshow;
    std::string at("*");
    model->Select(1, 1, 1, 1, resshow, at, NULL, NULL, 0, 0, 0, 0, 1);
    Point3<double> center;
    int atomCount = 0;
    for (MIIter<RESIDUE> res = model->GetResidues(); res; ++res)
    {
        for (int i = 0; i < res->atomCount(); i++)
        {
            MIAtom *a1 = res->atom(i);
            ++atomCount;
            center.x += a1->x();
            center.y += a1->y();
            center.z += a1->z();
        }
    }
    center.x /= atomCount;
    center.y /= atomCount;
    center.z /= atomCount;

    MIIter<RESIDUE> res = model->GetResidues();
    for (int i = 0; i < res->atomCount(); i++)
    {
        res->atom(i)->setType(0);
        res->atom(i)->setColor(model->getcolor(res, res->atom(i), 1, Colors::YELLOW, Colors::COLORC, at));
        CurrentAtoms.push_back(res->atom(i));
    }
    model->FixAtomicNumbers();
    // Should the bond orders be guessed here? For mol files, it messes up bond orders read from file
    //chemlib::GuessBondOrders(res, geomrefiner->dict.RefiBonds);
    UpdateGeom();
    UpdateButtons();
    // Should the bond orders be guessed here? For mol files, it messes up bond orders read from file
    //chemlib::GuessBondOrders(res, model->getBonds());

    Point3<double> pos;
    float maxDistance = 0.0;
    for (MIIter<RESIDUE> res = model->GetResidues(); res; ++res)
    {
        for (int i = 0; i < res->atomCount(); i++)
        {
            MIAtom *a1 = res->atom(i);
            pos.set(a1->x(), a1->y(), a1->z());
            float distance = center.distance(pos);
            if (distance > maxDistance)
            {
                maxDistance = distance;
            }
        }
    }
    // Add margin to distance
    maxDistance += 1.0;

    float fieldOfView = 40.0f;
    double cameraDistance = maxDistance * 2.0 / (2.0 *  tan(toRadians(fieldOfView) * 0.5));
    Vector3<float> eye(center.x, center.y, center.z + cameraDistance);
    camera->setEye(eye);
    camera->lookAt(Vector3<float>(center.x, center.y, center.z));
    frustum->setFieldOfView(fieldOfView);
    frustum->setPerspective(stereoView->isStereo());
    frustum->setFocalLength(cameraDistance);
    cameraMouseOrbitor->setDistanceToTarget(frustum->getFocalLength());


    frustum->setNearClipping(0.01f);
    frustum->setFarClipping(3.0f * cameraDistance);

    updateViewDependentSettings();

    scene->models = Models;
    scene->renderer->setFogEnabled(true);
    scene->renderer->setFogStart(0.9 * cameraDistance);
    scene->renderer->setFogEnd(2.0 * cameraDistance);

    scene->initializeForRender();
}

void DictEditCanvas::render()
{
    scene->renderer->setViewVector(camera->getViewVector());
    stereoView->render(*scene);
    //stereoView->render(*atomPickingRenderable);
    //stereoView->render(*bondPickingRenderable);
    //stereoView->render(*anglePickingRenderable);
    //stereoView->render(*annotationPickingRenderable);
}

bool IsAutoChiral(CHIRAL c)
{
    return c.flags == CHIRAL_AUTO;
}

void DictEditCanvas::AutoFindChirals(bool copyChiralClass)
{

    //This erase-remove combination clears from RefiChirals
    //all (previously) automatically inferred chirals, since
    //these will be regenerated only if they are still present
    geomrefiner->dict.RefiChirals.erase(
        remove_if(geomrefiner->dict.RefiChirals.begin(),
                  geomrefiner->dict.RefiChirals.end(),
                  IsAutoChiral),
        geomrefiner->dict.RefiChirals.end());

    if (geomrefiner->GetCurrentModel() == NULL)
    {
        return;
    }
    RESIDUE *res = geomrefiner->GetCurrentModel()->getResidues();
    std::string list = FindChiralCenters(res, geomrefiner->dict.RefiBonds, copyChiralClass);

    CHIRAL c;
    vector<MIAtom*> nabors;
    if (list.size())
    {
        for (int i = 0; i < res->atomCount(); i++)
        {

            //Check if this is a (tetrahedral) chiral center
            if (res->atom(i)->chiral_class() != CH_TETRAHEDRAL)
            {
                continue;
            }

            //Get the bonded neighbors of this center, and check that
            //there are at least three neighbors
            nabors.clear();
            chemlib::GetNabors(res->atom(i), geomrefiner->dict.RefiBonds, nabors);
            if (nabors.size() < 3)
            {
                continue;
            }

            //Check if there's already a chiral here (could be a user-deleted chiral)
            if (SearchChirals(res->atom(i)) != 0)
            {
                continue;
            }

            //If we reach here, add a new chiral
            c.center = res->atom(i);
            c.setAtom1(nabors[0]);
            c.setAtom2(nabors[1]);
            c.atom3 = nabors[2];
            c.flags = CHIRAL_AUTO;
            geomrefiner->dict.RefiChirals.push_back(c);
        }
    }
}

void DictEditCanvas::UpdateGeom()
{
    if (bondsdirty)
    {
        std::sort(geomrefiner->dict.RefiBonds.begin(), geomrefiner->dict.RefiBonds.end(), ordering);
        bondsdirty = false;
    }
    model->clearAtomLabels();
    size_t i;
    std::string s;

    RESIDUE *res;
    res = model->getResidues();
    if (res == NULL)
    {
        return;
    }
    MIAtom *a1, *a2;

    model->clearAnnotations();
    Models->ClearVus();

    for (int i = 0; i < res->atomCount(); i++)
    {
        model->labelAtom(res->atom(i), res);
        model->labelAtomStyle(res->atom(i), 5);
    }

    AutoFindChirals();

    if (!ShowBonds)
    {
        SetPick((Bond*) NULL);
    }
    else
    {
        for (i = 0; i < geomrefiner->dict.RefiBonds.size(); i++)
        {
            if (geomrefiner->dict.RefiBonds[i].tolerance < 0)
            {
                continue;
            }
            a1 = geomrefiner->dict.RefiBonds[i].getAtom1();
            a2 = geomrefiner->dict.RefiBonds[i].getAtom2();
            if (!ShowHydrogens && (MIAtom::MIIsHydrogen(a1) || MIAtom::MIIsHydrogen(a2)))
            {
                continue;
            }
            s = ::format("%0.3f", AtomDist(*a1, *a2));
            float offset = std::max(a1->getRadius(), a2->getRadius()) * scene->renderer->getRenderStyle().getStickPercent();
            model->addAnnotation(s, (a1->x()+a2->x())/2.0F + offset, (a1->y()+a2->y())/2.0F + offset, (a1->z()+a2->z())/2.0F + offset);
        }
    }

    if (!ShowTorsions)
    {
        SetPick((TORSION*) NULL);
    }
    else
    {
        DrawTorsions();
    }

    if (!ShowAngles)
    {
        SetPick((ANGLE*) NULL);
    }
    else
    {
        for (i = 0; i < geomrefiner->dict.RefiAngles.size(); i++)
        {
            if (geomrefiner->dict.RefiAngles[i].tolerance < 0)
            {
                continue;
            }
            a1 = geomrefiner->dict.RefiAngles[i].getAtom1();
            a2 = geomrefiner->dict.RefiAngles[i].atom3;
            if (!ShowHydrogens && (MIAtom::MIIsHydrogen(a1) || MIAtom::MIIsHydrogen(a2)))
            {
                continue;
            }
            Models->AddLine(a1->x(), a1->y(), a1->z(), a2->x(), a2->y(), a2->z(), 7);
            s = ::format("%0.2f", CalcAtomAngle(*a1, *geomrefiner->dict.RefiAngles[i].getAtom2(), *a2));
            model->addAnnotation(s, (a1->x()+a2->x())/2.0F, (a1->y()+a2->y())/2.0F, (a1->z()+a2->z())/2.0F);
        }
    }

    createPlanes();
    if (!ShowPlanes)
    {
        SetPick((PLANE*) NULL);
    }
    else
    {
        DrawPlanes();
    }

    if (ShowChirals)
    {
        DrawChirals();
    }
    ReDraw();
}

void DictEditCanvas::initializeGL()
{
    glClearColor(0.0f, 0.0f, 0.0f, 0.0f);
}

void DictEditCanvas::paintGL()
{
    if (working)
        return;

    if (!initialized)
    {
        initialized = true;
        initializeForRender();
        return;
    }

    working = true;
    render();
    glFlush();
    working = false;
}


void DictEditCanvas::keyPressEvent(QKeyEvent *event)
{
    event->accept();
    int key = event->key();
    switch (key)
    {
    case Qt::Key_Up:
    case Qt::Key_I:
        //    view->getCamera().zoom(1.1);
        break;
    case Qt::Key_Down:
    case Qt::Key_O:
        //    view->getCamera().zoom(0.9);
        break;
    case Qt::Key_Left:
        //    view->getCamera().moveRight(-0.05);
        break;
    case Qt::Key_Right:
        //    view->getCamera().moveRight(0.05);
        break;
    case Qt::Key_Space:
        geomrefiner->Refine();
        UpdateGeom();
        break;
    case Qt::Key_Backslash: {
        Application::instance()->toggleStereo();
        bool stereo = MIConfig::Instance()->GetProfileInt("View Parameters", "stereo", 0) != 0;
        stereoView->setStereo(stereo);
        if (stereoView->isStereo())
        {
            if (stereoView->isHardwareStereo())
            {
                Logger::log("hardware stereo on");
            }
            else if (stereoView->isCrossEyed())
            {
                Logger::log("cross-eyed stereo on");
            }
            else
            {
                Logger::log("wall-eyed stereo on");
            }
        }
        else
        {
            Logger::log("stereo off");
        }
        break;
    }
    case Qt::Key_A:
        stereoView->setEyeSeparation(stereoView->getEyeSeparation() - 0.1);
        Logger::log("eye separation: %f", stereoView->getEyeSeparation());
        break;
    case Qt::Key_S:
        stereoView->setEyeSeparation(stereoView->getEyeSeparation() + 0.1);
        Logger::log("eye separation: %f", stereoView->getEyeSeparation());
        break;
    case Qt::Key_Q:
        stereoView->setImageSeparation(stereoView->getImageSeparation() - 0.1);
        Logger::log("image separation: %f", stereoView->getImageSeparation());
        break;
    case Qt::Key_W:
        stereoView->setImageSeparation(stereoView->getImageSeparation() + 0.1);
        Logger::log("image separation: %f", stereoView->getImageSeparation());
        break;
    case Qt::Key_X:
        stereoView->setCrossEyed(!stereoView->isCrossEyed());
        if (stereoView->isCrossEyed())
        {
            Logger::log("cross-eyed stereo on");
        }
        else
        {
            Logger::log("wall-eyed stereo on");
        }
        break;
    }

    ReDraw();
}

void DictEditCanvas::OnTimer()
{
    beginDragAction();
}

void DictEditCanvas::mousePressEvent(QMouseEvent *evt)
{
    evt->accept();
    QPoint point = evt->pos();

    if (evt->button() & Qt::LeftButton)
    {
        mouseStart = point;

        if (evt->modifiers() & Qt::ShiftModifier)
        {
            isZooming = true;
        }
        else
        {
            isRotating = true;
        }
        // Delay start of action to allow for clicking
        mouseTimer->start(mouseClickThreshold);
    }

    if (evt->button() & Qt::MidButton)
    {
        mouseStart = point;

        isZooming = true;
        beginDragAction();
        reRender();
    }

    if (evt->button() & Qt::RightButton)
    {
        mouseStart = point;
        mouseDownTime.start();

        isPanning = true;
        // Delay start of action to allow for clicking
        mouseTimer->start(mouseClickThreshold);
    }
}


void DictEditCanvas::mouseMoveEvent(QMouseEvent *evt)
{
    evt->accept();
    QPoint point = evt->pos();

    if (isRotating || isZooming || isPanning)
    {
        int interactWidth = width();
        if (stereoView->isStereo() && !stereoView->isHardwareStereo())
        {
            interactWidth /= 2;
        }
        cameraMouseOrbitor->setControlSize(interactWidth, height());

        if (mouseStart.x() != point.x() || mouseStart.y() != point.y())
        {
            if (mouseTimer->isActive())
            {
                mouseTimer->stop();
                beginDragAction();
            }

            if (isRotating)
            {
                cameraMouseOrbitor->rotate(point.x(), point.y());
            }
            if (isZooming)
            {
                mouseZoomer->zoom(point.x(), point.y());
                updateViewDependentSettings();
            }
            if (isPanning)
            {
                cameraMouseTranslator->translate(point.x(), point.y());
            }
            reRender();
        }
    }
}

void DictEditCanvas::mouseReleaseEvent(QMouseEvent *evt)
{
    QPoint point = evt->pos();
    evt->accept();

    if (evt->button() == Qt::LeftButton)
    {
        if (mouseStart.x() == point.x() && mouseStart.y() == point.y())
        {
            handlePick(evt->pos());
        }
        endDragAction();
        reRender();
    }


    if (evt->button() == Qt::MidButton)
    {
        endDragAction();
        reRender();
    }

    if (evt->button() == Qt::RightButton)
    {
        if (point.x() == mouseStart.x() && point.y() == mouseStart.y())
        {
            if (mouseDownTime.elapsed() < mouseClickThreshold)
            {
                QPoint pos(point.x(), point.y());
                pos = mapToGlobal(pos);
                popupMenu->doExec(pos);
            }
        }

        endDragAction();
        reRender();

        if (!(evt->buttons() & (Qt::LeftButton | Qt::MidButton | Qt::RightButton)))
        {
            if (mouseTimer->isActive())
            {
                mouseTimer->stop();
            }
        }
    }
}


void DictEditCanvas::wheelEvent(QWheelEvent *evt)
{
    evt->accept();
    if (evt->delta() != 0)
    {
        float change = (float) evt->delta();
        frustum->setFieldOfView(frustum->getFieldOfView() + change);
        updateViewDependentSettings();
        reRender();
    }
}


void DictEditCanvas::beginDragAction()
{
    if (isRotating)
    {
        cameraMouseOrbitor->beginRotate(mouseStart.x(), mouseStart.y());
    }
    if (isZooming)
    {
        mouseZoomer->beginZoom(mouseStart.x(), mouseStart.y());
        updateViewDependentSettings();
    }
    if (isPanning)
    {
        cameraMouseTranslator->beginTranslate(mouseStart.x(), mouseStart.y());
    }
    reRender();
}

void DictEditCanvas::endDragAction()
{
    if (isRotating)
    {
        cameraMouseOrbitor->endRotate();
        isRotating = false;
    }
    if (isZooming)
    {
        mouseZoomer->endZoom();
        updateViewDependentSettings();
        isZooming = false;
    }
    if (isPanning)
    {
        cameraMouseTranslator->endTranslate();
        isPanning = false;
    }
    reRender();
}

bool DictEditCanvas::handlePick(const QPoint &point)
{
    bool somethingPicked = false;
    if (!AtomStack->empty() && AtomStack->PickClearBox(point.x(), height() - point.y()))
    {
        AtomStack->Clear();
        ReDraw();
        somethingPicked = true;
    }
    else if (!AtomStack->empty() && AtomStack->PickHideBox(point.x(), height() - point.y()))
    {
        AtomStack->ToggleMinMax();
        ReDraw();
        somethingPicked = true;
    }
    else if (!AtomStack->empty() && AtomStack->PickPopBox(point.x(), height() - point.y()))
    {
        AtomStack->Pop();
        ReDraw();
        somethingPicked = true;
    }
    else
    {
        vector<GLuint> ids = mousePicker->pick(point.x(), point.y(), frustum, annotationPickingRenderable);
        Annotation *annot = NULL;
        if (ids.size() > 0)
        {
            annot = annotationPickingRenderable->getAnnotation(ids[0]);
        }

        if (annot != NULL && planes.find(annot->m_id) != planes.end())
        {
            PLANE *picked = planes.find(annot->m_id)->second;
            if (pickedPlane == picked)
            {
                SetPick((PLANE*)NULL);
            }
            else
            {
                SetPick(picked);
            }

            UpdateGeom();
            somethingPicked = true;
        }
        else
        {
            ids = mousePicker->pick(point.x(), point.y(), frustum, atomPickingRenderable);
            MIAtom *atom = NULL;
            if (ids.size() > 0)
            {
                atom = atomPickingRenderable->getAtom(ids[0]);
            }
            if (atom != NULL)
            {
                Molecule *molecule = 0;
                RESIDUE *res = 0;
                std::list<Molecule*>::iterator molIter = Models->getMolecules().begin();
                while (molIter != Models->getMolecules().end())
                {
                    molecule = *molIter;
                    molIter++;
                    res = residue_from_atom(molecule->getResidues(), atom);
                    if (res != NULL)
                    {
                        break;
                    }
                }
                Models->SetPicked(molecule, res, atom);
            }
            ids = mousePicker->pick(point.x(), point.y(), frustum, bondPickingRenderable);
            Bond *bestbond = NULL;
            if (ids.size() > 0)
            {
                bestbond = bondPickingRenderable->getBond(ids[0]);
            }
            ids = mousePicker->pick(point.x(), point.y(), frustum, anglePickingRenderable);
            ANGLE *bestangle = NULL;
            if (ids.size() > 0)
            {
                bestangle = anglePickingRenderable->getAngle(ids[0]);
            }
            if (bestbond != NULL)
            {
                SetPick(bestbond);
                SetPick(SearchTorsions(bestbond->getAtom1(), bestbond->getAtom2()));
                somethingPicked = true;
            }
            if (bestangle != NULL)
            {
                SetPick(bestangle);
                somethingPicked = true;
            }
            if (atom != NULL)
            {
                AtomStack->Push(atom, Models->GetPickedResidue(), Models->GetPickedMolecule());
                if (SearchChirals(atom))
                {
                    SetPick(SearchChirals(atom));
                }
                somethingPicked = true;
            }
            if (atom != NULL && (Application::instance()->LabelPicks))
            {
                model->labelAtom(atom, Models->GetPickedResidue());
                model->labelAtomStyle(atom, 5);
            }

            PaletteChanged = true;
            UpdateGeom();
        }
    }
    UpdateButtons();
    return somethingPicked;
}

void DictEditCanvas::ReDraw()
{
    updateGL();
}

void DictEditCanvas::resizeGL(int width, int height)
{
    stereoView->setSize(width, height);
    updateViewDependentSettings();
}

void DictEditCanvas::updateViewDependentSettings()
{
    Viewport *viewport = stereoView->getViewport();
    frustum->updateFrustum(*viewport);
    float heightAtTarget = frustum->getFrustumHeight(frustum->getFocalLength());
    float glUnitsPerPixel = heightAtTarget / (float) viewport->getHeight();
    scene->glUnitsPerPixel = glUnitsPerPixel;
    scene->renderer->updateTextScale(glUnitsPerPixel);
    annotationPickingRenderable->updateTextScale(glUnitsPerPixel);
}
void DictEditCanvas::UpdateButtons()
{

    // enable buttons based on size of stack
    menuBar->Enable(ID_REMATOM_MENU, !AtomStack->empty());
    menuBar->Enable(ID_PLANEREMATOM_MENU, !AtomStack->empty() && pickedPlane != NULL);
    menuBar->Enable(ID_PLANEREMATOMS_MENU, AtomStack->size() > 1 && pickedPlane != NULL);
    menuBar->Enable(ID_PLANEREM_MENU, pickedPlane != NULL);
    menuBar->Enable(ID_PLANEADD_MENU, AtomStack->size() > 2);
    menuBar->Enable(ID_PLANEINCATOM_MENU, !AtomStack->empty() && pickedPlane != NULL);
    menuBar->Enable(ID_PLANEINCATOMS_MENU, AtomStack->size() > 1 && pickedPlane != NULL);
    menuBar->Enable(ID_BONDSETLEN_MENU, pickedBond != NULL);
    menuBar->Enable(ID_BONDCHGORDER_MENU, pickedBond != NULL);
    menuBar->Enable(ID_BONDREM_MENU, pickedBond != NULL);
    menuBar->Enable(ID_REMHYDROGENS_MENU, true);
    menuBar->Enable(ID_REMATOMS_MENU, !AtomStack->empty());
    menuBar->Enable(ID_CHANGEATOM_MENU, !AtomStack->empty());
    menuBar->Enable(ID_SETANGLE_MENU, pickedAngle != NULL);
    menuBar->Enable(ID_INVERTCHI_MENU, !AtomStack->empty());
    menuBar->Enable(ID_ADDCHI_MENU, !AtomStack->empty());
    menuBar->Enable(ID_REMCHI_MENU, pickedChiral != NULL && pickedChiral->flags != CHIRAL_DELETED);
    menuBar->Enable(ID_NEXTCONFORMER_MENU, conf < confs.NumberSets()-1);
    menuBar->Enable(ID_PREVCONFORMER_MENU, conf > 1);

    form->changeNameButton->setEnabled(model != NULL && model->getResidues() != NULL);
    form->exportButton->setEnabled(model != NULL && model->getResidues() != NULL);
    form->conformerSpinbox->setRange(1, confs.NumberSets()-1);
    form->conformerSpinbox->setValue(conf);
    form->confCountLabel->setText(::format("of %d", confs.NumberSets()-1).c_str());
}

void DictEditCanvas::on_changeNameButton_clicked()
{
    // change the 3-letter dictionary code of the residue
    if (!model->getResidues())
    {
        return;
    }
    QString s = QInputDialog::getText(this, "Rename Dictionary Entry", "Enter 3-letter name:", QLineEdit::Normal, model->getResidues()->type().c_str());
    if (s.isEmpty())
    {
        return;
    }
    if (s.size())
    {
        if (s.size() > 3)
        {
            QMessageBox::warning(this, "Invalid name", "Name cannot be longer than 3 letters");
            return;
        }
        // Check for collisions
        vector<std::string> names = GetDictResList(MIFitGeomRefiner());
        vector<std::string>::iterator name = names.begin();
        while (name != names.end())
        {
            if (std::string(name->c_str()) == s.toStdString())
            {
                if (QMessageBox::question(this, "Replace dictionary entry?",
                                          "That name already exists!\nDid you want to replace it?",
                                          QMessageBox::Yes | QMessageBox::No, QMessageBox::No) == QMessageBox::No)
                {

                    // parent->Raise();
                    return;
                }
                else
                {
                    break;
                }
            }
            name++;
        }

        // Finalize rename
        model->getResidues()->setType(s.toStdString());
        std::string t = ::format("Dictionary Editor - %s", s.toAscii().constData());
        parent->setWindowTitle(t.c_str());
        ReDraw();
    }
    // parent->Raise();
}

void DictEditCanvas::on_exportButton_clicked()
{
    OnExport(0);
}

void DictEditCanvas::OnExport(const char *optionalFilename)
{
    // This makes use of the FIOS object for writing out various formats.
    map<int, PLANE*>::iterator i, e;

    MIMolInfo mi;
    MIMolIO fio;
    i = planes.begin();
    e = planes.end();
    while (i != e)
    {
        mi.planes.push_back(*i->second);
        i++;
    }
    Molecule *model = Models->CurrentItem();

    if (model)
    {
        mi.res = model->getResidues();
        mi.bonds = geomrefiner->dict.RefiBonds;
        mi.angles = geomrefiner->dict.RefiAngles;
        mi.torsions = geomrefiner->dict.RefiTorsions;
        mi.chirals = geomrefiner->dict.RefiChirals;
        if (optionalFilename != NULL)
        {
            fio.Write(mi, std::string(optionalFilename));
        }
        else
        {
            fio.Write(mi, "");
        }
    }
    UpdateGeom();
}

void DictEditCanvas::OnRemoveAtomsFromPlane()
{
    MIAtom *a = NULL;
    RESIDUE *res = NULL;
    if (pickedPlane == NULL)
    {
        return;
    }
    // What atom are we removing
    if (pickedPlane->natoms - AtomStack->size() < 3)
    {
        if (QMessageBox::question(this, "Delete Plane?",
                                  "Deleting this many atoms implies deleting the plane. Continue?",
                                  QMessageBox::Yes | QMessageBox::No, QMessageBox::No) == QMessageBox::No)
        {

            return;
        }
        OnRemovePlane();
        return;
    }
    while (true)
    {
        AtomStack->Pop(a, res);
        if (a == NULL)
        {
            return;
        }
        removeAtomFromPlane(a, true);
    }
}

void DictEditCanvas::OnRemoveAtomFromPlane()
{
    MIAtom *a = NULL;
    RESIDUE *res = NULL;
    if (pickedPlane == NULL)
    {
        return;
    }
    // What atom are we removing
    AtomStack->Pop(a, res);
    if (a == NULL)
    {
        return;
    }
    removeAtomFromPlane(a, true);
}

bool DictEditCanvas::AtomInPlane(MIAtom *a)
{
    int i;
    if (a == NULL)
    {
        return false;
    }
    if (pickedPlane == NULL)
    {
        return false;
    }
    for (i = 0; i < pickedPlane->natoms; i++)
    {
        if (pickedPlane->atoms[i] == a)
        {
            return true;
        }
    }
    return false;
}

bool DictEditCanvas::removeAtomFromPlane(MIAtom *a, bool atomNotToBeDeleted)
{
    int i, j;
    MIAtom **atoms;
    if (a == NULL)
    {
        return false;
    }
    if (pickedPlane == NULL)
    {
        return false;
    }
    if (!AtomInPlane(a))
    {
        ReDraw();
        return false;
    }

    // Will we still have a plane?
    if (pickedPlane->natoms == 3)
    {
        if (atomNotToBeDeleted)
        {
            if (QMessageBox::question(this, "Delete Plane?",
                                      "Removing this atom implies deleting the plane. Continue?",
                                      QMessageBox::Yes | QMessageBox::No, QMessageBox::No) == QMessageBox::No)
            {
                return false;
            }
        }
        OnRemovePlane();
        return true;
    }
    if ((atoms = (MIAtom**) malloc((pickedPlane->natoms - 1) * sizeof(MIAtom*))) == NULL)
    {
        QMessageBox::warning(this, "Aborting remove", "Memory allocation failed! Aborting remove");
        return false;
    }

    // Copy over the data from the old plane to the new one.
    j = 0;
    for (i = 0; i < pickedPlane->natoms; i++)
    {
        if (a != pickedPlane->atoms[i])
        {
            atoms[j++] = pickedPlane->atoms[i];
        }
    }

    pickedPlane->natoms--;
    pickedPlane->atoms = atoms; //Update our internal list of planes
    // Final updates to the new plane
    lsqplane(*pickedPlane);
    UpdateGeom();
    return true;
}

void DictEditCanvas::DeletePlane(PLANE *p)
{
    map<int, PLANE*>::iterator i, e;
    vector<PLANE>::iterator k, d;
    k = geomrefiner->dict.RefiPlanes.begin();
    d = geomrefiner->dict.RefiPlanes.end();

    //Delete the plane from our internal structure
    free(p->atoms);
    p->atoms = NULL;
    p->natoms = 0;
    i = planes.begin();
    e = planes.end();
    for (; i != e; i++)
    {
        if (i->second == p)
        {
            planes.erase(i);
            break;
        }
    }
    //Delete the plane from the geomrefiner
    for (; k != d; k++)
    {
        if (k->atoms == p->atoms)
        {
            geomrefiner->dict.RefiPlanes.erase(k);
            break;
        }
    }
}

void DictEditCanvas::OnRemovePlane()
{
    if (pickedPlane == NULL)
    {
        return;
    }
    DeletePlane(pickedPlane);
    SetPick((PLANE*) NULL);
    UpdateGeom();
}

void DictEditCanvas::OnAddPlane()
{
    int i;
    float a[3];
    float b[3];
    MIAtom *atom;
    RESIDUE *res = 0;
    PLANE p;
    p.natoms = AtomStack->size();
    if (p.natoms < 3)
    {
        QMessageBox::warning(this, "Must have 3 atoms", "You must have at least 3 atoms selected to add a plane");
        return;
    }
    // Add the atoms to the plane
    if ((p.atoms = (MIAtom**) malloc(p.natoms * sizeof(MIAtom*))) != NULL)
    {
        for (i = 0; i < p.natoms; i++)
        {
            AtomStack->Pop(atom, res);
            p.atoms[i] = atom;
        }
        p.res = res;
    }
    if (p.natoms == 3)
    {
        // Setup vectors for getting the normal
        a[0] = p.atoms[1]->x() - p.atoms[0]->x();
        a[1] = p.atoms[1]->y() - p.atoms[0]->y();
        a[2] = p.atoms[1]->z() - p.atoms[0]->z();
        b[0] = p.atoms[2]->x() - p.atoms[0]->x();
        b[1] = p.atoms[2]->y() - p.atoms[0]->y();
        b[2] = p.atoms[2]->z() - p.atoms[0]->z();
        cross(a, b, p.vm);
        normvect(p.vm, p.vm); //normalize the normal
    }
    else
    {
        lsqplane(p);
    }
    //Add it to the geomrefiner as well
    geomrefiner->dict.RefiPlanes.push_back(p);
    createPlanes();
    UpdateGeom();
}

void Russtransform(float rt[3][3], float tr[3], float *x, float *y, float *z)
{
    float xp, yp, zp;
    xp = rt[0][0]*(*x)+rt[1][0]*(*y)+rt[2][0]*(*z);
    yp = rt[0][1]*(*x)+rt[1][1]*(*y)+rt[2][1]*(*z);
    zp = rt[0][2]*(*x)+rt[1][2]*(*y)+rt[2][2]*(*z);
    *x = xp + tr[0];
    *y = yp + tr[1];
    *z = zp + tr[2];
}

void DictEditCanvas::DrawTorsions()
{
    int color;
    float tran[3];
    float BOXSIZE = 0.125; //Optimizer will do constant propigation on this
    vector<TORSION>::iterator i, e;

    i = geomrefiner->dict.RefiTorsions.begin();
    e = geomrefiner->dict.RefiTorsions.end();
    for (; i != e; i++)
    {
        color = Colors::YELLOW;
        if (pickedTorsion == &*i)
        {
            color = Colors::RED;
        }
        //Where do I translate to
        tran[0] = (i->getAtom2()->x() + i->atom3->x())/2.0;
        tran[1] = (i->getAtom2()->y() + i->atom3->y())/2.0;
        tran[2] = (i->getAtom2()->z() + i->atom3->z())/2.0;

        //Draw this torsion
        //Bottom
        Models->AddLine(tran[0] - BOXSIZE, tran[1] - BOXSIZE, tran[2],
                        tran[0] + BOXSIZE, tran[1] - BOXSIZE, tran[2], color);
        //Top
        Models->AddLine(tran[0] - BOXSIZE, tran[1] + BOXSIZE, tran[2],
                        tran[0] + BOXSIZE, tran[1] + BOXSIZE, tran[2], color);
        //Left
        Models->AddLine(tran[0] - BOXSIZE, tran[1] - BOXSIZE, tran[2],
                        tran[0] - BOXSIZE, tran[1] + BOXSIZE, tran[2], color);
        //Right
        Models->AddLine(tran[0] + BOXSIZE, tran[1] - BOXSIZE, tran[2],
                        tran[0] + BOXSIZE, tran[1] + BOXSIZE, tran[2], color);
    }
}

void DictEditCanvas::DrawChirals()
{
    int color;
    float tran[3];
    float BOXSIZE = 0.25;
    vector<CHIRAL>::iterator i, e;
    //Are there any chirals?
    if (geomrefiner->dict.RefiChirals.size() == 0)
    {
        return;
    }
    i = geomrefiner->dict.RefiChirals.begin();
    e = geomrefiner->dict.RefiChirals.end();
    for (; i != e; i++)
    {
        color = Colors::YELLOW;
        if (i->flags == CHIRAL_DELETED)
        {
            continue;
        }
        if (&*i == pickedChiral)
        {
            color = Colors::RED;
        }
        tran[0] = i->center->x();
        tran[1] = i->center->y();
        tran[2] = i->center->z();

        Models->AddLine(tran[0] - BOXSIZE, tran[1] - BOXSIZE, tran[2],
                        tran[0], tran[1], tran[2] - BOXSIZE, color);
        Models->AddLine(tran[0] - BOXSIZE, tran[1] - BOXSIZE, tran[2],
                        tran[0], tran[1], tran[2] + BOXSIZE, color);
        Models->AddLine(tran[0] - BOXSIZE, tran[1] + BOXSIZE, tran[2],
                        tran[0], tran[1], tran[2] - BOXSIZE, color);
        Models->AddLine(tran[0] - BOXSIZE, tran[1] + BOXSIZE, tran[2],
                        tran[0], tran[1], tran[2] + BOXSIZE, color);
        Models->AddLine(tran[0] + BOXSIZE, tran[1] - BOXSIZE, tran[2],
                        tran[0], tran[1], tran[2] - BOXSIZE, color);
        Models->AddLine(tran[0] + BOXSIZE, tran[1] - BOXSIZE, tran[2],
                        tran[0], tran[1], tran[2] + BOXSIZE, color);
        Models->AddLine(tran[0] + BOXSIZE, tran[1] + BOXSIZE, tran[2],
                        tran[0], tran[1], tran[2] - BOXSIZE, color);
        Models->AddLine(tran[0] + BOXSIZE, tran[1] + BOXSIZE, tran[2],
                        tran[0], tran[1], tran[2] + BOXSIZE, color);
    }
}

void DictEditCanvas::DrawPlanes()
{
    int a, b, color;
    float rot[3][3], tran[3], x1, y1, z1, x2, y2, z2;
    float maxdist = 0.0, curr = 0.0, i, incr;
    std::string name;
    map<int, PLANE*>::iterator currplane;
    map<int, PLANE*>::iterator lastplane;
    currplane = planes.begin();
    lastplane = planes.end();

    while (currplane != lastplane)
    {
        color = Colors::GREEN;
        if (currplane->second == pickedPlane)
        {
            color = Colors::RED;
        }
        rot[0][0] = currplane->second->atoms[1]->x() - currplane->second->atoms[0]->x();
        rot[0][1] = currplane->second->atoms[1]->y() - currplane->second->atoms[0]->y();
        rot[0][2] = currplane->second->atoms[1]->z() - currplane->second->atoms[0]->z();
        rot[2][0] = currplane->second->vm[0];
        rot[2][1] = currplane->second->vm[1];
        rot[2][2] = currplane->second->vm[2];
        normvect(rot[0], rot[0]);
        cross(rot[2], rot[0], rot[1]);

        tran[0] = tran[1] = tran[2] = 0.0;
        for (a = 0; a < currplane->second->natoms; a++)
        {
            tran[0] += currplane->second->atoms[a]->x();
            tran[1] += currplane->second->atoms[a]->y();
            tran[2] += currplane->second->atoms[a]->z();
        }
        tran[0] /= currplane->second->natoms;
        tran[1] /= currplane->second->natoms;
        tran[2] /= currplane->second->natoms;
        // Get the max distance between two atoms
        for (a = 0; a < currplane->second->natoms; a++)
        {
            for (b = a + 1; b < currplane->second->natoms; b++)
            {
                curr = AtomDist(*currplane->second->atoms[a], *currplane->second->atoms[b]);
                if (curr > maxdist)
                {
                    maxdist = curr;
                }
            }
        }
        curr = ceil(maxdist / 2.0);
        if (currplane->second == pickedPlane)
        {
            incr = 0.5;
        }
        else
        {
            incr = 2.0;
        }
        i = -curr;
        while (i <= curr)
        {
            // Vertical Lines
            x1 = x2 = (float) i;
            y1 = -curr;
            y2 = curr;
            z1 = z2 = 0.0;
            Russtransform(rot, tran, &x1, &y1, &z1);
            Russtransform(rot, tran, &x2, &y2, &z2);
            Models->AddLine(x1, y1, z1, x2, y2, z2, color);
            i += incr;
        }
        // Horizontal Lines
        i = -curr;
        while (i <= curr)
        {
            // Vertical Lines
            x1 = -curr;
            x2 = curr;
            y1 = y2 = (float) i;
            z1 = z2 = 0.0;
            Russtransform(rot, tran, &x1, &y1, &z1);
            Russtransform(rot, tran, &x2, &y2, &z2);
            Models->AddLine(x1, y1, z1, x2, y2, z2, color);
            i += incr;
        }
        // Normal
        x1 = x2 = y1 = y2 = z1 = 0.0;
        z2 = 1.0;
        Russtransform(rot, tran, &x1, &y1, &z1);
        // Add the annotation for it
        name = ::format("Plane %d", currplane->first);
        Annotation *ann = new Annotation;
        ann->m_x = x1;
        ann->m_y = y1;
        ann->m_z = z1;
        ann->m_text = name;
        ann->m_type = Annotation::Plane;
        ann->m_id = currplane->first;
        model->addAnnotation(ann);
        Russtransform(rot, tran, &x2, &y2, &z2);
        Models->AddLine(x1, y1, z1, x2, y2, z2, color);
        currplane++;
    }
}

void DictEditCanvas::OnInvertCenter()
{
    // invert the chirality of the current atom
    int c;
    MIAtom *center;
    std::string myerror;
    std::string error;
    center = NULL;
    center = AtomStack->Pop();
    if (center == NULL)
    {
        return; //No atoms on stack

    }
    GeomSaver tmp;
    for (c = 1; c < confs.NumberSets(); c++)
    {
        confs.Restore(c);
        if (!chemlib::InvertChiralCenter(center, geomrefiner->dict.RefiBonds, error) )
        {
            std::string myerror = ::format("Error: %s", error.c_str());
            QMessageBox::warning(this, myerror.c_str(), myerror.c_str());
            confs.Restore(conf);
            UpdateGeom();
            return;
        }
        tmp.Save(CurrentAtoms, model);
    }

    //  A roundabout way of copying one geomsaver to another
    confs.Clear();
    for (c = 1; c < tmp.NumberSets(); ++c)
    {
        tmp.Restore(c);
        confs.Save(CurrentAtoms, model);
    }
    on_conformerSpinbox_valueChanged(conf);
}

void DictEditCanvas::OnSetAngle()
{
    int i, m;
    Bond *edge;
    double newangle, a, b, c;
    if (pickedAngle == NULL)
    {
        return;
    }
    // Ask the user what the angle should be
    std::string s;
    edge = SearchBonds(pickedAngle->getAtom1(), pickedAngle->getAtom2());
    if (!edge)
    {
        QMessageBox::warning(this, "Error", "Unable to find associated bonds for this angle");
        return;
    }
    a = edge->ideal_length;
    edge = SearchBonds(pickedAngle->getAtom2(), pickedAngle->atom3);
    if (!edge)
    {
        QMessageBox::warning(this, "Error", "Unable to find associated bonds for this angle");
        return;
    }
    b = edge->ideal_length;
    c = acos( ((pickedAngle->ideal_angle * pickedAngle->ideal_angle) - (a*a) - (b*b))/(-2*a*b)) * RAD2DEG;

    bool ok;
    double na = QInputDialog::getDouble(this, "Set Angle", "Enter new angle in degrees:", c, 0.0, 360.0, 1, &ok);
    if (!ok)
    {
        return;
    }
    newangle = na;

    c = sqrt( (a*a) + (b*b) - 2*a*b*cos(newangle * DEG2RAD));
    pickedAngle->ideal_angle = c;
    if (!needsrefining)
    {
        QMessageBox::warning(this, "Must Optimize", "The new angle has been set. To apply it finish making your changes to angles, bond lengths and order then hit Optimize");
        SetNeedsRefine();
    }
    m = geomrefiner->dict.RefiAngles.size();
    for (i = 0; i < m; i++)
    {
        if (&(geomrefiner->dict.RefiAngles[i]) == pickedAngle)
        {
            modifiedAngles[i] = true;
            break;
        }
    }
}

MIAtom *AngleIncludesBond(ANGLE *a, Bond *b)
{
    if (a->getAtom1() == b->getAtom1() && a->getAtom2() == b->getAtom2())
    {
        return a->getAtom1();
    }
    if (a->getAtom1() == b->getAtom2() && a->getAtom2() == b->getAtom1())
    {
        return a->getAtom1();
    }
    if (a->atom3 == b->getAtom1() && a->getAtom2() == b->getAtom2())
    {
        return a->atom3;
    }
    if (a->atom3 == b->getAtom2() && a->getAtom2() == b->getAtom1())
    {
        return a->atom3;
    }
    return NULL;
}

void DictEditCanvas::OnSetBondLength()
{
    // set a bond to a user-input value
    if (pickedBond == NULL)
    {
        return;
    }

    bool ok;
    double d = QInputDialog::getDouble(this, "Set Ideal Bond Length", "Enter new bond length:", pickedBond->ideal_length, 0.0, std::numeric_limits<double>::max(), 4, &ok);
    if (ok)
    {
        OnSetBondLength2(d);
        ReDraw();
    }
    // parent->Raise();
}

void DictEditCanvas::OnSetBondLength2(float newideal)
{
    MIAtom *atom;
    MIAtom dummy;
    MIAtom dummy2;
    ANGLE *ang;
    float d_start = pickedBond->ideal_length;
    pickedBond->ideal_length = newideal;
    // go through angles and set their ideal_length to preserve their angle
    for (size_t i = 0; i < geomrefiner->dict.RefiAngles.size(); i++)
    {
        if (modifiedAngles[i])
        {
            continue; //Don't update angles if usermodified
        }
        ang = &(geomrefiner->dict.RefiAngles[i]);
        if (ang->ideal_angle < 0.0)
        {
            continue;
        }
        atom = AngleIncludesBond(ang, pickedBond);
        if (atom)
        {
            float x1, y1, z1, x3, y3, z3;
            float x2 = ang->getAtom2()->x();
            float y2 = ang->getAtom2()->y();
            float z2 = ang->getAtom2()->z();
            if (atom == ang->getAtom1())
            {
                x1 = ang->getAtom1()->x();
                y1 = ang->getAtom1()->y();
                z1 = ang->getAtom1()->z();
                x3 = ang->atom3->x();
                y3 = ang->atom3->y();
                z3 = ang->atom3->z();
            }
            else
            {
                x1 = ang->atom3->x();
                y1 = ang->atom3->y();
                z1 = ang->atom3->z();
                x3 = ang->getAtom1()->x();
                y3 = ang->getAtom1()->y();
                z3 = ang->getAtom1()->z();
            }

            float diff = d_start-pickedBond->ideal_length;
            float dx = (x1-x2)/d_start;
            float dy = (y1-y2)/d_start;
            float dz = (z1-z2)/d_start;
            dx *= diff;
            dy *= diff;
            dz *= diff;
            x1 -= dx;
            y1 -= dy;
            z1 -= dz;
            dummy.setPosition(x1, y1, z1);
            dummy2.setPosition(x3, y3, z3);
            ang->ideal_angle = AtomDist(dummy, dummy2);

        }
    }
    if (!needsrefining)
    {
        QMessageBox::warning(this, "Optimize required", "The new bond length has been set.\nTo apply it finish making your changes to \nangles, bond lengths and order then hit Optimize");
        SetNeedsRefine();
    }
}

float DictEditCanvas::GetBondOrder(vector<Bond> &edges)
{
    float curval = 0.0, retval = 0.0;
    vector<Bond>::iterator b = edges.begin();
    while (b != edges.end())
    {
        curval = GetBondOrder(b->getOrder());
        if (curval == 0)
        {
            QMessageBox::warning(this, "Error", "Unknown bond type. Aborting Change Bond Order\n");
            return -1.0;
        }
        retval += GetBondOrder((&*b));
        b++;
    }
    return retval;
}

unsigned char DictEditCanvas::GetBondOrder(float c)
{
    if (c == 1.0f)
    {
        return (unsigned char) SINGLEBOND;
    }
    else if (c == 2.0f)
    {
        return (unsigned char) DOUBLEBOND;
    }
    else if (c == 3.0f)
    {
        return (unsigned char) TRIPLEBOND;
    }
    else if (c == 1.5f)
    {
        return (unsigned char) PARTIALDOUBLEBOND;
    }
    else if (c == 0.5f)
    {
        return (unsigned char) HYDROGENBOND;
    }
    else   //Assume it's order one, and assign it to be that too
    {
        return NORMALBOND;
    }
}

float DictEditCanvas::GetBondOrder(unsigned char c)
{
    switch (c)
    {
    case NORMALBOND:
        return (float) SINGLEBOND; // Normal bonds are unknown, assume their order 1
    case SINGLEBOND:
        return (float) SINGLEBOND;
    case DOUBLEBOND:
        return (float) DOUBLEBOND;
    case TRIPLEBOND:
        return (float) TRIPLEBOND;
    case PARTIALDOUBLEBOND:
        return 1.5; // so to speak...
    case HYDROGENBOND:
        return .5; // but these will ultimately be ignored
    case IONICBOND:
    case METALLIGANDBOND:
        return 0;
    default: //Assume it's order one, and assign it to be that too
        return 1.0;
    }
}

void DictEditCanvas::OnChangeBondOrder()
{
    int choicei; // choice index
    float oldorder = 0, neworder = 0, a1bo = 0, a2bo = 0;
    std::string order;
    vector<Bond> edges;

    MIAtom *a1, *a2;
    // Check to make sure a bond is selected
    if (pickedBond == NULL)
    {
        return;
    }
    a1 = pickedBond->getAtom1();
    a2 = pickedBond->getAtom2();
    // Ask the user what they want the bond order to be
    QStringList choices;
    choices << "Single" << "Double" << "Triple" << "Partial Double";
    QString item = QInputDialog::getItem(this, "Change bond order", "Change bond order to:", choices);
    choicei = choices.indexOf(item);
    neworder = GetBondOrder((unsigned char) (choicei + 1));
    oldorder = GetBondOrder(pickedBond); //What is bond order for my current bond

    // Check the current bond order
    // TODO: This is the extremely simple case for changing the bond order.
    // Linear search the bonds and increment as needed. This will probably need
    // to be split off into another function that will build a proper MIAtom
    // which would hold a list of bonds in it's own structure.
    SearchBonds(a1, edges);
    a1bo = GetBondOrder(edges); //Calc total bond order for those edges
    edges.clear();
    SearchBonds(a2, edges);
    a2bo = GetBondOrder(edges); //Calc bond order for atom 2
    // Are we increasing the bond order or decreasing it?
    if (neworder == oldorder)
    {
        return;
    }
    Bond *moleculeBond = SearchMoleculeBonds(a1, a2);
    // Increase: Make sure the atoms on each end have a bond free to give
    if (neworder > oldorder  && (neworder > (MaxNumBonds(a1) - a1bo + oldorder)
                                 || neworder > (MaxNumBonds(a2) - a2bo + oldorder)))
    {
        QMessageBox::warning(this, "Error", "Cannot change bond order, not enough free e-");
        return;
    }
    pickedBond->setOrder(GetBondOrder(neworder));
    if (moleculeBond != NULL)
    {
        moleculeBond->setOrder(GetBondOrder(neworder));
    }
    OnSetBondLength2(chemlib::IdealBondLength(*pickedBond)); // Only need the refiner and updater, not the dialog
    UpdateGeom();
}

void DictEditCanvas::OnDeleteBond(Bond *tokill)
{
    // delete a bond
    if (!tokill)
    {
        return;
    }
    vector<Bond>::iterator b = geomrefiner->dict.RefiBonds.begin();
    while (b != geomrefiner->dict.RefiBonds.end())
    {
        if ((&*b) == tokill)
        {
            model->BreakBond(b->getAtom1(), b->getAtom2());
            geomrefiner->dict.RefiBonds.erase(b);
            UpdateGeom();
            return;
        }
        b++;
    }
}

void DictEditCanvas::OnDeleteBond()
{
    // delete a bond
    OnDeleteBond(pickedBond);
    SetPick((Bond*) NULL);
    UpdateGeom();
}

void DictEditCanvas::OnRemoveHydrogens()
{
    AtomStack->Clear();
    Molecule *mol = (Molecule*)geomrefiner->GetCurrentModel();
    for (MIIter<RESIDUE> res = mol->GetResidues(); res; ++res)
    {
        for (int i = 0; i < res->atomCount(); i++)
        {
            MIAtom *atom = res->atom(i);
            if (MIAtom::MIIsHydrogen(atom))
            {
                AtomStack->Push(atom, res, mol);
            }
        }
    }
    OnRemoveAtoms();
}

void DictEditCanvas::OnRemoveAtoms()
{
    // delete all the atoms on the stack
    while (!AtomStack->empty())
    {
        OnRemoveAtom();
    }
}

//TODO: looks like we need to remove chirals centered on the removed
// atom.  Generally these will go away in UpdateGeom anyhow, but a
// user-added or user-deleted chiral will remain.
void DictEditCanvas::OnRemoveAtom()
{
    if (AtomStack->size() <= 0)
    {
        return;
    }
    // delete the top atom on the stack
    MIAtom *atom = AtomStack->Pop();
    Bond *b;
    ANGLE *a;
    PLANE *p;
    int i;
    //Have to reset loop several times because pointers to vectors change on
    //deletes
    for (i = 0; i < (int)geomrefiner->dict.RefiBonds.size(); i++)
    {
        b = &(geomrefiner->dict.RefiBonds[i]);
        if (atom == b->getAtom1()
            || atom == b->getAtom2())
        {
            OnDeleteBond(b);
            //if a delete happens this is going to be
            SetPick((Bond*) NULL);
            //invalidated for sure
            i--; //Redo the same index since the "array" in one less now. As a
            //result this get the "next" entry from the one we just deleted
        }
    }
    for (i = 0; i < (int)geomrefiner->dict.RefiAngles.size(); i++)
    {
        a = &(geomrefiner->dict.RefiAngles[i]);
        if (strcmp(atom->name(), a->getAtom1()->name()) == 0
            || strcmp(atom->name(), a->getAtom2()->name()) == 0
            || strcmp(atom->name(), a->atom3->name()) == 0)
        {
            geomrefiner->dict.RefiAngles.erase(geomrefiner->dict.RefiAngles.begin() +i);
            SetPick((ANGLE*) NULL);
            i--;
        }
    }
    // Go through the planes
    // This is trickier than the others because we're dealing with a map object.
    // If we delete a plane, reset to the beginning, DO NOT INCREMENT at the end
    // of the loop (otherwise what's the point of resetting to the start) and
    // run through the map of planes again.
    bool allow_incr = true;
    map<int, PLANE*>::iterator pl;
    pl = planes.begin();
    while (pl != planes.end())
    {
        p = pl->second;
        allow_incr = true;
        //Does this plan contain the atom I'm deleteing?
        for (int j = 0; j < p->natoms; j++)
        {
            if (strcmp(atom->name(), p->atoms[j]->name()) == 0)
            {
                SetPick(p);
                if (removeAtomFromPlane(atom, false))
                {
                    SetPick((PLANE*) NULL);
                    pl = planes.begin();
                    allow_incr = false;
                    break;
                }
            }
        }
        if (allow_incr)
        {
            pl++;
        }
    }
    vector<TORSION>::iterator rt = geomrefiner->dict.RefiTorsions.begin();
    vector<TORSION>::iterator e = geomrefiner->dict.RefiTorsions.end();

    while (rt != e)
    {
        TORSION *t = &*rt;
        if (strcmp(atom->name(), t->getAtom1()->name()) == 0
            || strcmp(atom->name(), t->getAtom2()->name()) == 0
            || strcmp(atom->name(), t->atom3->name()) == 0
            || strcmp(atom->name(), t->atom4->name()) == 0)
        {
            geomrefiner->dict.RefiTorsions.erase(rt);
            if (t == pickedTorsion)
            {
                pickedTorsion = NULL;
            }
            rt = geomrefiner->dict.RefiTorsions.begin();
            e = geomrefiner->dict.RefiTorsions.end();
        }
        else
        {
            rt++;
        }
    }
    {
        MIAtomList atoms;
        atoms.push_back(atom);
        AtomStack->atomsToBeDeleted(0, atoms); // mol is not known here
    }
    confs.Purge(atom);
    model->DeleteAtom(atom);
    UpdateGeom();
}

void DictEditCanvas::OnChangeAtomType()
{
    int choice;
    int tmnb; //temp max number bonds
    float anb; //atom total bond order
    const char *oldname;
    MIAtom *atom;
    RESIDUE *res;
    Bond *orig_picked;
    std::string newname;
    MIAtom tatom;
    vector<Bond> edges;
    vector<Bond>::iterator i, e;
    static const char *choices[] =
    {
        "",
        "H", "HE",
        "LI", "BE", "B", "C", "N", "O", "F", "NE",
        "NA", "MG", "AL", "SI", "P", "S", "CL", "AR",
        "K", "CA", "SC", "TI", "V", "CR", "MN", "FE", "CO",
        "NI", "CU", "ZN", "GA", "GE", "AS", "SE", "BR", "KR",
        "RB", "SR", "Y", "ZR", "NB", "MO", "TC", "RU", "RH",
        "PD", "AG", "CD", "IN", "SN", "SB", "TE", "I", "XE",
        "CS", "BA", "LA",
        "CE", "PR", "ND", "PM", "SM", "EU", "GD", "TB", "DY",
        "HO", "ER", "TM", "YB", "LU",
        "HF", "TA", "W", "RE", "OS", "IR",
        "PT", "AU", "HG", "TL", "PB", "BI", "PO", "AT", "RN",
        "FR", "RA", "AC",
        "TH", "PA", "U", "NP", "PU", "AM", "CM", "BK", "CF",
        "ES", "FM", "MD", "NO", "LR",
        "KH"
    };

    static QStringList choice_vec;
    if (!choice_vec.size())
    {
        for (unsigned int i = 0; i < 105; ++i)
        {
            choice_vec += choices[i];
        }
    }
    QString item = QInputDialog::getItem(this, "Change atom type to:", "Change Atom Type", choice_vec);
    choice = choice_vec.indexOf(item);
    if (choice < 1)
    {
        return;
    }
    AtomStack->Pop(atom, res);
    if (atom == NULL)
    {
        return;
    }
    // Get the current atom's total bond order and check that it is <= to what
    // the new atom type can support
    SearchBonds(atom, edges);
    anb = GetBondOrder(edges);
    // Only change the first part of the atom name, keep the numbers and the rest
    oldname = atom->name();
    while (!isdigit(*oldname) && *oldname != '\0')
    {
        oldname++;
    }
    newname = choice_vec.at(choice).toAscii().constData();
    newname += oldname;
    tatom.setName(newname.c_str());
    tatom.setAtomicnumber(choice);
    tmnb = MaxNumBonds(&tatom);
    if (anb > (float) tmnb)
    {
        QMessageBox::warning(this, "Error", "The old atom type has more edges than the\nnew atom type can support. Plese delete\nsome of the edges.");
        return;
    }

    //Passes checks, apply it to the real atom
    atom->setName(newname.c_str());
    atom->setColor(color_by_name(atom->name()));
    atom->setAtomicnumber(choice);
    if (edges.size() == 0)
    {
        printf("Atom has no bonds!\n");
    }
    else
    {
        i = edges.begin();
        e = edges.end();
        orig_picked = pickedBond;
        while (i != e)
        {
            SetPick(&*i);
            OnSetBondLength2(chemlib::IdealBondLength(*i));
            i++;
        }
        SetPick(orig_picked);
    }
    UpdateGeom();
}

void DictEditCanvas::on_showChirals_clicked()
{
    ShowChirals ^= 1;
    UpdateGeom();
}

void DictEditCanvas::on_showTorsions_clicked()
{
    ShowTorsions ^= 1;
    UpdateGeom();
}

void DictEditCanvas::on_showPlanes_clicked()
{
    ShowPlanes ^= 1;
    UpdateGeom();
}

void DictEditCanvas::on_showAtomLabels_clicked()
{
    scene->showAtomLabels = !scene->showAtomLabels;
    UpdateGeom();
}

void DictEditCanvas::on_showHydrogens_clicked()
{
    ShowHydrogens ^= 1;
    scene->renderer->setHideHydrogens(!ShowHydrogens);
    UpdateGeom();
}

void DictEditCanvas::on_showAngles_clicked()
{
    ShowAngles ^= 1;
    UpdateGeom();
}

void DictEditCanvas::on_showBonds_clicked()
{
    ShowBonds ^= 1;
    UpdateGeom();
}

TORSION*DictEditCanvas::SearchTorsions(MIAtom *a1, MIAtom *a2)
{
    vector<TORSION>::iterator i, e;
    i = geomrefiner->dict.RefiTorsions.begin();
    e = geomrefiner->dict.RefiTorsions.end();
    for (; i != e; i++)
    {
        if ( (i->getAtom2() == a1 || i->atom3 == a1)
             && (i->getAtom2() == a2 || i->atom3 == a2) )
        {
            return &*i;
        }
    }

    return NULL;
}

Bond*DictEditCanvas::SearchBonds(MIAtom *a1, MIAtom *a2)
{
    if (bondsdirty)
    {
        std::sort(geomrefiner->dict.RefiBonds.begin(), geomrefiner->dict.RefiBonds.end(), ordering);
        bondsdirty = false;
    }
    vector<Bond>::iterator i = geomrefiner->dict.RefiBonds.begin();
    vector<Bond>::iterator e = geomrefiner->dict.RefiBonds.end();
    while (i != e)
    {
        if ( (i->getAtom1() == a1 && i->getAtom2() == a2) || (i->getAtom1() == a2 && i->getAtom2() == a1) )
        {
            return &*i;
        }
        i++;
    }
    return NULL;
}

Bond*DictEditCanvas::SearchMoleculeBonds(MIAtom *a1, MIAtom *a2)
{
    vector<Bond>::iterator i = model->getBonds().begin();
    vector<Bond>::iterator e = model->getBonds().end();
    while (i != e)
    {
        if ( (i->getAtom1() == a1 && i->getAtom2() == a2) || (i->getAtom1() == a2 && i->getAtom2() == a1) )
        {
            return &*i;
        }
        i++;
    }
    return NULL;
}

void DictEditCanvas::SearchBonds(const MIAtom *a1, vector<Bond> &edges)
{
    vector<Bond>::iterator i = geomrefiner->dict.RefiBonds.begin();
    vector<Bond>::iterator e = geomrefiner->dict.RefiBonds.end();
    while (i != e)
    {
        if (i->getAtom1() == a1 || i->getAtom2() == a1)
        {
            edges.push_back(*i);
        }
        i++;
    }
}

PLANE*DictEditCanvas::SearchPlane(MIAtom *a)
{
    map<int, PLANE*>::iterator i, e;
    int j;
    i = planes.begin();
    e = planes.end();
    while (i != e)
    {
        for (j = 0; j < i->second->natoms; j++)
        {
            if (a == i->second->atoms[j])
            {
                return (i->second);
            }
        }
        i++;
    }
    return NULL;
}

void DictEditCanvas::SetNeedsRefine(bool n)
{
    std::string output;
    QFont font = form->optimizeButton->font();
    needsrefining = n;
    if (needsrefining)
    {
        output = ::format("*Optimize*");
        font.setBold(true);
    }
    else
    {
        output = ::format("Optimize");
        font.setBold(false);
    }

    form->optimizeButton->setText(output.c_str());
    form->optimizeButton->setFont(font);
}

void DictEditCanvas::on_optimizeButton_clicked()
{
    OnRefine(form->optimizeIterations->value());
}

void DictEditCanvas::OnRefine(int iterations)
{
    if (confs.NumberSets()-1 > 1)
    {
        if (QMessageBox::question(this, "Remove conformers?",
                                  "Optimize can only be performed on a single conformer.\nRemove conformers and continue?",
                                  QMessageBox::Yes | QMessageBox::No, QMessageBox::No) != QMessageBox::Yes)
        {
            // parent->Raise();
            return;
        }
    }
    //Set all modified angles to false, since this is going to take it into
    //account and we don't want them to be skipped if the user changes the
    //bondlength after doing a refine
    clearModifiedAngles();
    confs.Clear();
    for (int i = 0; i < iterations; i++)
    {
        geomrefiner->Refine();
    }
    map<int, PLANE*>::iterator pln, e;
    pln = planes.begin();
    e = planes.end();
    for (; pln != e; pln++)
    {
        lsqplane(*pln->second);
    }
    confs.Save(CurrentAtoms, model);
    on_conformerSpinbox_valueChanged(1);

    UpdateGeom();
    ReDraw();
    //After refining, we're not going to need it again.
    SetNeedsRefine(false);
}

void DictEditCanvas::GrowPlane(int size)
{
    MIAtom **newatoms;
    if (size < 1)
    {
        return;
    }
    if ((newatoms = (MIAtom**) malloc((pickedPlane->natoms + size) * sizeof(MIAtom*))) == NULL)
    {
        QMessageBox::warning(this, "Error", "Memory allocation failed! Aborting add atom to plane");
        return;
    }

    memcpy(newatoms, pickedPlane->atoms, pickedPlane->natoms * sizeof(MIAtom*));
    pickedPlane->natoms += size;
    pickedPlane->atoms = newatoms;
}

void DictEditCanvas::OnIncludeAtomInPlane()
{
    MIAtom *a = NULL;
    RESIDUE *res = NULL;
    if (pickedPlane == NULL)
    {
        return;
    }
    // What atom are we adding
    AtomStack->Pop(a, res);
    if (a == NULL)
    {
        return;
    }
    //sanity checks
    if (pickedPlane == NULL)
    {
        return;
    }
    if (AtomInPlane(a))
    {
        ReDraw();
        return;
    }
    GrowPlane(1);
    pickedPlane->atoms[pickedPlane->natoms - 1] = a;
    lsqplane(*pickedPlane);
    UpdateGeom();
}

void DictEditCanvas::OnIncludeAtomsInPlane()
{
    if (pickedPlane == NULL)
    {
        return;
    }
    MIAtom *a = NULL;
    RESIDUE *res = NULL;
    vector<MIAtom*> atoms;
    while (!AtomStack->empty() )
    {
        AtomStack->Pop(a, res);
        if (a == NULL)
        {
            QMessageBox::warning(this, "No atoms on stack", "An atom must be on the stack.");
            return;
        }
        if (AtomInPlane(a))
        {
            continue;
        }
        atoms.push_back(a);
    }
    if (atoms.size() < 1)
    {
        ReDraw();
        return;
    }
    GrowPlane(atoms.size());
    vector<MIAtom*>::iterator i, e = atoms.end();
    int idx = atoms.size();
    for (i = atoms.begin(); i != e; i++)
    {
        pickedPlane->atoms[pickedPlane->natoms - idx] = *i;
        idx--; //natoms is the total number, not the 0 index
    }
    lsqplane(*pickedPlane);
    UpdateGeom();
}

// reset to starting conditions
void DictEditCanvas::on_resetButton_clicked()
{
    model = (Molecule*)geomrefiner->GetCurrentModel();
    RESIDUE *dictres = geomrefiner->dict.GetDictResidue(model->getResidues()->type().c_str(), 0);
    // Effectivly the cancel function from Geomrefiner
    geomrefiner->dict.Clear();
    delete model;
    /*geomrefiner->nRefiRes = 0;
       geomrefiner->RefiRes = NULL;
       geomrefiner->GetCurrentModel() = NULL;*/
    // end of the cancel function

    // blank out the atom stack so that we won't have bad pointers hanging around
    AtomStack->Clear();

    // similar to the geomrefiner OnEdit function and our own init function
    RESIDUE *reslist = new RESIDUE(*dictres);
    Molecule *newmodel = new Molecule(reslist, "Dictionary", NULL, NULL, 0, MoleculeType::Other);
    geomrefiner->SetRefiRes(reslist, reslist, newmodel);
    model = (Molecule*)geomrefiner->GetCurrentModel();
    clearModifiedAngles();

    // Get planes that already exist so that we can display them
    createPlanes();
    SetNeedsRefine(false); //After a reset we don't need a refine

    // Force redo of setup
    makeCurrent();
    initializeGL();

    //Reset conformations
    confs.Clear();
    GetConfs(confs, geomrefiner->GetCurrentModel()->getResidues(), &geomrefiner->dict, geomrefiner->GetCurrentModel());
    on_conformerSpinbox_valueChanged(1);

    UpdateGeom();
}

//Updated to handle multiple confs KWB 11/8/2005

void DictEditCanvas::OnOk()
{
    geomrefiner->unlockRefineTarget();

    //This needs to push all the results in the the actual dictionary
    const char *restype = model->getResidues()->type().c_str();

    // Adding the residue in append mode will blast out the other entires for it
    // in the other dictionaries (i.e. the PlaneDict)
    model->getResidues()->setPrefBonds(geomrefiner->dict.RefiBonds);
    model->getResidues()->setPrefAngles(geomrefiner->dict.RefiAngles);
    RESIDUE *multiConf = ExpandConfs(model->getResidues(), confs);
    geomrefiner->dict.LoadRes(multiConf, true, true, false);
    FreeResidueList(multiConf);
    model->getResidues()->clearPrefBonds();
    model->getResidues()->clearPrefAngles();

    //Add all the planes back in
    if (planes.size() > 0)
    {
        map<int, PLANE*>::iterator i, e;
        i = planes.begin();
        e = planes.end();
        for (; i != e; i++)
        {
            if (geomrefiner->dict.AddPlane(*i->second) != 1)
            {
                Logger::log("Plane was NOT saved!!\n");
            }
        }
    }

    if (geomrefiner->dict.RefiTorsions.size() > 0)
    {
        vector<TORSION>::iterator i, e;
        i = geomrefiner->dict.RefiTorsions.begin();
        e = geomrefiner->dict.RefiTorsions.end();
        for (; i != e; i++)
        {
            if (geomrefiner->dict.AddTorsion(*i) != 1)
            {
                Logger::log("Torsion was NOT saved!!\n");
            }
        }
    }

    if (geomrefiner->dict.RefiChirals.size() > 0)
    {
        UpdateChirals();
        vector<CHIRAL>::iterator i, e;
        i = geomrefiner->dict.RefiChirals.begin();
        e = geomrefiner->dict.RefiChirals.end();
        for (; i != e; i++)
        {
            if (i->flags != CHIRAL_DELETED)
            {
                if (geomrefiner->dict.AddChiral(*i, restype) != 1)
                {
                    Logger::log("Chiral was NOT saved!!\n");
                }
                i->center->chiral_class(CH_NONE);
            }
        }
    }
}

void DictEditCanvas::OnRenameAtom()
{
    const char *on, *nn;  //oldname ptr and new name ptr
    RESIDUE *res;
    MIAtom *atom;
    if (AtomStack->empty())
    {
        return;
    }
    AtomStack->Pop(atom, res);
    if (atom == NULL)
    {
        return;
    }
    // User CANNOT edit the full name of the atom since it determines type
    // Therefore the user may only edit the "number" of the atom i.e. they can
    // rename C65 to C14 and that's okay

    QString newname = QInputDialog::getText(this, "New Atom Name", "Enter new atom name:", QLineEdit::Normal, atom->name());
    if (newname.isEmpty())
    {
        return;
    }

    //Check size < chemlib::MAXATOMNAME
    if (newname.size() > chemlib::MAXATOMNAME)
    {
        QMessageBox::warning(this, "Error", "Atom name too long!");
        return;
    }

    //Check if user "changed" atom type
    std::string oldname = Atomic_Name(atom->atomicnumber());
    MIStringTrim(oldname, false);
    MIStringTrim(oldname); //Trim from both sides
    on = oldname.c_str();
    nn = newname.toAscii().constData();
    while (*on != '\0' && *nn != '\0')
    {
        if (*on == *nn)
        {
            on++;
            nn++;
        }
        else if (isdigit(*on))
        {
            break; //Atom type not changed, break out of loop
        }
        else
        {
            QMessageBox::warning(this, "Error", "Error: you tried to change the atom type!");
            return;
            break;
        }
    }

    atom->setName(newname.toAscii().constData());
    UpdateGeom();
}

void DictEditCanvas::OnAddChiral()
{
    CHIRAL c;
    MIAtom *a = NULL;
    a = AtomStack->Pop();
    if (!a)
    {
        return;
    }
    a->chiral_class(CH_TETRAHEDRAL);
    if (pickedChiral != NULL) //There might be a deleted chiral there...
    {
        if (pickedChiral->flags != CHIRAL_DELETED)
        {
            return;
        }
        pickedChiral->flags ^= CHIRAL_DELETED;
    }
    else
    {
        c.center = a;
        c.flags = CHIRAL_USERADD;
        geomrefiner->dict.RefiChirals.push_back(c);
    }
    SetPick((CHIRAL*) NULL);
    UpdateGeom();
}

void DictEditCanvas::OnRemoveChiral()
{
    if (pickedChiral == NULL)
    {
        return;
    }
    if (pickedChiral->flags == CHIRAL_USERADD)
    {
        //Outright delete it
        vector<CHIRAL>::iterator i, e;
        i = geomrefiner->dict.RefiChirals.begin();
        e = geomrefiner->dict.RefiChirals.end();
        for (; i != e; i++)
        {
            if (&*i == pickedChiral)
            {
                geomrefiner->dict.RefiChirals.erase(i);
                break;
            }
        }
    }
    else if (pickedChiral->flags == CHIRAL_AUTO)
    {
        pickedChiral->flags = CHIRAL_DELETED;
        pickedChiral->center->chiral_class(CH_NONE);
    }

    SetPick((CHIRAL*) NULL);

    UpdateGeom();
}

CHIRAL*DictEditCanvas::SearchChirals(MIAtom *atom)
{
    vector<CHIRAL>::iterator i, e;
    //Are there any chirals?
    if (geomrefiner->dict.RefiChirals.size() == 0)
    {
        return NULL;
    }
    i = geomrefiner->dict.RefiChirals.begin();
    e = geomrefiner->dict.RefiChirals.end();
    for (; i != e; i++)
    {
        if (i->center == atom)
        {
            return &*i;
        }
    }
    return NULL;
}

void DictEditCanvas::UpdateChirals()
{
    if (geomrefiner->dict.RefiChirals.size() == 0)
    {
        return;
    }
    vector<CHIRAL>::iterator i, e;
    vector<Bond>::iterator k, d;
    vector<Bond> temp;
    MIAtom *ary[3];
    int ctr = 0;
    i = geomrefiner->dict.RefiChirals.begin();
    e = geomrefiner->dict.RefiChirals.end();
    for (; i != e; i++)
    {
        temp.clear();
        SearchBonds(i->center, temp);
        if (temp.size() != 3) //TODO: Replace this with a smarter test for chiral centers
        {
            Logger::log("Invalid Chiral center found\n");
            return;
        }
        k = temp.begin();
        d = temp.end();
        for (; k != d; k++)
        {
            if (ctr > 2)
            {
                break;
            }
            if (k->getAtom1() != i->center)
            {
                ary[ctr++] = k->getAtom1();
                continue;
            }
            else if (k->getAtom2() != i->center)
            {
                ary[ctr++] = k->getAtom2();
                continue;
            }
            else
            {
                Logger::log("Couldn't figure out who belongs to who! %d: 0x%x 0x%x 0x%x", ctr, i->center, k->getAtom1(), k->getAtom2());
            }
        }
        i->setAtom1(ary[0]);
        i->setAtom2(ary[1]);
        i->setAtom3(ary[2]);
    }
}

void DictEditCanvas::OnNextConformer()
{
    on_conformerSpinbox_valueChanged(conf+1);
}

void DictEditCanvas::OnPrevConformer()
{
    on_conformerSpinbox_valueChanged(conf-1);
}

void DictEditCanvas::on_conformerSpinbox_valueChanged(int conformer)
{
    if (conformer > 0 && conformer < confs.NumberSets())
    {
        conf = conformer;
        confs.Restore(conf);
        UpdateGeom();
        UpdateButtons();
    }
}

void DictEditCanvas::OnGenerateConformers()
{
    if (QMessageBox::question(this, "Replace or Append?",
                              "Replace existing conformer(s)?",
                              QMessageBox::Yes | QMessageBox::No, QMessageBox::No) == QMessageBox::No)
    {
        generateConformers(false);
    }
    else
    {
        generateConformers(true);
    }
}


void DictEditCanvas::generateConformers(bool replace)
{
    if (replace)
    {
        confs.Clear();
    }

    int nConf = conflib::GenerateEnsemble(model->getResidues(),
                                          model,
                                          geomrefiner->dict.RefiBonds,
                                          geomrefiner->dict.RefiTorsions,
                                          confs);
    on_conformerSpinbox_valueChanged(1);

    std::string message = ::format("Generated %d conformers", nConf);
    QMessageBox::warning(this, message.c_str(), message.c_str());
}


void DictEditCanvas::OnPopupMenu(const MIActionEvent &event)
{
    int id = event.GetId();
    if (id == ID_POPUP_LINES)
    {
        scene->renderer->setRenderStyle(RenderStyle::getDefaultLine());
    }
    else if (id == ID_POPUP_BALL_LINES)
    {
        scene->renderer->setRenderStyle(RenderStyle::getDefaultBallAndLine());
    }
    else if (id == ID_POPUP_STICKS)
    {
        scene->renderer->setRenderStyle(RenderStyle::getDefaultStick());
    }
    else if (id == ID_POPUP_BALL_STICKS)
    {
        scene->renderer->setRenderStyle(RenderStyle::getDefaultBallAndStick());
    }
    reRender();
}

void DictEditCanvas::clearModifiedAngles()
{
    if (geomrefiner == NULL)
    {
        modifiedAngles.clear();
    }
    else
    {
        modifiedAngles.resize(geomrefiner->dict.RefiAngles.size());
        fill(modifiedAngles.begin(), modifiedAngles.end(), false);
    }
}

void DictEditCanvas::createPlanes()
{
    planes.clear();
    // Get planes that already exist so that we can display them
    nextAnnotationId = 1; //Id must start at 1 to prevent off by one error
    vector<PLANE>::iterator iter = geomrefiner->dict.RefiPlanes.begin();
    vector<PLANE>::iterator end = geomrefiner->dict.RefiPlanes.end();
    while (iter != end)
    {
        PLANE &plane = *iter;
        ++iter;
        lsqplane(plane);
        planes.insert(pair<int, PLANE*>(nextAnnotationId, &plane));
        nextAnnotationId++;
    }
}

void DictEditCanvas::reRender()
{
    updateGL();
}

int DictEditCanvas::getEventX(QMouseEvent &e)
{
    int x = e.x();
    if (stereoView->isStereo() && !stereoView->isHardwareStereo())
    {
        if (x > width() / 2)
        {
            x -= width() / 2;
        }
    }
    return x;
}

int DictEditCanvas::getEventY(QMouseEvent &e)
{
    return e.y();
}

float DictEditCanvas::GetBondOrder(Bond *e)
{
    return GetBondOrder(e->getOrder());
}

void DictEditCanvas::SetPick(Bond *b)
{
    pickedBond = b;
    scene->pickedBond = pickedBond;
}

void DictEditCanvas::SetPick(ANGLE *a)
{
    pickedAngle = a;
    scene->pickedAngle = pickedAngle;
}

void DictEditCanvas::SetPick(PLANE *p)
{
    pickedPlane = p;
    scene->pickedPlane = pickedPlane;
}

void DictEditCanvas::SetPick(TORSION *t)
{
    pickedTorsion = t;
}

void DictEditCanvas::SetPick(CHIRAL *c)
{
    if (c == pickedChiral)
    {
        pickedChiral = NULL;
    }
    else
    {
        pickedChiral = c;
    }
}

chemlib::MIMoleculeBase*DictEditCanvas::getMolecule()
{
    return geomrefiner->GetCurrentModel();
}
