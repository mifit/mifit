#ifndef MIFIT_UI_DICTEDITCANVAS_H_
#define MIFIT_UI_DICTEDITCANVAS_H_

#include <map>

#include <chemlib/chemlib.h>
#include "core/corelib.h"

class DictEditAnglePickingRenderable;
class DictEditAnnotationPickingRenderable;
class DictEditAtomPickingRenderable;
class DictEditBondPickingRenderable;
class Stack;
class Displaylist;

namespace mi
{
    namespace opengl
    {
        class StereoView;
        class Camera;
        class Frustum;

        namespace interact
        {
            class FieldOfViewZoomCommand;
            class MouseArcBallOrbitor;
            class MousePicker;
            class MouseTranslator;
            class MouseZoomer;
        }
    }
}

class DictEditScene;
class GeomRefiner;

class QKeyEvent;
class QMenuBar;

#include <QGLWidget>
#include "MIMenu.h"
#include "MIEventHandler.h"

namespace Ui
{
    class DictEditorForm;
}

/**
 * Canvas for editing a dictionary entry.
 */
class DictEditDialog;
class DictEditCanvas : public QGLWidget, public MIEventHandler
{
    Q_OBJECT

    mi::opengl::StereoView *stereoView;
    mi::opengl::Camera *camera;
    mi::opengl::Frustum *frustum;
    DictEditScene *scene;
    mi::opengl::interact::MouseArcBallOrbitor *cameraMouseOrbitor;
    mi::opengl::interact::MouseTranslator *cameraMouseTranslator;
    mi::opengl::interact::FieldOfViewZoomCommand *zoomCommand;
    mi::opengl::interact::MouseZoomer *mouseZoomer;
    mi::opengl::interact::MousePicker *mousePicker;
    DictEditAtomPickingRenderable *atomPickingRenderable;
    DictEditBondPickingRenderable *bondPickingRenderable;
    DictEditAnglePickingRenderable *anglePickingRenderable;
    DictEditAnnotationPickingRenderable *annotationPickingRenderable;

    QMenuBar *menuBar;

    QAction *setAngleAction;
    QAction *changeAtomTypeAction;
    QAction *renameAtomAction;
    QAction *removeAtomAction;
    QAction *removeAtomsAction;
    QAction *removeHydrogensAction;

    QAction *setBondLengthAction;
    QAction *changeBondOrderAction;
    QAction *removeBondAction;

    QAction *invertChiralAction;
    QAction *addChiralAction;
    QAction *removeChiralAction;

    QAction *nextConformerAction;
    QAction *prevConformerAction;

    QAction *addPlaneAction;
    QAction *addPlaneAtomAction;
    QAction *addPlaneAtomsAction;
    QAction *removePlaneAtomAction;
    QAction *removePlaneAtomsAction;
    QAction *removePlaneAction;

    QMenu *popupMenu;

    /**
     * Number of milliseconds between down and up mouse events to be considered a click.
     */
    static const int mouseClickThreshold;

    /**
     * Track if the bonds need to be sorted before doing a search
     */
    bool bondsdirty;

    /**
     * Track for if the model needs to be refined
     */
    bool needsrefining;

    std::map<int, chemlib::PLANE*> planes;

    /**
     * Array of bools for if an angle has been modified
     */
    std::vector<bool> modifiedAngles;

    /**
     * The next available annotation id
     */
    int nextAnnotationId;

    /**
     * The token (index) of the currently displayed conformer
     */
    int conf;

    /**
     * Stack of atoms used by the picking operations. Created when view is created.
     */
    Stack *AtomStack;

    Displaylist *Models;
    QPoint mouseStart;
    GeomRefiner *geomrefiner;
    Molecule *model;
    chemlib::GeomSaver confs;
    void DrawPlanes();
    void DrawTorsions();
    void DrawChirals();
    void UpdateChirals();
    void UpdateGeom();
    void AutoFindChirals(bool copyChiralClass = false);
    bool ShowAngles;
    bool ShowHydrogens;
    bool ShowBonds;
    bool ShowPlanes;
    bool ShowTorsions;
    bool ShowChirals;
    bool working;
    bool initialized;

    QTimer *mouseTimer;
    bool isRotating;
    bool isZooming;
    bool isPanning;

public slots:
    void OnTimer();

public:
    DictEditDialog *parent;

    std::vector<chemlib::MIAtom*> CurrentAtoms;
    chemlib::Bond *pickedBond;
    chemlib::ANGLE *pickedAngle;
    chemlib::PLANE *pickedPlane;
    chemlib::TORSION *pickedTorsion;
    chemlib::CHIRAL *pickedChiral;
    bool AtomInPlane(chemlib::MIAtom *a);
    float GetBondOrder(unsigned char);
    unsigned char GetBondOrder(float);
    float GetBondOrder(std::vector<chemlib::Bond>&);
    float GetBondOrder(chemlib::Bond *e);
    void SetPick(chemlib::Bond *b);
    void SetPick(chemlib::ANGLE *a);
    void SetPick(chemlib::PLANE *p);
    void SetPick(chemlib::TORSION *t);
    void SetPick(chemlib::CHIRAL *c);
    chemlib::Bond *SearchBonds(chemlib::MIAtom *a1, chemlib::MIAtom *a2);
    chemlib::Bond *SearchMoleculeBonds(chemlib::MIAtom *a1, chemlib::MIAtom *a2);
    chemlib::PLANE *SearchPlane(chemlib::MIAtom *a);
    chemlib::TORSION *SearchTorsions(chemlib::MIAtom *a1, chemlib::MIAtom *a2);
    chemlib::CHIRAL *SearchChirals(chemlib::MIAtom *atom);
    void SearchBonds(const chemlib::MIAtom *a1, std::vector<chemlib::Bond> &edges);
    void DeletePlane(chemlib::PLANE *p);
    void SetNeedsRefine(bool n = true);
    bool handlePick(const QPoint &point);

    void beginDragAction();
    void endDragAction();

    void clearModifiedAngles();

    void initializeForRender();
    void render();

    void updateViewDependentSettings();

    void createPlanes();

    void reRender();

    int getEventX(QMouseEvent &e);
    int getEventY(QMouseEvent &e);

public:

    DictEditCanvas(DictEditDialog *parent);
    virtual ~DictEditCanvas();

    void createMenus();

    void initializeGL();
    void paintGL();
    void resizeGL(int w, int h);
    void keyPressEvent(QKeyEvent *e);
    void mousePressEvent(QMouseEvent *e);
    void mouseReleaseEvent(QMouseEvent *e);
    void mouseMoveEvent(QMouseEvent *e);
    void wheelEvent(QWheelEvent *e);

    void OnOk();
    void ReDraw();


private:
    void OnSetBondLength2(float);
    void OnExport(const char *optionalFilename = NULL);
    bool removeAtomFromPlane(chemlib::MIAtom *a, bool atomNotToBeDeleted);
    // grey in and out buttons based on stack, etc.
    void UpdateButtons();
    void OnDeleteBond(chemlib::Bond*);

    chemlib::MIMoleculeBase *getMolecule();
    void generateConformers(bool replace = false);
    void GrowPlane(int);
    void OnRefine(int);

private slots:
    void OnPopupMenu(QAction *action);

    // button events
    void OnGenerateConformers();
    void OnRemoveAtomFromPlane();
    void OnRemoveAtomsFromPlane();
    void OnRemovePlane();
    void OnAddPlane();
    void OnAddChiral();
    void OnRemoveChiral();
    void OnInvertCenter();
    void OnIncludeAtomInPlane();
    void OnIncludeAtomsInPlane();
    void OnSetAngle();
    void OnSetBondLength();
    void OnChangeBondOrder();
    void OnDeleteBond();
    void OnRemoveHydrogens();
    void OnRemoveAtoms();
    void OnRemoveAtom();
    void OnChangeAtomType();
    void on_showTorsions_clicked();
    void on_showChirals_clicked();
    void on_showPlanes_clicked();
    void on_showAngles_clicked();
    void on_showAtomLabels_clicked();
    void on_showHydrogens_clicked();
    void on_showBonds_clicked();
    void OnRenameAtom();
    void OnNextConformer();
    void OnPrevConformer();

    void on_exportButton_clicked();
    void on_changeNameButton_clicked();
    void on_optimizeButton_clicked();
    void on_resetButton_clicked();
    void on_conformerSpinbox_valueChanged(int conformerNumber);

private:
    Ui::DictEditorForm *form;
};


#endif // MIFIT_UI_DICTEDITCANVAS_H_
