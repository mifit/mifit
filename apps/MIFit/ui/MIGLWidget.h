#ifndef mifit_ui_MIGLWidget_h
#define mifit_ui_MIGLWidget_h

#include <cmath>
#include <QDeclarativeView>
#include <QMdiSubWindow>
#include <QTime>
#include <QVBoxLayout>

#include <math/Quaternion.h>
#include "core/corelib.h"
#include <chemlib/chemlib.h>
#include "macafxwin.h"
#include "SaveModel.h"

#include "RamaPlot.h"

#ifdef _WIN32
#include <time.h>
#endif


class EMap;
class Displaylist;
class MIGLWidget;

class CMolwViewAnnotationPickingRenderable;
class CMolwViewAtomPickingRenderable;
class CMolwViewBondPickingRenderable;
class CMolwViewSlabPickingRenderable;
class CMolwViewScene;

class GLRenderer;
class QActionGroup;
class QDeclarativeEngine;
class QGraphicsObject;
class QGraphicsSceneMouseEvent;
class QResizeEvent;
class InterpBox;
class MainItem;
class NavigatorCanvas;
namespace mi
{
    namespace opengl
    {
        class Camera;
        class Frustum;
        class Scene;
        class StereoView;
        class Viewport;
        namespace interact
        {
            class MousePicker;
            class TargetFeedback;
        }
    }
}

class MAP_POINT;


#define sign(a) (a < 0 ? -1 : 1)


class MIGLWidget : public QDeclarativeView
{
    Q_OBJECT

    //Display list contains the list of objects owned by this document
    Displaylist *Models;
    int newfile;

    typedef std::vector<CMapHeaderBase> MapHeaderList;
    /**
     * These map headers are read from XML session files, but since
     * the maps are not loaded when the XML is parsed, they must be
     * stored for later use when the maps are loaded.
     */
    MapHeaderList inputMapHeaders;


    std::string _filename;
    std::string _title;
    bool _modified;

    mi::opengl::StereoView *stereoView;
    CMolwViewScene *scene;
    mi::opengl::Frustum *frustum;
    mi::opengl::Camera *camera;
    mi::opengl::interact::TargetFeedback *targetFeedback;
    mi::opengl::interact::MousePicker *mousePicker;

    CMolwViewAnnotationPickingRenderable *annotationPickingRenderable;
    CMolwViewAtomPickingRenderable *atomPickingRenderable;
    CMolwViewBondPickingRenderable *bondPickingRenderable;
    CMolwViewSlabPickingRenderable *slabPickingRenderable;

    bool is_drawing;
    chemlib::Residue *PentamerStart;
    int batonposition;
    chemlib::MIAtom *batonatom;

    int m_timer;
    int ClockTick;
    int TimerOn;

    QTime time;

public:
    void timerEvent(QTimerEvent *e);
private:
    QWidget *fullScreenWidget;
    QVBoxLayout *fullScreenLayout;
    QMdiSubWindow *oldParent;
    QWidget *containerWidget;
    MainItem *mainItem;
    QGraphicsObject *stackItem;

    chemlib::Residue *focusres;
    bool focusresDeleted;
    bool DragStart;

    chemlib::GeomSaver savefits;
    std::vector<chemlib::MIAtom*> CurrentAtoms;
    std::vector<MAP_POINT> PList;
    int SurfaceCurrent;
    unsigned int FitToken;
    int Roll;
    int Rock;
    float RockRange;
    float RockAngle;
    float RockIncrement;
    float RollIncrement;
    bool TopView;
    bool MouseCaptured;
    bool Blink;
    int BlinkCounter;
    int BlinkTime;
    long MouseStillTime;
    long ToolTipInterval;
    bool MouseInWindow;
    bool AutoFit;

    CPoint mouse, mousestart;

    QMenu *popup_menu;
    QAction *fitResidueAction;
    QAction *rotateAction;
    QAction *translateAction;
    QAction *setTorsionAction;
    QAction *acceptFitAction;
    QAction *cancelFitAction;
    QAction *replaceAndFitAction;
    QAction *deleteResidueAction;
    QAction *deleteAtomAction;
    QAction *refineRegionAction;
    QAction *acceptRefineAction;
    QAction *cancelRefineAction;
    QAction *saveAction;
    QAction *fullscreenAction;
    QAction *clearLabelsAction;

    QAction *dotSurfVdwAction;
    QAction *dotSurfSolventExposedAction;
    QAction *dotSurfAtomSphereAction;
    QAction *dotSurfResidueAction;
    QAction *dotSurfResiduesAction;
    QAction *dotSurfAtomAction;
    QAction *dotSurfAtomsAction;
    QAction *dotSurfClearAction;

    QMenu *solidsurf_popup_menu;
    QAction *solidsurfActions[10];
    QActionGroup *solidsurfActionGroup;
    QActionGroup *solidsurfCommonActionGroup;

    QMenu *newSolidsurfMenu(bool include_selection_items = true);

    unsigned int solidsurf_current_surface;
    int SelectType;

    NavigatorCanvas *navigator;
    void createNavigatorWindow();
    unsigned int DoingRange;

    unsigned int ClusterSize;


    void prepareAtomPicking();
    void prepareBondPicking();

    void doSlabDrag(int x, int y, int dx, int dy);

    void findResidueAndMoleculeForAtom(chemlib::MIAtom *atom, chemlib::Residue* &res, Molecule* &mol);
    void doSetFocusResidue(chemlib::Residue *res);

    void StackExpandTop2Range();

    /**
     * Get rid of any pointers or objects pointing to model.
     * Prevents accessing bad pointers leading to a seg fault.
     */
    void Purge(Molecule *model);
    void PurgeChain(chemlib::Residue *chain);
    void Purge(chemlib::Residue *res);
    void Purge(chemlib::MIAtom *atom);

    /**
     * List of checkpoint saved models.
     */
    std::vector<SaveModel> SaveModels;

public:
    MIGLWidget();
    virtual ~MIGLWidget();

    void SetFilename(const std::string &fname)
    {
        _filename = fname;
    }
    std::string GetFilename() const
    {
        return _filename;
    }

    void SetTitle(const std::string &title);
    std::string GetTitle() const
    {
        return _title;
    }

    void SetDocumentSaved(bool saved)
    {
        _modified = !saved;
    }


    bool LoadScript(const char *path);
    virtual bool Close();
    virtual bool SaveAs();
    virtual bool IsModified();
    virtual void Modify(bool modify);

    //****************************************
    //		Miscellaneous functions
    //****************************************
    //@{
    // Create an EMap and load a density map into it.
    //@}
    bool LoadMap(const char *command, const char *path, const std::string &columnNames);
    void exportCurrentModel();

    //@{
    // Save the document into a session (.mlw file).
    //@}
    bool SaveDocument(const std::string &pathname);
    //@{
    // Returns true if the file is an XML .mlw file.
    // Otherwise it is assumed to
    // be the older Molw .mlw file format.
    //@}
    bool IsXMLDocument(const char *pathname);
    //@{
    // Returns true if the fie is XML and contains a <Molecule> tag.
    //@}
    bool IsXMLMolecule(io &file);
    //@{
    // Load an Molecule from an XML .mlw file.
    //@}
    bool LoadXMLMolecule(const char *pathname, int nmodel);
    //@{
    // Purge all references to an EMap.
    // To be executed before deleting an EMap.
    //@}
    void Purge(EMap*);

    //@{
    // Return a pointer to the display list.
    //@}
    Displaylist *GetDisplaylist()
    {
        return Models;
    }

    //@{
    // Return true if this a newly loaded file.
    //@}
    int NewFile()
    {
        return (newfile);
    }

    //@{
    // Clear the newfile flag.
    //@}
    void ResetNewFile()
    {
        newfile = 0;
    }

    //@{
    // Set the newfile flag on.
    //@}
    void SetNewFile()
    {
        newfile = 1;
    }

    //@{
    //
    //@}
    //	void OnSelchangeCombo1();
    //@{
    // Load a model from a PDB format file.
    //@}
    bool LoadPDBFile(const char *pathname);
    //@{
    // Load a dictionary residue.
    //@}
    chemlib::Residue *GetDictRes(const char *key, int confomer = 0);
    //@{
    // Load a dictionary residue.
    //@}
    chemlib::Residue *GetDictRes(const char single, int confomer = 0);
    //@{
    // Return a pointer to the dictionary residue linked list.
    //@}
    std::vector<std::string> GetDictResList();

    //@{
    // Save the document and view data to an XML .mlw file.
    //@}
    void Serialize(XMLArchive &ar); // overridden for document i/o
    //@{
    // Run on the creation of a new document.
    //@}
    virtual bool    OnNewDocument();
    //@{
    // Menu callback for File/Open.
    //@}
    virtual bool    OnOpenDocument(const std::string &pathname);
    //@{
    // Menu callback for File/Save.
    //@}

    virtual bool    OnSaveDocument(const std::string &pathname);

    // Prompts for saving if about to close a modified document. Returns true
    // if ok to close the document (may have saved in the meantime, or set
    // modified to false)
    virtual bool OnSaveModified();

    //@{
    // The current annotation.
    //@}
    Annotation *CurrentAnnotation;


    Molecule *PentamerFrom;
    Molecule *PentamerModel;

    /**
     * Full screen flag.
     */
    bool IsFullScreen;
    /**
     * Cluster residue being optimized
     */
    InterpBox *BoundingBox;
    /**
     * string containing the message printed at the bottom of the screen.
     */
    std::string message;
    /**
     * The viewpoint containing the current screen center, rotation and so forth.
     * Created when view is created.
     */
    ViewPoint *viewpoint;
    /**
     * Stack of atoms used by the picking operations.
     * Created when view is created.
     */
    Stack *AtomStack;
    /**
     * if true then the colors information is saved for an undo operation
     */
    bool SaveColors;
    /**
     * The current residue.
     * If NULL there is currently no focus.
     */
    const chemlib::Residue *getFocusResidue();
    /**
     * Set the current residue.
     */
    void setFocusResidue(chemlib::Residue *res, bool recenter = true, bool deleted = false);

    /**
     * The current model containing the focus residue.
     */
    Molecule *focusnode;
    /**
     * The residue being currently fit, if any.
     * NULL value indicate no residue being fit.
     */
    chemlib::Residue *fitres;
    /**
     * The model being currently fit, if any.
     * NULL value indicate no model being fit.
     */
    Molecule *fitmol;
    /**
     * An atom representing the current center of rotation for fitting operations.
     */
    chemlib::MIAtom fitcenter;
    /**
     * pointer to the atom from which the fitcenter was copied.
     */
    chemlib::MIAtom *fromcenter;
    /**
     * if true then we are torsioning a bond.
     */
    int Torsioning;
    /**
     * if true then the right mouse is used for translation.
     */
    int RightMouseMode;
    /**
     * flag to indicate the canvas should be updated.
     */
    int Update; // picture needs updating
    /**
     * flag to indicate the size is being changed.
     */
    int NewSize;
    /**
     * if true then show the gnomon.
     */
    bool ShowGnomon;
    bool showUnitCell;
    /**
     * if true then show the labels.
     */
    bool ShowLabels;
    /**
     * if true then show the contacts.
     */
    bool ShowContacts;
    /**
     * used to differntiate a single-click of the mouse from a double click.
     */
    bool doubleclicked;

    /**
     * true is we are fitting.
     */
    bool IsFitting();
    /**
     * true is we are torsioning around a bond.
     */
    bool IsTorsioning();
    /**
     * an arrow representing the bond being torsioned for drawing on screen
     * to let the user know which bond and which direction the torsion is.
     */
    std::vector<PLINE> TorsionArrow;
    /**
     * initialize a fitting operation.
     * Fills current atoms and some other things.
     */
    unsigned int SetupFit(const std::vector<chemlib::MIAtom*> &atoms, Molecule *model, char altloc = ' ');
    unsigned int SetupFit(chemlib::Residue *start, chemlib::Residue *end, Molecule *model, char altloc = ' ');

    /**
     * True if it is OK to do a fitting.
     * There are instances when setting up a fit can be bad, such as when
     * there is already a fitting going on.
     */
    bool FitOK();

    /**
     * Used by coloring function.
     * Saves the color to be used.
     */
    int WhenShownColor;
    /**
     * Used by coloring function. Saves the method to be used.
     */
    int WhenShownColorMethod;
    /**
     * Used by radius function.
     */
    int WhenShownRadius;

    /**
     * Activated when hovering over a residue to display its name.
     */
    void OnToolTip();

    /**
     * Called if the xfit compatible mouse mode is activated.
     */
    void MouseMoveXfit(unsigned short nFlags, CPoint point);

    /**
     * Return a pointer to the viewpoint.
     */
    ViewPoint *GetViewPoint();
    /**
     * Check whether the viewpoint has benn changed since the last draw.
     */
    bool ViewChanged();
    /**
     * Clear the atom stack.
     */
    void ClearAtomStack();
    /**
     * Read the atom stack from a file.
     * @param pathname the filename.
     * @param append if false the stack is overwritten, if true then the read items are appended to the stack.
     */
    void ReadStack(const char *pathname, bool append = true);
    /**
     * Write the stack to an archive file to save for later.
     */
    void WriteStack(CArchive &ar);

    /**
     * Make the canvas cover the entire screen.
     * @param full if true then make canvas full screen, if false restore the canvas size.
     */
    void FullScreen(bool full);
    /**
     * Redraw the canvas.
     */
    void ReDraw();
    /**
     * Redraw the canvas and the other windows associated with the view.
     */
    void ReDrawAll();
    /**
     * return the value of Update.
     */
    int GetUpdate();
    /**
     * reset the value of Update to false.
     */
    void ResetUpdate();
    /**
     * Apply the fit and make it permanent.
     */
    void ApplyFit();
    /**
     * Put the fit atoms back where they started but stay in fitting mode.
     */
    void ResetFit();

    /**
     * Ask the user if he wants to change the render mode.
     * Used when the user asks for an operation that he won't see in the other
     * rendering modes.
     */
    bool DoYouWantBallAndCylinder();
    /**
     * return whether we are currently drawing.
     * Prevents a callback from interupting drawing.
     */
    bool IsDrawing();


    /**
     * Called by the system when the canvas needs drawing.
     */
    void draw(QPainter *painter);

    /**
     * main frame will call this when this widget is activated.
     */
    void OnActivated();


    /**
     * Fit the torsion named.  The trsion name must be in the dictionary.
     */
    bool FitTorsion(const char *torsion_name);
    /**
     * Implement a right mouse drag.
     */
    void RightMouseDrag(CPoint d, float xang, float yang, float zang);

    bool DraggingRefiAtom;
    bool DontConfirmWater;
    void ColorModel(Molecule *model);

    void CheckCenter();

    bool AutoSave;
    /**
     * Parse and execute a script command.
     * If successful returns 1.
     * If error returns -1.
     * If command not found returns 0.
     */
    int ScriptCommand(const char *command);
    /**
     * Checkpoints the model for saving when the model will be changed too
     * much for a simple undo to recover.
     */
    bool SaveModelFile(Molecule *model, const char *description);

    /**
     * Return true if symmetry models should be shown as links (C-alpha only).
     */
    bool link_symm;

    /**
     * Find the next baton position in list.  Keypress ' ' when in the add MRK mode.
     */
    void OnNextBatonPosition();
    /**
     * Insert a MRK residue.
     * @param where the residue to insert before/after.
     * @param model the model into which the MRK is added.
     * @param after if true the resideu is added after where, otherwise before where.
     */
    void InsertMRK(chemlib::Residue *where, Molecule *model, bool after);
    void UpdateCurrent();

    /**
     * Mouse handler for a mouse movement event.
     */
    void OnMouseMove(unsigned short nFlags, CPoint point);
    /**
     * Mouse handler for a mouse left-button double-click event.
     */
    void OnLButtonDblClk(unsigned short nFlags, CPoint point);
    /**
     * Keyboard event handler.
     */
    bool OnKeyDown(unsigned int nChar, unsigned int nRepCnt, unsigned int nFlags);
    /**
     * handler for a window size event.
     */
    void resizeEvent(QResizeEvent *evt);

    virtual QSize sizeHint() const
    {
        return QSize(320, 200);
    }

    /**
     * Prompts for phase file and then calls mapLoadfromphsfile(string).
     */
    void mapLoadfromphsfile();
    /**
     * Loads the given phase file. If crystal information is not available
     * (file is not mtz file and model is not loaded), the user is prompted for
     * crystal information.
     */
    void mapLoadfromphsfile(const std::string &file);

    /**
     * Prompts for map file and then calls mapLoadfromfile(string).
     */
    void mapLoadfromfile();
    /**
     * Loads the given map file. If crystal information is not available
     * (model is not loaded), the user is prompted for crystal information.
     */
    void mapLoadfromfile(const std::string &file);
    void doMapContour(EMapBase *emap);

    void RemoveFileFromHistory(const std::string &pathname);
    void solidSurfaceCommand(int id, std::vector<Molecule*> &mols, std::vector<unsigned int> &selected);

    void deleteAtom();

    void moveTo(float x, float y, float z);

public slots:
    void connectToModel(Molecule *model);
    void modelAdded(Molecule *model);
    void currentModelChanged(Molecule *oldModel, Molecule *newModel);
    void connectToMap(EMapBase *map);
    void mapAdded(EMap *map);
    void annotationAdded(Molecule *model, Annotation *annotation);
    void annotationDeleted(Molecule *model);
    void annotationChanged(Annotation *annotation);
    void atomLabelChanged(Molecule *model, ATOMLABEL *label);
    void surfaceChanged(Molecule *model);
    void atomChanged(chemlib::MIMoleculeBase *model, chemlib::MIAtomList &atom);
    void atomsDeleted(chemlib::MIMoleculeBase *model);
    void modelChanged(chemlib::MIMoleculeBase *model);
    void mapContourLevelsChanged(EMapBase *map);
    void mapVisibilityChanged(EMapBase *map);

    void modelAtomsToBeDeleted(chemlib::MIMoleculeBase *model, const chemlib::MIAtomList &atoms);
    void modelResiduesToBeDeleted(chemlib::MIMoleculeBase *model, std::vector<chemlib::Residue*> &res);
    void moleculeToBeDeleted(chemlib::MIMoleculeBase *model);

    void symmetryToBeCleared(chemlib::MIMoleculeBase *mol);

public:
    void recenter(chemlib::Residue *residue, chemlib::MIAtom *atom);

    void select(Molecule *model, chemlib::Residue *residue, chemlib::MIAtom *atom, bool label = true);

signals:
    void focusResidueChanged(chemlib::Residue*);
public:

    void doRefresh();

    void acceptRefine();
    void cancelRefine();

    void handleKey_space(bool spaceKeyDown);

    using QGraphicsView::mousePressEvent;
    using QGraphicsView::mouseMoveEvent;
    using QGraphicsView::mouseReleaseEvent;
    using QGraphicsView::mouseDoubleClickEvent;
    void handleMousePress(QGraphicsSceneMouseEvent *e);
    void handleMouseRelease(QGraphicsSceneMouseEvent *e);
    void handleMouseMove(QGraphicsSceneMouseEvent *e);
    void handleMouseDoubleClick(QGraphicsSceneMouseEvent *e);
    void keyPressEvent(QKeyEvent *e);
    void enterEvent(QEvent *e);
    void leaveEvent(QEvent *e);

    void closeEvent(QCloseEvent *evt);


    /**
     * Replace the residue on the top of the stack with a new type.
     * @param type a residue type in the dictionary to replace with.
     */
    void replaceResidue(const char *residueType);
    void replaceResidueWithPrompt();
    void fitResidue();
    void refineConformer();
    void replaceAndFitResidue(const char *residueType);
    void replaceAndFitResidueWithPrompt();

    bool promptForReplaceResidue(std::string &value);

    void gotoXyzWithPrompt();
    void gotoXyz(float x, float y, float z);
    void addResidueWithPrompt();
    void findLigandFit();
    float densityForCurrentAtoms();
    void generateSymmAtoms();
    void clearSymmAtoms();
    void saveSymmAtoms();

    void updateRamachandranPlot();

    void clearFitTorsion();
    void cancelFitting();

    void polyAla();
    void polyAlaChain();
    void deleteResidueOnTopOfStack();
    void refineRange();
    void clearStack();
    void readLowerSequenceWithPrompt();
    void readLowerSequence(const std::string &file);
    void modelReplaceWithSequence();
    void fitResidueRange();

    void clearPentamer();

    void revertToSavedModelWithPrompt();
    void revertToSavedModel(int choice);

    void promptSecondaryStructureOptions(const std::string &type);

public:
    unsigned int CurrentSolidSurface()
    {
        return solidsurf_current_surface;
    }
    void SetCursor(unsigned int);
    void OpenAnyFile(const std::string &fname);

    const QActionGroup *solidSurfCommonActionGroup() const
    {
        return solidsurfCommonActionGroup;
    }

    QAction *solidSurfMenuAction() const;

public slots:

    void updatePopupMenu();
    void updateDotSurfMenu();
    void updateSolidsurfMenu();
    void solidSurfaceActionTriggered(QAction *action);

    void OnPrint();

    //@{
    // Menu callback for View/Clear Labels.
    //@}
    void OnEditClearlabels();
    //@{
    // Menu callback for View/Edit Labels...
    //@}
    void OnEditLabels();
    //@{
    // Callback for File/Add Model to Current Document.
    //@}
    void OnFileSave();
    void OnFileSaveAs();

    //@{
    // Update menu callback for Object/Add Annotation to Model.
    //@}
    void OnAnnotation();
    void OnMapAddFree();

    //@{
    // Update menu callback for Object/Move Annotation to Center.
    //@}
    void OnUpdateMoveAnnotationCenter(QAction *action);
    //@{
    // Menu callback for Object/Move Annotation to Center.
    //@}
    void OnMoveAnnotationCenter();
    //@{
    // Update menu callback for Object/MoveAnnotation to Picked Atom.
    //@}
    void OnUpdateMoveAnnotationAtom(QAction *action);
    //@{
    // Menu callback for Object/MoveAnnotation to Picked Atom.
    //@}
    void OnMoveAnnotationAtom();

    void OnUpdateExportModel(QAction *action);
    void OnExportModel();
    void OnUpdateAnnotation(QAction *action);
    //@{
    // Create a new, empty model.
    //@}
    void OnNewModel();

    void OnUpdateMapReindex(QAction *action);
    void OnMapReindex();
    void OnUpdateFitSplitFit(QAction *action);
    void OnFitSplitFit();
    void OnUpdateFitSplitTorsion(QAction *action);
    void OnFitSplitTorsion();
    void OnUpdateFitRange(QAction *action);
    void OnFitRange();
    void OnUpdateRefiResidue(QAction *action);
    void OnRefiResidue();
    void OnUpdateFullScreen(QAction *action);
    void OnFullScreen();
    void OnFindLigandDensity();
    void OnUpdateFindLigandDensity(QAction *action);
    void OnUpdateMapCenterVisibleDensity(QAction *action);
    void OnMapCenterVisibleDensity();
    void OnUpdateLabelEveryNth(QAction *action);
    void OnLabelEveryNth();
    void OnClearGeomAnnotations();
    void OnUpdateFindGeomErrors(QAction *action);
    void OnFindGeomErrors();
    /**
     * Handler for 'W' key which places an H2O at the poistion of the cursor
     * and then refines the position.
     */
    void OnAddWaterAtCursor();
    void OnViewLoad();
    void OnViewSave();
    /**
     * Update menu function for hiding the model.
     */
    void OnUpdateHideModel(QAction *action);

    void OnRefiLigandFit();
    void OnUpdateRefiLigandFit(QAction *action);

    /**
     * Update menu callback for Model/Add Waters.
     */
    void OnUpdateAddWater(QAction *action);
    /**
     * Menu callback for Model/Add Waters.
     */
    void OnAddWater();
    /**
     * Handler for the key ']' to center of the last residue in the chain or fragment.
     */
    void OnGotoCter();
    /**
     * Update menu Handler to center of the last residue in the chain or fragment.
     */
    void OnUpdateGotoCter(QAction *action);
    /**
     * Handler for the key '[' to center of the first residue in the chain or fragment.
     */
    void OnGotoNter();
    /**
     * Update menu Handler to center of the last residue in the chain or fragment.
     */
    void OnUpdateGotoNter(QAction *action);

    /**
     * Menu callback for Model/Build CB Range.
     */
    void OnBuildCBRange();
    /**
     * Update menu callback for Model/Checkpoint Model.
     */
    void OnUpdateCheckpointModel(QAction *action);
    /**
     * Menu callback for Model/Checkpoint Model.
     */
    void OnCheckpointModel();
    void OnUpdateAutoCheckpointModel(QAction *action);
    void OnAutoCheckpointModel();
    /**
     * Update menu callback for Model/Revert Model.
     */
    void OnUpdateRevertModel(QAction *action);
    /**
     * Menu callback for Model/Revert Model.
     */
    void OnRevertModel();
    /**
     * Update menu callback for Model/Replace with Sequence.
     */
    void OnUpdateReplaceSequence(QAction *action);
    /**
     * Menu callback for Model/Replace with Sequence.
     */
    void OnReplaceSequence();
    /**
     * Update menu callback for Sequence/Set Sequence Chain.
     */
    void OnUpdateSequencePositionChain(QAction *action);
    /**
     * Update menu callback for Sequence/Set Sequence Range.
     */
    void OnUpdateSequencePosition(QAction *action);
    /**
     * Menu callback for Sequence/Set Sequence Chain.
     */
    void OnSequencePositionChain();
    /**
     * Menu callback for Sequence/Set Sequence Range.
     */
    void OnSequencePosition();
    /**
     * Update menu callback for Sequence/Save Model Sequence.
     */
    void OnUpdateSequenceSaveModel(QAction *action);
    /**
     * Menu callback for Sequence/Save Model Sequence.
     */
    void OnSequenceSaveModel();
    void OnBuildMainchainRange();
    void OnUpdatePolyAlaChain(QAction *action);
    void OnPolyAlaChain();
    void OnUpdatePolyAla(QAction *action);
    void OnPolyAla();
    /**
     * Open a graph window with a phi-psi/Ramachandrn plot of the current model.
     */
    void OnRamachandranPlotShowAllowed();
    void OnUpdateRamachandranPlotShowAllowed(QAction *action);
    void OnUpdateFitSurfaceProbe(QAction *action);
    void OnUpdateFitSurfaceExtended(QAction *action);
    void OnUpdateFitSurfaceVdw(QAction *action);
    void OnUpdateFitSurfaceNone(QAction *action);
    /**
     * Menu callback for Fit/Surface Fit Atoms/Extended Surface.
     */
    void OnFitSurfaceExtended();
    /**
     * Menu callback for Fit/Surface Fit Atoms/Contact Surface.
     */
    void OnFitSurfaceProbe();
    /**
     * Menu callback for Fit/Surface Fit Atoms/Van der Waal Surface.
     */
    void OnFitSurfaceVdw();
    /**
     * Menu callback for Fit/Surface Fit Atoms/No Surface.
     */
    void OnFitSurfaceNone();
    void OnUpdateAddMarkBefore(QAction *action);
    /**
     * Menu callback for Model/Add MRK Before.
     */
    void OnAddMarkBefore();
    void OnUpdateAddMarkAfter(QAction *action);
    /**
     * Menu callback for Model/Add MRK After.
     */
    void OnAddMarkAfter();
    /**
     * Hide the model.
     */
    void OnHideModel();
    void OnUpdateFlipPeptide(QAction *action);
    /**
     * Menu callback for Fit/Fix Backbone/Flip Peptide.
     */
    void OnFlipPeptide();
    void OnUpdateMarkBeta(QAction *action);
    void OnUpdateMarkAlpha(QAction *action);
    /**
     * Set the refine target for this residue as beta sheet.
     */
    void OnMarkBeta();
    /**
     * Set the refine target for this residue as alpha helix.
     */
    void OnMarkAlpha();
    void OnUpdateReplaceAll(QAction *action);
    /**
     * Menu callback for Fit/Fix Backbone/Replace All 5.
     */
    void OnReplaceAll();
    void OnUpdateReplaceLast4(QAction *action);
    /**
     * Menu callback for Fit/Fix Backbone/Replace Last 4.
     */
    void OnReplaceLast4();
    void OnUpdateReplaceFirst4(QAction *action);
    /**
     * Menu callback for Fit/Fix Backbone/Replace First 4.
     */
    void OnReplaceFirst4();
    void OnUpdateReplaceMiddle3(QAction *action);
    /**
     * Menu callback for Fit/Fix Backbone/Replace Middle 3.
     */
    void OnReplaceMiddle3();
    void OnUpdateClearPentamer(QAction *action);
    /**
     * Menu callback for Fit/Fix Backbone/Clear Backbone Match.
     */
    void OnClearPentamer();
    void OnUpdateFindPentamer(QAction *action);
    /**
     * Menu callback for Fit/Fix Backbone/Suggest Backbone Match.
     */
    void OnFindPentamer();
    /**
     * Callback for Show/Ribbon Colors
     */
    void OnSetRibbonColors();
    void OnUpdateFitReplaceAndFit(QAction *action);
    /**
     * Callback for Model/Fit and Replace.
     */
    void OnFitReplaceAndFit();
    void OnUpdateRefiRigidBody(QAction *action);

    /**
     * Callback for Refi/Rigid Body Refine Current Atoms.
     */
    void OnRefiRigidBody();
    void OnUpdateFitRedo(QAction *action);
    void OnUpdateFitUndo(QAction *action);
    /**
     * Callback for Fit/Redo.  Goes back up the Undo Stack.
     */
    void OnFitRedo();
    /**
     * Callback for Fit/Undo.  Goes down up the Undo Stack.
     */
    void OnFitUndo();

    void OnUpdateRefiRedo(QAction *action);
    /**
     * Callback for Refi/Redo.
     */
    void OnRefiReDo();
    void OnUpdateRefiReset(QAction *action);
    void OnUpdateRefiCancel(QAction *action);
    void OnUpdateRefiAccept(QAction *action);
    /**
     * Callback for Refi/Reset.
     */
    void OnRefiReset();
    /**
     * Callback for Refi/Cancal.
     */
    void OnRefiCancel();
    /**
     * Callback for Refi/Accept.
     */
    void OnRefiAccept();
    void OnUpdateRenderLinesmooth(QAction *action);
    /**
     * Callback for Render/Smooth Lines.
     */
    void OnRenderLinesmooth();
    /**
     * Callback for Refi/Refine Molecule.
     */
    void OnRefiMolecule();
    /**
     * Callback for Refi/Refine Options.
     */
    void OnRefiOptions();
    /**
     * Callback for Refi/Undo.
     */
    void OnRefiUndo();
    /**
     * Callback for Refi/Refine Range.
     */
    void OnRefiRange();
    /**
     * Callback for Refi/Refine Region.  A region is the residue plus the one before and after.
     */
    void OnRefiRegion();
    void OnUpdateRefiMolecule(QAction *action);
    void OnUpdateRefiUndo(QAction *action);
    void OnUpdateRefiRange(QAction *action);
    void OnUpdateRefiRegion(QAction *action);
    void OnUpdateMapSwitch(QAction *action);
    /**
     * Callback for Map/Switch Map.  Cahnges teh current map to the one picked from a list.
     */
    void OnMapSwitch();
    void OnUpdateMapContourLevels(QAction *action);
    /**
     * Callback for Map/Contour.
     */
    void OnMapContourLevels();
    void OnUpdateMapFFT(QAction *action);
    /**
     * Callback for Map/FFT Phases.
     */
    void OnMapFFT();
    void OnUpdateMapSFCalc(QAction *action);
    /**
     * Callback for Map/Calculate Structure Factors.
     */
    void OnMapSFCalc();

    void OnContourListFile();

    void OnRenderTargetSize();

    /**
     * handler for a timer event.
     */
    void OnHardwareStereo();

    void OnStereoToggle();
    /**
     * Callback for decreasing perspective.
     */
    void OnDecreasePersp();
    /**
     * Callback for increasing perspective.
     */
    void OnIncreasePersp();
    /**
     * Callback for setting perspective to 0.
     */
    void OnViewOrthonormal();
    void OnUpdateDecreasePersp(QAction *action);
    /**
     * Callback for Geom/Angle.
     */
    void OnGeometryAngle();
    /**
     * Callback for Geom/Distance.
     */
    void OnGeometryDistance();
    /**
     * Callback for Geom/Dihedral.
     */
    void OnGeometryTorsion();
    /**
     * Callback for View/Slab In.
     */
    void OnViewSlabin();
    /**
     * Callback for View/Slab Out.
     */
    void OnViewSlabout();
    /**
     * Callback for View/Atom Stack.
     */
    void OnViewAtomstack();
    /**
     * Callback for View/Gnomon.
     */
    void OnViewGnomon();
    void OnViewUnitCell();
    /**
     * Callback for View/Animate/Roll.
     */
    void OnAnimateRoll();
    /**
     * Callback for View/Animate/Rock.
     */
    void OnAnimateRock();

    void setViewpointLineThickness(int thickness);

    void OnRenderingDepthcuedcolors();
    void OnRenderingDepthcuedlinewidth();

    void OnDimNonactiveModels();
    void OnAmountToDimNonactiveModels();
    void OnUpdateDimNonactiveModels(QAction *action);

    /**
     * Callback for View/Rotate View +90.
     */
    void OnViewRotatey90();
    /**
     * Callback for View/Rotate View -90.
     */
    void OnViewRotateyminus();
    /**
     * Not implemented.  Used to select parts of the model to be colored and shown.
     * Left over from Molw.
     */
    void OnEditSelectmodel();

    void OnUpdateAnimateRock(QAction *action);
    void OnUpdateAnimateRoll(QAction *action);
    void OnUpdateLinethicknessFour(QAction *action);
    void OnUpdateLinethicknessOne(QAction *action);
    void OnUpdateLinethicknessThree(QAction *action);
    void OnUpdateLinethicknessTwo(QAction *action);
    void OnUpdateRenderingDepthcuedcolors(QAction *action);
    void OnUpdateRenderingDepthcuedlinewidth(QAction *action);
    void OnUpdateHardwareStereo(QAction *action);
    void OnUpdateStereoToggle(QAction *action);
    void OnUpdateViewAtomstack(QAction *action);
    void OnUpdateViewGnomon(QAction *action);
    void OnUpdateViewUnitCell(QAction *action);
    void OnUpdateViewOrthonormal(QAction *action);
    /**
     * Callback for Render/Ball-and-cylinder.
     */
    void OnRenderingBallandstick();
    void OnUpdateRenderingBallandstick(QAction *action);
    /**
     * Go to a specific x ,y z coordinate.  Not implemented in this version.
     */
    void OnGotoGotoxyz();
    /**
     * Callback for View/Clipping Planes.  Not implemented in this version.  Needs
     * a dialog box written.
     */
    void OnViewClipplanes();
    /**
     * Callback for View/Edit Labels.  NEEDS IMPLEMENTING!
     */
    void OnViewLabels();
    void OnUpdateViewLabels(QAction *action);
    void OnVdwDotSurface();
    void OnViewClearmessage();
    // void OnObjectsModels();
    void OnFitResidue();
    void OnFitResidues();
    void OnUpdateFitResidues(QAction *action);
    void OnFitRotate();
    void OnFitTorsion();
    void OnUpdateFitTorsion(QAction *action);
    void OnFitTranslate();
    void OnUpdateFitTranslate(QAction *action);
    void OnFitCentermode();
    void OnUpdateFitCentermode(QAction *action);
    void OnFitCancel();
    void OnFitReset();
    void OnFitApply();
    void OnFitSingleatom();
    void OnFitAtoms();
    void OnUpdateFitAtoms(QAction *action);
    void OnFitSetuptorsion();
    void OnUpdateFitSetuptorsion(QAction *action);
    void OnFitCleartorsion();
    void OnFitNextConfomer();
    void OnUpdateFitNextConfomer(QAction *action);
    void OnUpdateFitCleartorsion(QAction *action);
    void OnUpdateFitApply(QAction *action);
    void OnUpdateFitCancel(QAction *action);
    void OnUpdateFitReset(QAction *action);
    void OnUpdateFitRotate(QAction *action);
    void OnFitReplacewith();
    void OnUpdateFitReplacewith(QAction *action);
    void OnObjectsAllatoms();
    void OnViewTopview();
    void OnUpdateViewTopview(QAction *action);
    void OnRenderSpacefilling();
    void OnUpdateRenderSpacefilling(QAction *action);
    void OnRenderSticks();
    void OnUpdateRenderSticks(QAction *action);
    void OnGeomBond();
    void OnUpdateGeomBond(QAction *action);
    void OnGeomUnbond();
    void OnUpdateGeomUnbond(QAction *action);
    void OnUpdateGeometryAngle(QAction *action);
    void OnUpdateGeometryDistance(QAction *action);
    void OnUpdateGeometryTorsion(QAction *action);
    void OnRenderBallsize();
    void OnUpdateRenderBallsize(QAction *action);
    void OnRenderBallandcylinder();
    void OnUpdateRenderBallandcylinder(QAction *action);
    void OnGotoZoomiin();
    void OnGotoZoomout();
    void OnGotoFittoscreen();
    void OnObjectClearstack();
    void OnEditCopy();
    void OnAnimateRockandrollparameters();
    void OnUpdateFitResidue(QAction *action);
    void OnUpdateFitSingleatom(QAction *action);
    void OnObjectShowresidue();
    void OnUpdateObjectShowresidue(QAction *action);
    void OnObjectShowsidechain();
    void OnUpdateObjectShowsidechain(QAction *action);
    void OnGeomNeighbours();
    void OnUpdateGeomNeighbours(QAction *action);
    void OnGeomClearneighbours();
    void OnGeomHbonds();
    void OnGeomClearhbonds();
    void OnObjectBackboneribbon();
    void OnObjectClearribbon();
    void OnSchematicSecondaryStructure();
    void OnRibbonSecondaryStructure();
    void OnDeleteSecondaryStructure();
    void OnSecondaryStructureOptionsTube();
    void OnSecondaryStructureOptionsSheet();
    void OnSecondaryStructureOptionsTurn();
    void OnSecondaryStructureOptionsRandom();
    void OnSecondaryStructureOptionsHelix();
    void OnGeomAddsinglehbond();
    void OnUpdateGeomAddsinglehbond(QAction *action);
    void OnViewContacts();
    void OnUpdateViewContacts(QAction *action);
    void OnFitLsqsuperpose();
    void OnFitFitmolecule();
    void OnUpdateFitFitmolecule(QAction *action);
    void OnFitDeleteresidue();
    void OnUpdateFitDeleteresidue(QAction *action);
    void OnFitInsertresidue();
    void OnUpdateFitInsertresidue(QAction *action);
    void OnFitRenameresidue();
    void OnUpdateFitRenameresidue(QAction *action);
    void OnGotoFitalltoscreen();
    void OnShowBackboneAsCATrace();
    void OnShowBackboneAsAtoms();
    void OnAnimateBlink();
    void OnUpdateAnimateBlink(QAction *action);
    void OnObjectSurfaceresidue();
    void OnUpdateObjectSurfaceresidue(QAction *action);
    void OnObjectSurfaceClearsurface();
    void OnObjectSurfaceresidues();
    void OnUpdateObjectSurfaceresidues(QAction *action);
    void OnUpdateObjectSurfaceatoms(QAction *action);
    void OnObjectSurfaceatom();
    void OnUpdateObjectSurfaceAtom(QAction *action);
    void OnObjectSurfaceAtoms();
    void OnObjectStackExpandtopallatomsinresidue();
    void OnUpdateObjectStackExpandtopallatomsinresidue(QAction *action);
    void OnUpdateObjectStackExpandtop2residues(QAction *action);
    void OnObjectStackExpandtop2residues();
    void OnObjectStackExpandtop2allatomsinrange();
    void OnUpdateObjectStackExpandtop2allatomsinrange(QAction *action);
    void OnObjectShowresidues();
    void OnUpdateObjectShowresidues(QAction *action);
    void OnObjectShowresiduerange();
    void OnObjectShowsidechainrange();
    void OnUpdateObjectShowresiduerange(QAction *action);
    void OnUpdateObjectShowsidechainrange(QAction *action);
    void OnObjectShowsidechains();
    void OnUpdateObjectShowsidechains(QAction *action);
    void OnMapLoadfromphsfile();
    void OnMapLoadfromfile();
    void OnMapContour();
    void OnUpdateMapContour(QAction *action);
    void OnShowShowwithinsphere();
    void OnObjectWhenshowncolor();
    void OnObjectRadiusisBall();
    void OnObjectRadiusisCpk();
    void OnObjectRadiusisCylinder();
    void OnUpdateObjectRadiusisBall(QAction *action);
    void OnUpdateObjectRadiusisCpk(QAction *action);
    void OnUpdateObjectRadiusisCylinder(QAction *action);
    void OnUpdateObjectWhenshowncolor(QAction *action);
    void OnObjectResiduerangeColor();
    void OnUpdateObjectResiduerangeColor(QAction *action);
    void OnObjectResiduerangeRadius();
    void OnObjectResidueColor();
    void OnObjectResidueRadius();
    void OnUpdateObjectResiduerangeRadius(QAction *action);
    void OnObjectResiduerangeTurnoff();
    void OnObjectResiduesColor();
    void OnUpdateObjectResiduesColor(QAction *action);
    void OnObjectResidueTurnoff();
    void OnObjectResiduesRadius();
    void OnUpdateObjectResiduerangeTurnoff(QAction *action);
    void OnUpdateObjectResiduesRadius(QAction *action);
    void OnObjectResiduesTurnoff();
    void OnUpdateObjectResiduesTurnoff(QAction *action);
    void OnObjectAtomColor();
    void OnObjectAtomRadius();
    void OnUpdateObjectAtomColor(QAction *action);
    void OnUpdateObjectAtomRadius(QAction *action);
    void OnObjectAtomsColor();
    void OnUpdateObjectAtomsColor(QAction *action);
    void OnObjectAtomsRadius();
    void OnUpdateObjectAtomsRadius(QAction *action);
    void OnShowColorallatoms();
    void OnUpdateObjectResidueColor(QAction *action);
    void OnUpdateObjectResidueRadius(QAction *action);
    void OnUpdateObjectResidueTurnoff(QAction *action);
    void OnShowUndocolorradius();
    void OnUpdateShowUndocolorradius(QAction *action);
    void OnObjectStackDeletetopitem();
    void OnViewUndo();
    void OnUpdateViewUndo(QAction *action);
    void OnShowPickedatomTurnon();
    void OnUpdateShowPickedatomTurnon(QAction *action);
    void OnShowPickedatomTurnoff();
    void OnUpdateShowPickedatomTurnoff(QAction *action);
    void OnShowAllpickedatomsTurnoff();
    void OnUpdateShowAllpickedatomsTurnoff(QAction *action);
    void OnShowAllpickedatomsTurnon();
    void OnUpdateShowAllpickedatomsTurnon(QAction *action);
    void OnShowRadiusmodel();
    void OnUpdateShowRadiusmodel(QAction *action);
    void OnUpdateShowColorallatoms(QAction *action);
    void OnUpdateSurfaceSolvent(QAction *action);
    void OnShowHideBackbone();
    void OnShowSidechainAtoms();
    void OnHideSidechainAtoms();
    void OnShowHidehydrogens();
    void OnUpdateShowHidehydrogens(QAction *action);
    void OnUpdateObjectSurfaceSpherearoundatom(QAction *action);
    void OnSurfaceSolvent();
    void OnObjectSurfaceSpherearoundatom();
    void OnSequenceEnter();
    void OnUpdateSequenceEnter(QAction *action);
    void OnSequenceRead();
    void OnUpdateSequenceRead(QAction *action);
    void OnSequenceSave();
    void OnUpdateSequenceSave(QAction *action);
    void OnSequenceInsertgap();
    void OnUpdateSequenceInsertgap(QAction *action);
    void OnSequenceDeletegap();
    void OnUpdateSequenceDeletegap(QAction *action);
    void OnSequenceInsertlowergap();
    void OnUpdateSequenceInsertlowergap(QAction *action);
    void OnSequenceDeletelowergap();
    void OnUpdateSequenceDeletelowergap(QAction *action);

    void OnInvertChiralCenter();
    void OnUpdateInvertChiralCenter(QAction *action);
    void OnViewChiralCenters();
    void OnUpdateViewChiralCenters(QAction *action);

    void OnSolidSurfaceCommand(QAction *action);

    void OnExportImage();

    void OnDeleteAtom();
    void OnUpdateDeleteAtom(QAction *action);

    void OnShowSymmAtomsAsAtoms();
    void OnShowSymmAtomsAsCATrace();
    void OnShowHideSymmAtoms();
    void OnShowSaveSymmAtoms();
    void OnUpdateShowSaveSymmAtoms(QAction *action);

};

inline bool MIGLWidget::IsFitting()
{
    return CurrentAtoms.size() > 0;
}

inline bool MIGLWidget::IsTorsioning()
{
    return ( Torsioning > 0 && CurrentAtoms.size() > 0 && fitmol);
}


inline ViewPoint*MIGLWidget::GetViewPoint()
{
    return viewpoint;
}

inline int MIGLWidget::GetUpdate()
{
    return Update;
}

inline void MIGLWidget::ResetUpdate()
{
    Update = 0;
}

inline bool MIGLWidget::IsDrawing()
{
    return is_drawing;
}

class ViewController : public QObject
{
    Q_OBJECT

    friend class MIGLWidget;

    ViewController() : QObject()
    {
    }

public:
    static ViewController *instance();

signals:
    void viewActivated(MIGLWidget *view);
    void viewDeactivated(MIGLWidget *view);
};

#endif // ifndef mifit_ui_MIGLWidget_h
