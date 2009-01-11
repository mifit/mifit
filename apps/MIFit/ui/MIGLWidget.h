#ifndef mifit_ui_MIGLWidget_h
#define mifit_ui_MIGLWidget_h

#include <math/Quaternion.h>
#include <boost/signal.hpp>
#include <cmath>

#include "corelib.h"
#include "chemlib.h"
#include "macafxwin.h"
#include "SaveModel.h"

#include "MIEventHandler.h"

#ifdef USE_QT_RAMAPLOT
#include "RamaPlot.h"
#endif

#ifdef _WIN32
#include <time.h>
#endif


class EMap;
class Displaylist;
class MIGLWidget;
class BatchJob;
class MIMenu;
class MIMenuBar;

class CMolwViewAnnotationPickingRenderable;
class CMolwViewAtomPickingRenderable;
class CMolwViewBondPickingRenderable;
class CMolwViewSlabPickingRenderable;
class CMolwViewScene;

class GLRenderer;
class QResizeEvent;
class InterpBox;
#ifdef USE_SEQNAV_WINDOW
class SequenceWindow;
#endif //USE_SEQ_WINDOW
#ifdef USE_NAV_WINDOW
class NavigatorCanvas;
#endif //USE_NAV_WINDOW
namespace mi {
namespace opengl {
class Camera;
class Frustum;
class Scene;
class StereoView;
class Viewport;
namespace interact {
class MousePicker;
class TargetFeedback;
}
}
}

class MAP_POINT;

#include <QGLWidget>
#include <QMdiSubWindow>
#include <QVBoxLayout>

#define sign(a) (a < 0 ? -1 : 1)


class MIGLWidget : public QGLWidget, public MIEventHandler {
  Q_OBJECT

  //Display list contains the list of objects owned by this document
  Displaylist* Models;
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
  bool        _modified;

  mi::opengl::StereoView* stereoView;
  CMolwViewScene* scene;
  mi::opengl::Frustum* frustum;
  mi::opengl::Camera* camera;
  mi::opengl::interact::TargetFeedback* targetFeedback;
  mi::opengl::interact::MousePicker* mousePicker;

  CMolwViewAnnotationPickingRenderable* annotationPickingRenderable;
  CMolwViewAtomPickingRenderable* atomPickingRenderable;
  CMolwViewBondPickingRenderable* bondPickingRenderable;
  CMolwViewSlabPickingRenderable* slabPickingRenderable;

  bool is_drawing;
  chemlib::RESIDUE* PentamerStart;
  int batonposition;
  chemlib::MIAtom* batonatom;

  int m_timer;
  int ClockTick;
  int TimerOn;
  public:
  void timerEvent(QTimerEvent *e);
  private:
  QWidget *fullScreenWidget;
  QVBoxLayout *fullScreenLayout;
  QMdiSubWindow *oldParent;
  QWidget *containerWidget;

  chemlib::RESIDUE* focusres;
  bool focusresDeleted;
  bool DragStart;
  bool FromMouse;

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

  QMenu* popup_menu;
  QAction* fitResidueAction;
  QAction* rotateAction;
  QAction* translateAction;
  QAction* setTorsionAction;
  QAction* acceptFitAction;
  QAction* cancelFitAction;
  QAction* replaceAndFitAction;
  QAction* deleteResidueAction;
  QAction* deleteAtomAction;
  QAction* refineRegionAction;
  QAction* bestConformerAction;
  QAction* acceptRefineAction;
  QAction* cancelRefineAction;
  QAction* saveAction;
  QAction* fullscreenAction;
  QAction* clearLabelsAction;

  QAction* dotSurfVdwAction;
  QAction* dotSurfSolventExposedAction;
  QAction* dotSurfAtomSphereAction;
  QAction* dotSurfResidueAction;
  QAction* dotSurfResiduesAction;
  QAction* dotSurfAtomAction;
  QAction* dotSurfAtomsAction;
  QAction* dotSurfClearAction;

  QMenu* solidsurf_popup_menu;
  QAction *solidsurfActions[10];
  QActionGroup* solidsurfActionGroup;
  QActionGroup* solidsurfCommonActionGroup;

  QMenu* newSolidsurfMenu(bool include_selection_items = true);
  QAction* solidsurf_menu_action(QMenu* menu, QActionGroup* group, int actionId, const char* text);

  unsigned int solidsurf_current_surface;
  int SelectType;

#ifdef USE_NAV_WINDOW
  NavigatorCanvas* navigator;
  void createNavigatorWindow();
#endif
#ifdef USE_SEQ_WINDOW
  void createSequenceWindow();
#endif

  unsigned int DoingRange;

  unsigned int ClusterSize;


  void prepareAtomPicking();
  void prepareBondPicking();

  void doSlabDrag(int x, int y, int dx, int dy);

  void findResidueAndMoleculeForAtom(chemlib::MIAtom* atom, chemlib::RESIDUE*& res, Molecule*& mol);
  void doSetFocusResidue(chemlib::RESIDUE* res);

  /**
   * Get rid of any pointers or objects pointing to model.
   * Prevents accessing bad pointers leading to a seg fault.
   */
  void Purge(Molecule* model);
  void PurgeChain(chemlib::RESIDUE* chain);
  void Purge(chemlib::RESIDUE* res);
  void Purge(chemlib::MIAtom* atom);

  /**
   * List of checkpoint saved models.
   */
  std::vector<SaveModel> SaveModels;

public:
  MIGLWidget();
  virtual ~MIGLWidget();

  void SetFilename(const std::string &fname) { _filename=fname; }
  std::string GetFilename() const { return _filename; }

  void SetTitle(const std::string &title);
  std::string GetTitle() const { return _title; }

  void SetDocumentSaved(bool saved) { _modified=!saved; }


  bool LoadScript(const char* path);
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
  bool LoadMap(const char* command, const char* path, const std::string& columnNames);
  void exportCurrentModel();

  //@{
  // Save the document into a session (.mlw file).
  //@}
  bool SaveDocument(const std::string& pathname);
  //@{
  // Returns true if the file is an XML .mlw file.
  // Otherwise it is assumed to
  // be the older Molw .mlw file format.
  //@}
  bool IsXMLDocument(const char* pathname);
  //@{
  // Returns true if the fie is XML and contains a <Molecule> tag.
  //@}
  bool IsXMLMolecule(io& file);
  //@{
  // Load an Molecule from an XML .mlw file.
  //@}
  bool LoadXMLMolecule(const char* pathname, int nmodel);
  //@{
  // Purge all references to an EMap.
  // To be executed before deleting an EMap.
  //@}
  void Purge(EMap*);

  //@{
  // Return a pointer to the display list.
  //@}
  Displaylist* GetDisplaylist() {
    return Models;
  }

  //@{
  // Return true if this a newly loaded file.
  //@}
  int NewFile() {
    return (newfile);
  }

  //@{
  // Clear the newfile flag.
  //@}
  void ResetNewFile() {
    newfile = 0;
  }

  //@{
  // Set the newfile flag on.
  //@}
  void SetNewFile() {
    newfile = 1;
  }

  //@{
  //
  //@}
  //	void OnSelchangeCombo1();
  //@{
  // Load a model from a PDB format file.
  //@}
  bool LoadPDBFile(const char* pathname);
  //@{
  // Load a dictionary residue.
  //@}
  chemlib::RESIDUE* GetDictRes(const char* key, int confomer = 0);
  //@{
  // Load a dictionary residue.
  //@}
  chemlib::RESIDUE* GetDictRes(const char single, int confomer = 0);
  //@{
  // Return a pointer to the dictionary residue linked list.
  //@}
  std::vector<std::string> GetDictResList();

  //@{
  // Save the document and view data to an XML .mlw file.
  //@}
  void Serialize(XMLArchive& ar);   // overridden for document i/o
  //@{
  // Run on the creation of a new document.
  //@}
  virtual bool    OnNewDocument();
  //@{
  // Menu callback for File/Open.
  //@}
  virtual bool    OnOpenDocument(const std::string& pathname);
  //@{
  // Menu callback for File/Save.
  //@}

  virtual bool    OnSaveDocument(const std::string& pathname);

  // Prompts for saving if about to close a modified document. Returns true
  // if ok to close the document (may have saved in the meantime, or set
  // modified to false)
  virtual bool OnSaveModified();

  //@{
  // The current annotation.
  //@}
  Annotation* CurrentAnnotation;


  Molecule* PentamerFrom;
  Molecule* PentamerModel;

  /**
   * Full screen flag.
   */
  bool IsFullScreen;
  /**
   * Cluster residue being optimized
   */
  InterpBox* BoundingBox;
  /**
   * string containing the message printed at the bottom of the screen.
   */
  std::string message;
  /**
   * The viewpoint containing the current screen center, rotation and so forth.
   * Created when view is created.
   */
  ViewPoint* viewpoint;
  /**
   * Stack of atoms used by the picking operations.
   * Created when view is created.
   */
  Stack* AtomStack;
  /**
   * if true then the colors information is saved for an undo operation
   */
  bool SaveColors;
  /**
   * The current residue.
   * If NULL there is currently no focus.
   */
  const chemlib::RESIDUE* getFocusResidue();
  /**
   * Set the current residue.
   */
  void setFocusResidue(chemlib::RESIDUE* res, bool recenter = true, bool deleted = false);

  /**
   * The current model containing the focus residue.
   */
  Molecule* focusnode;
  /**
   * The residue being currently fit, if any.
   * NULL value indicate no residue being fit.
   */
  chemlib::RESIDUE* fitres;
  /**
   * The model being currently fit, if any.
   * NULL value indicate no model being fit.
   */
  Molecule* fitmol;
  /**
   * An atom representing the current center of rotation for fitting operations.
   */
  chemlib::MIAtom fitcenter;
  /**
   * pointer to the atom from which the fitcenter was copied.
   */
  chemlib::MIAtom* fromcenter;
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
   * if true then show the stack.
   */
  bool ShowStack;
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
  unsigned int SetupFit(const std::vector<chemlib::MIAtom*>& atoms, Molecule* model, char altloc = ' ');
  unsigned int SetupFit(chemlib::RESIDUE* start, chemlib::RESIDUE* end, Molecule* model, char altloc = ' ');

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
  ViewPoint* GetViewPoint();
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
  void ReadStack(const char* pathname, bool append = true);
  /**
   * Write the stack to an archive file to save for later.
   */
  void WriteStack(CArchive& ar);

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


#ifdef USE_SEQ_WINDOW
  SequenceWindow* seqwin;
#endif //USE_SEQ_WINDOW

  /**
   * Called by the system when the canvas needs drawing.
   */
  virtual void paintGL();

  /**
   * main frame will call this when this widget is activated.
   */
  void OnActivated();


  /**
   * Fit the torsion named.  The trsion name must be in the dictionary.
   */
  bool FitTorsion(const char* torsion_name);
  /**
   * Implement a right mouse drag.
   */
  void RightMouseDrag(CPoint d, float xang, float yang, float zang);

  bool DraggingRefiAtom;
  bool DontConfirmWater;
  void ColorModel(Molecule* model);

  void CheckCenter();

  bool AutoSave;
  /**
   * Parse and execute a script command.
   * If successful returns 1.
   * If error returns -1.
   * If command not found returns 0.
   */
  int ScriptCommand(const char* command);
  /**
   * Checkpoints the model for saving when the model will be changed too
   * much for a simple undo to recover.
   */
  bool SaveModelFile(Molecule* model, const char* description);

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
  void InsertMRK(chemlib::RESIDUE* where, Molecule* model, bool after);
  void UpdateCurrent();

  /**
   * Mouse handler for a mouse movement event.
   */
  void OnMouseMove(unsigned short nFlags, CPoint point);
  /**
   * Mouse handler for a mouse left-button up event.
   */
  void OnLButtonUp(unsigned short nFlags, CPoint point);
  /**
   * Mouse handler for a mouse left-button down event.
   */
  void OnLButtonDown(unsigned short nFlags, CPoint point);
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

  virtual QSize sizeHint() const {
    return QSize(320, 200);
  }

  void OnRButtonUp(unsigned short nFlags, CPoint point);
  void OnRButtonDown(unsigned short nFlags, CPoint point);

  /**
   * Prompts for phase file and then calls mapLoadfromphsfile(string).
   */
  void mapLoadfromphsfile();
  /**
   * Loads the given phase file. If crystal information is not available
   * (file is not mtz file and model is not loaded), the user is prompted for
   * crystal information.
   */
  void mapLoadfromphsfile(const std::string& file);

  /**
   * Prompts for map file and then calls mapLoadfromfile(string).
   */
  void mapLoadfromfile();
  /**
   * Loads the given map file. If crystal information is not available
   * (model is not loaded), the user is prompted for crystal information.
   */
  void mapLoadfromfile(const std::string& file);
  void doMapContour(EMapBase* emap);

  void RemoveFileFromHistory(const std::string& pathname);
  void solidSurfaceCommand(int id, std::vector<Molecule*> &mols, std::vector<unsigned int> &selected);

  void deleteAtom();

  void moveTo(float x, float y, float z);

  static boost::signal1<void, MIGLWidget*> viewActivated;
  static boost::signal1<void, MIGLWidget*> viewDeactivated;


  void connectToModel(Molecule* model);
  void modelAdded(Molecule* model);
  void currentModelChanged(Molecule* oldModel, Molecule* newModel);
  void connectToMap(EMapBase* map);
  void mapAdded(EMapBase* map);
  void annotationAdded(Molecule* model, Annotation* annotation);
  void annotationDeleted(Molecule* model);
  void annotationChanged(Annotation* annotation);
  void atomLabelChanged(Molecule* model, ATOMLABEL* label);
  void surfaceChanged(Molecule* model);
  void atomChanged(chemlib::MIMoleculeBase* model, std::vector<chemlib::MIAtom*>& atom);
  void atomsDeleted(chemlib::MIMoleculeBase* model);
  void modelChanged(chemlib::MIMoleculeBase* model);
  void mapContourLevelsChanged(EMapBase* map);
  void mapVisibilityChanged(EMapBase* map);

  void modelAtomsToBeDeleted(chemlib::MIMoleculeBase* model, const std::vector<chemlib::MIAtom*> &atoms);
  void modelResiduesToBeDeleted(chemlib::MIMoleculeBase* model, std::vector<chemlib::RESIDUE*> &res);
  void moleculeToBeDeleted(chemlib::MIMoleculeBase *model);

  void symmetryToBeCleared(chemlib::MIMoleculeBase* mol);


  void recenter(chemlib::RESIDUE* residue, chemlib::MIAtom* atom);

  void select(Molecule* model, chemlib::RESIDUE* residue, chemlib::MIAtom* atom, bool label = true);

  boost::signal1<void, chemlib::RESIDUE*> focusResidueChanged;

  bool HandleHistory(MIData& data);
  bool PlaybackMouseHistory(MIData& data);

  void doRefresh();

  void acceptRefine();
  void cancelRefine();

  void handleKey_space(bool spaceKeyDown);

  void mousePressEvent(QMouseEvent *e);
  void mouseReleaseEvent(QMouseEvent *e);
  void mouseMoveEvent(QMouseEvent *e);
  void mouseDoubleClickEvent(QMouseEvent *e);
  void keyPressEvent(QKeyEvent *e);
  void enterEvent(QEvent *e);
  void leaveEvent(QEvent *e);

  /*!
   * \todo MIGLWidget::closeEvent does not get called. Close events only occur
   * for top-level windows.
   */
  void closeEvent(QCloseEvent *evt);


  /**
   * Replace the residue on the top of the stack with a new type.
   * @param type a residue type in the dictionary to replace with.
   */
  void replaceResidue(const char* residueType);
  void replaceResidueWithPrompt();
  void fitResidue();
  void refineConformer();
  void replaceAndFitResidue(const char* residueType);
  void replaceAndFitResidueWithPrompt();

  bool promptForReplaceResidue(std::string& value);

  void gotoXyzWithPrompt();
  void gotoXyz(float x, float y, float z);
  void addResidueWithPrompt();
  void findLigandFit();
  float densityForCurrentAtoms();
  void generateSymmAtoms();
  void clearSymmAtoms();
  void saveSymmAtoms();

#ifdef USE_QT_RAMAPLOT
  void updateRamachandranPlot();
#endif

  void clearFitTorsion();
  void cancelFitting();

  void polyAla();
  void polyAlaChain();
  void deleteResidueOnTopOfStack();
  void refineRange();
  void clearStack();
  void readLowerSequenceWithPrompt();
  void readLowerSequence(const std::string& file);
  void modelReplaceWithSequence();
  void fitResidueRange();

  void clearPentamer();

  void revertToSavedModelWithPrompt();
  void revertToSavedModel(int choice);

  void promptSecondaryStructureOptions(const std::string& type);

public:
  unsigned int CurrentSolidSurface() { return solidsurf_current_surface; }
  void SetCursor(unsigned int);
  void OpenAnyFile(const std::string &fname);

  const QActionGroup* solidSurfCommonActionGroup() const {
    return solidsurfCommonActionGroup;
  }

  QAction* solidSurfMenuAction() const;

public Q_SLOTS:

  void updatePopupMenu();
  void updateDotSurfMenu();
  void updateSolidsurfMenu();
  void solidSurfaceActionTriggered(QAction* action);

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
  void OnFileAddModel();
  void OnFileSave();
  void OnFileSaveAs();

  //@{
  // Update menu callback for Object/Add Annotation to Model.
  //@}
  void OnAnnotation();
  void OnUpdateMapAddFree(const MIUpdateEvent& pCmdUI);
  void OnMapAddFree();

  //@{
  // Update menu callback for Object/Move Annotation to Center.
  //@}
  void OnUpdateMoveAnnotationCenter(const MIUpdateEvent& pCmdUI);
  //@{
  // Menu callback for Object/Move Annotation to Center.
  //@}
  void OnMoveAnnotationCenter();
  //@{
  // Update menu callback for Object/MoveAnnotation to Picked Atom.
  //@}
  void OnUpdateMoveAnnotationAtom(const MIUpdateEvent& pCmdUI);
  //@{
  // Menu callback for Object/MoveAnnotation to Picked Atom.
  //@}
  void OnMoveAnnotationAtom();

  void OnUpdateExportModel(const MIUpdateEvent& pCmdUI);
  void OnExportModel();
  //@{
  // Update menu callback for Object/Edit Annotation
  //@}
  void OnUpdateEditAnnotation(const MIUpdateEvent& pCmdUI);
  //@{
  // Menu callback for Object/Edit Annotation
  //@}
  void OnEditAnnotation();
  void OnUpdateAnnotation(const MIUpdateEvent& pCmdUI);
  //@{
  // Create a new, empty model.
  //@}
  void OnNewModel();

  void OnUpdateMapReindex(const MIUpdateEvent& pCmdUI);
  void OnMapReindex();
  void OnUpdateFitSplitFit(const MIUpdateEvent& pCmdUI);
  void OnFitSplitFit();
  void OnUpdateFitSplitTorsion(const MIUpdateEvent& pCdUI);
  void OnFitSplitTorsion();
  void OnUpdateFitRange(const MIUpdateEvent& pCmdUI);
  void OnFitRange();
  void OnUpdateRefiResidue(const MIUpdateEvent& pCmdUI);
  void OnRefiResidue();
  void OnUpdateFullScreen(const MIUpdateEvent& pCmdUI);
  void OnFullScreen();
  void OnFindLigandDensity();
  void OnUpdateFindLigandDensity(const MIUpdateEvent& pCmdUI);
  void OnUpdateMapCenterVisibleDensity(const MIUpdateEvent& pCmdUI);
  void OnMapCenterVisibleDensity();
  void OnUpdateLabelEveryNth(const MIUpdateEvent& pCmdUI);
  void OnLabelEveryNth();
  void OnClearGeomAnnotations();
  void OnUpdateFindGeomErrors(const MIUpdateEvent& pCmdUI);
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
  void OnUpdateHideModel(const MIUpdateEvent& pCmdUI);

  void OnRefiLigandFit();
  void OnUpdateRefiLigandFit(const MIUpdateEvent &pCmdUI);

  /**
   * Update menu callback for Model/Add Waters.
   */
  void OnUpdateAddWater(const MIUpdateEvent& pCmdUI);
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
  void OnUpdateGotoCter(const MIUpdateEvent& pCmdUI);
  /**
   * Handler for the key '[' to center of the first residue in the chain or fragment.
   */
  void OnGotoNter();
  /**
   * Update menu Handler to center of the last residue in the chain or fragment.
   */
  void OnUpdateGotoNter(const MIUpdateEvent& pCmdUI);

  /**
   * Menu callback for Model/Build CB Range.
   */
  void OnBuildCBRange();
  /**
   * Update menu callback for Model/Checkpoint Model.
   */
  void OnUpdateCheckpointModel(const MIUpdateEvent& pCmdUI);
  /**
   * Menu callback for Model/Checkpoint Model.
   */
  void OnCheckpointModel();
  void OnUpdateAutoCheckpointModel(const MIUpdateEvent& pCmdUI);
  void OnAutoCheckpointModel();
  /**
   * Update menu callback for Model/Revert Model.
   */
  void OnUpdateRevertModel(const MIUpdateEvent& pCmdUI);
  /**
   * Menu callback for Model/Revert Model.
   */
  void OnRevertModel();
  /**
   * Update menu callback for Model/Replace with Sequence.
   */
  void OnUpdateReplaceSequence(const MIUpdateEvent& pCmdUI);
  /**
   * Menu callback for Model/Replace with Sequence.
   */
  void OnReplaceSequence();
  /**
   * Update menu callback for Sequence/Set Sequence Chain.
   */
  void OnUpdateSequencePositionChain(const MIUpdateEvent& pCmdUI);
  /**
   * Update menu callback for Sequence/Set Sequence Range.
   */
  void OnUpdateSequencePosition(const MIUpdateEvent& pCmdUI);
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
  void OnUpdateSequenceSaveModel(const MIUpdateEvent& pCmdUI);
  /**
   * Menu callback for Sequence/Save Model Sequence.
   */
  void OnSequenceSaveModel();
  void OnBuildMainchainRange();
  void OnUpdatePolyAlaChain(const MIUpdateEvent& pCmdUI);
  void OnPolyAlaChain();
  void OnUpdatePolyAla(const MIUpdateEvent& pCmdUI);
  void OnPolyAla();
  /**
   * Open a graph window with a phi-psi/Ramachandrn plot of the current model.
   */
#ifdef USE_QT_RAMAPLOT
  void OnRamachandranPlotShowAllowed();
  void OnUpdateRamachandranPlotShowAllowed(const MIUpdateEvent& pCmdUI);
#endif
  void OnUpdateFitSurfaceProbe(const MIUpdateEvent& pCmdUI);
  void OnUpdateFitSurfaceExtended(const MIUpdateEvent& pCmdUI);
  void OnUpdateFitSurfaceVdw(const MIUpdateEvent& pCmdUI);
  void OnUpdateFitSurfaceNone(const MIUpdateEvent& pCmdUI);
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
  void OnUpdateAddMarkBefore(const MIUpdateEvent& pCmdUI);
  /**
   * Menu callback for Model/Add MRK Before.
   */
  void OnAddMarkBefore();
  void OnUpdateAddMarkAfter(const MIUpdateEvent& pCmdUI);
  /**
   * Menu callback for Model/Add MRK After.
   */
  void OnAddMarkAfter();
  /**
   * Hide the model.
   */
  void OnHideModel();
  void OnUpdateFlipPeptide(const MIUpdateEvent& pCmdUI);
  /**
   * Menu callback for Fit/Fix Backbone/Flip Peptide.
   */
  void OnFlipPeptide();
  void OnUpdateMarkBeta(const MIUpdateEvent& pCmdUI);
  void OnUpdateMarkAlpha(const MIUpdateEvent& pCmdUI);
  /**
   * Set the refine target for this residue as beta sheet.
   */
  void OnMarkBeta();
  /**
   * Set the refine target for this residue as alpha helix.
   */
  void OnMarkAlpha();
  void OnUpdateReplaceAll(const MIUpdateEvent& pCmdUI);
  /**
   * Menu callback for Fit/Fix Backbone/Replace All 5.
   */
  void OnReplaceAll();
  void OnUpdateReplaceLast4(const MIUpdateEvent& pCmdUI);
  /**
   * Menu callback for Fit/Fix Backbone/Replace Last 4.
   */
  void OnReplaceLast4();
  void OnUpdateReplaceFirst4(const MIUpdateEvent& pCmdUI);
  /**
   * Menu callback for Fit/Fix Backbone/Replace First 4.
   */
  void OnReplaceFirst4();
  void OnUpdateReplaceMiddle3(const MIUpdateEvent& pCmdUI);
  /**
   * Menu callback for Fit/Fix Backbone/Replace Middle 3.
   */
  void OnReplaceMiddle3();
  void OnUpdateClearPentamer(const MIUpdateEvent& pCmdUI);
  /**
   * Menu callback for Fit/Fix Backbone/Clear Backbone Match.
   */
  void OnClearPentamer();
  void OnUpdateFindPentamer(const MIUpdateEvent& pCmdUI);
  /**
   * Menu callback for Fit/Fix Backbone/Suggest Backbone Match.
   */
  void OnFindPentamer();
  /**
   * Callback for Show/Ribbon Colors
   */
  void OnSetRibbonColors();
  void OnUpdateFitReplaceAndFit(const MIUpdateEvent& pCmdUI);
  /**
   * Callback for Model/Fit and Replace.
   */
  void OnFitReplaceAndFit();
  void OnUpdateRefiConfomer(const MIUpdateEvent& pCmdUI);
  /**
   * Callback for Refi/Find Best Conformer.
   */
  void OnRefiConformer();
  void OnUpdateRefiRigidBody(const MIUpdateEvent& pCmdUI);

  /**
   * Callback for Refi/Rigid Body Refine Current Atoms.
   */
  void OnRefiRigidBody();
  void OnUpdateFitRedo(const MIUpdateEvent& pCmdUI);
  void OnUpdateFitUndo(const MIUpdateEvent& pCmdUI);
  /**
   * Callback for Fit/Redo.  Goes back up the Undo Stack.
   */
  void OnFitRedo();
  /**
   * Callback for Fit/Undo.  Goes down up the Undo Stack.
   */
  void OnFitUndo();

  void OnUpdateRefiRedo(const MIUpdateEvent& pCmdUI);
  /**
   * Callback for Refi/Redo.
   */
  void OnRefiReDo();
  void OnUpdateRefiReset(const MIUpdateEvent& pCmdUI);
  void OnUpdateRefiCancel(const MIUpdateEvent& pCmdUI);
  void OnUpdateRefiAccept(const MIUpdateEvent& pCmdUI);
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
  void OnUpdateRenderLinesmooth(const MIUpdateEvent& p);
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
  void OnUpdateRefiMolecule(const MIUpdateEvent& pCmdUI);
  void OnUpdateRefiUndo(const MIUpdateEvent& pCmdUI);
  void OnUpdateRefiRange(const MIUpdateEvent& pCmdUI);
  void OnUpdateRefiRegion(const MIUpdateEvent& pCmdUI);
  void OnUpdateMapSwitch(const MIUpdateEvent& pCmdUI);
  /**
   * Callback for Map/Switch Map.  Cahnges teh current map to the one picked from a list.
   */
  void OnMapSwitch();
  void OnUpdateMapContourLevels(const MIUpdateEvent& pCmdUI);
  /**
   * Callback for Map/Contour.
   */
  void OnMapContourLevels();
  void OnUpdateMapFFT(const MIUpdateEvent& pCmdUI);
  /**
   * Callback for Map/FFT Phases.
   */
  void OnMapFFT();
  void OnUpdateMapSFCalc(const MIUpdateEvent& pCmdUI);
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
  void OnUpdateDecreasePersp(const MIUpdateEvent& pCmdUI);
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
  /**
   * Callback for Render/Line Thickness/1 pixel.
   */
  void OnLinethicknessOne();
  /**
   * Callback for Render/Line Thickness/2 pixel.
   */
  void OnLinethicknessTwo();
  /**
   * Callback for Render/Line Thickness/3 pixel.
   */
  void OnLinethicknessThree();
  /**
   * Callback for Render/Line Thickness/4 pixel.
   */
  void OnLinethicknessFour();
  void OnRenderingDepthcuedcolors();
  void OnRenderingDepthcuedlinewidth();

  void OnDimNonactiveModels();
  void OnAmountToDimNonactiveModels();
  void OnUpdateDimNonactiveModels(const MIUpdateEvent& pCmdUI);

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

  void OnUpdateAnimateRock(const MIUpdateEvent& pCmdUI);
  void OnUpdateAnimateRoll(const MIUpdateEvent& pCmdUI);
  void OnUpdateLinethicknessFour(const MIUpdateEvent& pCmdUI);
  void OnUpdateLinethicknessOne(const MIUpdateEvent& pCmdUI);
  void OnUpdateLinethicknessThree(const MIUpdateEvent& pCmdUI);
  void OnUpdateLinethicknessTwo(const MIUpdateEvent& pCmdUI);
  void OnUpdateRenderingDepthcuedcolors(const MIUpdateEvent& pCmdUI);
  void OnUpdateRenderingDepthcuedlinewidth(const MIUpdateEvent& pCmdUI);
  void OnUpdateHardwareStereo(const MIUpdateEvent& pCmdUI);
  void OnUpdateStereoToggle(const MIUpdateEvent& pCmdUI);
  void OnUpdateViewAtomstack(const MIUpdateEvent& pCmdUI);
  void OnUpdateViewGnomon(const MIUpdateEvent& pCmdUI);
  void OnUpdateViewUnitCell(const MIUpdateEvent& pCmdUI);
  void OnUpdateViewOrthonormal(const MIUpdateEvent& pCmdUI);
  /**
   * Callback for Render/Ball-and-cylinder.
   */
  void OnRenderingBallandstick();
  void OnUpdateRenderingBallandstick(const MIUpdateEvent& pCmdUI);
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
  void OnUpdateViewLabels(const MIUpdateEvent& pCmdUI);
  void OnVdwDotSurface();
  void OnViewClearmessage();
  // void OnObjectsModels();
  void OnFitResidue();
  void OnFitResidues();
  void OnUpdateFitResidues(const MIUpdateEvent& pCmdUI);
  void OnFitRotate();
  void OnFitTorsion();
  void OnUpdateFitTorsion(const MIUpdateEvent& pCmdUI);
  void OnFitTranslate();
  void OnUpdateFitTranslate(const MIUpdateEvent& pCmdUI);
  void OnFitCentermode();
  void OnUpdateFitCentermode(const MIUpdateEvent& pCmdUI);
  void OnFitCancel();
  void OnFitReset();
  void OnFitApply();
  void OnFitSingleatom();
  void OnFitAtoms();
  void OnUpdateFitAtoms(const MIUpdateEvent& pCmdUI);
  void OnFitSetuptorsion();
  void OnUpdateFitSetuptorsion(const MIUpdateEvent& pCmdUI);
  void OnFitCleartorsion();
  void OnFitNextConfomer();
  void OnUpdateFitNextConfomer(const MIUpdateEvent& pCmdUI);
  void OnUpdateFitCleartorsion(const MIUpdateEvent& pCmdUI);
  void OnUpdateFitApply(const MIUpdateEvent& pCmdUI);
  void OnUpdateFitCancel(const MIUpdateEvent& pCmdUI);
  void OnUpdateFitReset(const MIUpdateEvent& pCmdUI);
  void OnUpdateFitRotate(const MIUpdateEvent& pCmdUI);
  void OnFitReplacewith();
  void OnUpdateFitReplacewith(const MIUpdateEvent& pCmdUI);
  void OnObjectsAllatoms();
  void OnViewTopview();
  void OnUpdateViewTopview(const MIUpdateEvent& pCmdUI);
  void OnRenderSpacefilling();
  void OnUpdateRenderSpacefilling(const MIUpdateEvent& pCmdUI);
  void OnRenderSticks();
  void OnUpdateRenderSticks(const MIUpdateEvent& pCmdUI);
  void OnGeomBond();
  void OnUpdateGeomBond(const MIUpdateEvent& pCmdUI);
  void OnGeomUnbond();
  void OnUpdateGeomUnbond(const MIUpdateEvent& pCmdUI);
  void OnUpdateGeometryAngle(const MIUpdateEvent& pCmdUI);
  void OnUpdateGeometryDistance(const MIUpdateEvent& pCmdUI);
  void OnUpdateGeometryTorsion(const MIUpdateEvent& pCmdUI);
  void OnRenderBallsize();
  void OnUpdateRenderBallsize(const MIUpdateEvent& pCmdUI);
  void OnRenderBallandcylinder();
  void OnUpdateRenderBallandcylinder(const MIUpdateEvent& pCmdUI);
  void OnGotoZoomiin();
  void OnGotoZoomout();
  void OnGotoFittoscreen();
  void OnObjectClearstack();
  void OnEditCopy();
  void OnAnimateRockandrollparameters();
  void OnUpdateFitResidue(const MIUpdateEvent& pCmdUI);
  void OnUpdateFitSingleatom(const MIUpdateEvent& pCmdUI);
  void OnObjectShowresidue();
  void OnUpdateObjectShowresidue(const MIUpdateEvent& pCmdUI);
  void OnObjectShowsidechain();
  void OnUpdateObjectShowsidechain(const MIUpdateEvent& pCmdUI);
  void OnGeomNeighbours();
  void OnUpdateGeomNeighbours(const MIUpdateEvent& pCmdUI);
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
  void OnUpdateGeomAddsinglehbond(const MIUpdateEvent& pCmdUI);
  void OnViewContacts();
  void OnUpdateViewContacts(const MIUpdateEvent& pCmdUI);
  void OnFitLsqsuperpose();
  void OnFitFitmolecule();
  void OnUpdateFitFitmolecule(const MIUpdateEvent& pCmdUI);
  void OnFitDeleteresidue();
  void OnUpdateFitDeleteresidue(const MIUpdateEvent& pCmdUI);
  void OnFitInsertresidue();
  void OnUpdateFitInsertresidue(const MIUpdateEvent& pCmdUI);
  void OnFitRenameresidue();
  void OnUpdateFitRenameresidue(const MIUpdateEvent& pCmdUI);
  void OnGotoFitalltoscreen();
  void OnShowBackboneAsCATrace();
  void OnShowBackboneAsAtoms();
  void OnAnimateBlink();
  void OnUpdateAnimateBlink(const MIUpdateEvent& pCmdUI);
  void OnObjectSurfaceresidue();
  void OnUpdateObjectSurfaceresidue(const MIUpdateEvent& pCmdUI);
  void OnObjectSurfaceClearsurface();
  void OnObjectSurfaceresidues();
  void OnUpdateObjectSurfaceresidues(const MIUpdateEvent& pCmdUI);
  void OnUpdateObjectSurfaceatoms(const MIUpdateEvent& pCmdUI);
  void OnObjectSurfaceatom();
  void OnUpdateObjectSurfaceAtom(const MIUpdateEvent& pCmdUI);
  void OnObjectSurfaceAtoms();
  void OnObjectStackExpandtopallatomsinresidue();
  void OnUpdateObjectStackExpandtopallatomsinresidue(const MIUpdateEvent& pCmdUI);
  void OnUpdateObjectStackExpandtop2residues(const MIUpdateEvent& pCmdUI);
  void OnObjectStackExpandtop2residues();
  void OnObjectStackExpandtop2allatomsinrange();
  void OnUpdateObjectStackExpandtop2allatomsinrange(const MIUpdateEvent& pCmdUI);
  void OnObjectShowresidues();
  void OnUpdateObjectShowresidues(const MIUpdateEvent& pCmdUI);
  void OnObjectShowresiduerange();
  void OnObjectShowsidechainrange();
  void OnUpdateObjectShowresiduerange(const MIUpdateEvent& pCmdUI);
  void OnUpdateObjectShowsidechainrange(const MIUpdateEvent& pCmdUI);
  void OnObjectShowsidechains();
  void OnUpdateObjectShowsidechains(const MIUpdateEvent& pCmdUI);
  void OnMapLoadfromphsfile();
  void OnMapLoadfromfile();
  void OnMapContour();
  void OnUpdateMapContour(const MIUpdateEvent& pCmdUI);
  void OnShowShowwithinsphere();
  void OnObjectWhenshowncolor();
  void OnObjectRadiusisBall();
  void OnObjectRadiusisCpk();
  void OnObjectRadiusisCylinder();
  void OnUpdateObjectRadiusisBall(const MIUpdateEvent& pCmdUI);
  void OnUpdateObjectRadiusisCpk(const MIUpdateEvent& pCmdUI);
  void OnUpdateObjectRadiusisCylinder(const MIUpdateEvent& pCmdUI);
  void OnUpdateObjectWhenshowncolor(const MIUpdateEvent& pCmdUI);
  void OnObjectResiduerangeColor();
  void OnUpdateObjectResiduerangeColor(const MIUpdateEvent& pCmdUI);
  void OnObjectResiduerangeRadius();
  void OnObjectResidueColor();
  void OnObjectResidueRadius();
  void OnUpdateObjectResiduerangeRadius(const MIUpdateEvent& pCmdUI);
  void OnObjectResiduerangeTurnoff();
  void OnObjectResiduesColor();
  void OnUpdateObjectResiduesColor(const MIUpdateEvent& pCmdUI);
  void OnObjectResidueTurnoff();
  void OnObjectResiduesRadius();
  void OnUpdateObjectResiduerangeTurnoff(const MIUpdateEvent& pCmdUI);
  void OnUpdateObjectResiduesRadius(const MIUpdateEvent& pCmdUI);
  void OnObjectResiduesTurnoff();
  void OnUpdateObjectResiduesTurnoff(const MIUpdateEvent& pCmdUI);
  void OnObjectAtomColor();
  void OnObjectAtomRadius();
  void OnUpdateObjectAtomColor(const MIUpdateEvent& pCmdUI);
  void OnUpdateObjectAtomRadius(const MIUpdateEvent& pCmdUI);
  void OnObjectAtomsColor();
  void OnUpdateObjectAtomsColor(const MIUpdateEvent& pCmdUI);
  void OnObjectAtomsRadius();
  void OnUpdateObjectAtomsRadius(const MIUpdateEvent& pCmdUI);
  void OnShowColorallatoms();
  void OnUpdateObjectResidueColor(const MIUpdateEvent& pCmdUI);
  void OnUpdateObjectResidueRadius(const MIUpdateEvent& pCmdUI);
  void OnUpdateObjectResidueTurnoff(const MIUpdateEvent& pCmdUI);
  void OnShowUndocolorradius();
  void OnUpdateShowUndocolorradius(const MIUpdateEvent& pCmdUI);
  void OnObjectStackDeletetopitem();
  void OnViewUndo();
  void OnUpdateViewUndo(const MIUpdateEvent& pCmdUI);
  void OnShowPickedatomTurnon();
  void OnUpdateShowPickedatomTurnon(const MIUpdateEvent& pCmdUI);
  void OnShowPickedatomTurnoff();
  void OnUpdateShowPickedatomTurnoff(const MIUpdateEvent& pCmdUI);
  void OnShowAllpickedatomsTurnoff();
  void OnUpdateShowAllpickedatomsTurnoff(const MIUpdateEvent& pCmdUI);
  void OnShowAllpickedatomsTurnon();
  void OnUpdateShowAllpickedatomsTurnon(const MIUpdateEvent& pCmdUI);
  void OnShowRadiusmodel();
  void OnUpdateShowRadiusmodel(const MIUpdateEvent& pCmdUI);
  void OnUpdateShowColorallatoms(const MIUpdateEvent& pCmdUI);
  void OnUpdateSurfaceSolvent(const MIUpdateEvent& pCmdUI);
  void OnShowHideBackbone();
  void OnShowSidechainAtoms();
  void OnHideSidechainAtoms();
  void OnShowHidehydrogens();
  void OnUpdateShowHidehydrogens(const MIUpdateEvent& pCmdUI);
  void OnUpdateObjectSurfaceSpherearoundatom(const MIUpdateEvent& pCmdUI);
  void OnSurfaceSolvent();
  void OnObjectSurfaceSpherearoundatom();
  void OnSequenceEnter();
  void OnUpdateSequenceEnter(const MIUpdateEvent& pCmdUI);
  void OnSequenceRead();
  void OnUpdateSequenceRead(const MIUpdateEvent& pCmdUI);
  void OnSequenceSave();
  void OnUpdateSequenceSave(const MIUpdateEvent& pCmdUI);
  void OnSequenceInsertgap();
  void OnUpdateSequenceInsertgap(const MIUpdateEvent& pCmdUI);
  void OnSequenceDeletegap();
  void OnUpdateSequenceDeletegap(const MIUpdateEvent& pCmdUI);
  void OnSequenceInsertlowergap();
  void OnUpdateSequenceInsertlowergap(const MIUpdateEvent& pCmdUI);
  void OnSequenceDeletelowergap();
  void OnUpdateSequenceDeletelowergap(const MIUpdateEvent& pCmdUI);
#ifdef USE_ASPLOT
  void OnSitePlot();
  void OnUpdateSitePlot(const MIUpdateEvent& pCmdUI);
#endif

  void OnInvertChiralCenter();
  void OnUpdateInvertChiralCenter(const MIUpdateEvent& pCmdUI);
  void OnViewChiralCenters();
  void OnUpdateViewChiralCenters(const MIUpdateEvent& pCmdUI);

  void OnSolidSurfaceCommand(const MIActionEvent &evt);

  void OnExportImage();

  void OnDeleteAtom();
  void OnUpdateDeleteAtom(const MIUpdateEvent& evt);

  void OnShowSymmAtomsAsAtoms();
  void OnShowSymmAtomsAsCATrace();
  void OnShowHideSymmAtoms();
  void OnShowSaveSymmAtoms();
  void OnUpdateShowSaveSymmAtoms(const MIUpdateEvent& event);

};

inline bool MIGLWidget::IsFitting() {
  return CurrentAtoms.size() > 0;
}

inline bool MIGLWidget::IsTorsioning() {
  return ( Torsioning > 0 && CurrentAtoms.size() > 0 && fitmol);
}


inline ViewPoint* MIGLWidget::GetViewPoint() {
  return viewpoint;
}

inline int MIGLWidget::GetUpdate() {
  return Update;
}

inline void MIGLWidget::ResetUpdate() {
  Update = 0;
}

inline bool MIGLWidget::IsDrawing() {
  return is_drawing;
}


#endif
