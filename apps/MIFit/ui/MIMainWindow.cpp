#include <errno.h>
#include <sys/stat.h>
#include <sys/types.h>
#include <ctype.h> // for alanum
#include <QtGui>
#include <core/Version.h>
#include <chemlib/chemlib.h>
#include <chemlib/Residue.h>
#include <conflib/conflib.h>
#include <jobs/jobslib.h>
#include <util/utillib.h>
#include <core/corelib.h>
#include <ligand/ligandlib.h>
#include <map/maplib.h>
#include <nongui/nonguilib.h>

#include "JobsView.h"
#include "DictEditCanvas.h"
#include "DisplayView.h"
#include "EMap.h"
#include "Displaylist.h"
#include "ManageCrystals.h"
#include "ModelsView.h"
#include "MIMolIO.h"
#include "molw.h"
//#include "preferences/MapPreferencesPanel.h"
//#include "preferences/EnvironmentPreferencesPanel.h"
//#include "preferences/GeneralPreferencesPanel.h"
//#include "preferences/PreferencesDialog.h"
#include "PreferencesDialog.h"
#include "MIMainWindow.h"
#include "MIMenu.h"
#include "MIMenuBar.h"
#include "MIToolBar.h"
#include "MIGLWidget.h"
#include "Application.h"
#include "id.h"
#include "surf.h"
#include "MIMolIO.h"
#include "molw.h"
#include "GLOverviewCanvas.h"
#include "RamaPlot.h"
#include "DictEditDialog.h"
#include "MIEventHandlerMacros.h"
#include "GLFormatDialog.h"
#include "tools.h"
#include "GenericDataDialog.h"
#include "ui/SmilesDialog.h"
#include "ui/AtomColors.h"
#include "ui/BValueColors.h"
#include "CurrentMIGLWidgetAction.h"
#include "script/LocalSocketScript.h"

#ifdef _WIN32
#include <images/mifit_icon_32x32.xpm>
#else
#include <images/mifit_icon.xpm>
#endif

static bool GetOpenFilenames(std::vector<std::string> &fnames);
static void SplitPath(const std::string &origPath,
                      std::string *dirname,
                      std::string *fname,
                      std::string *ext);

using namespace chemlib;

template <typename T>
static T firstChildWithName(const QWidget *widget, const QString &name)
{
    return widget->findChildren<T>(name).first();
}

//  The sole purpose of this class is to allow us to dispatch menu/update
//  events to whichever mdi child is currently active, while still only
//  having one instance of each menu item, instead of a different copy of
//  each menu item for each MdiChild.
class ActiveMIGLWidgetFtor
    : public MIChildEventHandlerFtor
{
public:
    virtual MIChildEventHandlerFtor *CreateCopy()
    {
        return new ActiveMIGLWidgetFtor();
    }

    QObject*operator()()
    {
        return MIMainWindow::instance()->currentMIGLWidget();
    }
};

void MIMainWindow::initMenuHandlers()
{
    MIEventHandler *_dispatcher = this;
    ActiveMIGLWidgetFtor *ftor = new ActiveMIGLWidgetFtor();

//
    EVT_MENU(ID_EDIT_LABELS, MIGLWidget::OnEditLabels)

    EVT_MENU(ID_OBJECT_NEWMODEL, MIGLWidget::OnNewModel)

//geommenu
    EVT_MENU(ID_GEOMETRY_DISTANCE, MIGLWidget::OnGeometryDistance)
    EVT_UPDATE_UI(ID_GEOMETRY_DISTANCE, MIGLWidget::OnUpdateGeometryDistance)

    EVT_MENU(ID_GEOMETRY_ANGLE, MIGLWidget::OnGeometryAngle)
    EVT_UPDATE_UI(ID_GEOMETRY_ANGLE, MIGLWidget::OnUpdateGeometryAngle)

    EVT_MENU(ID_GEOMETRY_TORSION, MIGLWidget::OnGeometryTorsion)
    EVT_UPDATE_UI(ID_GEOMETRY_TORSION, MIGLWidget::OnUpdateGeometryTorsion)

    EVT_MENU(ID_GEOM_BOND, MIGLWidget::OnGeomBond)
    EVT_UPDATE_UI(ID_GEOM_BOND, MIGLWidget::OnUpdateGeomBond)

    EVT_MENU(ID_GEOM_UNBOND, MIGLWidget::OnGeomUnbond)
    EVT_UPDATE_UI(ID_GEOM_UNBOND, MIGLWidget::OnUpdateGeomUnbond)

    EVT_MENU(ID_GEOM_NEIGHBOURS, MIGLWidget::OnGeomNeighbours)
    EVT_UPDATE_UI(ID_GEOM_NEIGHBOURS, MIGLWidget::OnUpdateGeomNeighbours)

    EVT_MENU(ID_GEOM_CLEARNEIGHBOURS, MIGLWidget::OnGeomClearneighbours)
    EVT_MENU(ID_GEOM_HBONDS, MIGLWidget::OnGeomHbonds)
    EVT_MENU(ID_GEOM_CLEARHBONDS, MIGLWidget::OnGeomClearhbonds)

    EVT_MENU(ID_GEOM_ADDSINGLEHBOND, MIGLWidget::OnGeomAddsinglehbond)
    EVT_UPDATE_UI(ID_GEOM_ADDSINGLEHBOND, MIGLWidget::OnUpdateGeomAddsinglehbond)

    EVT_MENU(ID_GEOM_FINDGEOMERRORS, MIGLWidget::OnFindGeomErrors)
    EVT_UPDATE_UI(ID_GEOM_FINDGEOMERRORS, MIGLWidget::OnUpdateFindGeomErrors)

    EVT_MENU(ID_GEOM_CLEARGEOMERRORS, MIGLWidget::OnClearGeomAnnotations)

    //showmenu
    EVT_MENU(ID_OBJECTS_ALLATOMS, MIGLWidget::OnObjectsAllatoms)
    EVT_MENU(ID_OBJECT_BACKBONERIBBON, MIGLWidget::OnObjectBackboneribbon)
    EVT_MENU(ID_OBJECT_CLEARRIBBON, MIGLWidget::OnObjectClearribbon)

    //
    EVT_MENU(ID_RIBBONSECONDARYSTRUCTURE, MIGLWidget::OnRibbonSecondaryStructure)
    EVT_MENU(ID_SCHEMATICSECONDARYSTRUCTURE, MIGLWidget::OnSchematicSecondaryStructure)
    EVT_MENU(ID_DELETESECONDARYSTRUCTURE, MIGLWidget::OnDeleteSecondaryStructure)
    EVT_MENU(ID_SECONDARYSTRUCTUREOPTIONS_TUBE, MIGLWidget::OnSecondaryStructureOptionsTube)
    EVT_MENU(ID_SECONDARYSTRUCTUREOPTIONS_SHEET, MIGLWidget::OnSecondaryStructureOptionsSheet)
    EVT_MENU(ID_SECONDARYSTRUCTUREOPTIONS_TURN, MIGLWidget::OnSecondaryStructureOptionsTurn)
    EVT_MENU(ID_SECONDARYSTRUCTUREOPTIONS_RANDOM, MIGLWidget::OnSecondaryStructureOptionsRandom)
    EVT_MENU(ID_SECONDARYSTRUCTUREOPTIONS_HELIX, MIGLWidget::OnSecondaryStructureOptionsHelix)

//
    EVT_MENU(ID_OBJECTS_SECONDARYSTRUCTURE, MIGLWidget::OnSetRibbonColors)
    EVT_MENU(ID_SHOW_BACKBONECA, MIGLWidget::OnShowBackboneAsCATrace)
    EVT_MENU(ID_SHOW_BACKBONEATOMS, MIGLWidget::OnShowBackboneAsAtoms)
    EVT_MENU(ID_SHOW_HIDEBACKBONE, MIGLWidget::OnShowHideBackbone)
    EVT_MENU(ID_SHOW_SYMMATOMSASATOMS, MIGLWidget::OnShowSymmAtomsAsAtoms)
    EVT_MENU(ID_SHOW_SYMMATOMSASCA, MIGLWidget::OnShowSymmAtomsAsCATrace)

    EVT_MENU(ID_SHOW_SAVESYMMATOMS, MIGLWidget::OnShowSaveSymmAtoms)
    EVT_UPDATE_UI(ID_SHOW_SAVESYMMATOMS, MIGLWidget::OnUpdateShowSaveSymmAtoms)

    EVT_MENU(ID_SHOW_HIDESYMMATOMS, MIGLWidget::OnShowHideSymmAtoms)
    EVT_MENU(ID_SHOW_SIDECHAINATOMS, MIGLWidget::OnShowSidechainAtoms)
    EVT_MENU(ID_SHOW_HIDESIDECHAINATOMS, MIGLWidget::OnHideSidechainAtoms)

    EVT_MENU(ID_SHOW_HIDEHYDROGENS, MIGLWidget::OnShowHidehydrogens)
    EVT_UPDATE_UI(ID_SHOW_HIDEHYDROGENS, MIGLWidget::OnUpdateShowHidehydrogens)

    EVT_MENU(ID_SHOW_UNDOCOLORRADIUS, MIGLWidget::OnShowUndocolorradius)
    EVT_UPDATE_UI(ID_SHOW_UNDOCOLORRADIUS, MIGLWidget::OnUpdateShowUndocolorradius)

    EVT_MENU(ID_SHOW_LABELEVERYNTH, MIGLWidget::OnLabelEveryNth)
    EVT_UPDATE_UI(ID_SHOW_LABELEVERYNTH, MIGLWidget::OnUpdateLabelEveryNth)

    //viewmenu
    EVT_MENU(ID_VIEW_SLABIN, MIGLWidget::OnViewSlabin)
    EVT_MENU(ID_VIEW_SLABOUT, MIGLWidget::OnViewSlabout)

    EVT_MENU(ID_VIEW_ATOMSTACK, MIGLWidget::OnViewAtomstack)
    EVT_UPDATE_UI(ID_VIEW_ATOMSTACK, MIGLWidget::OnUpdateViewAtomstack)

    EVT_MENU(ID_VIEW_GNOMON, MIGLWidget::OnViewGnomon)
    EVT_UPDATE_UI(ID_VIEW_GNOMON, MIGLWidget::OnUpdateViewGnomon)

    EVT_MENU(ID_VIEW_UNITCELL, MIGLWidget::OnViewUnitCell)
    EVT_UPDATE_UI(ID_VIEW_UNITCELL, MIGLWidget::OnUpdateViewUnitCell)

    EVT_MENU(ID_VIEW_ROTATEY90, MIGLWidget::OnViewRotatey90)
    EVT_MENU(ID_VIEW_ROTATEYMINUS, MIGLWidget::OnViewRotateyminus)

    EVT_MENU(ID_VIEW_LABELS, MIGLWidget::OnViewLabels)
    EVT_UPDATE_UI(ID_VIEW_LABELS, MIGLWidget::OnUpdateViewLabels)

    EVT_MENU(ID_VIEW_CLEARMESSAGE, MIGLWidget::OnViewClearmessage)
    EVT_MENU(ID_GOTO_ZOOMIIN, MIGLWidget::OnGotoZoomiin)
    EVT_MENU(ID_GOTO_ZOOMOUT, MIGLWidget::OnGotoZoomout)

    EVT_MENU(ID_VIEW_TOPVIEW, MIGLWidget::OnViewTopview)
    EVT_UPDATE_UI(ID_VIEW_TOPVIEW, MIGLWidget::OnUpdateViewTopview)

    EVT_MENU(ID_GOTO_FITTOSCREEN, MIGLWidget::OnGotoFittoscreen)

    EVT_MENU(ID_VIEW_CONTACTS, MIGLWidget::OnViewContacts)
    EVT_UPDATE_UI(ID_VIEW_CONTACTS, MIGLWidget::OnUpdateViewContacts)

    EVT_MENU(ID_GOTO_FITALLTOSCREEN, MIGLWidget::OnGotoFitalltoscreen)
    EVT_MENU(ID_VIEW_UNDO, MIGLWidget::OnViewUndo)
    EVT_UPDATE_UI(ID_VIEW_UNDO, MIGLWidget::OnUpdateViewUndo)

    EVT_MENU(ID_INCREASE_PERSP, MIGLWidget::OnIncreasePersp)

    EVT_MENU(ID_DECREASE_PERSP, MIGLWidget::OnDecreasePersp)
    EVT_UPDATE_UI(ID_DECREASE_PERSP, MIGLWidget::OnUpdateDecreasePersp)

    EVT_MENU(ID_VIEW_ORTHONORMAL, MIGLWidget::OnViewOrthonormal)
    EVT_UPDATE_UI(ID_VIEW_ORTHONORMAL, MIGLWidget::OnUpdateViewOrthonormal)

    EVT_MENU(ID_VIEW_CLIPPLANES, MIGLWidget::OnViewClipplanes)
    EVT_MENU(ID_VIEW_SAVE, MIGLWidget::OnViewSave)
    EVT_MENU(ID_VIEW_LOAD, MIGLWidget::OnViewLoad)

    EVT_MENU(ID_VIEW_FULLSCREEN, MIGLWidget::OnFullScreen)
    EVT_UPDATE_UI(ID_VIEW_FULLSCREEN, MIGLWidget::OnUpdateFullScreen)

    //stackmenu
    EVT_MENU(ID_OBJECT_STACK_DELETETOPITEM, MIGLWidget::OnObjectStackDeletetopitem)

    EVT_MENU(ID_OBJECT_STACK_EXPANDTOPALLATOMSINRESIDUE, MIGLWidget::OnObjectStackExpandtopallatomsinresidue)
    EVT_UPDATE_UI(ID_OBJECT_STACK_EXPANDTOPALLATOMSINRESIDUE, MIGLWidget::OnUpdateObjectStackExpandtopallatomsinresidue)

    EVT_MENU(ID_OBJECT_CLEARSTACK, MIGLWidget::OnObjectClearstack)

    EVT_MENU(ID_OBJECT_STACK_EXPANDTOP2RESIDUES, MIGLWidget::OnObjectStackExpandtop2residues)
    EVT_UPDATE_UI(ID_OBJECT_STACK_EXPANDTOP2RESIDUES, MIGLWidget::OnUpdateObjectStackExpandtop2residues)

    EVT_MENU(ID_OBJECT_STACK_EXPANDTOP2ALLATOMSINRANGE, MIGLWidget::OnObjectStackExpandtop2allatomsinrange)
    EVT_UPDATE_UI(ID_OBJECT_STACK_EXPANDTOP2ALLATOMSINRANGE, MIGLWidget::OnUpdateObjectStackExpandtop2allatomsinrange)

//mapmenu
    EVT_MENU(ID_MAP_LOADMAP, MIGLWidget::OnMapLoadfromphsfile)
    EVT_MENU(ID_MAP_LOADMAPFILE, MIGLWidget::OnMapLoadfromfile)

    EVT_MENU(ID_MAP_CONTOUR, MIGLWidget::OnMapContour)
    EVT_UPDATE_UI(ID_MAP_CONTOUR, MIGLWidget::OnUpdateMapContour)

    EVT_MENU(ID_CONTOUR_LIST_FILE, MIGLWidget::OnContourListFile)

    EVT_MENU(ID_MAP_FFT, MIGLWidget::OnMapFFT)
    EVT_UPDATE_UI(ID_MAP_FFT, MIGLWidget::OnUpdateMapFFT)

    EVT_MENU(ID_MAP_SFCALC, MIGLWidget::OnMapSFCalc)
    EVT_UPDATE_UI(ID_MAP_SFCALC, MIGLWidget::OnUpdateMapSFCalc)

    EVT_MENU(ID_MAP_SWITCH, MIGLWidget::OnMapSwitch)
    EVT_UPDATE_UI(ID_MAP_SWITCH, MIGLWidget::OnUpdateMapSwitch)

    EVT_MENU(ID_MAP_CONTOURLEVELS, MIGLWidget::OnMapContourLevels)
    EVT_UPDATE_UI(ID_MAP_CONTOURLEVELS, MIGLWidget::OnUpdateMapContourLevels)

    EVT_MENU(ID_MAP_CENTERDENSITY, MIGLWidget::OnMapCenterVisibleDensity)
    EVT_UPDATE_UI(ID_MAP_CENTERDENSITY, MIGLWidget::OnUpdateMapCenterVisibleDensity)

    EVT_MENU(ID_FIND_DENSITY, MIGLWidget::OnFindLigandDensity)
    EVT_UPDATE_UI(ID_FIND_DENSITY, MIGLWidget::OnUpdateFindLigandDensity)

    EVT_MENU(ID_MAP_REINDEX, MIGLWidget::OnMapReindex)
    EVT_UPDATE_UI(ID_MAP_REINDEX, MIGLWidget::OnUpdateMapReindex)

    //fitmenuandmodelmenu
    EVT_MENU(ID_FIT_RESIDUE, MIGLWidget::OnFitResidue)
    EVT_UPDATE_UI(ID_FIT_RESIDUE, MIGLWidget::OnUpdateFitResidue)
    EVT_MENU(ID_FIT_RESIDUES, MIGLWidget::OnFitResidues)
    EVT_UPDATE_UI(ID_FIT_RESIDUES, MIGLWidget::OnUpdateFitResidues)
    EVT_MENU(ID_FIT_RANGE, MIGLWidget::OnFitRange)
    EVT_UPDATE_UI(ID_FIT_RANGE, MIGLWidget::OnUpdateFitRange)
    EVT_MENU(ID_FIT_SINGLEATOM, MIGLWidget::OnFitSingleatom)
    EVT_UPDATE_UI(ID_FIT_SINGLEATOM, MIGLWidget::OnUpdateFitSingleatom)
    EVT_MENU(ID_FIT_ATOMS, MIGLWidget::OnFitAtoms)
    EVT_UPDATE_UI(ID_FIT_ATOMS, MIGLWidget::OnUpdateFitAtoms)
    EVT_MENU(ID_FIT_FITMOLECULE, MIGLWidget::OnFitFitmolecule)
    EVT_UPDATE_UI(ID_FIT_FITMOLECULE, MIGLWidget::OnUpdateFitFitmolecule)
    EVT_MENU(ID_FIT_ROTATE, MIGLWidget::OnFitRotate)
    EVT_UPDATE_UI(ID_FIT_ROTATE, MIGLWidget::OnUpdateFitRotate)
    EVT_MENU(ID_FIT_TORSION, MIGLWidget::OnFitTorsion)
    EVT_UPDATE_UI(ID_FIT_TORSION, MIGLWidget::OnUpdateFitTorsion)
    EVT_MENU(ID_FIT_TRANSLATE, MIGLWidget::OnFitTranslate)
    EVT_UPDATE_UI(ID_FIT_TRANSLATE, MIGLWidget::OnUpdateFitTranslate)
    EVT_MENU(ID_FIT_CENTERMODE, MIGLWidget::OnFitCentermode)
    EVT_UPDATE_UI(ID_FIT_CENTERMODE, MIGLWidget::OnUpdateFitCentermode)
    EVT_MENU(ID_FIT_CANCEL, MIGLWidget::OnFitCancel)
    EVT_UPDATE_UI(ID_FIT_CANCEL, MIGLWidget::OnUpdateFitCancel)
    EVT_MENU(ID_FIT_RESET, MIGLWidget::OnFitReset)
    EVT_UPDATE_UI(ID_FIT_RESET, MIGLWidget::OnUpdateFitReset)
    EVT_MENU(ID_FIT_APPLY, MIGLWidget::OnFitApply)
    EVT_UPDATE_UI(ID_FIT_APPLY, MIGLWidget::OnUpdateFitApply)
    EVT_MENU(ID_FIT_UNDO, MIGLWidget::OnFitUndo)
    EVT_UPDATE_UI(ID_FIT_UNDO, MIGLWidget::OnUpdateFitUndo)
    EVT_MENU(ID_FIT_REDO, MIGLWidget::OnFitRedo)
    EVT_UPDATE_UI(ID_FIT_REDO, MIGLWidget::OnUpdateFitRedo)
    EVT_MENU(ID_FIT_SETUPTORSION, MIGLWidget::OnFitSetuptorsion)
    EVT_UPDATE_UI(ID_FIT_SETUPTORSION, MIGLWidget::OnUpdateFitSetuptorsion)
    EVT_MENU(ID_FIT_CLEARTORSION, MIGLWidget::OnFitCleartorsion)
    EVT_UPDATE_UI(ID_FIT_CLEARTORSION, MIGLWidget::OnUpdateFitCleartorsion)
    EVT_MENU(ID_FIT_REPLACEWITH, MIGLWidget::OnFitReplacewith)
    EVT_UPDATE_UI(ID_FIT_REPLACEWITH, MIGLWidget::OnUpdateFitReplacewith)
    EVT_MENU(ID_REPLACEANDFIT, MIGLWidget::OnFitReplaceAndFit)
    EVT_UPDATE_UI(ID_REPLACEANDFIT, MIGLWidget::OnUpdateFitReplaceAndFit)
    EVT_MENU(ID_FIT_DELETERESIDUE, MIGLWidget::OnFitDeleteresidue)
    EVT_UPDATE_UI(ID_FIT_DELETERESIDUE, MIGLWidget::OnUpdateFitDeleteresidue)
    EVT_MENU(ID_DELETEATOM, MIGLWidget::OnDeleteAtom)
    EVT_UPDATE_UI(ID_DELETEATOM, MIGLWidget::OnUpdateDeleteAtom)
    EVT_MENU(ID_FIT_RENAMERESIDUE, MIGLWidget::OnFitRenameresidue)
    EVT_UPDATE_UI(ID_FIT_RENAMERESIDUE, MIGLWidget::OnUpdateFitRenameresidue)
    EVT_MENU(ID_FIT_INSERTRESIDUE, MIGLWidget::OnFitInsertresidue)
    EVT_UPDATE_UI(ID_FIT_INSERTRESIDUE, MIGLWidget::OnUpdateFitInsertresidue)
    EVT_MENU(ID_FIT_LSQSUPERPOSE, MIGLWidget::OnFitLsqsuperpose)
    EVT_MENU(ID_FIT_NEXTCONFOMER, MIGLWidget::OnFitNextConfomer)
    EVT_UPDATE_UI(ID_FIT_NEXTCONFOMER, MIGLWidget::OnUpdateFitNextConfomer)

    EVT_MENU(ID_FIT_SURFVDW, MIGLWidget::OnFitSurfaceVdw)
    EVT_UPDATE_UI(ID_FIT_SURFVDW, MIGLWidget::OnUpdateFitSurfaceVdw)

    EVT_MENU(ID_FIT_SURFEXT, MIGLWidget::OnFitSurfaceExtended)
    EVT_UPDATE_UI(ID_FIT_SURFEXT, MIGLWidget::OnUpdateFitSurfaceExtended)

    EVT_MENU(ID_FIT_SURFPROBE, MIGLWidget::OnFitSurfaceProbe)
    EVT_UPDATE_UI(ID_FIT_SURFPROBE, MIGLWidget::OnUpdateFitSurfaceProbe)

    EVT_MENU(ID_FIT_SURFNONE, MIGLWidget::OnFitSurfaceNone)
    EVT_UPDATE_UI(ID_FIT_SURFNONE, MIGLWidget::OnUpdateFitSurfaceNone)

    EVT_MENU(ID_FIT_SPLITTORSION, MIGLWidget::OnFitSplitTorsion)
    EVT_UPDATE_UI(ID_FIT_SPLITTORSION, MIGLWidget::OnUpdateFitSplitTorsion)

    EVT_MENU(ID_FIT_SPLITFIT, MIGLWidget::OnFitSplitFit)
    EVT_UPDATE_UI(ID_FIT_SPLITFIT, MIGLWidget::OnUpdateFitSplitFit)

    //
    EVT_MENU(ID_FIT_POLYALA, MIGLWidget::OnPolyAla)
    EVT_UPDATE_UI(ID_FIT_POLYALA, MIGLWidget::OnUpdatePolyAla)

    EVT_MENU(ID_FIT_POLYALACHAIN, MIGLWidget::OnPolyAlaChain)
    EVT_UPDATE_UI(ID_FIT_POLYALACHAIN, MIGLWidget::OnUpdatePolyAlaChain)

    EVT_MENU(ID_FIT_MARKAFTER, MIGLWidget::OnAddMarkAfter)
    EVT_UPDATE_UI(ID_FIT_MARKAFTER, MIGLWidget::OnUpdateAddMarkAfter)

    EVT_MENU(ID_FIT_MARKBEFORE, MIGLWidget::OnAddMarkBefore)
    EVT_UPDATE_UI(ID_FIT_MARKBEFORE, MIGLWidget::OnUpdateAddMarkBefore)

    EVT_MENU(ID_FIT_MATCHPENTAMER, MIGLWidget::OnFindPentamer)
    EVT_UPDATE_UI(ID_FIT_MATCHPENTAMER, MIGLWidget::OnUpdateFindPentamer)

    EVT_MENU(ID_FIT_REPLACEFIRST4, MIGLWidget::OnReplaceFirst4)
    EVT_UPDATE_UI(ID_FIT_REPLACEFIRST4, MIGLWidget::OnUpdateReplaceFirst4)

    EVT_MENU(ID_FIT_REPLACEMIDDLE3, MIGLWidget::OnReplaceMiddle3)
    EVT_UPDATE_UI(ID_FIT_REPLACEMIDDLE3, MIGLWidget::OnUpdateReplaceMiddle3)

    EVT_MENU(ID_FIT_REPLACELAST4, MIGLWidget::OnReplaceLast4)
    EVT_UPDATE_UI(ID_FIT_REPLACELAST4, MIGLWidget::OnUpdateReplaceLast4)

    EVT_MENU(ID_FIT_REPLACEALL, MIGLWidget::OnReplaceAll)
    EVT_UPDATE_UI(ID_FIT_REPLACEALL, MIGLWidget::OnUpdateReplaceAll)

    EVT_MENU(ID_FIT_CLEARPENTAMER, MIGLWidget::OnClearPentamer)
    EVT_UPDATE_UI(ID_FIT_CLEARPENTAMER, MIGLWidget::OnUpdateClearPentamer)

    EVT_MENU(ID_FIT_MARK_ALPHA, MIGLWidget::OnMarkAlpha)
    EVT_UPDATE_UI(ID_FIT_MARK_ALPHA, MIGLWidget::OnUpdateMarkAlpha)

    EVT_MENU(ID_FIT_MARK_BETA, MIGLWidget::OnMarkBeta)
    EVT_UPDATE_UI(ID_FIT_MARK_BETA, MIGLWidget::OnUpdateMarkBeta)
    EVT_MENU(ID_FIT_FLIPPEPTIDE, MIGLWidget::OnFlipPeptide)
    EVT_UPDATE_UI(ID_FIT_FLIPPEPTIDE, MIGLWidget::OnUpdateFlipPeptide)
    EVT_MENU(ID_FIT_REPLACESEQUENCE, MIGLWidget::OnReplaceSequence)
    EVT_UPDATE_UI(ID_FIT_REPLACESEQUENCE, MIGLWidget::OnUpdateReplaceSequence)
    EVT_MENU(ID_MODEL_REVERT, MIGLWidget::OnRevertModel)
    EVT_UPDATE_UI(ID_MODEL_REVERT, MIGLWidget::OnUpdateRevertModel)
    EVT_MENU(ID_MODEL_AUTOCHECKPOINT, MIGLWidget::OnAutoCheckpointModel)
    EVT_UPDATE_UI(ID_MODEL_AUTOCHECKPOINT, MIGLWidget::OnUpdateAutoCheckpointModel)
    EVT_MENU(ID_MODEL_CHECKPOINT, MIGLWidget::OnCheckpointModel)
    EVT_UPDATE_UI(ID_MODEL_CHECKPOINT, MIGLWidget::OnUpdateCheckpointModel)
    EVT_MENU(ID_GOTO_NTER, MIGLWidget::OnGotoNter)
    EVT_UPDATE_UI(ID_GOTO_NTER, MIGLWidget::OnUpdateGotoNter)
    EVT_MENU(ID_GOTO_CTER, MIGLWidget::OnGotoCter)
    EVT_UPDATE_UI(ID_GOTO_CTER, MIGLWidget::OnUpdateGotoCter)
    EVT_MENU(ID_FIT_ADDWATER, MIGLWidget::OnAddWater)
    EVT_UPDATE_UI(ID_FIT_ADDWATER, MIGLWidget::OnUpdateAddWater)

//animemenu
    EVT_MENU(ID_ANIMATE_ROLL, MIGLWidget::OnAnimateRoll)
    EVT_UPDATE_UI(ID_ANIMATE_ROLL, MIGLWidget::OnUpdateAnimateRoll)
    EVT_MENU(ID_ANIMATE_ROCK, MIGLWidget::OnAnimateRock)
    EVT_UPDATE_UI(ID_ANIMATE_ROCK, MIGLWidget::OnUpdateAnimateRock)
    EVT_MENU(ID_ANIMATE_ROCKANDROLLPARAMETERS, MIGLWidget::OnAnimateRockandrollparameters)
    EVT_MENU(ID_ANIMATE_BLINK, MIGLWidget::OnAnimateBlink)
    EVT_UPDATE_UI(ID_ANIMATE_BLINK, MIGLWidget::OnUpdateAnimateBlink)

    EVT_MENU(ID_SEQU_SETCHAIN, MIGLWidget::OnSequencePositionChain)
    EVT_UPDATE_UI(ID_SEQU_SETCHAIN, MIGLWidget::OnUpdateSequencePositionChain)
    EVT_MENU(ID_SEQU_SETPOSITION, MIGLWidget::OnSequencePosition)
    EVT_UPDATE_UI(ID_SEQU_SETPOSITION, MIGLWidget::OnUpdateSequencePosition)
    EVT_MENU(ID_SEQUENCE_ENTER, MIGLWidget::OnSequenceEnter)
    EVT_UPDATE_UI(ID_SEQUENCE_ENTER, MIGLWidget::OnUpdateSequenceEnter)
    EVT_MENU(ID_SEQUENCE_READ, MIGLWidget::OnSequenceRead)
    EVT_UPDATE_UI(ID_SEQUENCE_READ, MIGLWidget::OnUpdateSequenceRead)
    EVT_MENU(ID_SEQUENCE_SAVE, MIGLWidget::OnSequenceSave)
    EVT_UPDATE_UI(ID_SEQUENCE_SAVE, MIGLWidget::OnUpdateSequenceSave)
    EVT_MENU(ID_SEQU_SAVEPROTEIN, MIGLWidget::OnSequenceSaveModel)
    EVT_UPDATE_UI(ID_SEQU_SAVEPROTEIN, MIGLWidget::OnUpdateSequenceSaveModel)
    EVT_MENU(ID_SEQUENCE_INSERTGAP, MIGLWidget::OnSequenceInsertgap)
    EVT_UPDATE_UI(ID_SEQUENCE_INSERTGAP, MIGLWidget::OnUpdateSequenceInsertgap)
    EVT_MENU(ID_SEQUENCE_DELETEGAP, MIGLWidget::OnSequenceDeletegap)
    EVT_UPDATE_UI(ID_SEQUENCE_DELETEGAP, MIGLWidget::OnUpdateSequenceDeletegap)
    EVT_MENU(ID_SEQUENCE_INSERTLOWERGAP, MIGLWidget::OnSequenceInsertlowergap)
    EVT_UPDATE_UI(ID_SEQUENCE_INSERTLOWERGAP, MIGLWidget::OnUpdateSequenceInsertlowergap)
    EVT_MENU(ID_SEQUENCE_DELETELOWERGAP, MIGLWidget::OnSequenceDeletelowergap)
    EVT_UPDATE_UI(ID_SEQUENCE_DELETELOWERGAP, MIGLWidget::OnUpdateSequenceDeletelowergap)
    EVT_MENU(ID_OBJECT_RADIUSIS_BALL, MIGLWidget::OnObjectRadiusisBall)
    EVT_MENU(ID_OBJECT_RADIUSIS_CPK, MIGLWidget::OnObjectRadiusisCpk)
    EVT_MENU(ID_OBJECT_RADIUSIS_CYLINDER, MIGLWidget::OnObjectRadiusisCylinder)
    EVT_UPDATE_UI(ID_OBJECT_RADIUSIS_BALL, MIGLWidget::OnUpdateObjectRadiusisBall)
    EVT_UPDATE_UI(ID_OBJECT_RADIUSIS_CPK, MIGLWidget::OnUpdateObjectRadiusisCpk)
    EVT_UPDATE_UI(ID_OBJECT_RADIUSIS_CYLINDER, MIGLWidget::OnUpdateObjectRadiusisCylinder)
//
    EVT_MENU(ID_OBJECT_RESIDUERANGE_COLOR, MIGLWidget::OnObjectResiduerangeColor)
    EVT_UPDATE_UI(ID_OBJECT_RESIDUERANGE_COLOR, MIGLWidget::OnUpdateObjectResiduerangeColor)
    EVT_MENU(ID_OBJECT_RESIDUERANGE_RADIUS, MIGLWidget::OnObjectResiduerangeRadius)
    EVT_UPDATE_UI(ID_OBJECT_RESIDUERANGE_RADIUS, MIGLWidget::OnUpdateObjectResiduerangeRadius)
    EVT_MENU(ID_OBJECT_RESIDUERANGE_TURNOFF, MIGLWidget::OnObjectResiduerangeTurnoff)
    EVT_UPDATE_UI(ID_OBJECT_RESIDUERANGE_TURNOFF, MIGLWidget::OnUpdateObjectResiduerangeTurnoff)
    EVT_MENU(ID_OBJECT_RESIDUES_COLOR, MIGLWidget::OnObjectResiduesColor)
    EVT_UPDATE_UI(ID_OBJECT_RESIDUES_COLOR, MIGLWidget::OnUpdateObjectResiduesColor)
    EVT_MENU(ID_OBJECT_RESIDUES_RADIUS, MIGLWidget::OnObjectResiduesRadius)
    EVT_UPDATE_UI(ID_OBJECT_RESIDUES_RADIUS, MIGLWidget::OnUpdateObjectResiduesRadius)
    EVT_MENU(ID_OBJECT_RESIDUES_TURNOFF, MIGLWidget::OnObjectResiduesTurnoff)
    EVT_UPDATE_UI(ID_OBJECT_RESIDUES_TURNOFF, MIGLWidget::OnUpdateObjectResiduesTurnoff)
    EVT_MENU(ID_HIDEMODEL, MIGLWidget::OnHideModel)
    EVT_UPDATE_UI(ID_HIDEMODEL, MIGLWidget::OnUpdateHideModel)
    EVT_MENU(ID_OBJECT_ATOM_COLOR, MIGLWidget::OnObjectAtomColor)
    EVT_UPDATE_UI(ID_OBJECT_ATOM_COLOR, MIGLWidget::OnUpdateObjectAtomColor)
    EVT_MENU(ID_OBJECT_ATOM_RADIUS, MIGLWidget::OnObjectAtomRadius)
    EVT_UPDATE_UI(ID_OBJECT_ATOMS_COLOR, MIGLWidget::OnUpdateObjectAtomsColor)
    EVT_UPDATE_UI(ID_OBJECT_ATOM_RADIUS, MIGLWidget::OnUpdateObjectAtomRadius)
    EVT_MENU(ID_OBJECT_ATOMS_COLOR, MIGLWidget::OnObjectAtomsColor)
    EVT_MENU(ID_OBJECT_ATOMS_RADIUS, MIGLWidget::OnObjectAtomsRadius)
    EVT_UPDATE_UI(ID_OBJECT_ATOMS_RADIUS, MIGLWidget::OnUpdateObjectAtomsRadius)
    EVT_MENU(ID_OBJECT_SHOWRESIDUE, MIGLWidget::OnObjectShowresidue)
    EVT_UPDATE_UI(ID_OBJECT_SHOWRESIDUE, MIGLWidget::OnUpdateObjectShowresidue)
    EVT_MENU(ID_OBJECT_SHOWSIDECHAIN, MIGLWidget::OnObjectShowsidechain)
    EVT_UPDATE_UI(ID_OBJECT_SHOWSIDECHAIN, MIGLWidget::OnUpdateObjectShowsidechain)
    EVT_MENU(ID_OBJECT_SHOWRESIDUES, MIGLWidget::OnObjectShowresidues)
    EVT_UPDATE_UI(ID_OBJECT_SHOWRESIDUES, MIGLWidget::OnUpdateObjectShowresidues)
    EVT_MENU(ID_OBJECT_SHOWSIDECHAINS, MIGLWidget::OnObjectShowsidechains)
    EVT_UPDATE_UI(ID_OBJECT_SHOWSIDECHAINS, MIGLWidget::OnUpdateObjectShowsidechains)
    EVT_MENU(ID_OBJECT_SHOWRESIDUERANGE, MIGLWidget::OnObjectShowresiduerange)
    EVT_UPDATE_UI(ID_OBJECT_SHOWRESIDUERANGE, MIGLWidget::OnUpdateObjectShowresiduerange)
    EVT_UPDATE_UI(ID_OBJECT_SHOWSIDECHAINRANGE, MIGLWidget::OnUpdateObjectShowsidechainrange)
    EVT_MENU(ID_OBJECT_SHOWSIDECHAINRANGE, MIGLWidget::OnObjectShowsidechainrange)
//
    EVT_MENU(ID_SHOW_UNDOCOLORRADIUS, MIGLWidget::OnShowUndocolorradius)
    EVT_UPDATE_UI(ID_SHOW_UNDOCOLORRADIUS, MIGLWidget::OnUpdateShowUndocolorradius)
//
    EVT_MENU(ID_SHOW_PICKEDATOM_TURNON, MIGLWidget::OnShowPickedatomTurnon)
    EVT_UPDATE_UI(ID_SHOW_PICKEDATOM_TURNON, MIGLWidget::OnUpdateShowPickedatomTurnon)
    EVT_MENU(ID_SHOW_PICKEDATOM_TURNOFF, MIGLWidget::OnShowPickedatomTurnoff)
    EVT_UPDATE_UI(ID_SHOW_PICKEDATOM_TURNOFF, MIGLWidget::OnUpdateShowPickedatomTurnoff)
    EVT_MENU(ID_SHOW_ALLPICKEDATOMS_TURNOFF, MIGLWidget::OnShowAllpickedatomsTurnoff)
    EVT_UPDATE_UI(ID_SHOW_ALLPICKEDATOMS_TURNOFF, MIGLWidget::OnUpdateShowAllpickedatomsTurnoff)
    EVT_MENU(ID_SHOW_ALLPICKEDATOMS_TURNON, MIGLWidget::OnShowAllpickedatomsTurnon)
    EVT_UPDATE_UI(ID_SHOW_ALLPICKEDATOMS_TURNON, MIGLWidget::OnUpdateShowAllpickedatomsTurnon)
    EVT_MENU(ID_SHOW_COLORALLATOMS, MIGLWidget::OnShowColorallatoms)
    EVT_UPDATE_UI(ID_SHOW_COLORALLATOMS, MIGLWidget::OnUpdateShowColorallatoms)
    EVT_MENU(ID_SHOW_RADIUSMODEL, MIGLWidget::OnShowRadiusmodel)
    EVT_UPDATE_UI(ID_SHOW_RADIUSMODEL, MIGLWidget::OnUpdateShowRadiusmodel)
//
    EVT_MENU(ID_OBJECT_RESIDUE_COLOR, MIGLWidget::OnObjectResidueColor)
    EVT_UPDATE_UI(ID_OBJECT_RESIDUE_COLOR, MIGLWidget::OnUpdateObjectResidueColor)
    EVT_MENU(ID_OBJECT_RESIDUE_RADIUS, MIGLWidget::OnObjectResidueRadius)
    EVT_UPDATE_UI(ID_OBJECT_RESIDUE_RADIUS, MIGLWidget::OnUpdateObjectResidueRadius)
    EVT_MENU(ID_OBJECT_RESIDUE_TURNOFF, MIGLWidget::OnObjectResidueTurnoff)
    EVT_UPDATE_UI(ID_OBJECT_RESIDUE_TURNOFF, MIGLWidget::OnUpdateObjectResidueTurnoff)
//
    EVT_MENU(ID_RAMAPLOT_ALLOWED, MIGLWidget::OnRamachandranPlotShowAllowed)
    EVT_UPDATE_UI(ID_RAMAPLOT_ALLOWED, MIGLWidget::OnUpdateRamachandranPlotShowAllowed)
    EVT_MENU(ID_INVERTCHIRAL, MIGLWidget::OnInvertChiralCenter)
    EVT_UPDATE_UI(ID_INVERTCHIRAL, MIGLWidget::OnUpdateInvertChiralCenter)
    EVT_MENU(ID_VIEWCHIRAL, MIGLWidget::OnViewChiralCenters)
    EVT_UPDATE_UI(ID_VIEWCHIRAL, MIGLWidget::OnUpdateViewChiralCenters)
//
    EVT_MENU(ID_SHOW_SHOWWITHINSPHERE, MIGLWidget::OnShowShowwithinsphere)
    EVT_MENU(ID_GOTO_GOTOXYZ, MIGLWidget::OnGotoGotoxyz)
//todo
    EVT_MENU(ID_EDIT_SELECTMODEL, MIGLWidget::OnEditSelectmodel)

    delete ftor;
}

MIMainWindow*MIMainWindow::_instance = NULL;

MIMainWindow::MIMainWindow()
    : MIEventHandler(this),
      modelsDock(NULL),
      displayDock(NULL),
      jobsDock(NULL),
      logDock(NULL),
      ramaDock(NULL),
      logWindow(NULL),
      solidSurfMenuAction(NULL),
      _scriptSocket(NULL)
{

    _instance = this;
    QApplication::setWindowIcon(QIcon(QPixmap(mifit_icon)));

    mdiArea = new QMdiArea;
    setCentralWidget(mdiArea);
    connect(mdiArea, SIGNAL(subWindowActivated(QMdiSubWindow*)),
            this, SLOT(updateMenus()));
    connect(mdiArea, SIGNAL(subWindowActivated(QMdiSubWindow*)),
            this, SLOT(childActivated(QMdiSubWindow*)));
    windowMapper = new QSignalMapper(this);
    connect(windowMapper, SIGNAL(mapped(QWidget*)),
            this, SLOT(setActiveSubWindow(QWidget*)));

    createActions();
    createMenus();
    createDockWindows();

    createToolBars();
    createStatusBar();
    updateMenus();

    setWindowTitle(tr("MIFit"));

    LeftFooter("MIFit: Ok", 2500);

    middleFooter = new QLabel(statusBar());
    statusBar()->addPermanentWidget(middleFooter, 0);

    rightFooter = new QLabel(statusBar());
    statusBar()->addPermanentWidget(rightFooter, 0);

    navigator = new GLOverviewCanvas(this);
    navigatorDock = AddAsDockWidget(navigator, "Navigator", Qt::BottomDockWidgetArea);

    memset(cursors, 0, imhCount*sizeof(QCursor*));

    QSettings *settings = MIGetQSettings(); // could use MIConfig (it's the same file), but this api is easier here
    restoreState(settings->value("WindowLayout", saveState()).toByteArray());
    setGeometry(settings->value("WindowGeometry", geometry()).toRect());

    tabifyDockWidget(modelsDock, displayDock);
    tabifyDockWidget(displayDock, jobsDock);
    modelsDock->raise();

    tabifyDockWidget(navigatorDock, ramaDock);
    navigatorDock->raise();

    std::string title = "MIFit ";
    title += MIFit_version;
    setWindowTitle(title.c_str());

    std::string s = ::format("Read in %d scattering factors",
                             MIMapInitializeScatteringFactorTables(Application::instance()->CrystalData.c_str(),
                                                                   Application::instance()->MolimageHome.c_str()));
    Log(s);
}


MIMainWindow::~MIMainWindow()
{
    for (int i = 0; i < imhCount; ++i)
        delete cursors[i];
    MIMapFreeScatteringFactorTables();
    delete JobManager;
}

static QWidget *dying_widget = 0;

void MIMainWindow::OnClose()
{
    dying_widget = currentMIGLWidget();
    mdiArea->closeActiveSubWindow();
    dying_widget = 0;
}

void MIMainWindow::childActivated(QMdiSubWindow *w)
{
    if (!w)
        return;
    MIGLWidget *child = currentMIGLWidget();
    if (child)
    {
        if (child==dying_widget)
            return;
        child->OnActivated();
        if (tool_bar)
            tool_bar->doUpdates();
    }
}

void MIMainWindow::saveLayout()
{
    QSettings *settings = MIGetQSettings(); // could use MIConfig (it's the same file), but this api is easier here
    settings->setValue("WindowLayout", saveState());
    settings->setValue("WindowGeometry", geometry());
}

void MIMainWindow::OnExit()
{
    saveLayout();
    qApp->closeAllWindows();
}

void MIMainWindow::OnAbout()
{
    std::string s = ::format("MIFit %s\n", MIFit_version);
    QMessageBox::about(this, "About MIFit", s.c_str());
}

void MIMainWindow::Log(const std::string &msg)
{
    if (logWindow != NULL)
    {
        logWindow->append(msg.c_str());
    }
}

void MIMainWindow::OnManageCrystals()
{
    ManageCrystals mcd(this);
    mcd.exec();
}

static QFileDialog *fileDialog = NULL;

static void initializeFileDialog()
{
    if (fileDialog == NULL)
    {
        fileDialog = new QFileDialog(0, "Open File(s)", "",
                                     "Recognized files (*.mlw *.pdb *.mtz *.phs *.fcf *.cif *.map *.py);;"
                                     "MIFit session files (*.mlw);;"
                                     "PDB files (*.pdb);;"
                                     "MTZ files (*.mtz);;"
                                     "Phase files (*.phs);;"
                                     "FCF files (*.fcf);;"
//        "Scalepack files (*.sca);;"
                                     "CIF files (*.cif);;"
//        "Reflection files (*.ref);;"
                                     "Map files (*.map);;"
                                     "Python scripts (*.py);;"
                                     "All files (*.*)");
        fileDialog->setFileMode(QFileDialog::ExistingFiles);
    }
}

void MIMainWindow::OnFileOpen()
{
    fileOpen();
}

void MIMainWindow::fileOpen()
{
    std::vector<std::string> fnames;
    if (!GetOpenFilenames(fnames))
        return;
    OpenFiles(fnames);
}


void MIMainWindow::OnFileOpenNew()
{
    std::vector<std::string> fnames;
    if (!GetOpenFilenames(fnames))
        return;
    newWindow();
    OpenFiles(fnames);
}

// unused in Qt version, currently preserved to facilitate merging
//
// CMolwDoc* MainFrame::newCMolwDoc(const wxString& path) {
//   wxDocTemplate* docTemplate = NULL;
//   wxNode* node = m_docManager->GetTemplates().GetFirst();
//   while (node) {
//     wxDocTemplate* t = (wxDocTemplate*) node->GetData();
//     if (t->GetDefaultExtension() == "mlw") {
//       docTemplate = t;
//       break;
//     }
//     node = node->GetNext();
//   }
//   if (docTemplate == NULL) {
//     return NULL;
//   }
//   CMolwDoc* doc = (CMolwDoc*) docTemplate->CreateDocument(path);
//   if (doc != NULL) {
//     doc->SetDocumentName(docTemplate->GetDocumentName());
//     doc->SetDocumentTemplate(docTemplate);
//     doc->Modify(false);
//     if (!doc->OnNewDocument()) {
//       doc->DeleteAllViews();
//     }
//   }
//   return doc;
// }

// void MainFrame::importFiles(const wxArrayString& files) {
//   if (files.GetCount() > 0) {
//     wxString pathname = files[0];

//     CMolwDoc* doc = newCMolwDoc(pathname);
//     if (!doc) {
//       return;
//     }
//     if (!doc->OnOpenDocument(pathname)) {
//       return;
//     }

//     // default name is made by changing extension to .mlw
//     wxString docname(pathname);
//     wxString path, name, ext;
//     wxSplitPath(docname, &path, &name, &ext);
//     name += ".mlw";

//     m_docManager->AddFileToHistory(pathname);
//     doc->SetFilename(name, true);
//     doc->SetTitle(wxFileNameFromPath(name));
//     for (unsigned int i = 1; i < files.GetCount(); ++i) {
//       pathname = files[i];
//       doc->LoadPDBFile(pathname);
//     }
//     doc->Modify(true);
//   }
// }

// void MainFrame::OnTimer(wxTimerEvent&) {
//   timer();
// }

// void MainFrame::timer() {
//   static unsigned int ndocs = 0;
//   wxList& m = m_docManager->GetDocuments();
//   if (m.GetCount() != ndocs) {
//     ndocs = m.GetCount();
//     if (ndocs == 0) {
//       EnableToolBar(0);
//     }
//   }
//   if (MIBusyManager::instance()->Busy()) {
//     MIBusyManager::instance()->SetWaitCursor();
//   }
// }


void MIMainWindow::OnGLFormat()
{
    GLFormatDialog dialog;
    if (currentMIGLWidget())
    {
        dialog.setCurrentFormat(currentMIGLWidget()->format());
    }
    dialog.exec();
}

void MIMainWindow::OnHelp()
{
    QUrl url("http://code.google.com/p/mifit");
    QDesktopServices::openUrl(url);
}

void MIMainWindow::OnPreferences()
{
    showPreferences();
}

void MIMainWindow::showPreferences(int page)
{
    PreferencesDialog dlg(this);
    dlg.setPage(page);
    if (dlg.exec() != QDialog::Accepted)
    {
        return;
    }
    dlg.savePreferences();
    if (Application::instance()->SetDictionary(Application::instance()->XFitDictSetting))
    {
        Application::instance()->LoadDictionary();
    }
    Application::instance()->SetEnv();
    Application::instance()->Write();
}

void MIMainWindow::OnDefineAtomColors()
{
    MIData data;
    AtomColors dlg(this);
    dlg.setWindowTitle("Define Atom Coloring Scheme");
    dlg.InitializeFromData(data);
    if (dlg.exec() != QDialog::Accepted)
    {
        return;
    }
    dlg.GetData(data);
    //  input: atomNames.strList, atomColors.strList
    data["atomNames"].strList = Colors::atomnames;
    data["atomColors"].strList = Colors::atomcolors;

    Colors::atomnames = data["atomNames"].strList;
    Colors::atomcolors = data["atomColors"].strList;

    MIGLWidget *mdoc = currentMIGLWidget();
    if (mdoc)
    {
        Molecule *model = mdoc->GetDisplaylist()->CurrentItem();
        mdoc->ColorModel(model);
        mdoc->ReDraw();
    }
}

void MIMainWindow::OnDefineBValueColors()
{

    MIData data;
    data["color1"].i = Colors::BValueColors[0];
    data["color2"].i = Colors::BValueColors[1];
    data["color3"].i = Colors::BValueColors[2];
    data["color4"].i = Colors::BValueColors[3];
    data["color5"].i = Colors::BValueColors[4];
    data["color6"].i = Colors::BValueColors[5];
    data["color7"].i = Colors::BValueColors[6];
    data["color8"].i = Colors::BValueColors[7];
    data["color9"].i = Colors::BValueColors[8];
    data["color10"].i = Colors::BValueColors[9];

    data["level1"].f = Colors::BValueRanges[0]/100.0f;
    data["level2"].f = Colors::BValueRanges[1]/100.0f;
    data["level3"].f = Colors::BValueRanges[2]/100.0f;
    data["level4"].f = Colors::BValueRanges[3]/100.0f;
    data["level5"].f = Colors::BValueRanges[4]/100.0f;
    data["level6"].f = Colors::BValueRanges[5]/100.0f;
    data["level7"].f = Colors::BValueRanges[6]/100.0f;
    data["level8"].f = Colors::BValueRanges[7]/100.0f;
    data["level9"].f = Colors::BValueRanges[8]/100.0f;
    data["level10"].f = Colors::BValueRanges[9]/100.0f;

    BValueColors dlg(this);
    dlg.setWindowTitle("Define Atom Coloring Scheme");
    dlg.InitializeFromData(data);
    if (dlg.exec() != QDialog::Accepted)
    {
        return;
    }
    dlg.GetData(data);

    Colors::BValueColors[0] = data["color1"].i;
    Colors::BValueColors[1] = data["color2"].i;
    Colors::BValueColors[2] = data["color3"].i;
    Colors::BValueColors[3] = data["color4"].i;
    Colors::BValueColors[4] = data["color5"].i;
    Colors::BValueColors[5] = data["color6"].i;
    Colors::BValueColors[6] = data["color7"].i;
    Colors::BValueColors[7] = data["color8"].i;
    Colors::BValueColors[8] = data["color9"].i;
    Colors::BValueColors[9] = data["color10"].i;

    Colors::BValueRanges[0] = static_cast<int>(data["level1"].f*100.0f);
    Colors::BValueRanges[1] = static_cast<int>(data["level2"].f*100.0f);
    Colors::BValueRanges[2] = static_cast<int>(data["level3"].f*100.0f);
    Colors::BValueRanges[3] = static_cast<int>(data["level4"].f*100.0f);
    Colors::BValueRanges[4] = static_cast<int>(data["level5"].f*100.0f);
    Colors::BValueRanges[5] = static_cast<int>(data["level6"].f*100.0f);
    Colors::BValueRanges[6] = static_cast<int>(data["level7"].f*100.0f);
    Colors::BValueRanges[7] = static_cast<int>(data["level8"].f*100.0f);
    Colors::BValueRanges[8] = static_cast<int>(data["level9"].f*100.0f);
    Colors::BValueRanges[9] = static_cast<int>(data["level10"].f*100.0f);

    if (data["save"].b)
    {
        Application::instance()->Write();
    }
}


/// NOTE: after this point, MIMainWindow, and MainFrame are roughly in sync

void MIMainWindow::OnLoadDictReplace()
{
    QString s = QFileDialog::getOpenFileName(this, "Choose a dictionary file", "",
                                          "Dictionary files (*.pdb,*.cif,*.ent);;All files (*.*)");
    if (!s.isEmpty())
    {
        MIFitDictionary()->LoadDictionary(s.toAscii().constData(), false, true);
    }
}

void MIMainWindow::ShowDictEditor(const char *type)
{

    // note this dialog is dynamcially allocated, and is deleted when closed.
    // this is so that history can function
    std::string s = ::format("Dictionary Editor - %s", type);
    DictEditDialog *dlg = new DictEditDialog(this);
    dlg->setWindowTitle(s.c_str());
    dlg->resize(800, 720);
    dlg->setModal(true);
    dlg->exec();
}

void MIMainWindow::OnLoadDictAppend()
{
    QString s = QFileDialog::getOpenFileName(this, "Choose a dictionary file", "",
                                          "Dictionary files (*.pdb,*.cif,*.ent);;All files (*.*)");
    if (!s.isEmpty())
    {
        MIFitDictionary()->LoadDictionary(s.toAscii().constData(), true);
    }
}


void MIMainWindow::OnEditDictResidue()
{
    std::vector<std::string> resList = MIFitDictionary()->GetDictResList();
    if (resList.empty())
    {
        Logger::message("No residues in dictionary to edit");
        return;
    }

    QStringList choices;
    foreach (const std::string &res, resList)
    {
        choices += res.c_str();
    }

    GenericDataDialog dlg(this);
    dlg.setWindowTitle("Edit Residue");
    dlg.addComboField("Select a residue to edit", choices, 0);
    if (dlg.exec() != QDialog::Accepted)
    {
        return;
    }
    std::string type = resList[dlg.value(0).toInt()];

    if (MIFitDictionary()->DictContains(type.c_str()))
    {
        if (MIFitGeomRefiner()->EditEntry(type.c_str()))
        {
            ShowDictEditor(type.c_str());
        }
    }
    else
    {
        std::string rplc_warn;
        rplc_warn += "Cannot find residue \"";
        rplc_warn += type.c_str();
        rplc_warn += "\" in the dictionary.";
        Logger::message(rplc_warn);
    }
}

void MIMainWindow::OnLoadLigMol()
{
    QString filename = QFileDialog::getOpenFileName(this, "Pick a MOL file", "", "MDL file formats (*.mol, *.sdf, *.sd);;MOL files (*.mol);;SD files (*.sdf, *.sd);;All files (*.*)");
    if (filename.isEmpty())
    {
        return;
    }
    std::string code;
    while (true)
    {
        QString str = QInputDialog::getText(this, "Ligand ID", "ID Code (3 letters)", QLineEdit::Normal, "UNK");
        if (str.isEmpty())
        {
            return;
        }
        code = str.toStdString();
        if (code.size() > 3
            || (code.size() > 0 && !isalnum(code[0]))
            || (code.size() > 1 && !isalnum(code[1]))
            || (code.size() > 2 && !isalnum(code[2])))
        {
            QMessageBox::critical(this, "Invalid ID", "ID Code must be alphanumeric and at most 3 letters long");
        }
        else
        {
            break;
        }
    }
    OnLoadLigand("*.mol", filename.toAscii().constData(), "", code.c_str());
    ShowDictEditor(code.c_str());
}

void MIMainWindow::OnLoadLigPdb()
{
    QString filename = QFileDialog::getOpenFileName(this, "Pick a PDB file", "", "PDB files (*.pdb);;All files (*.*)");
    if (filename.isEmpty())
    {
        return;
    }
    std::string code = OnLoadLigand("*.pdb", filename.toAscii().constData(), "", "");
    if (code.size())
    {
        ShowDictEditor(code.c_str());
    }

}

void MIMainWindow::OnLoadLigSmi()
{
    MIData data;
    data["mode"].radio = 0;
    data["mode"].radio_count = 3;
    data["filename"].str = "";
    data["code"].str = "LIG";
    data["smiles"].str = "";
    data["dbquery"].str = "";

    SmilesDialog dlg(this);
    dlg.setWindowTitle("Import Smiles String");
    if (dlg.exec() != QDialog::Accepted)
    {
        return;
    }
    dlg.GetResults(data);

    std::string res;
    switch (data["mode"].radio)
    {
    case 0:  // From File
        res = OnLoadLigand("*.smi", data["filename"].str.c_str(), "", data["code"].str.c_str());
        break;
    case 1:  // From String
        res = OnLoadLigand("*.smi", "", data["smiles"].str.c_str(), data["code"].str.c_str());
        break;
    case 2:  // From Db
    {
        if (Application::instance()->SmilesDbCommand.size())
        {
            QStringList args;
            QString cmd(Application::instance()->SmilesDbCommand.c_str());
            args.append(data["dbquery"].str.c_str());

            QProcess proc(this);
            proc.start(cmd, args);
            if (proc.waitForReadyRead(10000))
            {
                QByteArray qba = proc.readAll();
                std::string str(qba.constData());
                res = OnLoadLigand("*.smi", "", str, data["code"].str.c_str());
            }
            else
            {
                QMessageBox::critical(this, "Command failed", "Error running smiles database command.\nCheck the setting in File/Preferences.../Environment.");
                proc.kill();
                return;
            }
            proc.kill();

            break;
        }
    }
    default:
        return;
        break;
    }

    if (res==std::string(""))
        return;

    ShowDictEditor(data["code"].str.c_str());
}

void MIMainWindow::OnLoadLigCif()
{
    QString filename = QFileDialog::getOpenFileName(this, "Pick a CIF file", "", "CIF files (*.cif);;All files (*.*)");
    if (filename.isEmpty())
    {
        return;
    }
    std::string code = OnLoadLigand("*.cif", filename.toAscii().constData(), "", "");
    if (code.size())
    {
        ShowDictEditor(code.c_str());
    }
}

std::string MIMainWindow::OnLoadLigand(std::string wildcard, std::string filename, std::string smiles, std::string code)
{
    std::string smi_err;                       //Error message returned from smiles library
    std::string err_report;                    //Report of smiles error to user

    // get MI data
    MIMolInfo mi;

    // NOTE: each reading method (fio or LoadSmiles) allocates at least one
    // new RESIDUE (mi.res), which we pass on LigDictEntry, which takes
    // ownership of it, and deletes it when done.

    if (filename.size() == 0 && smiles.size() > 0)
    {
        //Load SMILES from string into dictionary
        if (!LoadSmiles(smiles, code, mi))
        {
            return "";
        }
    }
    else
    {
        MIMolIO fio;
        int sel = fio.getReaderIndex(wildcard.c_str());
        if (!fio.Read(mi, filename.c_str(), sel))
        {
            return "";
        }
    }


    LigDictEntry entry(mi.res);
    entry.bonds = mi.bonds;
    entry.angles = mi.angles;
    entry.torsions = mi.tordict;
    entry.planes = mi.planedict;
    entry.chirals = mi.chiralsdict;


    if (entry.res->atomCount() == 0)
    {
        QMessageBox::critical(this, "Error: No atoms in file", "No atoms were found in the file. Please check the file format.");
        return "";
    }

    //For MOLfile, set flag for whether it has 2D or 3D coordinates
    if (wildcard==std::string("*.mol"))
    {
        entry.res->setName1((entry.res->type() == "3D") ? '.' : 'M');
    }

    if (code.size() > 0)
    {
        entry.res->setType(code);
        entry.res->setName("1");
    }
    if (MIFitDictionary()->DictContains(entry.res->type()))
    {
        std::string rplc_warn = "The residue code you have entered is already in use.  ";
        rplc_warn += "Replace the existing residue \"";
        rplc_warn += entry.res->type().c_str();
        rplc_warn += "\" in the dictionary?";

        int replace = QMessageBox::question(this, "Replace residue?", rplc_warn.c_str(),
                                            QMessageBox::Yes | QMessageBox::No | QMessageBox::Cancel);
        if (replace == QMessageBox::No)
        {
            QString newcode = QInputDialog::getText(this, "Enter Three-Letter Code", "New Code:");
            if (newcode.isEmpty())
            {
                return "";
            }
            else
            {
                entry.res->setType(newcode.toStdString());
            }
            code = newcode.toStdString();
        }
        if (replace == QMessageBox::Cancel)
        {
            return "";
        }
    }

    if (DupeAtomNames(entry.res) != 0)
    {
        if (QMessageBox::question(this, "Duplicate Atom Names", "File contains duplicate atom names...Generate unique names?", QMessageBox::Yes | QMessageBox::No, QMessageBox::No) == QMessageBox::Yes)
        {
            renameResidueAtomsToUnique(entry.res);
        }
        else
        {
            return "";
        }
    }
    if (DupeAtomNames(entry.res) != 0)
    {
        QMessageBox::critical(this, "Error: Duplicate Atom Names", "File contains duplicate atom names...please edit with unique names!");
        return "";
    }

    LigPostProcessor lpp(entry, wildcard.c_str());
    try
    {
        lpp.Process();
    }
    catch (const char *ex)
    {
        std::string message("Unable to load ligand.\n");
        message += ex;
        QMessageBox::critical(this, "Unable to load ligand.", message.c_str());
        return "";
    }
    entry.res->setPrefBonds(entry.bonds);
    entry.res->setPrefAngles(entry.angles);


    unsigned int level;
    if (!MIFitDictionary()->DictHCheck(entry.res, level))
    {
        int ret = QMessageBox::question(this, "H Check", "There is no dictionary, or the number of hydrogens in the dictionary and the residue you are replacing"
                               "are different.\nDo you want to load the appropriate default dictionary?\n"
                               "\nNote: any changes to the current dictionary will be lost.", QMessageBox::Yes | QMessageBox::No | QMessageBox::Cancel);
        if (ret == QMessageBox::Yes)
        {
            std::string path, name, ext;
            SplitPath(Application::instance()->getDictionary(), &path, &name, &ext);
            path += "/";
            if (level == DictionaryHLevel::NoHydrogens)
            {
                path += "dict.noh.pdb";
            }
            if (level == DictionaryHLevel::Polar)
            {
                path += "dict.polarh.pdb";
            }
            if (level == DictionaryHLevel::All)
            {
                path += "dict.allh.pdb";
            }
            MIFitDictionary()->LoadDictionary(path.c_str(), false, true, level);
        }
        if (ret == QMessageBox::Cancel)
        {
            return "";
        }
    }


    //This code substitutes for the "LoadDictionary" function opt->dict.LoadRes(&entry.res, true, true);
    MIFitDictionary()->LoadRes(entry.res, true, true);
    for (std::vector<TORSDICT>::iterator tor = entry.torsions.begin(); tor != entry.torsions.end(); ++tor)
    {
        MIFitDictionary()->AddTorsion(*tor);
    }
    for (std::vector<PLANEDICT>::iterator pln = entry.planes.begin(); pln != entry.planes.end(); ++pln)
    {
        MIFitDictionary()->AddPlane(*pln);
    }
    for (std::vector<CHIRALDICT>::iterator chrl = entry.chirals.begin(); chrl != entry.chirals.end(); ++chrl)
    {
        MIFitDictionary()->AddChiral(*chrl);
    }

    if (MIFitGeomRefiner()->EditEntry(entry.res->type().c_str()))
    {
        std::string ret = entry.res->type().c_str();
        return ret;
    }
    return "";
}

bool MIMainWindow::LoadSmiles(std::string &smiles, std::string &tlcode, MIMolInfo &mi)
{
    std::string smistr(smiles.c_str());
    SMILES smi;
    if (!smi.Read(smistr, mi))
    {
        return false;
    }
    mi.res->setType(tlcode);
    Residue::fixnames(mi.res);
    return true;
}

void MIMainWindow::LeftFooter(const std::string &msg, int timeout)
{
    logWindow->append(msg.c_str());
    statusBar()->showMessage(msg.c_str(), timeout);
}

void MIMainWindow::MiddleFooter(const std::string &msg)
{
    middleFooter->setText(msg.c_str());
}

void MIMainWindow::RightFooter(const std::string &msg)
{
    rightFooter->setText(msg.c_str());
}

void MIMainWindow::OnBackgroundColor()
{
    Application::instance()->backgroundColor();
    MIGLWidget *view = currentMIGLWidget();
    if (view)
        view->ReDrawAll();
}


void MIMainWindow::OnSaveDict()
{
    Application::instance()->saveDict();
}

BatchJobManager *MIMainWindow::GetJobManager()
{
    return JobManager;
}
bool MIMainWindow::isJobLimit()
{
    return JobManager->numberOfRunningJobs() >= Application::instance()->concurrentJobLimit;
}

void MIMainWindow::UpdateToolBar()
{
    if (tool_bar)
        tool_bar->doUpdates();
}

void MIMainWindow::OnUpdateStereoToggle()
{
    bool stereo = MIConfig::Instance()->GetProfileInt("View Parameters", "stereo", 0) != 0;
    _stereoToggleAction->setChecked(stereo);
}

void MIMainWindow::OnUpdateHardwareStereo()
{
    bool hardwareStereo = MIConfig::Instance()->GetProfileInt("View Parameters", "hardwareStereo", 0) != 0;
    _hardwareStereoAction->setEnabled(Application::instance()->isHardwareStereoAvailable());
    _hardwareStereoAction->setChecked(hardwareStereo);
}

void MIMainWindow::OnStereoToggle()
{
    Application::instance()->toggleStereo();
}


void MIMainWindow::OnHardwareStereo()
{
    Application::instance()->toggleHardwareStereo();
}

void MIMainWindow::OnScript()
{
}

void MIMainWindow::OnSideChainTool()
{
    QPoint pos = QCursor::pos();
    side_menu->doExec(pos);
}

void MIMainWindow::OnHideTool()
{
    QPoint pos = QCursor::pos();
    hide_menu->doExec(pos);
}

void MIMainWindow::OnColorTool()
{
    QPoint pos = QCursor::pos();
    color_menu->doExec(pos);
}

void MIMainWindow::OnShowTool()
{
    QPoint pos = QCursor::pos();
    show_menu->doExec(pos);
}

MIMainWindow*MIMainWindow::instance()
{
    if (_instance)
        return _instance;
    new MIMainWindow(); // sets _instance
    return _instance;
}

QDockWidget*MIMainWindow::AddAsDockWidget(QWidget *w, const std::string &name, Qt::DockWidgetArea area)
{
    QDockWidget *dock = new QDockWidget(tr(name.c_str()), this);
    dock->setObjectName(name.c_str());
    w->setParent(dock);
    dock->setWidget(w);
    addDockWidget(area, dock);
    viewMenu->addAction(dock->toggleViewAction());
    return dock;
}

void MIMainWindow::createDockWindows()
{
    logWindow = new QTextEdit(logDock);
    logWindow->setReadOnly(true);
    //  logWindow->document()->setMaximumBlockCount(1000);
    logDock = AddAsDockWidget(logWindow, "Log", Qt::BottomDockWidgetArea);

    // add rama plot
    QWidget *graphwin = RamaPlotMgr::instance()->getGraphWin();
    ramaDock = AddAsDockWidget(graphwin, "Ramachandran plot", Qt::BottomDockWidgetArea);

    //add Models tree
    modelsView = new ModelsView(modelsDock);
    modelsDock = AddAsDockWidget(modelsView, "Models", Qt::LeftDockWidgetArea);

    //add Display tree
    QWidget *displayView = new DisplayView(displayDock);
    displayDock = AddAsDockWidget(displayView, "Display", Qt::LeftDockWidgetArea);

    // add Jobs tree
    JobsView *jobsView = new JobsView(jobsDock);
    jobsView->update(JobManager);
    jobsDock = AddAsDockWidget(jobsView, "Jobs", Qt::LeftDockWidgetArea);
}

static bool is_valid_extension(const char *ext)
{
    return (strcmp(ext, ".pdb") == 0
            || strcmp(ext, ".mlw") == 0
            || strcmp(ext, ".mtz") == 0
            || strcmp(ext, ".phs") == 0
            || strcmp(ext, ".fcf") == 0
            || // we no longer support these extensions b/c their readers are broken
               // || (strcmp(ext, ".sca") == 0)
               // || (strcmp(ext, ".ref") == 0)
            strcmp(ext, ".cif") == 0
            || strcmp(ext, ".map") == 0
            || strcmp(ext, ".py") == 0);
}


void MIMainWindow::OpenFiles(const std::vector<std::string> &files, bool newWin)
{
    if (newWin)
    {
        newWindow();
    }

    for (size_t i = 0; i<files.size(); ++i)
    {
        const std::string &file = files[i];
        const char *arg = file.c_str();

        if (arg[0] == '-')
        {
            // command line switch
            if (strcmp(arg, "-altTextRender") == 0)
            {
                MIConfig::Instance()->WriteProfileInt("View Parameters", "AlternateTextRendering", 1);
                OnExit();
            }
            continue;
        }

        // convert to absolute path
        QFileInfo fi(files[i].c_str());
        QString s = fi.absoluteFilePath();
        std::string ss = s.toStdString();
        arg = ss.c_str();
        const char *ext = file_extension(arg);

        if (!is_valid_extension(ext))
        {
            Log(::format("Unknown file type %s, ingoring.", arg));
            continue;
        }

        if (strcmp(ext, ".py") == 0)
        {
            runPythonScript(file);
            continue;
        }


        if (!currentMIGLWidget())
            newWindow();
        if (!currentMIGLWidget())
        {
            Log(::format("Error creating active MIGLWidget!"));
            continue;
        }

        setCurrentFile(arg);
        currentMIGLWidget()->OpenAnyFile(arg);
    }
}


static bool GetOpenFilenames(std::vector<std::string> &fnames)
{
    initializeFileDialog();
    // TODO use QFileDialog static functions to get native OS file dialog
    QString path = Application::instance()->latestFileBrowseDirectory("");
    fileDialog->setDirectory(path);
    if (fileDialog->exec())
    {
        Application::instance()->latestFileBrowseDirectory(fileDialog->directory().absolutePath());
        foreach (QString file, fileDialog->selectedFiles())
        {
            fnames.push_back(file.toStdString());
        }
        return true;
    }
    return false;
}

static void SplitPath(const std::string &origPath,
                      std::string *dirname,
                      std::string *fname,
                      std::string *ext)
{
    *dirname = std::string("");
    *fname = std::string("");
    *ext = std::string("");

    QFileInfo qfile(origPath.c_str());
    *dirname = qfile.path().toStdString();
    *fname = qfile.fileName().toStdString();
    *ext = qfile.suffix().toStdString();
}

void MIMainWindow::AfterInit()
{
    QStringList args = QApplication::arguments();
    std::vector<std::string> arglist;

#ifdef DEBUG
    if (!(args.count() > 1 && args.at(1) == QString("-recordoff")))
    {
        arglist.push_back("-record");
    }
#endif

    args.pop_front();
    foreach (QString str, args)
    {
        arglist.push_back(str.toStdString());
    }
    OpenFiles(arglist);

}


void MIMainWindow::OnMenuValidate()
{
    menu_bar->validateActions();
    //menu_bar->validateUpdates();
}

#include "../images/cursors/MW_CROSS.xpm"
#include "../images/cursors/MW_TRANS.xpm"
#include "../images/cursors/MW_ROTAT.xpm"
#include "../images/cursors/MW_ZROTA.xpm"
#include "../images/cursors/ZCURSOR.xpm"
#include "../images/cursors/MW_TORSI.xpm"
#include "../images/cursors/MW_CENTE.xpm"
#include "../images/cursors/MW_SCALE.xpm"
#include "../images/cursors/MW_SLAB.xpm"
#include "../images/cursors/slabDrag.xpm"
#include "../images/cursors/mw_wait1.xpm"
#include "../images/cursors/mw_wait2.xpm"
#include "../images/cursors/mw_wait3.xpm"
#include "../images/cursors/mw_wait4.xpm"
#include "../images/cursors/mw_wait5.xpm"
#include "../images/cursors/mw_wait6.xpm"

void MIMainWindow::SetCursor(int id, QWidget *w)
{
    if (!cursors[0])
    {
        cursors[imhCross] = new QCursor(QPixmap(MW_CROSS));
        cursors[imhTranslate] = new QCursor(QPixmap(MW_TRANS));
        cursors[imhRotate] = new QCursor(QPixmap(MW_ROTAT));
        cursors[imhZRotate] = new QCursor(QPixmap(MW_ZROTA));
        cursors[imhZCursor] = new QCursor(QPixmap(ZCURSOR));
        cursors[imhTorsion] = new QCursor(QPixmap(MW_TORSI));
        cursors[imhCenter] = new QCursor(QPixmap(MW_CENTE));
        cursors[imhScale] = new QCursor(QPixmap(MW_SCALE));
        cursors[imhSlab] = new QCursor(QPixmap(MW_SLAB));
        cursors[imhWait1] = new QCursor(QPixmap(mw_wait1));
        cursors[imhWait2] = new QCursor(QPixmap(mw_wait2));
        cursors[imhWait3] = new QCursor(QPixmap(mw_wait3));
        cursors[imhWait4] = new QCursor(QPixmap(mw_wait4));
        cursors[imhWait5] = new QCursor(QPixmap(mw_wait5));
        cursors[imhWait6] = new QCursor(QPixmap(mw_wait6));
        cursors[imhSlabDrag] = new QCursor(QPixmap(slabDrag));
    }

    if (!w)
        w = this;

    if (id >=0 && id <= (int)imhCount)
        w->setCursor(*cursors[id]);
    else
        w->unsetCursor();
}

void MIMainWindow::updateNavigator()
{
    navigator->updateGL();
}

void MIMainWindow::Debug(const std::string &msg)
{
#ifdef DEBUG
    Log(msg);
#endif
}

void MIMainWindow::closeEvent(QCloseEvent *event)
{
    saveLayout();
    mdiArea->closeAllSubWindows();
    if (currentMIGLWidget())
    {
        event->ignore();
    }
    else
    {
        event->accept();
    }
}

void MIMainWindow::newWindow()
{
    createMIGLWidget()->show();
}



void MIMainWindow::updateMenus()
{
    updateFileMenu();

    bool hasMIGLWidget = (currentMIGLWidget() != 0);
    closeAllAct->setEnabled(hasMIGLWidget);
    tileAct->setEnabled(hasMIGLWidget);
    cascadeAct->setEnabled(hasMIGLWidget);
    nextAct->setEnabled(hasMIGLWidget);
    previousAct->setEnabled(hasMIGLWidget);

    view_menu->setEnabled(hasMIGLWidget);
    show_menu->setEnabled(hasMIGLWidget);
    render_menu->setEnabled(hasMIGLWidget);
    model_menu->setEnabled(hasMIGLWidget);
    fit_menu->setEnabled(hasMIGLWidget);
    refi_menu->setEnabled(hasMIGLWidget);
    analyze_menu->setEnabled(hasMIGLWidget);
}

void MIMainWindow::updateWindowMenu()
{
    windowMenu->clear();
    windowMenu->addAction(closeAct);
    windowMenu->addAction(closeAllAct);
    windowMenu->addSeparator();
    windowMenu->addAction(tileAct);
    windowMenu->addAction(cascadeAct);
    windowMenu->addSeparator();
    windowMenu->addAction(nextAct);
    windowMenu->addAction(previousAct);

    QList<QMdiSubWindow*> windows = mdiArea->subWindowList();

    for (int i = 0; i < windows.size(); ++i)
    {
        MIGLWidget *child = dynamic_cast<MIGLWidget*>(windows.at(i)->widget());
        if (!child)
            continue;

        QString text;
        if (i == 0)
            windowMenu->addSeparator();
        if (i < 9)
        {
            text = tr("&%1 %2").arg(i + 1)
                   .arg(child->GetTitle().c_str());
        }
        else
        {
            text = tr("%1 %2").arg(i + 1)
                   .arg(child->GetTitle().c_str());
        }
        QAction *action  = windowMenu->addAction(text);
        action->setCheckable(true);
        action->setChecked(child == currentMIGLWidget());
        connect(action, SIGNAL(triggered()), windowMapper, SLOT(map()));
        windowMapper->setMapping(action, windows.at(i));
    }
}

QWidget*MIMainWindow::createMIGLWidget()
{
    QWidget *w = new MIGLWidget;
    w->setObjectName("MIGLWidget");
    mdiArea->addSubWindow(w);
    w->showMaximized();
    foreach (QKeySequence k, closeAction->shortcuts())
        grabShortcut(k, Qt::ApplicationShortcut);
    return w;
}

void MIMainWindow::createActions()
{

    closeAct = new QAction(tr("Cl&ose\tCtrl+W"), this);
    closeAct->setStatusTip(tr("Close the active window"));
    connect(closeAct, SIGNAL(triggered()),
            mdiArea, SLOT(closeActiveSubWindow()));

    closeAllAct = new QAction(tr("Close &All"), this);
    closeAllAct->setStatusTip(tr("Close all the windows"));
    connect(closeAllAct, SIGNAL(triggered()),
            mdiArea, SLOT(closeAllSubWindows()));
    ///////////////////////////////////////////////////////////

    for (int i = 0; i < MaxRecentFiles; ++i)
    {
        recentFileActs[i] = new QAction(this);
        recentFileActs[i]->setVisible(false);
        connect(recentFileActs[i], SIGNAL(triggered()),
                this, SLOT(openRecentFile()));
    }

    tileAct = new QAction(tr("&Tile"), this);
    tileAct->setStatusTip(tr("Tile the windows"));
    connect(tileAct, SIGNAL(triggered()), mdiArea, SLOT(tileSubWindows()));

    cascadeAct = new QAction(tr("&Cascade"), this);
    cascadeAct->setStatusTip(tr("Cascade the windows"));
    connect(cascadeAct, SIGNAL(triggered()), mdiArea, SLOT(cascadeSubWindows()));

    nextAct = new QAction(tr("Ne&xt"), this);
    nextAct->setStatusTip(tr("Move the focus to the next window"));
    connect(nextAct, SIGNAL(triggered()),
            mdiArea, SLOT(activateNextSubWindow()));

    previousAct = new QAction(tr("Pre&vious"), this);
    previousAct->setStatusTip(tr("Move the focus to the previous "
                                 "window"));
    connect(previousAct, SIGNAL(triggered()),
            mdiArea, SLOT(activatePreviousSubWindow()));

}

void MIMainWindow::updateFileMenu()
{
    bool hasDocument = currentMIGLWidget() != 0;
    fileSaveAction->setEnabled(hasDocument);
    fileSaveAsAction->setEnabled(hasDocument);
    bool canExport = hasDocument
                     && currentMIGLWidget()->GetDisplaylist()->CurrentItem() != NULL
                     && !MIBusyManager::instance()->Busy();
    exportModelAction->setEnabled(canExport);
    filePrintAction->setEnabled(hasDocument);
    copyCanvasAction->setEnabled(hasDocument);
    exportImageAction->setEnabled(hasDocument);

    closeAction->setEnabled(hasDocument);
}

void MIMainWindow::OnFileSave()
{
    if (currentMIGLWidget())
        currentMIGLWidget()->OnFileSave();
}

void MIMainWindow::OnFileSaveAs()
{
    if (currentMIGLWidget())
        currentMIGLWidget()->OnFileSaveAs();
}

void MIMainWindow::OnExportModel()
{
    if (currentMIGLWidget())
        currentMIGLWidget()->OnExportModel();
}

void MIMainWindow::OnPrint()
{
    if (currentMIGLWidget())
        currentMIGLWidget()->OnPrint();
}

void MIMainWindow::OnEditCopy()
{
    if (currentMIGLWidget())
        currentMIGLWidget()->OnEditCopy();
}

void MIMainWindow::OnExportImage()
{
    if (currentMIGLWidget())
        currentMIGLWidget()->OnExportImage();
}

void MIMainWindow::fill_file_menu(QMenu *file_menu)
{
    QAction *action;

    fileNewAction = file_menu->addAction(tr("New"), this, SLOT(OnNew()));
    fileNewAction->setToolTip(tr("Create a new document"));
    fileNewAction->setShortcut(tr("Ctrl+N"));

    fileOpenAction = file_menu->addAction(tr("Open models, data, maps, etc..."), this, SLOT(OnFileOpen()));
    fileOpenAction->setToolTip(tr("Open file(s) into the current document"));
    fileOpenAction->setShortcut(tr("Ctrl+O"));

    action = file_menu->addAction(tr("Open into new document..."), this, SLOT(OnFileOpenNew()));
    action->setToolTip(tr("Open file(s) into a new document"));
    action->setShortcut(tr("Ctrl+Shift+O"));

    closeAction = file_menu->addAction(tr("Close\tCtrl+W"), this, SLOT(OnClose()));

    fileSaveAction = file_menu->addAction(tr("Save Session"), this, SLOT(OnFileSave()));
    fileSaveAction->setToolTip(tr("Save the entire document to a file"));
    fileSaveAction->setShortcut(tr("Ctrl+S"));

    fileSaveAsAction = file_menu->addAction(tr("Save Session As..."), this, SLOT(OnFileSaveAs()));
    fileSaveAsAction->setToolTip(tr("Save the document under a new name"));
    fileSaveAsAction->setShortcut(tr("Ctrl+Shift+S"));

    exportModelAction = file_menu->addAction(tr("Export PDB..."), this, SLOT(OnExportModel()));
    exportModelAction->setToolTip(tr("Saves the current model as a PDB format file"));
    exportModelAction->setShortcut(tr("Ctrl+Alt+S"));

    filePrintAction = file_menu->addAction(tr("Print..."), this, SLOT(OnPrint()));
    filePrintAction->setToolTip(tr("Print the canvas"));

    file_menu->addSeparator();

    copyCanvasAction = file_menu->addAction(tr("Copy Canvas"), this, SLOT(OnEditCopy()));
    copyCanvasAction->setToolTip(tr("Copy the canvas to the clipboard to paste into other programs"));
    copyCanvasAction->setShortcut(tr("Ctrl+C"));

    exportImageAction = file_menu->addAction(tr("Export Image As..."), this, SLOT(OnExportImage()));
    exportImageAction->setToolTip(tr("Export view to a graphics file"));

    file_menu->addSeparator();

    action = file_menu->addAction(tr("Manage Crystals..."), this, SLOT(OnManageCrystals()));
    action->setToolTip(tr("Edit the crystal database"));

    file_menu->addSeparator();

    action = file_menu->addAction(tr("Preferences..."), this, SLOT(OnPreferences()));
    action->setToolTip(tr("Edit preferences"));

    action = file_menu->addAction(tr("Define Atom Colors..."), this, SLOT(OnDefineAtomColors()));
    action->setToolTip(tr("Define Atom Colors"));

    action = file_menu->addAction(tr("B-Value Colors Ranges..."), this, SLOT(OnDefineBValueColors()));
    action->setToolTip(tr("Define B-Value Coloring Ranges"));

    MIMainWindow::instance()->addRecentFileActions(file_menu);

    action = file_menu->addAction(tr("E&xit"), this, SLOT(OnExit()));
    action->setShortcut(tr("Alt+X"));

    connect(file_menu, SIGNAL(aboutToShow()), SLOT(updateFileMenu()));
}



void MIMainWindow::fill_view_menu(MIMenu *view_menu)
{
    view_menu->Append(ID_VIEW_UNDO, "&Undo Viewpoint", "Use this command to center objects", false);
    view_menu->Append(ID_VIEW_SAVE, "&Save Viewpoint...", "Use this command to save the viewpoint to a file", false);
    view_menu->Append(ID_VIEW_LOAD, "&Load Viewpoint...", "Use this command to load the viewpoint from a file", false);
    view_menu->AppendSeparator();
    view_menu->Append(ID_VIEW_TOPVIEW, "&Top view", "View from top and drag clipping planes", true);
    view_menu->Append(ID_VIEW_FULLSCREEN, "&Fullscreen\tESC", "Canvas fills the entire screen", true);

    _hardwareStereoAction = view_menu->addAction(tr("Use &Hardware Stereo"), this, SLOT(OnHardwareStereo()));
    _hardwareStereoAction->setCheckable(true);
    connect(view_menu, SIGNAL(aboutToShow()), this, SLOT(OnUpdateHardwareStereo()));

    _stereoToggleAction = view_menu->addAction(tr("S&tereo\t|"), this, SLOT(OnStereoToggle()));
    _stereoToggleAction->setCheckable(true);
    connect(view_menu, SIGNAL(aboutToShow()), this, SLOT(OnUpdateStereoToggle()));

    view_menu->AppendSeparator();
    view_menu->Append(ID_VIEW_SLABIN, "Slab In\tShift+I", "Decrease the distance beteween the front and back clipping planes", false);
    view_menu->Append(ID_VIEW_SLABOUT, "Slab Out\tShift+O", "Increase the distance beteween the front and back clipping planes", false);
    view_menu->Append(ID_VIEW_CLIPPLANES, "Set s&lab...", "Set the values of the front and back clipping planes", false);
    view_menu->Append(ID_GOTO_ZOOMIIN, "Zoom In\tI", "Zoom viewpoint in by 20%", false);
    view_menu->Append(ID_GOTO_ZOOMOUT, "Zoom Out\tO", "Zoom viewpoint in by 20%", false);
    view_menu->Append(ID_VIEW_ROTATEY90, "Rotate View &+90", "Use this command (and -90) to center objects", false);
    view_menu->Append(ID_VIEW_ROTATEYMINUS, "Rotate View &-90", "Use this command to center objects", false);
    view_menu->Append(ID_VIEW_ORTHONORMAL, "Orthonor&mal", "Set perspective to 0 (infinite viewpoint)", true);
    view_menu->Append(ID_INCREASE_PERSP, "I&ncrease Perspective", "Increase perspective +20%", false);
    view_menu->Append(ID_DECREASE_PERSP, "&Decrease Perspective", "Decrease perspective -20%", false);
    view_menu->AppendSeparator();
    view_menu->Append(ID_GOTO_FITTOSCREEN, "&Center model on screen", "Center molecule on screen", false);
    view_menu->Append(ID_GOTO_FITALLTOSCREEN, "&Center All Models On Screen", "Center all the molecules on the screen", false);
    view_menu->Append(ID_MAP_CENTERDENSITY, "Center Visible Density", "Center the density visible on the screen", false);
    view_menu->Append(ID_GOTO_GOTOXYZ, "&Go to x,y,z...", "Center view at a coordinate in Angstroms", false);
    view_menu->AppendSeparator();

    MIMenu *anime_menu = new MIMenu(view_menu->GetReceiver());
    anime_menu->Append(ID_ANIMATE_ROCK, "&Rock", "Rock the model back and forth", true);
    anime_menu->Append(ID_ANIMATE_ROLL, "Ro&ll", "Roll the model around the vertical", true);
    anime_menu->Append(ID_ANIMATE_BLINK, "&Blink", "Blink between two or more models", true);
    anime_menu->Append(ID_ANIMATE_ROCKANDROLLPARAMETERS, "Rock/Roll Rates...", "Set the rates for rock and roll", false);
    view_menu->Append(ID_ANIMEMENU, "&Animate", anime_menu);
}

void MIMainWindow::fill_surf_menu(QMenu *surf_menu)
{
    new CurrentMIGLWidgetAction("Surface Residue", "Calculate van der Waal dot surface around picked residue", surf_menu, SLOT(OnObjectSurfaceresidue()), SLOT(OnUpdateObjectSurfaceresidue(QAction*)));

    new CurrentMIGLWidgetAction("Surface Residues", "Surface all the residues on the stack", surf_menu, SLOT(OnObjectSurfaceresidues()), SLOT(OnUpdateObjectSurfaceresidues(QAction*)));

    new CurrentMIGLWidgetAction("Surface Atom", "Surface the last atom picked", surf_menu, SLOT(OnObjectSurfaceatom()), SLOT(OnUpdateObjectSurfaceAtom(QAction*)));

    new CurrentMIGLWidgetAction("Surface Atoms", "Surface all the atoms on the stack", surf_menu, SLOT(OnObjectSurfaceAtoms()), SLOT(OnUpdateObjectSurfaceatoms(QAction*)));

    new CurrentMIGLWidgetAction("van der Waal Surface", "Calculate van der Waal dot surface", surf_menu, SLOT(OnVdwDotSurface()));

    new CurrentMIGLWidgetAction("Solvent Exposed Surface", "Surface through the center of the solvent atoms touching the molecule", surf_menu, SLOT(OnSurfaceSolvent()), SLOT(OnUpdateSurfaceSolvent(QAction*)));

    new CurrentMIGLWidgetAction("Sphere around atom", "Calculate sphere around last picked atom", surf_menu, SLOT(OnObjectSurfaceSpherearoundatom()), SLOT(OnUpdateObjectSurfaceSpherearoundatom(QAction*)));

    surf_menu->addSeparator();
    new CurrentMIGLWidgetAction("Clear Surface", "Clears the surface dots", surf_menu, SLOT(OnObjectSurfaceClearsurface()));

}


void MIMainWindow::createMenus()
{
    menu_bar = new MIMenuBar(menuBar());

    //// Make a menubar
    MIMenu *file_menu = new MIMenu(*this);
    fill_file_menu(file_menu);
    updateRecentFileActions();



    MIMenu *stack_menu = new MIMenu(*this);
    stack_menu->Append(ID_VIEW_ATOMSTACK, "Show/Hide Atom Stack", "Toggle the atom stack visibility", true);
    stack_menu->AppendSeparator();
    stack_menu->Append(ID_OBJECT_STACK_DELETETOPITEM, "Clear &top item", "Clear the top item on the stack", false);
    stack_menu->Append(ID_OBJECT_CLEARSTACK, "&Clear stack", "Clear all items on the stack", false);
    stack_menu->AppendSeparator();
    stack_menu->Append(ID_OBJECT_STACK_EXPANDTOPALLATOMSINRESIDUE, "&Expand Top Residue", "Expand the top residue to put all of its atoms on the stack", false);
    stack_menu->Append(ID_OBJECT_STACK_EXPANDTOP2RESIDUES, "Expand &range to All Residue", "Expand stack to include all residues between top two", false);
    stack_menu->Append(ID_OBJECT_STACK_EXPANDTOP2ALLATOMSINRANGE, "Expand Range to All &Atoms", "Expand stack to include all atoms in all residues between top two", false);

    color_menu = new MIMenu(*this);
    color_menu->Append(ID_SHOW_COLORALLATOMS, "Color all atoms", "Colors the whole model with current color", false);
    color_menu->Append(ID_OBJECT_ATOM_COLOR, "Color last picked atom", "Colors the last picked atom with current color", false);
    color_menu->Append(ID_OBJECT_ATOMS_COLOR, "Color all picked atoms", "Colors all the atoms on the stack with current color", false);
    color_menu->AppendSeparator();
    color_menu->Append(ID_OBJECT_RESIDUE_COLOR, "Color last picked residue", "Colors the last picked residue with current color", false);
    color_menu->Append(ID_OBJECT_RESIDUES_COLOR, "Color all picked Residues", "Colors all the residues on the stack with current color", false);
    color_menu->Append(ID_OBJECT_RESIDUERANGE_COLOR, "Color residue range", "Colors the residue at the top of stack with current color", false);
    color_menu->AppendSeparator();
    color_menu->Append(ID_SHOW_UNDOCOLORRADIUS, "&Undo color", "Undo the last color or radius command", false);

    render_menu = new QMenu(this);
    connect(render_menu, SIGNAL(aboutToShow()), SLOT(updateRenderMenu()));

    QActionGroup *renderStyleActionGroup = new QActionGroup(render_menu);
    QAction *action;
    action = new CurrentMIGLWidgetAction("&Sticks", "Draw bonds as sticks", render_menu, SLOT(OnRenderSticks()), SLOT(OnUpdateRenderSticks(QAction*)));
    action->setActionGroup(renderStyleActionGroup);
    action->setCheckable(true);

    action = new CurrentMIGLWidgetAction("&Knob and Stick", "Draw atoms as knobs and bonds as stick", render_menu, SLOT(OnRenderingBallandstick()), SLOT(OnUpdateRenderingBallandstick(QAction*)));
    action->setActionGroup(renderStyleActionGroup);
    action->setCheckable(true);

    action = new CurrentMIGLWidgetAction("&Ball and Cylinder", "Draw bonds as ball and cylinders for bonds", render_menu, SLOT(OnRenderBallandcylinder()), SLOT(OnUpdateRenderBallandcylinder(QAction*)));
    action->setActionGroup(renderStyleActionGroup);
    action->setCheckable(true);

    action = new CurrentMIGLWidgetAction("&CPK", "Draw atoms as CPK/space filling", render_menu, SLOT(OnRenderSpacefilling()), SLOT(OnUpdateRenderSpacefilling(QAction*)));
    action->setActionGroup(renderStyleActionGroup);
    action->setCheckable(true);

    render_menu->addSeparator();
    color_menu->setTitle("Color");
    render_menu->addMenu(color_menu);
    render_menu->addAction("Set Background C&olor...", this, SLOT(OnBackgroundColor()));
    render_menu->addSeparator();
    new CurrentMIGLWidgetAction("Set Ball/Cylinder S&ize...", "Set the relative diameters of ball and cylinder", render_menu, SLOT(OnRenderBallsize()), SLOT(OnUpdateRenderBallsize(QAction*)));

    new CurrentMIGLWidgetAction("Set View Center Target S&ize...", "Set size of the marker in the view center", render_menu, SLOT(OnRenderTargetSize()));

    render_menu->addSeparator();
    action = new CurrentMIGLWidgetAction("Depthcue Colo&rs", "Depthcue using darker colors", render_menu, SLOT(OnRenderingDepthcuedcolors()), SLOT(OnUpdateRenderingDepthcuedcolors(QAction*)));
    action->setCheckable(true);

    action = new CurrentMIGLWidgetAction("Depthcue &Lines", "Depthcue with line thickness", render_menu, SLOT(OnRenderingDepthcuedlinewidth()), SLOT(OnUpdateRenderingDepthcuedlinewidth(QAction*)));
    action->setCheckable(true);

    action = new CurrentMIGLWidgetAction("S&mooth Lines", "Smooths the lines", render_menu, SLOT(OnRenderLinesmooth()), SLOT(OnUpdateRenderLinesmooth(QAction*)));
    action->setCheckable(true);

    renderLineThicknessMenu_ = render_menu->addMenu("Line &Thickness");
    QSignalMapper *lineMenuSignalMapper = new QSignalMapper(renderLineThicknessMenu_);
    connect(lineMenuSignalMapper, SIGNAL(mapped(int)), SLOT(setRenderLineThickness(int)));
    QAction *lineThickness1 = renderLineThicknessMenu_->addAction("&1 pixel", lineMenuSignalMapper, SLOT(map()));
    lineMenuSignalMapper->setMapping(lineThickness1, 1);
    action = renderLineThicknessMenu_->addAction("&2 pixel", lineMenuSignalMapper, SLOT(map()));
    lineMenuSignalMapper->setMapping(action, 2);
    action = renderLineThicknessMenu_->addAction("&3 pixel", lineMenuSignalMapper, SLOT(map()));
    lineMenuSignalMapper->setMapping(action, 3);
    action = renderLineThicknessMenu_->addAction("&4 pixel", lineMenuSignalMapper, SLOT(map()));
    lineMenuSignalMapper->setMapping(action, 4);
    QActionGroup *lineThicknessGroup = new QActionGroup(renderLineThicknessMenu_);
    foreach (QAction *a, renderLineThicknessMenu_->actions())
    {
        a->setCheckable(true);
        lineThicknessGroup->addAction(a);
    }
    lineThickness1->setChecked(true);

    action = new CurrentMIGLWidgetAction("Dim Non-active Models", "Dim non-active models", render_menu, SLOT(OnDimNonactiveModels()), SLOT(OnUpdateDimNonactiveModels(QAction*)));
    action->setCheckable(true);

    new CurrentMIGLWidgetAction("Set Amount to Dim Non-active Models", "Dim non-active models", render_menu, SLOT(OnAmountToDimNonactiveModels()));


    MIMenu *showres_menu = new MIMenu(*this);
    showres_menu->Append(ID_OBJECTS_ALLATOMS, "Show &all residues", "Show the whole model", false);
    showres_menu->Append(ID_OBJECT_SHOWRESIDUE, "Show &last picked residue", "Show the last picked residue", false);
    showres_menu->Append(ID_OBJECT_SHOWRESIDUES, "Show &all picked residues", "Show all the residues on the stack", false);
    showres_menu->Append(ID_OBJECT_SHOWRESIDUERANGE, "Show &residue range", "Show the residue at the top of stack", false);
    showres_menu->Append(ID_SHOW_SHOWWITHINSPHERE, "Show residues within &sphere...", "Show residue if it is witin a radius", false);
    showres_menu->AppendSeparator();
    showres_menu->Append(ID_HIDEMODEL, "&Hide all residues", "Hide the whole model", false);
    showres_menu->Append(ID_OBJECT_RESIDUE_TURNOFF, "Hide last picked residue", "Hide the last picked residue", false);
    showres_menu->Append(ID_OBJECT_RESIDUES_TURNOFF, "Hide all picked residues", "Hide all the residues on the stack", false);
    showres_menu->Append(ID_OBJECT_RESIDUERANGE_TURNOFF, "Hide residue range", "Hide the residue at the top of stack", false);

    side_menu = new MIMenu(*this);
    side_menu->Append(ID_SHOW_SIDECHAINATOMS, "Show sidechain &atoms", "Show sidechain atoms", false);
    side_menu->Append(ID_SHOW_HIDESIDECHAINATOMS, "&Hide sidechain atoms", "Hide sidechain atoms", false);
    side_menu->Append(ID_OBJECT_SHOWSIDECHAIN, "&Show sidechain of last picked", "Show the sidechain of the last picked residue", false);
    side_menu->Append(ID_OBJECT_SHOWSIDECHAINS, "S&how sidechains of all picked", "Show the sidechains of all the residues on the stack", false);
    side_menu->Append(ID_OBJECT_SHOWSIDECHAINRANGE, "Sh&ow sidechains of range", "Show the sidechains of the residue range at the top of stack", false);

    MIMenu *backbone_menu = new MIMenu(*this);
    backbone_menu->Append(ID_SHOW_BACKBONEATOMS, "Show backbone as &atoms", "Show backbone as atoms", false);
    backbone_menu->Append(ID_SHOW_BACKBONECA, "Show backbone as &CA trace", "Show backbone as CA trace", false);
    backbone_menu->Append(ID_SHOW_HIDEBACKBONE, "&Hide backbone", "Hide backbone", false);

    MIMenu *symmatoms_menu = new MIMenu(*this);
    symmatoms_menu->Append(ID_SHOW_SYMMATOMSASATOMS, "Show symmetry atoms as &atoms");
    symmatoms_menu->Append(ID_SHOW_SYMMATOMSASCA, "Show symmetry atoms as &CA trace");
    symmatoms_menu->Append(ID_SHOW_HIDESYMMATOMS, "&Hide symmetry atoms");
    symmatoms_menu->Append(ID_SHOW_SAVESYMMATOMS, "&Save symmetry atoms");

    canvas_menu = new MIMenu(*this);
    canvas_menu->Append(ID_OBJECT_STACKMENU, "S&tack", stack_menu);
    canvas_menu->Append(ID_VIEW_UNITCELL, "&Unit cell", "Toggle the unit cell visibility", true);
    canvas_menu->Append(ID_VIEW_CONTACTS, "&Contacts", "Toggle contacts on/off", true);
    canvas_menu->Append(ID_VIEW_GNOMON, "&Gnomon", "Toggle the axes visibility", true);
    canvas_menu->Append(ID_VIEW_CLEARMESSAGE, "Clea&r Message", "Use this command to clear the message at the screen bottom", false);

    MIMenu *labels_menu = new MIMenu(*this);
    labels_menu->Append(ID_VIEW_LABELS, "&Show/Hide Labels", "Toggle labels on and off", true);
    new CurrentMIGLWidgetAction("Cl&ear Labels",
                                "Clear the labels permanently - use view/labels to hide them temporarily",
                                labels_menu, SLOT(OnEditClearlabels()));

    new CurrentMIGLWidgetAction("La&bel Options...",
                                "Edit labels and set label options",
                                labels_menu, SLOT(OnEditLabels()));

    labels_menu->Append(ID_SHOW_LABELEVERYNTH, "Label E&very nth...", "Label every nth residue", false);

    MIMenu *secstruct_menu = new MIMenu(*this);
    secstruct_menu->Append(ID_OBJECT_BACKBONERIBBON, "&Make Ribbon", "Make a backbone ribbon", false);
    secstruct_menu->Append(ID_OBJECT_CLEARRIBBON, "C&lear Ribbon", "Clear the backbone ribbon", false);
    secstruct_menu->Append(ID_OBJECTS_SECONDARYSTRUCTURE, "Ri&bbon Colors...", "Set ribbon colours", false);
    secstruct_menu->AppendSeparator();
    secstruct_menu->Append(ID_RIBBONSECONDARYSTRUCTURE, "Show Tube Secondary Structure", "Display a solid tube for the secondary structure", false);
    secstruct_menu->Append(ID_SCHEMATICSECONDARYSTRUCTURE, "Show Schematic Secondary Structure", "Display a Schematic representation for the secondary structure", false);
    secstruct_menu->Append(ID_DELETESECONDARYSTRUCTURE, "Hide Secondary Structure", "Removes the display of the secondary structure", false);
    MIMenu *secstruct_options_menu = new MIMenu(*this);
    secstruct_menu->Append(0, "Options", secstruct_options_menu);
    secstruct_options_menu->Append(ID_SECONDARYSTRUCTUREOPTIONS_TUBE, "Tube", "Set secondary structure options", false);
    secstruct_options_menu->Append(ID_SECONDARYSTRUCTUREOPTIONS_SHEET, "Beta sheet", "Set secondary structure options", false);
    secstruct_options_menu->Append(ID_SECONDARYSTRUCTUREOPTIONS_TURN, "Turn", "Set secondary structure options", false);
    secstruct_options_menu->Append(ID_SECONDARYSTRUCTUREOPTIONS_RANDOM, "Random coil", "Set secondary structure options", false);
    secstruct_options_menu->Append(ID_SECONDARYSTRUCTUREOPTIONS_HELIX, "Helix", "Set secondary structure options", false);


    show_menu = new MIMenu(*this);
    connect(show_menu, SIGNAL(aboutToShow()), this, SLOT(updateShowMenu()));

    MIMenu *showatoms_menu = new MIMenu(*this);
    show_menu->Append(ID_SHOWATOMSMENU, "Atoms", showatoms_menu);
    showatoms_menu->Append(ID_OBJECTS_ALLATOMS, "Show &all atoms", "Show all the atoms in the model", false);
    showatoms_menu->Append(ID_SHOW_PICKEDATOM_TURNON, "Show last &picked atom", "Show the last picked atom", false);
    showatoms_menu->Append(ID_SHOW_ALLPICKEDATOMS_TURNON, "Show all p&icked atoms", "Show all the atoms on the stack", false);
    showatoms_menu->AppendSeparator();
    showatoms_menu->Append(ID_HIDEMODEL, "&Hide all atoms", "Hide the whole model", false);
    showatoms_menu->Append(ID_SHOW_PICKEDATOM_TURNOFF, "Hide last picked a&tom", "Hide the last picked atom", false);
    showatoms_menu->Append(ID_SHOW_ALLPICKEDATOMS_TURNOFF, "Hide all picked at&oms", "Hide all the atoms on the stack", false);
    showatoms_menu->Append(ID_SHOW_HIDEHYDROGENS, "Toggle h&ydrogens", "Toggle hydrogens on/off", true);


    show_menu->Append(ID_SHOWRESMENU, "Residues", showres_menu);
    show_menu->Append(ID_BACKBONEMENU, "Backbone", backbone_menu);
    show_menu->Append(ID_SIDEMENU, "Sidechains", side_menu);
    show_menu->Append(ID_SECONDARYSTRUCTUREMENU, "Secondary Structure", secstruct_menu);
    show_menu->Append(ID_LABELSMENU, "Labels", labels_menu);

    QMenu *annotation_menu = show_menu->addMenu("Annotation");
    new CurrentMIGLWidgetAction("A&dd annotation to model",
                                "Adds an annotation to the current model at the screen center",
                                annotation_menu,
                                SLOT(OnAnnotation()),
                                SLOT(OnUpdateAnnotation(QAction*)));

    new CurrentMIGLWidgetAction("Move annotation to &picked atom",
                                "Move the last picked Annotation to the last picked atom ccords",
                                annotation_menu,
                                SLOT(OnMoveAnnotationAtom()),
                                SLOT(OnUpdateMoveAnnotationAtom(QAction*)));

    new CurrentMIGLWidgetAction("Move annotation to c&enter",
                                "Move the last picked Annotation to the screen center",
                                annotation_menu,
                                SLOT(OnMoveAnnotationCenter()),
                                SLOT(OnUpdateMoveAnnotationCenter(QAction*)));

    show_menu->Append(ID_SYMMATOMSMENU, "Symmetry Atoms", symmatoms_menu);

    QMenu *surf_menu = show_menu->addMenu("&Dot Surface");
    fill_surf_menu(surf_menu);

    show_menu->Append(ID_CANVASMENU, "Canvas", canvas_menu);

    view_menu = new MIMenu(*this);
    fill_view_menu(view_menu);

    analyze_menu = new MIMenu(*this);
    analyze_menu->Append(ID_GEOMETRY_DISTANCE, "Measure &Distance", "Measure the distance between last two picked atoms", false);
    analyze_menu->Append(ID_GEOMETRY_ANGLE, "Measure &Angle", "Measure the angle of last three picked atoms", false);
    analyze_menu->Append(ID_GEOMETRY_TORSION, "Measure D&ihedral", "Measure the dihedral/torsion of last four picked atoms", false);
    analyze_menu->AppendSeparator();
    analyze_menu->AppendCheckItem(ID_RAMAPLOT_ALLOWED, "Detailed Ramachandran plot", "Show allowed regions in Ramachandran plot, too");
    analyze_menu->AppendSeparator();
    analyze_menu->Append(ID_GEOM_FINDGEOMERRORS, "A&nalyze Geometry Errors", "Find all the geometry errors above a threshold in the model", false);
    analyze_menu->Append(ID_GEOM_CLEARGEOMERRORS, "Clear &Geometry Errors", "Clear all the Annotation of geometry errors in the model", false);
    analyze_menu->AppendSeparator();
    analyze_menu->Append(ID_GEOM_ADDSINGLEHBOND, "Add &H-bond", "Build an H-bond between the last two picked atoms", false);
    analyze_menu->Append(ID_GEOM_HBONDS, "B&uild H-Bonds", "Build all the H-bonds in the model", false);
    analyze_menu->Append(ID_GEOM_CLEARHBONDS, "&Clear H-bonds", "Clear all the H-bonds in the model", false);
    analyze_menu->AppendSeparator();
    analyze_menu->Append(ID_GEOM_NEIGHBOURS, "Bui&ld Contacts", "Build all the Contacts in the model", false);
    analyze_menu->Append(ID_GEOM_CLEARNEIGHBOURS, "Cl&ear Contacts", "Clear all the Contacts in the model", false);

    // Refine menu actions which have shortcuts must be immediately updated
    // when refinement state changes. They can't wait for the update which occurs
    // when the menu displays.
    refi_menu = new QMenu(this);
    refineResidueAction = new CurrentMIGLWidgetAction("R&efine Residue\tCtrl+R", "Real-space refine the last picked residue", refi_menu, SLOT(OnRefiResidue()), SLOT(OnUpdateRefiResidue(QAction*)));

    connect(MIFitGeomRefiner(), SIGNAL(isRefiningChanged(bool)),
            this, SLOT(updateIsRefining(bool)));
    new CurrentMIGLWidgetAction("Re&fine Local Region", "Real-space refine the last picked residue and its 2 neighbours", refi_menu, SLOT(OnRefiRegion()), SLOT(OnUpdateRefiRegion(QAction*)));

    new CurrentMIGLWidgetAction("Ref&ine Range", "Real-space refine the last 2 picks and the intervening residues", refi_menu, SLOT(OnRefiRange()), SLOT(OnUpdateRefiRange(QAction*)));

    new CurrentMIGLWidgetAction("Refi&ne Molecule", "Real-space refine the entire molecule in all 6 dimensions", refi_menu, SLOT(OnRefiMolecule()), SLOT(OnUpdateRefiMolecule(QAction*)));

    refi_menu->addSeparator();
    new CurrentMIGLWidgetAction("Ri&gid-Body Refine Current Atoms", "Rigid Body Refine the current atoms (cyan color)", refi_menu, SLOT(OnRefiRigidBody()), SLOT(OnUpdateRefiRigidBody(QAction*)));

    new CurrentMIGLWidgetAction("Fin&d Ligand Fit and Conformer", "Search for the best conformer ligand fit and torsion angles", refi_menu, SLOT(OnRefiLigandFit()), SLOT(OnUpdateRefiLigandFit(QAction*)));

    refi_menu->addSeparator();
    acceptRefineAction = new CurrentMIGLWidgetAction("&Accept Refine\tCtrl+Shift+R", "Accept the refinement and finalize the atoms positions", refi_menu, SLOT(OnRefiAccept()), SLOT(OnUpdateRefiAccept(QAction*)));

    new CurrentMIGLWidgetAction("&Reset Refine", "Reset the current real-space refinement by putting atoms back where they were", refi_menu, SLOT(OnRefiReset()), SLOT(OnUpdateRefiReset(QAction*)));

    new CurrentMIGLWidgetAction("&Cancel Refine", "Cancel the current real-space refinement by putting atoms back where they were", refi_menu, SLOT(OnRefiCancel()), SLOT(OnUpdateRefiCancel(QAction*)));

    refi_menu->addSeparator();
    new CurrentMIGLWidgetAction("&Undo Refine", "", refi_menu, SLOT(OnRefiUndo()), SLOT(OnUpdateRefiUndo(QAction*)));

    new CurrentMIGLWidgetAction("Re&do Refine", "", refi_menu, SLOT(OnRefiReDo()), SLOT(OnUpdateRefiRedo(QAction*)));

    refi_menu->addSeparator();
    new CurrentMIGLWidgetAction("Refine &Options", "Real-space refine options", refi_menu, SLOT(OnRefiOptions()));

    MIMenu *disorder_menu = new MIMenu(*this);
    disorder_menu->Append(ID_FIT_SPLITTORSION, "Split at &Torsion", "Adds disorder (A and B parts) at a torsion angle", false);
    disorder_menu->Append(ID_FIT_SPLITFIT, "Split &Fit", "duplicates residues being fit making an A and B part", false);

    MIMenu *pentamer_menu = new MIMenu(*this);
    pentamer_menu->Append(ID_FIT_MATCHPENTAMER, "&Suggest Backbone Match", "Find a pentamer that best matches the CA and CB atoms (if present)", false);
    pentamer_menu->Append(ID_FIT_REPLACEMIDDLE3, "&Replace Middle 3", "Replace the backbone atoms of the middle 3 residues", false);
    pentamer_menu->Append(ID_FIT_REPLACEFIRST4, "R&eplace First 4", "Replace the backbone atoms of the first 4 residues", false);
    pentamer_menu->Append(ID_FIT_REPLACELAST4, "Re&place Last 4", "Replace the backbone atoms of the last 4 residues", false);
    pentamer_menu->Append(ID_FIT_REPLACEALL, "Rep&lace All 5", "Replace the backbone atoms of all 5 residues", false);
    //  Removed the Clear backbone match menu entry because it failed to update the tree properly, leading to a crash
    pentamer_menu->Append(ID_FIT_CLEARPENTAMER, "&Clear Backbone Match", "Clear the backbone pentamer from the screen (and memory)", false);
    pentamer_menu->Append(ID_FIT_FLIPPEPTIDE, "&Flip Peptide", "Flip petide bond to put carbonyl on opposite side", false);

    MIMenu *fitsurf_menu = new MIMenu(*this);
    fitsurf_menu->Append(ID_FIT_SURFVDW, "&Van Der Waal Surface", "Surface the current atoms (fit) with a van der Waal surface", true);
    fitsurf_menu->Append(ID_FIT_SURFEXT, "&Extended Surface", "Surface the current atoms (fit) with an extended (2x) surface", true);
    fitsurf_menu->Append(ID_FIT_SURFPROBE, "&Contact Surface", "Surface the current atoms (fit) wherever there is a contact", true);
    fitsurf_menu->Append(ID_FIT_SURFNONE, "&No Surface", "Turn off the surface for the current atoms", true);

    fit_menu = new MIMenu(*this);
    fit_menu->Append(ID_FIT_RESIDUE, "Fit R&esidue\tF", "Set the last picked residue for fitting", true);
    fit_menu->Append(ID_FIT_SINGLEATOM, "Fit A&tom", "Set the last picked atom for fitting", true);
    fit_menu->Append(ID_FIT_RESIDUES, "Fit Re&sidues", "Set all the picked residues for fitting", true);
    fit_menu->Append(ID_FIT_RANGE, "Fit Res&idue Range", "Fit the range marked by the top two residues on the stack", true);
    fit_menu->Append(ID_FIT_ATOMS, "Fit At&oms", "Set all the picked atoms for fitting", true);
    fit_menu->Append(ID_FIT_FITMOLECULE, "Fit &Molecule", "Set the selected molecule for fitting", true);
    fit_menu->AppendSeparator();
    fit_menu->Append(ID_FIT_APPLY, "&Accept Fit\t;", "Accept fitting, leave atoms at new positions");
    fit_menu->Append(ID_FIT_RESET, "&Reset Fit", "Put atoms back to start and continue fitting");
    fit_menu->Append(ID_FIT_CANCEL, "&Cancel Fit", "Stop fitting and put atoms back to starting positions");
    fit_menu->AppendSeparator();
    fit_menu->Append(ID_FIT_ROTATE, "&Rotate", "Set right mouse to rotate mode", true);
    fit_menu->Append(ID_FIT_TRANSLATE, "&Translate", "Set right mouse to translate mode", true);
    fit_menu->Append(ID_FIT_TORSION, "&Torsion", "Set right mouse to torsion mode", true);
    fit_menu->Append(ID_FIT_CENTERMODE, "C&enter", "Set right mouse to move screen center (default mode)", true);
    fit_menu->AppendSeparator();
    fit_menu->Append(ID_FIT_SETUPTORSION, "&Set Up Torsion", "Torsion about last two atoms picked, second end moves");
    fit_menu->Append(ID_FIT_CLEARTORSION, "C&lear Torsion", "Clear Torsion, leave in fitting mode");
    fit_menu->AppendSeparator();
    fit_menu->Append(ID_FIT_SURFMENU, "Surface &Fit Atoms", fitsurf_menu);
    //fit_menu->AppendSeparator();
    fit_menu->Append(ID_FIT_BACKBONE_MENU, "Fix &Backbone", pentamer_menu);
    //fit_menu->AppendSeparator();
    fit_menu->Append(ID_FIT_BACKBONE_MENU, "&Disorder", disorder_menu);
    fit_menu->AppendSeparator();
    fit_menu->Append(ID_FIT_UNDO, "&Undo Fit", "Undo Last fit");
    fit_menu->Append(ID_FIT_REDO, "&Redo Fit", "Undo the Undo Last fit");
    fit_menu->AppendSeparator();
    fit_menu->Append(ID_FIT_LSQSUPERPOSE, "S&uperimpose...", "Superimpose models by LSQ fitting");

    model_menu = new MIMenu(*this);
    model_menu->Append(ID_OBJECT_NEWMODEL, "&New Model...", "Adds a new blank model", false);
    model_menu->AppendSeparator();
    model_menu->Append(ID_FIT_INSERTRESIDUE, "&Add residue", "Add or insert a new residue after the last pick", false);
    model_menu->Append(ID_FIT_REPLACEWITH, "&Replace residue", "Replace a residue with another type", false);
    model_menu->Append(ID_REPLACEANDFIT, "Replace and &Fit\tR", "Replace a residue with another type and search confomers", false);
    model_menu->Append(ID_FIT_REPLACESEQUENCE, "Replace with &Sequence", "Replace a residue with the lower sequence", false);
    model_menu->Append(ID_FIT_NEXTCONFOMER, "Next &Conformer", "Replace a residue with another type", false);
    model_menu->Append(ID_FIT_DELETERESIDUE, "&Delete residue\tD", "Delete the last picked residue", false);
    model_menu->Append(ID_FIT_RENAMERESIDUE, "R&ename residue", "Rename(renumber) the last picked residue", false);
    model_menu->Append(ID_DELETEATOM, "&Delete atom\tDelete", "Delete the last picked atom", false);
    model_menu->Append(ID_FIT_ADDWATER, "Add &Waters...", "Sprinkle water throughout map", false);
    model_menu->AppendSeparator();
    model_menu->Append(ID_GEOM_BOND, "&Bond", "Force a bond between the last two picked atoms", false);
    model_menu->Append(ID_GEOM_UNBOND, "B&reak Bond", "Break the bond between the last two picked atoms", false);
    model_menu->AppendSeparator();
    model_menu->Append(ID_FIT_MARKBEFORE, "Add MRK &Before\t<", "Adds a MRK (C-alpha) before the focus residue for tracing the chain", false);
    model_menu->Append(ID_FIT_MARKAFTER, "Add MRK &After\t>", "Adds a MRK (C-alpha) after the focus residue for tracing the chain", false);
    model_menu->Append(ID_FIT_POLYALA, "Poly-Ala &Range", "Turns the model into poly-alanine between the two picks", false);
    model_menu->Append(ID_FIT_POLYALACHAIN, "Poly-Ala &Chain", "Turns the chain containing the last pick into poly-alinine", false);

    model_menu->Append(ID_GOTO_NTER, "Go To &N-terminus\t[", "Follows the current chain to its N-terminal end", false);
    model_menu->Append(ID_GOTO_CTER, "Go To &C-terminus\t]", "Follows the current chain to its C-terminal end", false);
    model_menu->AppendSeparator();
    model_menu->Append(ID_MODEL_CHECKPOINT, "&Checkpoint Model", "Checkpoints the model - use to save before", false);
    model_menu->Append(ID_MODEL_REVERT, "&Revert Model...", "Reverts the model to an earlier version saved", false);
    model_menu->Append(ID_MODEL_AUTOCHECKPOINT, "Auto-checkpoint Model", "Automatically checkpoints the model every few minutes", true);


    menu_bar->Append(file_menu, "&File");

    QMenu *di_menu = menuBar()->addMenu(tr("&Dictionary"));
    di_menu->addAction(tr("Load &New Dictionary..."), this, SLOT(OnLoadDictReplace()));
    di_menu->addAction(tr("Load & &Append Dictionary..."), this, SLOT(OnLoadDictAppend()));

    QMenu *ligImportMenu = di_menu->addMenu("&Import Ligand");
    ligImportMenu->addAction(tr("&CIF"), this, SLOT(OnLoadLigCif()));
    ligImportMenu->addAction(tr("&Mol"), this, SLOT(OnLoadLigMol()));
    ligImportMenu->addAction(tr("&Pdb"), this, SLOT(OnLoadLigPdb()));
    ligImportMenu->addAction(tr("&Smiles"), this, SLOT(OnLoadLigSmi()));

    di_menu->addAction(tr("&Save Dictionary..."), this, SLOT(OnSaveDict()));
    di_menu->addAction(tr("&Edit Residue..."), this, SLOT(OnEditDictResidue()));

    menu_bar->Append(view_menu, "&Viewpoint");
    menu_bar->Append(show_menu, "&Show");
    render_menu->setTitle("Re&nder");
    menu_bar->Append(render_menu);
    menu_bar->Append(model_menu, "&Model");
    menu_bar->Append(fit_menu, "F&it");
    refi_menu->setTitle("&Refine");
    menu_bar->Append(refi_menu);
    menu_bar->Append(analyze_menu, "&Analyze");

    JobManager = new BatchJobManager();
    _jobMenu = menuBar()->addMenu("&Job");
    QTimer::singleShot(0, this, SLOT(fillJobMenu()));


    hide_menu = new MIMenu(*this);
    hide_menu->Append(ID_HIDEMODEL, "Whole Model", "Hide the whole model", false);
    hide_menu->Append(ID_OBJECT_RESIDUE_TURNOFF, "Hide Last picked Residue", "Hide the last picked residue", false);
    hide_menu->Append(ID_OBJECT_RESIDUES_TURNOFF, "Hide All Picked Residues", "Hide all the residues on the stack", false);
    hide_menu->Append(ID_OBJECT_RESIDUERANGE_TURNOFF, "Hide Residue Range", "Hide the residue at the top of stack", false);
    hide_menu->Append(ID_SHOW_PICKEDATOM_TURNOFF, "Hide Last Picked Atom", "Hide the last picked atom", false);
    hide_menu->Append(ID_SHOW_ALLPICKEDATOMS_TURNOFF, "Hide All Picked Atoms", "Hide all the atoms on the stack", false);




    windowMenu = menuBar()->addMenu(tr("&Window"));
    updateWindowMenu();
    connect(windowMenu, SIGNAL(aboutToShow()), this, SLOT(updateWindowMenu()));

    viewMenu = menuBar()->addMenu(tr("View"));

    QMenu *help_menu = menuBar()->addMenu(tr("&Help"));
    help_menu->addAction(tr("&Help..."), this, SLOT(OnHelp()));
    help_menu->addAction(tr("&About..."), this, SLOT(OnAbout()));
    help_menu->addAction(tr("&OpenGL format..."), this, SLOT(OnGLFormat()));

    initMenuHandlers();
}


#include <images/annotate.xpm>
#include <images/apply.xpm>
#include <images/cancel.xpm>
#include <images/center.xpm>
#include <images/chekpoin.xpm>
#include <images/colortool.xpm>
#include <images/copy.xpm>
#include <images/decrper.xpm>
#include <images/fit.xpm>
#include <images/hidetool.xpm>
#include <images/incrper.xpm>
#include <images/new.xpm>
#include <images/open.xpm>
#include <images/print.xpm>
#include <images/rotate.xpm>
#include <images/roty90.xpm>
#include <images/rotym90.xpm>
#include <images/save.xpm>
#include <images/showtool.xpm>
#include <images/sidetool.xpm>
#include <images/slabin.xpm>
#include <images/slabout.xpm>
#include <images/topview.xpm>
#include <images/torsion.xpm>
#include <images/translate.xpm>
#include <images/zoomin.xpm>
#include <images/zoomout.xpm>
//#include <images/annotateDisabled.xpm>
//#include <images/applyDisabled.xpm>
//#include <images/cancelDisabled.xpm>
//#include <images/centerDisabled.xpm>
//#include <images/chekpoinDisabled.xpm>
//#include <images/colortoolDisabled.xpm>
//#include <images/copyDisabled.xpm>
//#include <images/decrperDisabled.xpm>
//#include <images/fitDisabled.xpm>
//#include <images/hidetoolDisabled.xpm>
//#include <images/incrperDisabled.xpm>
//#include <images/newDisabled.xpm>
//#include <images/openDisabled.xpm>
//#include <images/printDisabled.xpm>
//#include <images/rotateDisabled.xpm>
//#include <images/roty90Disabled.xpm>
//#include <images/rotym90Disabled.xpm>
//#include <images/saveDisabled.xpm>
//#include <images/showtoolDisabled.xpm>
//#include <images/sidetoolDisabled.xpm>
//#include <images/slabinDisabled.xpm>
//#include <images/slaboutDisabled.xpm>
//#include <images/topviewDisabled.xpm>
//#include <images/torsionDisabled.xpm>
//#include <images/translateDisabled.xpm>
//#include <images/zoominDisabled.xpm>
//#include <images/zoomoutDisabled.xpm>
void MIMainWindow::createToolBars()
{
    tool_bar = new MIToolBar(menu_bar, this);

    //TODO: use disabled icon specialization for QIcons?
    fileNewAction->setIcon(QIcon(QPixmap(new_xpm)));
    tool_bar->toolBar()->addAction(fileNewAction);
    fileOpenAction->setIcon(QIcon(QPixmap(open_xpm)));
    tool_bar->toolBar()->addAction(fileOpenAction);
    fileSaveAction->setIcon(QIcon(QPixmap(save_xpm)));
    tool_bar->toolBar()->addAction(fileSaveAction);
    filePrintAction->setIcon(QIcon(QPixmap(print_xpm)));
    tool_bar->toolBar()->addAction(filePrintAction);
    copyCanvasAction->setIcon(QIcon(QPixmap(copy_xpm)));
    tool_bar->toolBar()->addAction(copyCanvasAction);

    tool_bar->AddSeparator();
    tool_bar->AddTool(ID_OBJECT_ANNOTATION, annotate_xpm);
    tool_bar->AddTool(ID_MODEL_CHECKPOINT, chekpoin_xpm);
    tool_bar->AddTool(ID_VIEW_SLABIN, slabin_xpm);
    tool_bar->AddTool(ID_VIEW_SLABOUT, slabout_xpm);
    tool_bar->AddTool(ID_GOTO_ZOOMIIN, zoomin_xpm);
    tool_bar->AddTool(ID_GOTO_ZOOMOUT, zoomout_xpm);
    tool_bar->AddTool(ID_VIEW_TOPVIEW, topview_xpm);
    tool_bar->AddTool(ID_VIEW_ROTATEY90, roty90_xpm);
    tool_bar->AddTool(ID_VIEW_ROTATEYMINUS, rotym90_xpm);
    tool_bar->AddTool(ID_DECREASE_PERSP, decrper_xpm);
    tool_bar->AddTool(ID_INCREASE_PERSP, incrper_xpm);

    tool_bar->AddSeparator();
    tool_bar->AddTool(ID_FIT_RESIDUE, fit_xpm);
    tool_bar->AddTool(ID_FIT_APPLY, apply_xpm);
    tool_bar->AddTool(ID_FIT_CANCEL, cancel_xpm);
    tool_bar->AddTool(ID_FIT_TRANSLATE, translate_xpm);
    tool_bar->AddTool(ID_FIT_ROTATE, rotate_xpm);
    tool_bar->AddTool(ID_FIT_TORSION, torsion_xpm);
    tool_bar->AddTool(ID_FIT_CENTERMODE, center_xpm);

    QToolBar *toolBar = addToolBar("Display Tools");
    toolBar->setObjectName("Display Tools");
    toolBar->addAction(QIcon(colortool_xpm), tr("Color"), this, SLOT(OnColorTool()));
    toolBar->addAction(QIcon(showtool_xpm), tr("Show"), this, SLOT(OnShowTool()));
    toolBar->addAction(QIcon(sidetool_xpm), tr("Side Chain"), this, SLOT(OnSideChainTool()));
    toolBar->addAction(QIcon(hidetool_xpm), tr("Hide"), this, SLOT(OnHideTool()));

    tool_bar->show();
}


void MIMainWindow::createStatusBar()
{
    statusBar()->showMessage(tr("Ready"));
}

MIGLWidget*MIMainWindow::currentMIGLWidget()
{
    QMdiSubWindow *activeSubWindow = mdiArea->currentSubWindow();
    if (activeSubWindow)
    {
        return dynamic_cast<MIGLWidget*>(activeSubWindow->widget());
    }
    return 0;
}

QMdiArea*MIMainWindow::getMdiArea()
{
    return mdiArea;
}

void MIMainWindow::setActiveMIGLWidget(MIGLWidget *w)
{
    foreach (QMdiSubWindow *window, mdiArea->subWindowList())
    {
        if (window->widget() == w)
            mdiArea->setActiveSubWindow(window);
    }
}


QMdiSubWindow*MIMainWindow::findMIGLWidget(const QString &fileName)
{
    QString canonicalFilePath = QFileInfo(fileName).canonicalFilePath();

    foreach (QMdiSubWindow *window, mdiArea->subWindowList())
    {
        MIGLWidget *mdiChild = dynamic_cast<MIGLWidget*>(window->widget());
        if (mdiChild && QString(mdiChild->GetFilename().c_str()) == canonicalFilePath)
            return window;
    }
    return 0;
}

void MIMainWindow::switchLayoutDirection()
{
    if (layoutDirection() == Qt::LeftToRight)
        qApp->setLayoutDirection(Qt::RightToLeft);
    else
        qApp->setLayoutDirection(Qt::LeftToRight);
}

void MIMainWindow::setActiveSubWindow(QWidget *window)
{
    if (!window)
        return;
    mdiArea->setActiveSubWindow(qobject_cast<QMdiSubWindow*>(window));
}

void MIMainWindow::OnNew()
{
    newWindow();
}


void MIMainWindow::openRecentFile()
{
    QAction *action = qobject_cast<QAction*>(sender());
    if (action)
    {
        std::vector<std::string> fnames;
        fnames.push_back(action->data().toString().toStdString());
        OpenFiles(fnames);
    }
}

void MIMainWindow::addRecentFileActions(QMenu *fileMenu)
{
    fileMenu->addSeparator();
    for (int i = 0; i < MaxRecentFiles; ++i)
        fileMenu->addAction(recentFileActs[i]);
    fileMenu->addSeparator();
}

void MIMainWindow::setCurrentFile(const std::string &fname)
{
    QString fileName(fname.c_str());

    QSettings *settings = MIGetQSettings(); // could use MIConfig (it's the same file), but this api is easier here
    QStringList files = settings->value("recentFileList").toStringList();
    files.removeAll(fileName);
    files.prepend(fileName);
    while (files.size() > MaxRecentFiles)
        files.removeLast();
    settings->setValue("recentFileList", files);
    updateRecentFileActions();
}

void MIMainWindow::updateRecentFileActions()
{
    QSettings *settings = MIGetQSettings(); // could use MIConfig (it's the same file), but this api is easier here
    QStringList files = settings->value("recentFileList").toStringList();

    int numRecentFiles = qMin(files.size(), (int)MaxRecentFiles);

    for (int i = 0; i < numRecentFiles; ++i)
    {
        QString text = QString("&%1 %2").arg(i + 1).arg(files[i]);
        recentFileActs[i]->setText(text);
        recentFileActs[i]->setData(files[i]);
        recentFileActs[i]->setVisible(true);
    }
    for (int j = numRecentFiles; j < MaxRecentFiles; ++j)
        recentFileActs[j]->setVisible(false);
}

RamaPlotMgr*MIMainWindow::RamaPlotManager()
{
    return ramaPlotMgr;
}

void MIMainWindow::HasCurrentMIGLWidget(const MIUpdateEvent &pCmdUI)
{
    pCmdUI.Enable(currentMIGLWidget() != 0);
}

//some shortcuts for MIMainWindow::instance()->Log, etc
void MIMainWindowLog(const std::string &str)
{
    return MIMainWindow::instance()->Log(str);
}

void MIMainWindowDebug(const std::string &str)
{
    return MIMainWindow::instance()->Debug(str);
}

void MIMainWindowLeftFooter(const std::string &str, int timeout)
{
    MIMainWindow::instance()->LeftFooter(str, timeout);
}

void MIMainWindowMiddleFooter(const std::string &str)
{
    MIMainWindow::instance()->MiddleFooter(str);
}

void MIMainWindowRightFooter(const std::string &str)
{
    MIMainWindow::instance()->RightFooter(str);
}

void MIMainWindow::updateShowMenu()
{
    if (solidSurfMenuAction)
    {
        show_menu->removeAction(solidSurfMenuAction);
        solidSurfMenuAction = NULL;
    }
    MIGLWidget *view = currentMIGLWidget();
    if (view)
    {
        solidSurfMenuAction = view->solidSurfMenuAction();
    }
    else
    {
        solidSurfMenuAction = new QAction(this);
        solidSurfMenuAction->setText(tr("Solid surface"));
    }
    show_menu->insertAction(canvas_menu->menuAction(), solidSurfMenuAction);
}

void MIMainWindow::updateIsRefining(bool isRefining)
{
    refineResidueAction->setEnabled(!isRefining);
    acceptRefineAction->setEnabled(isRefining);
}

void MIMainWindow::setRenderLineThickness(int thickness)
{
    if (currentMIGLWidget())
        currentMIGLWidget()->setViewpointLineThickness(thickness);
}

void MIMainWindow::updateRenderMenu()
{
    bool lineThicknessEnabled = false;
    if (currentMIGLWidget())
    {
        ViewPoint *viewpoint = currentMIGLWidget()->viewpoint;
        int thickness = viewpoint->GetLineThickness();
        if (thickness > 0 && renderLineThicknessMenu_->actions().size() <= thickness)
        renderLineThicknessMenu_->actions().at(thickness-1)->setChecked(true);
        lineThicknessEnabled = viewpoint->GetBallandStick() == ViewPoint::BALLANDSTICK
                  || viewpoint->GetBallandStick() == ViewPoint::STICKS;
    }
    renderLineThicknessMenu_->setEnabled(lineThicknessEnabled);
}

void MIMainWindow::fillJobMenu()
{
    createLocalSocketScript();
    Tools::instance().FillToolsMenu(_jobMenu);
}

QMenu *MIMainWindow::jobMenu() const
{
    return _jobMenu;
}

void MIMainWindow::saveJobMenu()
{
    GetJobManager()->saveJobMenu(_jobMenu);
}

void MIMainWindow::createLocalSocketScript()
{
    _scriptSocket = new LocalSocketScript(this);
}

QString MIMainWindow::scriptPort() const
{
    if (_scriptSocket)
        return _scriptSocket->name();
    return QString::null;
}

void MIMainWindow::addJob(const QString &menuName, const QString &jobName, const QString &executable, const QStringList &arguments, const QString &workingDirectory)
{
    QAction *jobAction = GetJobManager()->customJobAction(menuName, jobName, executable, arguments, workingDirectory);
    jobMenu()->addAction(jobAction);
}

void MIMainWindow::runPythonScript(const std::string &file)
{
    BatchJob *pythonJob = GetJobManager()->CreateJob();
    pythonJob->setProgram(BatchJobManager::pythonExe());
    pythonJob->setArguments(file.c_str());

    pythonJob->StartJob();
}
