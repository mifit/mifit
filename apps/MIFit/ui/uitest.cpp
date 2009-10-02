#include <nongui/nonguilib.h>
#include <chemlib/chemlib.h>
#include <chemlib/RESIDUE_.h>
#include "core/corelib.h"


#include "MIMenuBar.h"
#include "MIMainWindow.h"
#include "MIGLWidget.h"

#include "Displaylist.h"
#include "uitest.h"
#include "ui/MIDialog.h"
#include "id.h"

#include "ViewSyncedPanel.h"

using namespace chemlib;

//#define MY_PRINTF printf
#define MY_PRINTF Logger::log


void PrintTabs(unsigned int level) {
  for (unsigned int i = 0; i < level; ++i) {
    MY_PRINTF("\t");
  }
}

static std::vector<unsigned int> SKIPPED_MENU_ITEMS;

// these menu items are excluded from the random tests
static void SetupSkippedItems(bool skip_new_model=true, bool skip_clear_stack=false)
{
  SKIPPED_MENU_ITEMS.clear();

  // FIXME: not yet converted to MIDialog
  SKIPPED_MENU_ITEMS.push_back(ID_FIT_LSQSUPERPOSE); //"Superimpose..."
  SKIPPED_MENU_ITEMS.push_back(ID_PREFERENCES); // "Preferences..."
  SKIPPED_MENU_ITEMS.push_back(ID_MANAGECRYSTALS); //"Manage Crystals..."
  SKIPPED_MENU_ITEMS.push_back(ID_OBJECT_DEFINECOLORS); // "Define Atom Colors..."
  SKIPPED_MENU_ITEMS.push_back(ID_BVALUECOLORS); //"B-Value Colors Ranges..."

  // these work, but are slow down rendering too much
  SKIPPED_MENU_ITEMS.push_back(ID_RENDER_SPACEFILLING); // "CPK"
  SKIPPED_MENU_ITEMS.push_back(ID_RENDERING_BALLANDSTICK); // "Knob and Stick"
  SKIPPED_MENU_ITEMS.push_back(ID_RENDER_BALLANDCYLINDER); // "Ball and Cylinder"

  // something about the event loop/timeout of these prevents any other command from
  // being called in debug mode
  SKIPPED_MENU_ITEMS.push_back(ID_ANIMATE_ROCK); // ROCK
  SKIPPED_MENU_ITEMS.push_back(ID_ANIMATE_ROLL); // ROLL

  // bad idea...
  SKIPPED_MENU_ITEMS.push_back(ID_PRINT); // "Print..."
  SKIPPED_MENU_ITEMS.push_back(DOCVIEW_HELP); //"Help..."

  // don't test debug items: could interfere with testing (and prevents history logging)
  SKIPPED_MENU_ITEMS.push_back(ID_DO_RANDOM_TEST); // "Debug: random tests"
  SKIPPED_MENU_ITEMS.push_back(ID_PLAY_HISTORY); // "Debug: play history"
  SKIPPED_MENU_ITEMS.push_back(ID_RECORD_HISTORY); // "Debug: record history"
  SKIPPED_MENU_ITEMS.push_back(ID_STOPRECORD_HISTORY) ; //"Debug: stop recording history"

  // makes following what's happening unpleasant
  SKIPPED_MENU_ITEMS.push_back(ID_VIEW_FULLSCREEN); //"Toggle Fullscreen"

  // well, duh
  SKIPPED_MENU_ITEMS.push_back(ID_FILE_EXIT); // "Exit"

  // Changes wxView, so testing stops working
  SKIPPED_MENU_ITEMS.push_back(ID_SITEPLOT); // "Export Active Site Plot"

  // possibly opens/closes windows.  
  // new windows don't have a model or map, so testing them is pointless
  // closed windows would be silly to test, too
  SKIPPED_MENU_ITEMS.push_back(ID_FILE_CLOSE); // "Close"
  if (skip_new_model) {
    SKIPPED_MENU_ITEMS.push_back(ID_FILE_NEW); //"New"
    SKIPPED_MENU_ITEMS.push_back(ID_OBJECT_NEWMODEL); // "New Model..."
    SKIPPED_MENU_ITEMS.push_back(ID_FILE_OPEN_NEW); // "Open into new document..."
  }
  if (skip_clear_stack)
    SKIPPED_MENU_ITEMS.push_back(ID_OBJECT_CLEARSTACK); //"Clear stack"
}


static void DoStackPush(MIGLWidget* view, MIAtom* atom, RESIDUE* res, Molecule* mol, bool recenter = false) {
}

static void PushResidueOrAtom(Molecule* mol, int resnum, MIGLWidget* view, bool random_atom = false, bool do_select=true) {
  RESIDUE* res = mol->getResidues();
  for (int i = 0; i < resnum; ++i) {
    res = res->next();
  }
  if (!random_atom) {
    if(do_select)
		view->select(mol, res, res->atom(0));
    DoStackPush(view, res->atom(0), res, mol);
  } else {
    MIAtom* a = res->atom(GetRand("Random atom", res->atomCount()));
	if (do_select)
		view->select(mol, res, a);
    DoStackPush(view, a, res, mol);
  }
}

static int GetResidueCount() {
  MIGLWidget* view = MIMainWindow::instance()->currentMIGLWidget();
  if (!view)
    return 0;
  Displaylist* dl = view->GetDisplaylist();
  Molecule* mol = dl->GetCurrentModel();
  if (!mol) {
    return 0;
  }

  return mol->getnresidues();
}

static bool ValidateMolecule() {
  MIGLWidget* view = MIMainWindow::instance()->currentMIGLWidget();
  if (!view)
    return true;
  Displaylist* dl = view->GetDisplaylist();
  Molecule* mol = dl->GetCurrentModel();
  if (!mol) {
    return true;
  }

  MIAtom_const_iter atom, endAtom;
  for (MIIter<RESIDUE> res = mol->GetResidues(); res; ++res) {
    RESIDUE& r = *res;
    const MIAtomList& atoms = r.atoms();
    endAtom = atoms.end();
    for (atom = atoms.begin(); atom != endAtom; ++atom) {
      MIAtom& a = **atom;
      if (MIIsNan(a.x()) || MIIsNan(a.y()) || MIIsNan(a.z()))
        return false;
    }
  }
  return true;
}

// limit the active model to 750 residues
void LimitGrowth() {
  MIGLWidget* view = MIMainWindow::instance()->currentMIGLWidget();
  if (!view)
    return;
  Displaylist* dl = view->GetDisplaylist();
  Molecule* mol = dl->GetCurrentModel();
  if (!mol) {
    return;
  }

  const int MAX_RES_COUNT=750;

  if (mol->getnresidues() > MAX_RES_COUNT) {

    MIIter<RESIDUE> res = mol->GetResidues();
    int count=0;
    for (; (bool)res && count < MAX_RES_COUNT; ++res, ++count) {}
    
    std::vector<RESIDUE*> deaders;
    for (; res; ++res)
      deaders.push_back(res);
    
    mol->DeleteResidues(deaders);
  }
}

static bool ValidateStackItem(const StackItem &item, Displaylist *dl) {
  MIAtom *a = item.atom;
  RESIDUE *r = item.residue;
  Molecule *m = item.molecule;
  if (!m || !r || !a)
    return false;
  
  bool found=false;
  for (Displaylist::ModelList::iterator modelIter = dl->begin();
       modelIter != dl->end(); ++modelIter) {
    if (m == *modelIter) {
      found=true;
      break;
    }
  }
  if (!found)
    return false;


  if (!m->Contains(r))
  {
    // check symmetry residues, too
    bool found=false;
    for (MIIter<RESIDUE> res=m->GetSymmResidues(); Residue::isValid(res); ++res) {
      if (r == res) {
        found=true;
        break;
      }
    }
    if (!found)
      return false;
  }

  for (int i=0; i < r->atomCount(); ++i)
    if (a == r->atom(i))
      return true;

  return false;
}


//validate that all items (if any) in the stack are in the display list, 
// and that the atom, residue and molecule in the stack are consistent
static bool ValidateStack(MIGLWidget *view) {
    // TODO remove
    return true;
}


static void CreateRandomStack() {
  // TODO remove
}



static void DoRandomEvents(MIMenuBar *mb, unsigned int iterations) {

    // TODO remove
}


static void ExhaustiveMenuTest(MIMenuBar* mb, unsigned int iterations) {

    // TODO remove
}


static bool IN_TEST_MODE=false;
bool IsInTestMode() {
  return IN_TEST_MODE;
}


void TestMode() {
    // TODO remove
}
