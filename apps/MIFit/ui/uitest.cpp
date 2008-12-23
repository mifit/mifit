#include "nonguilib.h"
#include "chemlib.h"
#include "RESIDUE_.h"
#include "corelib.h"


#include "MIMenuBar.h"
#include "MIMainWindow.h"
#include "MIGLWidget.h"

#include "Displaylist.h"
#include "uitest.h"
#include "MIDialog.h"
#include "MIHistory.h"
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

#ifdef USE_ASPLOT
  // Changes wxView, so testing stops working
  SKIPPED_MENU_ITEMS.push_back(ID_SITEPLOT); // "Export Active Site Plot"
#endif

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
  view->AtomStack->Push(atom, res, mol);
  Displaylist* dl = view->GetDisplaylist();
  MIGetHistory()->AddAtomPick("rand", dl, atom, res, mol, recenter);
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
    MIData values;
    values["command"].str="prune";
    MIGetHistory()->AddCommand(values);

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
  if (!view)
    return false;
  Displaylist* dl = view->GetDisplaylist();
  if (!dl)
    return false;

  Stack *stack=view->AtomStack;
  if (!stack->size() || !stack->StackChanged())
    return true;
  stack->ClearChanged();


  const Stack::DataContainer &data=stack->getData();
  Stack::DataContainer::const_iterator iter = data.begin();
  while (iter != data.end()) {
    StackItem item = *iter;
    if (!ValidateStackItem(item, dl))
      return false;
    ++iter;
  }

  return true;
}


static void CreateRandomStack() {
  MIGLWidget* view = MIMainWindow::instance()->currentMIGLWidget();
  if (!view)
    return;
  Displaylist* dl = view->GetDisplaylist();
  Molecule* mol = dl->GetCurrentModel();
  if (!mol) {
    return;
  }

  int nres = mol->getnresidues();
  if (!nres) {
    return;
  }

  view->AtomStack->Clear();
  long mode = GetRand("Random Stack Mode", 5);
  switch (mode) {
    case 0: // one residue
      {
        int resnum = GetRand("Random Residue", nres);
        PushResidueOrAtom(mol, resnum, view);
        MY_PRINTF("Pushed single residue onto stack\n");
        break;
      }
    case 1: // residue range
      {
        int resnum2 = GetRand("Random Residue 1", nres);
        int resnum1 = GetRand("Random Residue 2", resnum2);
        RESIDUE* res = mol->getResidues();
        for (int i = 0; i < resnum1; ++i) {
          if (res) {
            res = res->next();
          }
        }
        RESIDUE* res2 = res;
        for (int i = resnum1; i < resnum2; ++i) {
          if (res2) {
            res2 = res2->next();
          }
        }

        if (res && res2) {
          DoStackPush(view, res->atom(0), res, mol);
          DoStackPush(view, res2->atom(0), res2, mol);
          view->AtomStack->ExpandTop2AllAtoms();
          MIData values;
          values["command"].str = "rand";
          values["type"].str = "stack_range";
          MIGetHistory()->AddCommand(values);
        }
        MY_PRINTF("Pushed %d consecutive residues onto stack\n", resnum2-resnum1);
        break;
      }
    case 2:  // 4 atoms (for torsion testing)
      {
        int resnum = GetRand("Random Residue Torsion", nres);
        RESIDUE* res = mol->getResidues();
        for (int i = 0; i < resnum; ++i) {
          res = res->next();
        }
        MIAtom* N = atom_from_name("N", *res);
        MIAtom* CA = atom_from_name("CA", *res);
        MIAtom* C = atom_from_name("C", *res);
        MIAtom* O = atom_from_name("O", *res);
        if (N && CA && C && O) {
          DoStackPush(view, N, res, mol);
          DoStackPush(view, CA, res, mol);
          DoStackPush(view, C, res, mol);
          DoStackPush(view, O, res, mol);
          MY_PRINTF("Pushed torsion (consecutive) atoms onto stack\n");
        }
      }
      break;

    case 3: // random number of scattered residues
      {
        long push_num = GetRand("Random Residue Count", (nres/10));
        for (int i = 0; i < push_num; ++i) {
          int resnum = GetRand("Random Residue", nres);
          PushResidueOrAtom(mol, resnum, view, false, i==push_num-1);
        }
        MY_PRINTF("Pushed %d random residues onto stack\n", push_num);
        break;
      }
    case 4: // random number of scattered atoms
      {
        long push_num = GetRand("Random Atom Count", (nres/10));
        for (int i = 0; i < push_num; ++i) {
          int resnum = GetRand("Random Residue Number", nres);
          PushResidueOrAtom(mol, resnum, view, true, i==push_num-1);
        }
        MY_PRINTF("Pushed %d random atoms onto stack\n", push_num);
        break;
      }
  }
}



static void DoRandomEvents(MIMenuBar *mb, unsigned int iterations) {

  // first 3/4 of iterations don't allow creation of new models or new (empty) views.
  SetupSkippedItems(true);
  MIGLWidget* view = MIMainWindow::instance()->currentMIGLWidget();

  bool rebuilt=false;
  WaitCursor wait("Random testing");
  for (unsigned int i = 0; i < iterations; ++i) {
    if (wait.CheckForAbort()) {
      return;
    }

//#if 0
    if (!rebuilt &&  i > 3*iterations/4 ) {
      // rebuild skipped items, allowing creation of new models/empty windows
      rebuilt=true;
      SetupSkippedItems(false);
    }
//#endif

    // do one random menu item
    std::string s=::format("Iteration %d of %d",i,iterations);
    MIGetHistory()->AddComment(s);
    MIMainWindowMiddleFooter(s);
    mb->DoRandomItem(SKIPPED_MENU_ITEMS);
    MY_PRINTF("%s\n",s.c_str());

    if (!ValidateStack(view)) {
      Logger::log("Inconsistent stack!\n");
    }
    MIHistory::ProcessGUIEvents();

    if (GetRand("Create Random Stack?", 10) == 0) {
      CreateRandomStack();
    }

#ifdef MI_USE_TREE
    ViewSyncedPanel* panel = MIMainWindow::instance()->GetModelsTree();
    if (panel) {
      panel->RandomTest();
    }
#endif

    LimitGrowth();
    if (!ValidateStack(view)) {
      Logger::log("Inconsistent stack!\n");
    }
    //if (!ValidateMolecule()) 
    //Logger::log("Molecule is now invalid at iteration %d\n",i);

    MIHistory::ProcessGUIEvents();
  }
}


static void ExhaustiveMenuTest(MIMenuBar* mb, unsigned int iterations) {

  for (unsigned int j = 0; j < iterations; ++j) {
    for (unsigned int cycle=0; cycle < 4; ++cycle) {

      bool skip_new_model=j<iterations/2;
      bool skip_clear_stack;
      bool reverse=false;
      switch(cycle)
      {
        case 0:
          skip_clear_stack=true;
          reverse=false;
          break;
        case 1:
          skip_clear_stack=true;
          reverse=true;
          break;
        case 2:
          skip_clear_stack=false;
          reverse=false;
          break;
        case 3:
          skip_clear_stack=false;
          reverse=true;
          break;
      }

      //for the first half of the iterations, we disallow creating a new, empty model
      SetupSkippedItems(skip_new_model, skip_clear_stack);
      
      MIMainWindowMiddleFooter(::format("Iteration %d/%d, cycle %d/%d reverse: %d",j,iterations,cycle,4,reverse));
      mb->DoExhaustiveTest(SKIPPED_MENU_ITEMS,reverse);
      
      MIHistory::ProcessGUIEvents();
      if (GetRand("Create Random Stack?", 10) == 0) {
        CreateRandomStack();
      }
      MIHistory::ProcessGUIEvents();

    }
  }
}


static bool IN_TEST_MODE=false;
bool IsInTestMode() {
  return IN_TEST_MODE;
}


void TestMode() {
  unsigned int iterations = 10000;
  if (IN_TEST_MODE) {
    MY_PRINTF("Can't recursively do random debugging!\n");
    return;
  }

  MIGenericDialog dlg(0,"Random testing");
  MIData values;
  values["exhaustive"].b = false;
  values["iterations"].u = iterations;
  values["seed"].u = (unsigned int)time(0);

  dlg.order("exhaustive");
  dlg.order("iterations");
  dlg.order("seed");

  dlg.label("exhaustive","Exhaustive menu tests?");
  dlg.label("iterations","Number of iterations (10 max for exhaustive)");
  dlg.label("seed","Random seed");

  if (!dlg.GetResults(values)) {
    return;
  }

  bool exhaustive=values["exhaustive"].b;
  iterations=values["iterations"].u;
  unsigned int seed=values["seed"].u;

  if (exhaustive && iterations > 10)
    iterations=10;
  
#ifdef _WIN32
  srand(seed);
#else
  srandom(seed);
#endif

  // we put this after the dialog to swallow the dialog response in
  // playback mode
  if (MIGetHistory()->IsPlaying() && !MIGetHistory()->IsPlaybackPaused()) {
    return;
  }

  // NOTE: we might like to copy the current mihistory.mih file here into
  // another file, and then use that other file to as the history file to
  // record the random commands.  Unfortunately, because this program
  // changes cwd at will, we don't really know where mihistory.mih is, and
  // if we opened a new file now it might be in a different directory
  // anyway.
  //
  // In any case, the motivation for a copy was to avoid having to copy
  // mihistory.mih ourselves in the event of a crash from the random
  // tests.  We'll have to move the file anyway in the event of a
  // non-random-testing crash, so we'll be used to that.  Therefore, the
  // random commands will just go into the normal history log file.

  IN_TEST_MODE = true;
  MIDialog::SetRandom(true);
  if (exhaustive) {
    ExhaustiveMenuTest(MIMainWindow::instance()->getMenuBar(), iterations);
  } else {
    DoRandomEvents(MIMainWindow::instance()->getMenuBar(), iterations);
  }
  MIDialog::SetRandom(false);
  IN_TEST_MODE = false;


  MIMessageBox("FINISHED!","Finished all testing iterations!", MIDIALOG_ICON_INFORMATION);
}
