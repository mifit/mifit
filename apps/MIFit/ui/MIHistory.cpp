#include <cstdarg>

#include "nonguilib.h"
#include "utillib.h"
#include "Application.h"
#include "RESIDUE_.h"
#include "molw.h"
#include "MIHistory.h"
#include "Displaylist.h"

#include "ViewSyncedPanel.h"
#include "uitest.h"

#include "DictEditCanvas.h"
#include "DictEditDialog.h"

#include "MIMenuBar.h"
#include "MIMenu.h"
#include "MIMainWindow.h"
#include "MIGLWidget.h"
#include <QCoreApplication>
#include <QApplication>

using namespace chemlib;

//
//  History format
//
//  command{data}
//    where command is one of the following and {data} is
//    formatted MIData data
//
//  Command:  Data keys:
//  menu      menu=str(menu description string)
//  pick      data is has a special data["type"].str which is one of
//    clear
//    hide
//    pop
//    center
//    bond      atom1 res1 mol1 atom2 res2 mol2
//    atom      atom res mol recenter
//  dialog    data is MIData for the dialog
//  key       key flags
//  mouse     data is has a special data["type"].str which is one of
//    slab      frontclip backclip        (view slab)
//    zoom      zoom                      (view zoom)
//    rot       xang yang zang            (view rotate)
//    frot      xang yang zang            (fitting rotate)
//    tran      x y z                     (fitting translate)
//    tors      ang                       (fitting torsion)
//    pan       x y z                     (view translate)
//    move      x y width height          (for mouse tracking??)
//  rand      from random testing
//    stack_range                         (expand top 2 stack items to residue range)
//    atom      atom res mol recenter     (push item on stack)
//  tree      data is has a special data["type"].str which is one of
//    atoms              sel (hex-encoded selection for tree)
//    residues           sel
//    models             sel
//    atom_selected      sync atom res mol  
//    atom_activated     sync atom res mol  
//    residue_selected   sync atom res mol  
//    residue_activated  sync atom res mol  
//    model_selected     sync atom res mol  
//    chain_selected     sync atom res mol  
//    chain_activated    sync atom res mol  
//    emap_selected      map
//  seqw      data is has a special data["type"].str which is one of
//    select    atom res mol ss           (select residue)
//    change    atom res mol ss           (change ss type)
//    focus     atom res mol ss           (change focus)
//  rama      data is has a special data["type"].str which is one of
//    PHI       inc                       (alter phi of active fitting residue)
//    PSI       inc                       (alter psi of active fitting residue)
//    focus     atom res mol              (change focus)
//  frame     index name
//  dict  ok cancel

MIHistory::MIHistory()
  : HISTOUT(NULL), recording_paused(0), playback_paused(0) {
  _current_frame = 0;
  AddFrame(0);
}

MIHistory::~MIHistory() {
  for (size_t i = 0; i < HISTIN.size(); ++i) {
    fclose(HISTIN[i]);
  }
  StopRecording();
}

bool MIHistory::Record(const std::string& fname) {
  StopRecording();
  HISTOUT = fopen(fname.c_str(), "w");
  if (IsRecording()) {
    AddComment(std::string("MIFit version ") + MIFit_version);
  }
  return IsRecording();
}

bool MIHistory::Play(const std::string& fname) {
  FILE* fil = fopen(fname.c_str(), "r");
  if (!fil) {
    return false;
  }

  HISTIN.push_back(fil);
  PauseRecording();
  bool ret = PlaybackFile();
  if (HISTIN.size() != 0 && HISTIN.back() == fil) {
    StopPlaying();
  }
  ResumeRecording();
  return ret;
}

bool MIHistory::AddCommand(MIData& data) {
  if (!IsRecording() || IsRecordingPaused()) {
    return false;
  }
  if (IsPlaying()) { // don't add new stuff if we're playing back a line
    return true;
  }

  std::string str = CommandToString(data);
  WriteLine(str);
  return true;
}

static void searchAndReplace(std::string& str, const std::string &search, const std::string &replace) {
  std::string::size_type i = str.find(search);
  while (i != std::string::npos) {
    str.erase(i, search.length());
    str.insert(i, replace);
    i = str.find(search);
  }
}

bool MIHistory::AddComment(const std::string& str) {
  if (!IsRecording() || IsRecordingPaused()) {
    return false;
  }
  if (IsPlaying()) { // don't add new stuff if we're playing back a line
    return true;
  }
  std::string comment = str;
  searchAndReplace(comment, "\n", "\\n");
  searchAndReplace(comment, "\r", "\\r");
  WriteLine("--" + comment);
  return true;
}

void MIHistory::WriteLine(const std::string& line) {
  if (line.length() > 0) {
    fprintf(HISTOUT, "%s\n", line.c_str());
    fflush(HISTOUT);
  }
}


bool MIHistory::IsRecording() const {
  return HISTOUT != NULL;
}

bool MIHistory::StopRecording() {
  if (!IsRecording()) {
    return false;
  }
  fclose(HISTOUT);
  HISTOUT = NULL;
  return true;
}

bool MIHistory::IsPlaying() const {
  return HISTIN.size() != 0;
}

bool MIHistory::StopPlaying() {
  if (!IsPlaying()) {
    return false;
  }

  FILE* f = HISTIN.back();
  HISTIN.pop_back();

  fclose(f);
  return true;
}




bool MIHistory::AddMenuEvent(int id) {
  if (!IsRecording() || IsRecordingPaused()) {
    return false;
  }
  std::string str, name;
  void* frame = 0;

  // try all menu bars, then popup menus in order of frames
  //     first) current frame
  //    second) default frame
  //   finally) all other frames
  for (int vi = -1; vi < (int)_frames.size() && !str.size() ; ++vi) {
    frame = ( vi < 0 ? _current_frame : _frames[vi]);
    std::map<std::string, MIMenuBar*>::iterator i;
    for (i = _menu_bars[frame].begin(); i != _menu_bars[frame].end() && !str.size(); ++i) {
      name = i->first;
      str = i->second->Stringify(id);
    }

    std::map<std::string, MIMenu*>::iterator j;
    for (j = _menus[frame].begin(); j != _menus[frame].end() && !str.size(); ++j) {
      name = j->first;
      str = j->second->Stringify(id);
    }
  }

  if (str.size() == 0) {
    return false;
  }

  MIData values;
  values["command"].str = "menu";
  values["menu"].str = format("MENU: %s | %s", name.c_str(), str.c_str());
  return MIGetHistory()->AddCommand(values);
}

void MIHistory::AddMenuBar(MIMenuBar* mb, const std::string& name) {
  // printf("Adding menu bar '%s', %p to frame %p\n",name.c_str(),mb,_current_frame);
  _menu_bars[_current_frame][name] = mb;
}

void MIHistory::RemoveMenuBar(MIMenuBar* mb) {
  for (int vi = 0; vi < (int)_frames.size(); ++vi) {
    std::map<std::string, MIMenuBar*>& m = _menu_bars[_frames[vi]];
    std::map<std::string, MIMenuBar*>::iterator i;
    for (i = m.begin(); i != m.end();) {
      if (i->second == mb) {
        m.erase(i++);
      } else {
        ++i;
      }
    }
  }
}


void MIHistory::AddMenu(MIMenu* menu, const std::string& name) {
  // printf("Adding menu '%s', %p to frame %p\n",name.c_str(),menu,_current_frame);
  _menus[_current_frame][name] = menu;
}

void MIHistory::RemoveMenu(MIMenu* mb) {
  for (int vi = 0; vi < (int)_frames.size(); ++vi) {
    std::map<std::string, MIMenu*>& m = _menus[_frames[vi]];
    std::map<std::string, MIMenu*>::iterator i;
    for (i = m.begin(); i != m.end();) {
      if (i->second == mb) {
        m.erase(i++);
      } else {
        ++i;
      }
    }
  }
}


bool MIHistory::PlaybackMenuLine(const std::string& line) {
  char buf[1024];
  if (sscanf(line.c_str(), "MENU: %s", buf) == 0) {
    return false;
  }

  std::string::size_type pos = line.find_first_of('|');
  if (pos == std::string::npos) {
    return false;
  }


  // try all menu bars, then popup menus in order of frames
  //     first) current frame
  //    second) default frame
  //   finally) all other frames
  for (int vi = -1; vi < (int)_frames.size(); ++vi) {
    void* frame = ( vi < 0 ? _current_frame : _frames[vi]);

    // search menu bar(s)
    std::map<std::string, MIMenuBar*>& mb = _menu_bars[frame];
    std::map<std::string, MIMenuBar*>::iterator i = mb.find(buf);
    if (i != mb.end()) {
      QAction *act=i->second->findAction(&line[pos+2]);
      if (act) {
        act->trigger();
        return true;
      } 
    }

    // search popup menu(s)
    std::map<std::string, MIMenu*>& m = _menus[frame];
    std::map<std::string, MIMenu*>::iterator j = m.find(buf);
    if (j != m.end()) {
      QAction *act=j->second->findAction(&line[pos+2]);
      if (act) {
        act->trigger();
        return true;
      }
    }
  }

  Logger::log("Error in menu line: %s\n",line.c_str());
  return false;
}

void MIHistory::PushDialogResponse(const std::string &line) {
  _dialog_responses.push_back(line);
}

void MIHistory::ClearDialogResponses() {
  _dialog_responses.clear();
}

bool MIHistory::GetDialogResponse(std::string& line) {
  line.clear();

  if (_dialog_responses.size()!=0) {
    line=_dialog_responses.front();
    _dialog_responses.pop_front();
    return true;
  }

  if (!IsPlaying() || IsPlaybackPaused()) {
    return false;
  }

  char buf[4096];
  long pos = ftell(HISTIN.back());
  while (IsPlaying() // command might have closed file via StopPlaying
         && !feof(HISTIN.back())) {
    if (fgets(buf, 4096, HISTIN.back())) {
      if (strlen(buf) == 0) {
        continue;
      }
      // zap trailing \n
      if (buf[strlen(buf)-1] == '\n') {
        buf[strlen(buf)-1] = 0;
      }
      if (buf[strlen(buf)-1] == '\r') {
        buf[strlen(buf)-1] = 0;
      }
      if (IsComment(buf)) {
        continue;
      }
      if (strncmp("dialog", buf, 6) == 0 && strlen(buf) >= 6) {
        if (strncmp("dialog{}", buf, 8) == 0) {
          return false;
        }
        line = std::string(&buf[6]);
        return true;
      }

      // if we're here we got a non-comment, non DLOG line which was unexpected,
      // restore file position and return false
      fseek(HISTIN.back(), pos, SEEK_SET);
      return false;
    }
  }
  return false;
}

void MIHistory::ProcessGUIEvents()
{
  QCoreApplication::processEvents(QEventLoop::ExcludeUserInputEvents);
}

bool MIHistory::PlaybackFile() {
  char buf[4096];

  //unsigned int cnt=0;
  while (IsPlaying() // command might have closed file via StopPlaying
         && !feof(HISTIN.back())) {
    if (fgets(buf, 4096, HISTIN.back())) {
      if (strlen(buf) == 0) {
        continue;
      }
      // zap trailing \n
      if (buf[strlen(buf)-1] == '\n') {
        buf[strlen(buf)-1] = 0;
      }
      // and/or trailing \r (for windows/linux history file compatibility)
      if (buf[strlen(buf)-1] == '\r') {
        buf[strlen(buf)-1] = 0;
      }
      // printf("HISTORY(%05d): %s\n",++cnt,buf);

      bool result = PlaybackCommand(buf);
      if (!result) {
        return result;
      }
    }
    ProcessGUIEvents();
  }
  return true;
}

bool MIHistory::PlaybackCommand(const char* line) {

  Logger::debug("Playback> %s", line);
  if (IsComment(line)) {
    return true;
  }
  bool result = false;

  MIGLWidget* view = MIMainWindow::instance()->currentMIGLWidget();

  if (strncmp("menu", line, 4) == 0) {
    MIData data;
    StringToCommand(&line[4], data);
    result = PlaybackMenuLine(data["menu"].str);
  } else if (strncmp("dialog", line, 6) == 0) {
    // ignore, should not get this here
    Logger::debug("Error reading history file: unexpected dialog response.");
  } else if (strncmp("pick", line, 4) == 0) {
    MIData data;
    StringToCommand(&line[4], data);
    data["command"].str = "pick";
    if (view)
      result = view->HandleHistory(data);
  } else if (strncmp("key", line, 3) == 0) {
    MIData data;
    StringToCommand(&line[3], data);
    data["command"].str = "key";
    if (view)
      result = view->HandleHistory(data);
  } else if (strncmp("mouse", line, 5) == 0) {
    MIData data;
    StringToCommand(&line[5], data);
    if (view)
      result = view->PlaybackMouseHistory(data);
  } else if (strncmp("tree", line, 4) == 0) {
    MIData data;
    StringToCommand(&line[4], data);
    ViewSyncedPanel* panel = MIMainWindow::instance()->GetModelsTree();
    if (panel != NULL) {
      result = panel->HandleHistory(data);
    }
  } else if (strncmp("dict", line, 4) == 0) {
    MIFitGeomRefiner()->EditEntryCleanup(false);
    MIData data;
    StringToCommand(&line[4], data);

    Q_FOREACH (QWidget *w, QApplication::topLevelWidgets()) {
      DictEditDialog *dlg=dynamic_cast<DictEditDialog*>(w);
      if (dlg) {
        if (data["state"].str=="ok") {
          dlg->accept();
          result = true;
          break;
        } else {
          dlg->reject();
          result = true;
          break;
        }
      }
    }
  } else if (strncmp("frame", line, 5) == 0) {
    MIData data;
    StringToCommand(&line[4], data);
    unsigned int vnum = data["index"].u;
    vnum++;  // account for default frame (0)

    if (vnum < _frames.size()) {
      _current_frame = _frames[vnum];
      MIMainWindow::instance()->setActiveMIGLWidget((MIGLWidget*)_frames[vnum]);
    }
    result = true;
  } else if (strncmp("seqw", line, 4) == 0) {
    MIData data;
    StringToCommand(&line[4], data);
    unsigned int anum, rnum, mnum;
    Molecule* mol = 0;
    RESIDUE* res = 0;
    MIAtom* atom = 0;
    char c;
    if (view) {
      anum = data["atom"].u;
      rnum = data["res"].u;
      mnum = data["mol"].u;
      c = data["ss"].str[0];
      if (PickDeSerialize(view->GetDisplaylist(),
                          atom, res, mol, anum, rnum, mnum) || !res) {

        if (data["type"].str == "select") {
          view->select(mol, res, atom);
        } else if (data["type"].str == "change") {
          res->setSecstr(c);
        } else if (data["type"].str == "focus") {
          view->setFocusResidue(res);
        }
        PaletteChanged = true;
        view->ReDraw();
        result = true;
      }
    }
  } else if (strncmp("rama", line, 4) == 0) {
    MIData data;
    StringToCommand(&line[4], data);
    unsigned int anum, rnum, mnum;
    Molecule* mol = 0;
    RESIDUE* res = 0;
    MIAtom* atom = 0;
    float increment;
    if (view) {
      if (data["type"].str == "PHI" || data["type"].str == "PSI") {
        increment = data["inc"].f;
        view->FitTorsion(data["type"].str.c_str());
        view->fitmol->RotateTorsion(increment);
        view->UpdateCurrent();
        view->ReDraw();
        result = true;
      } else if (data["type"].str == "focus") {
        anum = data["atom"].u;
        rnum = data["res"].u;
        mnum = data["mol"].u;
        if (PickDeSerialize(view->GetDisplaylist(),
                            atom, res, mol, anum, rnum, mnum) && res) {
          view->setFocusResidue(res, true);
          result = true;
        }
      }
    }
  } else if (strncmp("prune", line, 5) == 0) {
    LimitGrowth();
  } else if (strncmp("rand", line, 4) == 0) {
    MIData data;
    StringToCommand(&line[4], data);
    if (data["type"].str == "stack_range") {
      if (view)
        view->AtomStack->ExpandTop2AllAtoms();
      result = true;
    } else {
      unsigned int anum, rnum, mnum, recenter;
      Molecule* mol = 0;
      RESIDUE* res = 0;
      MIAtom* atom = 0;
      anum = data["atom"].u;
      rnum = data["res"].u;
      mnum = data["mol"].u;
      recenter = data["recenter"].b;
      if (view) {
        if (PickDeSerialize(view->GetDisplaylist(),
                            atom, res, mol, anum, rnum, mnum)) {
          view->AtomStack->Push(atom, res, mol);
          result = true;
        }
      }
    }
  }
  Logger::debug("Playback> %s", (result ? "success" : "failure"));
  return result;
}

bool MIHistory::PickSerialize(Displaylist* dl,
                              MIAtom* atom,
                              RESIDUE* res,
                              Molecule* mol,
                              unsigned int& anum,
                              unsigned int& rnum,
                              unsigned int& mnum) {
  if (!IsRecording() || IsRecordingPaused()) {
    return false;
  }
  mnum = 0;
  rnum = 0;
  anum = 0;

  for (std::list<Molecule*>::iterator i = dl->begin(); i != dl->end(); ++i) {
    if (!mol || *i == mol) {
      mol = *i; // in case it's null, make sure we define it
      rnum = 0;

      MIIter<RESIDUE> r = mol->GetResidues();
      if (!r) {
        return true;  // there are no residues, so return mnum, 0, 0

      }
      for (; r; ++r) {
        if (!res || r == res) {
          if (!atom || r->atomCount() == 0) { // no atoms, so return mnum, rnum, 0
            return true;
          }
          for (anum = 0; anum < (unsigned int)r->atomCount(); ++anum) {
            if (atom == r->atom(anum)) {
              return true;
            }
          }
        }
        rnum++;
      }
    }
    mnum++;
  }
  return false;
}

bool MIHistory::PickDeSerialize(Displaylist* dl,
                                MIAtom*& atom,
                                RESIDUE*& res,
                                Molecule*& mol,
                                unsigned int anum,
                                unsigned int rnum,
                                unsigned int mnum) {
  unsigned int mcnt = 0, rcnt = 0;
  mol = 0;
  for (std::list<Molecule*>::iterator i = dl->begin(); i != dl->end(); ++i) {
    if (mcnt == mnum) {
      mol = *i;
      break;
    }
    ++mcnt;
  }
  if (!mol) {
    return false;
  }

  res = 0;
  for (MIIter<RESIDUE> r = mol->GetResidues(); r; ++r) {
    if (rcnt == rnum) {
      res = r;
      break;
    }
    ++rcnt;
  }
  if (!res) {
    return false;
  }

  if (anum >= (unsigned int)res->atomCount()) {
    return false;
  }
  atom = res->atom(anum);
  return true;
}

bool MIHistory::AddBondPick(const std::string& prefix,
                            Displaylist* dl,
                            MIAtom* atom1,
                            RESIDUE* res1,
                            Molecule* mol1,
                            MIAtom* atom2,
                            RESIDUE* res2,
                            Molecule* mol2) {
  if (!IsRecording() || IsRecordingPaused()) {
    return false;
  }
  unsigned int anum1, rnum1, mnum1, anum2, rnum2, mnum2;
  if (PickSerialize(dl, atom1, res1, mol1, anum1, rnum1, mnum1)
      && PickSerialize(dl, atom2, res2, mol2, anum2, rnum2, mnum2)) {
    MIData values;
    values["command"].str = prefix;
    values["type"].str = "bond";
    values["atom1"].u = anum1;
    values["res1"].u = rnum1;
    values["mol1"].u = mnum1;
    values["atom2"].u = anum2;
    values["res2"].u = rnum2;
    values["mol2"].u = mnum2;
    return MIGetHistory()->AddCommand(values);
  }
  return false;
}

bool MIHistory::AddAtomPick(const std::string& prefix,
                            Displaylist* dl,
                            MIAtom* atom,
                            RESIDUE* res,
                            Molecule* mol,
                            bool recenter) {
  if (!IsRecording() || IsRecordingPaused()) {
    return false;
  }
  unsigned int anum = 0, rnum = 0, mnum = 0;
  if (!PickSerialize(dl, atom, res, mol, anum, rnum, mnum)) {
    return false;
  }

  MIData values;
  values["command"].str = prefix;
  values["type"].str = "atom";
  values["atom"].u = anum;
  values["res"].u = rnum;
  values["mol"].u = mnum;
  values["recenter"].b = recenter;
  return MIGetHistory()->AddCommand(values);
}

bool MIHistory::AddFrame(void* v) {
  // printf("Adding frame %s: %p\n",name.c_str(),v);
  _current_frame = v;
  _frames.push_back(v);
  return true;
}

bool MIHistory::RemoveFrame(void* v) {
  if (!v) {
    return false;
  }

  // remove menu bars and popup menus associated with frame (client code
  // can still call RemoveMenu, but it won't have anything to do if this
  // is called first)
  _menu_bars.erase(v);
  _menus.erase(v);

  // reset current frame, if necessary
  if (_current_frame == v) {
    _current_frame = 0;
  }

  /// remove the frame and frame name
  std::vector<void*>::iterator frames = _frames.begin();
  for (; frames != _frames.end(); ++frames) {
    if (*frames == v) {
      _frames.erase(frames);
      return true;
    }
  }

  return false;
}

bool MIHistory::ActivateFrame(void* v) {
  if (!IsRecording() || IsRecordingPaused()) {
    return false;
  }
  if (v==_current_frame)
    return true;
  _current_frame = v;
  for (unsigned int i = 0; i < _frames.size(); ++i) {
    if (v == _frames[i]) {
      // printf("Activating view: %s %p\n",_view_names[i].c_str(),v);
      MIData data;
      data["command"].str = "frame";
      data["index"].u = i;
      return AddCommand(data);
    }
  }
  return false;
}

MIHistory* MIGetHistory() {
  static MIHistory* h = NULL;
  if (h == NULL) {
    h = new MIHistory();    
  }
  return h;
}

bool MIHistory::IsComment(const char* line) {
  return (strlen(line) >= 2) && (line[0] == '-') && (line[1] == '-');
}

std::string MIHistory::CommandToString(MIData& data) {
  std::string command = data["command"].str;
  std::string str;
  data.erase("command");
  
  return command+MIDataToString(data);
}

bool MIHistory::StringToCommand(const std::string& str, MIData& data) {
  return StringToMIData(str, data);
}
