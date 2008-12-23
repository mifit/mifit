#ifndef MI_HISTORY_H
#define MI_HISTORY_H

#include "corelib.h"
#include <map>
#include <deque>
#include "core/MIData.h"

class MIMenuBar;
class MIMenu;
class Displaylist;

class MIHistory {
  MIHistory();

  friend MIHistory* MIGetHistory();

public:
  ~MIHistory();
  
  bool Record(const std::string& fname);
  bool Play(const std::string& fname);
  bool AddComment(const std::string& str);
  bool AddCommand(MIData& data);
  static bool StringToCommand(const std::string& str, MIData& data);
  static std::string CommandToString(MIData& data);

  bool IsRecording() const;
  bool IsRecordingPaused() const {
    return recording_paused > 0;
  }

  void PauseRecording() {
    ++recording_paused;
  }

  void ResumeRecording() {
    if (recording_paused) {
      --recording_paused;
    }
  }

  bool StopRecording();

  bool IsPlaying() const;
  bool IsPlaybackPaused() const {
    return playback_paused > 0;
  }

  void PausePlayback() {
    ++playback_paused;
  }

  void ResumePlayback() {
    if (playback_paused) {
      --playback_paused;
    }
  }

  bool StopPlaying();


  void AddMenuBar(MIMenuBar* mb, const std::string& name);
  void RemoveMenuBar(MIMenuBar* mb);

  void AddMenu(MIMenu* menu, const std::string& name);
  void RemoveMenu(MIMenu* menu);

  bool AddFrame(void* v);
  bool RemoveFrame(void* v);
  bool ActivateFrame(void* v);

  bool AddMenuEvent(int id);

  bool AddAtomPick(const std::string& prefix,
                   Displaylist* dl,
                   chemlib::MIAtom* atom,
                   chemlib::RESIDUE* res,
                   Molecule* mol,
                   bool recenter = true);

  bool AddBondPick(const std::string& prefix,
                   Displaylist* dl,
                   chemlib::MIAtom* atom1,
                   chemlib::RESIDUE* res1,
                   Molecule* mol1,
                   chemlib::MIAtom* atom2,
                   chemlib::RESIDUE* res2,
                   Molecule* mol2);

  bool PickSerialize(Displaylist* dl,
                     chemlib::MIAtom* atom,
                     chemlib::RESIDUE* res,
                     Molecule* mol,
                     unsigned int& anum,
                     unsigned int& rnum,
                     unsigned int& mnum);
  bool PickDeSerialize(Displaylist* dl,
                       chemlib::MIAtom*& atom,
                       chemlib::RESIDUE*& res,
                       Molecule*& mol,
                       unsigned int anum,
                       unsigned int rnum,
                       unsigned int mnum);

  bool GetDialogResponse(std::string& line);
  void PushDialogResponse(const std::string &line);
  void ClearDialogResponses();

  static void ProcessGUIEvents();

  bool PlaybackCommand(const char* line);

private:
  bool PlaybackFile();
  bool PlaybackMenuLine(const std::string& buf);
  bool IsComment(const char* line);
  void WriteLine(const std::string& line);

  FILE* HISTOUT;
  std::vector<FILE*> HISTIN;
  unsigned int recording_paused;
  unsigned int playback_paused;
  std::map<void*, std::map<std::string, MIMenuBar*> > _menu_bars;
  std::map<void*, std::map<std::string, MIMenu*> > _menus;
  std::vector<void*> _frames;
  void* _current_frame;
  std::deque<std::string> _dialog_responses;
};

MIHistory* MIGetHistory();

#endif
