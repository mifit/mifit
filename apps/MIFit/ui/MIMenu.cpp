#include <util/utillib.h>

#include "MIMenu.h"
#include "MIEventHandler.h"
#include "MIMainWindow.h"

#include "uitest.h" // for IsInTestMode

#include <algorithm>

MIMenu::MIMenu(MIEventHandler &receiver, QWidget* parent)
  : QMenu("", parent), _receiver(receiver)
{
  connect(this, SIGNAL(aboutToShow()), this, SLOT(showHandler()));
}

MIMenu::~MIMenu() {}

QAction *MIMenu::Append(unsigned int event_id, const std::string &text, const std::string &tooltip, bool checkable)
{
  if (event_id == 0) {
    event_id = _idToAction.size();
  }
  //parse shortcut out of text and add as shortcut
  QString qtext(text.c_str());
  QStringList sl=qtext.split('\t');
  QKeySequence key;
  if (sl.count()==0)
    return 0;
  if (sl.count() > 1)
    key=QKeySequence(sl.at(1));

  // and now create the action
  QAction *action=addAction(sl.at(0), this, SLOT(triggered(bool)), key);
  _actionToId[action]=event_id;
  _idToAction[event_id]=action;
  action->setToolTip(tooltip.c_str());
  action->setCheckable(checkable);

  return action;
}

QAction *MIMenu::AppendCheckItem(unsigned int event_id, const std::string &menu_text, const std::string &tooltip) {
  return Append(event_id, menu_text, tooltip, true);
}

void MIMenu::Append(unsigned int, const std::string &menu_text, MIMenu* submenu) {
  submenu->setTitle(menu_text.c_str());
  addMenu(submenu);
  _submenus.push_back(submenu);
}


void MIMenu::AppendSeparator() {
  addSeparator();
}


void MIMenu::Check(unsigned int id, bool state) {
  if (_idToAction.find(id) == _idToAction.end()) {
    MIMainWindowDebug("MIMenu::Check: invalid id");
    return;
  }

  QAction *action=_idToAction[id];
  if (!action->isCheckable()) {
    MIMainWindowDebug("MIMenu::Check: id is not checkable");
    return;
  }
  action->setChecked(state);
}

bool MIMenu::isChecked(unsigned int id)
{
  if (_idToAction.find(id) == _idToAction.end()) {
    MIMainWindowDebug("MIMenu::idChecked: invalid id");
    return false;
  }

  QAction *action=_idToAction[id];
  if (!action->isCheckable()) {
    MIMainWindowDebug("MIMenu::isChecked: id is not checkable");
    return false;
  }
  return action->isChecked();
}



void MIMenu::Enable(unsigned int id, bool state) {
  if (_idToAction.find(id) == _idToAction.end()) {
    MIMainWindowDebug("MIMenu::Enable: invalid id");
    return;
  }

  QAction *action=_idToAction[id];
  action->setEnabled(state);
}

bool MIMenu::isEnabled(unsigned int id)
{
  if (!QMenu::isEnabled()) // entire menu disabled
    return false;

  if (_idToAction.find(id) == _idToAction.end()) {
    MIMainWindowDebug("MIMenu::isEnabled: invalid id");
    return false;
  }

  QAction *action=_idToAction[id];
  return action->isEnabled();
}

void MIMenu::ClearMenu() {
  clear();
  _submenus.clear();
  _actionToId.clear();
  _idToAction.clear();
}

void MIMenu::triggered(bool state)
{
  QAction *action=dynamic_cast<QAction*>(sender());
  if (!action) {
    MIMainWindowDebug("MIMenu::triggered: invalid sender.\n");
    return;
  }

  if (_actionToId.find(action) == _actionToId.end()) {
    MIMainWindowDebug("MIMenu::triggered invalid (unknown) action!\n");
    return;
  }

  unsigned int id=_actionToId[action];

  // send enhanced signal to registered receiver or main window
  MIActionEvent ae(action,id,state);
  _receiver.handleAction(ae);
}


void MIMenu::showHandler(){

  std::map<QAction*, unsigned int>::iterator i=_actionToId.begin();

  for (; i != _actionToId.end(); ++i) {
    unsigned int id=i->second;

    MIUpdateEvent evt(this, id);
    _receiver.doUpdate(evt);
  }
}


QAction *MIMenu::findAction(unsigned int id) {

  if (_idToAction.find(id) != _idToAction.end())
    return _idToAction[id];

  for (size_t i=0; i < _submenus.size(); ++i) {
    QAction *act=_submenus[i]->findAction(id);
    if (act)
      return act;
  }

  return 0;
}

/*
  internal: guesses a descriptive text from a text suited for a menu entry
 */
static QString qt_strippedText(QString s)
{
    s.remove( QString::fromLatin1("...") );
    int i = 0;
    while (i < s.size()) {
        ++i;
        if (s.at(i-1) != QLatin1Char('&'))
            continue;
        if (i < s.size() && s.at(i) == QLatin1Char('&'))
            ++i;
        s.remove(i-1,1);
    }
    return s.trimmed();
};

std::string MIMenu::Stringify(unsigned int id) {

  if (_idToAction.find(id) != _idToAction.end()) {
    // it's in this menu
    QAction *act=_idToAction[id];
    std::string this_title=qt_strippedText(title()).toStdString();
    std::string ret=qt_strippedText(act->text()).toStdString();
    return this_title + " | " + ret;
  }


  for (size_t i=0; i < _submenus.size(); ++i) {
    std::string substr=_submenus[i]->Stringify(id);
    if (substr.size()) {
      std::string this_title=qt_strippedText(title()).toStdString();
      std::string res=this_title + " | " + substr;
      return res;
    }
  }

  return "";
}


QAction *MIMenu::findAction(const std::string &line) {
  std::string::size_type pos = line.find_first_of('|');

  if (pos == std::string::npos) {
    // no '|' found, just search *menu* for item and return...it's not a submenu
    for (std::map<QAction*, unsigned int>::iterator i=_actionToId.begin();
         i!=_actionToId.end(); ++i) {
      QAction *act=i->first;
      QString t=qt_strippedText(act->text());
      std::string s=t.toStdString();
      if (s == line) {
        return act;
      }
    }
    return 0;
  }

  if (pos < 2) {
    return 0;
  }

  // search the submenus
  std::string label = line.substr(0, pos-1);
  for (size_t i = 0; i < _submenus.size(); ++i) {
    QString t=_submenus[i]->title();
    if (t.toStdString() == label) {
      return _submenus[i]->findAction(&line[pos+2]);
    }
  }
  return 0;
}

bool MIMenu::doUpdate(unsigned int id) {
  if (!QMenu::isEnabled())
    return true; // update was successful, entire menu is disabled

  if (_idToAction.find(id) != _idToAction.end())
  {
    MIUpdateEvent evt(this, id);
    _receiver.doUpdate(evt);
    return true;
  }

  for (size_t i=0; i < _submenus.size(); ++i) {
    bool done=_submenus[i]->doUpdate(id);
    if (done)
      return true;
  }

  return false;
}

bool MIMenu::validateAction(unsigned int id) {
  return _receiver.validateAction(id);
}

bool MIMenu::validateActions() {
  bool ret=true;
  for (size_t i=0; i < _submenus.size(); ++i) {
    if (!_submenus[i]->validateActions())
      ret=false;
  }

  QList<QAction *> acts = actions();
  int actionCount = acts.count();
  for (std::map<unsigned int, QAction*>::iterator i=_idToAction.begin();
       i!=_idToAction.end(); ++i) {
    unsigned int id=i->first;
    if (!validateAction(id))
    {
      QString t=title();
      std::string menu_name=t.toStdString();
      std::string action_name=i->second->text().toStdString();
      MIMainWindow::instance()->Debug(::format("Menu %s -> %s: unhandled action, id: %d", menu_name.c_str(), action_name.c_str(), id));
      ret=false;
    }
  }
  return ret;
}

bool MIMenu::validateUpdates() {
  bool ret=true;
  for (size_t i=0; i < _submenus.size(); ++i) {
    if (!_submenus[i]->validateUpdates())
      ret=false;
  }

  for (std::map<unsigned int, QAction*>::iterator i=_idToAction.begin();
       i!=_idToAction.end(); ++i) {
    unsigned int id=i->first;
    if (!validateUpdate(id))
    {
      std::string menu_name=title().toStdString();
      std::string action_name=i->second->text().toStdString();
      MIMainWindow::instance()->Debug(::format("Menu %s -> %s: unhandled update, id: %d", menu_name.c_str(), action_name.c_str(), id));
      ret=false;
    }
  }
  return ret;
}


bool MIMenu::validateUpdate(unsigned int id) {
  return _receiver.validateUpdate(id);
}


bool MIMenu::DoExhaustiveTest(const std::vector<unsigned int> &exclusion_list, bool reverse) {
  bool ret=true;


  // do all submenus
  if (!reverse) {
    for (unsigned int i=0; i < _submenus.size(); ++i) {
      if (!_submenus[i]->DoExhaustiveTest(exclusion_list))
        ret=false;
    }
    // do all items in this menu

    // critical: call all update code for this menu to make sure enabled state is correct
    showHandler();

    for (std::map<QAction*, unsigned int>::iterator i=_actionToId.begin();
         i!=_actionToId.end(); ++i) {

      if (std::find(exclusion_list.begin(), exclusion_list.end(), i->second) == exclusion_list.end()) {
        // i.e. not in exclusion list
        QAction *act=i->first;
        if (act->isEnabled()) {
          MIMainWindowRightFooter(Stringify(i->second));
          act->trigger();
          showHandler(); // critical: call all update code for this menu to make sure enabled state is correct
        }
      }
    }


  } else {
    // reversed
    for (unsigned int i=_submenus.size(); i > 0; --i) {
      if (!_submenus[i-1]->DoExhaustiveTest(exclusion_list))
        ret=false;
    }

    // critical: call all update code for this menu to make sure enabled state is correct
    showHandler();

    // do all items in this menu
    for (std::map<QAction*, unsigned int>::reverse_iterator i=_actionToId.rbegin();
         i!=_actionToId.rend(); ++i) {

      if (std::find(exclusion_list.begin(), exclusion_list.end(), i->second) == exclusion_list.end()) {
        // i.e. not in exclusion list
        QAction *act=i->first;
        if (act->isEnabled()) {
          MIMainWindowRightFooter(Stringify(i->second));
          act->trigger();
          showHandler(); // critical: call all update code for this menu to make sure enabled state is correct
        }
      }
    }
  }



  return ret;
}

#ifdef _WIN32
#define random rand
#endif
bool MIMenu::DoRandomItem(const std::vector<unsigned int> &exclusion_list) {

  unsigned long i=random()%2;
  if (_submenus.size() && (i==1 || _actionToId.size() == 0)) {
    //pick random submenu;
    i=random()%_submenus.size();
    return _submenus[i]->DoRandomItem(exclusion_list);
  }

  if (!_actionToId.size())
    return true; // nothing to do, not an error.

  // critical: call all update code for this menu to make sure enabled state is correct
  showHandler();

  // try 5 times to get a valid (non-excluded, enabled) entry in this menu
  for (unsigned int trial=0; trial<5; ++trial) {

    i=random()%_actionToId.size();
    std::map<QAction*, unsigned int>::iterator j=_actionToId.begin();
    for (unsigned tmp=0; tmp < i && j!=_actionToId.end(); ++tmp, ++j) {
      ;
    }

    if (std::find(exclusion_list.begin(), exclusion_list.end(), j->second) == exclusion_list.end()) {
      // i.e. not in exclusion list
      QAction *act=j->first;
      if (act->isEnabled()) {
        MIMainWindowRightFooter(Stringify(j->second));
        act->trigger();
        return true;
      }
    }
  }

  return false;
}

QAction *MIMenu::doExec(const QPoint &p) {
  if (IsInTestMode()) {
    std::vector<unsigned int> excluded; // nothing excluded from popup menus
    DoRandomItem(excluded);
    return 0;
  }

  return exec(p);
}
