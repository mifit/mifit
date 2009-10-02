#include <util/utillib.h>

#include "MIMenu.h"
#include "MIEventHandler.h"
#include "MIMainWindow.h"

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


QAction *MIMenu::doExec(const QPoint &p) {
  return exec(p);
}
