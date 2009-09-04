#include <nongui/nongui.h>
#include <util/utillib.h>
#include "MIMenuBar.h"
#include "MIMenu.h"

void MIMenuBar::Append(MIMenu *menu, const std::string &title)
{
  menu->setTitle(title.c_str());
  //menu->setParent(_mb);
  _mb->addMenu(menu);
  _menus.push_back(menu);
}

void MIMenuBar::AppendHidden(MIMenu *menu)
{
  _menus.push_back(menu);
}

void MIMenuBar::Enable(unsigned int id, bool state) {
  QAction *act=findAction(id);
  if (!act) {
    Logger::debug(::format("Can't find menu item %d for MIMenuBar::Enable",id));
    return;
  }
  act->setEnabled(state);
}


QAction* MIMenuBar::findAction(unsigned int id) {
  for (size_t i=0; i < _menus.size(); ++i) {
    QAction *act=_menus[i]->findAction(id);
    if (act)
      return act;
  }
  return 0;
}

bool MIMenuBar::doUpdate(unsigned int id) {
  for (size_t i=0; i < _menus.size(); ++i) {
    bool done=_menus[i]->doUpdate(id);
    if (done) 
      return true;
  }
  return false;
}

bool MIMenuBar::validateActions() {
  bool ret=true;
  for (size_t i=0; i < _menus.size(); ++i) {
    if (!_menus[i]->validateActions())
      ret=false;
  }
  return ret;
}

bool MIMenuBar::validateUpdates() {
  bool ret=true;
  for (size_t i=0; i < _menus.size(); ++i) {
    if (!_menus[i]->validateUpdates())
      ret=false;
  }
  return ret;
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

std::string MIMenuBar::Stringify(unsigned int id) {
  for (size_t i=0; i < _menus.size(); ++i) {
    std::string res=_menus[i]->Stringify(id);
    if (res.size())  {
      return res;
    }
  }
  return "";
}

#include <iostream>

QAction *MIMenuBar::findAction(const std::string &line) {
  std::string::size_type pos = line.find_first_of('|');
  if (pos < 2) {
    return 0;
  }

  std::string label = line.substr(0, pos-1);
  for (unsigned int i=0; i < _menus.size(); ++i) {
    QString t=qt_strippedText(_menus[i]->title());
    std::string s=t.toStdString();
    if (s == label)
      return _menus[i]->findAction(&line[line.find_first_of('|')+2]);
  }
  return 0;
}

bool MIMenuBar::DoExhaustiveTest(const std::vector<unsigned int> &exclusion_list, bool reverse) {
  bool ret=true;
  if (!reverse) {
    for (unsigned int i=0; i < _menus.size(); ++i) {
      if (!_menus[i]->DoExhaustiveTest(exclusion_list))
        ret=false;
    }
  } else {
    for (unsigned int i=_menus.size(); i > 0; --i) {
      if (!_menus[i-1]->DoExhaustiveTest(exclusion_list))
        ret=false;
    }
  }
  return ret;
}


#ifdef _WIN32
#define random rand
#endif
bool MIMenuBar::DoRandomItem(const std::vector<unsigned int> &exclusion_list) {
  if (_menus.size()==0)
    return false;

  long i=random()%_menus.size();
  return _menus[i]->DoRandomItem(exclusion_list);
}

