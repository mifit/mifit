#ifndef MIMenuBar_H
#define MIMenuBar_H

#include <QMenuBar>
#include <string>
#include <vector>

class QAction;
class MIMenu;

class MIMenuBar : public QMenuBar {
  public:
    friend class MIToolBar;

    MIMenuBar(QMenuBar *mb) : _mb(mb) {}
    void Append(MIMenu *menu, const std::string &title);

    void Enable(unsigned int id, bool state);

    // useful for toolbars, find the action for the given id
    QAction *findAction(unsigned int id);

    // find the action for the given menu string e.g. "Menu | Submenu | Item"
    QAction *findAction(const std::string &s);

    // return the menu string for the given item
    std::string Stringify(unsigned int id);

    bool doUpdate(unsigned int id);

    // returns true if all actions in all menus can be handled
    bool validateActions();

    // returns true if all updates in all menus can be handled
    bool validateUpdates();

    // only for hidden menu for toolbar actions which don't have parent menu item
    void AppendHidden(MIMenu *menu);

  private:
    QMenuBar *_mb;
    std::vector<MIMenu*> _menus;
};


#endif //MIMenuBar_H
