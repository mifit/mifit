#ifndef MIMENU_H
#define MIMENU_H

#include <string>
#include <map>
#include <vector>
#include <QMenu>

#include "MIEventHandler.h"

class MIMenu : public QMenu
{
    Q_OBJECT

public:
    // NOTE: if receiver is null, then events are sent to MIMainWindow::instance()
    MIMenu(MIEventHandler &receiver, QWidget *parent = 0);
    ~MIMenu();

    // append an item to the menu
    QAction *Append(unsigned int event_id, const std::string &menu_text,
                    const std::string &tooltip = "", bool checkable = false);

    // Note: same as append with checkable==true
    QAction *AppendCheckItem(unsigned int event_id, const std::string &menu_text,
                             const std::string &tooltip = "");

    // Note: event_id unused
    void Append(unsigned int /* event_id */, const std::string &menu_text,
                MIMenu *submenu);

    void AppendSeparator();

    void Check(unsigned int id, bool state);
    void Enable(unsigned int id, bool state);

    bool isChecked(unsigned int id);
    bool isEnabled(unsigned int id);

    MIEventHandler&GetReceiver()
    {
        return _receiver;
    }
    void ClearMenu();

    // useful for toolbars, find the action for the given id
    QAction *findAction(unsigned int id);

    // find the action for the given menu string e.g. "Menu | Submenu | Item"
    QAction *findAction(const std::string &s);

    // return the menu string for the given item
    std::string Stringify(unsigned int id);

    // call the receiver's update function
    bool doUpdate(unsigned int id);

    // returns true if the every menu item has a handler registered
    bool validateActions();

    // returns true if the every menu item has a handler registered
    bool validateUpdates();

    // do QMenu::exec(), respecting test mode
    QAction *doExec(const QPoint &p);

private slots:
    void triggered(bool);
    void showHandler();

private:
    // returns true if the given id has a handler registered
    bool validateAction(unsigned int id);

    // returns true if the given id has a handler registered
    bool validateUpdate(unsigned int id);

    std::map<QAction*, unsigned int> _actionToId;
    std::map<unsigned int, QAction*> _idToAction;
    MIEventHandler &_receiver;
    std::vector<MIMenu*> _submenus;
};
#endif // ifndef MIMENU_H
