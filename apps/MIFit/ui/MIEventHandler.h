#ifndef MI_EVENT_HANDLER
#define MI_EVENT_HANDLER

#include <QAction>
#include <string>

class QObject;
class MIMenu;

//  Helper classes to emulate a wxWidgets-like API for associating event ids with their handlers


// a menu (and eventually toolbar) action was triggered by the user
class MIActionEvent
{
public:
    MIActionEvent(QAction *action = 0, unsigned int event_id = 0, bool state = false) : _action(action), _event_id(event_id), _state(state)
    {
    }

    unsigned int GetId() const
    {
        return _event_id;
    }
    bool IsChecked() const
    {
        return _state;
    }

    // QAction *GetAction() const { return _action; } prob not necessary

private:
    QAction *_action;
    unsigned int _event_id;
    bool _state;
};


// a menu (and eventually toolbar) update is needed
class MIUpdateEvent
{
public:
    MIUpdateEvent(MIMenu *menu, unsigned int event_id) : _menu(menu), _event_id(event_id)
    {
    }

    unsigned int GetId() const
    {
        return _event_id;
    }

    //MIMenu* GetMenu() const { return _menu; } // prob not necessary

    void Enable(bool state) const;
    void Check(bool state) const;

    bool GetChecked() const;
    bool GetEnabled() const;


private:
    MIMenu *_menu;
    unsigned int _event_id;
};


// can use this class to register a different event receiver for events
// which should be directed to a child object
class MIChildEventHandlerFtor
{
public:
    virtual MIChildEventHandlerFtor *CreateCopy()
    {
        return 0;
    }
    virtual ~MIChildEventHandlerFtor()
    {
    }
    virtual QObject*operator()()
    {
        return 0;
    }
};

//
//  classes which want to process their own menu/update events should subclass this
//
class MIEventHandler
{
public:
    MIEventHandler(QObject *obj) : _obj(obj)
    {
    }
    virtual ~MIEventHandler();

    virtual void registerActionHandler(unsigned int id, const std::string &fn,
                                       MIChildEventHandlerFtor *child_ftor = 0);
    virtual void registerUpdateHandler(unsigned int id, const std::string &fn,
                                       MIChildEventHandlerFtor *child_ftor = 0);

    virtual void handleAction(const MIActionEvent &ae);
    virtual void doUpdate(const MIUpdateEvent &ue);

    // returns true if this event handler or its designated child can handle the given action id
    bool validateAction(unsigned int id);

    // returns true if this event handler or its designated child can handle the given update id
    bool validateUpdate(unsigned int id);

protected:
    typedef std::map<unsigned int, std::pair<std::string, MIChildEventHandlerFtor*> > EventMap;
    EventMap _actionEventMap;
    EventMap _updateEventMap;
    QObject *_obj;
};

#endif // MI_EVENT_HANDLER
