#include "MIEventHandler.h"
#include "MIMenu.h"
#include "MIMainWindow.h"

#include <util/utillib.h>

#include <utility>
#include <QMetaObject>


void MIUpdateEvent::Check(bool state) const
{
    if (_menu)
        _menu->Check(_event_id, state);
}

void MIUpdateEvent::Enable(bool state) const
{
    if (_menu)
        _menu->Enable(_event_id, state);
}

bool MIUpdateEvent::GetChecked() const
{
    if (_menu)
        return _menu->isChecked(_event_id);
    return false;
}

bool MIUpdateEvent::GetEnabled() const
{
    if (_menu)
        return _menu->isEnabled(_event_id);
    return false;
}


MIEventHandler::~MIEventHandler()
{
    // delete child handler ftors
    EventMap::iterator i = _actionEventMap.begin();
    for (; i!=_actionEventMap.end(); ++i)
    {
        delete i->second.second;
    }

    i = _updateEventMap.begin();
    for (; i!=_updateEventMap.end(); ++i)
    {
        delete i->second.second;
    }
}


void MIEventHandler::registerActionHandler(unsigned int id, const std::string &fn,
                                           MIChildEventHandlerFtor *child_ftor)
{
    std::string f = fn;

    std::string::size_type i = f.find_last_of(":");
    if (i != std::string::npos)
    {
        f = std::string(&fn[i+1]);
    }

    MIChildEventHandlerFtor *ftor = child_ftor;
    if (ftor)
        ftor = ftor->CreateCopy();

    _actionEventMap[id] = std::pair<std::string, MIChildEventHandlerFtor*>(f, ftor);
}

bool MIEventHandler::validateAction(unsigned int id)
{
    EventMap::iterator i = _actionEventMap.find(id);
    if (i==_actionEventMap.end())
        return false;

    std::string slotname = i->second.first + "()";
    slotname = QMetaObject::normalizedSignature(slotname.c_str()).constData();
    QObject *child = 0;
    if (i->second.second)
        child = (*i->second.second)();

    if (child && child->metaObject()->indexOfSlot(slotname.c_str()) != -1)
        return true;

    return (_obj->metaObject()->indexOfSlot(slotname.c_str()) != -1);
}

bool MIEventHandler::validateUpdate(unsigned int id)
{
    EventMap::iterator i = _updateEventMap.find(id);
    if (i==_actionEventMap.end())
        return false;

    std::string slotname = i->second.first + "(const MIUpdateEvent&)";
    slotname = QMetaObject::normalizedSignature(slotname.c_str()).constData();
    QObject *child = 0;
    if (i->second.second)
        child = (*i->second.second)();

    if (child && child->metaObject()->indexOfSlot(slotname.c_str()) != -1)
        return true;

    return (_obj->metaObject()->indexOfSlot(slotname.c_str()) != -1);
}


void MIEventHandler::handleAction(const MIActionEvent &evt)
{
    EventMap::iterator i = _actionEventMap.find(evt.GetId());
    // note, if there's simply no handler registered for this id, we return sliently, since validate will report that
    if (i==_actionEventMap.end())
        return;

    const char *slotname = i->second.first.c_str();
    QObject *child = 0;
    if (i->second.second)
        child = (*i->second.second)();

    // first, try to invoke on child
    if (child && (QMetaObject::invokeMethod(child, slotname, Qt::DirectConnection, Q_ARG(MIActionEvent, evt))
                  || QMetaObject::invokeMethod(child, slotname, Qt::DirectConnection)))
    {
        MIMainWindow::instance()->UpdateToolBar();
        return;
    }

    // try on _obj
    if (!QMetaObject::invokeMethod(_obj, slotname, Qt::DirectConnection, Q_ARG(MIActionEvent, evt)))
    {
        MIMainWindow::instance()->Debug(format("Failed to invoke %s::%s(const MIActionEvent&) event id %d instance %p",
                                               _obj->metaObject()->className(), slotname, evt.GetId(), _obj));
        if (child)
            MIMainWindow::instance()->Debug(format("Failed to invoke %s::%s(const MIActionEvent&) event id %d instance %p",
                                                   child->metaObject()->className(), slotname, evt.GetId(), child));
    }
    if (!QMetaObject::invokeMethod(_obj, slotname, Qt::DirectConnection))
    {
        MIMainWindow::instance()->Debug(format("Failed to invoke %s::%s() event id %d instance %p",
                                               _obj->metaObject()->className(), slotname, evt.GetId(), _obj));
        if (child)
            MIMainWindow::instance()->Debug(format("Failed to invoke %s::%s() event id %d instance %p",
                                                   child->metaObject()->className(), slotname, evt.GetId(), child));
    }


    MIMainWindow::instance()->UpdateToolBar();
}


void MIEventHandler::registerUpdateHandler(unsigned int id, const std::string &fn,
                                           MIChildEventHandlerFtor *child_ftor)

{
    std::string f = fn;
    std::string::size_type i = f.find_last_of(":");
    if (i != std::string::npos)
    {
        f = std::string(&fn[i+1]);
    }

    MIChildEventHandlerFtor *ftor = child_ftor;
    if (ftor)
        ftor = ftor->CreateCopy();

    _updateEventMap[id] = std::pair<std::string, MIChildEventHandlerFtor*>(f, ftor);
}

void MIEventHandler::doUpdate(const MIUpdateEvent &evt)
{
    EventMap::iterator i = _updateEventMap.find(evt.GetId());
    // note, if there's simply no handler registered for this id, we return sliently, since validate will report that
    if (i==_updateEventMap.end())
        return;

    QObject *child = 0;
    if (i->second.second)
        child = (*i->second.second)();

    const char *slotname = i->second.first.c_str();

    // first, try to invoke on child
    if (child && (QMetaObject::invokeMethod(child, slotname, Qt::DirectConnection, Q_ARG(MIUpdateEvent, evt))
                  || QMetaObject::invokeMethod(child, slotname, Qt::DirectConnection)))
        return;

    // try on _obj
    if (!QMetaObject::invokeMethod(_obj, slotname, Qt::DirectConnection, Q_ARG(MIUpdateEvent, evt)))
    {
        MIMainWindow::instance()->Debug(format("Failed to invoke %s::%s(const MIUpdateEvent&) event id %d instance %p, disabling item",
                                               _obj->metaObject()->className(), slotname, evt.GetId(), _obj));
        if (child)
            MIMainWindow::instance()->Debug(format("Failed to invoke %s::%s(const MIUpdateEvent&) event id %d instance %p, disabling item",
                                                   child->metaObject()->className(), slotname, evt.GetId(), child));

        if (!QMetaObject::invokeMethod(_obj, slotname, Qt::DirectConnection))
        {
            MIMainWindow::instance()->Debug(format("Failed to invoke %s::%s() event id %d instance %p, disabling item",
                                                   _obj->metaObject()->className(), slotname, evt.GetId(), _obj));
            if (child)
                MIMainWindow::instance()->Debug(format("Failed to invoke %s::%s() event id %d instance %p, disabling item",
                                                       child->metaObject()->className(), slotname, evt.GetId(), child));

            evt.Enable(false);
        }
    }
}


