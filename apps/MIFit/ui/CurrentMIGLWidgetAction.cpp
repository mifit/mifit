#include "CurrentMIGLWidgetAction.h"

#include "MIMainWindow.h"
#include "MIGLWidget.h"
#include <QDebug>
#include <QMenu>
#include <QMetaMethod>
#include <QMetaObject>

namespace
{
    int methodCode(const char *method)
    {
        return static_cast<int>(*method - '0');
    }

    void invokeSlot(QObject *object, const char *slot, QAction *arg = 0)
    {
        const QMetaObject *meta = object->metaObject();
        int index = meta->indexOfSlot(slot);
        if (index < 0)
            qWarning("CurrentMIGLWidgetAction: Unable to find slot %s on %s", slot, meta->className());
        else
        {
            QMetaMethod method = meta->method(index);
            bool invoked = false;
            if (arg)
                invoked = method.invoke(object, Qt::DirectConnection, Q_ARG(QAction*, arg));
            else
                invoked = method.invoke(object, Qt::DirectConnection);

            if (!invoked)
                qWarning("CurrentMIGLWidgetAction: Unable to invoke slot %s on %s",
                         method.signature(), meta->className());
        }
    }

}

CurrentMIGLWidgetAction::CurrentMIGLWidgetAction(const QString& text, const QString& statusTip, QMenu *menu,
                                                 const char *slot, const char *updateSlot)
    : QAction(text, menu), _slot(slot), _updateSlot(updateSlot)
{
    if (_slot)
    {
        if (methodCode(_slot) != QSLOT_CODE)
            qWarning("CurrentMIGLWidgetAction: Use the SLOT macro for %s", _slot);
        else
            ++_slot;
    }

    if (_updateSlot)
    {
        if (methodCode(_updateSlot) != QSLOT_CODE)
            qWarning("CurrentMIGLWidgetAction: Use the SLOT macro for %s", _updateSlot);
        else
            ++_updateSlot;
    }

    setStatusTip(statusTip);
    connect(this, SIGNAL(triggered()), SLOT(on_triggered()));

    menu->addAction(this);
    if (_updateSlot)
        connect(menu, SIGNAL(aboutToShow()), this, SLOT(update()));
}

void CurrentMIGLWidgetAction::on_triggered()
{
    if (_slot && MIMainWindow::instance()->currentMIGLWidget())
        invokeSlot(MIMainWindow::instance()->currentMIGLWidget(), _slot);
}

void CurrentMIGLWidgetAction::update()
{
    if (_updateSlot && MIMainWindow::instance()->currentMIGLWidget())
        invokeSlot(MIMainWindow::instance()->currentMIGLWidget(), _updateSlot, this);
}
