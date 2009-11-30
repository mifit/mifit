#include "ViewSyncedPanel.h"

#include <QStackedLayout>

#include "MIGLWidget.h"


ViewSyncedPanel::ViewSyncedPanel(QWidget *parent)
    : QWidget(parent)
{

    stackedLayout = new QStackedLayout(this);
    stackedLayout->setContentsMargins(0, 0, 0, 0);
    setLayout(stackedLayout);

    createDefaultNullPanel();
    setViewPanel(NULL);

    connect(ViewController::instance(), SIGNAL(viewActivated(MIGLWidget*)),
            this, SLOT(viewActivated(MIGLWidget*)));
    connect(ViewController::instance(), SIGNAL(viewDeactivated(MIGLWidget*)),
            this, SLOT(viewDeactivated(MIGLWidget*)));
}

ViewSyncedPanel::~ViewSyncedPanel()
{
}

void ViewSyncedPanel::createDefaultNullPanel()
{
    QWidget *panel = new QWidget(this);
    viewToPanel[NULL] = panel;
    stackedLayout->addWidget(panel);
}

void ViewSyncedPanel::createViewPanel(MIGLWidget *view)
{
    QWidget *panel = createPanelForView(view, this);
    viewToPanel[view] = panel;
    stackedLayout->addWidget(panel);
}

void ViewSyncedPanel::viewActivated(MIGLWidget *view)
{
    setViewPanel(view);
}

void ViewSyncedPanel::viewDeactivated(MIGLWidget *view)
{
    if (panelForViewExists(view))
    {
        QWidget *panel = viewToPanel[view];
        if (stackedLayout->currentWidget() == panel)
        {
            setViewPanel(NULL);
        }
        viewToPanel.erase(view);
        stackedLayout->removeWidget(panel);
        destroyContentsForView(view, panel);
        delete panel;
    }
}

void ViewSyncedPanel::setViewPanel(MIGLWidget *view)
{
    if (!panelForViewExists(view))
    {
        createViewPanel(view);
    }
    stackedLayout->setCurrentWidget(viewToPanel[view]);
}

bool ViewSyncedPanel::panelForViewExists(MIGLWidget *view)
{
    return viewToPanel.find(view) != viewToPanel.end();
}

