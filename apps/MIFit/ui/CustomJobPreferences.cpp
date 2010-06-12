#include "CustomJobPreferences.h"
#include "ui_CustomJobPreferences.h"
#include "MIMainWindow.h"
#include <QMenu>

CustomJobPreferences::CustomJobPreferences(QWidget *parent) :
    QWidget(parent),
    ui(new Ui::CustomJobPreferences)
{
    ui->setupUi(this);

    updateMenuList();
    updateRemoveButton();
}

CustomJobPreferences::~CustomJobPreferences()
{
    delete ui;
}

void CustomJobPreferences::savePreferences()
{
    QMenu *jobMenu = MIMainWindow::instance()->jobMenu();
    QList<QAction*> actions = jobMenu->actions();
    int listIndex = 0;
    int actionIndex = 0;

    // Skip to first separator
    while (actionIndex < actions.size())
    {
        if (actions.at(actionIndex)->isSeparator())
            break;
        ++actionIndex;
    }

    // Remove actions not found in list
    while (listIndex < ui->menuList->count()
        && actionIndex < actions.size())
    {
        QAction *action = actions.at(actionIndex);
        if (!action->isSeparator())
        {
            if (ui->menuList->item(listIndex)->text() != action->text())
                jobMenu->removeAction(action);
            else
                ++listIndex;
        }

        ++actionIndex;
    }
    while (actionIndex < actions.size())
    {
        QAction *action = actions.at(actionIndex);
        if (!action->isSeparator())
            jobMenu->removeAction(action);
        ++actionIndex;
    }

    MIMainWindow::instance()->saveJobMenu();
}

void CustomJobPreferences::updateMenuList()
{
    ui->menuList->clear();
    QMenu *jobMenu = MIMainWindow::instance()->jobMenu();
    QList<QAction*> actions = jobMenu->actions();
    QList<QAction*>::iterator actionIter = actions.begin();

    // Skip to first separator
    for (; actionIter != actions.end(); ++actionIter)
        if ((*actionIter)->isSeparator())
            break;

    for (; actionIter != actions.end(); ++actionIter)
    {
        QAction *action = *actionIter;
        if (action->isSeparator())
            continue;
        ui->menuList->addItem(action->text());
    }
}

void CustomJobPreferences::updateRemoveButton()
{
    ui->removeButton->setEnabled(ui->menuList->selectedItems().size());
}

void CustomJobPreferences::on_removeButton_clicked()
{
    QList<QListWidgetItem*> items = ui->menuList->selectedItems();
    qDeleteAll(items.begin(), items.end());
}
