#include "PreferencesDialog.h"
#include "GeneralPreferences.h"
#include "EnvironmentPreferences.h"
#include "ContourOptionsWidget.h"
#include "CustomJobPreferences.h"

PreferencesDialog::PreferencesDialog(QWidget *parent)
    : QDialog(parent)
{
    setupUi(this);

    for (int i = 0; i < stackedWidget->count(); ++i)
    {
        QWidget *w = stackedWidget->widget(i);
        stackedWidget->removeWidget(w);
        delete w;
    }

    generalPrefs = new GeneralPreferences(this);
    environmentPrefs = new EnvironmentPreferences(this);
    contourOptions = new ContourOptionsWidget(this, true);
    customJobPrefs = new CustomJobPreferences(this);

    stackedWidget->addWidget(generalPrefs);
    stackedWidget->addWidget(environmentPrefs);
    stackedWidget->addWidget(contourOptions);
    stackedWidget->addWidget(customJobPrefs);
}

void PreferencesDialog::setPage(int index)
{
    listWidget->setCurrentRow(index, QItemSelectionModel::SelectCurrent);
}

void PreferencesDialog::savePreferences()
{
    generalPrefs->savePreferences();
    environmentPrefs->savePreferences();
    contourOptions->savePreferences();
    customJobPrefs->savePreferences();
}
