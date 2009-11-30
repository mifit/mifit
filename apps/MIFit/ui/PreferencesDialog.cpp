#include "PreferencesDialog.h"
#include "GeneralPreferences.h"
#include "EnvironmentPreferences.h"
#include "ContourOptions.h"

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
    contourOptions = new ContourOptions(this, true);

    stackedWidget->addWidget(generalPrefs);
    stackedWidget->addWidget(environmentPrefs);
    stackedWidget->addWidget(contourOptions);
}


void PreferencesDialog::savePreferences()
{
    generalPrefs->savePreferences();
    environmentPrefs->savePreferences();
    contourOptions->savePreferences();
}
