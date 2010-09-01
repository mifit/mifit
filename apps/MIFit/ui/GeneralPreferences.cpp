#include "GeneralPreferences.h"
#include "core/corelib.h"
#include "ui/uilib.h"

#include <QSettings>

GeneralPreferences::GeneralPreferences(QWidget *parent)
    : QWidget(parent)
{
    setupUi(this);
    Application *app = Application::instance();
    QSettings settings;

    bool autoOpen = settings.value("openResultsOnJobFinished", false).toBool();
    autoOpenCheckBox->setChecked(autoOpen);
    incrementallyColorCheckBox->setChecked(app->incrementallyColorModels);
    xfitMouseCheckBox->setChecked(app->xfitMouseMode);
    saveOnCloseCheckBox->setChecked(app->onCloseSaveActiveModelToPdb);
    concurrentJobsSpinBox->setValue(app->concurrentJobLimit);

    bool breakByDiscontinuityPref = settings.value("Options/breakByDiscontinuity", true).toBool();
    bool breakByNonpeptidePref = settings.value("Options/breakByNonpeptide", false).toBool();
    breakOnDiscontinuityCheckBox->setChecked(breakByDiscontinuityPref);
    breakOnNonPeptideCheckBox->setChecked(breakByNonpeptidePref);
}

void GeneralPreferences::savePreferences()
{
    Application *app = Application::instance();

    QSettings settings;
    settings.setValue("openResultsOnJobFinished", autoOpenCheckBox->isChecked());

    app->incrementallyColorModels = incrementallyColorCheckBox->isChecked();
    app->xfitMouseMode = xfitMouseCheckBox->isChecked();
    app->onCloseSaveActiveModelToPdb = saveOnCloseCheckBox->isChecked();
    app->concurrentJobLimit = concurrentJobsSpinBox->value();

    settings.setValue("Options/breakByDiscontinuity", breakOnDiscontinuityCheckBox->isChecked());
    settings.setValue("Options/breakByNonpeptide", breakOnNonPeptideCheckBox->isChecked());
}
