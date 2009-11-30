#include "GeneralPreferences.h"
#include "core/corelib.h"
#include "ui/uilib.h"


GeneralPreferences::GeneralPreferences(QWidget *parent)
    : QWidget(parent)
{
    setupUi(this);
    Application *app = Application::instance();
    MIConfig *config = MIConfig::Instance();

    bool autoOpen;
    config->Read("openResultsOnJobFinished", &autoOpen, 0);
    autoOpenCheckBox->setChecked(autoOpen);
    incrementallyColorCheckBox->setChecked(app->incrementallyColorModels);
    xfitMouseCheckBox->setChecked(app->xfitMouseMode);
    saveOnCloseCheckBox->setChecked(app->onCloseSaveActiveModelToPdb);
    concurrentJobsSpinBox->setValue(app->concurrentJobLimit);

    bool breakByDiscontinuityPref;
    bool breakByNonpeptidePref;
    config->Read("Options/breakByDiscontinuity", &breakByDiscontinuityPref, true);
    config->Read("Options/breakByNonpeptide", &breakByNonpeptidePref, false);
    breakOnDiscontinuityCheckBox->setChecked(breakByDiscontinuityPref);
    breakOnNonPeptideCheckBox->setChecked(breakByNonpeptidePref);
}

void GeneralPreferences::savePreferences()
{
    Application *app = Application::instance();

    MIConfig *config = MIConfig::Instance();
    config->Write("openResultsOnJobFinished", autoOpenCheckBox->isChecked());

    app->incrementallyColorModels = incrementallyColorCheckBox->isChecked();
    app->xfitMouseMode = xfitMouseCheckBox->isChecked();
    app->onCloseSaveActiveModelToPdb = saveOnCloseCheckBox->isChecked();
    app->concurrentJobLimit = concurrentJobsSpinBox->value();

    config->Write("Options/breakByDiscontinuity", breakOnDiscontinuityCheckBox->isChecked());
    config->Write("Options/breakByNonpeptide", breakOnNonPeptideCheckBox->isChecked());
}
