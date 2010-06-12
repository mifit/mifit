#ifndef PREFERENCES_DIALOG_H
#define PREFERENCES_DIALOG_H

#include "ui_PreferencesDialog.h"

class GeneralPreferences;
class EnvironmentPreferences;
class ContourOptions;
class CustomJobPreferences;

class PreferencesDialog : public QDialog, public Ui::PreferencesDialog
{
    Q_OBJECT
public:
    PreferencesDialog(QWidget *parent = 0);
    void setPage(int index);
    void savePreferences();

    static const int CustomJobsPage = 3;

private:
    GeneralPreferences *generalPrefs;
    EnvironmentPreferences *environmentPrefs;
    ContourOptions *contourOptions;
    CustomJobPreferences *customJobPrefs;
};

#endif // ifndef PREFERENCES_DIALOG_H
