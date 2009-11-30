#ifndef PREFERENCES_DIALOG_H
#define PREFERENCES_DIALOG_H

#include "ui_PreferencesDialog.h"

class GeneralPreferences;
class EnvironmentPreferences;
class ContourOptions;

class PreferencesDialog : public QDialog, public Ui::PreferencesDialog
{
    Q_OBJECT
public:
    PreferencesDialog(QWidget *parent = 0);
    void savePreferences();
private:
    GeneralPreferences *generalPrefs;
    EnvironmentPreferences *environmentPrefs;
    ContourOptions *contourOptions;
};

#endif // ifndef PREFERENCES_DIALOG_H
