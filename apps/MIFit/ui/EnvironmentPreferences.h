#ifndef ENVIRONMENT_PREFERENCES_H
#define ENVIRONMENT_PREFERENCES_H

#include "core/corelib.h"

#include "ui_EnvironmentPreferences.h"

class EnvironmentPreferences : public QWidget, public Ui::EnvironmentPreferences
{
    Q_OBJECT

public:
    EnvironmentPreferences(QWidget *parent = 0);
    void savePreferences();

public slots:
    void on_crystalDataButton_pressed();
    void on_customDictionaryButton_pressed();
    void on_shelxHomeButton_pressed();
    void on_htmlBrowserButton_pressed();
    void on_checkpointDirectoryButton_pressed();
};

#endif // ifndef ENVIRONMENT_PREFERENCES_H
