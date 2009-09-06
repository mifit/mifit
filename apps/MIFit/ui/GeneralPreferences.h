#ifndef GENERAL_PREFERENCES_H
#define GENERAL_PREFERENCES_H

#include "core/corelib.h"

#include "ui_GeneralPreferences.h"

class GeneralPreferences : public QWidget, public Ui::GeneralPreferences
{
    Q_OBJECT

public:
    GeneralPreferences(QWidget *parent = 0);
    void savePreferences();

};

#endif
