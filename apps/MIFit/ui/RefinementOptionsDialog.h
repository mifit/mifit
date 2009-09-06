#ifndef REFINEMENTOPTIONSDIALOG_H
#define REFINEMENTOPTIONSDIALOG_H

#include "core/corelib.h" // for MIData

#include "ui_RefinementOptionsDialog.h"

class RefinementOptionsDialog : public QDialog, public Ui::RefinementOptionsDialog
{
    Q_OBJECT

public:
    RefinementOptionsDialog(const MIData &dat, QWidget *parent = 0);
    void GetResults(MIData &data);
};

#endif
