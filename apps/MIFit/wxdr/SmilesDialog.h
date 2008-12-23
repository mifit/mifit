#ifndef SMILESDIALOG_H
#define SMILESDIALOG_H

#include "corelib.h" // for MIData

#include "ui_SmilesDialog.h"

class SmilesDialog : public QDialog, public Ui::SmilesDialog
{
    Q_OBJECT

public:
    SmilesDialog(QWidget *parent = 0);
    void GetResults(MIData &data);

public Q_SLOTS:
    void on_browseButton_clicked();
};

#endif
