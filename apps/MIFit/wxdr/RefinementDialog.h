#ifndef REFINEMENT_DIALOG_H
#define REFINEMENT_DIALOG_H

#include <string>

#include "corelib.h"

#include "ui_RefinementDialog.h"

#include "MIQDialog.h"

class QAbstractButton;

class RefinementDialog : public MIQDialog, public Ui::RefinementDialog
{
    Q_OBJECT

public:
    RefinementDialog(QWidget *parent = 0);
    static void GetInitialData(MIData &dat);

    void InitializeFromData(const MIData &dat);
    bool GetData(MIData &data);

  private Q_SLOTS:
    void validateTimeout();

  private:
    QAbstractButton *_okButton;
};

#endif
