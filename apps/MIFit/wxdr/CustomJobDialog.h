#ifndef CUSTOMJOBDIALOG_H
#define CUSTOMJOBDIALOG_H

#include <string>

#include "core/corelib.h"

#include "ui_CustomJobDialog.h"
#include "MIQDialog.h"

class QAbstractButton;

class CustomJobDialog : public MIQDialog, public Ui::CustomJobDialog
{
    Q_OBJECT

public:
    CustomJobDialog(QWidget *parent = 0);
    static void GetInitialData(MIData &dat);

    void InitializeFromData(const MIData &dat);
    bool GetData(MIData &data);

  private Q_SLOTS:
    void validateTimeout();

  private:
    QAbstractButton *_okButton;
};

#endif
