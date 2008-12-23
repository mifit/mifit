#ifndef INTEGRATEDIALOG_H
#define INTEGRATEDIALOG_H

#include <string>

#include "corelib.h"

#include "ui_IntegrateDialog.h"
#include "MIQDialog.h"

class QAbstractButton;

class IntegrateDialog : public MIQDialog, public Ui::IntegrateDialog
{
    Q_OBJECT

public:
    IntegrateDialog(QWidget *parent = 0);
    static void GetInitialData(MIData &dat);

    void InitializeFromData(const MIData &dat);
    bool GetData(MIData &data);

  private Q_SLOTS:
    void on_removePushButton_clicked();
    void on_addPushButton_clicked();
    void validateTimeout();

  private:
    QAbstractButton *_okButton;
};

#endif
