#ifndef BINDNGRIND_DIALOG_H
#define BINDNGRIND_DIALOG_H

#include <string>

#include "core/corelib.h"

#include "ui_BindNGrind.h"

#include "MIQDialog.h"

class QAbstractButton;

class BindNGrind : public MIQDialog, public Ui::BindNGrind
{
    Q_OBJECT

public:
    BindNGrind(QWidget *parent = 0);
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
