#ifndef SADPHASING_H
#define SADPHASING_H

#include <string>

#include "core/corelib.h"

#include "ui_SadPhasing.h"
#include "MIQDialog.h"

class QAbstractButton;

class SadPhasing : public MIQDialog, public Ui::SadPhasing
{
    Q_OBJECT

public:
    SadPhasing(QWidget *parent = 0);
    static void GetInitialData(MIData &dat);

    void InitializeFromData(const MIData &dat);
    bool GetData(MIData &data);

  private slots:
    void validateTimeout();

  private:
    QAbstractButton *_okButton;
};

#endif
