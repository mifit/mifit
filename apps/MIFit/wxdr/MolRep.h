#ifndef MOLREP_H
#define MOLREP_H

#include <string>

#include "core/corelib.h"

#include "ui_MolRep.h"
#include "MIQDialog.h"

class QAbstractButton;

class MolRep : public MIQDialog, public Ui::MolRep
{
    Q_OBJECT

public:
    MolRep(QWidget *parent = 0);
    static void GetInitialData(MIData &dat);

    void InitializeFromData(const MIData &dat);
    bool GetData(MIData &data);

  private Q_SLOTS:
    void validateTimeout();

  private:
    QAbstractButton *_okButton;
};

#endif
