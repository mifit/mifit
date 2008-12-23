#ifndef COCRYSTALSUPERPOS_H
#define COCRYSTALSUPERPOS_H

#include <string>

#include "corelib.h"

#include "ui_CocrystalSuperPos.h"

#include "MIQDialog.h"

class QAbstractButton;

class CocrystalSuperPos : public MIQDialog, public Ui::CocrystalSuperPos
{
    Q_OBJECT

public:
    CocrystalSuperPos(QWidget *parent = 0);
    static void GetInitialData(MIData &dat);

    void InitializeFromData(const MIData &dat);
    bool GetData(MIData &data);

  private Q_SLOTS:
    void validateTimeout();

  private:
    QAbstractButton *_okButton;
};

#endif
