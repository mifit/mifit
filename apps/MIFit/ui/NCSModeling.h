#ifndef NCSMODELING_H
#define NCSMODELING_H

#include <string>

#include "core/corelib.h"

#include "ui_NCSModeling.h"
#include "MIQDialog.h"

class QAbstractButton;

class NCSModeling : public MIQDialog, public Ui::NCSModeling
{
    Q_OBJECT

public:
    NCSModeling(QWidget *parent = 0);
    static void GetInitialData(MIData &dat);

    void InitializeFromData(const MIData &dat);
    bool GetData(MIData &data);

  private Q_SLOTS:
    void validateTimeout();
    void togglePhaseOptions(bool state);

  private:
    QAbstractButton *_okButton;
};

#endif
