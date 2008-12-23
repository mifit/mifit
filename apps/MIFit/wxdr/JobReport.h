#ifndef JOBREPORT_H
#define JOBREPORT_H

#include <string>

#include "corelib.h"

#include "ui_JobReport.h"

#include "MIQDialog.h"

class QAbstractButton;

class JobReport : public MIQDialog, public Ui::JobReport
{
    Q_OBJECT

public:
    JobReport(QWidget *parent = 0);
    static void GetInitialData(MIData &dat);

    void InitializeFromData(const MIData &dat);
    bool GetData(MIData &data);

  private Q_SLOTS:
    void validateTimeout();

  private:
    QAbstractButton *_okButton;
};

#endif
