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

    void setJobName(const QString& jobName);
    void setProgram(const QString& program);
    void setArguments(const QString& arguments);
    enum ModelMode {
      CURRENT,
      FILE
    };
    void setModelMode(CustomJobDialog::ModelMode mode);
    void setWorkingDirectory(const QString& dir);
    void setModelFile(const QString& modelFile);
    void setDataFile(const QString& dataFile);

    QString jobName() const;
    QString program() const;
    QString arguments() const;
    ModelMode modelMode() const;
    QString workingDirectory() const;
    QString modelFile() const;
    QString dataFile() const;

  private slots:
    void validateTimeout();

  private:
    QAbstractButton *_okButton;
};

#endif
