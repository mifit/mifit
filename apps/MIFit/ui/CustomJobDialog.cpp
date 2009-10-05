#include "CustomJobDialog.h"
#include "MIBrowsePair.h"

#include "core/corelib.h"
#include <map/maplib.h>

#include <vector>
#include <string>

#include <QDialogButtonBox>
#include <QAbstractButton>
#include <QList>
#include <QTimer>
#include <QFileInfo>
#include <QFileDialog>


CustomJobDialog::CustomJobDialog(QWidget *parent) : MIQDialog(parent) {
  setupUi(this);

  new MIBrowsePair(commandPushButton, commandLineEdit,"Command file (*.* *)");
  new MIBrowsePair(workingDirPushButton, workingDirLineEdit,"", true);
  new MIBrowsePair(modelPushButton, modelLineEdit,"PDB file (*.pdb *.ent)");
  new MIBrowsePair(dataPushButton, dataLineEdit,"PDB file (*.pdb *.ent)");
  
  _okButton = buttonBox->button(QDialogButtonBox::Ok);

  QTimer *timer = new QTimer(this);
  connect(timer, SIGNAL(timeout()), this, SLOT(validateTimeout()));
  timer->start(100);
}

void CustomJobDialog::validateTimeout() {
  bool globalEnabled = true;

  bool thisEnabled;
  markEnabled(jobNameLineEdit, jobNameLineEdit->text().size()!=0, globalEnabled);

  thisEnabled = (QFileInfo(commandLineEdit->text()).exists());
  markEnabled(commandLineEdit, thisEnabled, globalEnabled);

  thisEnabled=QFileInfo(workingDirLineEdit->text()).exists() &&
    QFileInfo(workingDirLineEdit->text()).isDir();
  markEnabled(workingDirLineEdit, thisEnabled, globalEnabled);

  thisEnabled = (!fileModelRadioButton->isChecked() ||
                 QFileInfo(modelLineEdit->text()).exists());
  markEnabled(modelLineEdit, thisEnabled, globalEnabled);

  thisEnabled = (dataLineEdit->text().size()==0 ||
                 QFileInfo(dataLineEdit->text()).exists());
  markEnabled(dataLineEdit, thisEnabled, globalEnabled);

  if (_okButton)
    _okButton->setEnabled(globalEnabled);
}

void CustomJobDialog::setJobName(const QString& jobName)
{
    jobNameLineEdit->setText(jobName);

}

void CustomJobDialog::setProgram(const QString& program)
{
    commandLineEdit->setText(program);

}

void CustomJobDialog::setArguments(const QString& arguments)
{
    argumentsLineEdit->setText(arguments);

}

void CustomJobDialog::setModelMode(CustomJobDialog::ModelMode mode)
{
    if (mode == CURRENT) {
        currentModelRadioButton->setChecked(true);
    } else {
        fileModelRadioButton->setChecked(true);
    }
}

void CustomJobDialog::setWorkingDirectory(const QString& dir)
{
    workingDirLineEdit->setText(dir);
}

void CustomJobDialog::setModelFile(const QString& modelFile)
{
    modelLineEdit->setText(modelFile);
}

void CustomJobDialog::setDataFile(const QString& dataFile)
{
    dataLineEdit->setText(dataFile);
}

QString CustomJobDialog::jobName() const
{
    return jobNameLineEdit->text();
}

QString CustomJobDialog::program() const
{
    return commandLineEdit->text();
}

QString CustomJobDialog::arguments() const
{
    return argumentsLineEdit->text();
}

CustomJobDialog::ModelMode CustomJobDialog::modelMode() const
{
    if (currentModelRadioButton->isChecked())
        return CURRENT;
    else
        return FILE;
}

QString CustomJobDialog::workingDirectory() const
{
    return workingDirLineEdit->text();
}

QString CustomJobDialog::modelFile() const
{
    return modelLineEdit->text();
}

QString CustomJobDialog::dataFile() const
{
    return dataLineEdit->text();
}
