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
  
  _okButton = 0;
  QList<QAbstractButton*> buttons=buttonBox->buttons();
  for (int i=0; i < buttons.size(); ++i) {
    if (buttonBox->buttonRole(buttons[i]) == QDialogButtonBox::AcceptRole) {
      _okButton=buttons[i];
      break;
    }
  }

  QTimer *timer = new QTimer(this);
  connect(timer, SIGNAL(timeout()), this, SLOT(validateTimeout()));
  timer->start(100);
}

// create data array with all proper elements, but unset.
// this is for encapsulation purposes --- instead of having the knowlege
// about which fields will be set spread across multiple files, it's all
// contained here.

void CustomJobDialog::GetInitialData(MIData &data) { // static
  static int customJobNumber = 1;
  std::string str;

  MIConfig* config = MIConfig::Instance();
  data["jobName"].str = format("Custom job %d", customJobNumber++);

  config->Read("CustomJob/executable", str);
  data["executable"].str = str;
  config->Read("CustomJob/arguments", str);
  data["arguments"].str = str;
  bool b;
  config->Read("CustomJob/useCurrentModel", &b);
  data["useCurrentModel"].b = b;
  config->Read("CustomJob/workingDirectory", str);
  data["workingDirectory"].str = str;
  config->Read("CustomJob/modelFile", str);
  data["modelFile"].str = str;
  config->Read("CustomJob/dataFile", str);
  data["dataFile"].str = str;
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


void CustomJobDialog::InitializeFromData(const MIData &data)
{
  MIData dat=data;
  jobNameLineEdit->setText(dat["jobName"].str.c_str());
  commandLineEdit->setText(dat["executable"].str.c_str());
  argumentsLineEdit->setText(dat["arguments"].str.c_str());
  if (dat["useCurrentModel"].b) {
    currentModelRadioButton->setChecked(true);
  } else {
    fileModelRadioButton->setChecked(true);
  }
  workingDirLineEdit->setText(dat["workingDirectory"].str.c_str());
  modelLineEdit->setText(dat["modelFile"].str.c_str());
  dataLineEdit->setText(dat["dataFile"].str.c_str());
}

bool CustomJobDialog::GetData(MIData &data) {
  data["jobName"].str = jobNameLineEdit->text().toStdString();
  data["executable"].str = commandLineEdit->text().toStdString();
  data["arguments"].str = argumentsLineEdit->text().toStdString();
  data["useCurrentModel"].b = currentModelRadioButton->isChecked();
  data["workingDirectory"].str = workingDirLineEdit->text().toStdString();
  data["modelFile"].str = modelLineEdit->text().toStdString();
  data["dataFile"].str = dataLineEdit->text().toStdString();
  MIConfig* config = MIConfig::Instance();
  config->Write("CustomJob/executable", data["executable"].str);
  config->Write("CustomJob/arguments", data["arguments"].str);
  config->Write("CustomJob/useCurrentModel", data["useCurrentModel"].b);
  config->Write("CustomJob/workingDirectory", data["workingDirectory"].str);
  config->Write("CustomJob/modelFile", data["modelFile"].str);
  config->Write("CustomJob/dataFile", data["dataFile"].str);
  return true;
}
