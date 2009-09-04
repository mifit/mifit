#include "NCSModeling.h"
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


NCSModeling::NCSModeling(QWidget *parent) : MIQDialog(parent) {
  setupUi(this);

  new MIBrowsePair(workingDirPushButton, workingDirLineEdit,"", true);
  new MIBrowsePair(modelPushButton, modelLineEdit,"Exclude file (*.txt)");
  new MIBrowsePair(dataPushButton, dataLineEdit,"Data files (*.mtz)");
  new MIBrowsePair(maskAdditionsPushButton, maskAdditionsLineEdit,"Include file (*.txt)");
  
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

  connect(dataCheckBox, SIGNAL(toggled(bool)), 
          this, SLOT(togglePhaseOptions(bool)));
  togglePhaseOptions(false);
}


void NCSModeling::togglePhaseOptions(bool state) {
  phaseOptions->setEnabled(state);
  maskAdditionsCheckBox->setEnabled(state);
  maskAdditionsLineEdit->setEnabled(state);
  maskAdditionsPushButton->setEnabled(state);
}

// create data array with all proper elements, but unset.
// this is for encapsulation purposes --- instead of having the knowlege
// about which fields will be set spread across multiple files, it's all
// contained here.

void NCSModeling::GetInitialData(MIData &data) { // static
  data["workdir"].str="";
  data["model"].str = "";
  data["mtzdata"].str = "";
  data["maskadditions"].str="";
  data["chainid"].str;
  data["phasecalc"].radio;
  data["phasecalc"].radio_count=3;
}


void NCSModeling::validateTimeout() {
  bool globalEnabled = true;

  bool thisEnabled=QFileInfo(workingDirLineEdit->text()).exists() &&
    QFileInfo(workingDirLineEdit->text()).isDir();
  markEnabled(workingDirLineEdit, thisEnabled, globalEnabled);

  thisEnabled=(!modelCheckBox->isChecked() || 
               QFileInfo(modelLineEdit->text()).exists());
  markEnabled(modelLineEdit, thisEnabled, globalEnabled);

  thisEnabled = (!dataCheckBox->isChecked() || 
                 QFileInfo(dataLineEdit->text()).exists());
  markEnabled(dataLineEdit, thisEnabled, globalEnabled);

  thisEnabled = (!maskAdditionsCheckBox->isChecked() || 
                 QFileInfo(maskAdditionsLineEdit->text()).exists());
  markEnabled(maskAdditionsLineEdit, thisEnabled, globalEnabled);

  if (_okButton)
    _okButton->setEnabled(globalEnabled);
}


void NCSModeling::InitializeFromData(const MIData &)
{
  //FIXME: pre-populate chain/data?
}

bool NCSModeling::GetData(MIData &data) {
  data["workdir"].str = workingDirLineEdit->text().toStdString();

  data["model"].str = "";
  if (modelCheckBox->isChecked())
    data["model"].str = modelLineEdit->text().toStdString();

  data["mtzdata"].str = "";
  if (dataCheckBox->isChecked())
    data["mtzdata"].str = dataLineEdit->text().toStdString();


  data["chainid"].str=primaryChainIDLineEdit->text().toStdString();

  data["maskadditions"].str = "";
  if (maskAdditionsCheckBox->isChecked()) 
    data["maskadditions"].str=maskAdditionsLineEdit->text().toStdString();

  if (calculatedRadioButton->isChecked())
    data["phasecalc"].radio=0;
  else if (averageRadioButton->isChecked())
    data["phasecalc"].radio=1;
  else if (combinedRadioButton->isChecked())
    data["phasecalc"].radio=2;
  return true;
}
