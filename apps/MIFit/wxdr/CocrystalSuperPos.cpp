#include "CocrystalSuperPos.h"
#include "MIBrowsePair.h"

#include "corelib.h"
#include "maplib.h"

#include <vector>
#include <string>

#include <QDialogButtonBox>
#include <QAbstractButton>
#include <QList>
#include <QTimer>
#include <QFileInfo>
#include <QFileDialog>


CocrystalSuperPos::CocrystalSuperPos(QWidget *parent) : MIQDialog(parent) {
  setupUi(this);

  new MIBrowsePair(workingDirPushButton, workingDirLineEdit,"", true);
  new MIBrowsePair(structurePushButton, structureLineEdit,"", true);
  new MIBrowsePair(targetFilePushButton, targetFileLineEdit,"PDB file (*.pdb *.res *.ent *.*)");
  
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

void CocrystalSuperPos::GetInitialData(MIData &data) { // static
  data["workdir"].str="";
  data["pdbdir"].str = "";
  data["targetpdb"].str = "";
  data["x"].f = FLT_MIN;
  data["y"].f = FLT_MIN;
  data["z"].f = FLT_MIN;
}


void CocrystalSuperPos::validateTimeout() {
  bool globalEnabled = true;

  bool thisEnabled=QFileInfo(workingDirLineEdit->text()).exists() &&
    QFileInfo(workingDirLineEdit->text()).isDir();
  markEnabled(workingDirLineEdit, thisEnabled, globalEnabled);

  thisEnabled=QFileInfo(structureLineEdit->text()).exists() && 
    QFileInfo(structureLineEdit->text()).isDir();
  markEnabled(structureLineEdit, thisEnabled, globalEnabled);

  thisEnabled=QFileInfo(targetFileLineEdit->text()).exists();
  markEnabled(targetFileLineEdit, thisEnabled, globalEnabled);

  if (_okButton)
    _okButton->setEnabled(globalEnabled);
}


void CocrystalSuperPos::InitializeFromData(const MIData &)
{
}

bool CocrystalSuperPos::GetData(MIData &data) {
  data["workdir"].str=workingDirLineEdit->text().toStdString();
  data["pdbdir"].str = structureLineEdit->text().toStdString();
  data["targetpdb"].str = targetFileLineEdit->text().toStdString();
  data["x"].f = xSpinBox->value();
  data["y"].f = ySpinBox->value();
  data["z"].f = zSpinBox->value();
  return true;
}
