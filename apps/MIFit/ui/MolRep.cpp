#include "MolRep.h"
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


MolRep::MolRep(QWidget *parent) : MIQDialog(parent) {
  setupUi(this);

  new MIBrowsePair(workingDirPushButton, workingDirLineEdit,"", true);
  new MIBrowsePair(modelPushButton, modelLineEdit,"PDB file (*.pdb *.res *.ent *.*)");
  new MIBrowsePair(dataPushButton, dataLineEdit,"Data file (*.mtz *.ref *.sca)");
  new MIBrowsePair(fixedModelPushButton, fixedModelLineEdit,"PDB file (*.pdb *.res *.ent *.*)");
  
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

void MolRep::GetInitialData(MIData &data) { // static
  data["workdir"].str="";
  data["pdbfile"].str = "";
  data["mtzfile"].str = "";
  data["fixed_pdb"].str = "";

  data["multi_search"].b = false;
  data["match_pdbin"].b = false;
  data["engine"].str= "";


  std::vector<std::string> sg;
  MIGetSpacegroups(sg);
  data["spacegroup_no"].radio = 0; 
  data["spacegroup_no"].radio_count = sg.size();
  data["spacegroup_no"].radio_labels = sg;
  data["sg_search"].b = false;

  data["copies"].u = 1;
}


void MolRep::validateTimeout() {
  bool globalEnabled = true;

  bool thisEnabled=QFileInfo(workingDirLineEdit->text()).exists() &&
    QFileInfo(workingDirLineEdit->text()).isDir();
  markEnabled(workingDirLineEdit, thisEnabled, globalEnabled);

  thisEnabled = (QFileInfo(modelLineEdit->text()).exists());
  markEnabled(modelLineEdit, thisEnabled, globalEnabled);

  thisEnabled = (QFileInfo(dataLineEdit->text()).exists());
  markEnabled(dataLineEdit, thisEnabled, globalEnabled);

  thisEnabled = (!fixedModelCheckBox->isChecked() || 
                 QFileInfo(fixedModelLineEdit->text()).exists());
  markEnabled(fixedModelLineEdit, thisEnabled, globalEnabled);

  if (_okButton)
    _okButton->setEnabled(globalEnabled);
}


void MolRep::InitializeFromData(const MIData &data)
{
  MIData dat=data;
  // populate spacegroup combo box from dat;
  spaceGroupComboBox->clear();
  for (size_t i=0; i< dat["spacegroup_no"].radio_labels.size(); ++i) {
    spaceGroupComboBox->addItem(dat["spacegroup_no"].radio_labels[i].c_str());
  }
  spaceGroupComboBox->setCurrentIndex(0);
}

bool MolRep::GetData(MIData &data) {
  data["workdir"].str=workingDirLineEdit->text().toStdString();
  data["model"].str = modelLineEdit->text().toStdString();
  data["mtzfile"].str = dataLineEdit->text().toStdString();

  data["fixed_pdb"].str = "";
  if (fixedModelCheckBox->isChecked())
    data["fixed_pdb"].str = fixedModelLineEdit->text().toStdString();

  data["spacegroup_no"].radio = 0;
  if (spaceGroupRadioButton->isChecked()) 
    data["spacegroup_no"].radio = spaceGroupComboBox->currentIndex() + 1; // first spacegroup is 1, first index is 0

  data["multi_search"].b = searchMultipleCheckBox->isChecked();
  data["match_pdbin"].b = matchInputCheckBox->isChecked();

  data["engine"].str = "molrep";
  if (phaserRadioButton->isChecked())
    data["engine"].str = "phaser";
    
  data["copies"].u = copiesSpinBox->value();
  return true;
}
