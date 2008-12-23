#include "SadPhasing.h"
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


SadPhasing::SadPhasing(QWidget *parent) : MIQDialog(parent) {
  setupUi(this);

  new MIBrowsePair(workingDirectoryPushButton, workingDirectoryLineEdit,"", true);
  new MIBrowsePair(intensityDataPushButton, intensityDataLineEdit,"Intensity files (*.ref *.sca *.mtz)");
  new MIBrowsePair(sitesFromFilePushButton, sitesFromFileLineEdit,"PDB file (*.pdb *.res *.ent *.*)");
  
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

void SadPhasing::GetInitialData(MIData &data) { // static
  data["workdir"].str="";
  data["saddatafile"].str = "";

  data["sitefile"].str = "none";
  data["scatterer"].str = "";

  data["sitenumber"].u=UINT_MAX;
  data["ssnumber"].u=UINT_MAX;
  data["separation"].f=FLT_MIN;
  data["solventfraction"].f=FLT_MIN;
  data["siterefinemethod"].str="";
  data["bothhands"].b = false;

  std::vector<std::string> sg;
  MIGetSpacegroups(sg);
  data["spacegroup_no"].radio = 0; 
  data["spacegroup_no"].radio_count = sg.size();
  data["spacegroup_no"].radio_labels = sg;
  data["change_spacegroup"].b = false;
}


void SadPhasing::validateTimeout() {
  bool globalEnabled = true;

  bool thisEnabled=QFileInfo(workingDirectoryLineEdit->text()).exists() &&
    QFileInfo(workingDirectoryLineEdit->text()).isDir();
  markEnabled(workingDirectoryLineEdit, thisEnabled, globalEnabled);

  thisEnabled=QFileInfo(intensityDataLineEdit->text()).exists();
  markEnabled(intensityDataLineEdit, thisEnabled, globalEnabled);

  thisEnabled = (!sitesFromFileRadioButton->isChecked() || 
                 QFileInfo(sitesFromFileLineEdit->text()).exists());
  markEnabled(sitesFromFileLineEdit, thisEnabled, globalEnabled);

  thisEnabled = scattererTypeLineEdit->text().size();
  markEnabled(scattererTypeLineEdit, thisEnabled, globalEnabled);


  if (_okButton)
    _okButton->setEnabled(globalEnabled);
}


void SadPhasing::InitializeFromData(const MIData &data)
{
  MIData dat=data;
  // populate spacegroup combo box from dat;
  spaceGroupComboBox->clear();
  for (size_t i=0; i< dat["spacegroup_no"].radio_labels.size(); ++i) {
    spaceGroupComboBox->addItem(dat["spacegroup_no"].radio_labels[i].c_str());
  }
  spaceGroupComboBox->setCurrentIndex(0);
}

bool SadPhasing::GetData(MIData &data) {
  data["workdir"].str=workingDirectoryLineEdit->text().toStdString();
  data["saddatafile"].str = intensityDataLineEdit->text().toStdString();

  data["sitefile"].str = "none";
  if (sitesFromFileRadioButton->isChecked())
    data["sitefile"].str = sitesFromFileLineEdit->text().toStdString();

  data["scatterer"].str = scattererTypeLineEdit->text().toStdString();

  data["sitenumber"].u=numSitesSpinBox->value();
  data["ssnumber"].u=numDisulfidesSpinBox->value();
  data["separation"].f=minScattererSeparationSpinBox->value();
  data["solventfraction"].f=solventFractionSpinBox->value();
  if (bp3RadioButton->isChecked())
    data["siterefinemethod"].str="bp3";
  else
    data["siterefinemethod"].str="mlphare";
  data["bothhands"].b = phaseBothSiteEnantiomorphsCheckBox->isChecked();

  data["change_spacegroup"].b = changeSpaceGroupCheckBox->isChecked();
  data["spacegroup_no"].radio = spaceGroupComboBox->currentIndex() + 1; // first spacegroup is 1, first index is 0
  return true;
}
