#include "IntegrateDialog.h"
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


IntegrateDialog::IntegrateDialog(QWidget *parent) : MIQDialog(parent) {
  setupUi(this);

  connect(intensityDataBrowse, SIGNAL(clicked()),
          this, SLOT(selectIntensityData()));
  new MIBrowsePair(hardwareParametersPushButton, hardwareParametersLineEdit,"All Files (*.* *)");

  _okButton = buttonBox->button(QDialogButtonBox::Ok);

  QTimer *timer = new QTimer(this);
  connect(timer, SIGNAL(timeout()), this, SLOT(validateTimeout()));
  timer->start(100);
}

// create data array with all proper elements, but unset.
// this is for encapsulation purposes --- instead of having the knowlege
// about which fields will be set spread across multiple files, it's all
// contained here.

void IntegrateDialog::GetInitialData(MIData &data) { // static
  data["template_image"].str = "";
  data["detector_constants"].str = "";

  std::vector<std::string> sg;
  MIGetSpacegroups(sg);
  data["spacegroup_no"].radio = 0; 
  data["spacegroup_no"].radio_count = sg.size();
  data["spacegroup_no"].radio_labels = sg;

  data["first_image"].u = 0;
  data["last_image"].u = 0;

  data["integrate_resolution"].str="";
}

void IntegrateDialog::validateTimeout() {
  bool globalEnabled = true;

  bool thisEnabled = !intensityData->text().isEmpty();
  markEnabled(intensityData, thisEnabled, globalEnabled);

  thisEnabled = (!resetHardwareParametersCheckBox->isChecked() || 
                 QFileInfo(hardwareParametersLineEdit->text()).exists());
  markEnabled(hardwareParametersLineEdit, thisEnabled, globalEnabled);

  if (_okButton)
    _okButton->setEnabled(globalEnabled);
}


void IntegrateDialog::InitializeFromData(const MIData &data)
{
  MIData dat=data;
  // populate spacegroup combo box from dat;
  spaceGroupComboBox->clear();
  for (size_t i=0; i< dat["spacegroup_no"].radio_labels.size(); ++i) {
    spaceGroupComboBox->addItem(dat["spacegroup_no"].radio_labels[i].c_str());
  }
  spaceGroupComboBox->setCurrentIndex(0);
}

bool IntegrateDialog::GetData(MIData &data) {

  data["template_image"].str = intensityData->text().toStdString();

  data["detector_constants"].str = "";
  if (resetHardwareParametersCheckBox->isChecked())
    data["detector_constants"].str = hardwareParametersLineEdit->text().toStdString();
    
  data["first_image"].u = UINT_MAX;
  if (firstImageCheckBox->isChecked())
    data["first_image"].u = firstImageSpinBox->value();
  data["last_image"].u = UINT_MAX;
  if (lastImageCheckBox->isChecked())
    data["last_image"].u = lastImageSpinBox->value();

  data["integrate_resolution"].str="";
  if (resolutionRangeCheckBox->isChecked()) 
    data["integrate_resolution"].str=maxResSpinBox->text().toStdString() + " " + minResSpinBox->text().toStdString();
    
  
  data["spacegroup_no"].radio = spaceGroupComboBox->currentIndex() + 1; // first spacegroup is 1, first index is 0

  return true;
}

void IntegrateDialog::selectIntensityData() {
  QString str = QFileDialog::getExistingDirectory();
  if (str.isEmpty()) {
    return;
  }
  intensityData->setText(str);
}
