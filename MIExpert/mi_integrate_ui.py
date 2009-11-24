import sys, os, pickle, subprocess
from PyQt4 import QtCore, QtGui, uic


class IntegrateDialog(QtGui.QDialog):
    def __init__(self, parent=None):
        super(IntegrateDialog, self).__init__(parent)
        self.mifitDir = QtCore.QString()

        config = {
            #  data["template_image"].str = "";
            #  data["detector_constants"].str = "";
            #
            #  std::vector<std::string> sg;
            #  MIGetSpacegroups(sg);
            #  data["spacegroup_no"].radio = 0; 
            #  data["spacegroup_no"].radio_count = sg.size();
            #  data["spacegroup_no"].radio_labels = sg;
            #
            #  data["first_image"].u = 0;
            #  data["last_image"].u = 0;
            #
            #  data["integrate_resolution"].str="";
          }
        settings = QtCore.QSettings("MIFit", "MIExpert")
        appSettings = settings.value("integrate").toString()
        if not appSettings.isEmpty():
            config = pickle.loads(str(appSettings))

        uiFile = os.path.join(os.path.dirname(sys.argv[0]), "mi_integrate.ui")
        uic.loadUi(uiFile, self)

        #  connect(intensityDataBrowse, SIGNAL(clicked()),
        #          this, SLOT(selectIntensityData()));
        #  new MIBrowsePair(hardwareParametersPushButton, hardwareParametersLineEdit,"All Files (*.* *)");
        #
        #  _okButton = buttonBox->button(QDialogButtonBox::Ok);
        #
        #  QTimer *timer = new QTimer(this);
        #  connect(timer, SIGNAL(timeout()), this, SLOT(validateTimeout()));
        #  timer->start(100);

    def runJob(self):

        #settings = QtCore.QSettings("MIFit", "MIExpert")
        #settings.setValue("integrate", pickle.dumps(config))

        miexpert = os.path.join(os.path.dirname(sys.argv[0]), "MIExpert.py")
        args = [ sys.executable, miexpert, "integrate" ]

        #  QStringList args;
        #  args << MIExpertPy() << "integrate";
        #  args << "--template_image" << data["template_image"].str.c_str();
        #  if (data["detector_constants"].str.size()) {
        #    args << "--detector_constants" << buildAbsPath(data["detector_constants"].str.c_str());
        #  }
        #  args << "--spacegroup" << QString::number(data["spacegroup_no"].radio);
        #  if (data["first_image"].u != UINT_MAX) {
        #      args << "--first_image" << QString::number(data["first_image"].u);
        #  }
        #  if (data["last_image"].u != UINT_MAX) {
        #      args << "--last_image" << QString::number(data["last_image"].u);
        #  }
        #  if (data["integrate_resolution"].str.size()) {
        #    args << "--integrate_resolution" << data["integrate_resolution"].str.c_str();
        #  }

        result = subprocess.call(args)
        exit(result)

#void IntegrateDialog::validateTimeout() {
#  bool globalEnabled = true;
#
#  bool thisEnabled = !intensityData->text().isEmpty();
#  markEnabled(intensityData, thisEnabled, globalEnabled);
#
#  thisEnabled = (!resetHardwareParametersCheckBox->isChecked() || 
#                 QFileInfo(hardwareParametersLineEdit->text()).exists());
#  markEnabled(hardwareParametersLineEdit, thisEnabled, globalEnabled);
#
#  if (_okButton)
#    _okButton->setEnabled(globalEnabled);
#}
#
#
#void IntegrateDialog::InitializeFromData(const MIData &data)
#{
#  MIData dat=data;
#  // populate spacegroup combo box from dat;
#  spaceGroupComboBox->clear();
#  for (size_t i=0; i< dat["spacegroup_no"].radio_labels.size(); ++i) {
#    spaceGroupComboBox->addItem(dat["spacegroup_no"].radio_labels[i].c_str());
#  }
#  spaceGroupComboBox->setCurrentIndex(0);
#}
#
#bool IntegrateDialog::GetData(MIData &data) {
#
#  data["template_image"].str = intensityData->text().toStdString();
#
#  data["detector_constants"].str = "";
#  if (resetHardwareParametersCheckBox->isChecked())
#    data["detector_constants"].str = hardwareParametersLineEdit->text().toStdString();
#    
#  data["first_image"].u = UINT_MAX;
#  if (firstImageCheckBox->isChecked())
#    data["first_image"].u = firstImageSpinBox->value();
#  data["last_image"].u = UINT_MAX;
#  if (lastImageCheckBox->isChecked())
#    data["last_image"].u = lastImageSpinBox->value();
#
#  data["integrate_resolution"].str="";
#  if (resolutionRangeCheckBox->isChecked()) 
#    data["integrate_resolution"].str=maxResSpinBox->text().toStdString() + " " + minResSpinBox
#->text().toStdString();
#    
#  
#  data["spacegroup_no"].radio = spaceGroupComboBox->currentIndex() + 1; // first spacegroup is
# 1, first index is 0
#
#  return true;
#}
#
#void IntegrateDialog::selectIntensityData() {
#  QString str = QFileDialog::getExistingDirectory();
#  if (str.isEmpty()) {
#    return;
#  }
#  intensityData->setText(str);
#}


if __name__ == '__main__':

    app = QtGui.QApplication(sys.argv)

    dialog = IntegrateDialog()

    if 'MIFIT_DIR' in os.environ.keys():
        dialog.mifitDir = os.environ['MIFIT_DIR']

    if 'SHELX_DIR' in os.environ.keys():
        dialog.shelxDir = os.environ['SHELX_DIR']

    if dialog.exec_():
        dialog.runJob()
