import sys, os, pickle, subprocess
from PyQt4 import QtCore, QtGui, uic


class NCSModelingDialog(QtGui.QDialog):
    def __init__(self, parent=None):
        super(NCSModelingDialog, self).__init__(parent)
        self.mifitDir = QtCore.QString()

        config = {
            #  data["workdir"].str = "";
            #  data["mtzfile"].str = "";
            #  data["pdbfile"].str = "";
            #  data["libfile"].str = "none";
            #  data["seqfile"].str = "none";
            #  data["templatefile"].str = "none";
            #  data["datalogfile"].str = "none";
            #
            #  data["cif_write"].b = false; // mmCIF report
            #  data["map_write"].b = false;
            #  data["map_border"].f = FLT_MIN;
            #  data["text_report"].b = false;
            #  data["html_report"].b = false;
            #  data["hkl_write"].b = false; // mmCIF data
            #
            #  data["html_report"].b = false;
            #  data["report_title"].str = "";
            #  data["image1"].str = "";
            #
            #  data["rootname"].str = "";
          }
        settings = QtCore.QSettings("MIFit", "MIExpert")
        appSettings = settings.value("ncsmodeler").toString()
        if not appSettings.isEmpty():
            config = pickle.loads(str(appSettings))

        uiFile = os.path.join(os.path.dirname(sys.argv[0]), "mi_ncsmodeler.ui")
        uic.loadUi(uiFile, self)

        #  new MIBrowsePair(workingDirPushButton, workingDirLineEdit,"", true);
        #  new MIBrowsePair(modelPushButton, modelLineEdit,"Exclude file (*.txt)");
        #  new MIBrowsePair(dataPushButton, dataLineEdit,"Data files (*.mtz)");
        #  new MIBrowsePair(maskAdditionsPushButton, maskAdditionsLineEdit,"Include file (*.txt)");
        #  
        #  _okButton = 0;
        #  QList<QAbstractButton*> buttons=buttonBox->buttons();
        #  for (int i=0; i < buttons.size(); ++i) {
        #    if (buttonBox->buttonRole(buttons[i]) == QDialogButtonBox::AcceptRole) {
        #      _okButton=buttons[i];
        #      break;
        #    }
        #  }
        #
        #  QTimer *timer = new QTimer(this);
        #  connect(timer, SIGNAL(timeout()), this, SLOT(validateTimeout()));
        #  timer->start(100);
        #
        #  connect(dataCheckBox, SIGNAL(toggled(bool)), 
        #          this, SLOT(togglePhaseOptions(bool)));
        #  togglePhaseOptions(false);

#void NCSModeling::togglePhaseOptions(bool state) {
#  phaseOptions->setEnabled(state);
#  maskAdditionsCheckBox->setEnabled(state);
#  maskAdditionsLineEdit->setEnabled(state);
#  maskAdditionsPushButton->setEnabled(state);
#}
#
#// create data array with all proper elements, but unset.
#// this is for encapsulation purposes --- instead of having the knowlege
#// about which fields will be set spread across multiple files, it's all
#// contained here.
#
#void NCSModeling::GetInitialData(MIData &data) { // static
#  data["workdir"].str="";
#  data["model"].str = "";
#  data["mtzdata"].str = "";
#  data["maskadditions"].str="";
#  data["chainid"].str;
#  data["phasecalc"].radio;
#  data["phasecalc"].radio_count=3;
#}
#
#
#void NCSModeling::validateTimeout() {
#  bool globalEnabled = true;
#
#  bool thisEnabled=QFileInfo(workingDirLineEdit->text()).exists() &&
#    QFileInfo(workingDirLineEdit->text()).isDir();
#  markEnabled(workingDirLineEdit, thisEnabled, globalEnabled);
#
#  thisEnabled=(!modelCheckBox->isChecked() || 
#               QFileInfo(modelLineEdit->text()).exists());
#  markEnabled(modelLineEdit, thisEnabled, globalEnabled);
#
#  thisEnabled = (!dataCheckBox->isChecked() || 
#                 QFileInfo(dataLineEdit->text()).exists());
#  markEnabled(dataLineEdit, thisEnabled, globalEnabled);
#
#  thisEnabled = (!maskAdditionsCheckBox->isChecked() || 
#                 QFileInfo(maskAdditionsLineEdit->text()).exists());
#  markEnabled(maskAdditionsLineEdit, thisEnabled, globalEnabled);
#
#  if (_okButton)
#    _okButton->setEnabled(globalEnabled);
#}
#
#
#void NCSModeling::InitializeFromData(const MIData &)
#{
#  //FIXME: pre-populate chain/data?
#}
#
#bool NCSModeling::GetData(MIData &data) {
#  data["workdir"].str = workingDirLineEdit->text().toStdString();
#
#  data["model"].str = "";
#  if (modelCheckBox->isChecked())
#    data["model"].str = modelLineEdit->text().toStdString();
#
#  data["mtzdata"].str = "";
#  if (dataCheckBox->isChecked())
#    data["mtzdata"].str = dataLineEdit->text().toStdString();
#
#
#  data["chainid"].str=primaryChainIDLineEdit->text().toStdString();
#
#  data["maskadditions"].str = "";
#  if (maskAdditionsCheckBox->isChecked()) 
#    data["maskadditions"].str=maskAdditionsLineEdit->text().toStdString();
#
#  if (calculatedRadioButton->isChecked())
#    data["phasecalc"].radio=0;
#  else if (averageRadioButton->isChecked())
#    data["phasecalc"].radio=1;
#  else if (combinedRadioButton->isChecked())
#    data["phasecalc"].radio=2;
#  return true;
#}

    def runJob(self):

        #settings = QtCore.QSettings("MIFit", "MIExpert")
        #settings.setValue("ncsmodeler", pickle.dumps(config))

        miexpert = os.path.join(os.path.dirname(sys.argv[0]), "MIExpert.py")
        args = [ sys.executable, miexpert, "ncsmodeler" ]

        #  //Write a PDB of current model
        #  MIGLWidget* doc = MIMainWindow::instance()->currentMIGLWidget();
        #  if (!doc)
        #    return;
        #  Molecule* model = doc->GetDisplaylist()->GetCurrentModel();
        #  QDir workdir(data["workdir"].str.c_str());
        #  QString pdbout = buildAbsPath(workdir.absoluteFilePath(
        #          QString("%1_out.pdb").arg(QFileInfo(doc->GetTitle().c_str()).fileName())));
        #  model->SavePDBFile(pdbout.toAscii().constData());
        #
        #  QStringList args;
        #  args << MIExpertPy() << "ncsmodeler"
        #          << "--pdbfile" << pdbout
        #          << "--targetchain" << data["chainid"].str.c_str();
        #  if (!data["model"].str.empty()) {
        #    args << "--preserve_model" << buildAbsPath(data["model"].str.c_str());
        #  }
        #  if (!data["maskadditions"].str.empty()) {
        #    args << "--maskextras" << buildAbsPath(data["maskadditions"].str.c_str());
        #  }
        #  if (!data["mtzdata"].str.empty()) {
        #    args << "--mtzfile" << buildAbsPath(data["mtzdata"].str.c_str());
        #    switch (data["phasecalc"].radio) {
        #      case 0:
        #        break;
        #      case 1:
        #        args << "--phase_prob" << "yes";
        #        break;
        #      case 2:
        #        args << "--phase_comb" << "yes";
        #        break;
        #    }
        #  }

        result = subprocess.call(args)
        exit(result)


if __name__ == '__main__':

    app = QtGui.QApplication(sys.argv)

    dialog = NCSModelingDialog()

    if 'MIFIT_DIR' in os.environ.keys():
        dialog.mifitDir = os.environ['MIFIT_DIR']

    if 'SHELX_DIR' in os.environ.keys():
        dialog.shelxDir = os.environ['SHELX_DIR']

    if dialog.exec_():
        dialog.runJob()
