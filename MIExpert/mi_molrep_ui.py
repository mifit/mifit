import sys, os, pickle, subprocess
from PyQt4 import QtCore, QtGui, uic


class MolRepDialog(QtGui.QDialog):
    def __init__(self, parent=None):
        super(MolRepDialog, self).__init__(parent)
        self.mifitDir = QtCore.QString()

        config = {
            #  data["workdir"].str="";
            #  data["pdbfile"].str = "";
            #  data["mtzfile"].str = "";
            #  data["fixed_pdb"].str = "";
            #
            #  data["multi_search"].b = false;
            #  data["match_pdbin"].b = false;
            #  data["engine"].str= "";
            #
            #
            #  std::vector<std::string> sg;
            #  MIGetSpacegroups(sg);
            #  data["spacegroup_no"].radio = 0; 
            #  data["spacegroup_no"].radio_count = sg.size();
            #  data["spacegroup_no"].radio_labels = sg;
            #  data["sg_search"].b = false;
            #
            #  data["copies"].u = 1;
          }
        settings = QtCore.QSettings("MIFit", "MIExpert")
        appSettings = settings.value("molrep").toString()
        if not appSettings.isEmpty():
            config = pickle.loads(str(appSettings))

        uiFile = os.path.join(os.path.dirname(sys.argv[0]), "mi_molrep.ui")
        uic.loadUi(uiFile, self)

        #  new MIBrowsePair(workingDirPushButton, workingDirLineEdit,"", true);
        #  new MIBrowsePair(modelPushButton, modelLineEdit,"PDB file (*.pdb *.res *.ent *.*)");
        #  new MIBrowsePair(dataPushButton, dataLineEdit,"Data file (*.mtz *.ref *.sca)");
        #  new MIBrowsePair(fixedModelPushButton, fixedModelLineEdit,"PDB file (*.pdb *.res *.ent *.*)");
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

    def runJob(self):

        #settings = QtCore.QSettings("MIFit", "MIExpert")
        #settings.setValue("molrep", pickle.dumps(config))

        miexpert = os.path.join(os.path.dirname(sys.argv[0]), "MIExpert.py")
        args = [ sys.executable, miexpert, "molrep" ]

        #  QStringList args;
        #  args << MIExpertPy() << "molrep";
        #  args << "--engine" << data["engine"].str.c_str();
        #  if (data["spacegroup_no"].radio != 0) {
        #      args << "--spacegroup" << QString::number(data["spacegroup_no"].radio);
        #  } else if (data["sg_search"].b) {
        #    args << "--sg_search" << "yes";
        #  }
        #  args << "--pdbfile" << buildAbsPath(data["model"].str.c_str())
        #          << "--mtzfile" << buildAbsPath(data["mtzfile"].str.c_str())
        #          << "--workdir" << buildAbsPath(data["workdir"].str.c_str())
        #          << "--multi_search" << (data["multi_search"].b ? "yes": "no")
        #          << "--match_pdbin" << (data["match_pdbin"].b ? "yes": "no")
        #          << "--copies" << QString::number(data["copies"].u);
        #  if (!data["fixed_pdb"].str.empty())
        #    args << "--fixed_pdb" << buildAbsPath(data["fixed_pdb"].str.c_str());

        result = subprocess.call(args)
        exit(result)

#void MolRep::validateTimeout() {
#  bool globalEnabled = true;
#
#  bool thisEnabled=QFileInfo(workingDirLineEdit->text()).exists() &&
#    QFileInfo(workingDirLineEdit->text()).isDir();
#  markEnabled(workingDirLineEdit, thisEnabled, globalEnabled);
#
#  thisEnabled = (QFileInfo(modelLineEdit->text()).exists());
#  markEnabled(modelLineEdit, thisEnabled, globalEnabled);
#
#  thisEnabled = (QFileInfo(dataLineEdit->text()).exists());
#  markEnabled(dataLineEdit, thisEnabled, globalEnabled);
#
#  thisEnabled = (!fixedModelCheckBox->isChecked() || 
#                 QFileInfo(fixedModelLineEdit->text()).exists());
#  markEnabled(fixedModelLineEdit, thisEnabled, globalEnabled);
#
#  if (_okButton)
#    _okButton->setEnabled(globalEnabled);
#}
#
#
#void MolRep::InitializeFromData(const MIData &data)
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
#bool MolRep::GetData(MIData &data) {
#  data["workdir"].str=workingDirLineEdit->text().toStdString();
#  data["model"].str = modelLineEdit->text().toStdString();
#  data["mtzfile"].str = dataLineEdit->text().toStdString();
#
#  data["fixed_pdb"].str = "";
#  if (fixedModelCheckBox->isChecked())
#    data["fixed_pdb"].str = fixedModelLineEdit->text().toStdString();
#
#  data["spacegroup_no"].radio = 0;
#  if (spaceGroupRadioButton->isChecked()) 
#    data["spacegroup_no"].radio = spaceGroupComboBox->currentIndex() + 1; // first spacegroup is 1, first index is 0
#
#  data["multi_search"].b = searchMultipleCheckBox->isChecked();
#  data["match_pdbin"].b = matchInputCheckBox->isChecked();
#
#  data["engine"].str = "molrep";
#  if (phaserRadioButton->isChecked())
#    data["engine"].str = "phaser";
#    
#  data["copies"].u = copiesSpinBox->value();
#  return true;
#}

if __name__ == '__main__':

    app = QtGui.QApplication(sys.argv)

    dialog = MolRepDialog()

    if 'MIFIT_DIR' in os.environ.keys():
        dialog.mifitDir = os.environ['MIFIT_DIR']

    if 'SHELX_DIR' in os.environ.keys():
        dialog.shelxDir = os.environ['SHELX_DIR']

    if dialog.exec_():
        dialog.runJob()
