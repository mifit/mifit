import sys, os, pickle, subprocess
from PyQt4 import QtCore, QtGui, uic


class LigandOverlapDialog(QtGui.QDialog):
    def __init__(self, parent=None):
        super(LigandOverlapDialog, self).__init__(parent)
        self.mifitDir = QtCore.QString()

        config = {
            #  data["workdir"].str="";
            #  data["pdbdir"].str = "";
            #  data["targetpdb"].str = "";
            #  data["x"].f = FLT_MIN;
            #  data["y"].f = FLT_MIN;
            #  data["z"].f = FLT_MIN;
          }
        settings = QtCore.QSettings("MIFit", "MIExpert")
        appSettings = settings.value("ligandoverlap").toString()
        if not appSettings.isEmpty():
            config = pickle.loads(str(appSettings))

        uiFile = os.path.join(os.path.dirname(sys.argv[0]), "mi_ligandoverlap.ui")
        uic.loadUi(uiFile, self)

        #  new MIBrowsePair(workingDirPushButton, workingDirLineEdit,"", true);
        #  new MIBrowsePair(structurePushButton, structureLineEdit,"", true);
        #  new MIBrowsePair(targetFilePushButton, targetFileLineEdit,"PDB file (*.pdb *.res *.ent *.*)");
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
        #settings.setValue("ligandoverlap", pickle.dumps(config))

        miexpert = os.path.join(os.path.dirname(sys.argv[0]), "MIExpert.py")
        args = [ sys.executable, miexpert, "ligandoverlap" ]

        #  QStringList args;
        #  args << MIExpertPy() << "ligandoverlap";
        #  args << "--workdir" << buildAbsPath(data["workdir"].str.c_str());
        #  args << "--pdbdir" << buildAbsPath(data["pdbdir"].str.c_str());
        #  args << "--targetpdb" << buildAbsPath(data["targetpdb"].str.c_str());
        #  args << "--targetsite" << QString::number(data["x"].f, 'f', 3)
        #          << QString::number(data["y"].f, 'f', 3)
        #          << QString::number(data["z"].f, 'f', 3);

        result = subprocess.call(args)
        exit(result)

#
#void CocrystalSuperPos::validateTimeout() {
#  bool globalEnabled = true;
#
#  bool thisEnabled=QFileInfo(workingDirLineEdit->text()).exists() &&
#    QFileInfo(workingDirLineEdit->text()).isDir();
#  markEnabled(workingDirLineEdit, thisEnabled, globalEnabled);
#
#  thisEnabled=QFileInfo(structureLineEdit->text()).exists() && 
#    QFileInfo(structureLineEdit->text()).isDir();
#  markEnabled(structureLineEdit, thisEnabled, globalEnabled);
#
#  thisEnabled=QFileInfo(targetFileLineEdit->text()).exists();
#  markEnabled(targetFileLineEdit, thisEnabled, globalEnabled);
#
#  if (_okButton)
#    _okButton->setEnabled(globalEnabled);
#}
#
#bool CocrystalSuperPos::GetData(MIData &data) {
#  data["workdir"].str=workingDirLineEdit->text().toStdString();
#  data["pdbdir"].str = structureLineEdit->text().toStdString();
#  data["targetpdb"].str = targetFileLineEdit->text().toStdString();
#  data["x"].f = xSpinBox->value();
#  data["y"].f = ySpinBox->value();
#  data["z"].f = zSpinBox->value();
#  return true;
#}

if __name__ == '__main__':

    app = QtGui.QApplication(sys.argv)

    dialog = LigandOverlapDialog()

    if 'MIFIT_DIR' in os.environ.keys():
        dialog.mifitDir = os.environ['MIFIT_DIR']

    if 'SHELX_DIR' in os.environ.keys():
        dialog.shelxDir = os.environ['SHELX_DIR']

    if dialog.exec_():
        dialog.runJob()
