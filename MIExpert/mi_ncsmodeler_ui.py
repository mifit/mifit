import sys, os, pickle, subprocess
from PyQt4 import QtCore, QtGui, uic
import mifit

def buildAbsPath(path):
    return str(QtCore.QFileInfo(path).absoluteFilePath())

def boolYesNo(value):
    if value:
        return "yes"
    else:
        return "no"

def markEnabled(w, thisEnabled, globalEnabled):
    w.setAutoFillBackground(True)
    if not thisEnabled:
        pal = w.palette()
        pal.setColor(QtGui.QPalette.Normal, QtGui.QPalette.Base, QtGui.QColor(255, 255, 128))
        w.setPalette(pal)
    else:
        w.setPalette(QtGui.QPalette())

    return globalEnabled and thisEnabled


class NCSModelingDialog(QtGui.QDialog):
    def __init__(self, parent=None):
        super(NCSModelingDialog, self).__init__(parent)

        data = {
            "workdir" : '',
            "model" : '',
            "mtzdata" : '',
            "maskadditions" : '',
            "chainid" : '',
            "phasecalc" : 0,

          }
        settings = QtCore.QSettings("MIFit", "MIExpert")
        appSettings = settings.value("ncsmodeler").toString()
        if not appSettings.isEmpty():
            data = pickle.loads(str(appSettings))

        uiFile = os.path.join(os.path.dirname(sys.argv[0]), "mi_ncsmodeler.ui")
        uic.loadUi(uiFile, self)

        self.workingDirLineEdit.setText(data["workdir"])

        self.modelCheckBox.setChecked(False)
        if len(data["model"]) > 0:
            self.modelCheckBox.setChecked(True)
            self.modelLineEdit.setText(data["model"])

        self.dataCheckBox.setChecked(False)
        if len(data["mtzdata"]) > 0:
            self.dataCheckBox.setChecked(True)
            self.dataLineEdit.setText(data["mtzdata"])

        self.primaryChainIDLineEdit.setText(data["chainid"])

        self.maskAdditionsCheckBox.setChecked(False)
        if len(data["maskadditions"]) > 0:
            self.maskAdditionsCheckBox.setChecked(True)
            self.maskAdditionsLineEdit.setText(data["maskadditions"])

        if data["phasecalc"] == 1:
            self.averageRadioButton.setChecked(True)
        elif data["phasecalc"] == 2:
            self.combinedRadioButton.setChecked(True)
        else:
            self.calculatedRadioButton.setChecked(True)

        self.timer = QtCore.QTimer(self)
        self.connect(self.timer, QtCore.SIGNAL("timeout()"),
                     self, QtCore.SLOT("validateTimeout()"))
        self.timer.start(100)

        self.connect(self.dataCheckBox, QtCore.SIGNAL("toggled(bool)"),
                     self, QtCore.SLOT("togglePhaseOptions(bool)"))
        self.togglePhaseOptions(False)

    @QtCore.pyqtSlot()
    def on_workingDirPushButton_clicked(self):
        dir = QtGui.QFileDialog.getExistingDirectory(None, "", self.workingDirLineEdit.text())
        if not dir.isEmpty():
            self.workingDirLineEdit.setText(dir)

    def browseFile(self, lineEdit, filter):
        f = QtGui.QFileDialog.getOpenFileName(None, "Select a file", lineEdit.text(), filter)
        if not f.isEmpty():
            lineEdit.setText(f)

    @QtCore.pyqtSlot()
    def on_modelPushButton_clicked(self):
        self.browseFile(self.modelLineEdit, "Exclude file (*.txt)")

    @QtCore.pyqtSlot()
    def on_dataPushButton_clicked(self):
        self.browseFile(self.dataLineEdit, "Data files (*.mtz)")

    @QtCore.pyqtSlot()
    def on_maskAdditionsPushButton_clicked(self):
        self.browseFile(self.maskAdditionsLineEdit, "Include file (*.txt)")

    @QtCore.pyqtSlot(bool)
    def togglePhaseOptions(self, enable):
        self.phaseOptions.setEnabled(enable)
        self.maskAdditionsCheckBox.setEnabled(enable)
        self.maskAdditionsLineEdit.setEnabled(enable)
        self.maskAdditionsPushButton.setEnabled(enable)

    @QtCore.pyqtSlot()
    def validateTimeout(self):
        globalEnabled = True

        thisEnabled = QtCore.QFileInfo(self.workingDirLineEdit.text()).exists() and QtCore.QFileInfo(self.workingDirLineEdit.text()).isDir()
        globalEnabled = markEnabled(self.workingDirLineEdit, thisEnabled, globalEnabled)

        thisEnabled = (not self.modelCheckBox.isChecked() or QFileInfo(self.modelLineEdit.text()).exists())
        globalEnabled = markEnabled(self.modelLineEdit, thisEnabled, globalEnabled)

        thisEnabled = (not self.dataCheckBox.isChecked() or QtCore.QFileInfo(self.dataLineEdit.text()).exists())
        globalEnabled = markEnabled(self.dataLineEdit, thisEnabled, globalEnabled)

        thisEnabled = (not self.maskAdditionsCheckBox.isChecked() or QtCore.QFileInfo(self.maskAdditionsLineEdit.text()).exists())
        globalEnabled = markEnabled(self.maskAdditionsLineEdit, thisEnabled, globalEnabled)

        okButton = self.buttonBox.button(QtGui.QDialogButtonBox.Ok)
        if okButton:
            okButton.setEnabled(globalEnabled)


    def runJob(self):

        data = {
            "workdir" : '',
            "model" : '',
            "mtzdata" : '',
            "maskadditions" : '',
            "chainid" : '',
            "phasecalc" : 0,
          }
        data["workdir"] = str(self.workingDirLineEdit.text())

        data["model"] = ""
        if self.modelCheckBox.isChecked():
            data["model"] = str(self.modelLineEdit.text())

        data["mtzdata"] = ""
        if self.dataCheckBox.isChecked():
            data["mtzdata"] = str(self.dataLineEdit.text())


        data["chainid"] = str(self.primaryChainIDLineEdit.text())

        data["maskadditions"] = ""
        if self.maskAdditionsCheckBox.isChecked():
            data["maskadditions"] = str(self.maskAdditionsLineEdit.text())

        if self.averageRadioButton.isChecked():
            data["phasecalc"] = 1
        elif self.combinedRadioButton.isChecked():
            data["phasecalc"] = 2
        else:
            data["phasecalc"] = 0

        settings = QtCore.QSettings("MIFit", "MIExpert")
        settings.setValue("ncsmodeler", pickle.dumps(data))

        # Write a PDB of current model
        workdir = QtCore.QDir(data["workdir"])
        pdbout = str(buildAbsPath(workdir.absoluteFilePath(QtCore.QString("model.pdb"))));
        mifit.writeCurrentModel(pdbout)

        mifit.setJobWorkDir(str(QtCore.QFileInfo(data['workdir']).absoluteFilePath()))
        miexpert = os.path.join(os.path.dirname(sys.argv[0]), "MIExpert.py")
        args = [ sys.executable, miexpert, "ncsmodeler" ]

        args += [ "--pdbfile", pdbout,
                  "--targetchain", data["chainid"] ]
        if len(data["model"]) > 0:
            args += [ "--preserve_model", buildAbsPath(data["model"]) ]

        if len(data["maskadditions"]) > 0:
            args += [ "--maskextras", buildAbsPath(data["maskadditions"]) ]

        if len(data["mtzdata"]) > 0:
            args += [ "--mtzfile", buildAbsPath(data["mtzdata"]) ]
            if data["phasecalc"] == 1:
                args += [ "--phase_prob", "yes" ]
            elif data["phasecalc"] == 2:
                args += [ "--phase_comb", "yes" ]

        result = subprocess.call(args)
        exit(result)


if __name__ == '__main__':

    app = QtGui.QApplication(sys.argv)

    dialog = NCSModelingDialog()

    if dialog.exec_():
        dialog.runJob()
