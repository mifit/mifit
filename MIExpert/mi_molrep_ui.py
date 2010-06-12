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

class MolRepDialog(QtGui.QDialog):
    def __init__(self, parent=None):
        super(MolRepDialog, self).__init__(parent)

        data = {
           "workdir" : '',
           "pdbfile" : '',
           "mtzfile" : '',
           "fixed_pdb" : '',

           "multi_search" : False,
           "match_pdbin" : False,
           "engine" : '',

           "spacegroup_no" : 0,
           "sg_search" : False,

           "copies" : 1,
        }
        settings = QtCore.QSettings("MIFit", "MIExpert")
        appSettings = settings.value("molrep").toString()
        if not appSettings.isEmpty():
            data = pickle.loads(str(appSettings))

        uiFile = os.path.join(os.path.dirname(sys.argv[0]), "mi_molrep.ui")
        uic.loadUi(uiFile, self)

        spaceGroupList = mifit.spacegroupList()
        if spaceGroupList:
            self.spaceGroupComboBox.clear()
            for i in spaceGroupList:
                self.spaceGroupComboBox.addItem(i)
            self.spaceGroupComboBox.setCurrentIndex(-1)

        self.workingDirLineEdit.setText(data["workdir"])
        self.modelLineEdit.setText(data["pdbfile"])
        self.dataLineEdit.setText(data["mtzfile"])

        self.fixedModelCheckBox.setChecked(False)
        if len(data["fixed_pdb"]) != 0:
            self.fixedModelCheckBox.setChecked(True)
            self.fixedModelLineEdit.setText(data["fixed_pdb"])

        self.spaceGroupRadioButton.setChecked(False)
        if data["spacegroup_no"] > 0:
            self.spaceGroupRadioButton.setChecked(True)
            self.spaceGroupComboBox.setCurrentIndex(data["spacegroup_no"]-1)

        self.searchMultipleCheckBox.setChecked(data["multi_search"])
        self.matchInputCheckBox.setChecked(data["match_pdbin"])

        if data["engine"] == "phaser":
            self.phaserRadioButton.setChecked(True)
        else:
            self.molRepRadioButton.setChecked(True)

        self.copiesSpinBox.setValue(data["copies"])

        self.timer = QtCore.QTimer(self)
        self.connect(self.timer, QtCore.SIGNAL("timeout()"),
                     self, QtCore.SLOT("validateTimeout()"))
        self.timer.start(100)

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
        self.browseFile(self.modelLineEdit, "PDB file (*.pdb *.res *.ent)")

    @QtCore.pyqtSlot()
    def on_fixedModelPushButton_clicked(self):
        self.browseFile(self.fixedModelLineEdit, "PDB file (*.pdb *.res *.ent)")

    @QtCore.pyqtSlot()
    def on_dataPushButton_clicked(self):
        self.browseFile(self.dataLineEdit, "Data file (*.mtz *.ref *.sca)")

    def runJob(self):

        data = {
           "workdir" : '',
           "pdbfile" : '',
           "mtzfile" : '',
           "fixed_pdb" : '',

           "multi_search" : False,
           "match_pdbin" : False,
           "engine" : '',

           "spacegroup_no" : 1,
           "sg_search" : False,

           "copies" : 1,
        }

        data["workdir"] = str(self.workingDirLineEdit.text())
        data["pdbfile"] = str(self.modelLineEdit.text())
        data["mtzfile"] = str(self.dataLineEdit.text())

        if self.fixedModelCheckBox.isChecked():
            data["fixed_pdb"] = str(self.fixedModelLineEdit.text())

        data["spacegroup_no"] = 0
        if self.spaceGroupRadioButton.isChecked():
            data["spacegroup_no"] = self.spaceGroupComboBox.currentIndex() + 1

        data["multi_search"] = self.searchMultipleCheckBox.isChecked()
        data["match_pdbin"] = self.matchInputCheckBox.isChecked()

        data["engine"] = "molrep"
        if self.phaserRadioButton.isChecked():
            data["engine"] = "phaser"

        data["copies"] = self.copiesSpinBox.value()

        settings = QtCore.QSettings("MIFit", "MIExpert")
        settings.setValue("molrep", pickle.dumps(data))

        miexpert = os.path.join(os.path.dirname(sys.argv[0]), "MIExpert.py")
        args = [ sys.executable, miexpert, "molrep" ]

        args += [ "--engine", data["engine"] ]
        if data["spacegroup_no"] > 0:
            args += [ "--spacegroup", str(data["spacegroup_no"]) ]
        elif data["sg_search"]:
          args += [ "--sg_search", "yes" ]

        args += [ "--pdbfile", buildAbsPath(data["pdbfile"]),
                  "--mtzfile", buildAbsPath(data["mtzfile"]),
                  "--workdir", buildAbsPath(data["workdir"]),
                  "--multi_search", boolYesNo(data["multi_search"]),
                  "--match_pdbin", boolYesNo(data["match_pdbin"]),
                  "--copies", str(data["copies"]) ]
        if len(data["fixed_pdb"]) != 0:
            args += [ "--fixed_pdb", buildAbsPath(data["fixed_pdb"]) ]

        result = subprocess.call(args)
        exit(result)

    @QtCore.pyqtSlot()
    def validateTimeout(self):
        globalEnabled = True

        thisEnabled = QtCore.QFileInfo(self.workingDirLineEdit.text()).exists() and QtCore.QFileInfo(self.workingDirLineEdit.text()).isDir()
        globalEnabled = markEnabled(self.workingDirLineEdit, thisEnabled, globalEnabled)

        thisEnabled = QtCore.QFileInfo(self.modelLineEdit.text()).exists()
        globalEnabled = markEnabled(self.modelLineEdit, thisEnabled, globalEnabled)

        thisEnabled = QtCore.QFileInfo(self.dataLineEdit.text()).exists()
        globalEnabled = markEnabled(self.dataLineEdit, thisEnabled, globalEnabled)

        thisEnabled = (not self.fixedModelCheckBox.isChecked() or QtCore.QFileInfo(self.fixedModelLineEdit.text()).exists())
        globalEnabled = markEnabled(self.fixedModelLineEdit, thisEnabled, globalEnabled)

        okButton = self.buttonBox.button(QtGui.QDialogButtonBox.Ok)
        if okButton:
            okButton.setEnabled(globalEnabled)


if __name__ == '__main__':

    app = QtGui.QApplication(sys.argv)

    dialog = MolRepDialog()

    if dialog.exec_():
        dialog.runJob()
