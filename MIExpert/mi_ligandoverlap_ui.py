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


class LigandOverlapDialog(QtGui.QDialog):
    def __init__(self, parent=None):
        super(LigandOverlapDialog, self).__init__(parent)

        data = {
            "workdir" : '',
            "pdbdir" : '',
            "targetpdb" : '',
            "x" : 0,
            "y" : 0,
            "z" : 0,
        }
        settings = QtCore.QSettings("MIFit", "MIExpert")
        appSettings = settings.value("ligandoverlap").toString()
        if not appSettings.isEmpty():
            data = pickle.loads(str(appSettings))

        uiFile = os.path.join(os.path.dirname(sys.argv[0]), "mi_ligandoverlap.ui")
        uic.loadUi(uiFile, self)

        self.workingDirLineEdit.setText(data["workdir"])
        self.structureLineEdit.setText(data["pdbdir"])
        self.targetFileLineEdit.setText(data["targetpdb"])
        self.xSpinBox.setValue(data["x"])
        self.ySpinBox.setValue(data["y"])
        self.zSpinBox.setValue(data["z"])

        self.timer = QtCore.QTimer(self)
        self.connect(self.timer, QtCore.SIGNAL("timeout()"),
                     self, QtCore.SLOT("validateTimeout()"))
        self.timer.start(100)

    @QtCore.pyqtSlot()
    def on_workingDirPushButton_clicked(self):
        dir = QtGui.QFileDialog.getExistingDirectory(None, "", self.workingDirLineEdit.text())
        if not dir.isEmpty():
            self.workingDirLineEdit.setText(dir)

    @QtCore.pyqtSlot()
    def on_structurePushButton_clicked(self):
        dir = QtGui.QFileDialog.getExistingDirectory(None, "", self.structureLineEdit.text())
        if not dir.isEmpty():
            self.structureLineEdit.setText(dir)

    @QtCore.pyqtSlot()
    def on_targetFilePushButton_clicked(self):
        f = QtGui.QFileDialog.getOpenFileName(None, "Select a file", self.targetFileLineEdit.text(), "PDB file (*.pdb *.res *.ent *.*)")
        if not f.isEmpty():
            self.targetFileLineEdit.setText(f)

    def runJob(self):

        data = {
            "workdir" : '',
            "pdbdir" : '',
            "targetpdb" : '',
            "x" : 0,
            "y" : 0,
            "z" : 0,
        }
        data["workdir"] = str(self.workingDirLineEdit.text())
        data["pdbdir"] = str(self.structureLineEdit.text())
        data["targetpdb"] = str(self.targetFileLineEdit.text())
        data["x"] = self.xSpinBox.value()
        data["y"] = self.ySpinBox.value()
        data["z"] = self.zSpinBox.value()

        settings = QtCore.QSettings("MIFit", "MIExpert")
        settings.setValue("ligandoverlap", pickle.dumps(data))

        mifit.setJobWorkDir(str(QtCore.QFileInfo(data['workdir']).absoluteFilePath()))
        miexpert = os.path.join(os.path.dirname(sys.argv[0]), "MIExpert.py")
        args = [ sys.executable, miexpert, "ligandoverlap" ]

        targetsite_parameters = str(QtCore.QString.number(data["x"], 'f', 3)) + ' ' + \
                                str(QtCore.QString.number(data["y"], 'f', 3)) + ' ' + \
                                str(QtCore.QString.number(data["z"], 'f', 3))

        args += [ "--workdir", buildAbsPath(data["workdir"]),
                  "--pdbdir", buildAbsPath(data["pdbdir"]),
                  "--targetpdb", buildAbsPath(data["targetpdb"]),
                  "--targetsite", targetsite_parameters]

        result = subprocess.call(args)
        exit(result)

    @QtCore.pyqtSlot()
    def validateTimeout(self):
        globalEnabled = True

        thisEnabled = QtCore.QFileInfo(self.workingDirLineEdit.text()).exists() and QtCore.QFileInfo(self.workingDirLineEdit.text()).isDir()
        globalEnabled = markEnabled(self.workingDirLineEdit, thisEnabled, globalEnabled)

        thisEnabled = QtCore.QFileInfo(self.structureLineEdit.text()).exists() and QtCore.QFileInfo(self.structureLineEdit.text()).isDir()
        globalEnabled = markEnabled(self.structureLineEdit, thisEnabled, globalEnabled)

        thisEnabled = QtCore.QFileInfo(self.targetFileLineEdit.text()).exists()
        globalEnabled = markEnabled(self.targetFileLineEdit, thisEnabled, globalEnabled)

        okButton = self.buttonBox.button(QtGui.QDialogButtonBox.Ok)
        if okButton:
            okButton.setEnabled(globalEnabled)


if __name__ == '__main__':

    app = QtGui.QApplication(sys.argv)

    dialog = LigandOverlapDialog()

    if dialog.exec_():
        dialog.runJob()
