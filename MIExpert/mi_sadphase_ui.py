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

class SadPhasingDialog(QtGui.QDialog):
    def __init__(self, parent=None):
        super(SadPhasingDialog, self).__init__(parent)

        data = {
            "workdir" : '',
            "saddatafile" : '',
            "sitefile" : 'none',
            "scatterer" : '',
            "sitenumber" : -1,
            "solventfraction" : 0,
            "bothhands" : False,
            "spacegroup_no" : 1,
            "change_spacegroup" : False
        }
        settings = QtCore.QSettings("MIFit", "MIExpert")
        appSettings = settings.value("sadphase").toString()
        if not appSettings.isEmpty():
            data = pickle.loads(str(appSettings))

        uiFile = os.path.join(os.path.dirname(sys.argv[0]), "mi_sadphase.ui")
        uic.loadUi(uiFile, self)

        spaceGroupList = mifit.spacegroupList()
        if spaceGroupList:
            self.spaceGroupComboBox.clear()
            for i in spaceGroupList:
                self.spaceGroupComboBox.addItem(i)
            self.spaceGroupComboBox.setCurrentIndex(0)

        self.workingDirectoryLineEdit.setText(data["workdir"])
        self.intensityDataLineEdit.setText(data["saddatafile"])
        self.sitesFromFileLineEdit.setText(data["sitefile"])        

        data["scatterer"] = str(self.scattererTypeLineEdit.text())

        self.numSitesSpinBox.setValue(data["sitenumber"])
        self.solventFractionSpinBox.setValue(data["solventfraction"])
        self.phaseBothSiteEnantiomorphsCheckBox.setChecked(data["bothhands"])
        self.changeSpaceGroupCheckBox.setChecked(data["change_spacegroup"])
        self.spaceGroupComboBox.setCurrentIndex(data["spacegroup_no"]-1)

        self.timer = QtCore.QTimer(self)
        self.connect(self.timer, QtCore.SIGNAL("timeout()"),
                     self, QtCore.SLOT("validateTimeout()"))
        self.timer.start(100)


    @QtCore.pyqtSlot()
    def on_workingDirectoryPushButton_clicked(self):
        dir = QtGui.QFileDialog.getExistingDirectory(None, "", self.workingDirectoryLineEdit.text())
        if not dir.isEmpty():
            self.workingDirectoryLineEdit.setText(dir)

    def browseFile(self, lineEdit, filter):
        f = QtGui.QFileDialog.getOpenFileName(None, "Select a file", lineEdit.text(), filter)
        if not f.isEmpty():
            lineEdit.setText(f)

    @QtCore.pyqtSlot()
    def on_intensityDataPushButton_clicked(self):
        self.browseFile(self.intensityDataLineEdit, "Data files (*.ref *.sca *.mtz)")

    @QtCore.pyqtSlot()
    def on_sitesFromFilePushButton_clicked(self):
        self.browseFile(self.sitesFromFileLineEdit, "PDB files, SHELXD files (*.pdb *.res)")


    def runJob(self):

        data = {
            "workdir" : '',
            "saddatafile" : '',
            "sitefile" : 'none',
            "scatterer" : '',
            "sitenumber" : -1,
            "solventfraction" : 0,
            "bothhands" : False,
            "spacegroup_no" : 1,
            "change_spacegroup" : False
        }
        data["workdir"] = str(self.workingDirectoryLineEdit.text())
        data["saddatafile"] = str(self.intensityDataLineEdit.text())
        data["sitefile"] = str(self.sitesFromFileLineEdit.text())
        
        data["scatterer"] = str(self.scattererTypeLineEdit.text())
        data["sitenumber"] = self.numSitesSpinBox.value()
        data["solventfraction"] = self.solventFractionSpinBox.value()
        data["bothhands"] = self.phaseBothSiteEnantiomorphsCheckBox.isChecked()
        data["change_spacegroup"] = self.changeSpaceGroupCheckBox.isChecked()
        data["spacegroup_no"] = self.spaceGroupComboBox.currentIndex() + 1

        settings = QtCore.QSettings("MIFit", "MIExpert")
        settings.setValue("sadphase", pickle.dumps(data))

        mifit.setJobWorkDir(str(QtCore.QFileInfo(data['workdir']).absoluteFilePath()))
        miexpert = os.path.join(os.path.dirname(sys.argv[0]), "MIExpert.py")
        args = [ sys.executable, miexpert, "sadphase" ]

        args += [ "--workdir", buildAbsPath(data["workdir"]) ]
        args += [ "--saddatafile", buildAbsPath(data["saddatafile"]) ]
        args += [ "--sitefile", buildAbsPath(data["sitefile"]) ]
        args += [ "--scatterer", data["scatterer"] ]
        args += [ "--sitenumber", str(data["sitenumber"]) ]
        args += [ "--solventfraction", str(data["solventfraction"]) ]
        args += [ "--bothhands", boolYesNo(data["bothhands"]) ]
        if data["change_spacegroup"]:
            args += [ "--spacegroup_no", str(data["spacegroup_no"]) ]

        result = subprocess.call(args)
        exit(result)


    @QtCore.pyqtSlot()
    def validateTimeout(self):
        globalEnabled = True

        thisEnabled = QtCore.QFileInfo(self.workingDirectoryLineEdit.text()).exists() and QtCore.QFileInfo(self.workingDirectoryLineEdit.text()).isDir()
        globalEnabled = markEnabled(self.workingDirectoryLineEdit, thisEnabled, globalEnabled);

        thisEnabled = QtCore.QFileInfo(self.intensityDataLineEdit.text()).exists();
        globalEnabled = markEnabled(self.intensityDataLineEdit, thisEnabled, globalEnabled);

        thisEnabled = QtCore.QFileInfo(self.sitesFromFileLineEdit.text()).exists();
        globalEnabled = markEnabled(self.sitesFromFileLineEdit, thisEnabled, globalEnabled);

        thisEnabled = not self.scattererTypeLineEdit.text().isEmpty()
        globalEnabled = markEnabled(self.scattererTypeLineEdit, thisEnabled, globalEnabled);

        okButton = self.buttonBox.button(QtGui.QDialogButtonBox.Ok)
        if okButton:
            okButton.setEnabled(globalEnabled)


if __name__ == '__main__':

    app = QtGui.QApplication(sys.argv)

    dialog = SadPhasingDialog()

    if dialog.exec_():
        dialog.runJob()
