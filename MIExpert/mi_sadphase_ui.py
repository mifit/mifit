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
            "ssnumber" : -1,
            "separation" : 0,
            "solventfraction" : 0,
            "siterefinemethod" : '',
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
            sg = spaceGroupList.split(',')
            self.spaceGroupComboBox.clear()
            for i in sg:
                self.spaceGroupComboBox.addItem(i)
            self.spaceGroupComboBox.setCurrentIndex(0)

        self.workingDirectoryLineEdit.setText(data["workdir"])
        self.intensityDataLineEdit.setText(data["saddatafile"])

        if data["sitefile"] != "none":
            self.sitesFromFileRadioButton.setChecked(True)
            self.sitesFromFileLineEdit.setText(data["sitefile"])

        data["scatterer"] = str(self.scattererTypeLineEdit.text())

        self.numSitesSpinBox.setValue(data["sitenumber"])
        self.numDisulfidesSpinBox.setValue(data["ssnumber"])
        self.minScattererSeparationSpinBox.setValue(data["separation"])
        self.solventFractionSpinBox.setValue(data["solventfraction"])
        if data["siterefinemethod"] == "bp3":
            self.bp3RadioButton.isChecked(True)
        else:
            self.mlphareRadioButton.setChecked(True)
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
        self.browseFile(self.intensityDataLineEdit, "All files (*)")

    @QtCore.pyqtSlot()
    def on_sitesFromFilePushButton_clicked(self):
        self.browseFile(self.sitesFromFileLineEdit, "All files (*)")


    def runJob(self):

        data = {
            "workdir" : '',
            "saddatafile" : '',

            "sitefile" : 'none',
            "scatterer" : '',

            "sitenumber" : -1,
            "ssnumber" : -1,
            "separation" : 0,
            "solventfraction" : 0,
            "siterefinemethod" : '',
            "bothhands" : False,

            "spacegroup_no" : 1,
            "change_spacegroup" : False
        }
        data["workdir"] = str(self.workingDirectoryLineEdit.text())
        data["saddatafile"] = str(self.intensityDataLineEdit.text())

        data["sitefile"] = "none"
        if self.sitesFromFileRadioButton.isChecked():
          data["sitefile"] = str(self.sitesFromFileLineEdit.text())

        data["scatterer"] = str(self.scattererTypeLineEdit.text())

        data["sitenumber"] = self.numSitesSpinBox.value()
        data["ssnumber"] = self.numDisulfidesSpinBox.value()
        data["separation"] = self.minScattererSeparationSpinBox.value()
        data["solventfraction"] = self.solventFractionSpinBox.value()
        if self.bp3RadioButton.isChecked():
          data["siterefinemethod"] = "bp3"
        else:
          data["siterefinemethod"] = "mlphare"
        data["bothhands"] = self.phaseBothSiteEnantiomorphsCheckBox.isChecked()

        data["change_spacegroup"] = self.changeSpaceGroupCheckBox.isChecked()
        data["spacegroup_no"] = self.spaceGroupComboBox.currentIndex() + 1

        settings = QtCore.QSettings("MIFit", "MIExpert")
        settings.setValue("sadphase", pickle.dumps(data))

        mifit.setJobWorkDir(str(QtCore.QFileInfo(config['workdir']).absoluteFilePath()))
        miexpert = os.path.join(os.path.dirname(sys.argv[0]), "MIExpert.py")
        args = [ sys.executable, miexpert, "sadphase" ]

        args += [ "--molimagehome", buildAbsPath(mifit.mifitDir) ]
        args += [ "--workdir", buildAbsPath(data["workdir"]) ]
        args += [ "--saddatafile", buildAbsPath(data["saddatafile"]) ]
        args += [ "--sitefile", buildAbsPath(data["sitefile"]) ]
        args += [ "--scatterer", data["scatterer"] ]
        args += [ "--sitenumber", str(data["sitenumber"]) ]
        args += [ "--ssnumber", str(data["ssnumber"]) ]
        args += [ "--separation", str(data["separation"]) ]
        args += [ "--solventfraction", str(data["solventfraction"]) ]
        args += [ "--siterefinemethod", data["siterefinemethod"] ]
        if mifit.shelxDir != None:
            args += [ "--shelx_dir", buildAbsPath(mifit.shelxDir) ]
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

        thisEnabled = (not self.sitesFromFileRadioButton.isChecked() or QtCore.QFileInfo(self.sitesFromFileLineEdit.text()).exists())
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
