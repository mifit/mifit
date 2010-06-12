import sys, os, pickle, subprocess
from PyQt4 import QtCore, QtGui, uic
import mifit

def buildAbsPath(path):
    return str(QtCore.QFileInfo(path).absoluteFilePath())

def markEnabled(w, thisEnabled, globalEnabled):
    w.setAutoFillBackground(True)
    if not thisEnabled:
        pal = w.palette()
        pal.setColor(QtGui.QPalette.Normal, QtGui.QPalette.Base, QtGui.QColor(255, 255, 128))
        w.setPalette(pal)
    else:
        w.setPalette(QtGui.QPalette())

    return globalEnabled and thisEnabled


class IntegrateDialog(QtGui.QDialog):
    def __init__(self, parent=None):
        super(IntegrateDialog, self).__init__(parent)

        self.data = {
            "template_image" : "",
            "detector_constants" : "",

            "first_image" : -1,
            "last_image" : -1,

            "integrate_resolution" : "",

            "spacegroup_no" : 1,
          }
        settings = QtCore.QSettings("MIFit", "MIExpert")
        appSettings = settings.value("integrate").toString()
        if not appSettings.isEmpty():
            self.data = pickle.loads(str(appSettings))

        uiFile = os.path.join(os.path.dirname(sys.argv[0]), "mi_integrate.ui")
        uic.loadUi(uiFile, self)

        self.spaceGroupComboBox.clear()
        sgList = mifit.spacegroupList()
        if sgList:
            for i in sgList:
                self.spaceGroupComboBox.addItem(i)
            self.spaceGroupComboBox.setCurrentIndex(self.data["spacegroup_no"] - 1)


        self.intensityData.setText(self.data["template_image"])

        if len(self.data["detector_constants"]) != 0:
            self.resetHardwareParametersCheckBox.setChecked(True)
            self.hardwareParametersLineEdit.setText(self.data["detector_constants"])

        if self.data["first_image"] != -1:
            self.firstImageCheckBox.setChecked(True)
            self.firstImageSpinBox.setValue(self.data["first_image"])
        if self.data["last_image"] != -1:
            self.lastImageCheckBox.setChecked(True)
            self.lastImageSpinBox.setValue(self.data["last_image"])

        if len(self.data["integrate_resolution"]) != 0:
            self.resolutionRangeCheckBox.setChecked(True)
            res = QtCore.QString(self.data["integrate_resolution"]).split(" ")
            self.maxResSpinBox.setValue(float(str(res[0])))
            self.minResSpinBox.setValue(float(str(res[1])))

        self.spaceGroupComboBox.setCurrentIndex(self.data["spacegroup_no"] - 1)

        self.timer = QtCore.QTimer(self)
        self.connect(self.timer, QtCore.SIGNAL("timeout()"),
                     self, QtCore.SLOT("validateTimeout()"))
        self.timer.start(100)


    @QtCore.pyqtSlot()
    def on_intensityDataBrowse_clicked(self):
        dir = QtGui.QFileDialog.getExistingDirectory(None, "", self.intensityData.text())
        if not dir.isEmpty():
            self.intensityData.setText(dir)

    @QtCore.pyqtSlot()
    def on_hardwareParametersPushButton_clicked(self):
        f = QtGui.QFileDialog.getOpenFileName(self, "Select a file", self.hardwareParametersLineEdit.text(), "All Files (*.* *)")
        if not f.isEmpty():
            self.hardwareParametersLineEdit.setText(f)


    def runJob(self):

        self.data["template_image"] = str(self.intensityData.text())

        self.data["detector_constants"] = ""
        if self.resetHardwareParametersCheckBox.isChecked():
          self.data["detector_constants"] = str(self.hardwareParametersLineEdit.text());

        self.data["first_image"] = -1
        if self.firstImageCheckBox.isChecked():
          self.data["first_image"] = self.firstImageSpinBox.value()
        self.data["last_image"] = -1
        if self.lastImageCheckBox.isChecked():
          self.data["last_image"] = self.lastImageSpinBox.value()

        self.data["integrate_resolution"] = ""
        if self.resolutionRangeCheckBox.isChecked():
            self.data["integrate_resolution"] = str(self.maxResSpinBox.text()) + " " + str(self.minResSpinBox.text())


        self.data["spacegroup_no"] = self.spaceGroupComboBox.currentIndex() + 1

        settings = QtCore.QSettings("MIFit", "MIExpert")
        settings.setValue("integrate", pickle.dumps(self.data))

        miexpert = os.path.join(os.path.dirname(sys.argv[0]), "MIExpert.py")
        args = [ sys.executable, miexpert, "integrate" ]

        args += [ "--template_image", self.data["template_image"] ]
        if len(self.data["detector_constants"]) != 0:
            args += [ "--detector_constants", buildAbsPath(self.data["detector_constants"]) ]
        args += [ "--spacegroup", str(self.data["spacegroup_no"]) ]
        if self.data["first_image"] >= 0:
            args += [ "--first_image", str(dself.ata["first_image"]) ]
        if self.data["last_image"] >= 0:
            args += [ "--last_image", str(self.data["last_image"]) ]

        if len(self.data["integrate_resolution"]) != 0:
            args += [ "--integrate_resolution", self.data["integrate_resolution"] ]


        result = subprocess.call(args)
        exit(result)

    @QtCore.pyqtSlot()
    def validateTimeout(self):
        globalEnabled = True

        thisEnabled = not self.intensityData.text().isEmpty()
        globalEnabled = markEnabled(self.intensityData, thisEnabled, globalEnabled)

        thisEnabled = (not self.resetHardwareParametersCheckBox.isChecked() or
                         QtCore.QFileInfo(self.hardwareParametersLineEdit.text()).exists())
        globalEnabled = markEnabled(self.hardwareParametersLineEdit, thisEnabled, globalEnabled)

        okButton = self.buttonBox.button(QtGui.QDialogButtonBox.Ok)
        if okButton:
            okButton.setEnabled(globalEnabled)


if __name__ == '__main__':

    app = QtGui.QApplication(sys.argv)

    dialog = IntegrateDialog()

    if dialog.exec_():
        dialog.runJob()
