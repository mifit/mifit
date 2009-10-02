import sys, os
from PyQt4 import QtCore, QtGui, uic


class RefinementDialog(QtGui.QDialog):
    def __init__(self, parent=None):
        super(RefinementDialog, self).__init__(parent)
        self.mifitDir = QtCore.QString()
        self.shelxDir = QtCore.QString()

        uic.loadUi("mi_refine.ui", self)

        self.timer = QtCore.QTimer(self)
        self.connect(self.timer, QtCore.SIGNAL("timeout()"),
                     self, QtCore.SLOT("validate()"))
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
        self.browseFile(self.modelLineEdit, "PDB files (*.pdb *.ent)")

    @QtCore.pyqtSlot()
    def on_dataPushButton_clicked(self):
        self.browseFile(self.dataLineEdit, "Data files (*.mtz)")

    @QtCore.pyqtSlot()
    def on_tlsSpecificationPushButton_clicked(self):
        self.browseFile(self.tlsSpecificationLineEdit,"TLS files (*.tls)")

    @QtCore.pyqtSlot()
    def on_dictionaryPushButton_clicked(self):
        self.browseFile(self.dictionaryLineEdit, "Dictionary files (*.cif *.lib)")

    def markRequired(self, widget, required):
        widget.setAutoFillBackground(True)
        pal = self.palette()
        if required:
            pal.setColor(QtGui.QPalette.Normal, QtGui.QPalette.Base, QtGui.QColor(255,255,128))

        widget.setPalette(pal);

    @QtCore.pyqtSlot()
    def validate(self):
        somethingRequired = False

        thisRequired = not (QtCore.QFile.exists(self.workingDirLineEdit.text()) \
                       and QtCore.QFileInfo(self.workingDirLineEdit.text()).isDir())
        self.markRequired(self.workingDirLineEdit, thisRequired)
        somethingRequired = somethingRequired or thisRequired

        thisRequired = not QtCore.QFile.exists(self.modelLineEdit.text())
        self.markRequired(self.modelLineEdit, thisRequired);
        somethingRequired = somethingRequired or thisRequired

        thisRequired = not QtCore.QFile.exists(self.dataLineEdit.text())
        self.markRequired(self.dataLineEdit, thisRequired)
        somethingRequired = somethingRequired or thisRequired

        thisRequired = not (not self.tlsSpecificationCheckBox.isChecked() \
                       or QtCore.QFile.exists(self.tlsSpecificationLineEdit.text()))
        self.markRequired(self.tlsSpecificationLineEdit, thisRequired)
        somethingRequired = somethingRequired or thisRequired

        thisRequired = not (not self.dictionaryCheckBox.isChecked() \
                       or QtCore.QFile.exists(self.dictionaryLineEdit.text()))
        self.markRequired(self.dictionaryLineEdit, thisRequired)
        somethingRequired = somethingRequired or thisRequired

        okButton = self.buttonBox.button(QtGui.QDialogButtonBox.Ok)
        if okButton:
            okButton.setEnabled(not somethingRequired)

    def runJob(self):
        args = [ "python", "MIExpert.py", "refine" ]
        args += [ "--mifithome", str(QtCore.QFileInfo(self.mifitDir).absoluteFilePath()) ]
        workDir = str(QtCore.QFileInfo(self.workingDirLineEdit.text()).absoluteFilePath())
        args += [ "--workdir", workDir ]
        args += [ "--pdbfile", str(QtCore.QFileInfo(self.modelLineEdit.text()).absoluteFilePath()) ]
        args += [ "--mtzfile", str(QtCore.QFileInfo(self.dataLineEdit.text()).absoluteFilePath()) ]
        args += [ "--weight", str(QtCore.QString.number(self.weightSpinBox.value())) ]
        args += [ "--cycles", str(QtCore.QString.number(self.numCyclesSpinBox.value())) ]

        if not self.shelxDir.isEmpty():
           args += [ "--shelx_dir", str(QtCore.QFileInfo(self.shelxDir).absoluteFilePath()) ]

        if self.waterPickingCheckBox.isChecked():
            args += [ "--water_cycles", str(QtCore.QString.number(self.waterPickingCyclesSpinBox.value())) ]

        args += [ "--bref_type" ]
        if self.anisotropicRadioButton.isChecked():
            args += [ "anisotropic" ]
        else:
            args += [ "isotropic" ]

        args += [ "--engine" ]
        if self.rigidBodyRadioButton.isChecked():
            args += [ "rigid" ]
        elif self.shelxRadioButton.isChecked():
            args += [ "shelx" ]
        else:
            args += [ "refmac5" ]


        if self.maxResolutionCheckBox.isChecked():
            args += [ "--max_res", str(QtCore.QString.number(self.maxResolutionSpinBox.value())) ]

        if self.tlsSpecificationCheckBox.isChecked():
            args += [ "--tls_file", str(QtCore.QFileInfo(self.tlsSpecificationLineEdit.text()).absoluteFilePath()) ]

        if self.dictionaryCheckBox.isChecked():
            args += [ "--libfile", str(QtCore.QFileInfo(self.dictionaryLineEdit.text()).absoluteFilePath()) ]

        # TODO Communicate workdir back to MIFit

        # Execute python mi_refine.py args
        print args
        os.execvp("python", args)


if __name__ == '__main__':

    app = QtGui.QApplication(sys.argv)

    dialog = RefinementDialog()
    # TODO read options from MIFit
    dialog.mifitDir = QtCore.QString("mifitDir")
    dialog.shelxDir = QtCore.QString("shelxDir")
    if dialog.exec_():
        dialog.runJob()

