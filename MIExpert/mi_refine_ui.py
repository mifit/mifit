import sys, os, pickle, subprocess
from PyQt4 import QtCore, QtGui, uic
import mifit


class RefinementDialog(QtGui.QDialog):
    def __init__(self, parent=None):
        super(RefinementDialog, self).__init__(parent)

        config = {
            'workdir': '',
            'pdbfile': '',
            'mtzfile': '',
            'weight': 0.1,
            'cycles': 5,
            'water_cycles': False,
            'water_cycles_value': 0,
            'bref_type': 'isotropic',
            'engine': 'refmac5',
            'max_res': False,
            'max_res_value': 0,
            'tls_file': False,
            'tls_file_value': '',
            'libfile': False,
            'libfile_value': '',
            'shelx_dir': '',
        }
        settings = QtCore.QSettings("MIFit", "MIExpert")
        appSettings = settings.value("refine").toString()
        if not appSettings.isEmpty():
            config = pickle.loads(str(appSettings))

        uiFile = os.path.join(os.path.dirname(sys.argv[0]), "mi_refine.ui")
        uic.loadUi(uiFile, self)

        self.workingDirLineEdit.setText(config['workdir'])
        self.modelLineEdit.setText(config['pdbfile'])
        self.dataLineEdit.setText(config['mtzfile'])
        self.weightSpinBox.setValue(config['weight'])
        self.numCyclesSpinBox.setValue(config['cycles'])
        self.waterPickingCheckBox.setChecked(config['water_cycles'])
        self.waterPickingCyclesSpinBox.setValue(config['water_cycles_value'])

        if config['bref_type'] == "anisotropic":
            self.anisotropicRadioButton.setChecked(True)
        else:
            self.anisotropicRadioButton.setChecked(False)

        if config['engine'] == "rigid":
            self.rigidBodyRadioButton.setChecked(True)
        elif config['engine'] == "shelx":
            self.shelxRadioButton.setChecked(True)
        else:
            self.refmac5RadioButton.setChecked(True)

        self.maxResolutionCheckBox.setChecked(config['max_res'])
        self.maxResolutionSpinBox.setValue(config['max_res_value'])
        self.tlsSpecificationCheckBox.setChecked(config['tls_file'])
        self.tlsSpecificationLineEdit.setText(config['tls_file_value'])
        self.dictionaryCheckBox.setChecked(config['libfile'])
        self.dictionaryLineEdit.setText(config['libfile_value'])

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

        config = {}
        config['workdir'] = str(self.workingDirLineEdit.text())
        config['pdbfile'] = str(self.modelLineEdit.text())
        config['mtzfile'] = str(self.dataLineEdit.text())
        config['weight'] = self.weightSpinBox.value()
        config['cycles'] = self.numCyclesSpinBox.value()
        config['water_cycles'] = self.waterPickingCheckBox.isChecked()
        config['water_cycles_value'] = self.waterPickingCyclesSpinBox.value()
        if self.anisotropicRadioButton.isChecked():
            config['bref_type'] = "anisotropic"
        else:
            config['bref_type'] = "isotropic"
        if self.rigidBodyRadioButton.isChecked():
            config['engine'] = "rigid"
        elif self.shelxRadioButton.isChecked():
            config['engine'] = "shelx"
        else:
            config['engine'] = "refmac5"
        config['max_res'] = self.maxResolutionCheckBox.isChecked()
        config['max_res_value'] = self.maxResolutionSpinBox.value()
        config['tls_file'] = self.tlsSpecificationCheckBox.isChecked()
        config['tls_file_value'] = str(self.tlsSpecificationLineEdit.text())
        config['libfile'] = self.dictionaryCheckBox.isChecked()
        config['libfile_value'] = str(self.dictionaryLineEdit.text())

        settings = QtCore.QSettings("MIFit", "MIExpert")
        settings.setValue("refine", pickle.dumps(config))

        mifit.setJobWorkDir(str(QtCore.QFileInfo(config['workdir']).absoluteFilePath()))
        miexpert = os.path.join(os.path.dirname(sys.argv[0]), "MIExpert.py")
        args = [ sys.executable, miexpert, "refine" ]
        args += [ "--mifithome", str(QtCore.QFileInfo(mifit.mifitDir).absoluteFilePath()) ]
        args += [ "--workdir", str(QtCore.QFileInfo(config['workdir']).absoluteFilePath()) ]
        args += [ "--pdbfile", str(QtCore.QFileInfo(config['pdbfile']).absoluteFilePath()) ]
        args += [ "--mtzfile", str(QtCore.QFileInfo(config['mtzfile']).absoluteFilePath()) ]
        args += [ "--weight", str(config['weight']) ]
        args += [ "--cycles", str(config['cycles']) ]

        if mifit.shelxDir != None:
           args += [ "--shelx_dir", str(QtCore.QFileInfo(mifit.shelxDir).absoluteFilePath()) ]

        if config['water_cycles']:
            args += [ "--water_cycles", str(config['water_cycles_value']) ]

        args += [ "--bref_type", config['bref_type'] ]
        args += [ "--engine", config['engine'] ]

        if config['max_res']:
            args += [ "--max_res", str(config['max_res_value']) ]

        if config['tls_file']:
            args += [ "--tls_file", str(QtCore.QFileInfo(config['tls_file_value']).absoluteFilePath()) ]

        if config['libfile']:
            args += [ "--libfile", str(QtCore.QFileInfo(config['libfile_value']).absoluteFilePath()) ]

        result = subprocess.call(args)
        exit(result)


if __name__ == '__main__':

    app = QtGui.QApplication(sys.argv)

    dialog = RefinementDialog()

    if dialog.exec_():
        dialog.runJob()
