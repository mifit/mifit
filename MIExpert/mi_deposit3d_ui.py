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


class JobReportDialog(QtGui.QDialog):
    def __init__(self, parent=None):
        super(JobReportDialog, self).__init__(parent)

        config = {
            'workdir' : "",
            'mtzfile' : "",
            'pdbfile' : "",
            'libfile' : "none",
            'seqfile' : "none",
            'templatefile' : "none",
            'datalogfile' : "none",

            'cif_write' : False,
            'map_write' : False,
            'map_border' : 0,
            'text_report' : False,
            'hkl_write' : False,

            'html_report' : False,
            'report_title' : "",
            'image1' : "",
            'image2' : "",
            'image3' : "",

            'rootname' : ""
          }
        settings = QtCore.QSettings("MIFit", "MIExpert")
        appSettings = settings.value("deposit3d").toString()
        if not appSettings.isEmpty():
            config = pickle.loads(str(appSettings))

        uiFile = os.path.join(os.path.dirname(sys.argv[0]), "mi_deposit3d.ui")
        uic.loadUi(uiFile, self)

        self.workingDirLineEdit.setText(config['workdir'])
        self.modelLineEdit.setText(config['pdbfile'])
        self.dataLineEdit.setText(config['mtzfile'])


        self.fileRootnameCheckBox.setChecked(False)
        if len(config['rootname']) != 0:
            self.fileRootnameCheckBox.setChecked(True)
            self.fileRootnameLineEdit.setText(config['rootname'])

        self.dictionaryCheckBox.setChecked(False)
        if config['libfile'] != 'none':
            self.dictionaryCheckBox.setChecked(True)
            self.dictionaryLineEdit.setText(config['libfile'])

        self.sequenceCheckBox.setChecked(False)
        if config['seqfile'] != 'none':
            self.sequenceCheckBox.setChecked(True)
            self.sequenceLineEdit.setText(config['seqfile'])

        self.processingLogCheckBox.setChecked(False)
        if config['datalogfile'] != 'none':
            self.processingLogCheckBox.setChecked(True)
            self.processingLogLineEdit.setText(config['datalogfile'])

        self.annotationCheckBox.setChecked(False)
        if config['templatefile'] != 'none':
            self.annotationCheckBox.setChecked(True)
            self.annotationLineEdit.setText(config['templatefile'])

        self.mmCIFReportCheckBox.setChecked(config['cif_write'])

        self.ligandDensityMapCheckBox.setChecked(config['map_write'])
        self.borderSpinBox.setValue(config['map_border'])
        self.textReportCheckBox.setChecked(config['text_report'])
        self.mmCIFDataCheckBox.setChecked(config['hkl_write'])

        self.htmlReportGroupBox.setChecked(config['html_report'])

        self.reportTitleCheckBox.setChecked(False)
        if len(config['report_title']) != 0:
            self.reportTitleCheckBox.setChecked(True)
            self.reportTitleLineEdit.setText(config['report_title'])

        self.image1CheckBox.setChecked(False)
        if len(config['image1']) != 0:
            self.image1CheckBox.setChecked(True)
            self.image1LineEdit.setText(config['image1'])

        self.image2CheckBox.setChecked(False)
        if len(config['image2']) != 0:
            self.image2CheckBox.setChecked(True)
            self.image2LineEdit.setText(config['image2'])

        self.image3CheckBox.setChecked(False)
        if len(config['image3']) != 0:
            self.image3CheckBox.setChecked(True)
            self.image3LineEdit.setText(config['image3'])

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
        self.browseFile(self.dataLineEdit, "Data file (*.mtz)")

    @QtCore.pyqtSlot()
    def on_dictionaryPushButton_clicked(self):
        self.browseFile(self.dictionaryLineEdit, "Dictionary file (*.cif *.lib)")

    @QtCore.pyqtSlot()
    def on_sequencePushButton_clicked(self):
        self.browseFile(self.sequenceLineEdit, "Sequence file (*.faa *.seq *.fasta)")

    @QtCore.pyqtSlot()
    def on_processingLogPushButton_clicked(self):
        self.browseFile(self.processingLogLineEdit, "Processing log file (*.*)")

    @QtCore.pyqtSlot()
    def on_annotationPushButton_clicked(self):
        self.browseFile(self.annotationLineEdit, "Annotation file (*.*)")

    @QtCore.pyqtSlot()
    def on_image1PushButton_clicked(self):
        self.browseFile(self.image1LineEdit, "Image file (*.jpg *.gif *.tif *.jpeg *.png)")

    @QtCore.pyqtSlot()
    def on_image2PushButton_clicked(self):
        self.browseFile(self.image2LineEdit, "Image file (*.jpg *.gif *.tif *.jpeg *.png)")

    @QtCore.pyqtSlot()
    def on_image3PushButton_clicked(self):
        self.browseFile(self.image3LineEdit, "Image file (*.jpg *.gif *.tif *.jpeg *.png)")


    def runJob(self):

        data = {
            'workdir' : "",
            'mtzfile' : "",
            'pdbfile' : "",
            'libfile' : "none",
            'seqfile' : "none",
            'templatefile' : "none",
            'datalogfile' : "none",

            'cif_write' : False,
            'map_write' : False,
            'map_border' : 0,
            'text_report' : False,
            'hkl_write' : False,

            'html_report' : False,
            'report_title' : "",
            'image1' : "",
            'image2' : "",
            'image3' : "",

            'rootname' : ""
          }

        data["workdir"] = str(self.workingDirLineEdit.text())
        data["mtzfile"] = str(self.dataLineEdit.text())
        data["pdbfile"] = str(self.modelLineEdit.text())

        if self.dictionaryCheckBox.isChecked() and self.dictionaryLineEdit.text().size() != 0:
            data["libfile"] = str(self.dictionaryLineEdit.text())
        else:
            data["libfile"] = "none"

        if self.sequenceCheckBox.isChecked() and self.sequenceLineEdit.text().size() != 0:
            data["seqfile"] = str(self.sequenceLineEdit.text())
        else:
            data["seqfile"] = "none"

        if self.annotationCheckBox.isChecked() and self.annotationLineEdit.text().size() != 0:
            data["templatefile"] = str(self.annotationLineEdit.text())
        else:
            data["templatefile"] = "none"

        if self.processingLogCheckBox.isChecked() and self.processingLogLineEdit.text().size() != 0:
            data["datalogfile"] = str(self.processingLogLineEdit.text())
        else:
            data["datalogfile"] = "none"

        data["cif_write"] = self.mmCIFReportCheckBox.isChecked()
        data["map_write"] = self.ligandDensityMapCheckBox.isChecked()
        data["map_border"] = self.borderSpinBox.value()
        data["text_report"] = self.textReportCheckBox.isChecked()
        data["hkl_write"] = self.mmCIFDataCheckBox.isChecked()
        data["html_report"] = self.htmlReportGroupBox.isChecked()        

        data["report_title"] = ""
        if self.reportTitleCheckBox.isChecked():
            data["report_title"] = str(self.reportTitleLineEdit.text())

        data["image1"] = ""
        if self.image1CheckBox.isChecked():
            data["image1"] = str(self.image1LineEdit.text())

        data["image2"] = ""
        if self.image2CheckBox.isChecked():
            data["image2"] = str(self.image2LineEdit.text())

        data["image3"] = ""
        if self.image3CheckBox.isChecked():
            data["image3"] = str(self.image3LineEdit.text())

        data["rootname"] = ""
        if self.fileRootnameCheckBox.isChecked():
            data["rootname"] = str(self.fileRootnameLineEdit.text())

        settings = QtCore.QSettings("MIFit", "MIExpert")
        settings.setValue("deposit3d", pickle.dumps(data))

        mifit.setJobWorkDir(str(QtCore.QFileInfo(data['workdir']).absoluteFilePath()))
        miexpert = os.path.join(os.path.dirname(sys.argv[0]), "MIExpert.py")
        args = [ sys.executable, miexpert, "deposit3d" ]

        args += [ "--workdir", buildAbsPath(data["workdir"]),
                  "--mtzfile", buildAbsPath(data["mtzfile"]),
                  "--pdbfile", buildAbsPath(data["pdbfile"]),
                  "--molimagehome", buildAbsPath(mifit.mifitDir),
                  "--libfile", buildAbsPath(data["libfile"]),
                  "--seqfile", buildAbsPath(data["seqfile"]),
                  "--templatefile", buildAbsPath(data["templatefile"]),
                  "--datalogfile", buildAbsPath(data["datalogfile"]),
                ]
        args += [ "--cif_write", boolYesNo(data["cif_write"]) ]

        if data["map_write"]:
            args += [ "--map_write", "yes", "--map_border",  str(data["map_border"]) ]

        args += [ "--text_write",  boolYesNo(data["text_report"]),
                  "--html_write", boolYesNo(data["html_report"]),
                  "--hkl_write", boolYesNo(data["hkl_write"]) ]

        if data["html_report"]:
            if len(data["report_title"]) != 0:
               args += [ "--title", data["report_title"] ]
            if len(data["image1"]) != 0:
                args += [ "--image1", buildAbsPath(data["image1"]) ]
            if len(data["image2"]) != 0:
                args += [ "--image2", buildAbsPath(data["image2"]) ]
            if len(data["image3"]) != 0:
                args += [ "--image3", buildAbsPath(data["image3"]) ]
        if len(data["rootname"]) != 0:
            args += [ "--rootname", buildAbsPath(data["rootname"]) ]

        result = subprocess.call(args)
        exit(result)


    @QtCore.pyqtSlot()
    def validate(self):
        globalEnabled = True

        thisEnabled = (QtCore.QFileInfo(self.workingDirLineEdit.text()).exists() and
                          QtCore.QFileInfo(self.workingDirLineEdit.text()).isDir())
        globalEnabled = markEnabled(self.workingDirLineEdit, thisEnabled, globalEnabled)

        thisEnabled = QtCore.QFileInfo(self.modelLineEdit.text()).exists();
        globalEnabled = markEnabled(self.modelLineEdit, thisEnabled, globalEnabled);

        thisEnabled = QtCore.QFileInfo(self.dataLineEdit.text()).exists();
        globalEnabled = markEnabled(self.dataLineEdit, thisEnabled, globalEnabled);

        thisEnabled = (not self.dictionaryCheckBox.isChecked() or
                         QtCore.QFileInfo(self.dictionaryLineEdit.text()).exists());
        globalEnabled = markEnabled(self.dictionaryLineEdit, thisEnabled, globalEnabled);

        thisEnabled = (not self.sequenceCheckBox.isChecked() or
                         QtCore.QFileInfo(self.sequenceLineEdit.text()).exists());
        globalEnabled = markEnabled(self.sequenceLineEdit, thisEnabled, globalEnabled);

        thisEnabled = (not self.processingLogCheckBox.isChecked() or
                         QtCore.QFileInfo(self.processingLogLineEdit.text()).exists());
        globalEnabled = markEnabled(self.processingLogLineEdit, thisEnabled, globalEnabled);

        thisEnabled = (not self.annotationCheckBox.isChecked() or
                         QtCore.QFileInfo(self.annotationLineEdit.text()).exists());
        globalEnabled = markEnabled(self.annotationLineEdit, thisEnabled, globalEnabled);


        thisEnabled = (not self.image1CheckBox.isChecked() or
                         QtCore.QFileInfo(self.image1LineEdit.text()).exists());
        globalEnabled = markEnabled(self.image1LineEdit, thisEnabled, globalEnabled);

        thisEnabled = (not self.image2CheckBox.isChecked() or
                         QtCore.QFileInfo(self.image2LineEdit.text()).exists());
        globalEnabled = markEnabled(self.image2LineEdit, thisEnabled, globalEnabled);

        thisEnabled = (not self.image3CheckBox.isChecked() or
                         QtCore.QFileInfo(self.image3LineEdit.text()).exists());
        globalEnabled = markEnabled(self.image3LineEdit, thisEnabled, globalEnabled);


        thisEnabled = (not self.reportTitleCheckBox.isChecked() or
                         self.reportTitleLineEdit.text().size());
        globalEnabled = markEnabled(self.reportTitleLineEdit, thisEnabled, globalEnabled);

        thisEnabled = (not self.reportTitleCheckBox.isChecked() or
                         self.reportTitleLineEdit.text().size());
        globalEnabled = markEnabled(self.reportTitleLineEdit, thisEnabled, globalEnabled);

        thisEnabled = (not self.fileRootnameCheckBox.isChecked() or
                         self.fileRootnameLineEdit.text().size());
        globalEnabled = markEnabled(self.fileRootnameLineEdit, thisEnabled, globalEnabled);

        okButton = self.buttonBox.button(QtGui.QDialogButtonBox.Ok)
        if okButton:
            okButton.setEnabled(globalEnabled)


if __name__ == '__main__':

    app = QtGui.QApplication(sys.argv)

    dialog = JobReportDialog()

    if dialog.exec_():
        dialog.runJob()
