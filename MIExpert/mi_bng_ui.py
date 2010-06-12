import sys, os, pickle, subprocess
from PyQt4 import QtCore, QtGui, QtNetwork, uic
import mifit


class BindNGrindDialog(QtGui.QDialog):
    def __init__(self, spaceGroupList=None, parent=None):
        super(BindNGrindDialog, self).__init__(parent)

        config = {
            'hklin': [],
            'pdbin': '',
            'process_engine': 'none',
            'detector_constants': '',
            'spacegroup_no': 0,
            'reference_mtz': 'none',
            'multi_search': False,
            'libfile': 'none',
            'viewpoint': False,
            'viewpoint_method': '',
            'viewpoint_coordinates': [],
            'viewpoint_file': '',
            'bngsummary': '',           
            'place_ligand': False,
            'ligand_name': '',
            'chem_name': ''
        }
        settings = QtCore.QSettings("MIFit", "MIExpert")
        appSettings = settings.value("bng").toString()
        if not appSettings.isEmpty():
            config = pickle.loads(str(appSettings))

        uiFile = os.path.join(os.path.dirname(sys.argv[0]), "mi_bng.ui")
        uic.loadUi(uiFile, self)

        if spaceGroupList:
            self.spaceGroupComboBox.addItem('Unchanged')
            for i in spaceGroupList:
                self.spaceGroupComboBox.addItem(i)
            self.spaceGroupComboBox.setCurrentIndex(0)

        for item in config['hklin']:
            self.intensityDataListWidget.addItem(item)

        self.modelPDBLineEdit.setText(config['pdbin'])

        if config['process_engine'] == 'none':
            self.preProcessedRadioButton.setChecked(True)
        elif config['process_engine'] == 'dstartrek':
            self.dtrekRadioButton.setChecked(True)
        else:
            self.mosflmRadioButton.setChecked(True)
        
        if len(config['detector_constants']) > 0:
            self.detectorConstantsCheckBox.setChecked(True)
            self.detectorConstantsLineEdit.setText(config['detector_constants'])

        if len(config['reference_mtz']) > 0:
            self.referenceDataCheckBox.setChecked(True)
            self.referenceDataLineEdit.setText(config['detector_constants'])

        self.multipleModelsCheckBox.setChecked(config['multi_search'])

        if config['libfile'] == 'none':
            self.dictionaryCIFCheckBox.setChecked(True)
            self.dictionaryCIFLineEdit.setText(config['detector_constants'])

        if len(config['bngsummary']) > 0:
            self.htmlSummaryCheckBox.setChecked(True)
            self.htmlSummaryLineEdit.setText(config['bngsummary'])

        if config['viewpoint']:
            self.viewpointGroupBox.setChecked(True)
            if config['viewpoint_method'] == 'coordinates':
                self.coordinatesRadioButton.setChecked(True)
                self.xSpinBox.setText(config['viewpoint_coordinates'][0])
                self.ySpinBox.setText(config['viewpoint_coordinates'][1])
                self.zSpinBox.setText(config['viewpoint_coordinates'][2])
            elif config['viewpoint_method'] == 'mlw':
                self.sessionFileRadioButton.setChecked(True)
                self.sessionFileLineEdit.setText(config['viewpoint_file'])
            elif config['viewpoint_method'] == 'pdbview':
                self.fromPDBMarkerFileRadioButton.setChecked(True)
                self.fromPDBMarkerFileLineEdit.setText(config['viewpoint_file'])

        if config['place_ligand']:
            self.placeLigandPDBGroupBox.setChecked(True)

            if self.chemFileRadioButton.setChecked(True):
                self.chemFileLineEdit.setText(config['chem_name'])
            else:
                self.ligandFileRadioButton.setChecked(True)
                self.ligandFileLineEdit.setText(config['ligand_name'])

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
    def on_modelPDBPushButton_clicked(self):
        self.browseFile(self.modelPDBLineEdit, "PDB files (*.pdb *.ent)")

    @QtCore.pyqtSlot()
    def on_detectorConstantsPushButton_clicked(self):
        self.browseFile(self.detectorConstantsLineEdit, "Constants file (*.txt)")

    @QtCore.pyqtSlot()
    def on_referenceDataPushButton_clicked(self):
        self.browseFile(self.referenceDataLineEdit, "Data file (*.mtz)")

    @QtCore.pyqtSlot()
    def on_htmlSummaryPushButton_clicked(self):
        dir = QtGui.QFileDialog.getExistingDirectory(None, "", self.htmlSummaryLineEdit.text())
        if not dir.isEmpty():
            self.htmlSummaryLineEdit.setText(dir)

    @QtCore.pyqtSlot()
    def on_dictionaryCIFPushButton_clicked(self):
        self.browseFile(self.dictionaryCIFLineEdit, "Dictionary file (*.cif *.lib)")

    @QtCore.pyqtSlot()
    def on_sessionFilePushButton_clicked(self):
        self.browseFile(self.sessionFileLineEdit, "Session file (*.mlw)")

    @QtCore.pyqtSlot()
    def on_fromPDBMarkerFilePushButton_clicked(self):
        self.browseFile(self.fromPDBMarkerFileLineEdit, "PDB Files (*.pdb *.ent)")

    @QtCore.pyqtSlot()
    def on_ligandFilePushButton_clicked(self):
        self.browseFile(self.ligandFileLineEdit, "PDB Files (*.pdb *.ent)")

    @QtCore.pyqtSlot()
    def on_chemFilePushButton_clicked(self):
        self.browseFile(self.chemFileLineEdit, "Chemistry Files (*.sdf *.sd *.mol *.cif)")       

    @QtCore.pyqtSlot()
    def on_addPushButton_clicked(self):
        f = QtGui.QFileDialog.getOpenFileName(None, "Select a file", "", "Data files (*.ref *.sca *.mtz *.img *.osc);;"
            + "CCP4 files (*.mtz);;"
            + "Reflection files (*.ref *.sca);;"
            + "Image files (*.img *.osc);;"
            + "All files (*.*)")
        if not f.isEmpty():
            self.intensityDataListWidget.addItem(f)


    @QtCore.pyqtSlot()
    def on_removePushButton_clicked(self):
        isel = self.intensityDataListWidget.currentRow()
        if isel >= 0:
          item = self.intensityDataListWidget.takeItem(isel)
          item = None

    def markRequired(self, widget, required):
        widget.setAutoFillBackground(True)
        pal = self.palette()
        if required:
            pal.setColor(QtGui.QPalette.Normal, QtGui.QPalette.Base, QtGui.QColor(255,255,128))

        widget.setPalette(pal);

    @QtCore.pyqtSlot()
    def validate(self):
        somethingRequired = False

        thisRequired = self.intensityDataListWidget.count() <= 0
        self.markRequired(self.intensityDataListWidget, thisRequired)
        somethingRequired = somethingRequired or thisRequired

        thisRequired = not QtCore.QFile.exists(self.modelPDBLineEdit.text())
        self.markRequired(self.modelPDBLineEdit, thisRequired)
        somethingRequired = somethingRequired or thisRequired

        thisRequired = self.detectorConstantsCheckBox.isChecked() and not QtCore.QFile.exists(self.detectorConstantsLineEdit.text())
        self.markRequired(self.detectorConstantsLineEdit, thisRequired)
        somethingRequired = somethingRequired or thisRequired

        thisRequired = self.referenceDataCheckBox.isChecked() and not QtCore.QFile.exists(self.referenceDataLineEdit.text())
        self.markRequired(self.referenceDataLineEdit, thisRequired)
        somethingRequired = somethingRequired or thisRequired

        thisRequired = self.htmlSummaryCheckBox.isChecked() and not QtCore.QFile.exists(self.htmlSummaryLineEdit.text())
        self.markRequired(self.htmlSummaryLineEdit, thisRequired)
        somethingRequired = somethingRequired or thisRequired

        thisRequired = self.dictionaryCIFCheckBox.isChecked() and not QtCore.QFile.exists(self.dictionaryCIFLineEdit.text())
        self.markRequired(self.dictionaryCIFLineEdit, thisRequired)
        somethingRequired = somethingRequired or thisRequired

        if self.viewpointGroupBox.isChecked():
            thisRequired = self.sessionFileRadioButton.isChecked() and not QtCore.QFile.exists(self.sessionFileLineEdit.text())
            self.markRequired(self.sessionFileLineEdit, thisRequired)
            somethingRequired = somethingRequired or thisRequired

            thisRequired = self.fromPDBMarkerFileRadioButton.isChecked() and not QtCore.QFile.exists(self.fromPDBMarkerFileLineEdit.text())
            self.markRequired(self.fromPDBMarkerFileLineEdit, thisRequired)
            somethingRequired = somethingRequired or thisRequired

            if self.placeLigandPDBGroupBox.isChecked():
                thisRequired = self.ligandFileRadioButton.isChecked() and not QtCore.QFile.exists(self.ligandFileLineEdit.text())
                self.markRequired(self.ligandFileLineEdit, thisRequired)
                somethingRequired = somethingRequired or thisRequired

        okButton = self.buttonBox.button(QtGui.QDialogButtonBox.Ok)
        if okButton:
            okButton.setEnabled(not somethingRequired)


    def runJob(self):

        config = {
            'hklin': [],
            'pdbin': '',
            'process_engine': 'none',
            'detector_constants': '',
            'spacegroup_no': 0,
            'reference_mtz': 'none',
            'multi_search': False,
            'libfile': 'none',
            'viewpoint': False,
            'viewpoint_method': '',
            'viewpoint_coordinates': [],
            'viewpoint_file': '',
            'bngsummary': '',          
            'place_ligand': False,
            'ligand_name': '',
            'chem_name':''
        }

        config['hklin'] = []
        for i in range(0, self.intensityDataListWidget.count()):
            config['hklin'] += [ str(self.intensityDataListWidget.item(i).text()) ]

        config['pdbin'] = str(self.modelPDBLineEdit.text())

        if self.preProcessedRadioButton.isChecked():
            config['process_engine'] = 'none'
        elif self.dtrekRadioButton.isChecked():
            config['process_engine'] = 'dstartrek'
        elif self.mosflmRadioButton.isChecked():
            config['process_engine'] = 'mosflm'
        else:
            config['process_engine'] = 'none'

        if self.detectorConstantsCheckBox.isChecked():
            config['detector_constants'] = str(self.detectorConstantsLineEdit.text())

        if self.spaceGroupComboBox.currentIndex() != 0:
            config['spacegroup_no'] = self.spaceGroupComboBox.currentIndex()

        if self.referenceDataCheckBox.isChecked():
            config['reference_mtz'] = str(self.referenceDataLineEdit.text())

        config['multi_search'] = self.multipleModelsCheckBox.isChecked()

        if self.dictionaryCIFCheckBox.isChecked():
            config['libfile'] = str(self.dictionaryCIFLineEdit.text())

        if self.htmlSummaryCheckBox.isChecked():
            config['bngsummary'] = str(self.htmlSummaryLineEdit.text())

        if self.viewpointGroupBox.isChecked():
            if self.coordinatesRadioButton.isChecked():
                config['viewpoint_method'] = 'coordinates'
                config['viewpoint_coordinates'] = [
                    str(self.xSpinBox.text()),
                    ' ',
                    str(self.ySpinBox.text()),
                    ' ',
                    str(self.zSpinBox.text()) ]

            elif self.sessionFileRadioButton.isChecked():
                config['viewpoint_method'] = 'mlw'
                config['viewpoint_file'] = str(self.sessionFileLineEdit.text())

            elif self.fromPDBMarkerFileRadioButton.isChecked():
                config['viewpoint_method'] = 'pdbview'
                config['viewpoint_file'] = str(self.fromPDBMarkerFileLineEdit.text())

        if self.placeLigandPDBGroupBox.isChecked():
            config['place_ligand'] = True
            
            if self.chemFileRadioButton.isChecked():
                config['chem_name'] = str(self.chemFileLineEdit.text())                
                
            elif self.ligandFileRadioButton.isChecked():
                config['ligand_name'] = str(self.ligandFileLineEdit.text())


        settings = QtCore.QSettings("MIFit", "MIExpert")
        settings.setValue("bng", pickle.dumps(config))

        miexpert = os.path.join(os.path.dirname(sys.argv[0]), "MIExpert.py")
        args = [ sys.executable, miexpert, "bng" ]
        args += [ "--molimagehome", str(QtCore.QFileInfo(mifit.mifitDir).absoluteFilePath()) ]
        for hklin in config['hklin']:
            args += [ "--hklin", str(QtCore.QFileInfo(hklin).absoluteFilePath()) ]
        args += [ "--pdbin", str(QtCore.QFileInfo(config['pdbin']).absoluteFilePath()) ]
        args += [ "--process_engine", str(config['process_engine']) ]

        if len(config['detector_constants']) > 0:
            args += [ '--detector_constants', config['detector_constants'] ]

        if config['spacegroup_no'] != 0:
            args += [ '--spacegroup_no', config['spacegroup_no'] ]

        if len(config['reference_mtz']) > 0:
            args += [ '--reference_mtz', config['reference_mtz'] ]

        args += [ "--multi_search" ]
        if config['multi_search']:
            args += [ 'yes' ]
        else:
            args += [ 'no' ]

        args += [ "--libfile", str(QtCore.QFileInfo(config['libfile']).absoluteFilePath()) ]

        if len(config['viewpoint_method']) > 0:
        
            if config['viewpoint_method'] == 'coordinates':
                args += [ '--frag_center', config['viewpoint_coordinates'] ]
            
            if config['viewpoint_method'] == 'mlw':
                args += [ '--mlwfile', str(QtCore.QFileInfo(config['viewpoint_file']).absoluteFilePath()) ]

            if config['viewpoint_method'] == 'pdbview':
                args += [ '--pdbviewfile', str(QtCore.QFileInfo(config['viewpoint_file']).absoluteFilePath())]               

        if len(config['bngsummary']) > 0:
            args += [ '--bngsummary', str(QtCore.QFileInfo(config['bngsummary']).absoluteFilePath()) ]           

        args += [ '--mifit', 'no' ]

        if config['place_ligand'] and len(config['viewpoint_method']) > 0:
            if self.chemFileRadioButton.isChecked():
                args += [ '--chemfit', config['chem_name'] ]                
            else:
                args += [ '--fragfit', config['ligand_name'] ]


        result = subprocess.call(args)
        exit(result)


if __name__ == '__main__':

    app = QtGui.QApplication(sys.argv)

    sgList = mifit.spacegroupList()
 
    dialog = BindNGrindDialog(sgList)    

    if dialog.exec_():
        dialog.runJob()
