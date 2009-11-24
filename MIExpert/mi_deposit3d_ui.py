import sys, os, pickle, subprocess
from PyQt4 import QtCore, QtGui, uic


class JobReportDialog(QtGui.QDialog):
    def __init__(self, parent=None):
        super(JobReportDialog, self).__init__(parent)
        self.mifitDir = QtCore.QString()

        config = {
            #  data["workdir"].str = "";
            #  data["mtzfile"].str = "";
            #  data["pdbfile"].str = "";
            #  data["libfile"].str = "none";
            #  data["seqfile"].str = "none";
            #  data["templatefile"].str = "none";
            #  data["datalogfile"].str = "none";
            #
            #  data["cif_write"].b = false; // mmCIF report
            #  data["map_write"].b = false;
            #  data["map_border"].f = FLT_MIN;
            #  data["text_report"].b = false;
            #  data["html_report"].b = false;
            #  data["hkl_write"].b = false; // mmCIF data
            #
            #  data["html_report"].b = false;
            #  data["report_title"].str = "";
            #  data["image1"].str = "";
            #
            #  data["rootname"].str = "";
          }
        settings = QtCore.QSettings("MIFit", "MIExpert")
        appSettings = settings.value("deposit3d").toString()
        if not appSettings.isEmpty():
            config = pickle.loads(str(appSettings))

        uiFile = os.path.join(os.path.dirname(sys.argv[0]), "mi_deposit3d.ui")
        uic.loadUi(uiFile, self)

        #  new MIBrowsePair(workingDirPushButton, workingDirLineEdit, "", true);
        #  new MIBrowsePair(modelPushButton, modelLineEdit,"PDB Files (*.pdb *.ent)");
        #  new MIBrowsePair(dataPushButton, dataLineEdit, "Data file (*.mtz)");
        #  new MIBrowsePair(dictionaryPushButton, dictionaryLineEdit, "Dictionary file (*.cif *.lib)");
        #  new MIBrowsePair(sequencePushButton, sequenceLineEdit, "Sequence file (*.faa *.seq *.fasta)");
        #  new MIBrowsePair(processingLogPushButton, processingLogLineEdit, "Processing log file (*.*)");
        #  new MIBrowsePair(annotationPushButton, annotationLineEdit, "Annotation file (*.*)");
        #  new MIBrowsePair(image1PushButton, image1LineEdit, "Image file (*.jpg *.gif *.tif *.jpeg *.png)");
        #  new MIBrowsePair(image2PushButton, image2LineEdit, "Image file (*.jpg *.gif *.tif *.jpeg *.png)");
        #  new MIBrowsePair(image3PushButton, image3LineEdit, "Image file (*.jpg *.gif *.tif *.jpeg *.png)");
        #
        #  _okButton = 0;
        #  QList<QAbstractButton*> buttons=buttonBox->buttons();
        #  for (int i=0; i < buttons.size(); ++i) {
        #    if (buttonBox->buttonRole(buttons[i]) == QDialogButtonBox::AcceptRole) {
        #      _okButton=buttons[i];
        #      break;
        #    }
        #  }
        #
        #  QTimer *timer = new QTimer(this);
        #  connect(timer, SIGNAL(timeout()), this, SLOT(validateTimeout()));
        #  timer->start(100);

    def runJob(self):

        #settings = QtCore.QSettings("MIFit", "MIExpert")
        #settings.setValue("deposit3d", pickle.dumps(config))

        miexpert = os.path.join(os.path.dirname(sys.argv[0]), "MIExpert.py")
        args = [ sys.executable, miexpert, "deposit3d" ]

        #  QStringList args;
        #  args << MIExpertPy() << "deposit3d"
        #          << "--workdir" << buildAbsPath(data["workdir"].str.c_str())
        #          << "--mtzfile" << buildAbsPath(data["mtzfile"].str.c_str())
        #          << "--pdbfile" << buildAbsPath(data["pdbfile"].str.c_str())
        #          << "--molimagehome" << buildAbsPath(Application::instance()->MolimageHome.c_str())
        #          << "--libfile" << buildAbsPath(data["libfile"].str.c_str())
        #          << "--seqfile" << buildAbsPath(data["seqfile"].str.c_str())
        #          << "--templatefile" << buildAbsPath(data["templatefile"].str.c_str())
        #          << "--datalogfile" << buildAbsPath(data["datalogfile"].str.c_str())
        #          << "--cif_write" << (data["cif_write"].b ? "yes" : "no");
        #
        #  if (data["map_write"].b) {
        #    args << "--map_write" << "yes"
        #            << "--map_border" << QString::number(data["map_border"].f);
        #  }
        #  args << "--text_write" << (data["text_report"].b ? "yes" : "no")
        #          << "--html_write" << (data["html_report"].b ? "yes" : "no")
        #          << "--hkl_write" << (data["hkl_write"].b ? "yes" : "no");
        #  if (data["html_report"].b) {
        #    if (data["report_title"].str.empty()) {
        #      args << "--title" << data["report_title"].str.c_str();
        #    }
        #    if (!data["image1"].str.empty()) {
        #      args << "--image1" << buildAbsPath(data["image1"].str.c_str());
        #    }
        #    if (!data["image2"].str.empty()) {
        #      args << "--image2" << buildAbsPath(data["image2"].str.c_str());
        #    }
        #    if (!data["image3"].str.empty()) {
        #      args << "--image3" << buildAbsPath(data["image3"].str.c_str());
        #    }
        #  }
        #  if (!data["rootname"].str.empty()) {
        #    args << "--rootname" << buildAbsPath(data["rootname"].str.c_str());
        #  }

        result = subprocess.call(args)
        exit(result)

#bool JobReport::GetData(MIData &data) {
#  data["workdir"].str=workingDirLineEdit->text().toStdString();
#  data["mtzfile"].str = dataLineEdit->text().toStdString();
#  data["pdbfile"].str = modelLineEdit->text().toStdString();
#
#  if (dictionaryCheckBox->isChecked() && dictionaryLineEdit->text().size())
#    data["libfile"].str = dictionaryLineEdit->text().toStdString();
#  else
#    data["libfile"].str = "none";
#
#  if (sequenceCheckBox->isChecked() && sequenceLineEdit->text().size())
#    data["seqfile"].str = sequenceLineEdit->text().toStdString();
#  else
#    data["seqfile"].str = "none";
#
#  if (annotationCheckBox->isChecked() && annotationLineEdit->text().size())
#    data["templatefile"].str = annotationLineEdit->text().toStdString();
#  else
#    data["templatefile"].str = "none";
#
#  if (processingLogCheckBox->isChecked() && processingLogLineEdit->text().size())
#    data["datalogfile"].str = processingLogLineEdit->text().toStdString();
#  else
#    data["datalogfile"].str = "none";
#
#  data["cif_write"].b = mmCIFReportCheckBox->isChecked();
#  data["map_write"].b = ligandDensityMapCheckBox->isChecked();
#  data["map_border"].f = borderSpinBox->value();
#  data["text_report"].b = textReportCheckBox->isChecked();
#  data["hkl_write"].b = mmCIFDataCheckBox->isChecked();
#
#  data["report_title"].str = "";
#  if (reportTitleCheckBox->isChecked())
#    data["report_title"].str = reportTitleLineEdit->text().toStdString();
#
#  data["image1"].str = "";
#  if (image1CheckBox->isChecked())
#    data["image1"].str = image1LineEdit->text().toStdString();
#
#  data["image2"].str = "";
#  if (image2CheckBox->isChecked())
#    data["image2"].str = image2LineEdit->text().toStdString();
#
#  data["image3"].str = "";
#  if (image3CheckBox->isChecked())
#    data["image3"].str = image3LineEdit->text().toStdString();
#
#  data["rootname"].str = "";
#  if (fileRootnameCheckBox->isChecked())
#    data["rootname"].str = fileRootnameLineEdit->text().toStdString();
#  return true;
#}

#void JobReport::validateTimeout() {
#  bool globalEnabled = true;
#
#  bool thisEnabled;
#
#  thisEnabled = (QFileInfo(workingDirLineEdit->text()).exists() &&
#                  QFileInfo(workingDirLineEdit->text()).isDir());
#  markEnabled(workingDirLineEdit, thisEnabled, globalEnabled);
#
#  thisEnabled = QFileInfo(modelLineEdit->text()).exists();
#  markEnabled(modelLineEdit, thisEnabled, globalEnabled);
#
#  thisEnabled = QFileInfo(dataLineEdit->text()).exists();
#  markEnabled(dataLineEdit, thisEnabled, globalEnabled);
#
#  thisEnabled = (!dictionaryCheckBox->isChecked() || 
#                 QFileInfo(dictionaryLineEdit->text()).exists());
#  markEnabled(dictionaryLineEdit, thisEnabled, globalEnabled);
#
#  thisEnabled = (!sequenceCheckBox->isChecked() || 
#                 QFileInfo(sequenceLineEdit->text()).exists());
#  markEnabled(sequenceLineEdit, thisEnabled, globalEnabled);
#
#  thisEnabled = (!processingLogCheckBox->isChecked() || 
#                 QFileInfo(processingLogLineEdit->text()).exists());
#  markEnabled(processingLogLineEdit, thisEnabled, globalEnabled);
#
#  thisEnabled = (!annotationCheckBox->isChecked() || 
#                 QFileInfo(annotationLineEdit->text()).exists());
#  markEnabled(annotationLineEdit, thisEnabled, globalEnabled);
#
#
#  thisEnabled = (!image1CheckBox->isChecked() || 
#                 QFileInfo(image1LineEdit->text()).exists());
#  markEnabled(image1LineEdit, thisEnabled, globalEnabled);
#
#  thisEnabled = (!image2CheckBox->isChecked() || 
#                 QFileInfo(image2LineEdit->text()).exists());
#  markEnabled(image2LineEdit, thisEnabled, globalEnabled);
#
#  thisEnabled = (!image3CheckBox->isChecked() || 
#                 QFileInfo(image3LineEdit->text()).exists());
#  markEnabled(image3LineEdit, thisEnabled, globalEnabled);
#
#
#  thisEnabled = (!reportTitleCheckBox->isChecked() || 
#                 reportTitleLineEdit->text().size());
#  markEnabled(reportTitleLineEdit, thisEnabled, globalEnabled);
#
#  thisEnabled = (!reportTitleCheckBox->isChecked() || 
#                 reportTitleLineEdit->text().size());
#  markEnabled(reportTitleLineEdit, thisEnabled, globalEnabled);
#
#  thisEnabled = (!fileRootnameCheckBox->isChecked() || 
#                 fileRootnameLineEdit->text().size());
#  markEnabled(fileRootnameLineEdit, thisEnabled, globalEnabled);
#
#  if (_okButton)
#    _okButton->setEnabled(globalEnabled);
#}

if __name__ == '__main__':

    app = QtGui.QApplication(sys.argv)

    dialog = JobReportDialog()

    if 'MIFIT_DIR' in os.environ.keys():
        dialog.mifitDir = os.environ['MIFIT_DIR']

    if 'SHELX_DIR' in os.environ.keys():
        dialog.shelxDir = os.environ['SHELX_DIR']

    if dialog.exec_():
        dialog.runJob()
