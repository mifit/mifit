import sys, os, pickle, subprocess
from PyQt4 import QtCore, QtGui, uic


class SadPhasingDialog(QtGui.QDialog):
    def __init__(self, parent=None):
        super(SadPhasingDialog, self).__init__(parent)
        self.mifitDir = QtCore.QString()

        config = {
#        data["workdir"].str="";
#        data["saddatafile"].str = "";
#
#        data["sitefile"].str = "none";
#        data["scatterer"].str = "";
#
#        data["sitenumber"].u=UINT_MAX;
#        data["ssnumber"].u=UINT_MAX;
#        data["separation"].f=FLT_MIN;
#        data["solventfraction"].f=FLT_MIN;
#        data["siterefinemethod"].str="";
#        data["bothhands"].b = false;
#
#        std::vector<std::string> sg;
#        MIGetSpacegroups(sg);
#        data["spacegroup_no"].radio = 0;
#        data["spacegroup_no"].radio_count = sg.size();
#        data["spacegroup_no"].radio_labels = sg;
#        data["change_spacegroup"].b = false;
          }
        settings = QtCore.QSettings("MIFit", "MIExpert")
        appSettings = settings.value("sadphase").toString()
        if not appSettings.isEmpty():
            config = pickle.loads(str(appSettings))

        uiFile = os.path.join(os.path.dirname(sys.argv[0]), "mi_sadphase.ui")
        uic.loadUi(uiFile, self)


    def runJob(self):

        #settings = QtCore.QSettings("MIFit", "MIExpert")
        #settings.setValue("sadphase", pickle.dumps(config))

        miexpert = os.path.join(os.path.dirname(sys.argv[0]), "MIExpert.py")
        args = [ sys.executable, miexpert, "sadphase" ]


        #  QStringList args;
        #  args << MIExpertPy() << "sadphase";
        #  args << "--molimagehome \"" + buildAbsPath(Application::instance()->MolimageHome.c_str());
        #  args << "--workdir" << buildAbsPath(data["workdir"].str.c_str());
        #  args << "--saddatafile" << buildAbsPath(data["saddatafile"].str.c_str());
        #  args << "--sitefile" << buildAbsPath(data["sitefile"].str.c_str());
        #  args << "--scatterer" << data["scatterer"].str.c_str();
        #  args << "--sitenumber" << QString::number(data["sitenumber"].u);
        #  args << "--ssnumber" << QString::number(data["ssnumber"].u);
        #  args << "--separation" << QString::number(data["separation"].f);
        #  args << "--solventfraction" << QString::number(data["solventfraction"].f);
        #  args << "--siterefinemethod" << data["siterefinemethod"].str.c_str();
        #  if (!Application::instance()->ShelxHome.empty()) {
        #    args << "--shelx_dir" << buildAbsPath(Application::instance()->ShelxHome.c_str());
        #  }
        #  args << "--bothhands" << (data["bothhands"].b ? "yes" : "no");
        #  if (data["change_spacegroup"].b ) {
        #    args << "--spacegroup_no" << QString::number(data["spacegroup_no"].radio);
        #  }

        result = subprocess.call(args)
        exit(result)

#void SadPhasing::validateTimeout() {
#  bool globalEnabled = true;
#
#  bool thisEnabled=QFileInfo(workingDirectoryLineEdit->text()).exists() &&
#    QFileInfo(workingDirectoryLineEdit->text()).isDir();
#  markEnabled(workingDirectoryLineEdit, thisEnabled, globalEnabled);
#
#  thisEnabled=QFileInfo(intensityDataLineEdit->text()).exists();
#  markEnabled(intensityDataLineEdit, thisEnabled, globalEnabled);
#
#  thisEnabled = (!sitesFromFileRadioButton->isChecked() || 
#                 QFileInfo(sitesFromFileLineEdit->text()).exists());
#  markEnabled(sitesFromFileLineEdit, thisEnabled, globalEnabled);
#
#  thisEnabled = scattererTypeLineEdit->text().size();
#  markEnabled(scattererTypeLineEdit, thisEnabled, globalEnabled);
#
#
#  if (_okButton)
#    _okButton->setEnabled(globalEnabled);
#}
#
#
#void SadPhasing::InitializeFromData(const MIData &data)
#{
#  MIData dat=data;
#  // populate spacegroup combo box from dat;
#  spaceGroupComboBox->clear();
#  for (size_t i=0; i< dat["spacegroup_no"].radio_labels.size(); ++i) {
#    spaceGroupComboBox->addItem(dat["spacegroup_no"].radio_labels[i].c_str());
#  }
#  spaceGroupComboBox->setCurrentIndex(0);
#}
#
#bool SadPhasing::GetData(MIData &data) {
#  data["workdir"].str=workingDirectoryLineEdit->text().toStdString();
#  data["saddatafile"].str = intensityDataLineEdit->text().toStdString();
#
#  data["sitefile"].str = "none";
#  if (sitesFromFileRadioButton->isChecked())
#    data["sitefile"].str = sitesFromFileLineEdit->text().toStdString();
#
#  data["scatterer"].str = scattererTypeLineEdit->text().toStdString();
#
#  data["sitenumber"].u=numSitesSpinBox->value();
#  data["ssnumber"].u=numDisulfidesSpinBox->value();
#  data["separation"].f=minScattererSeparationSpinBox->value();
#  data["solventfraction"].f=solventFractionSpinBox->value();
#  if (bp3RadioButton->isChecked())
#    data["siterefinemethod"].str="bp3";
#  else
#    data["siterefinemethod"].str="mlphare";
#  data["bothhands"].b = phaseBothSiteEnantiomorphsCheckBox->isChecked();
#
#  data["change_spacegroup"].b = changeSpaceGroupCheckBox->isChecked();
#  data["spacegroup_no"].radio = spaceGroupComboBox->currentIndex() + 1; // first spacegroup is 1, first index is 0
#  return true;
#}

if __name__ == '__main__':

    app = QtGui.QApplication(sys.argv)

    dialog = SadPhasingDialog()

    if 'MIFIT_DIR' in os.environ.keys():
        dialog.mifitDir = os.environ['MIFIT_DIR']

    if 'SHELX_DIR' in os.environ.keys():
        dialog.shelxDir = os.environ['SHELX_DIR']

    if dialog.exec_():
        dialog.runJob()
