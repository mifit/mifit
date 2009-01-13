#include "BindNGrind.h"
#include "MIBrowsePair.h"

#include "corelib.h"
#include "maplib.h"
#include "uilib.h"
#include "MIDialog.h"

#include <vector>
#include <string>

#include <QDialogButtonBox>
#include <QAbstractButton>
#include <QList>
#include <QTimer>
#include <QFileInfo>


BindNGrind::BindNGrind(QWidget *parent) : MIQDialog(parent) {
  setupUi(this);

  new MIBrowsePair(modelPDBPushButton, modelPDBLineEdit,"PDB Files (*.pdb *.ent)");
  new MIBrowsePair(detectorConstantsPushButton, detectorConstantsLineEdit,"Constants file (*.txt)");
  new MIBrowsePair(referenceDataPushButton, referenceDataLineEdit, "Data file (*.mtz)");
  new MIBrowsePair(htmlSummaryPushButton, htmlSummaryLineEdit, "", true);
  new MIBrowsePair(dictionaryCIFPushButton, dictionaryCIFLineEdit, "Dictionary file (*.cif *.lib)");
  new MIBrowsePair(sessionFilePushButton, sessionFileLineEdit, "Session file (*.mlw)");
  new MIBrowsePair(fromPDBMarkerFilePushButton, fromPDBMarkerFileLineEdit, "PDB File (*.pdb *.ent)");
  new MIBrowsePair(ligandFilePushButton, ligandFileLineEdit, "PDB Files (*.pdb *.ent)");
  
  _okButton = 0;
  QList<QAbstractButton*> buttons=buttonBox->buttons();
  for (int i=0; i < buttons.size(); ++i) {
    if (buttonBox->buttonRole(buttons[i]) == QDialogButtonBox::AcceptRole) {
      _okButton=buttons[i];
      break;
    }
  }

  QTimer *timer = new QTimer(this);
  connect(timer, SIGNAL(timeout()), this, SLOT(validateTimeout()));
  timer->start(100);
}

// create data array with all proper elements, but unset.
// this is for encapsulation purposes --- instead of having the knowlege
// about which fields will be set spread across multiple files, it's all
// contained here.

void BindNGrind::GetInitialData(MIData &data) { // static
  data["hklin"].strList = std::vector<std::string>();
  data["pdbin"].str = "";
  data["process_engine"].str = "none";
  data["arpwarpmap"].b = false;
  data["detector_constants"].str = "";

  std::vector<std::string> sg;
  sg.push_back("Unchanged");
  MIGetSpacegroups(sg);

  data["spacegroup_no"].radio = 0; 
  data["spacegroup_no"].radio_count = sg.size();
  data["spacegroup_no"].radio_labels = sg;

  data["reference_mtz"].str = "none";
  data["multi_search"].b = false;
  data["libfile"].str = "none";
  data["viewpoint_method"].str="";
  data["bngsummary"].str="";
  data["place_ligand"].b = false;
  data["ligand_name"].str = "";
  data["ligand_from_dictionary"].b = false;
}


void BindNGrind::validateTimeout() {
  bool globalEnabled = true;

  bool thisEnabled=intensityDataListWidget->count() > 0;
  markEnabled(intensityDataListWidget, thisEnabled, globalEnabled);

  thisEnabled = QFileInfo(modelPDBLineEdit->text()).exists();
  markEnabled(modelPDBLineEdit, thisEnabled, globalEnabled);

  thisEnabled = (!detectorConstantsCheckBox->isChecked() || 
                 QFileInfo(detectorConstantsLineEdit->text()).exists());
  markEnabled(detectorConstantsLineEdit, thisEnabled, globalEnabled);

  thisEnabled = (!referenceDataCheckBox->isChecked() || 
                 QFileInfo(referenceDataLineEdit->text()).exists());
  markEnabled(referenceDataLineEdit, thisEnabled, globalEnabled);

  thisEnabled = (!htmlSummaryCheckBox->isChecked() || 
                 (QFileInfo(htmlSummaryLineEdit->text()).exists() &&
                  QFileInfo(htmlSummaryLineEdit->text()).isDir()));
  markEnabled(htmlSummaryLineEdit, thisEnabled, globalEnabled);

  thisEnabled = (!dictionaryCIFCheckBox->isChecked() || 
                 QFileInfo(dictionaryCIFLineEdit->text()).exists());
  markEnabled(dictionaryCIFLineEdit, thisEnabled, globalEnabled);

  thisEnabled = (!dictionaryCIFCheckBox->isChecked() || 
                 QFileInfo(dictionaryCIFLineEdit->text()).exists());
  markEnabled(dictionaryCIFLineEdit, thisEnabled, globalEnabled);


  if (viewpointGroupBox->isChecked()) {
    thisEnabled = (!sessionFileRadioButton->isChecked() || 
                   QFileInfo(sessionFileLineEdit->text()).exists());
    markEnabled(sessionFileLineEdit, thisEnabled, globalEnabled);

    thisEnabled = (!fromPDBMarkerFileRadioButton->isChecked() || 
                   QFileInfo(fromPDBMarkerFileLineEdit->text()).exists());
    markEnabled(fromPDBMarkerFileLineEdit, thisEnabled, globalEnabled);

    if (placeLigandPDBGroupBox->isChecked()) {
      thisEnabled = (!ligandFileCheckBox->isChecked() || 
                     QFileInfo(ligandFileLineEdit->text()).exists());
      markEnabled(ligandFileLineEdit, thisEnabled, globalEnabled);
    }
  }


  if (_okButton)
    _okButton->setEnabled(globalEnabled);
}


void BindNGrind::InitializeFromData(const MIData &data)
{
  MIData dat=data;
  // populate spacegroup combo box from dat;
  spaceGroupComboBox->clear();
  for (size_t i=0; i< dat["spacegroup_no"].radio_labels.size(); ++i) {
    spaceGroupComboBox->addItem(dat["spacegroup_no"].radio_labels[i].c_str());
  }
  spaceGroupComboBox->setCurrentIndex(0);

  // populate dictionary combobox
  dictionaryComboBox->clear();
  std::vector<std::string> choices = MIFitDictionary()->GetDictResList();
  for (size_t i=0; i< choices.size(); ++i) {
    dictionaryComboBox->addItem(choices[i].c_str());
  }
  dictionaryComboBox->setCurrentIndex(0);
}

bool BindNGrind::GetData(MIData &data) {

  //convert intensityDataListWidget contents to strlist;
  std::vector<std::string> hklin;
  for (int i=0;i<intensityDataListWidget->count(); ++i) {
    hklin.push_back(intensityDataListWidget->item(i)->text().toStdString());
  }
  data["hklin"].strList = hklin;


  data["pdbin"].str = modelPDBLineEdit->text().toStdString();
  if (preProcessedRadioButton->isChecked())
    data["process_engine"].str = "none";
  else if (dtrekRadioButton->isChecked())
    data["process_engine"].str = "dstartrek";
  else
    data["process_engine"].str = "mosflm";

  data["arpwarpmap"].b = arpWarpCheckBox->isChecked();
  if (detectorConstantsCheckBox->isChecked() && detectorConstantsLineEdit->text().size())
    data["detector_constants"].str = detectorConstantsLineEdit->text().toStdString();
  else
    data["detector_constants"].str = "";
  data["spacegroup_no"].radio = 
    (spaceGroupComboBox->currentIndex() == -1 ? 0 : spaceGroupComboBox->currentIndex());

  if (referenceDataCheckBox->isChecked() && referenceDataLineEdit->text().size())
    data["reference_mtz"].str = referenceDataLineEdit->text().toStdString();
  else
    data["reference_mtz"].str = "none";
  data["multi_search"].b = multipleModelsCheckBox->isChecked();

  if (dictionaryCIFCheckBox->isChecked() && dictionaryCIFLineEdit->text().size())
    data["libfile"].str = dictionaryCIFLineEdit->text().toStdString();
  else
    data["libfile"].str = "none";

  data["bngsummary"].str="";
  if (htmlSummaryCheckBox->isChecked() && htmlSummaryLineEdit->text().size()) {
    data["bngsummary"].str=htmlSummaryLineEdit->text().toStdString();
  }


  data["viewpoint_method"].str="";
  if (viewpointGroupBox->isChecked()) {
    if (coordinatesRadioButton->isChecked()) {
      data["viewpoint_method"].str = " --frag_center \"" +
        xSpinBox->text().toStdString() + " " +
        ySpinBox->text().toStdString() + " " +
        zSpinBox->text().toStdString() + " " + "\"";
    } else if (sessionFileRadioButton->isChecked()) {
      data["viewpoint_method"].str =
        " --mlwfile \"" + sessionFileLineEdit->text().toStdString() + "\"";
    } else  {
      data["viewpoint_method"].str =
        " --pdbviewfile \"" + fromPDBMarkerFileLineEdit->text().toStdString() + "\"";
    }
  }

  data["place_ligand"].b = false;
  data["ligand_from_dictionary"].b = false;
  if (placeLigandPDBGroupBox->isChecked()) {
    data["place_ligand"].b = true;
    if (ligandFileCheckBox->isChecked()) {
      data["ligand_name"].str = ligandFileLineEdit->text().toStdString();
    } else if(dictionaryRadioButton->isChecked()) {
      data["ligand_name"].str = dictionaryComboBox->currentText().toStdString();
      data["ligand_from_dictionary"].b = true;
    } else {
      data["ligand_name"].str="all";
    }
  }
  return true;
}

void BindNGrind::on_removePushButton_clicked() {
  int isel = intensityDataListWidget->currentRow();
  if (isel >= 0) {
    QListWidgetItem *item=intensityDataListWidget->takeItem(isel);
    delete item;
  }
}

void BindNGrind::on_addPushButton_clicked() {
  static std::string path = "";
  std::string file = MIFileSelector("Select intensity file", path, "", "",
    "Data files (*.ref, *.sca, *.mtz, *.img, *.osc)|*.ref;*.sca;*.mtz;*.img;*.osc"
    "|CCP4 files (*.mtz)|*.mtz"
    "|Reflection files (*.ref, *.sca)|*.ref;*.sca"
    "|Image files (*.img, *.osc)|*.img;*.osc"
    "|All files (*.*)|*.*", 0, 0);
  if (file.size() == 0) {
    return;
  }
  intensityDataListWidget->addItem(file.c_str());
}
