#include "RefinementDialog.h"
#include "MIBrowsePair.h"

#include "core/corelib.h"
#include <map/maplib.h>
#include "ui/uilib.h"
#include "MIDialog.h"

#include <vector>
#include <string>

#include <QDialogButtonBox>
#include <QAbstractButton>
#include <QList>
#include <QTimer>
#include <QFileInfo>


RefinementDialog::RefinementDialog(QWidget *parent) : MIQDialog(parent) {
  setupUi(this);

  new MIBrowsePair(workingDirPushButton, workingDirLineEdit, "", true);
  new MIBrowsePair(modelPushButton, modelLineEdit,"PDB Files (*.pdb *.ent)");
  new MIBrowsePair(dataPushButton, dataLineEdit, "Data file (*.mtz)");
  new MIBrowsePair(tlsSpecificationPushButton, tlsSpecificationLineEdit,"Select a TLS file (*.tls)");
  new MIBrowsePair(dictionaryPushButton, dictionaryLineEdit, "Dictionary file (*.cif *.lib)");

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

void RefinementDialog::GetInitialData(MIData &data) { // static
  data["workdir"].str = "";
  data["pdbfile"].str = "";
  data["mtzfile"].str = "";
  data["weight"].f = FLT_MIN;
  data["cycles"].u = UINT_MAX;

  data["water_cycles"].u = UINT_MAX;
  data["build_cycles"].u = UINT_MAX;
  data["bref_type"].str = "";
  data["engine"].str = "";

  data["use_max_res"].b = false;
  data["max_res"].f = FLT_MIN;

  data["libfile"].str = "";
  data["tls_file"].str = "";
}


void RefinementDialog::validateTimeout() {
  bool globalEnabled = true;

  bool thisEnabled;

  thisEnabled = (QFileInfo(workingDirLineEdit->text()).exists() &&
                 QFileInfo(workingDirLineEdit->text()).isDir());
  markEnabled(workingDirLineEdit, thisEnabled, globalEnabled);

  thisEnabled = QFileInfo(modelLineEdit->text()).exists();
  markEnabled(modelLineEdit, thisEnabled, globalEnabled);

  thisEnabled = QFileInfo(dataLineEdit->text()).exists();
  markEnabled(dataLineEdit, thisEnabled, globalEnabled);

  thisEnabled = (!tlsSpecificationCheckBox->isChecked() || 
                 QFileInfo(tlsSpecificationLineEdit->text()).exists());
  markEnabled(tlsSpecificationLineEdit, thisEnabled, globalEnabled);

  thisEnabled = (!dictionaryCheckBox->isChecked() || 
                 QFileInfo(dictionaryLineEdit->text()).exists());
  markEnabled(dictionaryLineEdit, thisEnabled, globalEnabled);

  if (_okButton)
    _okButton->setEnabled(globalEnabled);
}


void RefinementDialog::InitializeFromData(const MIData &)
{
}

bool RefinementDialog::GetData(MIData &data) {
  data["workdir"].str = workingDirLineEdit->text().toStdString();
  data["pdbfile"].str = modelLineEdit->text().toStdString();
  data["mtzfile"].str = dataLineEdit->text().toStdString();
  data["weight"].f = weightSpinBox->value();
  data["cycles"].u = numCyclesSpinBox->value();

  data["water_cycles"].u = UINT_MAX;
  if (waterPickingCheckBox->isChecked())
    data["water_cycles"].u = waterPickingCyclesSpinBox->value();

  data["build_cycles"].u = UINT_MAX;

  data["bref_type"].str = "isotropic";
  if (anisotropicRadioButton->isChecked())
    data["bref_type"].str = "anisotropic";
    
  data["engine"].str = "refmac5";
  if (rigidBodyRadioButton->isChecked())
    data["engine"].str = "rigid";
  if (shelxRadioButton->isChecked())
    data["engine"].str = "shelx";
    
  data["use_max_res"].b = maxResolutionCheckBox->isChecked();
  data["max_res"].f = maxResolutionSpinBox->value();

  data["tls_file"].str = "";
  if (tlsSpecificationCheckBox->isChecked())
    data["tls_file"].str = tlsSpecificationLineEdit->text().toStdString();
    
  data["libfile"].str = "";
  if (dictionaryCheckBox->isChecked())
    data["libfile"].str = dictionaryLineEdit->text().toStdString();
  return true;
}

