#include "corelib.h"
#include "SmilesDialog.h"

#include <QFileDialog>

void SmilesDialog::on_browseButton_clicked() {
  fileLineEdit->setText(QFileDialog::getOpenFileName(this,"Select smiles file","","Smiles files (*.smi *.ism *.can)"));
}

SmilesDialog::SmilesDialog(QWidget *parent)
  : QDialog(parent)
{
  setupUi(this);
}


void SmilesDialog::GetResults(MIData &data) {
  data["mode"].radio = 0;
  if (fileRadioButton->isChecked()) {
    data["mode"].radio=0;
  }
  else if (smilesRadioButton->isChecked()) {
    data["mode"].radio=1;
  } else if (databaseRadioButton->isChecked()) {
    data["mode"].radio=2;
  }

  
  data["filename"].str = fileLineEdit->text().toStdString();
  data["smiles"].str = smilesLineEdit->text().toStdString();
  data["dbquery"].str = databaseLineEdit->text().toStdString();
  data["code"].str = idLineEdit->text().toStdString();
}
