#include "core/corelib.h"
#include "SmilesDialog.h"

#include <QFileDialog>

void SmilesDialog::on_browseButton_clicked()
{
    fileLineEdit->setText(QFileDialog::getOpenFileName(this, "Select smiles file", "", "Smiles files (*.smi *.ism *.can)"));
}

SmilesDialog::SmilesDialog(QWidget *parent)
    : QDialog(parent)
{
    setupUi(this);
}


void SmilesDialog::GetResults(Data &data)
{
    data.mode = 0;
    if (fileRadioButton->isChecked())
    {
        data.mode = 0;
    }
    else if (smilesRadioButton->isChecked())
    {
        data.mode = 1;
    }
    else if (databaseRadioButton->isChecked())
    {
        data.mode = 2;
    }


    data.filename = fileLineEdit->text().toStdString();
    data.smiles = smilesLineEdit->text().toStdString();
    data.dbquery = databaseLineEdit->text().toStdString();
    data.code = idLineEdit->text().toStdString();
}
