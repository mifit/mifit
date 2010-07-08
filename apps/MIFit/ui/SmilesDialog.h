#ifndef SMILESDIALOG_H
#define SMILESDIALOG_H

#include "ui_SmilesDialog.h"
#include <string>

class SmilesDialog : public QDialog, public Ui::SmilesDialog
{
    Q_OBJECT

public:

    struct Data
    {
        int mode;
        std::string filename;
        std::string code;
        std::string smiles;
        std::string dbquery;
    };

    SmilesDialog(QWidget *parent = 0);
    void GetResults(Data &data);

public slots:
    void on_browseButton_clicked();
};

#endif // SMILESDIALOG_H
