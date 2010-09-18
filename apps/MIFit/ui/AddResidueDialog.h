#ifndef ADDRESIDUEDIALOG_H
#define ADDRESIDUEDIALOG_H

#include <QDialog>
#include <string>
#include <vector>

namespace Ui {
    class AddResidueDialog;
}

class AddResidueDialog : public QDialog
{
    Q_OBJECT

public:

    enum ModelLocation
    {
        AfterSelectedResidue,
        BeforeSelectedResidue,
        StartOfModel,
        EndOfModel
    };

    enum ResidueOffset
    {
        ScreenCenter,
        BestFit,
        AlphaHelix,
        BetaSheet
    };

    AddResidueDialog(const std::vector<std::string> &resList, QWidget *parent = 0);
    virtual ~AddResidueDialog();

    QString residueType() const;
    char chainId() const;
    ModelLocation modelLocation() const;
    ResidueOffset residueOffset() const;

private:
    Ui::AddResidueDialog *ui;

private slots:
    void loadSettings();
    void saveSettings();
};

#endif // ADDRESIDUEDIALOG_H
