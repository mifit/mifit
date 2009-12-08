#ifndef SELECT_CRYSTAL_DIALOG_H
#define SELECT_CRYSTAL_DIALOG_H

#include <string>

#include "core/MIData.h"
#include "ui_SelectCrystal.h"

class SelectCrystal : public QDialog, public Ui::SelectCrystal
{
    Q_OBJECT

public:
    SelectCrystal(const std::string &info, QWidget *parent = 0);
    const std::string getLabel();

    static bool doSelectCrystal(MIData& data);

private:
    void scanCrystalsDirectory();

private slots:
    void on_crystalListWidget_itemClicked(QListWidgetItem *item);
    void on_selectCrystalRadioButton_clicked();
    void on_specifyParametersRadioButton_clicked();
    void cellParamsChanged(const QString &unitCell, const QString &spaceGroup);
    void on_unitCellLineEdit_textEdited(const QString &str);
    void on_spaceGroupLineEdit_textEdited(const QString &str);
};

#endif // ifndef SELECT_CRYSTAL_DIALOG_H
