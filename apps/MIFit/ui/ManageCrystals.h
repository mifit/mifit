#ifndef MANAGE_CRYSTALS_DIALOG_H
#define MANAGE_CRYSTALS_DIALOG_H

#include <string>

#include "ui_ManageCrystals.h"

class CMapHeader;

class ManageCrystals : public QDialog, public Ui::ManageCrystals
{
    Q_OBJECT

public:
    ManageCrystals(QWidget *parent = 0);
    ~ManageCrystals();

  private slots:
    void on_crystalListBox_itemClicked(QListWidgetItem *);
    void on_newCrystalButton_clicked();
    void on_copyCrystalButton_clicked();
    void on_deleteCrystalButton_clicked();
    void on_findSpaceGroupButton_clicked();
    void on_resetButton_clicked();
    void on_applyButton_clicked();
    void on_crystalName_textEdited(const QString&);
    void on_title_textEdited(const QString&);
    void on_unitCell_textEdited(const QString&);
    void on_spaceGroup_textEdited(const QString&);
    void on_symmOps_textChanged();
    void on_ncsOps_textChanged();

    void dialog_finished(int);

  private:
    void scanCrystalsDirectory();
    void updateDetails();
    void updateSpacegroupDetails();
    void updateButtons();
    void lookupSpaceGroup();
    void loadSelectedCrystal();
    void saveCrystal(bool deleteOldCrystal = true);
    void saveAs(const std::string& name);
    bool crystalFileExists(std::string& crystal);
    void deleteCrystal(std::string& crystal);
    void clearDetails();
    void promptSaveIfModified();

    CMapHeader *mapHeader;
    bool modified;
    bool fieldModified;
};

#endif
