#include "ManageCrystals.h"

#include "ui/uilib.h"
#include <map/maplib.h>
#include "MIDialog.h"

#include <QDir>
#include <QFileInfo>
#include <QMessageBox>


ManageCrystals::ManageCrystals(QWidget *parent)
    : QDialog(parent)
{
    mapHeader = new CMapHeader;
    modified = false;
    fieldModified = false;

    setupUi(this);

    updateButtons();
    scanCrystalsDirectory();

    connect(this, SIGNAL(finished(int)), this, SLOT(dialog_finished(int)));
}

ManageCrystals::~ManageCrystals()
{
    delete mapHeader;
}

void ManageCrystals::scanCrystalsDirectory()
{
    std::string crystal_data = Application::instance()->GetCrystalData();
    QDir crystal_dir(crystal_data.c_str());
    if (!crystal_dir.exists())
        return;
    crystalListBox->clear();
    QFileInfoList file_list = crystal_dir.entryInfoList();
    for (int i = 0; i<file_list.count(); ++i)
    {
        QFileInfo &fi = file_list[i];
        QString fn = fi.fileName();
        if (fn[0] != '.' && fn != QString("crystals"))
        {
            crystalListBox->addItem(fn);
        }
    }
    if (crystalListBox->count())
        crystalListBox->setCurrentRow(0);
}

void ManageCrystals::on_crystalListBox_itemClicked(QListWidgetItem*)
{
    promptSaveIfModified();
    loadSelectedCrystal();
}

void ManageCrystals::updateDetails()
{
    crystalName->setText(mapHeader->crystal_name.c_str());
    title->setText(MIAfterFirst(mapHeader->title, ' ').c_str());
    std::string unitCellString = ::format("%0.2f %0.2f %0.2f %0.2f %0.2f %0.2f", mapHeader->a, mapHeader->b, mapHeader->c, mapHeader->alpha, mapHeader->beta, mapHeader->gamma);
    unitCell->setText(unitCellString.c_str());
    updateSpacegroupDetails();

    ncsOps->clear();
    if (mapHeader->nNCRSymmops > 0)
    {
        std::string fmt;
        std::string ncsOpsString;
        for (int i = 0; i < mapHeader->nNCRSymmops; i++)
        {
            fmt = ::format("%0.4f %0.4f %0.4f\n%0.4f %0.4f %0.4f\n%0.4f %0.4f %0.4f\n%0.2f %0.2f %0.2f\n",
                           mapHeader->NCRSymmops[i][0], mapHeader->NCRSymmops[i][1], mapHeader->NCRSymmops[i][2],
                           mapHeader->NCRSymmops[i][3], mapHeader->NCRSymmops[i][4], mapHeader->NCRSymmops[i][5],
                           mapHeader->NCRSymmops[i][6], mapHeader->NCRSymmops[i][7], mapHeader->NCRSymmops[i][8],
                           mapHeader->NCRSymmops[i][9], mapHeader->NCRSymmops[i][10], mapHeader->NCRSymmops[i][11]);
            ncsOpsString += fmt;
        }
        ncsOps->setText(ncsOpsString.c_str());
    }
}

void ManageCrystals::updateSpacegroupDetails()
{
    std::string spacegroupString;
    if (mapHeader->spgpname.length() > 1 && mapHeader->spgpno > 0)
    {
        spacegroupString = ::format("%s %d", mapHeader->spgpname.c_str(), mapHeader->spgpno);
    }
    else
    {
        spacegroupString = ::format("%s", mapHeader->spgpname.c_str());
    }
    spaceGroup->setText(spacegroupString.c_str());

    symmOps->clear();
    if (mapHeader->nsym > 0)
    {
        for (int i = 0; i < mapHeader->nsym; i++)
        {
            symmOps->append(mapHeader->SymopsString[i].c_str());
            //symmOps->append("\n");
        }
    }
}

void ManageCrystals::on_applyButton_clicked()
{
    saveCrystal();
}

void ManageCrystals::on_resetButton_clicked()
{
    loadSelectedCrystal();
}

void ManageCrystals::on_findSpaceGroupButton_clicked()
{
    lookupSpaceGroup();
}



void ManageCrystals::on_crystalName_textEdited(const QString&)
{
    updateButtons();
}
void ManageCrystals::on_title_textEdited(const QString&)
{
    updateButtons();
}
void ManageCrystals::on_unitCell_textEdited(const QString&)
{
    updateButtons();
}
void ManageCrystals::on_spaceGroup_textEdited(const QString&)
{
    updateButtons();
}
void ManageCrystals::on_symmOps_textChanged()
{
    updateButtons();
}
void ManageCrystals::on_ncsOps_textChanged()
{
    updateButtons();
}

void ManageCrystals::updateButtons()
{
    fieldModified = crystalName->isModified()
                    || title->isModified()
                    || unitCell->isModified()
                    || spaceGroup->isModified()
                    || ncsOps->document()->isModified();
    applyButton->setEnabled(modified || fieldModified);
    resetButton->setEnabled(modified || fieldModified);
    findSpaceGroupButton->setEnabled(spaceGroup->isModified());
    bool crystalSelected = (crystalListBox->currentItem() && crystalListBox->currentItem()->text().toStdString().size() > 0);
    copyCrystalButton->setEnabled(crystalSelected);
    deleteCrystalButton->setEnabled(crystalSelected);

}

void ManageCrystals::lookupSpaceGroup()
{
    if (!spaceGroup->isModified())
    {
        return;
    }
    std::string s = spaceGroup->text().toStdString();
    s = MIBeforeFirst(s, ' ');
    int l;
    if (sscanf(s.c_str(), "%d", &l) == 1)
    {
        if (l > 0 && l <= 230)
        {
            mapHeader->spgpno = l;
            mapHeader->SetSymmOps();
            modified = true;
            updateDetails();
        }
    }
    else
    {
        if (mapHeader->FindSpacegroup(s.c_str()))
        {
            modified = true;
            updateDetails();
        }
    }
}

void ManageCrystals::loadSelectedCrystal()
{
    if (!crystalListBox->currentItem())
        return;
    std::string selection = crystalListBox->currentItem()->text().toStdString();
    mapHeader->LoadCrystal(selection.c_str());
    modified = false;
    updateDetails();
    updateButtons();
}

void ManageCrystals::saveCrystal(bool deleteOldCrystal)
{
    std::string oldCrystalName = mapHeader->crystal_name.c_str();
    std::string newCrystalName = crystalName->text().toStdString();
    mapHeader->crystal_name = newCrystalName;
    mapHeader->title = "name " + title->text().toStdString();
    sscanf(unitCell->text().toStdString().c_str(), "%f%f%f%f%f%f", &mapHeader->a, &mapHeader->b, &mapHeader->c, &mapHeader->alpha, &mapHeader->beta, &mapHeader->gamma);
    lookupSpaceGroup();
    mapHeader->nNCRSymmops = 0;
    std::string ncsstr = ncsOps->toPlainText().toStdString();
    std::vector<std::string> results;
    MIStringSplit(ncsstr, "\n", results);
    for (unsigned int i = 0; i < results.size(); i += 4)
    {
        std::string &ncrline1 = results[i];
        std::string &ncrline2 = results[i+1];
        std::string &ncrline3 = results[i+2];
        std::string &ncrline4 = results[i+3];
        float s1, s2, s3, s4, s5, s6, s7, s8, s9, s10, s11, s12;
        if (sscanf(ncrline1.c_str(), "%f%f%f", &s1, &s2, &s3) == 3
            && sscanf(ncrline2.c_str(), "%f%f%f", &s4, &s5, &s6) == 3
            && sscanf(ncrline3.c_str(), "%f%f%f", &s7, &s8, &s9) == 3
            && sscanf(ncrline4.c_str(), "%f%f%f", &s10, &s11, &s12) == 3)
        {

            mapHeader->NCRSymmops[mapHeader->nNCRSymmops][0] = s1;
            mapHeader->NCRSymmops[mapHeader->nNCRSymmops][1] = s2;
            mapHeader->NCRSymmops[mapHeader->nNCRSymmops][2] = s3;
            mapHeader->NCRSymmops[mapHeader->nNCRSymmops][3] = s4;
            mapHeader->NCRSymmops[mapHeader->nNCRSymmops][4] = s5;
            mapHeader->NCRSymmops[mapHeader->nNCRSymmops][5] = s6;
            mapHeader->NCRSymmops[mapHeader->nNCRSymmops][6] = s7;
            mapHeader->NCRSymmops[mapHeader->nNCRSymmops][7] = s8;
            mapHeader->NCRSymmops[mapHeader->nNCRSymmops][8] = s9;
            mapHeader->NCRSymmops[mapHeader->nNCRSymmops][9] = s10;
            mapHeader->NCRSymmops[mapHeader->nNCRSymmops][10] = s11;
            mapHeader->NCRSymmops[mapHeader->nNCRSymmops][11] = s12;
            mapHeader->nNCRSymmops++;
        }
    }
    mapHeader->SaveCrystal(mapHeader->crystal_name);
    if (deleteOldCrystal && oldCrystalName != newCrystalName)
    {
        deleteCrystal(oldCrystalName);
        // find newCrystalName in crystalListBox and select it
        QList<QListWidgetItem*> items = crystalListBox->findItems(newCrystalName.c_str(), Qt::MatchFixedString);
        if (items.count() != 0)
            crystalListBox->setCurrentItem(items[0]);
        loadSelectedCrystal();
    }
    modified = false;
    updateDetails();
    updateButtons();
}

void ManageCrystals::on_newCrystalButton_clicked()
{
    clearDetails();
    saveAs("new");
}

void ManageCrystals::on_copyCrystalButton_clicked()
{
    std::string crystalCopy = mapHeader->crystal_name.c_str();
    if (strncmp(MIAfterLast(crystalCopy, '_').c_str(), "copy", 4)==0)
    {
        crystalCopy = MIBeforeLast(crystalCopy, '_');
    }
    crystalCopy += "_copy";
    saveAs(crystalCopy);
}

void ManageCrystals::saveAs(const std::string &name)
{
    std::string crystal = name;
    std::string fmt = name + "%d";
    int count = 1;
    while (crystalFileExists(crystal))
    {
        crystal = ::format(fmt.c_str(), count);
        ++count;
    }
    crystalName->setText(crystal.c_str());
    saveCrystal(false);
    scanCrystalsDirectory();
    QList<QListWidgetItem*> items = crystalListBox->findItems(crystal.c_str(), Qt::MatchFixedString);
    if (items.count() != 0)
        crystalListBox->setCurrentItem(items[0]);
    loadSelectedCrystal();
}

bool ManageCrystals::crystalFileExists(std::string &crystal)
{
    QFileInfo file(QDir(Application::instance()->GetCrystalData().c_str()), crystal.c_str());
    return file.exists();
}

void ManageCrystals::on_deleteCrystalButton_clicked()
{
    if (!crystalListBox->currentItem())
        return;
    std::string selection = crystalListBox->currentItem()->text().toStdString();
    deleteCrystal(selection);
    clearDetails();
    loadSelectedCrystal();
}

void ManageCrystals::deleteCrystal(std::string &crystal)
{
    if (crystalFileExists(crystal))
    {
        QFileInfo finfo(QDir(Application::instance()->GetCrystalData().c_str()), crystal.c_str());
        QFile file(finfo.absoluteFilePath());
        file.remove();
    }
    scanCrystalsDirectory();
}

void ManageCrystals::clearDetails()
{
    CMapHeader blank;
    *mapHeader = blank;
    crystalName->clear();
    title->clear();
    unitCell->clear();
    spaceGroup->clear();
    symmOps->clear();
    ncsOps->clear();
}

void ManageCrystals::dialog_finished(int)
{
    promptSaveIfModified();
}

void ManageCrystals::promptSaveIfModified()
{
    updateButtons();
    if (modified || fieldModified)
    {
        if (QMessageBox::question(this, "Save crystal?", "The current crystal has been modifed.\nShould the changes be saved?", QMessageBox::Yes | QMessageBox::No) == QMessageBox::Yes)
        {
            saveCrystal();
        }
    }
}


