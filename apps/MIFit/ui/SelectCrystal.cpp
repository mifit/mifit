#include "SelectCrystal.h"

#include "ui/uilib.h"

#include <QDir>
#include <QFileInfo>

void SelectCrystal::scanCrystalsDirectory()
{
    std::string crystal_data = Application::instance()->GetCrystalData();
    QDir crystal_dir(crystal_data.c_str());
    if (!crystal_dir.exists())
        return;
    QFileInfoList file_list = crystal_dir.entryInfoList();
    for (int i = 0; i<file_list.count(); ++i)
    {
        QFileInfo &fi = file_list[i];
        QString fn = fi.fileName();
        if (fn[0] != '.' && fn != QString("crystals"))
        {
            crystalListWidget->addItem(fn);
        }
    }
}

void SelectCrystal::on_crystalListWidget_itemClicked(QListWidgetItem *item)
{
    CMapHeader mh;
    mh.LoadCrystal(item->text().toStdString().c_str());
    valueLineEdit->setText(mh.Label().c_str());
}

void SelectCrystal::on_selectCrystalRadioButton_clicked()
{
    crystalListWidget->setEnabled(true);
    unitCellLineEdit->setEnabled(false);
    spaceGroupLineEdit->setEnabled(false);

    if (crystalListWidget->selectedItems().count()==0
        && crystalListWidget->count()!=0)
    {
        crystalListWidget->setCurrentRow(0);
        on_crystalListWidget_itemClicked(crystalListWidget->currentItem());
    }
}

void SelectCrystal::on_specifyParametersRadioButton_clicked()
{
    crystalListWidget->setEnabled(false);
    unitCellLineEdit->setEnabled(true);
    spaceGroupLineEdit->setEnabled(true);
}

void SelectCrystal::cellParamsChanged(const QString &unitCell, const QString &spaceGroup)
{
    CMapHeader mh;
    mh.a = 100.0f;
    mh.b = 100.0f;
    mh.c = 100.0f;
    mh.alpha = 90.0f;
    mh.beta = 90.0f;
    mh.gamma = 90.0f;
    sscanf(unitCell.toStdString().c_str(), "%f%f%f%f%f%f", &mh.a, &mh.b, &mh.c, &mh.alpha, &mh.beta, &mh.gamma);

    mh.spgpno = 1;
    std::string s = spaceGroup.toStdString();
    int sg;
    if (sscanf(s.c_str(), "%d", &sg) == 1)
    {
        if (sg > 0 &&  sg<= 230)
        {
            mh.spgpno = sg;
        }
    }
    else
    {
        mh.FindSpacegroup(s.c_str());
    }

    mh.SetSymmOps();
    valueLineEdit->setText(mh.Label().c_str());
}


void SelectCrystal::on_unitCellLineEdit_textEdited(const QString &str)
{
    cellParamsChanged(str, spaceGroupLineEdit->text());
}

void SelectCrystal::on_spaceGroupLineEdit_textEdited(const QString &str)
{
    cellParamsChanged(unitCellLineEdit->text(), str);
}


SelectCrystal::SelectCrystal(const std::string &info, QWidget *parent)
    : QDialog(parent)
{
    setupUi(this);

    crystalListWidget->setEnabled(false);
    scanCrystalsDirectory();

    CMapHeaderBase mh(info);
    char buf[128];
    sprintf(buf, "%0.2f %0.2f %0.2f %0.1f %0.1f %0.1f", mh.a, mh.b, mh.c, mh.alpha, mh.beta, mh.gamma);
    unitCellLineEdit->setText(buf);
    spaceGroupLineEdit->setText(mh.spgpname.c_str());
    cellParamsChanged(unitCellLineEdit->text(), spaceGroupLineEdit->text());
}

const std::string SelectCrystal::getLabel()
{
    return valueLineEdit->text().toStdString();
}
