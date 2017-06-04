#include "AddResidueDialog.h"
#include "ui_AddResidueDialog.h"

#include <QtCore/QSettings>

AddResidueDialog::AddResidueDialog(const std::vector<std::string> &resList, QWidget *parent) :
    QDialog(parent),
    ui(new Ui::AddResidueDialog)
{
    ui->setupUi(this);
    setWindowTitle("Add Residue");

    foreach (const std::string &res, resList)
        ui->residueType->addItem(res.c_str());

    loadSettings();
    connect(this, SIGNAL(accepted()), this, SLOT(saveSettings()));
}

AddResidueDialog::~AddResidueDialog()
{
    delete ui;
}

QString AddResidueDialog::residueType() const
{
    return ui->residueType->currentText();
}

char AddResidueDialog::chainId() const
{
    if (!ui->chainId->text().isEmpty())
        return ui->chainId->text().at(0).toLatin1();
    return '\0';
}

AddResidueDialog::ModelLocation AddResidueDialog::modelLocation() const
{
    if (ui->beforeSelected->isChecked())
        return BeforeSelectedResidue;
    else if (ui->afterSelected->isChecked())
        return AfterSelectedResidue;
    else if (ui->startOfModel->isChecked())
        return StartOfModel;
    else if (ui->endOfModel->isChecked())
        return EndOfModel;
    return BeforeSelectedResidue;
}

AddResidueDialog::ResidueOffset AddResidueDialog::residueOffset() const
{
    if (ui->screenCenter->isChecked())
        return ScreenCenter;
    else if (ui->bestFit->isChecked())
        return BestFit;
    else if (ui->alphaHelix->isChecked())
        return AlphaHelix;
    else if (ui->betaSheet->isChecked())
        return BetaSheet;
    return ScreenCenter;
}

namespace
{
    const QString GroupKey = "AddResidueDialog";
    const QString ResidueTypeKey = "residueType";
    const QString ModelLocationKey = "modelLocation";
    const QString ResidueOffsetKey = "residueOffset";

} // anonymous namespace

void AddResidueDialog::loadSettings()
{
    QSettings settings;
    settings.beginGroup(GroupKey);
    QString residueType = settings.value(ResidueTypeKey).toString();
    ModelLocation modelLocation = static_cast<ModelLocation>(settings.value(ModelLocationKey, BeforeSelectedResidue).toInt());
    ResidueOffset residueOffset = static_cast<ResidueOffset>(settings.value(ResidueOffsetKey, ScreenCenter).toInt());
    settings.endGroup();

    ui->residueType->setCurrentIndex(ui->residueType->findText(residueType));
    switch (modelLocation)
    {
    default:
    case BeforeSelectedResidue:
        ui->beforeSelected->setChecked(true);
        break;
    case AfterSelectedResidue:
        ui->afterSelected->setChecked(true);
        break;
    case StartOfModel:
        ui->startOfModel->setChecked(true);
        break;
    case EndOfModel:
        ui->endOfModel->setChecked(true);
        break;
    }

    switch (residueOffset)
    {
    default:
    case ScreenCenter:
        ui->screenCenter->setChecked(true);
        break;
    case BestFit:
        ui->bestFit->setChecked(true);
        break;
    case AlphaHelix:
        ui->alphaHelix->setChecked(true);
        break;
    case BetaSheet:
        ui->betaSheet->setChecked(true);
        break;
    }
}

void AddResidueDialog::saveSettings()
{
    QSettings settings;
    settings.beginGroup(GroupKey);
    settings.setValue(ResidueTypeKey, residueType());
    settings.setValue(ModelLocationKey, modelLocation());
    settings.setValue(ResidueOffsetKey, residueOffset());
    settings.endGroup();
}
