#include "core/corelib.h"
#include "RefinementOptionsDialog.h"

#include <QFileDialog>

RefinementOptionsDialog::RefinementOptionsDialog(const Data &data, QWidget *parent)
    : QDialog(parent)
{
    setupUi(this);

    bondSlider->setValue(data.bondWeight);
    angleSlider->setValue(data.angleWeight);
    planeSlider->setValue(data.planeWeight);
    torsionSlider->setValue(data.torsionWeight);
    bumpSlider->setValue(data.bumpWeight);
    mapSlider->setValue(data.mapWeight);

    restrainCAsCheckBox->setChecked(data.constrainCA);
    restrainEndsCheckBox->setChecked(data.constrainEnds);
    verboseMessagesCheckBox->setChecked(data.verbose);
    refineWhileFittingCheckBox->setChecked(data.refineWhileFit);

    bondSpinBox->setValue(data.sigmaBond);
    angleSpinBox->setValue(data.sigmaAngle);
    torsionSpinBox->setValue(data.sigmaTorsion);
    planeSpinBox->setValue(data.sigmaPlane);
    bumpSpinBox->setValue(data.sigmaBump);
}


void RefinementOptionsDialog::GetResults(Data &data)
{
    data.bondWeight = bondSlider->value();
    data.angleWeight = angleSlider->value();
    data.planeWeight = planeSlider->value();
    data.torsionWeight = torsionSlider->value();
    data.bumpWeight = bumpSlider->value();
    data.mapWeight = mapSlider->value();

    data.constrainCA = restrainCAsCheckBox->isChecked();
    data.constrainEnds = restrainEndsCheckBox->isChecked();
    data.verbose = verboseMessagesCheckBox->isChecked();
    data.refineWhileFit = refineWhileFittingCheckBox->isChecked();

    data.sigmaBond = bondSpinBox->value();
    data.sigmaAngle = angleSpinBox->value();
    data.sigmaTorsion = torsionSpinBox->value();
    data.sigmaPlane = planeSpinBox->value();
    data.sigmaBump = bumpSpinBox->value();
}
