#include "corelib.h"
#include "RefinementOptionsDialog.h"

#include <QFileDialog>

RefinementOptionsDialog::RefinementOptionsDialog(const MIData &dat, QWidget *parent)
  : QDialog(parent)
{
  setupUi(this);

  MIData data=dat;

  //variables: BondWeight.i AngleWeight.i PlaneWeight.i TorsionWeight.i
  //           BumpWeight.i MapWeight.i ConstrainCA.b ConstrainEnds.b
  //           Verbose.b RefineWhileFit.b SigmaBond.f SigmaAngle.f
  //           SigmaTorsion.f SigmaPlane.f SigmaBump.f

  bondSlider->setValue(data["BondWeight"].i);
  angleSlider->setValue(data["AngleWeight"].i);
  planeSlider->setValue(data["PlaneWeight"].i);
  torsionSlider->setValue(data["TorsionWeight"].i);
  bumpSlider->setValue(data["BumpWeight"].i);
  mapSlider->setValue(data["MapWeight"].i);

  restrainCAsCheckBox->setChecked(data["ConstrainCA"].b);
  restrainEndsCheckBox->setChecked(data["ConstrainEnds"].b);
  verboseMessagesCheckBox->setChecked(data["Verbose"].b);
  refineWhileFittingCheckBox->setChecked(data["RefineWhileFit"].b);

  bondSpinBox->setValue(data["SigmaBond"].f);
  angleSpinBox->setValue(data["SigmaAngle"].f);
  torsionSpinBox->setValue(data["SigmaTorsion"].f);
  planeSpinBox->setValue(data["SigmaPlane"].f);
  bumpSpinBox->setValue(data["SigmaBump"].f);
}


void RefinementOptionsDialog::GetResults(MIData &data) {
  data["BondWeight"].i=bondSlider->value();
  data["AngleWeight"].i=angleSlider->value();
  data["PlaneWeight"].i=planeSlider->value();
  data["TorsionWeight"].i=torsionSlider->value();
  data["BumpWeight"].i=bumpSlider->value();
  data["MapWeight"].i=mapSlider->value();

  data["ConstrainCA"].b=restrainCAsCheckBox->isChecked();
  data["ConstrainEnds"].b=restrainEndsCheckBox->isChecked();
  data["Verbose"].b=verboseMessagesCheckBox->isChecked();
  data["RefineWhileFit"].b=refineWhileFittingCheckBox->isChecked();

  data["SigmaBond"].f=bondSpinBox->value();
  data["SigmaAngle"].f=angleSpinBox->value();
  data["SigmaTorsion"].f=torsionSpinBox->value();
  data["SigmaPlane"].f=planeSpinBox->value();
  data["SigmaBump"].f=bumpSpinBox->value();
}
