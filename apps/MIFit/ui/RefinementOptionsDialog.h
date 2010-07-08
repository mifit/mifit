#ifndef REFINEMENTOPTIONSDIALOG_H
#define REFINEMENTOPTIONSDIALOG_H

#include "ui_RefinementOptionsDialog.h"

class RefinementOptionsDialog : public QDialog, public Ui::RefinementOptionsDialog
{
    Q_OBJECT

public:

    struct Data
    {
        int bondWeight;
        int angleWeight;
        int planeWeight;
        int torsionWeight;
        int bumpWeight;
        int mapWeight;

        bool constrainCA;
        bool constrainEnds;
        bool verbose;
        bool refineWhileFit;

        float sigmaBond;
        float sigmaAngle;
        float sigmaPlane;
        float sigmaTorsion;
        float sigmaBump;

    };

    RefinementOptionsDialog(const Data &dat, QWidget *parent = 0);
    void GetResults(Data &data);
};

#endif
