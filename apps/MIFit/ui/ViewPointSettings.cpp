#include "ViewPointSettings.h"


const int ViewPointSettings::STICKS = 0;
const int ViewPointSettings::BALLANDSTICK = 1;
const int ViewPointSettings::CPK = 2;
const int ViewPointSettings::BALLANDCYLINDER = 3;

ViewPointSettings::ViewPointSettings()
    : depthcue_color(true),
      depthcue_width(true),
      dimNonactiveModels(true),
      amountToDimNonactiveModels(50),
      ballandstick(0),
      ballsize(10),
      cylindersize(33),
      lineThickness(1)
{
}


void ViewPointSettings::setDepthCuedLineWidth(bool on)
{
    depthcue_width = on;
}

bool ViewPointSettings::isDepthCuedLineWidth()
{
    return depthcue_width;
}

void ViewPointSettings::setDepthCuedColors(bool on)
{
    depthcue_color = on;
}

bool ViewPointSettings::isDepthCuedColors()
{
    return depthcue_color;
}

void ViewPointSettings::setDimNonactiveModels(bool on)
{
    dimNonactiveModels = on;
}

bool ViewPointSettings::isDimNonactiveModels()
{
    return dimNonactiveModels;
}

void ViewPointSettings::setAmountToDimNonactiveModels(float percent)
{
    amountToDimNonactiveModels = percent;
}

float ViewPointSettings::getAmountToDimNonactiveModels()
{
    return amountToDimNonactiveModels;
}

void ViewPointSettings::SetBallandStick()
{
    ballandstick = BALLANDSTICK;
}

void ViewPointSettings::SetBallandStick(int s)
{
    ballandstick = s;
}

void ViewPointSettings::SetBallandCylinder()
{
    ballandstick = BALLANDCYLINDER;
}

void ViewPointSettings::SetSpaceFilling()
{
    ballandstick = CPK;
}

void ViewPointSettings::SetSticks()
{
    ballandstick = STICKS;
}

int ViewPointSettings::GetBallandStick()
{
    return ballandstick;
}

int ViewPointSettings::GetBallSize()
{
    return ballsize;
}

void ViewPointSettings::SetBallSize(int b)
{
    ballsize = b;
}

int ViewPointSettings::GetCylinderSize()
{
    return cylindersize;
}

void ViewPointSettings::SetCylinderSize(int b)
{
    cylindersize = b;
}


int ViewPointSettings::GetLineThickness()
{
    return lineThickness;
}

void ViewPointSettings::SetLineThickness(int w)
{
    lineThickness = w;
}

