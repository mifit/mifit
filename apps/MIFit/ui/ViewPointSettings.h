#ifndef VIEWPOINTSETTINGS_H
#define VIEWPOINTSETTINGS_H

class ViewPointSettings
{
public:
    ViewPointSettings();

    //definitions for ballandstick
    static const int STICKS;
    static const int BALLANDSTICK;
    static const int CPK;
    static const int BALLANDCYLINDER;

    void setDepthCuedLineWidth(bool on);
    bool isDepthCuedLineWidth();
    void setDepthCuedColors(bool on);
    bool isDepthCuedColors();
    void setDimNonactiveModels(bool on);
    bool isDimNonactiveModels();
    void setAmountToDimNonactiveModels(float percent);
    float getAmountToDimNonactiveModels();
    int GetLineThickness();
    void SetLineThickness(int w);
    void SetBallandStick();
    void SetBallandStick(int s);
    void SetBallandCylinder();
    void SetSpaceFilling();
    void SetSticks();
    int GetBallandStick();
    int GetBallSize();
    void SetBallSize(int b);
    int GetCylinderSize();
    void SetCylinderSize(int b);

private:
    bool depthcue_color;
    bool depthcue_width;
    bool dimNonactiveModels;
    float amountToDimNonactiveModels;
    int ballandstick;
    int ballsize;
    int cylindersize;
    int lineThickness;

};

#endif // VIEWPOINTSETTINGS_H
