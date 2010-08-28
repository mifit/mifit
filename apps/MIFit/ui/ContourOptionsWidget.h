#ifndef CONTOUR_OPTIONS_WIDGET_H
#define CONTOUR_OPTIONS_WIDGET_H

#include "core/corelib.h"
#include "MapSettings.h"

#include "ui_ContourOptionsWidget.h"

const unsigned int NUM_DEFINED_STYLES = 6;

class ContourOptionsWidget : public QWidget, public Ui::ContourOptionsWidget
{
    Q_OBJECT

public:
    ContourOptionsWidget(QWidget *parent = 0, bool prefsMode = false);
    void InitializeFromData(const MapSettingsBase &dat);
    void GetData(MapSettingsBase &dat);
    void savePreferences();

public slots:
    void on_maxRadiusSpinBox_valueChanged(int i);
    void colorButtonPressed();
    void on_mapStylesComboBox_currentIndexChanged(int styleNum);
    void on_revertPushButton_clicked();

private:
    void setColor(QToolButton *valueControl, int c);
    std::map<QWidget*, int> colors;
    MapSettingsBase _styleSettings[NUM_DEFINED_STYLES];
    bool _preferencesMode;
    int _currentStyle;
};

#endif // CONTOUR_OPTIONS_WIDGET_H
