#include "ContourOptionsWidget.h"
#include "ui/uilib.h"
#include "ui/MIColorPickerDlg.h"

ContourOptionsWidget::ContourOptionsWidget(QWidget *parent, bool prefsMode)
    : QWidget(parent)
{
    setupUi(this);

    _currentStyle = 0;
    connect(colorButton_0, SIGNAL(clicked()), this, SLOT(colorButtonPressed()));
    connect(colorButton_1, SIGNAL(clicked()), this, SLOT(colorButtonPressed()));
    connect(colorButton_2, SIGNAL(clicked()), this, SLOT(colorButtonPressed()));
    connect(colorButton_3, SIGNAL(clicked()), this, SLOT(colorButtonPressed()));
    connect(colorButton_4, SIGNAL(clicked()), this, SLOT(colorButtonPressed()));

    if (!prefsMode)
    {
        revertPushButton->setVisible(false);
    }
    else
    {
        for (unsigned int i = 0; i<NUM_DEFINED_STYLES; ++i)
        {
            GetMapSettingsForStyle(i, _styleSettings[i]);
        }
        InitializeFromData(_styleSettings[0]);
    }
    _preferencesMode = prefsMode;
}

void ContourOptionsWidget::setColor(QToolButton *valueControl, int c)
{
    colors[valueControl] = c;

    MIPalette *m_palette = Application::instance()->GetLPpal();
    int ci = PaletteIndex(c);
    QColor color(m_palette->colors[ci].red,
                 m_palette->colors[ci].green,
                 m_palette->colors[ci].blue);

    // color the background with the color
    QPalette palette = valueControl->palette();
    palette.setColor(QPalette::Window, color);
    valueControl->setPalette(palette);

    // for mac, we have to set a solid-color pixmap too, b/c the
    // background palette color is ignored however, if we *just* use a
    // pixmap, there's no (straightforward) way to retreive the color.
    QPixmap p(32, 32);
    p.fill(color);
    QIcon icon(p);
    valueControl->setIcon(icon);
}


void ContourOptionsWidget::InitializeFromData(const MapSettingsBase &data)
{
    maxRadiusSpinBox->setValue(data.m_radiusmax);
    radiusSlider->setMaximum(data.m_radiusmax);

    radiusSpinBox->setValue((int)data.Radius);
    radiusSlider->setValue((int)data.Radius);

    if (data.mapmin != 0.0f)
    {
        int imin = (int)data.mapmin;

        contourSpinBox_0->setMinimum(imin);
        contourSlider_0->setMinimum(imin);
        minLabel_0->setNum(imin);

        contourSpinBox_1->setMinimum(imin);
        contourSlider_1->setMinimum(imin);
        minLabel_1->setNum(imin);

        contourSpinBox_2->setMinimum(imin);
        contourSlider_2->setMinimum(imin);
        minLabel_2->setNum(imin);

        contourSpinBox_3->setMinimum(imin);
        contourSlider_3->setMinimum(imin);
        minLabel_3->setNum(imin);

        contourSpinBox_4->setMinimum(imin);
        contourSlider_4->setMinimum(imin);
        minLabel_4->setNum(imin);
    }

    if (data.mapmax != 0.0f)
    {
        int imax = (int)data.mapmax;

        contourSpinBox_0->setMaximum(imax);
        contourSlider_0->setMaximum(imax);
        maxLabel_0->setNum(imax);

        contourSpinBox_1->setMaximum(imax);
        contourSlider_1->setMaximum(imax);
        maxLabel_1->setNum(imax);

        contourSpinBox_2->setMaximum(imax);
        contourSlider_2->setMaximum(imax);
        maxLabel_2->setNum(imax);

        contourSpinBox_3->setMaximum(imax);
        contourSlider_3->setMaximum(imax);
        maxLabel_3->setNum(imax);

        contourSpinBox_4->setMaximum(imax);
        contourSlider_4->setMaximum(imax);
        maxLabel_4->setNum(imax);
    }

    switch (data.ContourMethod)
    {
    case 0:
        boxMethodRadioButton->setChecked(true);
        break;
    case 1:
        sphereMethodRadioButton->setChecked(true);
        break;
    case 2:
        blobMethodRadioButton->setChecked(true);
        break;
    }

    blobDistanceSpinBox->setValue(data.BlobRadius);

    lineWidthSpinBox->setValue(data.maplinewidth);

    // contour 1
    contourSpinBox_0->setValue(data.MapLevel[0]);
    contourSlider_0->setValue(data.MapLevel[0]);
    showBox_0->setChecked(data.MapLevelOn[0]);
    setColor(colorButton_0, data.MapLevelColor[0]);

    // contour 2
    contourSpinBox_1->setValue(data.MapLevel[1]);
    contourSlider_1->setValue(data.MapLevel[1]);
    showBox_1->setChecked(data.MapLevelOn[1]);
    setColor(colorButton_1, data.MapLevelColor[1]);

    // contour 3
    contourSpinBox_2->setValue(data.MapLevel[2]);
    contourSlider_2->setValue(data.MapLevel[2]);
    showBox_2->setChecked(data.MapLevelOn[2]);
    setColor(colorButton_2, data.MapLevelColor[2]);

    // contour 4
    contourSpinBox_3->setValue(data.MapLevel[3]);
    contourSlider_3->setValue(data.MapLevel[3]);
    showBox_3->setChecked(data.MapLevelOn[3]);
    setColor(colorButton_3, data.MapLevelColor[3]);

    // contour 5
    contourSpinBox_4->setValue(data.MapLevel[4]);
    contourSlider_4->setValue(data.MapLevel[4]);
    showBox_4->setChecked(data.MapLevelOn[4]);
    setColor(colorButton_4, data.MapLevelColor[4]);
}

void ContourOptionsWidget::GetData(MapSettingsBase &data)
{
    data.m_radiusmax = maxRadiusSpinBox->value();

    data.Radius = (float)radiusSpinBox->value();

    data.mapmin = (float)contourSlider_0->minimum();
    data.mapmax = (float)contourSlider_0->maximum();

    if (boxMethodRadioButton->isChecked())
        data.ContourMethod = 0;
    if (sphereMethodRadioButton->isChecked())
        data.ContourMethod = 1;
    if (blobMethodRadioButton->isChecked())
        data.ContourMethod = 2;

    data.BlobRadius = blobDistanceSpinBox->value();
    data.maplinewidth = lineWidthSpinBox->value();

    // contour 1
    data.MapLevel[0] = (float)contourSlider_0->value();
    data.MapLevelOn[0] = showBox_0->isChecked();
    data.MapLevelColor[0] = colors[colorButton_0];

    // contour 2
    data.MapLevel[1] = (float)contourSlider_1->value();
    data.MapLevelOn[1] = showBox_1->isChecked();
    data.MapLevelColor[1] = colors[colorButton_1];

    // contour 3
    data.MapLevel[2] = (float)contourSlider_2->value();
    data.MapLevelOn[2] = showBox_2->isChecked();
    data.MapLevelColor[2] = colors[colorButton_2];

    // contour 4
    data.MapLevel[3] = (float)contourSlider_3->value();
    data.MapLevelOn[3] = showBox_3->isChecked();
    data.MapLevelColor[3] = colors[colorButton_3];

    // contour 5
    data.MapLevel[4] = (float)contourSlider_4->value();
    data.MapLevelOn[4] = showBox_4->isChecked();
    data.MapLevelColor[4] = colors[colorButton_4];
}

void ContourOptionsWidget::on_maxRadiusSpinBox_valueChanged(int i)
{
    radiusSpinBox->setMaximum(i);
    radiusSlider->setMaximum(i);
}

void ContourOptionsWidget::on_revertPushButton_clicked()
{
    for (unsigned int i = 0; i<NUM_DEFINED_STYLES; ++i)
    {
        MapSettingsBase s;
        GetMapSettingsForStyle(i, s, true);
        _styleSettings[i] = s;
    }

    InitializeFromData(_styleSettings[_currentStyle]);
}


void ContourOptionsWidget::colorButtonPressed()
{
    //find relevant control
    QToolButton *control = (QToolButton*)sender();

    int ci = MIColorPickerDlg::getColor(this, colors[control], "Choose color");
    setColor(control, ci);
}

void ContourOptionsWidget::on_mapStylesComboBox_currentIndexChanged(int styleNum)
{
    MapSettingsBase data;
    if (!_preferencesMode)
    {
        GetMapSettingsForStyle(styleNum, data);
    }
    else
    {
        GetData(data);
        _styleSettings[_currentStyle] = data; // save current settings
        data = _styleSettings[styleNum]; // get new settings
    }
    InitializeFromData(data);
    _currentStyle = styleNum;
}


void ContourOptionsWidget::savePreferences()
{
    if (!_preferencesMode)
        return;
    GetData(_styleSettings[_currentStyle]);

    for (unsigned int i = 0; i < NUM_DEFINED_STYLES; ++i)
    {
        MapSettings::saveStyle(i, _styleSettings[i]);
    }
}
