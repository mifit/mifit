#include "ContourOptions.h"
#include "MIDialog.h"
#include "ui/uilib.h"

ContourOptions::ContourOptions(QWidget *parent, bool prefsMode)
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
            MapSettings s;
            GetMapSettingsForStyle(i, s);
            MapSettingsToData(s, _styleSettings[i]);
        }
        InitializeFromData(_styleSettings[0]);
    }
    _preferencesMode = prefsMode;
}

void ContourOptions::setColor(QToolButton *valueControl, int c)
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


void ContourOptions::InitializeFromData(const MIData &dat)
{
    MIData data = dat;

    if (data.find("radiusMax") != data.end())
    {
        maxRadiusSpinBox->setValue(data["radiusMax"].i);
        radiusSlider->setMaximum(data["radiusMax"].i);
    }

    if (data.find("radius") != data.end())
    {
        radiusSpinBox->setValue((int)data["radius"].f);
        radiusSlider->setValue((int)data["radius"].f);
    }

    if (data.find("mapmin") != data.end())
    {
        int imin = (int)data["mapmin"].f;

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

    if (data.find("mapmax") != data.end())
    {
        int imax = (int)data["mapmax"].f;

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

    if (data.find("contourMethod") != data.end())
    {
        switch (data["contourMethod"].radio)
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
    }

    if (data.find("blobRadius") != data.end())
    {
        blobDistanceSpinBox->setValue(data["blobRadius"].f);
    }

    if (data.find("lineWidth") != data.end())
    {
        lineWidthSpinBox->setValue(data["lineWidth"].f);
    }

    // contour 1
    if (data.find("MapLevel0") != data.end())
    {
        contourSpinBox_0->setValue(data["MapLevel0"].f);
        contourSlider_0->setValue(data["MapLevel0"].f);
    }
    if (data.find("MapLevelOn0") != data.end())
    {
        showBox_0->setChecked(data["MapLevelOn0"].b);
    }
    if (data.find("MapLevelColor0") != data.end())
    {
        setColor(colorButton_0, data["MapLevelColor0"].i);
    }

    // contour 2
    if (data.find("MapLevel1") != data.end())
    {
        contourSpinBox_1->setValue(data["MapLevel1"].f);
        contourSlider_1->setValue(data["MapLevel1"].f);
    }
    if (data.find("MapLevelOn1") != data.end())
    {
        showBox_1->setChecked(data["MapLevelOn1"].b);
    }
    if (data.find("MapLevelColor1") != data.end())
    {
        setColor(colorButton_1, data["MapLevelColor1"].i);
    }

    // contour 3
    if (data.find("MapLevel2") != data.end())
    {
        contourSpinBox_2->setValue(data["MapLevel2"].f);
        contourSlider_2->setValue(data["MapLevel2"].f);
    }
    if (data.find("MapLevelOn2") != data.end())
    {
        showBox_2->setChecked(data["MapLevelOn2"].b);
    }
    if (data.find("MapLevelColor2") != data.end())
    {
        setColor(colorButton_2, data["MapLevelColor2"].i);
    }

    // contour 4
    if (data.find("MapLevel3") != data.end())
    {
        contourSpinBox_3->setValue(data["MapLevel3"].f);
        contourSlider_3->setValue(data["MapLevel3"].f);
    }
    if (data.find("MapLevelOn3") != data.end())
    {
        showBox_3->setChecked(data["MapLevelOn3"].b);
    }
    if (data.find("MapLevelColor3") != data.end())
    {
        setColor(colorButton_3, data["MapLevelColor3"].i);
    }

    // contour 5
    if (data.find("MapLevel4") != data.end())
    {
        contourSpinBox_4->setValue(data["MapLevel4"].f);
        contourSlider_4->setValue(data["MapLevel4"].f);
    }
    if (data.find("MapLevelOn4") != data.end())
    {
        showBox_4->setChecked(data["MapLevelOn4"].b);
    }
    if (data.find("MapLevelColor4") != data.end())
    {
        setColor(colorButton_4, data["MapLevelColor4"].i);
    }
}

void ContourOptions::GetData(MIData &data)
{
    data["radiusMax"].i = maxRadiusSpinBox->value();

    data["radius"].f = (float)radiusSpinBox->value();

    data["mapmin"].f = (float)contourSlider_0->minimum();
    data["mapmax"].f = (float)contourSlider_0->maximum();

    if (boxMethodRadioButton->isChecked())
        data["contourMethod"].radio = 0;
    if (sphereMethodRadioButton->isChecked())
        data["contourMethod"].radio = 1;
    if (blobMethodRadioButton->isChecked())
        data["contourMethod"].radio = 2;

    data["blobRadius"].f = blobDistanceSpinBox->value();
    data["lineWidth"].f = lineWidthSpinBox->value();

    // contour 1
    data["MapLevel0"].f = (float)contourSlider_0->value();
    data["MapLevelOn0"].b = showBox_0->isChecked();
    data["MapLevelColor0"].i = colors[colorButton_0];

    // contour 2
    data["MapLevel1"].f = (float)contourSlider_1->value();
    data["MapLevelOn1"].b = showBox_1->isChecked();
    data["MapLevelColor1"].i = colors[colorButton_1];

    // contour 3
    data["MapLevel2"].f = (float)contourSlider_2->value();
    data["MapLevelOn2"].b = showBox_2->isChecked();
    data["MapLevelColor2"].i = colors[colorButton_2];

    // contour 4
    data["MapLevel3"].f = (float)contourSlider_3->value();
    data["MapLevelOn3"].b = showBox_3->isChecked();
    data["MapLevelColor3"].i = colors[colorButton_3];

    // contour 5
    data["MapLevel4"].f = (float)contourSlider_4->value();
    data["MapLevelOn4"].b = showBox_4->isChecked();
    data["MapLevelColor4"].i = colors[colorButton_4];
}

void ContourOptions::on_maxRadiusSpinBox_valueChanged(int i)
{
    radiusSpinBox->setMaximum(i);
    radiusSlider->setMaximum(i);
}

void ContourOptions::on_revertPushButton_clicked()
{
    for (unsigned int i = 0; i<NUM_DEFINED_STYLES; ++i)
    {
        MapSettings s;
        GetMapSettingsForStyle(i, s, true);
        MapSettingsToData(s, _styleSettings[i]);
    }

    InitializeFromData(_styleSettings[_currentStyle]);
}


void ContourOptions::colorButtonPressed()
{
    //find relevant control
    QToolButton *control = (QToolButton*)sender();

    int ci = MIColorChooser(colors[control], "Choose color");
    setColor(control, ci);
}

void ContourOptions::on_mapStylesComboBox_currentIndexChanged(int styleNum)
{
    MIData data;
    if (!_preferencesMode)
    {
        MapSettings s;
        GetMapSettingsForStyle(styleNum, s);
        MapSettingsToData(s, data);
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


void ContourOptions::savePreferences()
{
    if (!_preferencesMode)
        return;
    GetData(_styleSettings[_currentStyle]);

    for (unsigned int i = 0; i < NUM_DEFINED_STYLES; ++i)
    {
        MapSettings s;
        float mapmin, mapmax;
        DataToMapSettings(_styleSettings[i], s, mapmin, mapmax);
        s.saveStyle(i);
    }
}
