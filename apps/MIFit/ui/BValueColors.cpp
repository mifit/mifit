#include "BValueColors.h"
#include "MIDialog.h"
#include "ui/uilib.h"

BValueColors::BValueColors(QWidget *parent)
    : QDialog(parent)
{
    setupUi(this);
}

void BValueColors::setColor(QToolButton *valueControl, int c)
{
    colors[valueControl] = c;

    MIPalette *m_palette = Application::instance()->GetLPpal();
    int pi = PaletteIndex(c);
    QColor color(m_palette->colors[pi].red,
                 m_palette->colors[pi].green,
                 m_palette->colors[pi].blue);

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


void BValueColors::colorButtonPressed(QToolButton *control)
{
    int ci = MIColorChooser(colors[control], "Choose color");
    setColor(control, ci);
}


void BValueColors::GetData(MIData &data)
{
    data["level1"].f = (float)levelSpinBox_1->value();
    data["level2"].f = (float)levelSpinBox_2->value();
    data["level3"].f = (float)levelSpinBox_3->value();
    data["level4"].f = (float)levelSpinBox_4->value();
    data["level5"].f = (float)levelSpinBox_5->value();
    data["level6"].f = (float)levelSpinBox_6->value();
    data["level7"].f = (float)levelSpinBox_7->value();
    data["level8"].f = (float)levelSpinBox_8->value();
    data["level9"].f = (float)levelSpinBox_9->value();
    data["level10"].f = (float)levelSpinBox_10->value();

    data["color1"].i = colors[colorToolButton1];
    data["color2"].i = colors[colorToolButton2];
    data["color3"].i = colors[colorToolButton3];
    data["color4"].i = colors[colorToolButton4];
    data["color5"].i = colors[colorToolButton5];
    data["color6"].i = colors[colorToolButton6];
    data["color7"].i = colors[colorToolButton7];
    data["color8"].i = colors[colorToolButton8];
    data["color9"].i = colors[colorToolButton9];
    data["color10"].i = colors[colorToolButton10];

    data["save"].b = saveCheckBox->isChecked();
}

void BValueColors::InitializeFromData(const MIData &dat)
{
    MIData data = dat;

    levelSpinBox_1->setValue(data["level1"].f);
    levelSpinBox_2->setValue(data["level2"].f);
    levelSpinBox_3->setValue(data["level3"].f);
    levelSpinBox_4->setValue(data["level4"].f);
    levelSpinBox_5->setValue(data["level5"].f);
    levelSpinBox_6->setValue(data["level6"].f);
    levelSpinBox_7->setValue(data["level7"].f);
    levelSpinBox_8->setValue(data["level8"].f);
    levelSpinBox_9->setValue(data["level9"].f);
    levelSpinBox_10->setValue(data["level10"].f);

    setColor(colorToolButton1, data["color1"].i);
    setColor(colorToolButton2, data["color2"].i);
    setColor(colorToolButton3, data["color3"].i);
    setColor(colorToolButton4, data["color4"].i);
    setColor(colorToolButton5, data["color5"].i);
    setColor(colorToolButton6, data["color6"].i);
    setColor(colorToolButton7, data["color7"].i);
    setColor(colorToolButton8, data["color8"].i);
    setColor(colorToolButton9, data["color9"].i);
    setColor(colorToolButton10, data["color10"].i);
}

