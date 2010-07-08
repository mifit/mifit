#include "BValueColors.h"
#include "ui/uilib.h"
#include "ui/MIColorPickerDlg.h"

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
    int ci = MIColorPickerDlg::getColor(this, colors[control], "Choose color");
    setColor(control, ci);
}


void BValueColors::GetData(ColorValues &data)
{
    data.BValueRanges[0] = (float)levelSpinBox_1->value();
    data.BValueRanges[1] = (float)levelSpinBox_2->value();
    data.BValueRanges[2] = (float)levelSpinBox_3->value();
    data.BValueRanges[3] = (float)levelSpinBox_4->value();
    data.BValueRanges[4] = (float)levelSpinBox_5->value();
    data.BValueRanges[5] = (float)levelSpinBox_6->value();
    data.BValueRanges[6] = (float)levelSpinBox_7->value();
    data.BValueRanges[7] = (float)levelSpinBox_8->value();
    data.BValueRanges[8] = (float)levelSpinBox_9->value();
    data.BValueRanges[9] = (float)levelSpinBox_10->value();

    data.BValueColors[0] = colors[colorToolButton1];
    data.BValueColors[1] = colors[colorToolButton2];
    data.BValueColors[2] = colors[colorToolButton3];
    data.BValueColors[3] = colors[colorToolButton4];
    data.BValueColors[4] = colors[colorToolButton5];
    data.BValueColors[5] = colors[colorToolButton6];
    data.BValueColors[6] = colors[colorToolButton7];
    data.BValueColors[7] = colors[colorToolButton8];
    data.BValueColors[8] = colors[colorToolButton9];
    data.BValueColors[9] = colors[colorToolButton10];

    data.save = saveCheckBox->isChecked();
}

void BValueColors::InitializeFromData(const ColorValues &dat)
{
    levelSpinBox_1->setValue(dat.BValueRanges[0]);
    levelSpinBox_2->setValue(dat.BValueRanges[1]);
    levelSpinBox_3->setValue(dat.BValueRanges[2]);
    levelSpinBox_4->setValue(dat.BValueRanges[3]);
    levelSpinBox_5->setValue(dat.BValueRanges[4]);
    levelSpinBox_6->setValue(dat.BValueRanges[5]);
    levelSpinBox_7->setValue(dat.BValueRanges[6]);
    levelSpinBox_8->setValue(dat.BValueRanges[7]);
    levelSpinBox_9->setValue(dat.BValueRanges[8]);
    levelSpinBox_10->setValue(dat.BValueRanges[9]);

    setColor(colorToolButton1, dat.BValueColors[0]);
    setColor(colorToolButton2, dat.BValueColors[1]);
    setColor(colorToolButton3, dat.BValueColors[2]);
    setColor(colorToolButton4, dat.BValueColors[3]);
    setColor(colorToolButton5, dat.BValueColors[4]);
    setColor(colorToolButton6, dat.BValueColors[5]);
    setColor(colorToolButton7, dat.BValueColors[6]);
    setColor(colorToolButton8, dat.BValueColors[7]);
    setColor(colorToolButton9, dat.BValueColors[8]);
    setColor(colorToolButton10, dat.BValueColors[9]);
}

