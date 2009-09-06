#ifndef BVALUE_COLORS_DIALOG_H
#define BVALUE_COLORS_DIALOG_H

#include <string>

#include "core/corelib.h"

#include "ui_BValueColors.h"

class BValueColors : public QDialog, public Ui::BValueColors
{
    Q_OBJECT

public:
    BValueColors(QWidget *parent = 0);
    void InitializeFromData(const MIData &dat);
    void GetData(MIData &dat);

  private slots:
    void on_colorToolButton1_clicked() { colorButtonPressed((QToolButton*)sender()); }
    void on_colorToolButton2_clicked() { colorButtonPressed((QToolButton*)sender()); }
    void on_colorToolButton3_clicked() { colorButtonPressed((QToolButton*)sender()); }
    void on_colorToolButton4_clicked() { colorButtonPressed((QToolButton*)sender()); }
    void on_colorToolButton5_clicked() { colorButtonPressed((QToolButton*)sender()); }
    void on_colorToolButton6_clicked() { colorButtonPressed((QToolButton*)sender()); }
    void on_colorToolButton7_clicked() { colorButtonPressed((QToolButton*)sender()); }
    void on_colorToolButton8_clicked() { colorButtonPressed((QToolButton*)sender()); }
    void on_colorToolButton9_clicked() { colorButtonPressed((QToolButton*)sender()); }
    void on_colorToolButton10_clicked() { colorButtonPressed((QToolButton*)sender()); }

  private:
    std::map<QToolButton *, int> colors;
    void colorButtonPressed(QToolButton *w);
    void setColor(QToolButton *b, int i);
};

#endif
