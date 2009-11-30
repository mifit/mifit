#include "MIQDialog.h"

MIQDialog::MIQDialog(QWidget *parent)
    : QDialog(parent)
{
}

void MIQDialog::markEnabled(QWidget *w, bool thisEnabled, bool &globalEnabled)
{
    w->setAutoFillBackground(true);
    QPalette pal = palette();
    if (!thisEnabled)
    {
        pal.setColor(QPalette::Normal, QPalette::Base, QColor(255, 255, 128));
    }
    w->setPalette(pal);

    globalEnabled &= thisEnabled;
}
