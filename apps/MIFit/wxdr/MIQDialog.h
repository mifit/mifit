#ifndef MIQDIALOG_H
#define MIQDIALOG_H

#include <QDialog>

class MIQDialog : public QDialog
{
  Q_OBJECT
  public:
    MIQDialog(QWidget *parent);

  protected:
    void markEnabled(QWidget *w, bool thisEnabled, bool &globalEnabled);
};


#endif
