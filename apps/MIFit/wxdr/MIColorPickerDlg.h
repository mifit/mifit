#ifndef MICOLORPICKER_H
#define MICOLORPICKER_H

#include <QDialog>
#include <QToolButton>
#include <QColor>
#include <QPalette>
#include <QPixmap>

class MIColorToolButton : public QToolButton
{
    Q_OBJECT
    Q_PROPERTY(QColor color READ color WRITE setColor)
      public:

    MIColorToolButton(QWidget *parent = 0)
      : QToolButton(parent), col(QColor()) { setColor(col); }

    MIColorToolButton(const QColor &c, QWidget *parent = 0)
      : QToolButton(parent), col(c) { setColor(col); }

    QColor color() const
    {
      return col;
    }

    void setColor(const QColor &cc)
    {
      col = cc;
      QPixmap p(32,32);
      p.fill(col);
      QIcon icon(p);
      setIcon(icon);

      update();
    }

    QSize sizeHint() const
    {
      return QSize(32, 32);
    }

  private:
    QColor col;
};


class MIColorPickerDlg : public QDialog
{
  Q_OBJECT

  public:
    MIColorPickerDlg(QWidget *parent = 0, int selected = 0);

    unsigned int GetResult() const { return result; }

 public Q_SLOTS:
    virtual void colorPicked(int i);

  private:
    unsigned int result;
};

#endif
