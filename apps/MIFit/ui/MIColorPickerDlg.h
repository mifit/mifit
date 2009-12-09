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
public :

        MIColorToolButton(QWidget *parent = 0)
            : QToolButton(parent), col(QColor())
    {
        setColor(col);
    }

    MIColorToolButton(const QColor &c, QWidget *parent = 0)
        : QToolButton(parent), col(c)
    {
        setColor(col);
    }

    QColor color() const
    {
        return col;
    }

    void setColor(const QColor &cc)
    {
        col = cc;
        QPixmap p(32, 32);
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

    static int getColor(QWidget *parent = 0, int selected = 0, const QString& title = QString("Choose color"))
    {
        MIColorPickerDlg dlg(parent, selected, title);
        if (dlg.exec() == QDialog::Accepted)
        {
            return dlg.color();
        }
        return -1;
    }

    MIColorPickerDlg(QWidget *parent = 0, int selected = 0, const QString& title = QString("Choose color"));

    int color() const
    {
        return result;
    }

    unsigned int GetResult() const
    {
        return result;
    }

public slots:
    virtual void colorPicked(int i);

private:
    int result;
};

#endif // ifndef MICOLORPICKER_H
