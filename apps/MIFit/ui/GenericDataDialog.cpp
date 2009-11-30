#include "GenericDataDialog.h"

#include <QApplication>
#include <QGridLayout>
#include <QLabel>
#include <QSpinBox>
#include <QLineEdit>
#include <QCheckBox>
#include <QDialogButtonBox>
#include <QToolButton>
#include <QColorDialog>
#include <QComboBox>
#include <QFormLayout>

#include <nongui/nongui.h>
#include <util/utillib.h>
#include "core/corelib.h"
#include "ui/uilib.h"
#include "MIDialog.h"

namespace
{
void setButtonColor(QAbstractButton *button, const QColor &color)
{
    QPalette palette = button->palette();
    palette.setColor(QPalette::Window, color);
    button->setPalette(palette);

    QPixmap p(32, 32);
    p.fill(color);
    QIcon icon(p);
    button->setIcon(icon);
}


class ColorIndexButton
    : public QToolButton
{
    Q_OBJECT
public:
    ColorIndexButton(QWidget *parent = 0)
        : QToolButton(parent)
    {
    }

    int colorIndex;
};
}

#include "GenericDataDialog.moc"

GenericDataDialog::GenericDataDialog(QWidget *parent, Qt::WindowFlags f)
    : QDialog(parent, f)
{
    QVBoxLayout *layout = new QVBoxLayout;

    formLayout_ = new QFormLayout;
    layout->addLayout(formLayout_, 1);

    QDialogButtonBox *buttons = new QDialogButtonBox;
    buttons->addButton(QDialogButtonBox::Ok);
    buttons->addButton(QDialogButtonBox::Cancel);
    layout->addWidget(buttons);
    QObject::connect(buttons, SIGNAL(accepted()), this, SLOT(accept()));
    QObject::connect(buttons, SIGNAL(rejected()), this, SLOT(reject()));

    setLayout(layout);
}

void GenericDataDialog::addBoolField(const QString &label, bool value)
{
    QCheckBox *widget = new QCheckBox(this);
    widget->setChecked(value);
    fieldTypes_ += Bool;
    values_ += value;
    formLayout_->addRow(label, widget);
}

void GenericDataDialog::addIntField(const QString &label, int value)
{
    QLineEdit *widget = new QLineEdit(this);
    widget->setText(QString::number(value));
    fieldTypes_ += Int;
    values_ += value;
    formLayout_->addRow(label, widget);
}

void GenericDataDialog::addUIntField(const QString &label, unsigned int value)
{
    QLineEdit *widget = new QLineEdit(this);
    widget->setText(QString::number(value));
    fieldTypes_ += UInt;
    values_ += value;
    formLayout_->addRow(label, widget);
}

void GenericDataDialog::addDoubleField(const QString &label, double value)
{
    QLineEdit *widget = new QLineEdit(this);
    widget->setText(QString::number(value, 'f'));
    fieldTypes_ += Double;
    values_ += value;
    formLayout_->addRow(label, widget);
}

void GenericDataDialog::addStringField(const QString &label, const QString &value)
{
    QLineEdit *widget = new QLineEdit(this);
    widget->setText(value);
    fieldTypes_ += String;
    values_ += value;
    formLayout_->addRow(label, widget);
}

void GenericDataDialog::addComboField(const QString &label, const QStringList &choices, int index)
{
    QComboBox *widget = new QComboBox(this);
    widget->addItems(choices);
    widget->setCurrentIndex(index);
    fieldTypes_ += Combo;
    values_ += index;
    formLayout_->addRow(label, widget);
}

void GenericDataDialog::addColorField(const QString &label, const QColor &color)
{
    QToolButton *widget = new QToolButton(this);
    connect(widget, SIGNAL(clicked()), SLOT(colorButtonPressed()));
    setButtonColor(widget, color);
    fieldTypes_ += Color;
    values_ += color;
    formLayout_->addRow(label, widget);
}

void GenericDataDialog::addColorIndexField(const QString &label, int colorIndex)
{
    MIPalette *palette = Application::instance()->GetLPpal();
    int ci = PaletteIndex(colorIndex);
    QColor color(palette->colors[ci].red,
                 palette->colors[ci].green,
                 palette->colors[ci].blue);

    ColorIndexButton *widget = new ColorIndexButton(this);
    widget->colorIndex = colorIndex;
    connect(widget, SIGNAL(clicked()), SLOT(colorIndexButtonPressed()));
    setButtonColor(widget, color);
    fieldTypes_ += ColorIndex;
    values_ += colorIndex;
    formLayout_->addRow(label, widget);
}

QVariant GenericDataDialog::value(int index) const
{
    return values_.at(index);
}

void GenericDataDialog::accept()
{
    for (int i = 0; i < formLayout_->count(); ++i)
    {
        QLayoutItem *item = formLayout_->itemAt(i, QFormLayout::FieldRole);
        if (item)
        {
            QWidget *widget = item->widget();
            switch (fieldTypes_.at(i))
            {
            case Bool:
            {
                QCheckBox *checkBox = static_cast<QCheckBox*>(widget);
                values_[i] = checkBox->isChecked();
            }
            break;
            case Int:
            {
                QLineEdit *lineEdit = static_cast<QLineEdit*>(widget);
                values_[i] = lineEdit->text().toInt();
            }
            break;
            case UInt:
            {
                QLineEdit *lineEdit = static_cast<QLineEdit*>(widget);
                values_[i] = lineEdit->text().toUInt();
            }
            break;
            case Double:
            {
                QLineEdit *lineEdit = static_cast<QLineEdit*>(widget);
                values_[i] = lineEdit->text().toDouble();
            }
            break;
            case String:
            {
                QLineEdit *lineEdit = static_cast<QLineEdit*>(widget);
                values_[i] = lineEdit->text();
            }
            break;
            case Combo:
            {
                QComboBox *comboBox = static_cast<QComboBox*>(widget);
                values_[i] = comboBox->currentIndex();
            }
            break;
            case Color:
            {
                QToolButton *button = static_cast<QToolButton*>(widget);
                values_[i] = button->palette().color(QPalette::Normal, QPalette::Window);
            }
            break;
            case ColorIndex:
            {
                ColorIndexButton *button = static_cast<ColorIndexButton*>(widget);
                values_[i] = button->colorIndex;
            }
            }
        }
    }
    QDialog::accept();
}

void GenericDataDialog::colorIndexButtonPressed()
{
    if (sender()->isWidgetType())
    {
        ColorIndexButton *widget = qobject_cast<ColorIndexButton*>(sender());
        if (widget)
        {
            widget->colorIndex = static_cast<unsigned char>(MIColorChooser(widget->colorIndex));

            MIPalette *palette = Application::instance()->GetLPpal();
            int ci = PaletteIndex(widget->colorIndex);
            QColor color(palette->colors[ci].red,
                         palette->colors[ci].green,
                         palette->colors[ci].blue);

            setButtonColor(widget, color);
        }
    }
}

void GenericDataDialog::colorButtonPressed()
{
    if (sender()->isWidgetType())
    {
        QToolButton *widget = qobject_cast<QToolButton*>(sender());
        if (widget)
        {
            QColor initColor = widget->palette().color(QPalette::Normal, QPalette::Window);
            QColor color = QColorDialog::getColor(initColor, this);
            if (color.isValid())
                setButtonColor(widget, color);
        }
    }
}
