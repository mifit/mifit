#ifndef mifit_ui_GenericDataDialog_h
#define mifit_ui_GenericDataDialog_h

#include <QList>
#include <QDialog>
#include <QVariant>

class QFormLayout;

class GenericDataDialog : public QDialog
{
    Q_OBJECT

public:

    GenericDataDialog(QWidget *parent = 0, Qt::WindowFlags f = 0);

    void addBoolField(const QString &label, bool value);
    void addIntField(const QString &label, int value);
    void addUIntField(const QString &label, unsigned int value);
    void addDoubleField(const QString &label, double value);
    void addStringField(const QString &label, const QString &value);
    void addComboField(const QString &label, const QStringList &choices, int index);
    void addColorField(const QString &label, const QColor &color);
    void addColorIndexField(const QString &label, int colorIndex);

    QVariant value(int index) const;

public slots:
    virtual void accept();

private:
    enum FieldType
    {
        Bool,
        Int,
        UInt,
        Double,
        String,
        Combo,
        Color,
        ColorIndex
    };

    QFormLayout *formLayout_;
    QList<FieldType> fieldTypes_;
    QList<QVariant> values_;

private slots:
    void colorButtonPressed();
    void colorIndexButtonPressed();

};

#endif // ifndef mifit_ui_GenericDataDialog_h

