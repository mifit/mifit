#ifndef CUSTOMJOBPREFERENCES_H
#define CUSTOMJOBPREFERENCES_H

#include <QWidget>

namespace Ui {
    class CustomJobPreferences;
}

class CustomJobPreferences : public QWidget
{
    Q_OBJECT

public:
    explicit CustomJobPreferences(QWidget *parent = 0);
    ~CustomJobPreferences();

    void savePreferences();

private:
    Ui::CustomJobPreferences *ui;

private slots:
    void on_removeButton_clicked();
    void updateMenuList();
    void updateRemoveButton();

};

#endif // CUSTOMJOBPREFERENCES_H
