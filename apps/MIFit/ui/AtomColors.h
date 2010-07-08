#ifndef __AtomColors_H__
#define __AtomColors_H__

#include "core/corelib.h"
#include "ui_AtomColors.h"
#include <string>
#include <vector>

class AtomColors : public QDialog, public Ui::AtomColors
{
    Q_OBJECT

public:
    AtomColors(QWidget *parent = 0);
    void GetData(std::vector<std::string> &atomNames,
                 std::vector<std::string> &atomColors, bool &save);

private slots:
    void on_deleteTypePushButton_clicked();
    void on_editNamePushButton_clicked();
    void on_addTypePushButton_clicked();
    void on_resetPushButton_clicked();
    void on_colorToolButton_clicked();
    void on_listWidget_currentRowChanged(int);

private:
    void populateList(const std::vector<std::string> &atomNames,
                      const std::vector<std::string> &atomColors);
    void setColor(int i);

    std::vector<std::string> atomNames;
    std::vector<std::string> atomColors;
    std::vector<std::string> currentAtomNames;
    std::vector<std::string> currentAtomColors;
};

#endif // ifndef __AtomColors_H__
