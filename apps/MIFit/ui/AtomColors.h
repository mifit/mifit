#ifndef __AtomColors_H__
#define __AtomColors_H__

#include "core/corelib.h"
#include "ui_AtomColors.h"

// variables: atomNames.strList, atomColors.strList,
class AtomColors : public QDialog, public Ui::AtomColors
{
    Q_OBJECT

public:
    AtomColors(QWidget *parent = 0);
    void InitializeFromData(const MIData &dat);
    void GetData(MIData &dat);

private Q_SLOTS:
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
    MIData data;
};

#endif
