#ifndef MI_TOOLS_H
#define MI_TOOLS_H

#include <QObject>
#include <QList>

class QAction;
class QMenu;

class Tools : public QObject {
Q_OBJECT
private:
  void CIFConvertlib(const char*);        // Generic function for running mi_convertlib

  Tools();

  QList<QAction*> actions;
  QList<QAction*> docActions;

public:     //Event handles for the tools menu
  static Tools& instance();

  static bool VerifyMIExpert();
  static bool VerifyCCP4();
  void FillToolsMenu(QMenu*);

private slots:
  void OnBindNGrind();
  void OnRunTestJob();

  void OnCIF2Shellx();
  void OnCIF2CNS();
  void OnMolRep();
  void OnRefmacRestraints();
  void OnRefine();
  void OnJobReport();
  void OnCoCrystalSuperPos();
  void OnSadPhasing();
  void OnNCSModeling();
  void OnCustom();
  void OnIntegrate();
  void OnUpdateForJobLimit();
};

#endif
