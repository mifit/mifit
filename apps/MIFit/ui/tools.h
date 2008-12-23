/*
 * Author: Paul Collins
 * Date: 8-17-2005
 * Description: This is just a seperate class to seperate the handling
 *		of the Tools menu. This is more for organizational purposes than anything
 *		else
 */
#ifndef MI_TOOLS_H
#define MI_TOOLS_H

#ifdef MI_USE_JOBS

#include <QObject>
#include "MIEventHandler.h"
#include "MIMenu.h"

class Tools : public QObject, public MIEventHandler {
Q_OBJECT
private:
  void CIFConvertlib(const char*);        // Generic function for running mi_convertlib
  static Tools *_instance;
  Tools();

public:     //Event handles for the tools menu
  static Tools *instance();

  static bool VerifyMIExpert();
  static bool VerifyCCP4();
  void FillToolsMenu(MIMenu*, bool havedoc = false);

private Q_SLOTS:
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
  void OnUpdateForJobLimit(const MIUpdateEvent &);
};

#endif // MI_USE_JOBS

#endif
