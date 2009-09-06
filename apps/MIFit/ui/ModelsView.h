#ifndef mifit_ui_ModelsView_h
#define mifit_ui_ModelsView_h

#include <map>
#include <vector>

#include "core/MIData.h"

#include "ViewSyncedPanel.h"

class QLineEdit;
class QPushButton;
class QToolButton;
class QResizeEvent;

class ModelsTree;
class ResiduesTree;
class AtomsTree;

class ModelsView : public ViewSyncedPanel {
  Q_OBJECT
  typedef std::map<QLineEdit*, ModelsTree*> LineEditToModelsTreeMap;
  static LineEditToModelsTreeMap lineEditToModelsTree;
  typedef std::map<QWidget*, QLineEdit*> PanelToLineEditMap;
  static PanelToLineEditMap panelToLineEdit;
  typedef std::map<QWidget*, QToolButton*> PanelToToolButtonMap;
  static PanelToToolButtonMap panelToToolButton;

  typedef std::map<QPushButton*, ModelsTree*> ButtonCtrlToModelsTreeMap;
  static ButtonCtrlToModelsTreeMap buttonCtrlToModelsTree;
  typedef std::map<QPushButton*, QLineEdit*> ButtonCtrlToLineEditMap;
  static ButtonCtrlToLineEditMap buttonCtrlToLineEdit;
  typedef std::map<QWidget*, QPushButton*> PanelToButtonCtrlMap;
  static PanelToButtonCtrlMap panelToButtonCtrl;

private slots:
  void OnGoToResidueReturnPressed();
  void OnSplitterChanged(int, int);
  void updateSyncView(bool state);

protected:

  void resizeEvent(QResizeEvent *event);

  virtual QWidget* createPanelForView(MIGLWidget* view, QWidget* parent);
  virtual void destroyContentsForView(MIGLWidget* view, QWidget* panel);

  ModelsTree* GetCurrentModelsTree() const;
  ResiduesTree* GetCurrentResiduesTree() const;
  AtomsTree* GetCurrentAtomsTree() const;

public:

  ModelsView(QWidget* parent);
  virtual ~ModelsView();

  bool HandleHistory(MIData& data);
  bool RandomTest();
};

#endif
