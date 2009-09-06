#ifndef mifit_wxdr_MIDataDialog_h
#define mifit_wxdr_MIDataDialog_h

#include <QDialog>
#include <QGridLayout>
#include <string>
#include <map>
#include "core/MIData.h"

class MIDataDialog : public QDialog {
  Q_OBJECT
  
  enum ControlType {
    TEXT,
    RADIO,
    INT,
    UNSIGNEDINT,
    SHORT,
    FLOAT,
    DOUBLE,
    BOOLEAN,
    COLOR,
    COLORINDEX
  };
    
  MIData* data;
  typedef std::map<std::string, std::pair<ControlType, QWidget*> > DataControlMap;
  DataControlMap dataControls;
  typedef std::vector<std::string> OrderList;
  OrderList order_;
  typedef std::map<std::string, std::string> LabelsMap;
  LabelsMap labels;

  std::map<QWidget*,unsigned char> colorIndexMap;

  QGridLayout* layout;

  void addControl(const std::string &key, const MIDatum &value, int row);
  
public:

  MIDataDialog(QWidget* parent = 0, Qt::WindowFlags f = 0); 
  void setMIData(MIData* data);

  void order(const std::string& key);
  void label(const std::string& key, const std::string& label);


public Q_SLOTS:
  void accepted();
  void colorButtonPressed();
  void colorIndexButtonPressed();
};

#endif

