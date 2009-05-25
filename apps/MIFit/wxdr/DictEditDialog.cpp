#include "DictEditDialog.h"
#include "uilib.h"


#include <QGridLayout>
#include <QDialogButtonBox>
#include <QMenuBar>
#include <QFrame>

DictEditDialog::DictEditDialog(QWidget* parent) : QDialog(parent) {

  int row=0;
  
  QGridLayout *_layout = new QGridLayout(this);
  
  _menuBar=new QMenuBar(this);
#ifndef __APPLE__
  _layout->addWidget(_menuBar,row++,0,1,2);
#endif

  _frame=new QFrame(this);
  _layout->addWidget(_frame,row,0);

  canvas=new DictEditCanvas(this);
  _layout->addWidget(canvas,row++,1);

  QDialogButtonBox *buttonBox = 
    new QDialogButtonBox(QDialogButtonBox::Ok |
                         QDialogButtonBox::Cancel,
                         Qt::Horizontal,
                         this);
  _layout->addWidget(buttonBox,row,0,1,2);
  
  connect(buttonBox, SIGNAL(accepted()), this, SLOT(accept()));
  connect(buttonBox, SIGNAL(rejected()), this, SLOT(reject()));
  
  setModal(true);
  setAttribute(Qt::WA_DeleteOnClose, true);
}

QWidget *DictEditDialog::getFrame() { return _frame; }
QMenuBar *DictEditDialog::getMenuBar() { return _menuBar; }


void DictEditDialog::accept() {
  // do something here then call
  // the underlying dialog box OnOK function
  canvas->OnOk();

  MIFitGeomRefiner()->EditEntryCleanup(true);
  MIData values;
  values["command"].str = "dict";
  values["state"].str = "ok";
  MIGetHistory()->AddCommand(values);

  QDialog::accept();
}

void DictEditDialog::reject() {

  MIFitGeomRefiner()->EditEntryCleanup(false);
  MIData values;
  values["command"].str = "dict";
  values["state"].str = "cancel";
  MIGetHistory()->AddCommand(values);

  QDialog::reject();
}
