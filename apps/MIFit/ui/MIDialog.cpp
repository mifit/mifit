#include <QMessageBox>
#include <QInputDialog>
#include <QColorDialog>
#include <QTextEdit>
#include <QDialogButtonBox>
#include <QLabel>
#include <QMdiArea>
#include <QFontDialog>
#include <QFileDialog>

#include <map/maplib.h>

#include "ui/MIMainWindow.h"

#include "MIDialog.h"

// local dialog class definitions
#include "MIDataDialog.h"
#include "MIColorPickerDlg.h"
#include "SelectCrystal.h"
#include "SmilesDialog.h"
#include "ContourOptions.h"
#include "RefinementOptionsDialog.h"
#include "PhaseFileLoadDialog.h"
#include "BValueColors.h"
#include "AtomColors.h"
#include "LSQFitDialog.h"
#include "IntegrateDialog.h"
#include "SadPhasing.h"
#include "NCSModeling.h"
#include "CocrystalSuperPos.h"
#include "MolRep.h"
#include "CustomJobDialog.h"
#include "JobReport.h"

MIDialog::MIDialog(QWidget* parent, const std::string& name) : _qparent(parent), _name(name) {
  // Note:
  //
  // we can't use MIMainWindow directly as a dialog parent, b/c if it is,
  // the active window changes to the dialog instead of the MIGLWidget,
  // (MIMainWindow gets subWindowActivated(0) signal), and that can bork
  // up handling of dialog results.  It appears to be safe to use the
  // MdiArea as the parent. FMH.

  if (!_qparent || _qparent == MIMainWindow::instance())
    _qparent=MIMainWindow::instance()->currentMIGLWidget();
  if (!_qparent)
    _qparent=MIMainWindow::instance()->getMdiArea();
}

MIDialog::~MIDialog() {
}

void ValidateData(const MIData& data)
{
  MIDataConstIter i=data.begin();
  for (; i!= data.end(); ++i) {
    MIDatum datum=i->second;
    if (datum.radio!=UINT_MAX && datum.radio_count==UINT_MAX) {
#ifdef DEBUG
      Logger::message("Programmer error: radio set, but radio_count not set!");
#endif
    }
  }
}

bool MIDialog::GetResults(MIData& data) {

  bool ret;

  ValidateData(data);
  ret = PromptForResults(data);

  return ret;
}

//
// Generic
//
MIGenericDialog::MIGenericDialog(QWidget* parent, const std::string& name)
  : MIDialog(parent, name) {
  _qdlg=new MIDataDialog(_qparent);
  _qdlg->setWindowTitle(name.c_str());
}


void MIGenericDialog::order(const std::string& key) {
  ((MIDataDialog*)_qdlg)->order(key);
}

void MIGenericDialog::label(const std::string& key, const std::string& label){
  ((MIDataDialog*)_qdlg)->label(key,label);
}

bool MIGenericDialog::PromptForResults(MIData& data) {
  ((MIDataDialog*)_qdlg)->setMIData(&data);
  if (_qdlg->exec() != QDialog::Accepted)
    return false;
  return true;
}

MIGenericDialog::~MIGenericDialog() {
  delete _qdlg;
}



int MIGetSingleChoiceIndex(const std::string& message,
                           const std::string& caption,
                           const std::vector<std::string>& choices,
                           QWidget *parent) {

  MIGenericDialog dialog(parent, caption);
  MIData data;
  data["choice"].radio = 0;
  data["choice"].radio_count = choices.size();
  data["choice"].radio_labels = choices;
  dialog.label("choice",message);
  if (!dialog.GetResults(data)) {
    return -1;
  }

  return data["choice"].radio;
}


//
// Message
//

MIMessageDialog::MIMessageDialog(QWidget* parent, const std::string& name, const std::string& message, int style)
  : MIDialog(parent, name), _message(message) {
  init(message, style);
}


void MIMessageDialog::init(const std::string& message, int style)
{
  _message=message;
  _style=style;
}

bool MIMessageDialog::PromptForResults(MIData& data) {
  int result=QMessageBox::Cancel;

  QString msg(_message.c_str());
#ifdef __APPLE__
  // on apple the title is ignored so we move it to the msg
  if (_name.size())
    msg = QString(_name.c_str()) + QString(":\n\n") + msg;
#endif

  if (_style & MIDIALOG_YES_NO) {
    QMessageBox::StandardButtons sb=QMessageBox::Yes | QMessageBox::No;
    QMessageBox::StandardButton deft=QMessageBox::Yes;

    if (_style & MIDIALOG_CANCEL)
      sb |= QMessageBox::Cancel;
    if (_style & MIDIALOG_NO_DEFAULT)
      deft=QMessageBox::No;

    result=QMessageBox::question(_qparent, _name.c_str(), msg, sb, deft);
  } else if (_style & MIDIALOG_ICON_WARNING) {
    result=QMessageBox::warning(_qparent, _name.c_str(), msg);

  } else if (_style & MIDIALOG_ICON_ERROR) {
    result=QMessageBox::critical(_qparent, _name.c_str(), msg);

  } else if (_style & MIDIALOG_ICON_INFORMATION) {
    result=QMessageBox::information(_qparent, _name.c_str(), msg);
  }

  switch (result) {
    case QMessageBox::No:
      data["val"].radio = 0;
      break;

    case QMessageBox::Ok:
    case QMessageBox::Yes:
      data["val"].radio = 1;
      break;

    default: // cancel
      data["val"].radio = 2;
      break;
  }
  return true;
}

int MIMessageBox(const std::string& message, const std::string& caption, int style, QWidget* parent) {
  std::string name = caption;
//  if (name.size() == 0) {
//    name = message;
//  }

  MIMessageDialog dlg(parent, name, message, style);
  MIData data;
  data["val"].radio = 0;
  data["val"].radio_count = 3;
  if (!dlg.GetResults(data)) {
    return MI_CANCEL;
  }

  switch (data["val"].radio) {
    case 0: return MI_NO; break;
    case 1: return MI_YES; break;
    case 2: return MI_CANCEL; break;
  }
  return 0;
}

bool MIYesNo(const std::string& message, const std::string& caption)
{
  return MIMessageBox(message, caption, MIDIALOG_YES_NO) == MI_YES;
}


//
// GetInteger
//
MIGetIntegerDialog::MIGetIntegerDialog(QWidget* parent, const std::string& name, const std::string& message)
  : MIDialog(parent, name), _message(message) {
  _min = INT_MIN;
  _max = INT_MAX;
}

bool MIGetIntegerDialog::GetValue(int deft, int& result, int min, int max) {
  _min=min;
  _max=max;

  MIData data;
  data["val"].i = deft;
  if (!GetResults(data)) {
    return false;
  }
  if (data["val"].i < _min || data["val"].i > _max) {
    return false;
  }
  result = data["val"].i;
  return true;
}

bool MIGetIntegerDialog::PromptForResults(MIData& data) {

  bool ok;
  int res=QInputDialog::getInteger(_qparent, _name.c_str(), _message.c_str(), data["val"].i,
                                   _min,_max,1,&ok);

  if (!ok) {
    return false;
  }
  data["val"].i = res;
  return true;
}

//
// GetUnsignedInteger
//
MIGetUnsignedIntegerDialog::MIGetUnsignedIntegerDialog(QWidget* parent, const std::string& name, const std::string& message)
  : MIDialog(parent, name), _message(message) {
  _max = UINT_MAX;
}

bool MIGetUnsignedIntegerDialog::GetValue(unsigned int deft, unsigned int& result, unsigned int min, unsigned int max) {
  _max=max;
  _min=min;

  MIData data;
  data["val"].u = deft;
  if (!GetResults(data)) {
    return false;
  }
  if (data["val"].u < _min || data["val"].u > _max) {
    return false;
  }
  result = data["val"].u;
  return true;
}

bool MIGetUnsignedIntegerDialog::PromptForResults(MIData& data) {

  bool ok;
  if (_max > INT_MAX)
    _max=INT_MAX;

  int res=QInputDialog::getInteger(_qparent, _name.c_str(), _message.c_str(), (int)(data["val"].u), (int)_min, (int)_max, 1, &ok);

  if (!ok) {
    return false;
  }
  data["val"].u = (unsigned int)res;
  return true;
}

//
// GetFloat
//
MIGetFloatDialog::MIGetFloatDialog(QWidget* parent, const std::string& name, const std::string& message)
  : MIDialog(parent, name), _message(message) {
  _min = FLT_MIN;
  _max = FLT_MAX;
}

bool MIGetFloatDialog::GetValue(float deft, float& result, float min, float max) {
  _max=max;
  _min=min;

  MIData data;
  data["val"].f = deft;
  if (!GetResults(data)) {
    return false;
  }
  if (data["val"].f < _min || data["val"].f > _max) {
    return false;
  }
  result = data["val"].f;
  return true;
}

bool MIGetFloatDialog::PromptForResults(MIData& data) {

  bool ok;
  double res=QInputDialog::getDouble(_qparent, _name.c_str(), _message.c_str(), (double)(data["val"].f), (double)_min, (double)_max, 3, &ok);

  if (!ok) {
    return false;
  }
  data["val"].f = (float)res;
  return true;
}

//
// GetString
//
MIGetStringDialog::MIGetStringDialog(QWidget* parent, const std::string& name, const std::string& message, bool multiline)
  : MIDialog(parent, name), _message(message),_multiline(multiline) {
}

bool MIGetStringDialog::GetValue(const std::string &deft, std::string &result) {
  MIData data;
  data["val"].str = deft;
  if (!GetResults(data)) {
    return false;
  }
  result = data["val"].str;
  return true;
}

bool MIGetStringDialog::PromptForResults(MIData& data) {

  if (_multiline) {
    QDialog dlg(_qparent);
    dlg.setWindowTitle(_name.c_str());
    dlg.setModal(true);
    dlg.setSizeGripEnabled(true);
    QLabel *msg=new QLabel(_message.c_str(),&dlg);
    QTextEdit *textEdit = new QTextEdit(data["val"].str.c_str(), &dlg);
    QDialogButtonBox *bb=new QDialogButtonBox(QDialogButtonBox::Ok | QDialogButtonBox::Cancel, Qt::Horizontal, &dlg);
    dlg.connect(bb, SIGNAL(accepted()), &dlg, SLOT(accept()));
    dlg.connect(bb, SIGNAL(rejected()), &dlg, SLOT(reject()));

    QVBoxLayout *mainLayout = new QVBoxLayout;
    mainLayout->addWidget(msg);
    mainLayout->addWidget(textEdit);
    mainLayout->addWidget(bb);
    dlg.setLayout(mainLayout);

    if (dlg.exec() != QDialog::Accepted) {
      return false;
    }
    data["val"].str = textEdit->toPlainText().toStdString();
    return true;
  }


  bool ok;
  QString res=QInputDialog::getText(_qparent, _name.c_str(), _message.c_str(), QLineEdit::Normal, data["val"].str.c_str(), &ok);

  if (!ok) {
    return false;
  }
  data["val"].str = res.toStdString();
  return true;
}


//
// ColorPrompt
//
MIColorPromptDialog::MIColorPromptDialog(QWidget* parent, const std::string& name)
  : MIDialog(parent, name) {
}


bool MIColorPromptDialog::PromptForResults(MIData& data) {
  QColor initColor=QColor((int)data["red"].u,(int)data["green"].u,(int)data["blue"].u);
  QColor color=QColorDialog::getColor(initColor, _qparent);
  if (!color.isValid()) {
    return false;
  }
  data["red"].u = (unsigned int)color.red();
  data["green"].u = (unsigned int)color.green();
  data["blue"].u = (unsigned int)color.blue();
  return true;
}


PaletteColor MIGetColorFromUser(QWidget* parent, const PaletteColor& deft) {
  MIColorPromptDialog dlg(parent, "Color dialog");
  MIData data;
  data["red"].u = (unsigned int)deft.red;
  data["green"].u = (unsigned int)deft.green;
  data["blue"].u = (unsigned int)deft.blue;
  if (!dlg.GetResults(data)) {
    return deft;
  }
  PaletteColor c;
  c.red=(unsigned char)data["red"].u;
  c.green=(unsigned char)data["green"].u;
  c.blue=(unsigned char)data["blue"].u;
  return c;

}

//
// ColorPalette
//
MIColorPaletteDialog::MIColorPaletteDialog(QWidget* parent, const std::string& name)
  : MIDialog(parent, name) {
  //GetColorFunc(_dlg, true, true);
}

bool MIColorPaletteDialog::PromptForResults(MIData& data) {
  MIColorPickerDlg dlg(_qparent,data["color"].radio);
  dlg.exec();
  data["color"].radio = dlg.GetResult();
  return true;
}

int MIColorChooser(int colorstart, const std::string name) {
  MIColorPaletteDialog dlg(0, name.c_str());
  MIData data;
  data["color"].radio = colorstart;
  data["color"].radio_count = Colors_NUMBERCOLORS;
  if (!dlg.GetResults(data)) {
    return -1;
  }
  return data["color"].radio;
}



//
// FontDialog
//
MIFontDialog::MIFontDialog(QWidget* parent, const std::string& name)
  : MIDialog(parent, name) {
}

bool MIFontDialog::PromptForResults(MIData &data)
{
  QFont f;
  f.fromString(data["fontDesc"].str.c_str());

  bool ok;
  f=QFontDialog::getFont(&ok,f,_qparent);
  if (!ok)
    return false;
  data["fontDesc"].str = f.toString().toStdString();
  return true;
}

bool MIGetFontFromUser(std::string &fontDesc) {
  MIData data;
  data["fontDesc"].str = fontDesc;
  MIFontDialog dlg(0,"Font selection");
  if (!dlg.GetResults(data)) {
    return false;
  }
  fontDesc=data["fontDesc"].str;
  return true;
}


//
// SelectDirectory
//
MISelectDirectoryDialog::MISelectDirectoryDialog(QWidget* parent,
                                                 const std::string& name,
                                                 const std::string& deft)
  : MIDialog(parent, name),_deft(deft) {
}

bool MISelectDirectoryDialog::PromptForResults(MIData& data) {
  QString str=QFileDialog::getExistingDirectory(_qparent, _name.c_str(),_deft.c_str());

  if (str.isEmpty()) {
    return false;
  }
  data["dir"].str = str.toStdString();
  return true;
}


//
// File
//
MIFileDialog::MIFileDialog(QWidget* parent, const std::string& message,
                           const std::string& deftDir,
                           const std::string& deftFile,
                           const std::string& filter,
                           unsigned int mode)
  : MIDialog(parent, message), _deftDir(deftDir), _deftFile(deftFile),
    _filter(filter), _mode(mode) {
}

static void stringSplit(std::string str,
                        const std::string &delim,
                        std::vector<std::string>& results) {
  unsigned int cutAt;
  while ( (cutAt = str.find_first_of(delim)) != str.npos) {
    if (cutAt > 0) {
      results.push_back(str.substr(0, cutAt));
    }
    str = str.substr(cutAt+1);
  }
  if (str.length() > 0) {
    results.push_back(str);
  }
}


bool MIFileDialog::PromptForResults(MIData& data) {
  std::string path=_deftDir;
  if (_deftFile.size()) {
    path = path+ "/" + _deftFile;
  }
  if (data["path"].str.size() && data["path"].str != MIDatum::INVALID_STRING) {
    path=data["path"].str;
  }

  // build filter string/ vector;
  QString filter;
  std::vector<std::string> filterList;
  MIStringReplace(_filter, ",", " ");
  if (_filter.size()) {
    std::vector<std::string> results;
    stringSplit(_filter,"|",results);
    for (unsigned int i=0;i<results.size(); i += 2) {
      filterList.push_back(results[i]);
      if (i > 0)
        filter += ";;";
      filter += QString(results[i].c_str());
    }
  }

 	QString selectedFilter;
  QString response;
  QStringList responseList;

  std::vector<std::string>& pathlist = data["pathList"].strList;
  pathlist.clear();

  path = Application::instance()->latestFileBrowseDirectory(path.c_str()).toStdString();

  if (_mode == MI_SAVE_MODE) {
    response=QFileDialog::getSaveFileName(_qparent, _name.c_str(), path.c_str(), filter, &selectedFilter);
    if (response.isEmpty())
      return false;
    data["path"].str = response.toStdString();
    pathlist.push_back(response.toStdString());
  } else if (_mode==MI_OPEN_MODE) {
    response=QFileDialog::getOpenFileName(_qparent, _name.c_str(), path.c_str(), filter, &selectedFilter);
    if (response.isEmpty())
      return false;
    data["path"].str = response.toStdString();
    pathlist.push_back(response.toStdString());
  }

  QFileInfo fileInfo(data["path"].str.c_str());
  Application::instance()->latestFileBrowseDirectory(fileInfo.absolutePath());

  // set selected filter index
  std::string selFilter=selectedFilter.toStdString();
  for (size_t i=0; i<filterList.size(); ++i)
  {
    if (filterList[i] == selFilter) {
      data["filterIndex"].radio = i;
      break;
    }
  }

  return true;
}


//
// wraps MIFileDialog
//
std::string MIFileSelector(const std::string& title,
                           const std::string& defaultDirString,
                           const std::string& defaultFileNameString,
                           const std::string& defaultExtension,
                           const std::string& filter,
                           unsigned int flags,
                           QWidget* parent) {

    std::string dir = Application::instance()->latestFileBrowseDirectory(defaultDirString.c_str()).toStdString();

  std::string filter2 = filter;
  if (defaultExtension.size() && !filter.size() ) {
    filter2 = std::string("*.") + defaultExtension;
  }

  MIFileDialog fileDialog(parent, title, dir,
                          defaultFileNameString, filter2,
                          flags);

  MIData data;
  data["path"].str = "";
  data["pathList"].strList=std::vector<std::string>();

  // if filter is of form "All files (*)|*|..." set correct filter index
//   if (defaultExtension.size() != 0 && filter2.find_first_of('|') != std::string::npos) {
//     int filterIndex = 0;

//     wxArrayString descriptions, filters;
//     // don't care about errors, handled already by wxFileDialog
//     (void)wxParseCommonDialogsFilter(filter2.c_str(), descriptions, filters);
//     for (size_t n = 0; n < filters.GetCount(); n++) {
//       if (filters[n].Contains(defaultExtension.c_str())) {
//         filterIndex = n;
//         break;
//       }
//     }

//     if (filterIndex > 0) {
//       data["filterIndex"].radio = filterIndex;
//       data["filterIndex"].radio_count = filters.GetCount();
//     }
//   }

  std::string filename;
  if (fileDialog.GetResults(data)) {
    filename = data["path"].str;
  }
  QFileInfo fileInfo(filename.c_str());
  Application::instance()->latestFileBrowseDirectory(fileInfo.absolutePath());
  return filename;
}


//
// SelectCrystal
//
MISelectCrystalDialog::MISelectCrystalDialog(QWidget* parent, const std::string& name)
  : MIDialog(parent, name) {
}

bool MISelectCrystalDialog::PromptForResults(MIData& data) {
  SelectCrystal dlg(data["info"].str,_qparent);
  if (dlg.exec() != QDialog::Accepted) {
    return false;
  }
  data["info"].str = dlg.getLabel();
  return true;
}


//
// SmilesImport
//
MISmilesImportDialog::MISmilesImportDialog(QWidget* parent, const std::string& name)
  : MIDialog(parent, name) {
}

bool MISmilesImportDialog::PromptForResults(MIData& data) {
  SmilesDialog dlg(_qparent);
  if (dlg.exec() != QDialog::Accepted) {
    return false;
  }
  dlg.GetResults(data);
  return true;
}


//
// Contour
//
MIContourDialog::MIContourDialog(QWidget* parent, const std::string& name)
  : MIDialog(parent, name) {
}

bool MIContourDialog::PromptForResults(MIData& data) {
  QDialog dlg(_qparent);

  dlg.setWindowTitle(_name.c_str());
  dlg.setModal(true);
  dlg.setSizeGripEnabled(true);

  ContourOptions *co=new ContourOptions(&dlg,false);
  QDialogButtonBox *bb=new QDialogButtonBox(QDialogButtonBox::Ok | QDialogButtonBox::Cancel, Qt::Horizontal, &dlg);
  dlg.connect(bb, SIGNAL(accepted()), &dlg, SLOT(accept()));
  dlg.connect(bb, SIGNAL(rejected()), &dlg, SLOT(reject()));

  QVBoxLayout *mainLayout = new QVBoxLayout;
  mainLayout->addWidget(co);
  mainLayout->addWidget(bb);
  dlg.setLayout(mainLayout);

  co->InitializeFromData(data);
  if (dlg.exec() != QDialog::Accepted) {
    return false;
  }
  co->GetData(data);
  return true;
}

//
// BValue Colors
//
MIBValueColorsDialog::MIBValueColorsDialog(QWidget* parent, const std::string& name)
  : MIDialog(parent, name) {
}

bool MIBValueColorsDialog::PromptForResults(MIData& data) {
  BValueColors dlg(_qparent);
  dlg.InitializeFromData(data);
  if (dlg.exec() != QDialog::Accepted) {
    return false;
  }
  dlg.GetData(data);
  return true;
}

//
// AtomColors
//
MIAtomColorsDialog::MIAtomColorsDialog(QWidget* parent, const std::string& name)
  : MIDialog(parent, name) {
}

bool MIAtomColorsDialog::PromptForResults(MIData& data) {
  AtomColors dlg(_qparent);
  dlg.InitializeFromData(data);
  if (dlg.exec() != QDialog::Accepted) {
    return false;
  }
  dlg.GetData(data);
  return true;
}

//
// LSQFitDialog
//
MILSQFitDialog::MILSQFitDialog(QWidget* parent, const std::string& name)
  : MIDialog(parent, name) {
}

bool MILSQFitDialog::PromptForResults(MIData& data) {
  static LSQFitDialog dlg(_qparent);
  dlg.InitializeFromData(data);
  if (dlg.exec() != QDialog::Accepted) {
    return false;
  }
  dlg.GetData(data);
  return true;
}

//
// IntegrateDialogDialog
//
MIIntegrateDialog::MIIntegrateDialog(QWidget* parent, const std::string& name)
  : MIDialog(parent, name) {
}

void MIIntegrateDialog::GetInitialData(MIData& data) {
  IntegrateDialog::GetInitialData(data);
}

bool MIIntegrateDialog::PromptForResults(MIData& data) {
  static IntegrateDialog dlg(_qparent);
  dlg.InitializeFromData(data);
  if (dlg.exec() != QDialog::Accepted) {
    return false;
  }
  return dlg.GetData(data);
}


//
// SadPhasingDialog
//
MISadPhasingDialog::MISadPhasingDialog(QWidget* parent, const std::string& name)
  : MIDialog(parent, name) {
}

void MISadPhasingDialog::GetInitialData(MIData& data) {
  SadPhasing::GetInitialData(data);
}

bool MISadPhasingDialog::PromptForResults(MIData& data) {
  static SadPhasing dlg(_qparent);
  dlg.InitializeFromData(data);
  if (dlg.exec() != QDialog::Accepted) {
    return false;
  }
  return dlg.GetData(data);
}


//
// NCSModelingDialog
//
MINCSModelingDialog::MINCSModelingDialog(QWidget* parent, const std::string& name)
  : MIDialog(parent, name) {
}

void MINCSModelingDialog::GetInitialData(MIData& data) {
  NCSModeling::GetInitialData(data);
}

bool MINCSModelingDialog::PromptForResults(MIData& data) {
  static NCSModeling dlg(_qparent);
  dlg.InitializeFromData(data);
  if (dlg.exec() != QDialog::Accepted) {
    return false;
  }
  return dlg.GetData(data);
}



//
// CocrystalSuperpositionDialog
//
MICocrystalSuperpositionDialog::MICocrystalSuperpositionDialog(QWidget* parent, const std::string& name)
  : MIDialog(parent, name) {
}

void MICocrystalSuperpositionDialog::GetInitialData(MIData& data) {
  CocrystalSuperPos::GetInitialData(data);
}

bool MICocrystalSuperpositionDialog::PromptForResults(MIData& data) {
  static CocrystalSuperPos dlg(_qparent);
  dlg.InitializeFromData(data);
  if (dlg.exec() != QDialog::Accepted) {
    return false;
  }
  return dlg.GetData(data);
}



//
// MolRepDialog
//
MIMolRepDialog::MIMolRepDialog(QWidget* parent, const std::string& name)
  : MIDialog(parent, name) {
}

void MIMolRepDialog::GetInitialData(MIData& data) {
  MolRep::GetInitialData(data);
}

bool MIMolRepDialog::PromptForResults(MIData& data) {
  static MolRep dlg(_qparent);
  dlg.InitializeFromData(data);
  if (dlg.exec() != QDialog::Accepted) {
    return false;
  }
  return dlg.GetData(data);
}


//
// JobReportDialog
//
MIJobReportDialog::MIJobReportDialog(QWidget* parent, const std::string& name)
  : MIDialog(parent, name) {
}

void MIJobReportDialog::GetInitialData(MIData& data) {
  JobReport::GetInitialData(data);
}

bool MIJobReportDialog::PromptForResults(MIData& data) {
  static JobReport dlg(_qparent);
  dlg.InitializeFromData(data);
  if (dlg.exec() != QDialog::Accepted) {
    return false;
  }
  return dlg.GetData(data);
}




