#include "OpenJobResults.h"

#include <util/utillib.h>
#include "ui/uilib.h"
#include "ui/MIMainWindow.h"

#include "ui/MIDialog.h"

#include <QFileInfo>
#include <QDateTime>
#include <QDir>


#include <vector>
#include <algorithm>


static bool later(const std::pair<std::string, QDateTime>& f1, const std::pair<std::string, QDateTime>& f2) {
  return f1 > f2;
}

static void sortFileNameByTime(std::vector<std::string>& names) {
  std::vector<std::pair<std::string, QDateTime> > fileTimes;
  std::vector<std::string>::iterator iter;
  for (iter = names.begin(); iter != names.end(); ++iter) {
    QFileInfo fn(iter->c_str());
    QDateTime time = fn.lastModified();
    fileTimes.push_back(std::make_pair(*iter, time));
  }
  std::sort(fileTimes.begin(), fileTimes.end(), later);

  names.clear();
  std::vector<std::pair<std::string, QDateTime> >::iterator iter2;
  for (iter2 = fileTimes.begin(); iter2 != fileTimes.end(); ++iter2) {
    names.push_back(iter2->first);
  }
}


void OpenJobResults::prompt(const std::string& workdir, const std::string& jobName) {
  QDir dir(workdir.c_str());
  if (!dir.exists()) {
    return;
  }

  std::vector<std::string> mlwFiles;
  std::vector<std::string> pdbFiles;
  std::vector<std::string> mtzFiles;

  QStringList fileList=dir.entryList(QDir::Files | QDir::Readable | QDir::NoDotAndDotDot);
  for (int i = 0; i < fileList.count(); ++i) {
    std::string file = fileList[i].toStdString();
    QFileInfo fi(file.c_str());

    if (fi.suffix() == QString("mlw")) {
      mlwFiles.push_back(file);
    } else if (fi.suffix() == QString("pdb")) {
      pdbFiles.push_back(file);
    } else if (fi.suffix() == QString("mtz")) {
      mtzFiles.push_back(file);
    }
  }

  sortFileNameByTime(mlwFiles);
  sortFileNameByTime(pdbFiles);
  sortFileNameByTime(mtzFiles);


  if (!(mlwFiles.size() > 0 || pdbFiles.size() > 0 || mtzFiles.size() > 0)) {
    std::string s = format("Sorry, unable to open any results because\nno session, PDB, or MTZ files were found in job working directory:\n%s", workdir.c_str());
    MIMessageBox(s.c_str(), "Open Results unsuccessful", MIDIALOG_ICON_WARNING);
    return;
  }


  std::string title=::format("Job %s Finished",jobName.size() ? jobName.c_str(): "");
  MIGenericDialog dlg(0,title);
  MIData data;

  if (mlwFiles.size() > 0) {
    mlwFiles.push_back("None");
    data["mlw"].radio = 0;
    data["mlw"].radio_count = mlwFiles.size();
    data["mlw"].radio_labels = mlwFiles;
    dlg.order("mlw");
    dlg.label("mlw","Load session:");
  }
  if (mtzFiles.size() > 0) {
    mtzFiles.push_back("None");
    data["mtz"].radio = 0;
    data["mtz"].radio_count = mtzFiles.size();
    data["mtz"].radio_labels = mtzFiles;
    dlg.order("mtz");
    dlg.label("mtz","Load data:");
  }
  if (pdbFiles.size() > 0) {
    pdbFiles.push_back("None");
    data["pdb"].radio = 0;
    data["pdb"].radio_count = pdbFiles.size();
    data["pdb"].radio_labels = pdbFiles;
    dlg.order("pdb");
    dlg.label("pdb","Load PDB:");
  }
  if (!dlg.GetResults(data))
    return;

  std::vector<std::string> files;
  bool newWin=true;
  if (mlwFiles.size()) {
    std::string res=data["mlw"].radio_labels[data["mlw"].radio];
    if ( res != std::string("None")) {
      files.push_back(res);
      newWin=false;
    }
  }
  if (pdbFiles.size()) {
    std::string res=data["pdb"].radio_labels[data["pdb"].radio];
    if ( res != std::string("None")) {
      files.push_back(res);
    }
  }
  if (mtzFiles.size()) {
    std::string res=data["mtz"].radio_labels[data["mtz"].radio];
    if (res != std::string("None")) {
      files.push_back(res);
    }
  }
  MIMainWindow::instance()->OpenFiles(files, newWin);
}
