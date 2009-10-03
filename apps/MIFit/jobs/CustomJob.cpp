#include <QDir>


#include <nongui/nonguilib.h>


#include <util/utillib.h>

#include "ui/uilib.h"
#include "CustomJob.h"
#include "ui/MIMainWindow.h"
#include "ui/MIGLWidget.h"

#include <vector>



CustomJob::CustomJob() {
}

std::string CustomJob::Info() {
  return BatchJob::Info();
}


bool CustomJob::StartJob() {
  typedef std::map<std::string, std::string> SubstitutionMap;
  SubstitutionMap subs;
  subs["\\$DATA"] = settings["dataFile"].str;

  if (settings["useCurrentModel"].b) {
    MIGLWidget *doc = MIMainWindow::instance()->currentMIGLWidget();
    if (doc != NULL) {
      Molecule* model = doc->GetDisplaylist()->GetCurrentModel();
      if (model) {
        std::string modelFile = format("%s%cmifit_%ld.pdb", settings["workingDirectory"].str.c_str(), QDir::separator().toAscii(), jobId_);
        model->SavePDBFile(modelFile.c_str());
        subs["\\$MODEL"] = modelFile.c_str();
      }
    }
  } else {
    subs["\\$MODEL"] = settings["modelFile"].str;
  }

  std::string args = settings["arguments"].str;
  SubstitutionMap::iterator iter = subs.begin();
  while (iter != subs.end()) {
    std::string pattern = iter->first;
    QString second(iter->second.c_str());
    second.replace("\\\\","\\\\\\\\");
    QString qargs(args.c_str());
    qargs.replace(pattern.c_str(),iter->second.c_str());
    args=qargs.toStdString();
    ++iter;
  }

  setJobDir(settings["workingDirectory"].str.c_str());
  LogFile = QString("%1%2%3_%4.log").arg(jobDir).arg(QDir::separator())
            .arg(settings["jobName"].str.c_str()).arg(jobId_);

  program_ = settings["executable"].str.c_str();
  arguments_ = parseArgs(args.c_str());

  return BatchJob::StartJob();
}

void CustomJob::doJobFinished() {
  BatchJob::doJobFinished();
}

