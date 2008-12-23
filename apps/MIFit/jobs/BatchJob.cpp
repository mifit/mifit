#include <QDir>
#include <QFile>
#include <QProcess>
#include <QDialog>
#include <QTextBrowser>
#include <QDialogButtonBox>
#include <QVBoxLayout>

#include "nonguilib.h"
#include "utillib.h"
#include "uilib.h"
#include "MIDialog.h"

#include "BatchJob.h"
#include "OpenJobResults.h"

#ifdef _WIN32
#include <process.h>
#include <time.h>
#endif



//#include "LogView.h"

//NOTE: this could stand to be substantially revised to take better
//      advantage of Qt's QProcess to manage job state, but for now, we'll
//      just do the minimal port

using namespace std;

#ifdef _WIN32
#define strncasecmp strnicmp
#endif

BatchJob::BatchJob()
  : Cleaned(false), pid(0), running(false), completed(false), success(false), 
     CommandFileFp(NULL), m_doc(NULL) {

  setJobId();
  setJobDir("");
}

BatchJob::BatchJob(const std::string& dir)
  : Cleaned(false), pid(0), running(false), completed(false), success(false), 
     CommandFileFp(NULL), m_doc(NULL) {

  setJobId();

  setJobDir(dir.c_str());

  openCommandFile();
}

void BatchJob::openCommandFile() {
  std::string savedir = QDir::current().absolutePath().toStdString();
  QDir::setCurrent(jobDir.c_str());

  UpdateScript = format("mifit%lu_update.mlw", JobId);
  FinishedFile = format("mifit%lu_finished.txt", JobId);
#ifndef _WIN32
  CommandFile = format("mifit%lu_commands.sh", JobId);
#else
  CommandFile = format("mifit%lu_commands.bat", JobId);
#endif
  CommandFileFp = fopen(CommandFile.c_str(), "w");
  if (CommandFileFp == NULL) {
    QDir::setCurrent(savedir.c_str());
    return;
  }
#ifndef _WIN32
  fprintf(CommandFileFp, "#!/bin/sh\n");
#endif
  fprintf(CommandFileFp, "cd \"%s\"\n", QDir::toNativeSeparators(jobDir.c_str()).toStdString().c_str());

  QDir::setCurrent(savedir.c_str());
}

void BatchJob::setJobDir(const char* dir) {
  jobDir = dir;
  if (jobDir.size() == 0) {
    jobDir = QDir::current().absolutePath().toStdString();
  }
}

void BatchJob::CleanUp() {
  std::string savedir = QDir::current().absolutePath().toStdString();
  QDir::setCurrent(jobDir.c_str());

  if (CommandFileFp) {
    fclose(CommandFileFp);
    CommandFileFp = NULL;
  }
  QFile::remove(UpdateScript.c_str());
  QFile::remove(CommandFile.c_str());
  QFile::remove(FinishedFile.c_str());
  if (CleanupList.size() > 0) {
    for (unsigned int i = 0; i < CleanupList.size(); i++) {
      QFile::remove(CleanupList[i].c_str());
    }
  }
  Cleaned = true;

  QDir::setCurrent(savedir.c_str());
}

BatchJob::~BatchJob() {
  if (!Cleaned) {
    CleanUp();
  }
}

void BatchJob::setJobId() {
  JobId = getpid()*1000;
  time_t t;
  time(&t);
  JobId += abs((int)(t%1000));
}

bool BatchJob::StartJob() {
  std::string savedir = QDir::current().absolutePath().toStdString();
  QDir::setCurrent(jobDir.c_str());

  if (CommandFileFp) {
#ifndef _WIN32
    fprintf(CommandFileFp,
      "if [ $? -eq  0 ]\n"
      "then\n"
      "echo SUCCESS > %s\n"
      "else\n"
      "echo FAILED > %s\n"
      "fi\n\n",
      FinishedFile.c_str(), FinishedFile.c_str());
#else
    fprintf(CommandFileFp, "\nif %%ERRORLEVEL%% EQU 0 (\n"
      "echo SUCCESS > %s\n) else (\n"
      "echo FAILED > %s\n)\n\n", FinishedFile.c_str(), FinishedFile.c_str());
#endif
    fclose(CommandFileFp);
    CommandFileFp = NULL;
#ifndef _WIN32
    std::string s;
    QStringList args;
    args.append(CommandFile.c_str());
    pid = QProcess::startDetached("/bin/sh",args);
    if (pid == 0) {
      running = false;
      setSuccess(false);
      Logger::log("Command failed!");
      QDir::setCurrent(savedir.c_str());
      return false;
    }
#else
    pid = QProcess::startDetached(CommandFile.c_str());
    if (pid == 0) {
      running = false;
      setSuccess(false);
      Logger::log("Command failed!");
      QDir::setCurrent(savedir.c_str());
      return false;
    }
#endif

    running = true;
    jobChanged(this);
  } else {
    running = false;
    setSuccess(false);
  }
  QDir::setCurrent(savedir.c_str());
  return running;
}

bool BatchJob::IsRunning() {
  if (running) {
    HasEnded();
  }
  return running;
}

bool BatchJob::HasEnded() {
  if (!running) {
    return true;
  }

  std::string savedir = QDir::current().absolutePath().toStdString();
  QDir::setCurrent(jobDir.c_str());

  FILE* fp = fopen(FinishedFile.c_str(), "r");
  char buf[100];
  if (fp) {
    // we should make this more sophisticated!
    fgets(buf, sizeof buf, fp);
    if (strlen(buf) < 6) {
      // if too small then file still not closed - try again later
      fclose(fp);
      QDir::setCurrent(savedir.c_str());
      return false;
    }
    running = false;
    if (strncasecmp(buf, "SUCCESS", 7) == 0) {
      setSuccess(true);
    } else {
      setSuccess(false);
    }
    fclose(fp);

    doJobFinished();
    QDir::setCurrent(savedir.c_str());
    return true;
  }
  QDir::setCurrent(savedir.c_str());
  return false;
}

void BatchJob::doJobFinished() {
  std::string savedir = QDir::current().absolutePath().toStdString();
  QDir::setCurrent(jobDir.c_str());

  // run the update script
  if (m_doc && success) {
    m_doc->LoadScript(UpdateScript.c_str());
  }
  std::string jobName;
  if (settings["jobName"].str != MIDatum::INVALID_STRING) {
    jobName = format("%s (%d)", settings["jobName"].str.c_str(), JobId);
  } else {
    jobName = format("%d", JobId);
  }
  MIMessageBox(jobName, "MIFit Job Finished");
  if (!m_doc && success) {
    bool openResultsOnJobFinished = false;
    MIConfig::Instance()->Read("openResultsOnJobFinished", &openResultsOnJobFinished, false);
    if (openResultsOnJobFinished) {
      openResults();
    }
  }
  QDir::setCurrent(savedir.c_str());
}

bool BatchJob::WriteCommand(const std::string& buf) {
  std::string command;
  command += buf;
  if (CommandFileFp) {
    fprintf(CommandFileFp, "%s", command.c_str());
    if (command[command.size()-1] != '\n') {
      fprintf(CommandFileFp, "\n");
    }
    return true;
  } else {
    return false;
  }
}

bool BatchJob::AbortJob() {
  return false;
}

void BatchJob::AddtoCleanup(const std::string& file) {
  //canonicalize file name for QFile first.

  CleanupList.push_back(file);
}

std::string BatchJob::Info() {
  std::string s;
  if (settings["jobName"].str != MIDatum::INVALID_STRING) {
    s += "Job name: " + settings["jobName"].str + "\n"; 
  }
  s += format("Job id: %ld\n"
              "Log file: %s\n"
              "Command file: %s\n"
              "Finished file: %s\n"
              "Update script: %s\n"
              "Job directory: %s\n"
              "Running: %d\n"
              "Success: %d\n",
              JobId, LogFile.c_str(), CommandFile.c_str(), FinishedFile.c_str(),
              UpdateScript.c_str(), jobDir.c_str(), (int)running, (int)success);
  if (settings["workingDirectory"].str != MIDatum::INVALID_STRING) {
    s += "Working directory: " + settings["workingDirectory"].str + "\n"; 
  }
  return s;
}

void BatchJob::ShowLog() {
  std::string savedir = QDir::current().absolutePath().toStdString();
  QDir::setCurrent(jobDir.c_str());

  QDialog dlg(m_doc);
  dlg.setWindowTitle(LogFile.c_str());
  dlg.setModal(true);
  dlg.setSizeGripEnabled(true);

  QTextBrowser *browse = new QTextBrowser(&dlg);
  browse->setSource(QUrl::fromLocalFile(LogFile.c_str()));
  QDialogButtonBox *bb=new QDialogButtonBox(QDialogButtonBox::Ok, Qt::Horizontal, &dlg);
  dlg.connect(bb, SIGNAL(accepted()), &dlg, SLOT(accept()));
  
  QVBoxLayout *mainLayout = new QVBoxLayout;
  mainLayout->addWidget(browse);
  mainLayout->addWidget(bb);
  dlg.setLayout(mainLayout);
  dlg.exec();

  QDir::setCurrent(savedir.c_str());
}

void BatchJob::SetDocument(MIGLWidget* doc) {
  m_doc = doc;
}

MIGLWidget* BatchJob::GetDocument() {
  return m_doc;
}

void BatchJob::setSuccess(bool value) {
  completed = true;
  success = value;
  jobChanged(this);
}

std::string BatchJob::getJobDir() const {
  return jobDir;
}

void BatchJob::setSettings(const MIData& jobSettings) {
  settings = jobSettings;
  jobChanged(this);
}

MIData& BatchJob::getSettings() {
  return settings;
}

void BatchJob::openResults() {
  std::string workDir;
  if (settings["workingDirectory"].str != MIDatum::INVALID_STRING) {
    workDir = settings["workingDirectory"].str; 
  } else {
    workDir = jobDir;
  }
  std::string jobName;
  if (settings["jobName"].str != MIDatum::INVALID_STRING) {
    jobName = settings["jobName"].str;
  }
  OpenJobResults::prompt(workDir, jobName);
}
