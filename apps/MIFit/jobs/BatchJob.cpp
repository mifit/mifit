#include <QDir>
#include <QFile>
#include <QProcess>
#include <QDialog>
#include <QMessageBox>
#include <QTextBrowser>
#include <QDialogButtonBox>
#include <QVBoxLayout>

#include <nongui/nonguilib.h>
#include <util/utillib.h>
#include "ui/uilib.h"
#include "ui/MIDialog.h"

#include "BatchJob.h"
#include "OpenJobResults.h"

#ifdef _WIN32
#include <process.h>
#include <time.h>
#endif

using namespace std;

#ifdef _WIN32
#define strncasecmp strnicmp
#endif

BatchJob::BatchJob()
    : process(NULL)
{
  setJobId();
  setJobDir("");
}

BatchJob::BatchJob(const std::string& dir)
    : process(NULL)
{
  setJobId();
  setJobDir(dir.c_str());
}

void BatchJob::setJobDir(const char* dir) {
  jobDir = dir;
  if (jobDir.size() == 0) {
    jobDir = QDir::current().absolutePath();
  }
}

BatchJob::~BatchJob() {
    delete process;
}

void BatchJob::setJobId() {
  jobId_ = getpid()*1000;
  time_t t;
  time(&t);
  jobId_ += abs((int)(t%1000));
}

bool BatchJob::StartJob() {

    process = new QProcess(this);
    process->setWorkingDirectory(jobDir);
    process->setProcessChannelMode(QProcess::MergedChannels);
    process->setStandardOutputFile(LogFile);
    connect(process, SIGNAL(finished(int)),
            this, SLOT(doJobFinished()));
    connect(process, SIGNAL(finished(int)),
            this, SLOT(signalJobChanged()));

    process->start(program_, arguments_);
    process->closeWriteChannel();

    bool started = process->waitForStarted();
    jobChanged(this);
    return started;

}

void BatchJob::signalJobChanged() {
    jobChanged(this);
}

bool BatchJob::isRunning() {
    return process && (process->state() == QProcess::Running
                       || process->state() == QProcess::Starting);
}


void BatchJob::doJobFinished() {
std::string jobName;
  if (settings["jobName"].str != MIDatum::INVALID_STRING) {
    jobName = format("%s (%d)", settings["jobName"].str.c_str(), jobId_);
  } else {
    jobName = format("%d", jobId_);
  }
  QMessageBox::information(MIMainWindow::instance(), jobName.c_str(), "MIFit Job Finished");
}

void BatchJob::AbortJob() {
    process->kill();
}

std::string BatchJob::Info() {
  std::string s;
  if (settings["jobName"].str != MIDatum::INVALID_STRING) {
    s += "Job name: " + settings["jobName"].str + "\n"; 
  }
  s += format("Job id: %ld\n"
              "Log file: %s\n"
              "Job directory: %s\n"
              "Running: %d\n"
              "Success: %d\n",
              jobId_, LogFile.toAscii().constData(),
              jobDir.toAscii().constData(), (int)isRunning(), (int)isSuccess());
  if (settings["workingDirectory"].str != MIDatum::INVALID_STRING) {
    s += "Working directory: " + settings["workingDirectory"].str + "\n"; 
  }
  return s;
}

void BatchJob::ShowLog() {

    QDialog dlg(MIMainWindow::instance());
    dlg.setWindowTitle(LogFile);
    dlg.setModal(true);
    dlg.setSizeGripEnabled(true);
    QVBoxLayout* mainLayout = new QVBoxLayout;
    dlg.setLayout(mainLayout);

    QTextBrowser *browse = new QTextBrowser(&dlg);
    mainLayout->addWidget(browse);

    QDialogButtonBox* bb = new QDialogButtonBox(QDialogButtonBox::Ok, Qt::Horizontal, &dlg);
    mainLayout->addWidget(bb);
    dlg.connect(bb, SIGNAL(accepted()), &dlg, SLOT(accept()));

    QFile logFileObj(LogFile);
    if ( !logFileObj.open(QFile::ReadOnly | QFile::Text) ) {
        browse->setPlainText(QObject::tr("ERROR: Log file %1 not found!")
                             .arg(LogFile));
    } else {
        QTextStream logStream(&logFileObj);
        browse->setPlainText(logStream.readAll());
    }

    dlg.exec();
}

QString BatchJob::getJobDir() const {
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
  QString workDir;
  if (settings["workingDirectory"].str != MIDatum::INVALID_STRING) {
    workDir = settings["workingDirectory"].str.c_str();
  } else {
    workDir = jobDir;
  }
  std::string jobName;
  if (settings["jobName"].str != MIDatum::INVALID_STRING) {
    jobName = settings["jobName"].str;
  }
  OpenJobResults::prompt(workDir.toStdString(), jobName);
}

QStringList BatchJob::parseArgs(const QString &program)
{
    QStringList args;
    QString tmp;
    int quoteCount = 0;
    bool inQuote = false;

    // Tokens can be surrounded by double quotes "hello world".
    // Three consecutive double quotes represent
    // the quote character itself.
    for (int i = 0; i < program.size(); ++i) {
        if (program.at(i) == QLatin1Char('"')) {
            ++quoteCount;
            if (quoteCount == 3) {
                // third consecutive quote
                quoteCount = 0;
                tmp += program.at(i);
            }
            continue;
        }
        if (quoteCount) {
            if (quoteCount == 1)
                inQuote = !inQuote;
            quoteCount = 0;
        }
        if (!inQuote && program.at(i).isSpace()) {
            if (!tmp.isEmpty()) {
                args += tmp;
                tmp.clear();
            }
        } else {
            tmp += program.at(i);
        }
    }
    if (!tmp.isEmpty())
        args += tmp;

    return args;
}

void BatchJob::setCommandLine(const QString& command)
{
    arguments_ = parseArgs(command);
    program_ = arguments_.front();
    arguments_.pop_front();
}
