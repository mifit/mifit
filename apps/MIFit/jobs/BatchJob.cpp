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
    : workingDirectory_(QDir::current().absolutePath()), process(NULL)
{
  setJobId();
}

BatchJob::BatchJob(const QString& dir)
    : workingDirectory_(dir), process(NULL)
{
  setJobId();
}

void BatchJob::setWorkingDirectory(const QString& dir) {
  workingDirectory_ = dir;
  if (workingDirectory_.isEmpty()) {
    workingDirectory_ = QDir::current().absolutePath();
  }
}

BatchJob::~BatchJob() {
    if (QFile::exists(LogFile)) {
        QFile::remove(LogFile);
    }
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
    process->setWorkingDirectory(workingDirectory_);
    process->setProcessChannelMode(QProcess::MergedChannels);
    process->setStandardOutputFile(LogFile);
    connect(process, SIGNAL(finished(int)),
            this, SLOT(doJobFinished()));
    connect(process, SIGNAL(finished(int)),
            this, SLOT(signalJobChanged()));

    process->start(program_, arguments_);
    process->closeWriteChannel();

    bool started = process->waitForStarted();
    if (!started) {
        QMessageBox::warning(NULL, "Job Error", QString("Unable to start job %1").arg(jobId_));
        delete process;
        process = NULL;
    }
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
  QString message;
  if (!jobName_.isEmpty())
      message = QString("%1 finished (job %2)").arg(jobName_).arg(jobId_);
  else
      message = QString("Job %2 finished").arg(jobName_).arg(jobId_);
  QMessageBox::information(MIMainWindow::instance(), "MIFit Job Finished", message);
}

void BatchJob::AbortJob() {
    process->kill();
}

std::string BatchJob::Info() {
  return format("Job name: %s\n"
              "Job id: %ld\n"
              "Program: %s\n"
              "Arguments: \"%s\"\n"
              "Log file: %s\n"
              "Job directory: %s\n"
              "Running: %s\n"
              "Success: %s\n",
              jobName_.toAscii().constData(),
              jobId_,
              program_.toAscii().constData(),
              arguments_.join("\" \"").toAscii().constData(),
              LogFile.toAscii().constData(),
              workingDirectory_.toAscii().constData(),
              isRunning() ? "true" : "false",
              isSuccess() ? "true" : "false");
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

QString BatchJob::workingDirectory() const {
  return workingDirectory_;
}

void BatchJob::openResults() {
  OpenJobResults::prompt(workingDirectory_.toStdString(), jobName_.toStdString());
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

void BatchJob::setArguments(const QString& arguments)
{
    arguments_ = parseArgs(arguments);
}

void BatchJob::setCommandLine(const QString& command)
{
    arguments_ = parseArgs(command);
    program_ = arguments_.front();
    arguments_.pop_front();
}
