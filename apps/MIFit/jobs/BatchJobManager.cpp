#include <nongui/nonguilib.h>
#include <QAction>
#include <QDir>
#include <QFile>
#include <QFileDialog>
#include <QMenu>
#include <QSettings>
#include "ui/Application.h"
#include "BatchJob.h"
#include "BatchJobManager.h"
#include "core/corelib.h"

using namespace std;


static QString pythonExe()
{
    static QString pythonExePath;
    if (pythonExePath.isEmpty() || !QFile::exists(pythonExePath))
    {
        QSettings *settings = MIGetQSettings();
        pythonExePath = settings->value("pythonExe").toString();
    }
    if (pythonExePath.isEmpty() || !QFile::exists(pythonExePath))
    {
#ifdef Q_OS_WIN32
        QString separator = ";";
        QString exe = "python.exe";
        QString filters = "Programs (*.exe);;All files (*.*)";
#else
        QString separator = ":";
        QString exe = "python";
        QString filters = "All files (*)";
#endif
        QString pathEnv = getenv("PATH");
        QStringList paths = pathEnv.split(separator);
        foreach (QString p, paths)
        {
            QDir dir(p);
            if (dir.exists(exe))
            {
                pythonExePath = dir.absoluteFilePath(exe);
                break;
            }
        }
        if (pythonExePath.isEmpty())
        {
            QString fileName = QFileDialog::getOpenFileName(NULL, "Select Python Executable",
                                                            "/", filters);
            if (!fileName.isEmpty())
                pythonExePath = fileName;
        }

        if (!pythonExePath.isEmpty() && QFile::exists(pythonExePath))
        {
            QSettings *settings = MIGetQSettings();
            settings->setValue("pythonExe", pythonExePath);
        }
    }

    return pythonExePath;
}

struct CustomJob
{
    QString menuName;
    QString jobName;
    QString executable;
    QStringList arguments;
    QString workingDirectory;

    CustomJob()
    {
    }

    CustomJob(const QString &menuName, const QString &jobName, const QString &executable, const QStringList &arguments, const QString &workingDirectory)
        : menuName(menuName),
          jobName(jobName),
          executable(executable),
          arguments(arguments),
          workingDirectory(workingDirectory)
    {
    }

    friend QDataStream &operator<<(QDataStream &out, const CustomJob &job);
    friend QDataStream &operator>>(QDataStream &in, CustomJob &job);
};

QDataStream &operator<<(QDataStream &out, const CustomJob &job)
{
    out << job.menuName << job.jobName << job.executable << job.arguments
            << job.workingDirectory;
    return out;
}

QDataStream &operator>>(QDataStream &in, CustomJob &job)
{
    in >> job.menuName >> job.jobName >> job.executable >> job.arguments
            >> job.workingDirectory;
    return in;
}

Q_DECLARE_METATYPE(CustomJob);


BatchJobManager::BatchJobManager()
{
    qRegisterMetaType<CustomJob>("CustomJob");
    qRegisterMetaTypeStreamOperators<CustomJob>("CustomJob");
}

BatchJobManager::~BatchJobManager()
{
    vector<BatchJob*>::iterator i, e;
    i = JobList.begin();
    e = JobList.end();
    for (; i != e; i++)
    {
        delete *i;
    }
}

QAction *BatchJobManager::customJobAction(const QString &menuName, const QString &jobName,
                                         const QString &executable, const QStringList &arguments,
                                         const QString &workingDirectory)
{
    CustomJob customJob(menuName, jobName, executable, arguments, workingDirectory);
    QAction *jobAction = new QAction(menuName, this);
    jobAction->setData(QVariant::fromValue(customJob));
    connect(jobAction, SIGNAL(triggered()), this, SLOT(handleCustomJobAction()));

    QSettings *settings = MIGetQSettings();
    settings->beginWriteArray("customJobs");
    settings->setArrayIndex(_customJobIndex);
    settings->setValue("job", QVariant::fromValue(customJob));
    ++_customJobIndex;
    settings->endArray();

    return jobAction;
}

void BatchJobManager::setupJobMenu(QMenu *menu)
{
    QSettings *settings = MIGetQSettings();
    _customJobIndex = settings->beginReadArray("customJobs");
    for (int i = 0; i < _customJobIndex; ++i)
    {
        settings->setArrayIndex(i);
        CustomJob job = settings->value("job").value<CustomJob>();
        QAction *jobAction = new QAction(job.menuName, this);
        jobAction->setData(QVariant::fromValue(job));
        connect(jobAction, SIGNAL(triggered()), this, SLOT(handleCustomJobAction()));
        menu->addAction(jobAction);
    }
    settings->endArray();
}

void BatchJobManager::saveJobMenu(QMenu *menu)
{
    QSettings *settings = MIGetQSettings();
    settings->beginWriteArray("customJobs");
    int index = 0;
    foreach (QAction *action, menu->actions())
    {
        if (action->isSeparator())
            continue;
        if (!action->data().isValid() || !action->data().canConvert<CustomJob>())
            continue;

        settings->setArrayIndex(index);
        settings->setValue("job", action->data());
        ++index;
    }
    settings->endArray();
}

void BatchJobManager::handleCustomJobAction()
{
    QAction *action = static_cast<QAction*>(sender());
    CustomJob customJob = action->data().value<CustomJob>();
    BatchJob *job = CreateJob();
    job->setJobName(customJob.jobName);
    if (customJob.executable == "python")
    {
        QString python = pythonExe();
        if (python.isEmpty())
            return;
        job->setProgram(python);
    }
    else
    {
        job->setProgram(customJob.executable);
    }
    job->setArguments(customJob.arguments);
    job->setWorkingDirectory(customJob.workingDirectory);
    job->StartJob();
}

BatchJob *BatchJobManager::CreateJob()
{
    BatchJob *job = new BatchJob;
    JobList.push_back(job);
    jobAdded(job);
    return job;
}

void BatchJobManager::CleanSucc()
{
    int i, size;
    BatchJob *job;
    vector<BatchJob*> tokill;
    size = JobList.size();
    for (i = 0; i < size; i++)
    {
        job = *(JobList.begin()+i);
        if (!job->isRunning() && job->isSuccess())
        {
            tokill.push_back(job);
        }
    }
    size = tokill.size();
    for (i = 0; i < size; i++)
    {
        DeleteJob(tokill[i]);
    }
}

void BatchJobManager::CleanAll()
{
    int i, size;
    BatchJob *job;
    vector<BatchJob*> tokill;
    size = JobList.size();
    for (i = 0; i < size; i++)
    {
        job = *(JobList.begin()+i);
        if (!job->isRunning())
        {
            tokill.push_back(job);
        }
    }
    size = tokill.size();
    for (i = 0; i < size; i++)
    {
        DeleteJob(tokill[i]);
    }
}

bool BatchJobManager::DeleteJob(BatchJob *job)
{
    if (!job)
        return false;
    if (job->isRunning())
    {
        Logger::message("You cannot delete a job while it is still running");
        return false;
    }
    for (size_t i = 0; i < JobList.size(); i++)
    {
        if (job == JobList[i])
        {
            delete *(JobList.begin()+i);
            JobList.erase(JobList.begin()+i);
            jobDeleted(job);
            return true;
        }
    }
    return false;
}

void BatchJobManager::ShowLogFile(BatchJob *b)
{
    if (b)
        b->ShowLog();
}

int BatchJobManager::numberOfRunningJobs()
{
    int count = 0;
    std::vector<BatchJob*>::iterator iter = JobList.begin();
    while (iter != JobList.end())
    {
        BatchJob *job = *iter;
        ++iter;
        if (job && job->isRunning())
        {
            ++count;
        }
    }
    return count;
}

