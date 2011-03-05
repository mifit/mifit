#include <QAction>
#include <QDir>
#include <QFile>
#include <QFileDialog>
#include <QMenu>
#include <QMessageBox>
#include <QSettings>
#include "ui/Application.h"
#include "BatchJob.h"
#include "BatchJobManager.h"
#include "core/corelib.h"

using namespace std;

namespace
{
    bool testPyQt(QString pythonPath)
    {
        QProcess python;
        python.start(pythonPath, QStringList() << "-");
        if (python.waitForStarted())
        {
            python.write("from PyQt4 import QtCore\n");
            python.write("print 'PyQt', QtCore.PYQT_VERSION_STR\n");
            python.closeWriteChannel();
            if (python.waitForFinished())
            {
                QString versions = python.readAll();
                if (versions.contains("PyQt 4.7"))
                    return true;
            }
        }
        return false;
    }

    QString checkPyQtPath(QString dir)
    {
#ifdef Q_OS_WIN32
        const QString exe = "python.exe";
#else
        const QString exe = "python";
#endif
        QDir d(dir);
        if (d.exists(exe))
        {
            QString exePath = d.absoluteFilePath(exe);
            if (testPyQt(exePath))
                return exePath;
        }
        return QString::null;
    }

    QString findPyQt()
    {
        QString pythonExePath;
#ifdef Q_OS_WIN32
        const QString separator = ";";
#else
        const QString separator = ":";
#endif

        QString pathEnv = getenv("PATH");
        QStringList paths = pathEnv.split(separator);
        foreach (QString p, paths)
        {
            QString exePath = checkPyQtPath(p);
            if (!exePath.isEmpty())
                return exePath;
        }
        return QString::null;
    }

} // anonymous namespace

QString BatchJobManager::pythonExe()
{
    static QString pythonExePath;
    if (pythonExePath.isEmpty() || !QFile::exists(pythonExePath))
    {
        QSettings settings;
        pythonExePath = settings.value("pythonExe").toString();
    }
    if (pythonExePath.isEmpty() || !QFile::exists(pythonExePath))
    {
        pythonExePath = findPyQt();

#ifdef Q_OS_WIN32
        QString filters = "Programs (*.exe);;All files (*.*)";
#else
        const QString filters = "All files (*)";
#endif
        while (pythonExePath.isEmpty())
        {
            QString fileName = QFileDialog::getOpenFileName(NULL, "Select Python Executable",
                                                            "/", filters);
            if (!testPyQt(fileName))
            {
                if (QMessageBox::warning(0, "Invalid Python",
                                         "MIFit requires Python with PyQt4 installed\n"
                                         "  Python: http://python.org/\n"
                                         "  PyQt4: http://www.riverbankcomputing.co.uk/software/pyqt",
                                         QMessageBox::Retry | QMessageBox::Abort) != QMessageBox::Retry)
                    break;
            }
            else
                pythonExePath = fileName;
        }

        if (!pythonExePath.isEmpty() && QFile::exists(pythonExePath))
        {
            QSettings settings;
            settings.setValue("pythonExe", pythonExePath);
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

    QSettings settings;
    settings.beginWriteArray("customJobs");
    settings.setArrayIndex(_customJobIndex);
    settings.setValue("job", QVariant::fromValue(customJob));
    ++_customJobIndex;
    settings.endArray();

    return jobAction;
}

void BatchJobManager::setupJobMenu(QMenu *menu)
{
    QSettings settings;
    _customJobIndex = settings.beginReadArray("customJobs");
    for (int i = 0; i < static_cast<int>(_customJobIndex); ++i)
    {
        settings.setArrayIndex(i);
        CustomJob job = settings.value("job").value<CustomJob>();
        QAction *jobAction = new QAction(job.menuName, this);
        jobAction->setData(QVariant::fromValue(job));
        connect(jobAction, SIGNAL(triggered()), this, SLOT(handleCustomJobAction()));
        menu->addAction(jobAction);
    }
    settings.endArray();
}

void BatchJobManager::saveJobMenu(QMenu *menu)
{
    QSettings settings;
    settings.beginWriteArray("customJobs");
    int index = 0;
    foreach (QAction *action, menu->actions())
    {
        if (action->isSeparator())
            continue;
        if (!action->data().isValid() || !action->data().canConvert<CustomJob>())
            continue;

        settings.setArrayIndex(index);
        settings.setValue("job", action->data());
        ++index;
    }
    settings.endArray();
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

