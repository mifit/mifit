#ifndef mifit_jobs_BatchJob_h
#define mifit_jobs_BatchJob_h

#include <string>
#include <QObject>
#include <QProcess>
#include <QString>
#include <QStringList>
#include <core/MIData.h>

class BatchJob : public QObject
{
    Q_OBJECT

public:

    BatchJob();

    /**
     * @param working_directory - the directory to run the command in - if blank, uses current directory
     */
    BatchJob(const QString& working_directory);

    virtual ~BatchJob();

    void ShowLog();

    virtual QString Info();

    void setProgram(const QString& program) {
        program_ = program;
    }

    void setArguments(const QStringList& arguments) {
        arguments_ = arguments;
    }

    void setArguments(const QString& arguments);

    void setCommandLine(const QString& command);

    /**
     * Start the job running
     */
    virtual bool StartJob();

    /**
   * Aborts the job by sending a kill to the operating system
   */
    void AbortJob();

    /**
   * returns true if the job is still running
   */
    bool isRunning();

    /**
   * returns true if the job is a success
   */
    bool isSuccess() {
        return process && process->exitCode() == 0;
    }

    /**
   * A unique job id
   */
    unsigned long jobId() const {
        return jobId_;
    }

    QString workingDirectory() const;
    void setWorkingDirectory(const QString& dir);

    QString jobName() const {
        return jobName_;
    }

    void setJobName(const QString& jobName) {
        jobName_ = jobName;
    }

    void openResults();

signals:
    void jobChanged(BatchJob*);


protected:
    void setJobId();

    unsigned long jobId_;
    QString jobName_;
    QString program_;
    QStringList arguments_;
    QString workingDirectory_;
    QString logFile;

    QProcess* process;

    static QStringList parseArgs(const QString &program);

protected slots:
    virtual void doJobFinished();
    void signalJobChanged();

};

#endif
