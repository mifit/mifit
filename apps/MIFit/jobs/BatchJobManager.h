#ifndef mifit_jobs_BatchJobManager_h
#define mifit_jobs_BatchJobManager_h

#include <QObject>
#include <QTimer>
#include <string>
#include <vector>

class BatchJob;
class QAction;
class QMenu;

/**
 * Runs a batch job in the background.
 * It can start a job but the user has to pick the output up manually
 * and there is no way to track jobs currently.
 */
class BatchJobManager : public QObject
{
    Q_OBJECT

    std::vector<BatchJob*> JobList;
    uint _customJobIndex;

public:
    void ShowLogFile(BatchJob *b);
    /**
     * Create a new job and return a pointer to that job.
     *  Owner is usually this.
     */
    BatchJob *CreateJob();

    /**
     *  Constructor - should be just one copy per program.
     */
    BatchJobManager();
    ~BatchJobManager();

    QAction *customJobAction(const QString &menuName, const QString &jobName,
                             const QString &executable, const QStringList &arguments,
                             const QString &workingDirectory);

    /**
     * Returns a pointer to the ith job in the list (indexed from 0).
     */
    BatchJob *GetJob(int i)
    {
        return JobList[i];
    }

    /**
     * Deletes a job in the list given its pointer.
     * returns true if it finds the job in the lsit and destroys it, else false.
     */
    bool DeleteJob(BatchJob *p_job);

    /**
     * Returns the current JobList
     */
    std::vector<BatchJob*> *GetJobList()
    {
        return &JobList;
    }

    /**
     * Clean the Successfully completed jobs from the TreeMenu
     */
    void CleanSucc();

    /**
     * Clean ALL the completed jobs from the TreeMenu
     */
    void CleanAll();

    int numberOfRunningJobs();

    void setupJobsMenu(QMenu *menu);

signals:
    void jobAdded(BatchJob*);
    void jobDeleted(BatchJob*);

private slots:
    void handleCustomJobAction();

};

#endif // ifndef mifit_jobs_BatchJobManager_h
