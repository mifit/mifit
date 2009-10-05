#ifndef mifit_jobs_BatchJobManager_h
#define mifit_jobs_BatchJobManager_h

#include <QObject>
#include <QTimer>
#include <string>
#include <vector>

class BatchJob;
class CustomJob;

/**
 * Runs a batch job in the background.
 * It can start a job but the user has to pick the output up manually
 * and there is no way to track jobs currently.
 */
class BatchJobManager : public QObject {
  Q_OBJECT

  std::vector<BatchJob*> JobList;
  std::string workdir;

public:
  void ShowLogFile(BatchJob* b);
  /**
   * Create a new job and return a pointer to that job.
   *  Owner is usually this.
   */
  BatchJob* CreateJob();

  /**
   *  Constructor - should be just one copy per program.
   */
  BatchJobManager();
  ~BatchJobManager();
  /**
   * Returns a pointer to the ith job in the list (indexed from 0).
   */
  BatchJob* GetJob(int i) {
    return JobList[i];
  }

  /**
   * Deletes a job in the list given its pointer.
   * returns true if it finds the job in the lsit and destroys it, else false.
   */
  bool DeleteJob(BatchJob* p_job);

  /**
   * Returns the working directory
   */
  const char* GetWorkDirectory() {
    return workdir.c_str();
  }

  /**
   * Returns the current JobList
   */
  std::vector<BatchJob*>* GetJobList() {
    return &JobList;
  }

  /**
   * Sets the working directory - the working directory will default to the current directory unless this function is used to change it
   */
  void SetWorkDirectory(const char* dir) {
    workdir = dir;
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

signals:
  void jobAdded(BatchJob*);
  void jobDeleted(BatchJob*);

};

#endif
