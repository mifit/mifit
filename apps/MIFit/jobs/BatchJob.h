#ifndef mifit_jobs_BatchJob_h
#define mifit_jobs_BatchJob_h

#include <string>
#include <boost/signal.hpp>

#include "core/corelib.h"

class MIGLWidget;

class BatchJob {
public:
  void ShowLog();
  virtual std::string Info();
  std::string FinishedFile;
  std::string LogFile;
  bool Cleaned;
  /**
   * A string for a command to be run after the job is over
   *  such as "/bin/rm file.tmp"
   */
  void AddtoCleanup(const std::string&);
  /**
   * Aborts the job by sending a kill to the operating system
   */
  bool AbortJob();
  /**
   * Add this commsnd to the job batch file
   */
  bool WriteCommand(const std::string&);
  /**
   * Start the job running
   */
  virtual bool StartJob();
  /**
   * returns true if the job is still running
   */
  bool IsRunning();
  /**
   * Handle some cleanup stuff for an individual job (i.e. remove the script
   * files), but don't actually delete the object. That for iff it's removed from
   * the treeview
   */
  void CleanUp();

  /**
   * Returns whether the job has been run and is finished.
   */
  bool isCompleted() {
    return completed;
  }

  /**
   * returns true if the job is a success
   */
  bool RanOK() {
    return success;
  }

  /**
   * A unique job id
   */
  unsigned long JobId;
  /**
   * A string containing the filename of a script to be called upon the job ending.
   *  The script is in internal script format to load the files.
   */
  std::string UpdateScript;
  /**
   * The name of the command file
   */
  std::string CommandFile;
  /**
   * The constructor
   *    working_directory - the directory to run the command in - max be blank
   */
  BatchJob(const std::string& working_directory);
  /**
   * The destructor.
   */
  virtual ~BatchJob();

  void SetDocument(MIGLWidget* doc);
  MIGLWidget *GetDocument();
  
  std::string getJobDir() const;
  void setSettings(const MIData& jobSettings);
  MIData& getSettings();

  boost::signal1<void, BatchJob*> jobChanged;

  void openResults();
  
protected:

  BatchJob();

  void setJobId();
  void setJobDir(const char* dir);
  void setSuccess(bool value);
  void openCommandFile();
  virtual void doJobFinished();

  MIData settings;
  long pid;
  bool running;
  bool completed;
  bool success;
  //FILE * UpdateScriptFp;
  FILE* CommandFileFp;
  std::string jobDir;
  std::vector<std::string> CleanupList;
  MIGLWidget* m_doc;
  /**
   * checks to see if job has ended and the return status
   */
  bool HasEnded();
};

#endif
