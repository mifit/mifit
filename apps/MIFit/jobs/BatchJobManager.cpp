#include "nonguilib.h"
#include "BatchJobManager.h"

#include "BatchJob.h"
#include "CustomJob.h"
#include "corelib.h"
#include <boost/bind.hpp>

using namespace std;

BatchJobManager::BatchJobManager() {
  timer = new QTimer(this);
  connect(timer, SIGNAL(timeout()), this, SLOT(OnTimer()));
}

BatchJobManager::~BatchJobManager() {
  delete timer;
  vector<BatchJob*>::iterator i, e;
  i = JobList.begin();
  e = JobList.end();
  for (; i != e; i++) {
    delete *i;
  }
}

BatchJob* BatchJobManager::CreateJob() {
  BatchJob* job = new BatchJob(workdir.c_str());
  JobList.push_back(job);
  jobAdded(job);
  job->jobChanged.connect(boost::bind(&BatchJobManager::jobChanged, this, _1));
  return job;
}

CustomJob* BatchJobManager::CreateCustomJob() {
  CustomJob* job = new CustomJob;
  JobList.push_back(job);
  jobAdded(job);
  job->jobChanged.connect(boost::bind(&BatchJobManager::jobChanged, this, _1));
  return job;
}

void BatchJobManager::jobChanged(BatchJob*) {
  checkAllJobs();
}

void BatchJobManager::CleanSucc() {
  int i, size;
  BatchJob* job;
  vector<BatchJob*> tokill;
  size = JobList.size();
  for (i = 0; i < size; i++) {
    job = *(JobList.begin()+i);
    if (!job->IsRunning() && job->RanOK()) {
      tokill.push_back(job);
    }
  }
  size = tokill.size();
  for (i = 0; i < size; i++) {
    DeleteJob(tokill[i]);
  }
}

void BatchJobManager::CleanAll() {
  int i, size;
  BatchJob* job;
  vector<BatchJob*> tokill;
  size = JobList.size();
  for (i = 0; i < size; i++) {
    job = *(JobList.begin()+i);
    if (!job->IsRunning()) {
      tokill.push_back(job);
    }
  }
  size = tokill.size();
  for (i = 0; i < size; i++) {
    DeleteJob(tokill[i]);
  }
}

bool BatchJobManager::DeleteJob(BatchJob* job) {
  //todo - add code to detect if job running
  //and abort the job before deleting
  // or refuse to delete job when running?
  if (!job)
    return false;
  if (job->IsRunning()) {
    Logger::message("You cannot delete a job while it is still running");
    return false;
  }
  for (size_t i = 0; i < JobList.size(); i++) {
    if (job == JobList[i]) {
      delete *(JobList.begin()+i);
      JobList.erase(JobList.begin()+i);
      jobDeleted(job);
      return true;
    }
  }
  return false;
}

void BatchJobManager::DetachJob(BatchJob* job) {
  for (size_t i = 0; i < JobList.size(); i++) {
    if (job == JobList[i]) {
      delete *(JobList.begin()+i);
      JobList.erase(JobList.begin()+i);
      jobDeleted(job);
    }
  }
}

void BatchJobManager::ShowLogFile(BatchJob* b) {
  if (b)
    b->ShowLog();
}

void BatchJobManager::stopTimer() {
  if (timer != NULL) {
    timer->stop();
  }
}

void BatchJobManager::startTimer() {
  if (timer->isActive()) {
    return;
  }
  timer->start(3000);
}

void BatchJobManager::OnTimer() {
  if (JobList.size() != 0) {
    checkAllJobs();
  }
}

void BatchJobManager::checkAllJobs() {
  std::vector<BatchJob*>::iterator jobIter = JobList.begin();
  bool running = false;
  while (jobIter != JobList.end()) {
    BatchJob* job = *jobIter;
    ++jobIter;
    if (job && job->IsRunning()) {
      running = true;
    }
  }
  if (running) {
    startTimer();
  } else {
    stopTimer();
  }
}

int BatchJobManager::numberOfRunningJobs() {
  int count = 0;
  std::vector<BatchJob*>::iterator iter = JobList.begin();
  while (iter != JobList.end()) {
    BatchJob* job = *iter;
    ++iter;
    if (job && job->IsRunning()) {
      ++count;
    }
  }
  return count;
}

