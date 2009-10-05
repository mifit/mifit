#include <nongui/nonguilib.h>
#include "BatchJobManager.h"

#include "BatchJob.h"
#include "core/corelib.h"

using namespace std;

BatchJobManager::BatchJobManager() {
}

BatchJobManager::~BatchJobManager() {
  vector<BatchJob*>::iterator i, e;
  i = JobList.begin();
  e = JobList.end();
  for (; i != e; i++) {
    delete *i;
  }
}

BatchJob* BatchJobManager::CreateJob() {
  BatchJob* job = new BatchJob;
  JobList.push_back(job);
  jobAdded(job);
  return job;
}

void BatchJobManager::CleanSucc() {
  int i, size;
  BatchJob* job;
  vector<BatchJob*> tokill;
  size = JobList.size();
  for (i = 0; i < size; i++) {
    job = *(JobList.begin()+i);
    if (!job->isRunning() && job->isSuccess()) {
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
    if (!job->isRunning()) {
      tokill.push_back(job);
    }
  }
  size = tokill.size();
  for (i = 0; i < size; i++) {
    DeleteJob(tokill[i]);
  }
}

bool BatchJobManager::DeleteJob(BatchJob* job) {
  if (!job)
    return false;
  if (job->isRunning()) {
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

void BatchJobManager::ShowLogFile(BatchJob* b) {
  if (b)
    b->ShowLog();
}

int BatchJobManager::numberOfRunningJobs() {
  int count = 0;
  std::vector<BatchJob*>::iterator iter = JobList.begin();
  while (iter != JobList.end()) {
    BatchJob* job = *iter;
    ++iter;
    if (job && job->isRunning()) {
      ++count;
    }
  }
  return count;
}

