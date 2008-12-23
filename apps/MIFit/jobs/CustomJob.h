#ifndef mifit_jobs_CustomJob_h
#define mifit_jobs_CustomJob_h

#include <string>
#include "BatchJob.h"
#include "core/MIData.h"

class CustomJob : public BatchJob {

  std::string command;

  virtual void doJobFinished();

public:
  CustomJob();

  virtual std::string Info();

  bool StartJob();

};

#endif
