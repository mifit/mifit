#ifndef mifit_jobs_TestJob_h
#define mifit_jobs_TestJob_h

#include "BatchJob.h"

class TestJob
{
public:
    TestJob();
    ~TestJob();
    void StartJob(BatchJob *job);
};

#endif
