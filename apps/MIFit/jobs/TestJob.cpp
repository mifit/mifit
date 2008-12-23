#include "nonguilib.h"
#include "TestJob.h"

TestJob::TestJob() {
}

TestJob::~TestJob() {

}

void TestJob::StartJob(BatchJob* job) {
  std::string s;

  //s.Printf("New Test Job #%ld\n", job->JobId);
  Logger::log(s);
  printf("%s", s.c_str());

  // Begin writing the script
#ifdef WIN32
  s = "PAUSE\n";
#else
  s = "echo blarg > testoutput.txt\n";
  job->WriteCommand(s);
#ifdef __APPLE__
  s = "/bin/sleep 10\n";
#else
  s = "/usr/bin/sleep 10\n";
#endif
#endif
  if (job->WriteCommand(s.c_str()) == false) {
    Logger::log("Batch job failed to initialize - no job started");
    return;
  }
  job->StartJob();
}

