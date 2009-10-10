#include <nongui/nonguilib.h>
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
#ifdef __APPLE__
  s = "/bin/sleep 10\n";
#else
  s = "/usr/bin/sleep 10\n";
#endif
#endif
  job->StartJob();
}

