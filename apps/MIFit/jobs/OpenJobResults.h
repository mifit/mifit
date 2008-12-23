#ifndef mifit_jobs_OpenJobResults_h
#define mifit_jobs_OpenJobResults_h

#include <string>

class OpenJobResults {
public:
  static void prompt(const std::string& workdir, const std::string& jobName = "");
};

#endif
