#ifndef _WIN32
#include <unistd.h>
#include <pwd.h>
#include <sys/param.h>
#include <sys/types.h>
#include <sys/stat.h>
#endif

#include <string.h>
#include <limits.h>

#include "nonguilib.h"
#include "maplib.h"
#include "corelib.h"

#include "macafxwin.h"
#include "Application.h"
#include "MIDialog.h"

WaitCursor::WaitCursor(const char* op) {
  op_name = op;
  MIBusyManager::instance()->SetLabel(op);
  MIBusyManager::instance()->SetBusy(true);
  MIBusyManager::instance()->StartWaitCursor();
}

WaitCursor::~WaitCursor() {
  MIBusyManager::instance()->SetBusy(false);
  MIBusyManager::instance()->StopWaitCursor();
}

bool WaitCursor::CheckForAbort() {
  return MIBusyManager::instance()->CheckForAbort();
}
