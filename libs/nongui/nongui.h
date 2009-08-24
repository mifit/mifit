#ifndef MI_NONGUI_H
#define MI_NONGUI_H

#ifdef _MSC_VER
// Disable the warning about "conversion ... possible loss of data"
// until use of size_t is implemented
#pragma warning(disable:4267)
#pragma warning(disable:4244)
#endif

//
// the classes and functions defined in this library are intended to be
// used from non-gui applications
//
// gui applications may wish to include this header file, but *not* link
// against libnongui and provide their own implementations of these
// functions
//

#include <string>

class Logger {
public:
  static void log(const char* format, ...);
  static void debug(const char* format, ...);

  static void message(const char* format, ...);
  static void footer(const char* format, ...);


  static void log(const std::string& str) {
    log(str.c_str());
  }

  static void debug(const std::string& str) {
    debug(str.c_str());
  }

  static void message(const std::string& str) {
    message(str.c_str());
  }

  static void footer(const std::string& str) {
    footer(str.c_str());
  }

};

//@{
// Creates an hourglass cursor to tell the user to be patient.
// Declare as an automatic variable at the the beginning of a long subroutine.
// When the subroutine exits, the WaitCursor object will automatically
// and the cursor restored to its previous state.
//@}
class WaitCursor {
  std::string op_name;
public:
  bool CheckForAbort();
  WaitCursor(const char* op);
  ~WaitCursor();
};



#endif
