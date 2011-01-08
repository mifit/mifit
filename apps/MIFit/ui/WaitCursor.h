#ifndef mifit_ui_WaitCursor_h
#define mifit_ui_WaitCursor_h

#include <string>

//@{
// Creates an hourglass cursor to tell the user to be patient.
// Declare as an automatic variable at the the beginning of a long subroutine.
// When the subroutine exits, the WaitCursor object will automatically
// and the cursor restored to its previous state.
//@}
class WaitCursor
{
    std::string op_name;
public:
    bool CheckForAbort();
    WaitCursor(const char *op);
    ~WaitCursor();
};

#endif // mifit_ui_WaitCursor_h
