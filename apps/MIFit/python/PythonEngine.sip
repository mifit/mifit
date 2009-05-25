%Module PythonEngine

class PythonEngine
{

%TypeHeaderCode
#include <PythonEngine.h>
%End

public:
    void write(const char * text);
    void writeln(const char * text);
    void flush();

private:
    PythonEngine(const PythonEngine &);
};
