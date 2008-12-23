#include "PythonEngine.h"

#include <boost/python.hpp>

#include <QFile>
#include <QTextStream>

using namespace std;
using namespace boost::python;

extern "C" void initmi();

BOOST_PYTHON_MODULE(PythonEngine) {
  class_<PythonEngine, boost::noncopyable>("PythonEngine", no_init)
    .def("write", &PythonEngine::write);
}

QString PythonEngine::ps1(">>> ");
QString PythonEngine::ps2("... ");

PythonEngine::PythonEngine()
  : main_dict(), currentPrompt(ps1), started(false) {
}

void PythonEngine::start() {
  started = true;
  Py_Initialize();

  // Initialize bindings to this class
  initPythonEngine();

  initmi();

  try {
    // Redirect output and error messages to this object for
    // signaling via message signal
    object sys = import("sys");
    sys.attr("stdout") = boost::python::ptr(this);
    sys.attr("stderr") = boost::python::ptr(this);

    main_dict = import("__main__").attr("__dict__");

    exec("import code\n"
        "console = code.InteractiveConsole()",
        main_dict, main_dict);

  } catch(error_already_set const &) {
    PyErr_Print();
  }

  command("import mi");
}

PythonEngine::~PythonEngine() {
  if (started) {
    Py_Finalize();
  }
}

void PythonEngine::command(const QString& text) {
  try {
    string s = "console.push('''" + text.toStdString() + "''')";
    object result = eval(s.c_str(), main_dict, main_dict);
    if (extract<bool>(result)) {
      currentPrompt = ps2;
    } else {
      currentPrompt = ps1;
    }
  } catch(error_already_set const &) {
    PyErr_Print();
  }
  prompt(currentPrompt);
}

void PythonEngine::write(const char* text) {
  message(QString(text));
}
