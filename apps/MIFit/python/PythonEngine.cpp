#include "PythonEngine.h"

#include <Python.h>
#include <sip.h>
#include "sipAPIPythonEngine.h"
extern "C" void initPythonEngine();

#include <qfile.h>
#include <qstringlist.h>
#include <qtextstream.h>

QString PythonEngine::ps1(">>> ");
QString PythonEngine::ps2("... ");

PythonEngine::PythonEngine()
  : currentPrompt(ps1), started(false) {
}

void PythonEngine::start() {
  started = true;
  char programName[] = "MIFit";
  Py_SetProgramName(programName);
  Py_Initialize();

  // Initialize bindings to this class
  initPythonEngine();

  try {
    // Get a reference to the main module.
    PyObject* main_module = PyImport_AddModule("__main__");
    if (!main_module) {
        throw PyErr_Occurred();
    }

    // Get the main module's dictionary
    main_dict = PyModule_GetDict(main_module);
    if (!main_dict) {
        throw PyErr_Occurred();
    }

    PyObject* result = PyRun_String("import code\n"
                                    "console = code.InteractiveConsole()",
                                    Py_file_input, main_dict, main_dict);
    if (!result) {
        throw PyErr_Occurred();
    }

    PyObject* pythonEngine_module = PyImport_AddModule("PythonEngine");
    if (!pythonEngine_module) {
        throw PyErr_Occurred();
    }

    PyObject* pythonEngine_dict = PyModule_GetDict(pythonEngine_module);
    if (!pythonEngine_dict) {
        throw PyErr_Occurred();
    }

    PyObject* pythonEngine = sipConvertFromInstance(this, sipClass_PythonEngine, 0);
    if (!pythonEngine) {
        throw PyErr_Occurred();
    }

    int err = PyDict_SetItemString(pythonEngine_dict, "instance", pythonEngine);
    if (err) {
        throw PyErr_Occurred();
    }

    // Redirect output and error messages to this object for
    // signaling via message signal
    result = PyRun_String("import sys, PythonEngine\n"
                            "sys.stdout = PythonEngine.instance\n"
                            "sys.stderr = PythonEngine.instance",
                            Py_file_input, main_dict, main_dict);
    if (!result) {
        throw PyErr_Occurred();
    }

  } catch (PyObject* err) {
    PyErr_Print();
  }

  command("import PythonEngine, PyQt4");
}

PythonEngine::~PythonEngine() {
  if (started) {
    Py_Finalize();
  }
}

void PythonEngine::command(const QString& text) {
  QString s = QString("console.push('''") + text + "''')";
  PyObject * result = PyRun_String(s.toAscii(), Py_eval_input, main_dict, main_dict);
  if (result) {
    if (result == Py_True) {
      currentPrompt = ps2;
    } else {
      currentPrompt = ps1;
    }
  } else {
    PyErr_Print();
  }
  prompt(currentPrompt);
}

void PythonEngine::write(const char* text) {
  message(QString(text));
}

void PythonEngine::writeln(const char* text) {
  message(QString(text) + "\n");
}

void PythonEngine::flush(void) {
}
