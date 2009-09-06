#ifndef PythonEngine_h
#define PythonEngine_h

#include <Python.h>
#include <qobject.h>

class PythonEngine : public QObject {

  Q_OBJECT

  static QString ps1;
  static QString ps2;

  PyObject* main_dict;
  QString currentPrompt;
  bool started;

  PythonEngine(const PythonEngine&);
  PythonEngine& operator=(const PythonEngine&);

public:

  PythonEngine();
  ~PythonEngine();

  void start();

  void write(const char* text);
  void writeln(const char* text);
  void flush(void);

public Q_SLOTS:
  void command(const QString& text);

Q_SIGNALS:
  void message(const QString& text);
  void prompt(const QString& prompt);
};

#endif
