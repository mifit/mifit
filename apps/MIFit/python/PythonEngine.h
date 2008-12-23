#ifndef PythonEngine_h
#define PythonEngine_h

#include <boost/python.hpp>
#include <QObject>

class PythonEngine : public QObject {

  Q_OBJECT

  static QString ps1;
  static QString ps2;

  boost::python::object main_dict;
  QString currentPrompt;
  bool started;

  PythonEngine(const PythonEngine&);
  PythonEngine& operator=(const PythonEngine&);

public:

  PythonEngine();
  ~PythonEngine();

  void start();

  void write(const char* text);

private Q_SLOTS:
  void command(const QString& text);

Q_SIGNALS:
  void message(const QString& text);
  void prompt(const QString& prompt);
};

#endif
