#include "ui/Application.h"

#include <QApplication>
#include <QTimer>
#include "../ui/MIMainWindow.h"

int main(int argc, char **argv)
{
  QApplication app(argc, argv);
  Application::instance(); // make sure App data object is initialized

  MIMainWindow *mw=MIMainWindow::instance();
  QTimer::singleShot(500, mw, SLOT(AfterInit()));

  mw->show();
  int retval=app.exec();
  delete Application::instance();
  return retval;
}
