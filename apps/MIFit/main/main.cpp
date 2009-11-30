#include "ui/Application.h"

#include <iostream>
#include <QApplication>
#include <QTimer>
#include "ui/MIMainWindow.h"

int main(int argc, char **argv)
{
    int result = -1;
    QApplication app(argc, argv);
#ifdef Q_OS_LINUX
    QApplication::setStyle("Plastique");
#endif
    Application::instance();

    try
    {
        MIMainWindow *mw = MIMainWindow::instance();
        QTimer::singleShot(500, mw, SLOT(AfterInit()));

        mw->show();

        result = app.exec();
        delete mw;
    }
    catch (...)
    {
        std::cerr << "MIFit: unknown exception\n";
        result = 1;
    }
    delete Application::instance();
    return result;
}
