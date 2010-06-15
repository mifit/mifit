#include "ui/Application.h"

#include <iostream>
#include <QApplication>
#include <QTimer>
#include "ui/MIMainWindow.h"

int main(int argc, char **argv)
{
    int result = -1;
    Application app(argc, argv);
#ifdef Q_OS_LINUX
    QApplication::setStyle("Plastique");
#endif

    try
    {
        MIMainWindow *mw = MIMainWindow::instance();
        QTimer::singleShot(0, &app, SLOT(AfterInit()));
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
    return result;
}
