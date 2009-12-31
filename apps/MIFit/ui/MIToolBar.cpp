#include <nongui/nonguilib.h>
#include <util/utillib.h>

#include "MIToolBar.h"

#include <QToolBar>
#include <QMenuBar>

#include "MIMenuBar.h"
#include "MIMainWindow.h"

MIToolBar::MIToolBar(MIMenuBar *mb, QMainWindow *parent)
    : _mb(mb)
{
    if (!parent)
        parent = MIMainWindow::instance();
    _tb = parent->addToolBar("MIFit tools");
    _tb->setObjectName("MIFit tools");
    _tb->setIconSize(QSize(20, 20));
}

void MIToolBar::AddTool(unsigned int id, const char **xpm_data)
{
    QAction *act = findAction(id);
    if (!act)
    {
        std::string foo = ::format("Action %d not found in AddTool", id);
        Logger::debug(foo.c_str());
        return;
    }

    act->setIcon(QIcon(QPixmap(xpm_data)));
    _tb->addAction(act);
    _ids.push_back(id);
}


void MIToolBar::AddSeparator()
{
    _tb->addSeparator();
}

void MIToolBar::show()
{
    doUpdates();
    _tb->show();
}


QAction*MIToolBar::findAction(unsigned int id)
{
    return _mb->findAction(id);
}

void MIToolBar::doUpdates()
{
    for (size_t i = 0; i < _ids.size(); ++i)
    {
        _mb->doUpdate(_ids[i]);
    }
}

QToolBar *MIToolBar::toolBar()
{
    return _tb;
}
