#ifndef MITOOLBAR_H
#define MITOOLBAR_H

#include <vector>

class MIMenuBar;
class QToolBar;
class QAction;
class QMainWindow;

class MIToolBar
{
public:
    MIToolBar(MIMenuBar *mb, QMainWindow *parent = 0);
    ~MIToolBar()
    {
    }

    void AddTool(unsigned int id, const char **xpm_data);
    void AddSeparator();
    void show();

    void doUpdates();

private:
    QAction *findAction(unsigned int id);

    MIMenuBar *_mb;  // not owned by this class
    QToolBar *_tb;   // not owned by this class
    std::vector<unsigned int> _ids;
};


#endif // MITOOLBAR_H
