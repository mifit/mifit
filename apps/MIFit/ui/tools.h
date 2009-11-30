#ifndef MI_TOOLS_H
#define MI_TOOLS_H

#include <QObject>
#include <QList>

class QAction;
class QMenu;

class Tools : public QObject
{
    Q_OBJECT

public:
    static Tools&instance();

    void FillToolsMenu(QMenu*);

private:
    Tools();

    QList<QAction*> actions;

private slots:
    void OnCustom();
    void OnUpdateForJobLimit();

};

#endif // ifndef MI_TOOLS_H
