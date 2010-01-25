#ifndef CURRENTMIGLWIDGETACTION_H
#define CURRENTMIGLWIDGETACTION_H

#include <QAction>

class CurrentMIGLWidgetAction : public QAction
{
    Q_OBJECT
public:
    CurrentMIGLWidgetAction(const QString& text, const QString& statusTip, QMenu *menu,
                            const char *slot, const char *updateSlot = 0);

public slots:
    void update();

private slots:
    void on_triggered();

private:
    const char *_slot;
    const char *_updateSlot;
};

#endif // CURRENTMIGLWIDGETACTION_H
