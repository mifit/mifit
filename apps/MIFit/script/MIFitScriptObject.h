#ifndef script_MIFitScriptObject_h
#define script_MIFitScriptObject_h

#include <QObject>
#include <QStringList>
class QMenu;
class QScriptEngine;

class MIFitScriptObject : public QObject
{
    Q_OBJECT
    Q_PROPERTY( QString version READ version )

public:
    MIFitScriptObject(QScriptEngine* engine, QObject* parent = 0);

    void setJobMenu(QMenu *jobMenu);

public slots:
    QString version();
    QString directory();
    bool writeCurrentModel(const QString& file);
    QStringList dictionaryResidueList();
    QStringList spacegroupList();
    void addJob(const QString& menuName, const QString& jobName, const QString& executable, const QStringList& arguments);

private:
    QScriptEngine *engine;
    QMenu *jobMenu;
};

#endif
