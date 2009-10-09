#ifndef script_LocalSocketScript_h
#define script_LocalSocketScript_h

#include <QObject>
#include <QString>
class QLocalServer;
class QScriptEngine;

class LocalSocketScript : public QObject
{
    Q_OBJECT

public:
    LocalSocketScript(QObject* parent = 0);

    QString name() const;

private:
    QString name_;
    QLocalServer* localServer;
    QScriptEngine* engine;

private slots:
    void handleConnection();
};

#endif
