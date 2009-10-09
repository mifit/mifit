#ifndef script_MIFitScriptObject_h
#define script_MIFitScriptObject_h

#include <QObject>
class QScriptEngine;

class MIFitScriptObject : public QObject
{
    Q_OBJECT
    Q_PROPERTY( QString version READ version )

public:
    MIFitScriptObject(QScriptEngine* engine, QObject* parent = 0);

public slots:
    QString version();
    bool writeCurrentModel(const QString& file);

private:
    QScriptEngine* engine;
};

#endif
