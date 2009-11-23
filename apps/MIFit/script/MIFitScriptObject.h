#ifndef script_MIFitScriptObject_h
#define script_MIFitScriptObject_h

#include <QObject>
#include <QStringList>
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
    QStringList dictionaryResidueList();
    QStringList spacegroupList();

private:
    QScriptEngine* engine;
};

#endif
