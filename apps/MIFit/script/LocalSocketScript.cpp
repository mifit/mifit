#include "LocalSocketScript.h"
#include <nongui/nonguilib.h>
#include <QBuffer>
#include <QDataStream>
#include <QLocalServer>
#include <QLocalSocket>
#include <QScriptEngine>
#include <QStringList>
#include "MIFitScriptObject.h"

#ifdef _WIN32
#include <process.h>
#endif

LocalSocketScript::LocalSocketScript(QObject* parent)
    : QObject(parent)
{
    engine = new QScriptEngine(this);
    QObject* mifitObject = new MIFitScriptObject(engine, this);
    QScriptValue objectValue = engine->newQObject(mifitObject);
    engine->globalObject().setProperty("mifit", objectValue);

    uint id = static_cast<uint>(getpid());
    name_ = "MIFit-" + QString::number(id, 16);

    localServer = new QLocalServer(this);
    localServer->listen(name_);
    connect(localServer, SIGNAL(newConnection()),
            this, SLOT(handleConnection()));

    Logger::log(("Local scripting socket: " + name_).toStdString());
}

QString LocalSocketScript::name() const
{
    return name_;
}

void LocalSocketScript::handleConnection()
{
    QLocalSocket* connection = localServer->nextPendingConnection();
    connect(connection, SIGNAL(disconnected()),
            connection, SLOT(deleteLater()));
    if (connection) {
        QString script;
        bool moreInput = true;
        QDataStream stream(connection);
        stream.setVersion(QDataStream::Qt_4_0);
        while (moreInput) {
            if (!connection->waitForReadyRead())
                return;
            if (stream.atEnd())
                return;
            QString str;
            stream >> str;
            int i = str.indexOf('\b');
            if (i >= 0) {
                str = str.mid(0, i);
                moreInput = false;
            }
            script += str;
        }
        Logger::log(("script: " + script).toStdString());

        QString result;
        QScriptValue scriptResult = engine->evaluate(script, name_);
        if (engine->hasUncaughtException()) {
            QScriptValue exception = engine->uncaughtException();
            int lineNumber = engine->uncaughtExceptionLineNumber();
            QStringList backtrace = engine->uncaughtExceptionBacktrace();
            result = QString("Exception %1 on line %2\n\t%3")
                              .arg(exception.toString()).arg(lineNumber)
                              .arg(backtrace.join("\n\t"));
        } else {
            result = scriptResult.toString();
        }

        Logger::log(("script result: " + result).toStdString());

        QByteArray data;
        QDataStream outStream(&data, QIODevice::WriteOnly);
        outStream.setVersion(QDataStream::Qt_4_0);
        outStream << static_cast<qint64>(0);
        outStream << result;
        outStream.device()->seek(0);
        outStream << static_cast<qint64>(data.size() - sizeof(qint64));
        connection->write(data);

        connection->disconnectFromServer();
        Logger::log("script done");
    }
}
