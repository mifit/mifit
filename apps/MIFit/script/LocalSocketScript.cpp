#include "LocalSocketScript.h"
#include <QBuffer>
#include <QDataStream>
#include <QLocalServer>
#include <QLocalSocket>
#include <QScriptEngine>
#include <QStringList>
#include <ui/Logger.h>
#include "MIFitScriptObject.h"

#ifdef _WIN32
#include <process.h>
#endif

LocalSocketScript::LocalSocketScript(QObject *parent)
    : QObject(parent),
      engine(0)
{
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
    if (!engine)
    {
        engine = new QScriptEngine(this);
        QObject *mifitObject = new MIFitScriptObject(engine, this);
        QScriptValue objectValue = engine->newQObject(mifitObject);
        engine->globalObject().setProperty("mifit", objectValue);
    }

    QLocalSocket *connection = localServer->nextPendingConnection();
    connect(connection, SIGNAL(disconnected()),
            connection, SLOT(deleteLater()));
    if (connection)
    {
        QString script;
        bool moreInput = true;
        QDataStream stream(connection);
        stream.setVersion(QDataStream::Qt_4_5);
        quint32 dataSize = 0;
        while (moreInput)
        {
            if (!connection->waitForReadyRead())
                break;
            if (dataSize == 0)
            {
                if (connection->bytesAvailable() < sizeof(quint32))
                    continue;
                stream >> dataSize;
            }
            if (connection->bytesAvailable() < dataSize)
                continue;
            stream >> script;
            moreInput = false;
        }
        Logger::debug(("script: " + script).toStdString());

        QString result;
        QScriptValue scriptResult = engine->evaluate(script, name_);
        if (engine->hasUncaughtException())
        {
            QScriptValue exception = engine->uncaughtException();
            int lineNumber = engine->uncaughtExceptionLineNumber();
            QStringList backtrace = engine->uncaughtExceptionBacktrace();
            result = QString("Exception %1 on line %2\n\t%3")
                     .arg(exception.toString()).arg(lineNumber)
                     .arg(backtrace.join("\n\t"));
        }
        else
        {
            result = scriptResult.toString();
        }

        Logger::debug(("script result: " + result).toStdString());

        QByteArray data;
        QDataStream out(&data, QIODevice::WriteOnly);
        out.setVersion(QDataStream::Qt_4_5);

        out << quint32(0) << result;

        out.device()->seek(0);
        out << quint32(data.size() - sizeof(quint32));
        connection->write(data);
        connection->flush();
        connection->waitForDisconnected();
        connection->close();
        Logger::debug("script done");
    }
}
