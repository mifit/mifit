#include "MIFitScriptObject.h"
#include <core/Version.h>
#include <QDebug>
#include <QFileInfo>
#include <QScriptEngine>
#include <ui/uilib.h>

MIFitScriptObject::MIFitScriptObject(QScriptEngine* engine, QObject* parent)
    : QObject(parent), engine(engine)
{
}

QString MIFitScriptObject::version()
{
    return MIFit_version;
}

bool MIFitScriptObject::writeCurrentModel(const QString& file)
{
    if (file == "test.pdb")
        engine->currentContext()->throwError("invalid file");

    MIGLWidget* doc = MIMainWindow::instance()->currentMIGLWidget();
    if (!doc) {
        engine->currentContext()->throwError("no current document");
        return false;
    }
    Molecule* model = doc->GetDisplaylist()->GetCurrentModel();
    if (!model) {
        engine->currentContext()->throwError("no current model");
        return false;
    }
    QFileInfo fileInfo(file);
    return model->SavePDBFile(fileInfo.absoluteFilePath().toAscii().constData());
}
