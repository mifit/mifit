#include "tools.h"

#include <cstdarg>
#include <chemlib/chemlib.h>
#include <chemlib/RESIDUE_.h>
#include <conflib/conflib.h>
#include <map/maplib.h>
#include <nongui/nonguilib.h>
#include <math/mathlib.h>
#include <QDebug>
#include <QFile>
#include <QProcess>
#include <QString>
#include <QFileDialog>
#include <QFileInfo>
#include <QDir>
#include <QMenu>
#include <QMessageBox>
#include <QSettings>
#include <QScriptEngine>
#include <script/MIFitScriptObject.h>
#include <util/utillib.h>
#include <vector>

#include "Application.h"
#include "CustomJobDialog.h"
#include "Displaylist.h"
#include "EMap.h"
#include "jobs/jobslib.h"
#include "id.h"
#include "macafxwin.h"
#include "MIEventHandlerMacros.h"
#include "MIGLWidget.h"
#include "MIMainWindow.h"
#include "molw.h"
#include <script/LocalSocketScript.h>

using namespace chemlib;

void Tools::OnCustom()
{
    static CustomJobDialog dlg(MIMainWindow::instance());

    static int customJobNumber = 1;

    QString jobName(QString("Custom job %1").arg(customJobNumber++));

    QSettings *settings = MIGetQSettings();
    settings->beginGroup("CustomJob");
    QString program = settings->value("executable").toString();
    QString arguments = settings->value("arguments").toString();
    QString workingDirectory = settings->value("workingDirectory", QDir::currentPath()).toString();
    bool useCurrentModel = settings->value("useCurrentModel", true).toBool();
    QString modelFile = settings->value("modelFile").toString();
    QString dataFile = settings->value("dataFile").toString();
    settings->endGroup();

    MIGLWidget *doc = MIMainWindow::instance()->currentMIGLWidget();
    if (doc != NULL)
    {
        EMap *map = doc->GetDisplaylist()->GetCurrentMap();
        if (map != NULL)
        {
            dataFile = map->pathName.c_str();
        }
        if (modelFile.isEmpty())
        {
            Molecule *model = doc->GetDisplaylist()->CurrentItem();
            if (model != NULL)
            {
                modelFile = model->pathname.c_str();
            }
        }
    }
    if (workingDirectory.isEmpty())
    {
        workingDirectory = QDir::currentPath();
    }

    dlg.setJobName(jobName);
    dlg.setProgram(program);
    dlg.setArguments(arguments);
    dlg.setModelMode(useCurrentModel ? CustomJobDialog::CURRENT : CustomJobDialog::FILE);
    dlg.setWorkingDirectory(workingDirectory);
    dlg.setModelFile(modelFile);
    dlg.setDataFile(dataFile);

    if (dlg.exec() != QDialog::Accepted)
    {
        return;
    }

    jobName = dlg.jobName();
    program = dlg.program();
    arguments = dlg.arguments();
    useCurrentModel = dlg.modelMode() == CustomJobDialog::CURRENT;
    workingDirectory = dlg.workingDirectory();
    modelFile = dlg.modelFile();
    dataFile = dlg.dataFile();

    BatchJob *job = MIMainWindow::instance()->GetJobManager()->CreateJob();

    typedef std::map<QString, QString> SubstitutionMap;
    SubstitutionMap subs;
    subs["DATA"] = dataFile;

    QDir dir(workingDirectory);
    job->setWorkingDirectory(dir.absolutePath());

    if (useCurrentModel)
    {
        MIGLWidget *doc = MIMainWindow::instance()->currentMIGLWidget();
        if (doc != NULL)
        {
            Molecule *model = doc->GetDisplaylist()->GetCurrentModel();
            if (model)
            {
                QString modelFile = dir.absoluteFilePath(QString("mifit_%1.pdb").arg(job->jobId()));
                model->SavePDBFile(modelFile.toAscii().constData());
                subs["MODEL"] = modelFile;
            }
        }
    }
    else
    {
        subs["MODEL"] = modelFile;
    }

    QString args(arguments);
    SubstitutionMap::iterator iter = subs.begin();
    args.replace("$$", "\b");
    while (iter != subs.end())
    {
        args.replace("$" + iter->first, iter->second);
        ++iter;
    }
    args.replace("\b", "$");

    job->setJobName(jobName);
    job->setProgram(program);
    job->setArguments(args);

    job->StartJob();

    settings->beginGroup("CustomJob");
    settings->setValue("executable", program);
    settings->setValue("arguments", arguments);
    settings->setValue("useCurrentModel", useCurrentModel);
    settings->setValue("workingDirectory", workingDirectory);
    settings->setValue("modelFile", modelFile);
    settings->setValue("dataFile", dataFile);
    settings->endGroup();

}

void Tools::FillToolsMenu(QMenu *parent)
{

    connect(parent, SIGNAL(aboutToShow()),
            this, SLOT(OnUpdateForJobLimit()));

    parent->addAction("Run Custom Job", this, SLOT(OnCustom()));
    parent->addSeparator();

    static LocalSocketScript scriptSocket;
    std::auto_ptr<QScriptEngine> engine(new QScriptEngine);
    MIFitScriptObject *mifitObject = new MIFitScriptObject(engine.get(), this);
    mifitObject->setJobMenu(parent);
    mifitObject->setScriptPort(scriptSocket.name());
    QScriptValue objectValue = engine->newQObject(mifitObject);
    engine->globalObject().setProperty("mifit", objectValue);

    QFile jobsJsFile(Application::instance()->GetMolimageHome().c_str() + QString("/jobs.js"));
    if (jobsJsFile.open(QIODevice::ReadOnly))
    {
        QString script = jobsJsFile.readAll();
        QScriptValue scriptResult = engine->evaluate(script, "jobs.js");
        if (engine->hasUncaughtException())
        {
            QScriptValue exception = engine->uncaughtException();
            int lineNumber = engine->uncaughtExceptionLineNumber();
            QStringList backtrace = engine->uncaughtExceptionBacktrace();
            QString result = QString("Exception %1 on line %2\n\t%3")
                             .arg(exception.toString()).arg(lineNumber)
                             .arg(backtrace.join("\n\t"));
            Logger::log(result.toStdString());
        }
    }

    actions = parent->actions();
}

void Tools::OnUpdateForJobLimit()
{
    bool enable = !MIMainWindow::instance()->isJobLimit();
    foreach (QAction* act, actions)
    {
        act->setEnabled(enable);
    }

}

Tools::Tools()
    : QObject(0)
{
}


Tools&Tools::instance()
{
    static Tools _instance;
    return _instance;
}
