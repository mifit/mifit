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
#include "ui/MIDialog.h"
#include <script/LocalSocketScript.h>

using namespace chemlib;

// convert path to system-appropriate absolute path string:
static QString buildAbsPath(const QString& path) {
  if (path.isEmpty() || path == "none")
    return path;

  return QDir::toNativeSeparators(QFileInfo(path).absoluteFilePath());
}

static QString MIExpertPy()
{
    return QDir::toNativeSeparators(
            QString("%1/MIExpert/MIExpert.py")
            .arg(Application::instance()->GetMolimageHome().c_str()));
}

static QString MIExpertScript(const QString& name)
{
    return QDir::toNativeSeparators(
            QString("%1/MIExpert/%2")
            .arg(Application::instance()->GetMolimageHome().c_str())
            .arg(name));
}

static QString pythonExe()
{
    static QString pythonExePath;
    if (pythonExePath.isEmpty() || !QFile::exists(pythonExePath)) {
        QSettings* settings = MIGetQSettings();
        pythonExePath = settings->value("pythonExe").toString();
    }
    if (pythonExePath.isEmpty() || !QFile::exists(pythonExePath)) {
#ifdef Q_OS_WIN32
        QString separator = ";";
        QString exe = "python.exe";
        QString filters = "Programs (*.exe);;All files (*.*)";
#else
        QString separator = ":";
        QString exe = "python";
        QString filters = "All files (*)";
#endif
        QString pathEnv = getenv("PATH");
        QStringList paths = pathEnv.split(separator);
        foreach (QString p, paths) {
            QDir dir(p);
            if (dir.exists(exe)) {
                pythonExePath = dir.absoluteFilePath(exe);
                break;
            }
        }
        if (pythonExePath.isEmpty()) {
            QString fileName = QFileDialog::getOpenFileName(NULL, "Select Python Executable",
                                                            "/", filters);
            if (!fileName.isEmpty())
                pythonExePath = fileName;
        }

        if (!pythonExePath.isEmpty() && QFile::exists(pythonExePath)) {
            QSettings* settings = MIGetQSettings();
            settings->setValue("pythonExe", pythonExePath);
        }
    }

    return pythonExePath;
}

bool Tools::VerifyMIExpert() {
  if (QFile(MIExpertPy()).exists()) {
    return true;
  }
  QMessageBox::critical(0, "Error", "Cannot find MIExpert");
  return false;
}

bool Tools::VerifyCCP4() {
  static bool firsttime = true;
  static bool result = false;
  if (!firsttime && result) {
    return result;
  }
  firsttime = false;

  QByteArray pdbsetOutput;
  QProcess pdbsetProcess;
  pdbsetProcess.start("pdbset");
  pdbsetProcess.closeWriteChannel();

  pdbsetProcess.waitForFinished(2000);

  if (pdbsetProcess.exitStatus() != QProcess::NormalExit) {
    pdbsetProcess.kill();
#ifdef DEBUG
    QMessageBox::warning(MIMainWindow::instance(), "Error", "Cannot find CCP4\n(Unable to run pdbset)");
    result = true;
    return true;
#else
    QMessageBox::critical(MIMainWindow::instance(), "Error", "Cannot find CCP4\n(Unable to run pdbset)");
    result = false;
    return false;
#endif
  }

  QString outputText(pdbsetProcess.readAllStandardOutput());
  if (outputText.indexOf("PDBSET") == -1)  {
    pdbsetProcess.kill();
#ifdef DEBUG
    QMessageBox::warning(MIMainWindow::instance(), "Error", "Cannot find CCP4\n(Unable to run pdbset)");
    result = true;
    return true;
#else
    QMessageBox::critical(MIMainWindow::instance(), "Error", "Cannot find CCP4\n(Unable to identify output as from pdbset)");
    result = false;
    return false;
#endif
  }
  pdbsetProcess.kill();
  result = true;
  return true;
}

void Tools::OnBindNGrind() {
}

void Tools::CIFConvertlib(const char* format)
{
}

void Tools::OnCIF2Shellx() {
}

void Tools::OnCIF2CNS() {
}

void Tools::OnMolRep() {
}

void Tools::OnRefmacRestraints() {
  if (!VerifyMIExpert() || !VerifyCCP4()) {
    return;
  }
  QString python = pythonExe();
  if (python.isEmpty())
      return;

  QString filename = Application::getOpenFileName(0, "Choose a PDB file", "PDB files (*.pdb);;All files (*.*)");

  if (filename.isEmpty()) {
    return;
  }
  BatchJob* job = MIMainWindow::instance()->GetJobManager()->CreateJob();
  job->setJobName("Refmac Restraints");
  QFileInfo workdir(filename);
  QStringList args;
  args << MIExpertPy() << "restraints"
          << "--pdbfile" << filename
          << "--workdir" << workdir.absolutePath();
  job->setProgram(python);
  job->setArguments(args);
  job->StartJob();
}

void Tools::OnRefine() {
}

void Tools::OnJobReport() {
}

void Tools::OnCoCrystalSuperPos() {
  if (!VerifyMIExpert() || !VerifyCCP4()) {
    return;
  }
  QString python = pythonExe();
  if (python.isEmpty())
      return;

  static MICocrystalSuperpositionDialog dlg(MIMainWindow::instance(), "Cocrystal superposition");
  MIData data;
  dlg.GetInitialData(data);
  if (!dlg.GetResults(data)) {
    return;
  }

  QStringList args;
  args << MIExpertPy() << "ligandoverlap";
  args << "--workdir" << buildAbsPath(data["workdir"].str.c_str());
  args << "--pdbdir" << buildAbsPath(data["pdbdir"].str.c_str());
  args << "--targetpdb" << buildAbsPath(data["targetpdb"].str.c_str());
  args << "--targetsite" << QString::number(data["x"].f, 'f', 3)
          << QString::number(data["y"].f, 'f', 3)
          << QString::number(data["z"].f, 'f', 3);

  BatchJob* job = MIMainWindow::instance()->GetJobManager()->CreateJob();
  job->setJobName("Cocrystal Superpos");
  job->setProgram(python);
  job->setArguments(args);
  job->setWorkingDirectory(data["workdir"].str.c_str());
  job->StartJob();
}

void Tools::OnSadPhasing() {
}

void Tools::OnNCSModeling() {
}

void Tools::OnCustom()
{
    static CustomJobDialog dlg(MIMainWindow::instance());

    static int customJobNumber = 1;

    QString jobName(QString("Custom job %1").arg(customJobNumber++));

    QSettings* settings = MIGetQSettings();
    settings->beginGroup("CustomJob");
    QString program = settings->value("executable").toString();
    QString arguments = settings->value("arguments").toString();
    QString workingDirectory = settings->value("workingDirectory", QDir::currentPath()).toString();
    bool useCurrentModel = settings->value("useCurrentModel", true).toBool();
    QString modelFile = settings->value("modelFile").toString();
    QString dataFile = settings->value("dataFile").toString();
    settings->endGroup();

    MIGLWidget *doc = MIMainWindow::instance()->currentMIGLWidget();
    if (doc != NULL) {
        EMap* map = doc->GetDisplaylist()->GetCurrentMap();
        if (map != NULL) {
            dataFile = map->pathName.c_str();
        }
        if (modelFile.isEmpty()) {
            Molecule* model = doc->GetDisplaylist()->CurrentItem();
            if (model != NULL) {
                modelFile = model->pathname.c_str();
            }
        }
    }
    if (workingDirectory.isEmpty()) {
        workingDirectory = QDir::currentPath();
    }

    dlg.setJobName(jobName);
    dlg.setProgram(program);
    dlg.setArguments(arguments);
    dlg.setModelMode(useCurrentModel ? CustomJobDialog::CURRENT : CustomJobDialog::FILE);
    dlg.setWorkingDirectory(workingDirectory);
    dlg.setModelFile(modelFile);
    dlg.setDataFile(dataFile);

    if (dlg.exec() != QDialog::Accepted) {
        return;
    }

    jobName = dlg.jobName();
    program = dlg.program();
    arguments = dlg.arguments();
    useCurrentModel = dlg.modelMode() == CustomJobDialog::CURRENT;
    workingDirectory = dlg.workingDirectory();
    modelFile = dlg.modelFile();
    dataFile = dlg.dataFile();

    BatchJob* job = MIMainWindow::instance()->GetJobManager()->CreateJob();

    typedef std::map<QString, QString> SubstitutionMap;
    SubstitutionMap subs;
    subs["DATA"] = dataFile;

    QDir dir(workingDirectory);
    job->setWorkingDirectory(dir.absolutePath());

    if (useCurrentModel) {
        MIGLWidget *doc = MIMainWindow::instance()->currentMIGLWidget();
        if (doc != NULL) {
            Molecule* model = doc->GetDisplaylist()->GetCurrentModel();
            if (model) {
                QString modelFile = dir.absoluteFilePath(QString("mifit_%1.pdb").arg(job->jobId()));
                model->SavePDBFile(modelFile.toAscii().constData());
                subs["MODEL"] = modelFile;
            }
        }
    } else {
        subs["MODEL"] = modelFile;
    }

    QString args(arguments);
    SubstitutionMap::iterator iter = subs.begin();
    args.replace("$$", "\b");
    while (iter != subs.end()) {
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

void Tools::FillToolsMenu(QMenu* parent) {

    connect(parent, SIGNAL(aboutToShow()),
            this, SLOT(OnUpdateForJobLimit()));

    actions += parent->addAction("Run Custom Job", this, SLOT(OnCustom()));
    parent->addSeparator();

    MIFitScriptObject* mifitObject = new MIFitScriptObject(engine, this);
    mifitObject->setJobMenu(parent);
    QScriptValue objectValue = engine->newQObject(mifitObject);
    engine->globalObject().setProperty("mifit", objectValue);

    QFile jobsJsFile(Application::instance()->GetMolimageHome().c_str() + QString("/jobs.js"));
    if (jobsJsFile.open(QIODevice::ReadOnly)) {
        QString script = jobsJsFile.readAll();
        QScriptValue scriptResult = engine->evaluate(script, "jobs.js");
        if (engine->hasUncaughtException()) {
            QScriptValue exception = engine->uncaughtException();
            int lineNumber = engine->uncaughtExceptionLineNumber();
            QStringList backtrace = engine->uncaughtExceptionBacktrace();
            QString result = QString("Exception %1 on line %2\n\t%3")
                             .arg(exception.toString()).arg(lineNumber)
                             .arg(backtrace.join("\n\t"));
            Logger::log(result.toStdString());
        }
    }
}

void Tools::OnUpdateForJobLimit() {
    bool enable = !MIMainWindow::instance()->isJobLimit();
    foreach (QAction* act, actions) {
        act->setEnabled(enable);
    }
    bool havedoc = MIMainWindow::instance()->currentMIGLWidget() != NULL;
    foreach (QAction* act, docActions) {
        act->setEnabled(enable && havedoc);
    }

}


void Tools::OnIntegrate() {
}

Tools::Tools() : QObject(0)
{
    engine = new QScriptEngine(this);
    QObject* mifitObject = new MIFitScriptObject(engine, this);
    QScriptValue objectValue = engine->newQObject(mifitObject);
    engine->globalObject().setProperty("mifit", objectValue);
}


Tools& Tools::instance() {
    static Tools _instance;
    return _instance;
}
