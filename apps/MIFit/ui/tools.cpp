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
  if (!VerifyMIExpert() || !VerifyCCP4()) {
    return;
  }
  QString python = pythonExe();
  if (python.isEmpty())
      return;

  BatchJob* job;

  static MIBindNGrindDialog dlg(MIMainWindow::instance(), "Cocrystal Solution");
  MIData data;
  dlg.GetInitialData(data);
  if (!dlg.GetResults(data)) {
    return;
  }

  QStringList args;
  args << MIExpertPy() << "bng";

  for (unsigned int i=0; i < data["hklin"].strList.size(); ++i ) {
    args << "--hklin" << buildAbsPath(data["hklin"].strList[i].c_str());
    args << "--workdir" << "none";
  }
  args << "--pdbin" << buildAbsPath(data["pdbin"].str.c_str());
  args << "--process_engine" << data["process_engine"].str.c_str();

  args << "--arpwarpmap" << (data["arpwarpmap"].b ? "yes" : "no");

  if (data["detector_constants"].str.size()) {
    args << "--detector_constants" << data["detector_constants"].str.c_str();
  }

  if (data["spacegroup_no"].radio != 0) {
      args << "--spacegroup_no" << QString::number(data["spacegroup_no"].radio);
  }
  args << "--reference_mtz" << buildAbsPath(data["reference_mtz"].str.c_str());

  args << "--multi_search" << (data["multi_search"].b ? "yes" : "no");

  args << "--libfile" << buildAbsPath(data["libfile"].str.c_str());

  args << data["viewpoint_method"].str.c_str();
  if (data["bngsummary"].str.size()) {
      args << "--bngsummary" << buildAbsPath(data["bngsummary"].str.c_str());
  }
  args << "--mifit" << "no";

  job = MIMainWindow::instance()->GetJobManager()->CreateJob();
  job->setJobName("Cocrystal Solution");

  QString lig;
  /* Handle Place Ligand checkbox */
  if (data["place_ligand"].b && data["viewpoint_method"].str.size()) {
    if (data["ligand_from_dictionary"].b) {
      RESIDUE* dictres = MIFitDictionary()->GetDictResidue(data["ligand_name"].str.c_str(), 0);
      RESIDUE* reslist = new RESIDUE(*dictres);
      Molecule model(reslist, "Dictionary", NULL, NULL, 0, MoleculeType::Other);
      lig = "dictlig.pdb";
      FILE* ligOut = fopen(lig.toAscii().constData(), "w");
      chemlib::PDB fpdb;
      MIMolInfo mi;
      mi.res = model.getResidues();
      mi.bonds = model.getBonds();
      fpdb.Write(ligOut, mi);
      fclose(ligOut);
    } else {
      lig = data["ligand_name"].str.c_str();
    }
    args << "--fragfit" << lig;
  }
  args << "--molimagehome" << buildAbsPath(Application::instance()->MolimageHome.c_str());

  QFileInfo hklin(data["hklin"].strList[0].c_str());
  QDir workdir = hklin.absoluteDir();
  job->setProgram(python);
  job->setArguments(args);
  job->StartJob();
}

void Tools::CIFConvertlib(const char* format)
{
  QString python = pythonExe();
  if (python.isEmpty())
      return;

  BatchJob* job = MIMainWindow::instance()->GetJobManager()->CreateJob();

  job->setProgram(python);
  QStringList args;
  args << MIExpertScript("mi_convertlib_ui.py")
          << "--refprogram" << format;
  job->setArguments(args);
  job->setWorkingDirectory(Application::instance()->latestFileBrowseDirectory(""));

  job->StartJob();
}

void Tools::OnCIF2Shellx() {
  CIFConvertlib("shelx");
}

void Tools::OnCIF2CNS() {
  CIFConvertlib("cns");
}

void Tools::OnMolRep() {
  if (!VerifyMIExpert() || !VerifyCCP4()) {
    return;
  }
  QString python = pythonExe();
  if (python.isEmpty())
      return;

  static MIMolRepDialog dlg(MIMainWindow::instance(), "Molecular Replacement");
  MIData data;
  dlg.GetInitialData(data);
  if (!dlg.GetResults(data)) {
    return;
  }

  QStringList args;
  args << MIExpertPy() << "molrep";
  args << "--engine" << data["engine"].str.c_str();
  if (data["spacegroup_no"].radio != 0) {
      args << "--spacegroup" << QString::number(data["spacegroup_no"].radio);
  } else if (data["sg_search"].b) {
    args << "--sg_search" << "yes";
  }
  args << "--pdbfile" << buildAbsPath(data["model"].str.c_str())
          << "--mtzfile" << buildAbsPath(data["mtzfile"].str.c_str())
          << "--workdir" << buildAbsPath(data["workdir"].str.c_str())
          << "--multi_search" << (data["multi_search"].b ? "yes": "no")
          << "--match_pdbin" << (data["match_pdbin"].b ? "yes": "no")
          << "--copies" << QString::number(data["copies"].u);
  if (!data["fixed_pdb"].str.empty())
    args << "--fixed_pdb" << buildAbsPath(data["fixed_pdb"].str.c_str());

  BatchJob* job = MIMainWindow::instance()->GetJobManager()->CreateJob();
  job->setJobName("Molrep");
  job->setProgram(python);
  job->setArguments(args);
  job->setWorkingDirectory(data["workdir"].str.c_str());
  job->StartJob();
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
  if (!VerifyMIExpert() || !VerifyCCP4()) {
    return;
  }
  QString python = pythonExe();
  if (python.isEmpty())
      return;

  static MIRefinementDialog dlg(MIMainWindow::instance(), "Refinement");
  MIData data;
  dlg.GetInitialData(data);
  if (!dlg.GetResults(data)) {
    return;
  }

  QStringList args;
  args << MIExpertPy() << "refine"
          << "--mifithome" << buildAbsPath(Application::instance()->MolimageHome.c_str())
          << "--workdir" << buildAbsPath(data["workdir"].str.c_str())
          << "--pdbfile" << buildAbsPath(data["pdbfile"].str.c_str())
          << "--mtzfile" << buildAbsPath(data["mtzfile"].str.c_str())
          << "--weight" << QString::number(data["weight"].f)
          << "--cycles" << QString::number(data["cycles"].u);
  if (!Application::instance()->ShelxHome.empty()) {
    args << "--shelx_dir" << buildAbsPath(Application::instance()->ShelxHome.c_str());
  }

  if (data["water_cycles"].u != UINT_MAX) {
    args << "--water_cycles" << QString::number(data["water_cycles"].u);
  }
  if (data["build_cycles"].u != UINT_MAX) {
    args << "--build_cycles" << QString::number(data["build_cycles"].u);
  }

  args << "--bref_type" << data["bref_type"].str.c_str();
  args << "--engine" << data["engine"].str.c_str();

  if (data["use_max_res"].b) {
    args << " --max_res " << QString::number(data["max_res"].f);
  }
  if (!data["libfile"].str.empty()) {
    args << "--libfile" << buildAbsPath(data["libfile"].str.c_str());
  }
  if (!data["tls_file"].str.empty()) {
    args << "--tls_file" << buildAbsPath(data["tls_file"].str.c_str());
  }

  BatchJob* job = MIMainWindow::instance()->GetJobManager()->CreateJob();
  job->setJobName("Refinement");
  job->setProgram(python);
  job->setArguments(args);
  job->setWorkingDirectory(data["workdir"].str.c_str());
  job->StartJob();
}

void Tools::OnJobReport() {
  if (!VerifyMIExpert() || !VerifyCCP4()) {
    return;
  }
  QString python = pythonExe();
  if (python.isEmpty())
      return;

  static MIJobReportDialog dlg(MIMainWindow::instance(), "Job Report");
  MIData data;
  dlg.GetInitialData(data);
  if (!dlg.GetResults(data)) {
    return;
  }

  QStringList args;
  args << MIExpertPy() << "deposit3d"
          << "--workdir" << buildAbsPath(data["workdir"].str.c_str())
          << "--mtzfile" << buildAbsPath(data["mtzfile"].str.c_str())
          << "--pdbfile" << buildAbsPath(data["pdbfile"].str.c_str())
          << "--molimagehome" << buildAbsPath(Application::instance()->MolimageHome.c_str())
          << "--libfile" << buildAbsPath(data["libfile"].str.c_str())
          << "--seqfile" << buildAbsPath(data["seqfile"].str.c_str())
          << "--templatefile" << buildAbsPath(data["templatefile"].str.c_str())
          << "--datalogfile" << buildAbsPath(data["datalogfile"].str.c_str())
          << "--cif_write" << (data["cif_write"].b ? "yes" : "no");

  if (data["map_write"].b) {
    args << "--map_write" << "yes"
            << "--map_border" << QString::number(data["map_border"].f);
  }
  args << "--text_write" << (data["text_report"].b ? "yes" : "no")
          << "--html_write" << (data["html_report"].b ? "yes" : "no")
          << "--hkl_write" << (data["hkl_write"].b ? "yes" : "no");
  if (data["html_report"].b) {
    if (data["report_title"].str.empty()) {
      args << "--title" << data["report_title"].str.c_str();
    }
    if (!data["image1"].str.empty()) {
      args << "--image1" << buildAbsPath(data["image1"].str.c_str());
    }
    if (!data["image2"].str.empty()) {
      args << "--image2" << buildAbsPath(data["image2"].str.c_str());
    }
    if (!data["image3"].str.empty()) {
      args << "--image3" << buildAbsPath(data["image3"].str.c_str());
    }
  }
  if (!data["rootname"].str.empty()) {
    args << "--rootname" << buildAbsPath(data["rootname"].str.c_str());
  }

  BatchJob* job = MIMainWindow::instance()->GetJobManager()->CreateJob();
  job->setJobName("Report");
  job->setProgram(python);
  job->setArguments(args);
  job->setWorkingDirectory(data["workdir"].str.c_str());
  job->StartJob();
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
  if (!VerifyMIExpert() || !VerifyCCP4()) {
    return;
  }
  QString python = pythonExe();
  if (python.isEmpty())
      return;

  static MISadPhasingDialog dlg(MIMainWindow::instance(), "SAD Phasing");
  MIData data;
  dlg.GetInitialData(data);
  if (!dlg.GetResults(data)) {
    return;
  }

  QStringList args;
  args << MIExpertPy() << "sadphase";
  args << "--molimagehome \"" + buildAbsPath(Application::instance()->MolimageHome.c_str());
  args << "--workdir" << buildAbsPath(data["workdir"].str.c_str());
  args << "--saddatafile" << buildAbsPath(data["saddatafile"].str.c_str());
  args << "--sitefile" << buildAbsPath(data["sitefile"].str.c_str());
  args << "--scatterer" << data["scatterer"].str.c_str();
  args << "--sitenumber" << QString::number(data["sitenumber"].u);
  args << "--ssnumber" << QString::number(data["ssnumber"].u);
  args << "--separation" << QString::number(data["separation"].f);
  args << "--solventfraction" << QString::number(data["solventfraction"].f);
  args << "--siterefinemethod" << data["siterefinemethod"].str.c_str();
  if (!Application::instance()->ShelxHome.empty()) {
    args << "--shelx_dir" << buildAbsPath(Application::instance()->ShelxHome.c_str());
  }
  args << "--bothhands" << (data["bothhands"].b ? "yes" : "no");
  if (data["change_spacegroup"].b ) {
    args << "--spacegroup_no" << QString::number(data["spacegroup_no"].radio);
  }

  BatchJob* job = MIMainWindow::instance()->GetJobManager()->CreateJob();
  job->setJobName("SAD phasing");
  job->setProgram(python);
  job->setArguments(args);
  job->setWorkingDirectory(data["workdir"].str.c_str());
  job->StartJob();
}

void Tools::OnNCSModeling() {
  if (!VerifyMIExpert() || !VerifyCCP4()) {
    return;
  }
  QString python = pythonExe();
  if (python.isEmpty())
      return;

  static MINCSModelingDialog dlg(MIMainWindow::instance(), "NCS Modeling");
  MIData data;
  dlg.GetInitialData(data);
  if (!dlg.GetResults(data)) {
    return;
  }

  //Write a PDB of current model
  MIGLWidget* doc = MIMainWindow::instance()->currentMIGLWidget();
  if (!doc)
    return;
  Molecule* model = doc->GetDisplaylist()->GetCurrentModel();
  QDir workdir(data["workdir"].str.c_str());
  QString pdbout = buildAbsPath(workdir.absoluteFilePath(
          QString("%1_out.pdb").arg(QFileInfo(doc->GetTitle().c_str()).fileName())));
  model->SavePDBFile(pdbout.toAscii().constData());

  QStringList args;
  args << MIExpertPy() << "ncsmodeler"
          << "--pdbfile" << pdbout
          << "--targetchain" << data["chainid"].str.c_str();
  if (!data["model"].str.empty()) {
    args << "--preserve_model" << buildAbsPath(data["model"].str.c_str());
  }
  if (!data["maskadditions"].str.empty()) {
    args << "--maskextras" << buildAbsPath(data["maskadditions"].str.c_str());
  }
  if (!data["mtzdata"].str.empty()) {
    args << "--mtzfile" << buildAbsPath(data["mtzdata"].str.c_str());
    switch (data["phasecalc"].radio) {
      case 0:
        break;
      case 1:
        args << "--phase_prob" << "yes";
        break;
      case 2:
        args << "--phase_comb" << "yes";
        break;
    }
  }

  BatchJob* job = MIMainWindow::instance()->GetJobManager()->CreateJob();
  job->setJobName("NCS Modeling");
  job->setProgram(python);
  job->setArguments(args);
  job->setWorkingDirectory(workdir.absolutePath());
  job->StartJob();
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

    actions += parent->addAction("&Integrate with d*TREK", this, SLOT(OnIntegrate()));
    actions += parent->addAction("&SAD Phasing", this, SLOT(OnSadPhasing()));

    QMenu* convcif_menu = new QMenu("&Convert CIF to", parent);
    convcif_menu->addAction("&SHELX", this, SLOT(OnCIF2Shellx()));
    convcif_menu->addAction("&CNS / CNX", this, SLOT(OnCIF2CNS()));
    actions += parent->addMenu(convcif_menu);

    actions += parent->addAction("S&et Refmac5 restraints", this, SLOT(OnRefmacRestraints()));
    actions += parent->addAction("&Molecular Replacement", this, SLOT(OnMolRep()));
    actions += parent->addAction("&Refinement", this, SLOT(OnRefine()));
    actions += parent->addAction("C&ocrystal Solution", this, SLOT(OnBindNGrind()));
    actions += parent->addAction("Re&port", this, SLOT(OnJobReport()));

    actions += parent->addAction("Run Custom Job", this, SLOT(OnCustom()));

    parent->addSeparator();
    docActions += parent->addAction("&NCS Modeling", this, SLOT(OnNCSModeling()));
    docActions += parent->addAction("Cocr&ystal Superposition", this, SLOT(OnCoCrystalSuperPos()));
#ifdef DEBUG
    QAction* action = parent->addAction("Run &Test Job...", this, SLOT(OnRunTestJob()));
    action->setStatusTip("Run a test job that waits 5 seconds then terminates");
    docActions += action;
#endif

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
  if (!VerifyMIExpert() || !VerifyCCP4()) {
    return;
  }
  QString python = pythonExe();
  if (python.isEmpty())
      return;

  static MIIntegrateDialog dlg(MIMainWindow::instance(), "Integrate");
  MIData data;
  dlg.GetInitialData(data);
  if (!dlg.GetResults(data)) {
    return;
  }

  BatchJob* job = MIMainWindow::instance()->GetJobManager()->CreateJob();
  job->setJobName("Integrate");

  QStringList args;
  args << MIExpertPy() << "integrate";
  args << "--template_image" << data["template_image"].str.c_str();
  if (data["detector_constants"].str.size()) {
    args << "--detector_constants" << buildAbsPath(data["detector_constants"].str.c_str());
  }
  args << "--spacegroup" << QString::number(data["spacegroup_no"].radio);
  if (data["first_image"].u != UINT_MAX) {
      args << "--first_image" << QString::number(data["first_image"].u);
  }
  if (data["last_image"].u != UINT_MAX) {
      args << "--last_image" << QString::number(data["last_image"].u);
  }
  if (data["integrate_resolution"].str.size()) {
    args << "--integrate_resolution" << data["integrate_resolution"].str.c_str();
  }

  job->setProgram(python);
  job->setArguments(args);

  job->StartJob();
}

void Tools::OnRunTestJob() {
  BatchJob* job =   MIMainWindow::instance()->GetJobManager()->CreateJob();
  try {
    TestJob testJob;
    testJob.StartJob(job);
  }
  catch (...) {
    Logger::message("Job Failed!!");
    return;
  }
  WaitCursor* wait = new WaitCursor("Waiting for test job to finish");
  while (wait->CheckForAbort() == false && job->isRunning()) {
#ifdef _WIN32
    Sleep(100);
#else
    usleep(100000);
#endif
  }
  delete wait;
}


Tools::Tools() : QObject(0) {
}

Tools& Tools::instance() {
    static Tools _instance;
    return _instance;
}
