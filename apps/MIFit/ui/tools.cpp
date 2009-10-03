#include "tools.h"

#include <cstdarg>
#include <chemlib/chemlib.h>
#include <chemlib/RESIDUE_.h>
#include <conflib/conflib.h>
#include <map/maplib.h>
#include <nongui/nonguilib.h>
#include <math/mathlib.h>
#include <QFile>
#include <QProcess>
#include <QString>
#include <QFileInfo>
#include <QDir>
#include <QMenu>
#include <QMessageBox>
#include <util/utillib.h>
#include <vector>

#include "Application.h"
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

  //NOTE: mtzdump used to be the test program, but with the port to Qt
  //      it's been replaced by pdbset.  mtzdump does not terminate without
  //      additional input, but pdbset does, which makes it easier to test.
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

  BatchJob* job;
  QString log;

  static MIBindNGrindDialog dlg(MIMainWindow::instance(), "Cocrystal Solution");
  MIData data;
  dlg.GetInitialData(data);
  if (!dlg.GetResults(data)) {
    return;
  }

  QStringList args;
  args << MIExpertPy() << "bng";
  // Write config file
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
  MIData settings;
  settings["jobName"].str = "Cocrystal Solution";
  job->setSettings(settings);

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

  QString s = QString("New Cocrystal Solution job #%1").arg(job->jobId());
  Logger::log(s.toStdString());
  QFileInfo workdir(data["hklin"].strList[0].c_str());
  log = buildAbsPath(QString("%1/mifit%2.log").arg(workdir.absoluteFilePath()).arg(job->jobId()));
  job->setProgram("python");
  job->setArguments(args);
  job->setLogFile(log);
  if (data["bngsummary"].str.size()) {
    QDir dir(data["bngsummary"].str.c_str());
    if (!dir.exists()) {
        QMessageBox::critical(MIMainWindow::instance(), "Error", "HTML Summary Directory doesn't exist!\nAborting");
      return;
    }
    QStringList filters;
    filters << "bng_jobsummary*.htm";
    QFileInfoList fi=dir.entryInfoList(filters);
    int count=fi.count() +1;
    std::string htmlName=::format("bng_jobsummary_%d.htm", count);

    // TODO open browser with results
  }
  job->StartJob();
}

void Tools::CIFConvertlib(const char* format) {
  std::string s, log;
  if (!VerifyMIExpert() || !VerifyCCP4()) {
    return;
  }
  std::string filename = MIFileSelector("Choose a CIF file", MIMainWindow::instance()->GetJobManager()->GetWorkDirectory(), "", "", "CIF files (*.cif)|*.cif|All files (*.*)|*.*", 0, 0);
  if (!filename.size()) {
    return;
  }

  BatchJob* job = MIMainWindow::instance()->GetJobManager()->CreateJob();
  s=::format("New CIF restraints job #%ld", job->jobId());
  Logger::log(s);

  QFileInfo workdir(filename.c_str());
  log=::format("%s%cmifit%ld.log", workdir.absolutePath().toStdString().c_str(),QDir::separator().toAscii(), job->jobId());

  s=::format("python \"%s\" convertlib --cif \"%s\" --workdir \"%s\" --refprogram \"%s\"",
             MIExpertPy().toAscii().constData(),
             filename.c_str(),
             workdir.absolutePath().toStdString().c_str(),
             format);
  job->setCommandLine(s.c_str());
  job->setLogFile(log.c_str());
  job->StartJob();
}

void Tools::OnCIF2Shellx() {
  CIFConvertlib("shelx");
}

void Tools::OnCIF2CNS() {
  CIFConvertlib("cns");
}

void Tools::OnMolRep() {
  BatchJob* job = NULL;
  std::string jobtxt;
  std::string args;
  std::string log;

  // added to save working directory between commands and to use current working directory - dem
  if (!VerifyMIExpert() || !VerifyCCP4()) {
    return;
  }

  static MIMolRepDialog dlg(MIMainWindow::instance(), "Molecular Replacement");
  MIData data;
  dlg.GetInitialData(data);
  if (!dlg.GetResults(data)) {
    return;
  }

  args += " --engine " + data["engine"].str;
  if (data["spacegroup_no"].radio != 0) {
    args += " --spacegroup " + ::format("%d",data["spacegroup_no"].radio);
  } else if (data["sg_search"].b) {
    args += " --sg_search yes ";
  }
  args +=" --pdbfile \"" + buildAbsPath(data["model"].str.c_str()).toStdString()  + "\"" +
    " --mtzfile \"" + buildAbsPath(data["mtzfile"].str.c_str()).toStdString()  + "\"" +
    " --workdir \"" + buildAbsPath(data["workdir"].str.c_str()).toStdString() + "\"" +
    " --multi_search " + (data["multi_search"].b ? "yes": "no") +
    " --match_pdbin " + (data["match_pdbin"].b ? "yes": "no") +
    " --copies " + ::format("%d",data["copies"].u);
  if (data["fixed_pdb"].str.size())
    args += " --fixed_pdb \"" + buildAbsPath(data["fixed_pdb"].str.c_str()).toStdString()  + "\" ";

  job = MIMainWindow::instance()->GetJobManager()->CreateJob();
  MIData settings;
  settings["workingDirectory"].str = data["workdir"].str;
  settings["jobName"].str = "Molrep";
  job->setSettings(settings);
  log = buildAbsPath(QString("%1/mifit%2.log").arg(data["workdir"].str.c_str()).arg(job->jobId())).toStdString();
  jobtxt=::format("python \"%s\" molrep %s", MIExpertPy().toAscii().constData(), args.c_str());
  job->setLogFile(log.c_str());
  job->setCommandLine(jobtxt.c_str());
  job->StartJob();
}

void Tools::OnRefmacRestraints() {
  std::string s, log;
  if (!VerifyMIExpert() || !VerifyCCP4()) {
    return;
  }
  std::string filename = MIFileSelector("Choose a PDB file", MIMainWindow::instance()->GetJobManager()->GetWorkDirectory(), "", "", "PDB files (*.pdb)|*.pdb|All files (*.*)|*.*", 0, 0);
  if (!filename.size()) {
    return;
  }
  BatchJob* job = MIMainWindow::instance()->GetJobManager()->CreateJob();
  MIData settings;
  settings["jobName"].str = "Refmac Restraints";
  job->setSettings(settings);
  s=::format("New CIF to ShellX job #%ld", job->jobId());
  Logger::log(s);
  QFileInfo workdir(filename.c_str());
  log=::format("%s%cmifit%ld.log", workdir.absolutePath().toStdString().c_str(),QDir::separator().toAscii(), job->jobId());

  s=::format("python \"%s\" restraints --pdbfile \"%s\" --workdir \"%s\"",
             MIExpertPy().toAscii().constData(), filename.c_str(), workdir.absolutePath().toStdString().c_str());
  job->setCommandLine(s.c_str());
  job->setLogFile(log.c_str());
  job->StartJob();
}

void Tools::OnRefine() {
  static std::string workdir;
  BatchJob* job = NULL;
  std::string jobtxt, log, args;
  if (!VerifyMIExpert() || !VerifyCCP4()) {
    return;
  }

  static MIRefinementDialog dlg(MIMainWindow::instance(), "Refinement");
  MIData data;
  dlg.GetInitialData(data);
  if (!dlg.GetResults(data)) {
    return;
  }

  //Get the required options
  args +=" --mifithome \"" + buildAbsPath(Application::instance()->MolimageHome.c_str()).toStdString() + "\"" +
    " --workdir \"" + buildAbsPath(data["workdir"].str.c_str()).toStdString() + "\"" +
    " --pdbfile \"" + buildAbsPath(data["pdbfile"].str.c_str()).toStdString() + "\"" +
    " --mtzfile \"" + buildAbsPath(data["mtzfile"].str.c_str()).toStdString() + "\"" +
    " --weight " + ::format("%f",data["weight"].f) +
    " --cycles " + ::format("%d",data["cycles"].u);
  if (!Application::instance()->ShelxHome.empty()) {
    args += " --shelx_dir \"" + buildAbsPath(Application::instance()->ShelxHome.c_str()).toStdString() + "\"";
  }

  if (data["water_cycles"].u != UINT_MAX) {
    args += " --water_cycles " + ::format("%d",data["water_cycles"].u);
  }
  if (data["build_cycles"].u != UINT_MAX) {
    args += " --build_cycles " + ::format("%d",data["build_cycles"].u);
  }

  args += " --bref_type " + data["bref_type"].str;
  args += " --engine " + data["engine"].str;


  //Grab the optional ones
  if (data["use_max_res"].b) {
    args += " --max_res " + ::format("%f",data["max_res"].f);
  }
  if (data["libfile"].str.size()) {
    args += " --libfile \"" + buildAbsPath(data["libfile"].str.c_str()).toStdString() + "\"";
  }
  if (data["tls_file"].str.size()) {
    args += " --tls_file \"" + buildAbsPath(data["tls_file"].str.c_str()).toStdString() + "\"";
  }

  job = MIMainWindow::instance()->GetJobManager()->CreateJob();
  MIData settings;
  settings["workingDirectory"].str = data["workdir"].str;
  settings["jobName"].str = "Refinement";
  job->setSettings(settings);
  log=buildAbsPath(::format("%smifit%ld.log", (data["workdir"].str+QDir::separator().toAscii()).c_str(), job->jobId()).c_str()).toStdString();
  jobtxt=::format("python \"%s\" refine %s", MIExpertPy().toAscii().constData(), args.c_str());
  job->setLogFile(log.c_str());
  job->setCommandLine(jobtxt.c_str());
  //Find highest number error
  //refine_#_errors.txt
  std::string currentName, errorName, tempName;
  int ctrCurrent = 0, ctrOld = 0;
  QDir dir(data["workdir"].str.c_str());
  if (!dir.exists()) {
    //TODO: This should actually pop a dialog to the user too...
    return;
  }
  QStringList filters;
  filters << "refine_*.pdb";
  QFileInfoList fi=dir.entryInfoList(filters);
  if (fi.count()) {
    for (int i=0; i < fi.count(); i++) {
      std::string base=fi[i].baseName().toStdString();
      if (sscanf(base.c_str(),"refine_%d",&ctrCurrent)==1 &&
          ctrCurrent > ctrOld) {
        ctrOld=ctrCurrent;
      }
    }
    ctrOld++;
    errorName=::format("refine_%d_errors.txt", (int)ctrOld);
  } else {
    errorName = "refine_1_errors.txt";
  }

  if (errorName.size()) {
      // TODO open error file in browser
  }
  job->StartJob();
}

void Tools::OnJobReport() {
  BatchJob* job = NULL;
  std::string jobtxt, log, temp, args;
  if (!VerifyMIExpert() || !VerifyCCP4()) {
    return;
  }

  static MIJobReportDialog dlg(MIMainWindow::instance(), "Job Report");
  MIData data;
  dlg.GetInitialData(data);
  if (!dlg.GetResults(data)) {
    return;
  }

  // Required things
  args += " --workdir \"" + buildAbsPath(data["workdir"].str.c_str()).toStdString() + "\"";
  args += " --mtzfile \"" + buildAbsPath(data["mtzfile"].str.c_str()).toStdString() + "\"";
  args += " --pdbfile \"" + buildAbsPath(data["pdbfile"].str.c_str()).toStdString() + "\"";
  args += " --molimagehome \"" + buildAbsPath(Application::instance()->MolimageHome.c_str()).toStdString() + "\"";
  args += " --libfile \"" + buildAbsPath(data["libfile"].str.c_str()).toStdString() + "\"";
  args += " --seqfile \"" + buildAbsPath(data["seqfile"].str.c_str()).toStdString() + "\"";
  args += " --templatefile \"" + buildAbsPath(data["templatefile"].str.c_str()).toStdString() + "\"";
  args += " --datalogfile \"" + buildAbsPath(data["datalogfile"].str.c_str()).toStdString() + "\"";

  if (data["cif_write"].b) {
    args += " --cif_write yes";
  } else {
    args += " --cif_write no";
  }
  if (data["map_write"].b) {
    args += " --map_write yes";
    args += " --map_border=" + ::format("%f",data["map_border"].f);
  }
  if (data["text_report"].b) {
    args += " --text_write yes";
  } else {
    args += " --text_write no";
  }
  if (data["html_report"].b) {
    args += " --html_write yes";
  } else {
    args += " --html_write no";
  }
  if (data["hkl_write"].b) {
    args += " --hkl_write yes";
  } else {
    args += " --hkl_write no";
  }
  if (data["html_report"].b) {
    if (data["report_title"].str.size()) {
      args += " --title \"" + data["report_title"].str + "\"";
    }
    if (data["image1"].str.size()) {
      args += " --image1 \"" + buildAbsPath(data["image1"].str.c_str()).toStdString() + "\"";
    }
    if (data["image2"].str.size()) {
      args += " --image2 \"" + buildAbsPath(data["image2"].str.c_str()).toStdString() + "\"";
    }
    if (data["image3"].str.size()) {
      args += " --image3 \"" + buildAbsPath(data["image3"].str.c_str()).toStdString() + "\"";
    }
  }
  if (data["rootname"].str.size()) {
    args += " --rootname \"" + buildAbsPath(data["rootname"].str.c_str()).toStdString() + "\"";
  }

  // Create the job object
  job = MIMainWindow::instance()->GetJobManager()->CreateJob();
  MIData settings;
  settings["workingDirectory"].str = data["workdir"].str.c_str();
  settings["jobName"].str = "Report";
  job->setSettings(settings);
  log=buildAbsPath(::format("%smifit%ld.log", (data["workdir"].str+QDir::separator().toAscii()).c_str(), job->jobId()).c_str()).toStdString();
  jobtxt=::format("python \"%s\" deposit3d %s", MIExpertPy().toAscii().constData(), args.c_str());

  job->setLogFile(log.c_str());
  job->setCommandLine(jobtxt.c_str());
  if (data["html_report"].b) {
    std::string htmlName;
    if (data["rootname"].str.size()) {
      htmlName=::format("%s.htm", data["rootname"].str.c_str());
    } else {
      htmlName=::format("pdbdeposit.htm");
    }
    // TODO open report in browser
  }
  job->StartJob();
}

void Tools::OnCoCrystalSuperPos() {
  BatchJob* job = NULL;
  std::string jobtxt, log, temp, args;
  if (!VerifyMIExpert() || !VerifyCCP4()) {
    return;
  }

  static MICocrystalSuperpositionDialog dlg(MIMainWindow::instance(), "Cocrystal superposition");
  MIData data;
  dlg.GetInitialData(data);
  if (!dlg.GetResults(data)) {
    return;
  }

  // Required things
  args += " --workdir \"" + buildAbsPath(data["workdir"].str.c_str()).toStdString() + "\"";
  args += " --pdbdir \"" + buildAbsPath(data["pdbdir"].str.c_str()).toStdString() + "\"";
  args += " --targetpdb \"" + buildAbsPath(data["targetpdb"].str.c_str()).toStdString() + "\"";
  jobtxt=::format(" --targetsite \"%0.3f %0.3f %0.3f\"", data["x"].f,data["y"].f,data["z"].f);
  args += jobtxt;

  // Create the job object
  job = MIMainWindow::instance()->GetJobManager()->CreateJob();
  MIData settings;
  settings["workingDirectory"].str = data["workdir"].str.c_str();
  settings["jobName"].str = "Cocrystal Superpos";
  job->setSettings(settings);
  log=buildAbsPath(::format("%smifit%ld.log", (data["workdir"].str+QDir::separator().toAscii()).c_str(), job->jobId()).c_str()).toStdString();
  jobtxt=::format("python \"%s\" ligandoverlap %s", MIExpertPy().toAscii().constData(), args.c_str());
  job->setLogFile(log.c_str());
  job->setCommandLine(jobtxt.c_str());
  job->StartJob();
}

void Tools::OnSadPhasing() {
  BatchJob* job = NULL;
  std::string jobtxt, log, temp, args;
  if (!VerifyMIExpert() || !VerifyCCP4()) {
    return;
  }

  static MISadPhasingDialog dlg(MIMainWindow::instance(), "SAD Phasing");
  MIData data;
  dlg.GetInitialData(data);
  if (!dlg.GetResults(data)) {
    return;
  }
  // Required things
  args += " --molimagehome \"" + buildAbsPath(Application::instance()->MolimageHome.c_str()).toStdString() + "\"";
  args += " --workdir \"" + buildAbsPath(data["workdir"].str.c_str()).toStdString() + "\"";
  args += " --saddatafile \"" + buildAbsPath(data["saddatafile"].str.c_str()).toStdString() + "\"";
  args += " --sitefile \"" + buildAbsPath(data["sitefile"].str.c_str()).toStdString() + "\"";
  args += " --scatterer \"" + data["scatterer"].str + "\"";
  args += " --sitenumber " + ::format("%d",data["sitenumber"].u);
  args += " --ssnumber " + ::format("%d",data["ssnumber"].u);
  args += " --separation " + ::format("%f",data["separation"].f);
  args += " --solventfraction " + ::format("%f",data["solventfraction"].f);
  args += " --siterefinemethod \"" + data["siterefinemethod"].str  + "\"";
  if (!Application::instance()->ShelxHome.empty()) {
    args += " --shelx_dir \"" + buildAbsPath(Application::instance()->ShelxHome.c_str()).toStdString() + "\"";
  }
  if (data["bothhands"].b) {
    args += " --bothhands yes";
  } else {
    args += " --bothhands no";
  }
  if (data["change_spacegroup"].b ) {
    args += " --spacegroup_no " + ::format("%d",data["spacegroup_no"].radio);
  }

  // Create the job object
  job = MIMainWindow::instance()->GetJobManager()->CreateJob();
  MIData settings;
  settings["workingDirectory"].str = data["workdir"].str;
  settings["jobName"].str = "SAD phasing";
  job->setSettings(settings);
  log=buildAbsPath(::format("%smifit%ld.log", (data["workdir"].str + QDir::separator().toAscii()).c_str(), job->jobId()).c_str()).toStdString();
  jobtxt=::format("python \"%s\" sadphase %s", MIExpertPy().toAscii().constData(), args.c_str());
  job->setLogFile(log.c_str());
  job->setCommandLine(jobtxt.c_str());
  std::string htmlName = "mi_phase_summary.html";
  // TODO open summary in browser
  job->StartJob();
  return;
}

void Tools::OnNCSModeling() {
  BatchJob* job = NULL;
  std::string jobtxt, log, temp, pdbout;
  if (!VerifyMIExpert() || !VerifyCCP4()) {
    return;
  }

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
  pdbout=buildAbsPath(::format("%s%c%s_out.pdb", data["workdir"].str.c_str(), QDir::separator().toAscii(), QFileInfo(doc->GetTitle().c_str()).fileName().toStdString().c_str()).c_str()).toStdString();
  model->SavePDBFile(pdbout.c_str());
  //Build command line
  std::string cmd;
  cmd=::format("ncsmodeler --pdbfile \"%s\" --targetchain \"%s\"", buildAbsPath(pdbout.c_str()).toStdString().c_str(), data["chainid"].str.c_str());
  if (data["model"].str.size()) {
    cmd += " --preserve_model \"" + buildAbsPath(data["model"].str.c_str()).toStdString() + "\"";
  }
  if (data["maskadditions"].str.size()) {
    cmd += " --maskextras \"" + buildAbsPath(data["maskadditions"].str.c_str()).toStdString() + "\"";
  }
  if (data["mtzdata"].str.size()) {
    cmd += " --mtzfile \"" + buildAbsPath(data["mtzdata"].str.c_str()).toStdString() + "\"";
    switch (data["phasecalc"].radio) {
      case 0:
        break;
      case 1:
        cmd += " --phase_prob yes";
        break;
      case 2:
        cmd += " --phase_comb yes";
        break;
    }
  }
  //Copied from superposition script
  // Create the job object
  job = MIMainWindow::instance()->GetJobManager()->CreateJob();
  MIData settings;
  settings["workingDirectory"].str = data["workdir"].str;
  settings["jobName"].str = "NCS Modeling";
  job->setSettings(settings);

  log=buildAbsPath(::format("%smifit%ld.log", (data["workdir"].str+QDir::separator().toAscii()).c_str(), job->jobId()).c_str()).toStdString();
  jobtxt=::format("python \"%s\" %s", MIExpertPy().toAscii().constData(), cmd.c_str());
  job->setLogFile(log.c_str());
  job->setCommandLine(jobtxt.c_str());
  pdbout=buildAbsPath(::format("%s%cncs_%s_out.pdb", data["workdir"].str.c_str(), QDir::separator().toAscii(), QFileInfo(doc->GetTitle().c_str()).fileName().toStdString().c_str()).c_str()).toStdString();

  job->StartJob();
}

void Tools::OnCustom() {
  MICustomJobDialog dlg(0, "Custom Job");
  MIData data;
  dlg.GetInitialData(data);

  MIGLWidget *doc = MIMainWindow::instance()->currentMIGLWidget();
  if (doc != NULL) {
    EMap* map = doc->GetDisplaylist()->GetCurrentMap();
    if (map != NULL) {
      data["dataFile"].str = map->pathName;
    }
    if (data["modelFile"].str.size() == 0) {
      Molecule* model = doc->GetDisplaylist()->CurrentItem();
      if (model != NULL) {
        data["modelFile"].str = model->pathname;
      }
    }
  }
  if (data["workingDirectory"].str.size() == 0) {
    data["workingDirectory"].str = QDir::currentPath().toStdString();
  }
  if (!dlg.GetResults(data)) {
    return;
  }
  CustomJob* job = MIMainWindow::instance()->GetJobManager()->CreateCustomJob();
  job->setSettings(data);
  job->StartJob();
}

void Tools::FillToolsMenu(QMenu* parent) {

    connect(parent, SIGNAL(aboutToShow()),
            this, SLOT(OnUpdateForJobLimit()));

    actions += parent->addAction("&Integrate with d*TREK", this, SLOT(OnIntegrate()));
    actions += parent->addAction("&SAD Phasing", this, SLOT(OnSadPhasing()));

    QMenu* convcif_menu = new QMenu("&Convert CIF to", parent);
    convcif_menu->addAction("&SHELX", this, SLOT(OnCIF2Shellx()));
    convcif_menu->addAction("&CNS / CNX", this, SLOT(OnCIF2CNS()));

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
  static std::string workdir("");
  BatchJob* job = NULL;
  std::string jobtxt, log, temp;
  if (!VerifyMIExpert() || !VerifyCCP4()) {
    return;
  }

  static MIIntegrateDialog dlg(MIMainWindow::instance(), "Integrate");
  MIData data;
  dlg.GetInitialData(data);
  if (!dlg.GetResults(data)) {
    return;
  }

  job = MIMainWindow::instance()->GetJobManager()->CreateJob();
  MIData settings;
  settings["jobName"].str = "Integrate";
  job->setSettings(settings);

  std::string commonArguments;
  if (data["detector_constants"].str.size()) {
    commonArguments += " --detector_constants=\""+ buildAbsPath(data["detector_constants"].str.c_str()).toStdString() + "\"";
  }
  commonArguments += " --spacegroup=\""+ ::format("%d",data["spacegroup_no"].radio) + "\"";
  if (data["first_image"].u != UINT_MAX) {
    commonArguments += " --first_image=\""+ ::format("%d",data["first_image"].u)+ "\"";
  }
  if (data["last_image"].u != UINT_MAX) {
    commonArguments += " --last_image=\""+ ::format("%d",data["last_image"].u)+ "\"";
  }
  if (data["integrate_resolution"].str.size()) {
    commonArguments += " --integrate_resolution=\""+ data["integrate_resolution"].str + "\"";
  }

  log=buildAbsPath(::format("mifit%ld.log", job->jobId()).c_str()).toStdString();

  std::vector<std::string> templateList = data["template_image"].strList;
  std::vector<std::string>::iterator iter = templateList.begin();
  bool firstCommand = true;
  while (iter != templateList.end()) {
    std::string templateImage = *iter;
    ++iter;
    std::string cmd("integrate");
    cmd += " --template_image=\"" + templateImage + "\"";
    cmd += commonArguments;
    std::string redirect(">>");
    if (firstCommand) {
      redirect = ">";
      firstCommand = false;
    }
    jobtxt=::format("python \"%s\" %s", MIExpertPy().toAscii().constData(), cmd.c_str());
    // TODO rework to run in one command
    job->setCommandLine(jobtxt.c_str());
  }

  job->setLogFile(log.c_str());
  job->StartJob();
}

void Tools::OnRunTestJob() {
  // set the working directory
#ifdef WIN32
  MIMainWindow::instance()->GetJobManager()->SetWorkDirectory("C:\\temp");
#endif
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
