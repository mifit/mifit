#include <vector>

#include "nonguilib.h"
#include "mathlib.h"
#include "chemlib.h"
#include "RESIDUE_.h"
#include "conflib.h"
#include "maplib.h"
#include "utillib.h"
#include "jobslib.h"

#include "tools.h"
#include "macafxwin.h"
#include "id.h"
#include "Application.h"
#include "uitest.h"
#include "MIDialog.h"
#include "MIMainWindow.h"
#include "MIEventHandlerMacros.h"
#include "MIGLWidget.h"
#include "Displaylist.h"
#include "EMap.h"
#include "molw.h"

using namespace chemlib;

#include <stdarg.h>


#include <QFile>
#include <QProcess>
#include <QString>
#include <QFileInfo>
#include <QDir>

Tools *Tools::_instance = NULL;

// convert path to system-appropriate string:
static std::string buildPath(const char* fmt, ...) {
  va_list args;
  va_start(args, fmt);
  std::string path = format_arg_list(fmt, args);
  va_end(args);

  path=QDir::toNativeSeparators(path.c_str()).toStdString();
  return path;
}

// convert path to system-appropriate string:
static std::string buildPath(const std::string &path) {
  if (path.size()==0)
    return path;
  if (path==std::string("none"))
    return "none";

  return QDir::toNativeSeparators(path.c_str()).toStdString();
}

// convert path to system-appropriate absolute path string:
static std::string buildAbsPath(const std::string &path) {
  if (path.size()==0)
    return path;
  if (path==std::string("none"))
    return "none";

  return QDir::toNativeSeparators(
    QDir().absoluteFilePath(path.c_str())).toStdString();
}

// convert path to system-appropriate absolute path string:
static std::string buildAbsPath(const char* fmt, ...) {
  va_list args;
  va_start(args, fmt);
  std::string path = format_arg_list(fmt, args);
  va_end(args);

  return buildAbsPath(path);
}



static const char *MIExpertCommand()
{
  static std::string cmd; // static so that returned char * will stay valid

  // must be re-assigned on each call b/c MolimageHome might change
  cmd=Application::instance()->GetMolimageHome() + "/MIExpert/MIExpert.py";
  cmd=QDir::toNativeSeparators(cmd.c_str()).toStdString();
  return cmd.c_str();
}

bool Tools::VerifyMIExpert() {
  if (QFile(MIExpertCommand()).exists()) {
    return true;
  }
  MIMessageBox("Cannot find MIExpert", "Error", MIDIALOG_ICON_ERROR);
  return false;
}

bool Tools::VerifyCCP4() {
  static bool firsttime=true;
  static bool result=false;
  if (!firsttime) {
    return result;
  }
  firsttime=false;

  //NOTE: mtzdump used to be the test program, but with the port to Qt
  //      it's been replaced by pdbset.  mtzdump does not terminate without
  //      additional input, but pdbset does, which makes it easier to test.
  QByteArray pdbsetOutput;
  QProcess pdbsetProcess;
  pdbsetProcess.start("pdbset");

  pdbsetProcess.waitForFinished(2000);

  if (pdbsetProcess.exitStatus() != QProcess::NormalExit) {
    pdbsetProcess.kill();
    MIMessageBox("Cannot find CCP4\n(Unable to run pdbset)", "Error", MIDIALOG_ICON_ERROR);
    result=false;
    return false;
  }

  QString outputText(pdbsetProcess.readAllStandardOutput());
  if (outputText.indexOf("PDBSET") == -1)  {
    pdbsetProcess.kill();
    MIMessageBox("Cannot find CCP4\n(Unable to identify output as from pdbset)", "Error", MIDIALOG_ICON_ERROR);
    result=false;
    return false;
  }
  pdbsetProcess.kill();
  result=true;
  return true;
}

void Tools::OnBindNGrind() {
  BatchJob* job;
  std::string s, log, args;
  if (!VerifyMIExpert() || !VerifyCCP4()) {
    return;
  }

  MIBindNGrindDialog dlg(0, "Cocrystal Solution");
  MIData data;
  dlg.GetInitialData(data);
  if (!dlg.GetResults(data)) {
    return;
  }

  // Write config file
  for (unsigned int i=0; i < data["hklin"].strList.size(); ++i ) {
    args += " --hklin \"" + buildAbsPath(data["hklin"].strList[i]) + "\"";
    args += " --workdir none";
  }
  args += " --pdbin \"" + buildAbsPath(data["pdbin"].str) + "\"";
  args += " --process_engine " + data["process_engine"].str;

  if (data["arpwarpmap"].b) {
    args += " --arpwarpmap yes";
  } else {
    args += " --arpwarpmap no";
  }

  if (data["detector_constants"].str.size()) {
    args += " --detector_constants \""+ data["detector_constants"].str + "\"";
  }

  if (data["spacegroup_no"].radio != 0) {
    args += " --spacegroup_no " + ::format("%d",data["spacegroup_no"].radio);
  }
  args += " --reference_mtz \"" + buildAbsPath(data["reference_mtz"].str) + "\"";

  if (!data["multi_search"].b) {
    args += " --multi_search no";
  } else {
    args += " --multi_search yes";
  }

  args += " --libfile \"" + buildAbsPath(data["libfile"].str) + "\"";

  args += data["viewpoint_method"].str;
  if (data["bngsummary"].str.size()) {
    args += " --bngsummary \""+ buildAbsPath(data["bngsummary"].str) + "\"";
  }
  args += " --mifit no ";

  job = MIMainWindow::instance()->GetJobManager()->CreateJob();
  MIData settings;
  settings["jobName"].str = "Cocrystal Solution";
  job->setSettings(settings);

  std::string lig("");
  /* Handle Place Ligand checkbox */
  if (data["place_ligand"].b && data["viewpoint_method"].str.size()) {
    if (data["ligand_from_dictionary"].b) {
      RESIDUE* dictres = MIFitDictionary()->GetDictResidue(data["ligand_name"].str.c_str(), 0);
      RESIDUE* reslist = new RESIDUE(*dictres);
      Molecule model(reslist, "Dictionary", NULL, NULL, 0, MoleculeType::Other);
      lig = "dictlig.pdb";
      FILE* ligOut = fopen(lig.c_str(), "w");
      chemlib::PDB fpdb;
      MIMolInfo mi;
      mi.res = model.getResidues();
      mi.bonds = model.getBonds();
      fpdb.Write(ligOut, mi);
      fclose(ligOut);
      // delete reslist; // NOTE: this would be a double delete, as the Molecule dtor deletes this
      job->AddtoCleanup(lig.c_str());
    } else {
      lig=data["ligand_name"].str;
    }
    args += " --fragfit \""+ lig + "\"";
  }
  args += " --molimagehome \""+ buildAbsPath(Application::instance()->MolimageHome.c_str()) + "\"";

  s=::format("New Cocrystal Solution job #%ld", job->JobId);
  Logger::log(s);
  s=::format("mi_runbng.txt");
  job->AddtoCleanup(s.c_str());
  QFileInfo workdir(data["hklin"].strList[0].c_str());
  log=buildAbsPath("%smifit%ld.log", (workdir.absolutePath().toStdString()+QDir::separator().toAscii()).c_str(), job->JobId);
  s=::format("\"%s\" bng %s > \"%s\"", MIExpertCommand(), args.c_str(), log.c_str() );
  job->LogFile = log;
  job->AddtoCleanup(log.c_str());
  if (!job->WriteCommand(s.c_str())) {
    Logger::message("Job failed to initialize!");
    MIMainWindow::instance()->GetJobManager()->DeleteJob(job);
  }
  if (data["bngsummary"].str.size()) {
    QDir dir(data["bngsummary"].str.c_str());
    if (!dir.exists()) {
      MIMessageBox("HTML Summary Directory doesn't exist!\nAborting");
      return;
    }
    QStringList filters;
    filters << "bng_jobsummary*.htm";
    QFileInfoList fi=dir.entryInfoList(filters);
    int count=fi.count() +1;
    std::string htmlName=::format("bng_jobsummary_%d.htm", count);


#ifndef _WIN32 //Unix pathnames are needed
    s=::format("%s \"%s/%s\"", buildPath(Application::instance()->HTMLBrowser).c_str(),
               buildAbsPath(data["bngsummary"].str).c_str(),
               buildPath(htmlName.c_str()));
#else
    s=::format("start \"%s\" \"%s\\%s\"",
               buildPath(Application::instance()->HTMLBrowser).c_str(),
               buildAbsPath(data["bngsummary"].str).c_str(),
               buildPath(htmlName.c_str()));
#endif
    job->WriteCommand(s.c_str());
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
  s=::format("New CIF restraints job #%ld", job->JobId);
  Logger::log(s);

  QFileInfo workdir(filename.c_str());
  log=::format("%s%cmifit%ld.log", workdir.absolutePath().toStdString().c_str(),QDir::separator().toAscii(), job->JobId);

  s=::format("\"%s\" convertlib --cif \"%s\" --workdir \"%s\" --refprogram \"%s\" > \"%s\"",
             MIExpertCommand(),
             filename.c_str(),
             workdir.absolutePath().toStdString().c_str(),
             format,
             log.c_str() );
  job->LogFile = log;
  job->AddtoCleanup(log.c_str());
  if (!job->WriteCommand(s.c_str())) {
    Logger::message("Job failed to initialize!");
    MIMainWindow::instance()->GetJobManager()->DeleteJob(job);
  }
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
  std::string multi = "no";
  std::string match = "no";
  std::string fixed = "none";
  std::string jobtxt;
  std::string args;
  std::string log;

  // added to save working directory between commands and to use current working directory - dem
  if (!VerifyMIExpert() || !VerifyCCP4()) {
    return;
  }

  MIMolRepDialog dlg(0, "Molecular Replacement");
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
  args +=" --pdbfile \"" + buildAbsPath(data["model"].str)  + "\"" +
    " --mtzfile \"" + buildAbsPath(data["mtzfile"].str)  + "\"" +
    " --workdir \"" + buildAbsPath(data["workdir"].str) + "\"" +
    " --multi_search " + (data["multi_search"].b ? "yes": "no") +
    " --match_pdbin " + (data["match_pdbin"].b ? "yes": "no") +
    " --copies " + ::format("%d",data["copies"].u);
  if (data["fixed_pdb"].str.size())
    args += " --fixed_pdb \"" + buildAbsPath(data["fixed_pdb"].str)  + "\" ";

  job = MIMainWindow::instance()->GetJobManager()->CreateJob();
  MIData settings;
  settings["workingDirectory"].str = data["workdir"].str;
  settings["jobName"].str = "Molrep";
  job->setSettings(settings);
  log=buildAbsPath("%s%cmifit%ld.log", data["workdir"].str.c_str(),QDir::separator().toAscii(), job->JobId);
  jobtxt=::format("\"%s\" molrep %s > \"%s\"", MIExpertCommand(), args.c_str(), log.c_str());
  job->LogFile = log;
  job->AddtoCleanup(log.c_str());
  if (!job->WriteCommand(jobtxt.c_str())) {
    Logger::message("Job failed to initialize!");
    MIMainWindow::instance()->GetJobManager()->DeleteJob(job);
  }
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
  s=::format("New CIF to ShellX job #%ld", job->JobId);
  Logger::log(s);
  QFileInfo workdir(filename.c_str());
  log=::format("%s%cmifit%ld.log", workdir.absolutePath().toStdString().c_str(),QDir::separator().toAscii(), job->JobId);

  s=::format("\"%s\" restraints --pdbfile \"%s\" --workdir \"%s\" > \"%s\"",
             MIExpertCommand(), filename.c_str(), workdir.absolutePath().toStdString().c_str(), log.c_str() );
  job->LogFile = log;
  job->AddtoCleanup(log.c_str());
  if (!job->WriteCommand(s.c_str())) {
    Logger::message("Job failed to initialize!");
    MIMainWindow::instance()->GetJobManager()->DeleteJob(job);
  }
  job->StartJob();
}

void Tools::OnRefine() {
  static std::string workdir;
  BatchJob* job = NULL;
  std::string jobtxt, log, args;
  if (!VerifyMIExpert() || !VerifyCCP4()) {
    return;
  }

  MIRefinementDialog dlg(0, "Refinement");
  MIData data;
  dlg.GetInitialData(data);
  if (!dlg.GetResults(data)) {
    return;
  }

  //Get the required options
  args +=" --mifithome \"" + buildAbsPath(Application::instance()->MolimageHome.c_str()) + "\"" +
    " --workdir \"" + buildAbsPath(data["workdir"].str) + "\"" +
    " --pdbfile \"" + buildAbsPath(data["pdbfile"].str) + "\"" +
    " --mtzfile \"" + buildAbsPath(data["mtzfile"].str) + "\"" +
    " --weight " + ::format("%f",data["weight"].f) +
    " --cycles " + ::format("%d",data["cycles"].u);
  if (!Application::instance()->ShelxHome.empty()) {
    args += " --shelx_dir \"" + buildAbsPath(Application::instance()->ShelxHome) + "\"";
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
    args += " --libfile \"" + buildAbsPath(data["libfile"].str) + "\"";
  }
  if (data["tls_file"].str.size()) {
    args += " --tls_file \"" + buildAbsPath(data["tls_file"].str) + "\"";
  }

  job = MIMainWindow::instance()->GetJobManager()->CreateJob();
  MIData settings;
  settings["workingDirectory"].str = data["workdir"].str;
  settings["jobName"].str = "Refinement";
  job->setSettings(settings);
  log=buildAbsPath("%smifit%ld.log", (data["workdir"].str+QDir::separator().toAscii()).c_str(), job->JobId);
  jobtxt=::format("\"%s\" refine %s > \"%s\"", MIExpertCommand(), args.c_str(), log.c_str() );
  job->LogFile = log;
  job->AddtoCleanup(log.c_str());
  if (!job->WriteCommand(jobtxt.c_str())) {
    Logger::message("Job failed to initialize!");
    MIMainWindow::instance()->GetJobManager()->DeleteJob(job);
  }
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
#ifndef _WIN32 //Unix pathnames are needed
    jobtxt=::format("%s \"%s/%s\"",
                    buildPath(Application::instance()->HTMLBrowser).c_str(),
                    buildAbsPath(data["workdir"].str).c_str(),
                    buildPath(errorName).c_str());
#else
    jobtxt=::format("start \"%s\" \"%s\\%s\"",
                    buildPath(Application::instance()->HTMLBrowser).c_str(),
                    buildAbsPath(data["workdir"].str).c_str(),
                    buildPath(errorName).c_str());
#endif
    job->WriteCommand(jobtxt.c_str());
  }
  job->StartJob();
}

void Tools::OnJobReport() {
  BatchJob* job = NULL;
  std::string jobtxt, log, temp, args;
  if (!VerifyMIExpert() || !VerifyCCP4()) {
    return;
  }

  MIJobReportDialog dlg(0, "Job Report");
  MIData data;
  dlg.GetInitialData(data);
  if (!dlg.GetResults(data)) {
    return;
  }

  // Required things
  args += " --workdir \"" + buildAbsPath(data["workdir"].str) + "\"";
  args += " --mtzfile \"" + buildAbsPath(data["mtzfile"].str) + "\"";
  args += " --pdbfile \"" + buildAbsPath(data["pdbfile"].str) + "\"";
  args += " --molimagehome \"" + buildAbsPath(Application::instance()->MolimageHome.c_str()) + "\"";
  args += " --libfile \"" + buildAbsPath(data["libfile"].str) + "\"";
  args += " --seqfile \"" + buildAbsPath(data["seqfile"].str) + "\"";
  args += " --templatefile \"" + buildAbsPath(data["templatefile"].str) + "\"";
  args += " --datalogfile \"" + buildAbsPath(data["datalogfile"].str) + "\"";

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
      args += " --image1 \"" + buildAbsPath(data["image1"].str) + "\"";
    }
    if (data["image2"].str.size()) {
      args += " --image2 \"" + buildAbsPath(data["image2"].str) + "\"";
    }
    if (data["image3"].str.size()) {
      args += " --image3 \"" + buildAbsPath(data["image3"].str) + "\"";
    }
  }
  if (data["rootname"].str.size()) {
    args += " --rootname \"" + buildAbsPath(data["rootname"].str) + "\"";
  }

  // Create the job object
  job = MIMainWindow::instance()->GetJobManager()->CreateJob();
  MIData settings;
  settings["workingDirectory"].str = data["workdir"].str.c_str();
  settings["jobName"].str = "Report";
  job->setSettings(settings);
  log=buildAbsPath("%smifit%ld.log", (data["workdir"].str+QDir::separator().toAscii()).c_str(), job->JobId);
  jobtxt=::format("\"%s\" deposit3d %s > \"%s\"", MIExpertCommand(), args.c_str(), log.c_str());

  job->LogFile = log;
  job->AddtoCleanup(log.c_str());
  if (!job->WriteCommand(jobtxt.c_str())) {
    Logger::message("Job failed to initialize!");
    MIMainWindow::instance()->GetJobManager()->DeleteJob(job);
  }
  if (data["html_report"].b) {
    std::string htmlName;
    if (data["rootname"].str.size()) {
      htmlName=::format("%s.htm", data["rootname"].str.c_str());
    } else {
      htmlName=::format("pdbdeposit.htm");
    }
#ifndef _WIN32 //Unix pathnames are needed
    jobtxt=::format("%s \"%s/%s\"",
                    buildPath(Application::instance()->HTMLBrowser).c_str(),
                    buildAbsPath(data["workdir"].str).c_str(),
                    buildPath(htmlName).c_str());
#else
    jobtxt=::format("start \"%s\" \"%s\\%s\"",
                    buildPath(Application::instance()->HTMLBrowser).c_str(),
                    buildAbsPath(data["workdir"].str).c_str(),
                    buildPath(htmlName).c_str());
#endif
    job->WriteCommand(jobtxt.c_str());
  }
  job->StartJob();
}

void Tools::OnCoCrystalSuperPos() {
  BatchJob* job = NULL;
  std::string jobtxt, log, temp, args;
  if (!VerifyMIExpert() || !VerifyCCP4()) {
    return;
  }

  MICocrystalSuperpositionDialog dlg(0, "Cocrystal superposition");
  MIData data;
  dlg.GetInitialData(data);
  if (!dlg.GetResults(data)) {
    return;
  }

  // Required things
  args += " --workdir \"" + buildAbsPath(data["workdir"].str) + "\"";
  args += " --pdbdir \"" + buildAbsPath(data["pdbdir"].str) + "\"";
  args += " --targetpdb \"" + buildAbsPath(data["targetpdb"].str) + "\"";
  jobtxt=::format(" --targetsite \"%0.3f %0.3f %0.3f\"", data["x"].f,data["y"].f,data["z"].f);
  args += jobtxt;

  // Create the job object
  job = MIMainWindow::instance()->GetJobManager()->CreateJob();
  MIData settings;
  settings["workingDirectory"].str = data["workdir"].str.c_str();
  settings["jobName"].str = "Cocrystal Superpos";
  job->setSettings(settings);
  log=buildAbsPath("%smifit%ld.log", (data["workdir"].str+QDir::separator().toAscii()).c_str(), job->JobId);
  jobtxt=::format("\"%s\" ligandoverlap %s > \"%s\"", MIExpertCommand(), args.c_str(), log.c_str());
  job->LogFile = log;
  job->AddtoCleanup(log.c_str());
  if (!job->WriteCommand(jobtxt.c_str())) {
    Logger::message("Job failed to initialize!");
    MIMainWindow::instance()->GetJobManager()->DeleteJob(job);
  }

#ifndef _WIN32 //Unix pathnames are needed
  jobtxt=::format("echo loadpdb 0 %s/allligands.pdb > \"%s\"",
                  buildAbsPath(data["workdir"].str).c_str(),
                  job->UpdateScript.c_str());
  job->WriteCommand(jobtxt.c_str());
#else
  jobtxt=::format("echo loadpdb 0 %s\\allligands.pdb > \"%s\"",
                  buildAbsPath(data["workdir"].str).c_str(),
                  job->UpdateScript.c_str());
  job->WriteCommand(jobtxt.c_str());
#endif
  job->SetDocument(MIMainWindow::instance()->currentMIGLWidget());
  job->StartJob();
}

void Tools::OnSadPhasing() {
  BatchJob* job = NULL;
  std::string jobtxt, log, temp, args;
  if (!VerifyMIExpert() || !VerifyCCP4()) {
    return;
  }

  MISadPhasingDialog dlg(0, "SAD Phasing");
  MIData data;
  dlg.GetInitialData(data);
  if (!dlg.GetResults(data)) {
    return;
  }
  // Required things
  args += " --molimagehome \"" + buildAbsPath(Application::instance()->MolimageHome.c_str()) + "\"";
  args += " --workdir \"" + buildAbsPath(data["workdir"].str) + "\"";
  args += " --saddatafile \"" + buildAbsPath(data["saddatafile"].str) + "\"";
  args += " --sitefile \"" + buildAbsPath(data["sitefile"].str) + "\"";
  args += " --scatterer \"" + data["scatterer"].str + "\"";
  args += " --sitenumber " + ::format("%d",data["sitenumber"].u);
  args += " --ssnumber " + ::format("%d",data["ssnumber"].u);
  args += " --separation " + ::format("%f",data["separation"].f);
  args += " --solventfraction " + ::format("%f",data["solventfraction"].f);
  args += " --siterefinemethod \"" + data["siterefinemethod"].str  + "\"";
  if (!Application::instance()->ShelxHome.empty()) {
    args += " --shelx_dir \"" + buildAbsPath(Application::instance()->ShelxHome) + "\"";
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
  log=buildAbsPath("%smifit%ld.log", (data["workdir"].str + QDir::separator().toAscii()).c_str(), job->JobId);
  jobtxt=::format("\"%s\" sadphase %s> \"%s\"", MIExpertCommand(), args.c_str(), log.c_str());
  job->LogFile = log;
  job->AddtoCleanup(log.c_str());
  if (!job->WriteCommand(jobtxt.c_str())) {
    Logger::message("Job failed to initialize!");
    MIMainWindow::instance()->GetJobManager()->DeleteJob(job);
  }
  std::string htmlName = "mi_phase_summary.html";
#ifndef _WIN32 //Unix pathnames are needed
  jobtxt=::format("%s \"%s/%s\"",
                  buildPath(Application::instance()->HTMLBrowser).c_str(),
                  buildAbsPath(data["workdir"].str).c_str(),
                  buildPath(htmlName).c_str());
#else
  jobtxt=::format("start \"%s\" \"%s\\%s\"",
                  buildPath(Application::instance()->HTMLBrowser).c_str(),
                  buildAbsPath(data["workdir"].str).c_str(),
                  buildPath(htmlName).c_str());
#endif
  job->WriteCommand(jobtxt.c_str());
  job->StartJob();
  return;
}

void Tools::OnNCSModeling() {
  BatchJob* job = NULL;
  std::string jobtxt, log, temp, pdbout;
  if (!VerifyMIExpert() || !VerifyCCP4()) {
    return;
  }

  MINCSModelingDialog dlg(0, "NCS Modeling");
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
  pdbout=buildAbsPath("%s%c%s_out.pdb", data["workdir"].str.c_str(), QDir::separator().toAscii(), QFileInfo(doc->GetTitle().c_str()).fileName().toStdString().c_str());
  model->SavePDBFile(pdbout.c_str());
  //Build command line
  std::string cmd;
  cmd=::format("ncsmodeler --pdbfile \"%s\" --targetchain \"%s\"", buildAbsPath(pdbout).c_str(), data["chainid"].str.c_str());
  if (data["model"].str.size()) {
    cmd += " --preserve_model \"" + buildAbsPath(data["model"].str) + "\"";
  }
  if (data["maskadditions"].str.size()) {
    cmd += " --maskextras \"" + buildAbsPath(data["maskadditions"].str) + "\"";
  }
  if (data["mtzdata"].str.size()) {
    cmd += " --mtzfile \"" + buildAbsPath(data["mtzdata"].str) + "\"";
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

  log=buildAbsPath("%smifit%ld.log", (data["workdir"].str+QDir::separator().toAscii()).c_str(), job->JobId);
  jobtxt=::format("\"%s\" %s > \"%s\"", MIExpertCommand(), cmd.c_str(), log.c_str());
  job->LogFile = log;
  job->AddtoCleanup(log.c_str());
  if (!job->WriteCommand(jobtxt.c_str())) {
    Logger::message("Job failed to initialize!");
    MIMainWindow::instance()->GetJobManager()->DeleteJob(job);
  }
  pdbout=buildAbsPath("%s%cncs_%s_out.pdb", data["workdir"].str.c_str(), QDir::separator().toAscii(), QFileInfo(doc->GetTitle().c_str()).fileName().toStdString().c_str());

  jobtxt=::format("echo loadpdb 0 %s > \"%s\"", pdbout.c_str(), job->UpdateScript.c_str());
  job->WriteCommand(jobtxt.c_str());
  job->SetDocument(MIMainWindow::instance()->currentMIGLWidget());
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

void Tools::FillToolsMenu(MIMenu* parent, bool havedoc) {

  MIMenu* convcif_menu = new MIMenu(*this);
  parent->Append(ID_TOOLS_INTEGRATE, "&Integrate with d*TREK");
  parent->Append(ID_TOOLS_SAD, "&SAD Phasing");
  parent->Append(0, "&Convert CIF to", convcif_menu);
  convcif_menu->Append(ID_TOOLS_CIF2SHELLX, "&SHELX");
  convcif_menu->Append(ID_TOOLS_CIF2CNS, "&CNS / CNX");
  parent->Append(ID_TOOLS_REFMACRES, "S&et Refmac5 restraints");
  parent->Append(ID_TOOLS_MOLREP, "&Molecular Replacement");
  parent->Append(ID_TOOLS_REFINE, "&Refinement");
  parent->Append(ID_TOOLS_BNG, "C&ocrystal Solution");
  parent->Append(ID_TOOLS_REPORT, "Re&port");

  parent->Append(ID_TOOLS_CUSTOM, "Run Custom Job");

  // Things that are in the job menu
  if (havedoc) {
    parent->AppendSeparator();
    parent->Append(ID_TOOLS_NCS, "&NCS Modeling");
    parent->Append(ID_TOOLS_SUPER, "Cocr&ystal Superposition");
#ifdef DEBUG
    parent->Append(ID_RUN_TESTJOB, "Run &Test Job...", "Run a test job that waits 5 seconds then terminates");
#endif
  }
}

void Tools::OnUpdateForJobLimit(const MIUpdateEvent& pCmdUI) {
  pCmdUI.Enable(!MIMainWindow::instance()->isJobLimit());

  // disable menu options which depend on model window if not present
  bool havedoc=MIMainWindow::instance()->currentMIGLWidget() != NULL;
  switch (pCmdUI.GetId()) {
    case ID_TOOLS_NCS:
    case ID_TOOLS_SUPER:
    case ID_RUN_TESTJOB:
      pCmdUI.Enable(havedoc);
    default:
      break;
  }
}

void Tools::OnIntegrate() {
  static std::string workdir("");
  BatchJob* job = NULL;
  std::string jobtxt, log, temp;
  if (!VerifyMIExpert() || !VerifyCCP4()) {
    return;
  }

  MIIntegrateDialog dlg(0, "Integrate");
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
    commonArguments += " --detector_constants=\""+ buildAbsPath(data["detector_constants"].str) + "\"";
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

  log=buildAbsPath("mifit%ld.log", job->JobId);

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
    jobtxt=::format("\"%s\" %s %s \"%s\"", MIExpertCommand(), cmd.c_str(), redirect.c_str(), log.c_str() );
    if (!job->WriteCommand(jobtxt.c_str())) {
      Logger::message("Job failed to initialize!");
      MIMainWindow::instance()->GetJobManager()->DeleteJob(job);
    }
  }

  job->LogFile = log;
  job->AddtoCleanup(log.c_str());
  job->StartJob();
}

void Tools::OnRunTestJob() {
  // set the working directory
#ifdef WIN32
  MIMainWindow::instance()->GetJobManager()->SetWorkDirectory("C:\\temp");
#endif
  BatchJob* job =   MIMainWindow::instance()->GetJobManager()->CreateJob();
  job->SetDocument(MIMainWindow::instance()->currentMIGLWidget());
  try {
    TestJob testJob;
    testJob.StartJob(job);
  }
  catch (...) {
    Logger::message("Job Failed!!");
    return;
  }
  WaitCursor* wait = new WaitCursor("Waiting for test job to finish");
  while (wait->CheckForAbort() == false && job->IsRunning()) {
#ifdef _WIN32
    Sleep(100);
#else
    usleep(100000);
#endif
  }
  delete wait;
}


Tools::Tools() : QObject(0),MIEventHandler(this) {
  _instance = this;

BEGIN_EVENT_TABLE(this,none)
EVT_MENU(ID_TOOLS_CIF2SHELLX, Tools::OnCIF2Shellx)
EVT_MENU(ID_TOOLS_CIF2CNS, Tools::OnCIF2CNS)
EVT_MENU(ID_TOOLS_MOLREP, Tools::OnMolRep)
EVT_UPDATE_UI(ID_TOOLS_MOLREP, Tools::OnUpdateForJobLimit)
EVT_MENU(ID_TOOLS_REFMACRES, Tools::OnRefmacRestraints)
EVT_UPDATE_UI(ID_TOOLS_REFMACRES, Tools::OnUpdateForJobLimit)
EVT_MENU(ID_TOOLS_REFINE, Tools::OnRefine)
EVT_UPDATE_UI(ID_TOOLS_REFINE, Tools::OnUpdateForJobLimit)
EVT_MENU(ID_TOOLS_BNG, Tools::OnBindNGrind)
EVT_UPDATE_UI(ID_TOOLS_BNG, Tools::OnUpdateForJobLimit)
EVT_MENU(ID_TOOLS_REPORT, Tools::OnJobReport)
EVT_UPDATE_UI(ID_TOOLS_REPORT, Tools::OnUpdateForJobLimit)
EVT_MENU(ID_TOOLS_SUPER, Tools::OnCoCrystalSuperPos)
EVT_UPDATE_UI(ID_TOOLS_SUPER, Tools::OnUpdateForJobLimit)
EVT_MENU(ID_TOOLS_INTEGRATE, Tools::OnIntegrate)
EVT_UPDATE_UI(ID_TOOLS_INTEGRATE, Tools::OnUpdateForJobLimit)
EVT_MENU(ID_TOOLS_SAD, Tools::OnSadPhasing)
EVT_UPDATE_UI(ID_TOOLS_SAD, Tools::OnUpdateForJobLimit)
EVT_MENU(ID_TOOLS_NCS, Tools::OnNCSModeling)
EVT_UPDATE_UI(ID_TOOLS_NCS, Tools::OnUpdateForJobLimit)
EVT_MENU(ID_TOOLS_CUSTOM, Tools::OnCustom)
EVT_UPDATE_UI(ID_TOOLS_CUSTOM, Tools::OnUpdateForJobLimit)
EVT_MENU(ID_RUN_TESTJOB, Tools::OnRunTestJob)
EVT_UPDATE_UI(ID_RUN_TESTJOB, Tools::OnUpdateForJobLimit)
END_EVENT_TABLE()

}

Tools *Tools::instance() {
  if (_instance)
    return _instance;
  new Tools(); // sets _instance
  return _instance;
}
