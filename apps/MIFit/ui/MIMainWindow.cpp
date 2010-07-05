#include <errno.h>
#include <sys/stat.h>
#include <sys/types.h>
#include <ctype.h> // for alanum
#include <QtGui>
#include <core/Version.h>
#include <chemlib/chemlib.h>
#include <chemlib/Residue.h>
#include <conflib/conflib.h>
#include <jobs/jobslib.h>
#include <util/utillib.h>
#include <core/corelib.h>
#include <ligand/ligandlib.h>
#include <map/maplib.h>
#include <nongui/nonguilib.h>

#include "JobsView.h"
#include "DictEditCanvas.h"
#include "DisplayView.h"
#include "EMap.h"
#include "Displaylist.h"
#include "ManageCrystals.h"
#include "ModelsView.h"
#include "MIMolIO.h"
#include "molw.h"
//#include "preferences/MapPreferencesPanel.h"
//#include "preferences/EnvironmentPreferencesPanel.h"
//#include "preferences/GeneralPreferencesPanel.h"
//#include "preferences/PreferencesDialog.h"
#include "PreferencesDialog.h"
#include "MIMainWindow.h"
#include "MIGLWidget.h"
#include "Application.h"
#include "id.h"
#include "surf.h"
#include "MIMolIO.h"
#include "molw.h"
#include "GLOverviewCanvas.h"
#include "RamaPlot.h"
#include "DictEditDialog.h"
#include "GLFormatDialog.h"
#include "tools.h"
#include "GenericDataDialog.h"
#include "ui/SmilesDialog.h"
#include "ui/AtomColors.h"
#include "ui/BValueColors.h"
#include "CurrentMIGLWidgetAction.h"
#include "script/LocalSocketScript.h"

#ifdef _WIN32
#include <images/mifit_icon_32x32.xpm>
#else
#include <images/mifit_icon.xpm>
#endif

#include <images/annotate.xpm>
#include <images/apply.xpm>
#include <images/cancel.xpm>
#include <images/center.xpm>
#include <images/chekpoin.xpm>
#include <images/colortool.xpm>
#include <images/copy.xpm>
#include <images/decrper.xpm>
#include <images/fit.xpm>
#include <images/hidetool.xpm>
#include <images/incrper.xpm>
#include <images/new.xpm>
#include <images/open.xpm>
#include <images/print.xpm>
#include <images/rotate.xpm>
#include <images/roty90.xpm>
#include <images/rotym90.xpm>
#include <images/save.xpm>
#include <images/showtool.xpm>
#include <images/sidetool.xpm>
#include <images/slabin.xpm>
#include <images/slabout.xpm>
#include <images/topview.xpm>
#include <images/torsion.xpm>
#include <images/translate.xpm>
#include <images/zoomin.xpm>
#include <images/zoomout.xpm>
//#include <images/annotateDisabled.xpm>
//#include <images/applyDisabled.xpm>
//#include <images/cancelDisabled.xpm>
//#include <images/centerDisabled.xpm>
//#include <images/chekpoinDisabled.xpm>
//#include <images/colortoolDisabled.xpm>
//#include <images/copyDisabled.xpm>
//#include <images/decrperDisabled.xpm>
//#include <images/fitDisabled.xpm>
//#include <images/hidetoolDisabled.xpm>
//#include <images/incrperDisabled.xpm>
//#include <images/newDisabled.xpm>
//#include <images/openDisabled.xpm>
//#include <images/printDisabled.xpm>
//#include <images/rotateDisabled.xpm>
//#include <images/roty90Disabled.xpm>
//#include <images/rotym90Disabled.xpm>
//#include <images/saveDisabled.xpm>
//#include <images/showtoolDisabled.xpm>
//#include <images/sidetoolDisabled.xpm>
//#include <images/slabinDisabled.xpm>
//#include <images/slaboutDisabled.xpm>
//#include <images/topviewDisabled.xpm>
//#include <images/torsionDisabled.xpm>
//#include <images/translateDisabled.xpm>
//#include <images/zoominDisabled.xpm>
//#include <images/zoomoutDisabled.xpm>

static bool GetOpenFilenames(std::vector<std::string> &fnames);
static void SplitPath(const std::string &origPath,
                      std::string *dirname,
                      std::string *fname,
                      std::string *ext);

using namespace chemlib;

template <typename T>
static T firstChildWithName(const QWidget *widget, const QString &name)
{
    return widget->findChildren<T>(name).first();
}

MIMainWindow*MIMainWindow::_instance = NULL;

MIMainWindow::MIMainWindow()
    : modelsDock(NULL),
      displayDock(NULL),
      jobsDock(NULL),
      logDock(NULL),
      ramaDock(NULL),
      logWindow(NULL),
      solidSurfMenuAction(NULL),
      _scriptSocket(NULL)
{

    _instance = this;
    QApplication::setWindowIcon(QIcon(QPixmap(mifit_icon)));

    mdiArea = new QMdiArea;
    setCentralWidget(mdiArea);
    connect(mdiArea, SIGNAL(subWindowActivated(QMdiSubWindow*)),
            this, SLOT(updateMenus()));
    connect(mdiArea, SIGNAL(subWindowActivated(QMdiSubWindow*)),
            this, SLOT(childActivated(QMdiSubWindow*)));
    windowMapper = new QSignalMapper(this);
    connect(windowMapper, SIGNAL(mapped(QWidget*)),
            this, SLOT(setActiveSubWindow(QWidget*)));

    tool_bar = addToolBar("MIFit tools");
    tool_bar->setObjectName("MIFit tools");
    tool_bar->setIconSize(QSize(20, 20));

    createActions();
    createMenus();
    createDockWindows();
    createStatusBar();
    updateMenus();

    setWindowTitle(tr("MIFit"));

    LeftFooter("MIFit: Ok", 2500);

    middleFooter = new QLabel(statusBar());
    statusBar()->addPermanentWidget(middleFooter, 0);

    rightFooter = new QLabel(statusBar());
    statusBar()->addPermanentWidget(rightFooter, 0);

    navigator = new GLOverviewCanvas(this);
    navigatorDock = AddAsDockWidget(navigator, "Navigator", Qt::BottomDockWidgetArea);

    memset(cursors, 0, imhCount*sizeof(QCursor*));

    QSettings *settings = MIGetQSettings(); // could use MIConfig (it's the same file), but this api is easier here
    restoreState(settings->value("WindowLayout", saveState()).toByteArray());
    setGeometry(settings->value("WindowGeometry", geometry()).toRect());

    tabifyDockWidget(modelsDock, displayDock);
    tabifyDockWidget(displayDock, jobsDock);
    modelsDock->raise();

    tabifyDockWidget(navigatorDock, ramaDock);
    navigatorDock->raise();

    std::string title = "MIFit ";
    title += MIFit_version;
    setWindowTitle(title.c_str());

    std::string s = ::format("Read in %d scattering factors",
                             MIMapInitializeScatteringFactorTables(Application::instance()->CrystalData.c_str(),
                                                                   Application::instance()->MolimageHome.c_str()));
    Log(s);
}


MIMainWindow::~MIMainWindow()
{
    for (int i = 0; i < imhCount; ++i)
        delete cursors[i];
    MIMapFreeScatteringFactorTables();
    delete JobManager;
}

static QWidget *dying_widget = 0;

void MIMainWindow::OnClose()
{
    dying_widget = currentMIGLWidget();
    mdiArea->closeActiveSubWindow();
    dying_widget = 0;
}

void MIMainWindow::childActivated(QMdiSubWindow *w)
{
    if (!w)
        return;
    MIGLWidget *child = currentMIGLWidget();
    if (child)
    {
        if (child==dying_widget)
            return;
        child->OnActivated();
    }
}

void MIMainWindow::saveLayout()
{
    QSettings *settings = MIGetQSettings(); // could use MIConfig (it's the same file), but this api is easier here
    settings->setValue("WindowLayout", saveState());
    settings->setValue("WindowGeometry", geometry());
}

void MIMainWindow::OnExit()
{
    saveLayout();
    qApp->closeAllWindows();
}

void MIMainWindow::OnAbout()
{
    std::string s = ::format("MIFit %s\n", MIFit_version);
    QMessageBox::about(this, "About MIFit", s.c_str());
}

void MIMainWindow::Log(const std::string &msg)
{
    if (logWindow != NULL)
    {
        logWindow->append(msg.c_str());
    }
}

void MIMainWindow::OnManageCrystals()
{
    ManageCrystals mcd(this);
    mcd.exec();
}

static QFileDialog *fileDialog = NULL;

static void initializeFileDialog()
{
    if (fileDialog == NULL)
    {
        fileDialog = new QFileDialog(0, "Open File(s)", "",
                                     "Recognized files (*.mlw *.pdb *.mtz *.phs *.fcf *.cif *.map *.py);;"
                                     "MIFit session files (*.mlw);;"
                                     "PDB files (*.pdb);;"
                                     "MTZ files (*.mtz);;"
                                     "Phase files (*.phs);;"
                                     "FCF files (*.fcf);;"
//        "Scalepack files (*.sca);;"
                                     "CIF files (*.cif);;"
//        "Reflection files (*.ref);;"
                                     "Map files (*.map);;"
                                     "Python scripts (*.py);;"
                                     "All files (*.*)");
        fileDialog->setFileMode(QFileDialog::ExistingFiles);
    }
}

void MIMainWindow::OnFileOpen()
{
    fileOpen();
}

void MIMainWindow::fileOpen()
{
    std::vector<std::string> fnames;
    if (!GetOpenFilenames(fnames))
        return;
    OpenFiles(fnames);
}


void MIMainWindow::OnFileOpenNew()
{
    std::vector<std::string> fnames;
    if (!GetOpenFilenames(fnames))
        return;
    newWindow();
    OpenFiles(fnames);
}

// unused in Qt version, currently preserved to facilitate merging
//
// CMolwDoc* MainFrame::newCMolwDoc(const wxString& path) {
//   wxDocTemplate* docTemplate = NULL;
//   wxNode* node = m_docManager->GetTemplates().GetFirst();
//   while (node) {
//     wxDocTemplate* t = (wxDocTemplate*) node->GetData();
//     if (t->GetDefaultExtension() == "mlw") {
//       docTemplate = t;
//       break;
//     }
//     node = node->GetNext();
//   }
//   if (docTemplate == NULL) {
//     return NULL;
//   }
//   CMolwDoc* doc = (CMolwDoc*) docTemplate->CreateDocument(path);
//   if (doc != NULL) {
//     doc->SetDocumentName(docTemplate->GetDocumentName());
//     doc->SetDocumentTemplate(docTemplate);
//     doc->Modify(false);
//     if (!doc->OnNewDocument()) {
//       doc->DeleteAllViews();
//     }
//   }
//   return doc;
// }

// void MainFrame::importFiles(const wxArrayString& files) {
//   if (files.GetCount() > 0) {
//     wxString pathname = files[0];

//     CMolwDoc* doc = newCMolwDoc(pathname);
//     if (!doc) {
//       return;
//     }
//     if (!doc->OnOpenDocument(pathname)) {
//       return;
//     }

//     // default name is made by changing extension to .mlw
//     wxString docname(pathname);
//     wxString path, name, ext;
//     wxSplitPath(docname, &path, &name, &ext);
//     name += ".mlw";

//     m_docManager->AddFileToHistory(pathname);
//     doc->SetFilename(name, true);
//     doc->SetTitle(wxFileNameFromPath(name));
//     for (unsigned int i = 1; i < files.GetCount(); ++i) {
//       pathname = files[i];
//       doc->LoadPDBFile(pathname);
//     }
//     doc->Modify(true);
//   }
// }

// void MainFrame::OnTimer(wxTimerEvent&) {
//   timer();
// }

// void MainFrame::timer() {
//   static unsigned int ndocs = 0;
//   wxList& m = m_docManager->GetDocuments();
//   if (m.GetCount() != ndocs) {
//     ndocs = m.GetCount();
//     if (ndocs == 0) {
//       EnableToolBar(0);
//     }
//   }
//   if (MIBusyManager::instance()->Busy()) {
//     MIBusyManager::instance()->SetWaitCursor();
//   }
// }


void MIMainWindow::OnGLFormat()
{
    GLFormatDialog dialog;
    if (currentMIGLWidget())
    {
        dialog.setCurrentFormat(currentMIGLWidget()->format());
    }
    dialog.exec();
}

void MIMainWindow::OnHelp()
{
    QUrl url("http://code.google.com/p/mifit");
    QDesktopServices::openUrl(url);
}

void MIMainWindow::OnPreferences()
{
    showPreferences();
}

void MIMainWindow::showPreferences(int page)
{
    PreferencesDialog dlg(this);
    dlg.setPage(page);
    if (dlg.exec() != QDialog::Accepted)
    {
        return;
    }
    dlg.savePreferences();
    if (Application::instance()->SetDictionary(Application::instance()->XFitDictSetting))
    {
        Application::instance()->LoadDictionary();
    }
    Application::instance()->SetEnv();
    Application::instance()->Write();
}

void MIMainWindow::OnDefineAtomColors()
{
    MIData data;
    AtomColors dlg(this);
    dlg.setWindowTitle("Define Atom Coloring Scheme");
    dlg.InitializeFromData(data);
    if (dlg.exec() != QDialog::Accepted)
    {
        return;
    }
    dlg.GetData(data);
    //  input: atomNames.strList, atomColors.strList
    data["atomNames"].strList = Colors::atomnames;
    data["atomColors"].strList = Colors::atomcolors;

    Colors::atomnames = data["atomNames"].strList;
    Colors::atomcolors = data["atomColors"].strList;

    MIGLWidget *mdoc = currentMIGLWidget();
    if (mdoc)
    {
        Molecule *model = mdoc->GetDisplaylist()->CurrentItem();
        mdoc->ColorModel(model);
        mdoc->ReDraw();
    }
}

void MIMainWindow::OnDefineBValueColors()
{

    MIData data;
    data["color1"].i = Colors::BValueColors[0];
    data["color2"].i = Colors::BValueColors[1];
    data["color3"].i = Colors::BValueColors[2];
    data["color4"].i = Colors::BValueColors[3];
    data["color5"].i = Colors::BValueColors[4];
    data["color6"].i = Colors::BValueColors[5];
    data["color7"].i = Colors::BValueColors[6];
    data["color8"].i = Colors::BValueColors[7];
    data["color9"].i = Colors::BValueColors[8];
    data["color10"].i = Colors::BValueColors[9];

    data["level1"].f = Colors::BValueRanges[0]/100.0f;
    data["level2"].f = Colors::BValueRanges[1]/100.0f;
    data["level3"].f = Colors::BValueRanges[2]/100.0f;
    data["level4"].f = Colors::BValueRanges[3]/100.0f;
    data["level5"].f = Colors::BValueRanges[4]/100.0f;
    data["level6"].f = Colors::BValueRanges[5]/100.0f;
    data["level7"].f = Colors::BValueRanges[6]/100.0f;
    data["level8"].f = Colors::BValueRanges[7]/100.0f;
    data["level9"].f = Colors::BValueRanges[8]/100.0f;
    data["level10"].f = Colors::BValueRanges[9]/100.0f;

    BValueColors dlg(this);
    dlg.setWindowTitle("Define Atom Coloring Scheme");
    dlg.InitializeFromData(data);
    if (dlg.exec() != QDialog::Accepted)
    {
        return;
    }
    dlg.GetData(data);

    Colors::BValueColors[0] = data["color1"].i;
    Colors::BValueColors[1] = data["color2"].i;
    Colors::BValueColors[2] = data["color3"].i;
    Colors::BValueColors[3] = data["color4"].i;
    Colors::BValueColors[4] = data["color5"].i;
    Colors::BValueColors[5] = data["color6"].i;
    Colors::BValueColors[6] = data["color7"].i;
    Colors::BValueColors[7] = data["color8"].i;
    Colors::BValueColors[8] = data["color9"].i;
    Colors::BValueColors[9] = data["color10"].i;

    Colors::BValueRanges[0] = static_cast<int>(data["level1"].f*100.0f);
    Colors::BValueRanges[1] = static_cast<int>(data["level2"].f*100.0f);
    Colors::BValueRanges[2] = static_cast<int>(data["level3"].f*100.0f);
    Colors::BValueRanges[3] = static_cast<int>(data["level4"].f*100.0f);
    Colors::BValueRanges[4] = static_cast<int>(data["level5"].f*100.0f);
    Colors::BValueRanges[5] = static_cast<int>(data["level6"].f*100.0f);
    Colors::BValueRanges[6] = static_cast<int>(data["level7"].f*100.0f);
    Colors::BValueRanges[7] = static_cast<int>(data["level8"].f*100.0f);
    Colors::BValueRanges[8] = static_cast<int>(data["level9"].f*100.0f);
    Colors::BValueRanges[9] = static_cast<int>(data["level10"].f*100.0f);

    if (data["save"].b)
    {
        Application::instance()->Write();
    }
}


/// NOTE: after this point, MIMainWindow, and MainFrame are roughly in sync

void MIMainWindow::OnLoadDictReplace()
{
    QString s = QFileDialog::getOpenFileName(this, "Choose a dictionary file", "",
                                          "Dictionary files (*.pdb,*.cif,*.ent);;All files (*.*)");
    if (!s.isEmpty())
    {
        MIFitDictionary()->LoadDictionary(s.toAscii().constData(), false, true);
    }
}

void MIMainWindow::ShowDictEditor(const char *type)
{

    // note this dialog is dynamcially allocated, and is deleted when closed.
    // this is so that history can function
    std::string s = ::format("Dictionary Editor - %s", type);
    DictEditDialog *dlg = new DictEditDialog(this);
    dlg->setWindowTitle(s.c_str());
    dlg->resize(800, 720);
    dlg->setModal(true);
    dlg->exec();
}

void MIMainWindow::OnLoadDictAppend()
{
    QString s = QFileDialog::getOpenFileName(this, "Choose a dictionary file", "",
                                          "Dictionary files (*.pdb,*.cif,*.ent);;All files (*.*)");
    if (!s.isEmpty())
    {
        MIFitDictionary()->LoadDictionary(s.toAscii().constData(), true);
    }
}


void MIMainWindow::OnEditDictResidue()
{
    std::vector<std::string> resList = MIFitDictionary()->GetDictResList();
    if (resList.empty())
    {
        Logger::message("No residues in dictionary to edit");
        return;
    }

    QStringList choices;
    foreach (const std::string &res, resList)
    {
        choices += res.c_str();
    }

    GenericDataDialog dlg(this);
    dlg.setWindowTitle("Edit Residue");
    dlg.addComboField("Select a residue to edit", choices, 0);
    if (dlg.exec() != QDialog::Accepted)
    {
        return;
    }
    std::string type = resList[dlg.value(0).toInt()];

    if (MIFitDictionary()->DictContains(type.c_str()))
    {
        if (MIFitGeomRefiner()->EditEntry(type.c_str()))
        {
            ShowDictEditor(type.c_str());
        }
    }
    else
    {
        std::string rplc_warn;
        rplc_warn += "Cannot find residue \"";
        rplc_warn += type.c_str();
        rplc_warn += "\" in the dictionary.";
        Logger::message(rplc_warn);
    }
}

void MIMainWindow::OnLoadLigMol()
{
    QString filename = QFileDialog::getOpenFileName(this, "Pick a MOL file", "", "MDL file formats (*.mol, *.sdf, *.sd);;MOL files (*.mol);;SD files (*.sdf, *.sd);;All files (*.*)");
    if (filename.isEmpty())
    {
        return;
    }
    std::string code;
    while (true)
    {
        QString str = QInputDialog::getText(this, "Ligand ID", "ID Code (3 letters)", QLineEdit::Normal, "UNK");
        if (str.isEmpty())
        {
            return;
        }
        code = str.toStdString();
        if (code.size() > 3
            || (code.size() > 0 && !isalnum(code[0]))
            || (code.size() > 1 && !isalnum(code[1]))
            || (code.size() > 2 && !isalnum(code[2])))
        {
            QMessageBox::critical(this, "Invalid ID", "ID Code must be alphanumeric and at most 3 letters long");
        }
        else
        {
            break;
        }
    }
    OnLoadLigand("*.mol", filename.toAscii().constData(), "", code.c_str());
    ShowDictEditor(code.c_str());
}

void MIMainWindow::OnLoadLigPdb()
{
    QString filename = QFileDialog::getOpenFileName(this, "Pick a PDB file", "", "PDB files (*.pdb);;All files (*.*)");
    if (filename.isEmpty())
    {
        return;
    }
    std::string code = OnLoadLigand("*.pdb", filename.toAscii().constData(), "", "");
    if (code.size())
    {
        ShowDictEditor(code.c_str());
    }

}

void MIMainWindow::OnLoadLigSmi()
{
    MIData data;
    data["mode"].radio = 0;
    data["mode"].radio_count = 3;
    data["filename"].str = "";
    data["code"].str = "LIG";
    data["smiles"].str = "";
    data["dbquery"].str = "";

    SmilesDialog dlg(this);
    dlg.setWindowTitle("Import Smiles String");
    if (dlg.exec() != QDialog::Accepted)
    {
        return;
    }
    dlg.GetResults(data);

    std::string res;
    switch (data["mode"].radio)
    {
    case 0:  // From File
        res = OnLoadLigand("*.smi", data["filename"].str.c_str(), "", data["code"].str.c_str());
        break;
    case 1:  // From String
        res = OnLoadLigand("*.smi", "", data["smiles"].str.c_str(), data["code"].str.c_str());
        break;
    case 2:  // From Db
    {
        if (Application::instance()->SmilesDbCommand.size())
        {
            QStringList args;
            QString cmd(Application::instance()->SmilesDbCommand.c_str());
            args.append(data["dbquery"].str.c_str());

            QProcess proc(this);
            proc.start(cmd, args);
            if (proc.waitForReadyRead(10000))
            {
                QByteArray qba = proc.readAll();
                std::string str(qba.constData());
                res = OnLoadLigand("*.smi", "", str, data["code"].str.c_str());
            }
            else
            {
                QMessageBox::critical(this, "Command failed", "Error running smiles database command.\nCheck the setting in File/Preferences.../Environment.");
                proc.kill();
                return;
            }
            proc.kill();

            break;
        }
    }
    default:
        return;
        break;
    }

    if (res==std::string(""))
        return;

    ShowDictEditor(data["code"].str.c_str());
}

void MIMainWindow::OnLoadLigCif()
{
    QString filename = QFileDialog::getOpenFileName(this, "Pick a CIF file", "", "CIF files (*.cif);;All files (*.*)");
    if (filename.isEmpty())
    {
        return;
    }
    std::string code = OnLoadLigand("*.cif", filename.toAscii().constData(), "", "");
    if (code.size())
    {
        ShowDictEditor(code.c_str());
    }
}

std::string MIMainWindow::OnLoadLigand(std::string wildcard, std::string filename, std::string smiles, std::string code)
{
    std::string smi_err;                       //Error message returned from smiles library
    std::string err_report;                    //Report of smiles error to user

    // get MI data
    MIMolInfo mi;

    // NOTE: each reading method (fio or LoadSmiles) allocates at least one
    // new RESIDUE (mi.res), which we pass on LigDictEntry, which takes
    // ownership of it, and deletes it when done.

    if (filename.size() == 0 && smiles.size() > 0)
    {
        //Load SMILES from string into dictionary
        if (!LoadSmiles(smiles, code, mi))
        {
            return "";
        }
    }
    else
    {
        MIMolIO fio;
        int sel = fio.getReaderIndex(wildcard.c_str());
        if (!fio.Read(mi, filename.c_str(), sel))
        {
            return "";
        }
    }


    LigDictEntry entry(mi.res);
    entry.bonds = mi.bonds;
    entry.angles = mi.angles;
    entry.torsions = mi.tordict;
    entry.planes = mi.planedict;
    entry.chirals = mi.chiralsdict;


    if (entry.res->atomCount() == 0)
    {
        QMessageBox::critical(this, "Error: No atoms in file", "No atoms were found in the file. Please check the file format.");
        return "";
    }

    //For MOLfile, set flag for whether it has 2D or 3D coordinates
    if (wildcard==std::string("*.mol"))
    {
        entry.res->setName1((entry.res->type() == "3D") ? '.' : 'M');
    }

    if (code.size() > 0)
    {
        entry.res->setType(code);
        entry.res->setName("1");
    }
    if (MIFitDictionary()->DictContains(entry.res->type()))
    {
        std::string rplc_warn = "The residue code you have entered is already in use.  ";
        rplc_warn += "Replace the existing residue \"";
        rplc_warn += entry.res->type().c_str();
        rplc_warn += "\" in the dictionary?";

        int replace = QMessageBox::question(this, "Replace residue?", rplc_warn.c_str(),
                                            QMessageBox::Yes | QMessageBox::No | QMessageBox::Cancel);
        if (replace == QMessageBox::No)
        {
            QString newcode = QInputDialog::getText(this, "Enter Three-Letter Code", "New Code:");
            if (newcode.isEmpty())
            {
                return "";
            }
            else
            {
                entry.res->setType(newcode.toStdString());
            }
            code = newcode.toStdString();
        }
        if (replace == QMessageBox::Cancel)
        {
            return "";
        }
    }

    if (DupeAtomNames(entry.res) != 0)
    {
        if (QMessageBox::question(this, "Duplicate Atom Names", "File contains duplicate atom names...Generate unique names?", QMessageBox::Yes | QMessageBox::No, QMessageBox::No) == QMessageBox::Yes)
        {
            renameResidueAtomsToUnique(entry.res);
        }
        else
        {
            return "";
        }
    }
    if (DupeAtomNames(entry.res) != 0)
    {
        QMessageBox::critical(this, "Error: Duplicate Atom Names", "File contains duplicate atom names...please edit with unique names!");
        return "";
    }

    LigPostProcessor lpp(entry, wildcard.c_str());
    try
    {
        lpp.Process();
    }
    catch (const char *ex)
    {
        std::string message("Unable to load ligand.\n");
        message += ex;
        QMessageBox::critical(this, "Unable to load ligand.", message.c_str());
        return "";
    }
    entry.res->setPrefBonds(entry.bonds);
    entry.res->setPrefAngles(entry.angles);


    unsigned int level;
    if (!MIFitDictionary()->DictHCheck(entry.res, level))
    {
        int ret = QMessageBox::question(this, "H Check", "There is no dictionary, or the number of hydrogens in the dictionary and the residue you are replacing"
                               "are different.\nDo you want to load the appropriate default dictionary?\n"
                               "\nNote: any changes to the current dictionary will be lost.", QMessageBox::Yes | QMessageBox::No | QMessageBox::Cancel);
        if (ret == QMessageBox::Yes)
        {
            std::string path, name, ext;
            SplitPath(Application::instance()->getDictionary(), &path, &name, &ext);
            path += "/";
            if (level == DictionaryHLevel::NoHydrogens)
            {
                path += "dict.noh.pdb";
            }
            if (level == DictionaryHLevel::Polar)
            {
                path += "dict.polarh.pdb";
            }
            if (level == DictionaryHLevel::All)
            {
                path += "dict.allh.pdb";
            }
            MIFitDictionary()->LoadDictionary(path.c_str(), false, true, level);
        }
        if (ret == QMessageBox::Cancel)
        {
            return "";
        }
    }


    //This code substitutes for the "LoadDictionary" function opt->dict.LoadRes(&entry.res, true, true);
    MIFitDictionary()->LoadRes(entry.res, true, true);
    for (std::vector<TORSDICT>::iterator tor = entry.torsions.begin(); tor != entry.torsions.end(); ++tor)
    {
        MIFitDictionary()->AddTorsion(*tor);
    }
    for (std::vector<PLANEDICT>::iterator pln = entry.planes.begin(); pln != entry.planes.end(); ++pln)
    {
        MIFitDictionary()->AddPlane(*pln);
    }
    for (std::vector<CHIRALDICT>::iterator chrl = entry.chirals.begin(); chrl != entry.chirals.end(); ++chrl)
    {
        MIFitDictionary()->AddChiral(*chrl);
    }

    if (MIFitGeomRefiner()->EditEntry(entry.res->type().c_str()))
    {
        std::string ret = entry.res->type().c_str();
        return ret;
    }
    return "";
}

bool MIMainWindow::LoadSmiles(std::string &smiles, std::string &tlcode, MIMolInfo &mi)
{
    std::string smistr(smiles.c_str());
    SMILES smi;
    if (!smi.Read(smistr, mi))
    {
        return false;
    }
    mi.res->setType(tlcode);
    Residue::fixnames(mi.res);
    return true;
}

void MIMainWindow::LeftFooter(const std::string &msg, int timeout)
{
    logWindow->append(msg.c_str());
    statusBar()->showMessage(msg.c_str(), timeout);
}

void MIMainWindow::MiddleFooter(const std::string &msg)
{
    middleFooter->setText(msg.c_str());
}

void MIMainWindow::RightFooter(const std::string &msg)
{
    rightFooter->setText(msg.c_str());
}

void MIMainWindow::OnBackgroundColor()
{
    Application::instance()->backgroundColor();
    MIGLWidget *view = currentMIGLWidget();
    if (view)
        view->ReDrawAll();
}


void MIMainWindow::OnSaveDict()
{
    Application::instance()->saveDict();
}

BatchJobManager *MIMainWindow::GetJobManager()
{
    return JobManager;
}
bool MIMainWindow::isJobLimit()
{
    return JobManager->numberOfRunningJobs() >= Application::instance()->concurrentJobLimit;
}

void MIMainWindow::UpdateToolBar()
{
    foreach (QAction *action, tool_bar->actions())
    {
        CurrentMIGLWidgetAction *miAction = dynamic_cast<CurrentMIGLWidgetAction*>(action);
        if (miAction)
            miAction->update();
    }
}

void MIMainWindow::OnUpdateStereoToggle()
{
    bool stereo = MIConfig::Instance()->GetProfileInt("View Parameters", "stereo", 0) != 0;
    _stereoToggleAction->setChecked(stereo);
}

void MIMainWindow::OnUpdateHardwareStereo()
{
    bool hardwareStereo = MIConfig::Instance()->GetProfileInt("View Parameters", "hardwareStereo", 0) != 0;
    _hardwareStereoAction->setEnabled(Application::instance()->isHardwareStereoAvailable());
    _hardwareStereoAction->setChecked(hardwareStereo);
}

void MIMainWindow::OnStereoToggle()
{
    Application::instance()->toggleStereo();
}


void MIMainWindow::OnHardwareStereo()
{
    Application::instance()->toggleHardwareStereo();
}

void MIMainWindow::OnScript()
{
}

void MIMainWindow::OnSideChainTool()
{
    QPoint pos = QCursor::pos();
    side_menu->exec(pos);
}

void MIMainWindow::OnHideTool()
{
    QPoint pos = QCursor::pos();
    hide_menu->exec(pos);
}

void MIMainWindow::OnColorTool()
{
    QPoint pos = QCursor::pos();
    color_menu->exec(pos);
}

void MIMainWindow::OnShowTool()
{
    QPoint pos = QCursor::pos();
    show_menu->exec(pos);
}

MIMainWindow*MIMainWindow::instance()
{
    if (_instance)
        return _instance;
    new MIMainWindow(); // sets _instance
    return _instance;
}

QDockWidget*MIMainWindow::AddAsDockWidget(QWidget *w, const std::string &name, Qt::DockWidgetArea area)
{
    QDockWidget *dock = new QDockWidget(tr(name.c_str()), this);
    dock->setObjectName(name.c_str());
    w->setParent(dock);
    dock->setWidget(w);
    addDockWidget(area, dock);
    viewMenu->addAction(dock->toggleViewAction());
    return dock;
}

void MIMainWindow::createDockWindows()
{
    logWindow = new QTextEdit(logDock);
    logWindow->setReadOnly(true);
    //  logWindow->document()->setMaximumBlockCount(1000);
    logDock = AddAsDockWidget(logWindow, "Log", Qt::BottomDockWidgetArea);

    // add rama plot
    QWidget *graphwin = RamaPlotMgr::instance()->getGraphWin();
    ramaDock = AddAsDockWidget(graphwin, "Ramachandran plot", Qt::BottomDockWidgetArea);

    //add Models tree
    modelsView = new ModelsView(modelsDock);
    modelsDock = AddAsDockWidget(modelsView, "Models", Qt::LeftDockWidgetArea);

    //add Display tree
    QWidget *displayView = new DisplayView(displayDock);
    displayDock = AddAsDockWidget(displayView, "Display", Qt::LeftDockWidgetArea);

    // add Jobs tree
    JobsView *jobsView = new JobsView(jobsDock);
    jobsView->update(JobManager);
    jobsDock = AddAsDockWidget(jobsView, "Jobs", Qt::LeftDockWidgetArea);
}

static bool is_valid_extension(const char *ext)
{
    return (strcmp(ext, ".pdb") == 0
            || strcmp(ext, ".mlw") == 0
            || strcmp(ext, ".mtz") == 0
            || strcmp(ext, ".phs") == 0
            || strcmp(ext, ".fcf") == 0
            || // we no longer support these extensions b/c their readers are broken
               // || (strcmp(ext, ".sca") == 0)
               // || (strcmp(ext, ".ref") == 0)
            strcmp(ext, ".cif") == 0
            || strcmp(ext, ".map") == 0
            || strcmp(ext, ".py") == 0);
}


void MIMainWindow::OpenFiles(const std::vector<std::string> &files, bool newWin)
{
    if (newWin)
    {
        newWindow();
    }

    for (size_t i = 0; i<files.size(); ++i)
    {
        const std::string &file = files[i];
        const char *arg = file.c_str();

        if (arg[0] == '-')
        {
            // command line switch
            if (strcmp(arg, "-altTextRender") == 0)
            {
                MIConfig::Instance()->WriteProfileInt("View Parameters", "AlternateTextRendering", 1);
                OnExit();
            }
            continue;
        }

        // convert to absolute path
        QFileInfo fi(files[i].c_str());
        QString s = fi.absoluteFilePath();
        std::string ss = s.toStdString();
        arg = ss.c_str();
        const char *ext = file_extension(arg);

        if (!is_valid_extension(ext))
        {
            Log(::format("Unknown file type %s, ingoring.", arg));
            continue;
        }

        if (strcmp(ext, ".py") == 0)
        {
            runPythonScript(file);
            continue;
        }


        if (!currentMIGLWidget())
            newWindow();
        if (!currentMIGLWidget())
        {
            Log(::format("Error creating active MIGLWidget!"));
            continue;
        }

        setCurrentFile(arg);
        currentMIGLWidget()->OpenAnyFile(arg);
    }
}


static bool GetOpenFilenames(std::vector<std::string> &fnames)
{
    initializeFileDialog();
    // TODO use QFileDialog static functions to get native OS file dialog
    QString path = Application::instance()->latestFileBrowseDirectory("");
    fileDialog->setDirectory(path);
    if (fileDialog->exec())
    {
        Application::instance()->latestFileBrowseDirectory(fileDialog->directory().absolutePath());
        foreach (QString file, fileDialog->selectedFiles())
        {
            fnames.push_back(file.toStdString());
        }
        return true;
    }
    return false;
}

static void SplitPath(const std::string &origPath,
                      std::string *dirname,
                      std::string *fname,
                      std::string *ext)
{
    *dirname = std::string("");
    *fname = std::string("");
    *ext = std::string("");

    QFileInfo qfile(origPath.c_str());
    *dirname = qfile.path().toStdString();
    *fname = qfile.fileName().toStdString();
    *ext = qfile.suffix().toStdString();
}

void MIMainWindow::AfterInit()
{
    QStringList args = QApplication::arguments();
    std::vector<std::string> arglist;

#ifdef DEBUG
    if (!(args.count() > 1 && args.at(1) == QString("-recordoff")))
    {
        arglist.push_back("-record");
    }
#endif

    args.pop_front();
    foreach (QString str, args)
    {
        arglist.push_back(str.toStdString());
    }
    OpenFiles(arglist);

}


void MIMainWindow::OnMenuValidate()
{
//    menu_bar->validateActions();
    //menu_bar->validateUpdates();
}

#include "../images/cursors/MW_CROSS.xpm"
#include "../images/cursors/MW_TRANS.xpm"
#include "../images/cursors/MW_ROTAT.xpm"
#include "../images/cursors/MW_ZROTA.xpm"
#include "../images/cursors/ZCURSOR.xpm"
#include "../images/cursors/MW_TORSI.xpm"
#include "../images/cursors/MW_CENTE.xpm"
#include "../images/cursors/MW_SCALE.xpm"
#include "../images/cursors/MW_SLAB.xpm"
#include "../images/cursors/slabDrag.xpm"
#include "../images/cursors/mw_wait1.xpm"
#include "../images/cursors/mw_wait2.xpm"
#include "../images/cursors/mw_wait3.xpm"
#include "../images/cursors/mw_wait4.xpm"
#include "../images/cursors/mw_wait5.xpm"
#include "../images/cursors/mw_wait6.xpm"

void MIMainWindow::SetCursor(int id, QWidget *w)
{
    if (!cursors[0])
    {
        cursors[imhCross] = new QCursor(QPixmap(MW_CROSS));
        cursors[imhTranslate] = new QCursor(QPixmap(MW_TRANS));
        cursors[imhRotate] = new QCursor(QPixmap(MW_ROTAT));
        cursors[imhZRotate] = new QCursor(QPixmap(MW_ZROTA));
        cursors[imhZCursor] = new QCursor(QPixmap(ZCURSOR));
        cursors[imhTorsion] = new QCursor(QPixmap(MW_TORSI));
        cursors[imhCenter] = new QCursor(QPixmap(MW_CENTE));
        cursors[imhScale] = new QCursor(QPixmap(MW_SCALE));
        cursors[imhSlab] = new QCursor(QPixmap(MW_SLAB));
        cursors[imhWait1] = new QCursor(QPixmap(mw_wait1));
        cursors[imhWait2] = new QCursor(QPixmap(mw_wait2));
        cursors[imhWait3] = new QCursor(QPixmap(mw_wait3));
        cursors[imhWait4] = new QCursor(QPixmap(mw_wait4));
        cursors[imhWait5] = new QCursor(QPixmap(mw_wait5));
        cursors[imhWait6] = new QCursor(QPixmap(mw_wait6));
        cursors[imhSlabDrag] = new QCursor(QPixmap(slabDrag));
    }

    if (!w)
        w = this;

    if (id >=0 && id <= (int)imhCount)
        w->setCursor(*cursors[id]);
    else
        w->unsetCursor();
}

void MIMainWindow::updateNavigator()
{
    navigator->updateGL();
}

void MIMainWindow::Debug(const std::string &msg)
{
#ifdef DEBUG
    Log(msg);
#endif
}

void MIMainWindow::closeEvent(QCloseEvent *event)
{
    saveLayout();
    mdiArea->closeAllSubWindows();
    if (currentMIGLWidget())
    {
        event->ignore();
    }
    else
    {
        event->accept();
    }
}

void MIMainWindow::newWindow()
{
    createMIGLWidget()->show();
}



void MIMainWindow::updateMenus()
{
    updateFileMenu();

    bool hasMIGLWidget = (currentMIGLWidget() != 0);
    closeAllAct->setEnabled(hasMIGLWidget);
    tileAct->setEnabled(hasMIGLWidget);
    cascadeAct->setEnabled(hasMIGLWidget);
    nextAct->setEnabled(hasMIGLWidget);
    previousAct->setEnabled(hasMIGLWidget);

    view_menu->setEnabled(hasMIGLWidget);
    show_menu->setEnabled(hasMIGLWidget);
    render_menu->setEnabled(hasMIGLWidget);
    model_menu->setEnabled(hasMIGLWidget);
    fit_menu->setEnabled(hasMIGLWidget);
    refi_menu->setEnabled(hasMIGLWidget);
    analyze_menu->setEnabled(hasMIGLWidget);
}

void MIMainWindow::updateWindowMenu()
{
    windowMenu->clear();
    windowMenu->addAction(closeAct);
    windowMenu->addAction(closeAllAct);
    windowMenu->addSeparator();
    windowMenu->addAction(tileAct);
    windowMenu->addAction(cascadeAct);
    windowMenu->addSeparator();
    windowMenu->addAction(nextAct);
    windowMenu->addAction(previousAct);

    QList<QMdiSubWindow*> windows = mdiArea->subWindowList();

    for (int i = 0; i < windows.size(); ++i)
    {
        MIGLWidget *child = dynamic_cast<MIGLWidget*>(windows.at(i)->widget());
        if (!child)
            continue;

        QString text;
        if (i == 0)
            windowMenu->addSeparator();
        if (i < 9)
        {
            text = tr("&%1 %2").arg(i + 1)
                   .arg(child->GetTitle().c_str());
        }
        else
        {
            text = tr("%1 %2").arg(i + 1)
                   .arg(child->GetTitle().c_str());
        }
        QAction *action  = windowMenu->addAction(text);
        action->setCheckable(true);
        action->setChecked(child == currentMIGLWidget());
        connect(action, SIGNAL(triggered()), windowMapper, SLOT(map()));
        windowMapper->setMapping(action, windows.at(i));
    }
}

QWidget*MIMainWindow::createMIGLWidget()
{
    QWidget *w = new MIGLWidget;
    w->setObjectName("MIGLWidget");
    mdiArea->addSubWindow(w);
    w->showMaximized();
    foreach (QKeySequence k, closeAction->shortcuts())
        grabShortcut(k, Qt::ApplicationShortcut);
    return w;
}

void MIMainWindow::createActions()
{

    closeAct = new QAction(tr("Cl&ose\tCtrl+W"), this);
    closeAct->setStatusTip(tr("Close the active window"));
    connect(closeAct, SIGNAL(triggered()),
            mdiArea, SLOT(closeActiveSubWindow()));

    closeAllAct = new QAction(tr("Close &All"), this);
    closeAllAct->setStatusTip(tr("Close all the windows"));
    connect(closeAllAct, SIGNAL(triggered()),
            mdiArea, SLOT(closeAllSubWindows()));
    ///////////////////////////////////////////////////////////

    for (int i = 0; i < MaxRecentFiles; ++i)
    {
        recentFileActs[i] = new QAction(this);
        recentFileActs[i]->setVisible(false);
        connect(recentFileActs[i], SIGNAL(triggered()),
                this, SLOT(openRecentFile()));
    }

    tileAct = new QAction(tr("&Tile"), this);
    tileAct->setStatusTip(tr("Tile the windows"));
    connect(tileAct, SIGNAL(triggered()), mdiArea, SLOT(tileSubWindows()));

    cascadeAct = new QAction(tr("&Cascade"), this);
    cascadeAct->setStatusTip(tr("Cascade the windows"));
    connect(cascadeAct, SIGNAL(triggered()), mdiArea, SLOT(cascadeSubWindows()));

    nextAct = new QAction(tr("Ne&xt"), this);
    nextAct->setStatusTip(tr("Move the focus to the next window"));
    connect(nextAct, SIGNAL(triggered()),
            mdiArea, SLOT(activateNextSubWindow()));

    previousAct = new QAction(tr("Pre&vious"), this);
    previousAct->setStatusTip(tr("Move the focus to the previous "
                                 "window"));
    connect(previousAct, SIGNAL(triggered()),
            mdiArea, SLOT(activatePreviousSubWindow()));

}

void MIMainWindow::updateFileMenu()
{
    bool hasDocument = currentMIGLWidget() != 0;
    fileSaveAction->setEnabled(hasDocument);
    fileSaveAsAction->setEnabled(hasDocument);
    bool canExport = hasDocument
                     && currentMIGLWidget()->GetDisplaylist()->CurrentItem() != NULL
                     && !MIBusyManager::instance()->Busy();
    exportModelAction->setEnabled(canExport);
    filePrintAction->setEnabled(hasDocument);
    copyCanvasAction->setEnabled(hasDocument);
    exportImageAction->setEnabled(hasDocument);

    closeAction->setEnabled(hasDocument);
}

void MIMainWindow::OnFileSave()
{
    if (currentMIGLWidget())
        currentMIGLWidget()->OnFileSave();
}

void MIMainWindow::OnFileSaveAs()
{
    if (currentMIGLWidget())
        currentMIGLWidget()->OnFileSaveAs();
}

void MIMainWindow::OnExportModel()
{
    if (currentMIGLWidget())
        currentMIGLWidget()->OnExportModel();
}

void MIMainWindow::OnPrint()
{
    if (currentMIGLWidget())
        currentMIGLWidget()->OnPrint();
}

void MIMainWindow::OnEditCopy()
{
    if (currentMIGLWidget())
        currentMIGLWidget()->OnEditCopy();
}

void MIMainWindow::OnExportImage()
{
    if (currentMIGLWidget())
        currentMIGLWidget()->OnExportImage();
}


void MIMainWindow::fill_surf_menu(QMenu *surf_menu)
{
    new CurrentMIGLWidgetAction("Surface Residue", "Calculate van der Waal dot surface around picked residue", surf_menu, SLOT(OnObjectSurfaceresidue()), SLOT(OnUpdateObjectSurfaceresidue(QAction*)));

    new CurrentMIGLWidgetAction("Surface Residues", "Surface all the residues on the stack", surf_menu, SLOT(OnObjectSurfaceresidues()), SLOT(OnUpdateObjectSurfaceresidues(QAction*)));

    new CurrentMIGLWidgetAction("Surface Atom", "Surface the last atom picked", surf_menu, SLOT(OnObjectSurfaceatom()), SLOT(OnUpdateObjectSurfaceAtom(QAction*)));

    new CurrentMIGLWidgetAction("Surface Atoms", "Surface all the atoms on the stack", surf_menu, SLOT(OnObjectSurfaceAtoms()), SLOT(OnUpdateObjectSurfaceatoms(QAction*)));

    new CurrentMIGLWidgetAction("van der Waal Surface", "Calculate van der Waal dot surface", surf_menu, SLOT(OnVdwDotSurface()));

    new CurrentMIGLWidgetAction("Solvent Exposed Surface", "Surface through the center of the solvent atoms touching the molecule", surf_menu, SLOT(OnSurfaceSolvent()), SLOT(OnUpdateSurfaceSolvent(QAction*)));

    new CurrentMIGLWidgetAction("Sphere around atom", "Calculate sphere around last picked atom", surf_menu, SLOT(OnObjectSurfaceSpherearoundatom()), SLOT(OnUpdateObjectSurfaceSpherearoundatom(QAction*)));

    surf_menu->addSeparator();
    new CurrentMIGLWidgetAction("Clear Surface", "Clears the surface dots", surf_menu, SLOT(OnObjectSurfaceClearsurface()));

}


void MIMainWindow::createMenus()
{
    CurrentMIGLWidgetAction *miGLWidgetAction;

    //// Make a menubar
    QMenu *file_menu = menuBar()->addMenu("&File");

    QList<QAction *> toolBarActions;
    QAction *action;

    fileNewAction = file_menu->addAction(tr("New"), this, SLOT(OnNew()));
    fileNewAction->setToolTip(tr("Create a new document"));
    fileNewAction->setShortcut(tr("Ctrl+N"));

    fileOpenAction = file_menu->addAction(tr("Open models, data, maps, etc..."), this, SLOT(OnFileOpen()));
    fileOpenAction->setToolTip(tr("Open file(s) into the current document"));
    fileOpenAction->setShortcut(tr("Ctrl+O"));

    action = file_menu->addAction(tr("Open into new document..."), this, SLOT(OnFileOpenNew()));
    action->setToolTip(tr("Open file(s) into a new document"));
    action->setShortcut(tr("Ctrl+Shift+O"));

    closeAction = file_menu->addAction(tr("Close\tCtrl+W"), this, SLOT(OnClose()));

    fileSaveAction = file_menu->addAction(tr("Save Session"), this, SLOT(OnFileSave()));
    fileSaveAction->setToolTip(tr("Save the entire document to a file"));
    fileSaveAction->setShortcut(tr("Ctrl+S"));

    fileSaveAsAction = file_menu->addAction(tr("Save Session As..."), this, SLOT(OnFileSaveAs()));
    fileSaveAsAction->setToolTip(tr("Save the document under a new name"));
    fileSaveAsAction->setShortcut(tr("Ctrl+Shift+S"));

    exportModelAction = file_menu->addAction(tr("Export PDB..."), this, SLOT(OnExportModel()));
    exportModelAction->setToolTip(tr("Saves the current model as a PDB format file"));
    exportModelAction->setShortcut(tr("Ctrl+Alt+S"));

    filePrintAction = file_menu->addAction(tr("Print..."), this, SLOT(OnPrint()));
    filePrintAction->setToolTip(tr("Print the canvas"));

    file_menu->addSeparator();

    copyCanvasAction = file_menu->addAction(tr("Copy Canvas"), this, SLOT(OnEditCopy()));
    copyCanvasAction->setToolTip(tr("Copy the canvas to the clipboard to paste into other programs"));
    copyCanvasAction->setShortcut(tr("Ctrl+C"));

    exportImageAction = file_menu->addAction(tr("Export Image As..."), this, SLOT(OnExportImage()));
    exportImageAction->setToolTip(tr("Export view to a graphics file"));

    file_menu->addSeparator();

    action = file_menu->addAction(tr("Manage Crystals..."), this, SLOT(OnManageCrystals()));
    action->setToolTip(tr("Edit the crystal database"));

    file_menu->addSeparator();

    action = file_menu->addAction(tr("Preferences..."), this, SLOT(OnPreferences()));
    action->setToolTip(tr("Edit preferences"));

    action = file_menu->addAction(tr("Define Atom Colors..."), this, SLOT(OnDefineAtomColors()));
    action->setToolTip(tr("Define Atom Colors"));

    action = file_menu->addAction(tr("B-Value Colors Ranges..."), this, SLOT(OnDefineBValueColors()));
    action->setToolTip(tr("Define B-Value Coloring Ranges"));

    MIMainWindow::instance()->addRecentFileActions(file_menu);

    action = file_menu->addAction(tr("E&xit"), this, SLOT(OnExit()));
    action->setShortcut(tr("Alt+X"));

    connect(file_menu, SIGNAL(aboutToShow()), SLOT(updateFileMenu()));

    updateRecentFileActions();

    //TODO: use disabled icon specialization for QIcons?
    fileNewAction->setIcon(QIcon(new_xpm));
    toolBarActions += fileNewAction;
    fileOpenAction->setIcon(QIcon(open_xpm));
    toolBarActions += fileOpenAction;
    fileSaveAction->setIcon(QIcon(save_xpm));
    toolBarActions += fileSaveAction;
    filePrintAction->setIcon(QIcon(print_xpm));
    toolBarActions += filePrintAction;
    copyCanvasAction->setIcon(QIcon(copy_xpm));
    toolBarActions += copyCanvasAction;
    toolBarActions += 0;

    QToolBar *toolBar = addToolBar("Display Tools");
    toolBar->setObjectName("Display Tools");
    toolBar->addAction(QIcon(colortool_xpm), tr("Color"), this, SLOT(OnColorTool()));
    toolBar->addAction(QIcon(showtool_xpm), tr("Show"), this, SLOT(OnShowTool()));
    toolBar->addAction(QIcon(sidetool_xpm), tr("Side Chain"), this, SLOT(OnSideChainTool()));
    toolBar->addAction(QIcon(hidetool_xpm), tr("Hide"), this, SLOT(OnHideTool()));

    QMenu *di_menu = menuBar()->addMenu(tr("&Dictionary"));
    view_menu = menuBar()->addMenu(tr("&Viewpoint"));
    show_menu = menuBar()->addMenu(tr("&Show"));
    render_menu = menuBar()->addMenu(tr("Re&nder"));
    model_menu = menuBar()->addMenu(tr("&Model"));
    fit_menu = menuBar()->addMenu(tr("F&it"));
    refi_menu = menuBar()->addMenu(tr("&Refine"));
    analyze_menu = menuBar()->addMenu(tr("&Analyze"));
    _jobMenu = menuBar()->addMenu(tr("&Job"));

    connect(show_menu, SIGNAL(aboutToShow()), this, SLOT(updateShowMenu()));

    canvas_menu = show_menu->addMenu("Canvas");

    QMenu *stack_menu = canvas_menu->addMenu("S&tack");
    miGLWidgetAction = new CurrentMIGLWidgetAction("Show/Hide Atom Stack", "Toggle the atom stack visibility", stack_menu, SLOT(OnViewAtomstack()), SLOT(OnUpdateViewAtomstack(QAction*)));
    miGLWidgetAction->setCheckable(true);

    stack_menu->addSeparator();
    new CurrentMIGLWidgetAction("Clear &top item", "Clear the top item on the stack", stack_menu, SLOT(OnObjectStackDeletetopitem()));
    new CurrentMIGLWidgetAction("&Clear stack", "Clear all items on the stack", stack_menu, SLOT(OnObjectClearstack()));
    stack_menu->addSeparator();
    new CurrentMIGLWidgetAction("&Expand Top Residue", "Expand the top residue to put all of its atoms on the stack", stack_menu, SLOT(OnObjectStackExpandtopallatomsinresidue()), SLOT(OnUpdateObjectStackExpandtopallatomsinresidue(QAction*)));
    new CurrentMIGLWidgetAction("Expand &range to All Residue", "Expand stack to include all residues between top two", stack_menu, SLOT(OnObjectStackExpandtop2residues()), SLOT(OnUpdateObjectStackExpandtop2residues(QAction*)));
    new CurrentMIGLWidgetAction("Expand Range to All &Atoms", "Expand stack to include all atoms in all residues between top two", stack_menu, SLOT(OnObjectStackExpandtop2allatomsinrange()), SLOT(OnUpdateObjectStackExpandtop2allatomsinrange(QAction*)));

    color_menu = new QMenu(this);
    new CurrentMIGLWidgetAction("Color all atoms", "Colors the whole model with current color", color_menu, SLOT(OnShowColorallatoms()), SLOT(OnUpdateShowColorallatoms(QAction*)));
    new CurrentMIGLWidgetAction("Color last picked atom", "Colors the last picked atom with current color", color_menu, SLOT(OnObjectAtomColor()), SLOT(OnUpdateObjectAtomColor(QAction*)));

    new CurrentMIGLWidgetAction("Color all picked atoms", "Colors all the atoms on the stack with current color", color_menu, SLOT(OnObjectAtomsColor()));

    color_menu->addSeparator();
    new CurrentMIGLWidgetAction("Color last picked residue", "Colors the last picked residue with current color", color_menu, SLOT(OnObjectResidueColor()), SLOT(OnUpdateObjectResidueColor(QAction*)));
    new CurrentMIGLWidgetAction("Color all picked Residues", "Colors all the residues on the stack with current color", color_menu, SLOT(OnObjectResiduesColor()), SLOT(OnUpdateObjectResiduesColor(QAction*)));

    new CurrentMIGLWidgetAction("Color residue range", "Colors the residue at the top of stack with current color", color_menu, SLOT(OnObjectResiduerangeColor()), SLOT(OnUpdateObjectResiduerangeColor(QAction*)));

    color_menu->addSeparator();
    new CurrentMIGLWidgetAction("&Undo color", "Undo the last color or radius command", color_menu, SLOT(OnShowUndocolorradius()), SLOT(OnUpdateShowUndocolorradius(QAction*)));

    connect(render_menu, SIGNAL(aboutToShow()), SLOT(updateRenderMenu()));

    QActionGroup *renderStyleActionGroup = new QActionGroup(render_menu);

    action = new CurrentMIGLWidgetAction("&Sticks", "Draw bonds as sticks", render_menu, SLOT(OnRenderSticks()), SLOT(OnUpdateRenderSticks(QAction*)));
    action->setActionGroup(renderStyleActionGroup);
    action->setCheckable(true);

    action = new CurrentMIGLWidgetAction("&Knob and Stick", "Draw atoms as knobs and bonds as stick", render_menu, SLOT(OnRenderingBallandstick()), SLOT(OnUpdateRenderingBallandstick(QAction*)));
    action->setActionGroup(renderStyleActionGroup);
    action->setCheckable(true);

    action = new CurrentMIGLWidgetAction("&Ball and Cylinder", "Draw bonds as ball and cylinders for bonds", render_menu, SLOT(OnRenderBallandcylinder()), SLOT(OnUpdateRenderBallandcylinder(QAction*)));
    action->setActionGroup(renderStyleActionGroup);
    action->setCheckable(true);

    action = new CurrentMIGLWidgetAction("&CPK", "Draw atoms as CPK/space filling", render_menu, SLOT(OnRenderSpacefilling()), SLOT(OnUpdateRenderSpacefilling(QAction*)));
    action->setActionGroup(renderStyleActionGroup);
    action->setCheckable(true);

    render_menu->addSeparator();
    color_menu->setTitle("Color");
    render_menu->addMenu(color_menu);
    render_menu->addAction("Set Background C&olor...", this, SLOT(OnBackgroundColor()));
    render_menu->addSeparator();
    new CurrentMIGLWidgetAction("Set Ball/Cylinder S&ize...", "Set the relative diameters of ball and cylinder", render_menu, SLOT(OnRenderBallsize()), SLOT(OnUpdateRenderBallsize(QAction*)));

    new CurrentMIGLWidgetAction("Set View Center Target S&ize...", "Set size of the marker in the view center", render_menu, SLOT(OnRenderTargetSize()));

    render_menu->addSeparator();
    action = new CurrentMIGLWidgetAction("Depthcue Colo&rs", "Depthcue using darker colors", render_menu, SLOT(OnRenderingDepthcuedcolors()), SLOT(OnUpdateRenderingDepthcuedcolors(QAction*)));
    action->setCheckable(true);

    action = new CurrentMIGLWidgetAction("Depthcue &Lines", "Depthcue with line thickness", render_menu, SLOT(OnRenderingDepthcuedlinewidth()), SLOT(OnUpdateRenderingDepthcuedlinewidth(QAction*)));
    action->setCheckable(true);

    action = new CurrentMIGLWidgetAction("S&mooth Lines", "Smooths the lines", render_menu, SLOT(OnRenderLinesmooth()), SLOT(OnUpdateRenderLinesmooth(QAction*)));
    action->setCheckable(true);

    renderLineThicknessMenu_ = render_menu->addMenu("Line &Thickness");
    QSignalMapper *lineMenuSignalMapper = new QSignalMapper(renderLineThicknessMenu_);
    connect(lineMenuSignalMapper, SIGNAL(mapped(int)), SLOT(setRenderLineThickness(int)));
    QAction *lineThickness1 = renderLineThicknessMenu_->addAction("&1 pixel", lineMenuSignalMapper, SLOT(map()));
    lineMenuSignalMapper->setMapping(lineThickness1, 1);
    action = renderLineThicknessMenu_->addAction("&2 pixel", lineMenuSignalMapper, SLOT(map()));
    lineMenuSignalMapper->setMapping(action, 2);
    action = renderLineThicknessMenu_->addAction("&3 pixel", lineMenuSignalMapper, SLOT(map()));
    lineMenuSignalMapper->setMapping(action, 3);
    action = renderLineThicknessMenu_->addAction("&4 pixel", lineMenuSignalMapper, SLOT(map()));
    lineMenuSignalMapper->setMapping(action, 4);
    QActionGroup *lineThicknessGroup = new QActionGroup(renderLineThicknessMenu_);
    foreach (QAction *a, renderLineThicknessMenu_->actions())
    {
        a->setCheckable(true);
        lineThicknessGroup->addAction(a);
    }
    lineThickness1->setChecked(true);

    action = new CurrentMIGLWidgetAction("Dim Non-active Models", "Dim non-active models", render_menu, SLOT(OnDimNonactiveModels()), SLOT(OnUpdateDimNonactiveModels(QAction*)));
    action->setCheckable(true);

    new CurrentMIGLWidgetAction("Set Amount to Dim Non-active Models", "Dim non-active models", render_menu, SLOT(OnAmountToDimNonactiveModels()));


    QMenu *showres_menu = show_menu->addMenu("Residues");

    new CurrentMIGLWidgetAction("Show &all residues",
                                "Show the whole model",
                                showres_menu,
                                SLOT(OnObjectsAllatoms()));
    new CurrentMIGLWidgetAction("Show &last picked residue",
                                "Show the last picked residue",
                                showres_menu,
                                SLOT(OnObjectShowresidue()),
                                SLOT(OnUpdateObjectShowresidue(QAction*)));
    new CurrentMIGLWidgetAction("Show &all picked residues",
                                "Show all the residues on the stack",
                                showres_menu,
                                SLOT(OnObjectShowresidues()),
                                SLOT(OnUpdateObjectShowresidues(QAction*)));
    new CurrentMIGLWidgetAction("Show &residue range",
                                "Show the residue at the top of stack",
                                showres_menu,
                                SLOT(OnObjectShowresiduerange()),
                                SLOT(OnUpdateObjectShowresiduerange(QAction*)));
    new CurrentMIGLWidgetAction("Show residues within &sphere...",
                                "Show residue if it is witin a radius",
                                showres_menu,
                                SLOT(OnShowShowwithinsphere()));
    showres_menu->addSeparator();
    new CurrentMIGLWidgetAction("&Hide all residues",
                                "Hide the whole model",
                                showres_menu,
                                SLOT(OnHideModel()),
                                SLOT(OnUpdateHideModel(QAction*)));
    new CurrentMIGLWidgetAction("Hide last picked residue",
                                "Hide the last picked residue",
                                showres_menu,
                                SLOT(OnObjectResidueTurnoff()),
                                SLOT(OnUpdateObjectResidueTurnoff(QAction*)));
    new CurrentMIGLWidgetAction("Hide all picked residues",
                                "Hide all the residues on the stack",
                                showres_menu,
                                SLOT(OnObjectResiduesTurnoff()),
                                SLOT(OnUpdateObjectResiduesTurnoff(QAction*)));
    new CurrentMIGLWidgetAction("Hide residue range",
                                "Hide the residue at the top of stack",
                                showres_menu,
                                SLOT(OnObjectResiduerangeTurnoff()),
                                SLOT(OnUpdateObjectResiduerangeTurnoff(QAction*)));

    side_menu = show_menu->addMenu("Sidechains");
    new CurrentMIGLWidgetAction("Show sidechain &atoms", "Show sidechain atoms", side_menu, SLOT(OnShowSidechainAtoms()));

    new CurrentMIGLWidgetAction("&Hide sidechain atoms", "Hide sidechain atoms", side_menu, SLOT(OnHideSidechainAtoms()));

    new CurrentMIGLWidgetAction("&Show sidechain of last picked", "Show the sidechain of the last picked residue", side_menu, SLOT(OnObjectShowsidechain()), SLOT(OnUpdateObjectShowsidechain(QAction*)));
    new CurrentMIGLWidgetAction("S&how sidechains of all picked", "Show the sidechains of all the residues on the stack", side_menu, SLOT(OnObjectShowsidechains()), SLOT(OnUpdateObjectShowsidechains(QAction*)));
    new CurrentMIGLWidgetAction("Sh&ow sidechains of range", "Show the sidechains of the residue range at the top of stack", side_menu, SLOT(OnObjectShowsidechainrange()), SLOT(OnUpdateObjectShowsidechainrange(QAction*)));

    QMenu *backbone_menu = show_menu->addMenu("Backbone");
    new CurrentMIGLWidgetAction("Show backbone as &atoms", "Show backbone as atoms", backbone_menu, SLOT(OnShowBackboneAsAtoms()));

    new CurrentMIGLWidgetAction("Show backbone as &CA trace", "Show backbone as CA trace", backbone_menu, SLOT(OnShowBackboneAsCATrace()));

    new CurrentMIGLWidgetAction("&Hide backbone", "Hide backbone", backbone_menu, SLOT(OnShowHideBackbone()));

    QMenu *symmatoms_menu = show_menu->addMenu("Symmetry Atoms");
    new CurrentMIGLWidgetAction("Show symmetry atoms as &atoms", "", symmatoms_menu, SLOT(OnShowSymmAtomsAsAtoms()));

    new CurrentMIGLWidgetAction("Show symmetry atoms as &CA trace", "", symmatoms_menu, SLOT(OnShowSymmAtomsAsCATrace()));

    new CurrentMIGLWidgetAction("&Hide symmetry atoms", "", symmatoms_menu, SLOT(OnShowHideSymmAtoms()));

    new CurrentMIGLWidgetAction("&Save symmetry atoms", "", symmatoms_menu, SLOT(OnShowSaveSymmAtoms()), SLOT(OnUpdateShowSaveSymmAtoms(QAction*)));


    miGLWidgetAction = new CurrentMIGLWidgetAction("&Unit cell", "Toggle the unit cell visibility", canvas_menu, SLOT(OnViewUnitCell()), SLOT(OnUpdateViewUnitCell(QAction*)));
    miGLWidgetAction->setCheckable(true);

    miGLWidgetAction = new CurrentMIGLWidgetAction("&Contacts", "Toggle contacts on/off", canvas_menu, SLOT(OnViewContacts()), SLOT(OnUpdateViewContacts(QAction*)));
    miGLWidgetAction->setCheckable(true);

    miGLWidgetAction = new CurrentMIGLWidgetAction("&Gnomon", "Toggle the axes visibility", canvas_menu, SLOT(OnViewGnomon()), SLOT(OnUpdateViewGnomon(QAction*)));
    miGLWidgetAction->setCheckable(true);
    new CurrentMIGLWidgetAction("Clea&r Message", "Use this command to clear the message at the screen bottom", canvas_menu, SLOT(OnViewClearmessage()));

    QMenu *labels_menu = show_menu->addMenu("Labels");
    miGLWidgetAction = new CurrentMIGLWidgetAction("&Show/Hide Labels", "Toggle labels on and off", labels_menu, SLOT(OnViewLabels()), SLOT(OnUpdateViewLabels(QAction*)));
    miGLWidgetAction->setCheckable(true);

    new CurrentMIGLWidgetAction("Cl&ear Labels",
                                "Clear the labels permanently - use view/labels to hide them temporarily",
                                labels_menu, SLOT(OnEditClearlabels()));

    new CurrentMIGLWidgetAction("La&bel Options...",
                                "Edit labels and set label options",
                                labels_menu, SLOT(OnEditLabels()));

    new CurrentMIGLWidgetAction("Label E&very nth...", "Label every nth residue", labels_menu, SLOT(OnLabelEveryNth()), SLOT(OnUpdateLabelEveryNth(QAction*)));

    QMenu *secstruct_menu = show_menu->addMenu("Secondary Structure");
    new CurrentMIGLWidgetAction("&Make Ribbon", "Make a backbone ribbon", secstruct_menu, SLOT(OnObjectBackboneribbon()));

    new CurrentMIGLWidgetAction("C&lear Ribbon", "Clear the backbone ribbon", secstruct_menu, SLOT(OnObjectClearribbon()));

    new CurrentMIGLWidgetAction("Ri&bbon Colors...", "Set ribbon colours", secstruct_menu, SLOT(OnSetRibbonColors()));

    secstruct_menu->addSeparator();
    new CurrentMIGLWidgetAction("Show Tube Secondary Structure", "Display a solid tube for the secondary structure", secstruct_menu, SLOT(OnRibbonSecondaryStructure()));

    new CurrentMIGLWidgetAction("Show Schematic Secondary Structure", "Display a Schematic representation for the secondary structure", secstruct_menu, SLOT(OnSchematicSecondaryStructure()));

    new CurrentMIGLWidgetAction("Hide Secondary Structure", "Removes the display of the secondary structure", secstruct_menu, SLOT(OnDeleteSecondaryStructure()));

    QMenu *secstruct_options_menu = secstruct_menu->addMenu("Options");
    new CurrentMIGLWidgetAction("Tube", "Set secondary structure options", secstruct_options_menu, SLOT(OnSecondaryStructureOptionsTube()));

    new CurrentMIGLWidgetAction("Beta sheet", "Set secondary structure options", secstruct_options_menu, SLOT(OnSecondaryStructureOptionsSheet()));

    new CurrentMIGLWidgetAction("Turn", "Set secondary structure options", secstruct_options_menu, SLOT(OnSecondaryStructureOptionsTurn()));

    new CurrentMIGLWidgetAction("Random coil", "Set secondary structure options", secstruct_options_menu, SLOT(OnSecondaryStructureOptionsRandom()));

    new CurrentMIGLWidgetAction("Helix", "Set secondary structure options", secstruct_options_menu, SLOT(OnSecondaryStructureOptionsHelix()));



    QMenu *showatoms_menu = show_menu->addMenu("Atoms");
    new CurrentMIGLWidgetAction("Show &all atoms", "Show all the atoms in the model", showatoms_menu, SLOT(OnObjectsAllatoms()));
    new CurrentMIGLWidgetAction("Show last &picked atom", "Show the last picked atom", showatoms_menu, SLOT(OnShowPickedatomTurnon()), SLOT(OnUpdateShowPickedatomTurnon(QAction*)));
    new CurrentMIGLWidgetAction("Show all p&icked atoms", "Show all the atoms on the stack", showatoms_menu, SLOT(OnShowAllpickedatomsTurnon()), SLOT(OnUpdateShowAllpickedatomsTurnon(QAction*)));
    showatoms_menu->addSeparator();
    new CurrentMIGLWidgetAction("&Hide all atoms", "Hide the whole model", showatoms_menu, SLOT(OnHideModel()), SLOT(OnUpdateHideModel(QAction*)));

    new CurrentMIGLWidgetAction("Hide last picked a&tom", "Hide the last picked atom", showatoms_menu, SLOT(OnShowPickedatomTurnoff()), SLOT(OnUpdateShowPickedatomTurnoff(QAction*)));
    new CurrentMIGLWidgetAction("Hide all picked at&oms", "Hide all the atoms on the stack", showatoms_menu, SLOT(OnShowAllpickedatomsTurnoff()), SLOT(OnUpdateShowAllpickedatomsTurnoff(QAction*)));
    miGLWidgetAction = new CurrentMIGLWidgetAction("Toggle h&ydrogens", "Toggle hydrogens on/off", showatoms_menu, SLOT(OnShowHidehydrogens()), SLOT(OnUpdateShowHidehydrogens(QAction*)));
    miGLWidgetAction->setCheckable(true);


    QMenu *annotation_menu = show_menu->addMenu("Annotation");
    miGLWidgetAction = new CurrentMIGLWidgetAction("A&dd annotation to model",
                                "Adds an annotation to the current model at the screen center",
                                annotation_menu,
                                SLOT(OnAnnotation()),
                                SLOT(OnUpdateAnnotation(QAction*)));
    miGLWidgetAction->setIcon(QIcon(annotate_xpm));
    toolBarActions += miGLWidgetAction;

    new CurrentMIGLWidgetAction("Move annotation to &picked atom",
                                "Move the last picked Annotation to the last picked atom ccords",
                                annotation_menu,
                                SLOT(OnMoveAnnotationAtom()),
                                SLOT(OnUpdateMoveAnnotationAtom(QAction*)));

    new CurrentMIGLWidgetAction("Move annotation to c&enter",
                                "Move the last picked Annotation to the screen center",
                                annotation_menu,
                                SLOT(OnMoveAnnotationCenter()),
                                SLOT(OnUpdateMoveAnnotationCenter(QAction*)));

    QMenu *surf_menu = show_menu->addMenu("&Dot Surface");
    fill_surf_menu(surf_menu);

    new CurrentMIGLWidgetAction("&Undo Viewpoint", "Use this command to center objects", view_menu, SLOT(OnViewUndo()), SLOT(OnUpdateViewUndo(QAction*)));
    new CurrentMIGLWidgetAction("&Save Viewpoint...", "Use this command to save the viewpoint to a file", view_menu, SLOT(OnViewSave()));
    new CurrentMIGLWidgetAction("&Load Viewpoint...", "Use this command to load the viewpoint from a file", view_menu, SLOT(OnViewLoad()));
    view_menu->addSeparator();
    miGLWidgetAction = new CurrentMIGLWidgetAction("&Top view", "View from top and drag clipping planes", view_menu, SLOT(OnViewTopview()), SLOT(OnUpdateViewTopview(QAction*)));
    miGLWidgetAction->setCheckable(true);
    miGLWidgetAction->setIcon(QIcon(topview_xpm));
    toolBarActions += miGLWidgetAction;
    miGLWidgetAction = new CurrentMIGLWidgetAction("&Fullscreen\tESC", "Canvas fills the entire screen", view_menu, SLOT(OnFullScreen()), SLOT(OnUpdateFullScreen(QAction*)));
    miGLWidgetAction->setCheckable(true);

    _hardwareStereoAction = view_menu->addAction(tr("Use &Hardware Stereo"), this, SLOT(OnHardwareStereo()));
    _hardwareStereoAction->setCheckable(true);
    connect(view_menu, SIGNAL(aboutToShow()), this, SLOT(OnUpdateHardwareStereo()));

    _stereoToggleAction = view_menu->addAction(tr("S&tereo\t|"), this, SLOT(OnStereoToggle()));
    _stereoToggleAction->setCheckable(true);
    connect(view_menu, SIGNAL(aboutToShow()), this, SLOT(OnUpdateStereoToggle()));

    view_menu->addSeparator();
    miGLWidgetAction = new CurrentMIGLWidgetAction("Slab In\tShift+I", "Decrease the distance beteween the front and back clipping planes", view_menu, SLOT(OnViewSlabin()));
    miGLWidgetAction->setIcon(QIcon(slabin_xpm));
    toolBarActions.insert(toolBarActions.size()-1, miGLWidgetAction);

    miGLWidgetAction = new CurrentMIGLWidgetAction("Slab Out\tShift+O", "Increase the distance beteween the front and back clipping planes", view_menu, SLOT(OnViewSlabout()));
    miGLWidgetAction->setIcon(QIcon(slabout_xpm));
    toolBarActions.insert(toolBarActions.size()-1, miGLWidgetAction);

    new CurrentMIGLWidgetAction("Set s&lab...", "Set the values of the front and back clipping planes", view_menu, SLOT(OnViewClipplanes()));
    miGLWidgetAction = new CurrentMIGLWidgetAction("Zoom In\tI", "Zoom viewpoint in by 20%", view_menu, SLOT(OnGotoZoomiin()));
    miGLWidgetAction->setIcon(QIcon(zoomin_xpm));
    toolBarActions.insert(toolBarActions.size()-1, miGLWidgetAction);
    miGLWidgetAction = new CurrentMIGLWidgetAction("Zoom Out\tO", "Zoom viewpoint in by 20%", view_menu, SLOT(OnGotoZoomout()));
    miGLWidgetAction->setIcon(QIcon(zoomout_xpm));
    toolBarActions.insert(toolBarActions.size()-1, miGLWidgetAction);
    miGLWidgetAction = new CurrentMIGLWidgetAction("Rotate View &+90", "Use this command (and -90) to center objects", view_menu, SLOT(OnViewRotatey90()));
    miGLWidgetAction->setIcon(QIcon(roty90_xpm));
    toolBarActions += miGLWidgetAction;

    miGLWidgetAction = new CurrentMIGLWidgetAction("Rotate View &-90", "Use this command to center objects", view_menu, SLOT(OnViewRotateyminus()));
    miGLWidgetAction->setIcon(QIcon(rotym90_xpm));
    toolBarActions += miGLWidgetAction;
    miGLWidgetAction = new CurrentMIGLWidgetAction("Orthonor&mal", "Set perspective to 0 (infinite viewpoint)", view_menu, SLOT(OnViewOrthonormal()), SLOT(OnUpdateViewOrthonormal(QAction*)));
    miGLWidgetAction->setCheckable(true);
    miGLWidgetAction = new CurrentMIGLWidgetAction("I&ncrease Perspective", "Increase perspective +20%", view_menu, SLOT(OnIncreasePersp()));
    miGLWidgetAction->setIcon(QIcon(incrper_xpm));
    toolBarActions += miGLWidgetAction;
    miGLWidgetAction = new CurrentMIGLWidgetAction("&Decrease Perspective", "Decrease perspective -20%", view_menu, SLOT(OnDecreasePersp()), SLOT(OnUpdateDecreasePersp(QAction*)));
    miGLWidgetAction->setIcon(QIcon(decrper_xpm));
    toolBarActions.insert(toolBarActions.size()-1, miGLWidgetAction);
    toolBarActions += 0;

    view_menu->addSeparator();
    new CurrentMIGLWidgetAction("&Center model on screen", "Center molecule on screen", view_menu, SLOT(OnGotoFittoscreen()));
    new CurrentMIGLWidgetAction("&Center All Models On Screen", "Center all the molecules on the screen", view_menu, SLOT(OnGotoFitalltoscreen()));
    new CurrentMIGLWidgetAction("Center Visible Density", "Center the density visible on the screen", view_menu, SLOT(OnMapCenterVisibleDensity()), SLOT(OnUpdateMapCenterVisibleDensity(QAction*)));

    new CurrentMIGLWidgetAction("&Go to x,y,z...", "Center view at a coordinate in Angstroms", view_menu, SLOT(OnGotoGotoxyz()));
    view_menu->addSeparator();

    QMenu *anime_menu = view_menu->addMenu("&Animate");
    miGLWidgetAction = new CurrentMIGLWidgetAction("&Rock", "Rock the model back and forth", anime_menu, SLOT(OnAnimateRock()), SLOT(OnUpdateAnimateRock(QAction*)));
    miGLWidgetAction->setCheckable(true);
    miGLWidgetAction = new CurrentMIGLWidgetAction("Ro&ll", "Roll the model around the vertical", anime_menu, SLOT(OnAnimateRoll()), SLOT(OnUpdateAnimateRoll(QAction*)));
    miGLWidgetAction->setCheckable(true);
    miGLWidgetAction = new CurrentMIGLWidgetAction("&Blink", "Blink between two or more models", anime_menu, SLOT(OnAnimateBlink()), SLOT(OnUpdateAnimateBlink(QAction*)));
    miGLWidgetAction->setCheckable(true);
    new CurrentMIGLWidgetAction("Rock/Roll Rates...", "Set the rates for rock and roll", anime_menu, SLOT(OnAnimateRockandrollparameters()));

    new CurrentMIGLWidgetAction("Measure &Distance",
                                "Measure the distance between last two picked atoms",
                                analyze_menu,
                                SLOT(OnGeometryDistance()),
                                SLOT(OnUpdateGeometryDistance(QAction*)));

    new CurrentMIGLWidgetAction("Measure &Angle",
                                "Measure the angle of last three picked atoms",
                                analyze_menu,
                                SLOT(OnGeometryAngle()),
                                SLOT(OnUpdateGeometryAngle(QAction*)));

    new CurrentMIGLWidgetAction("Measure D&ihedral",
                                "Measure the dihedral/torsion of last four picked atoms",
                                analyze_menu,
                                SLOT(OnGeometryTorsion()),
                                SLOT(OnUpdateGeometryTorsion(QAction*)));
    analyze_menu->addSeparator();
    miGLWidgetAction = new CurrentMIGLWidgetAction("Detailed Ramachandran plot", "Show allowed regions in Ramachandran plot, too", analyze_menu, SLOT(OnRamachandranPlotShowAllowed()), SLOT(OnUpdateRamachandranPlotShowAllowed(QAction*)));
    miGLWidgetAction->setCheckable(true);
    analyze_menu->addSeparator();
    new CurrentMIGLWidgetAction("A&nalyze Geometry Errors",
                                "Find all the geometry errors above a threshold in the model",
                                analyze_menu,
                                SLOT(OnFindGeomErrors()),
                                SLOT(OnUpdateFindGeomErrors(QAction*)));

    new CurrentMIGLWidgetAction("Clear &Geometry Errors",
                                "Clear all the Annotation of geometry errors in the model",
                                analyze_menu,
                                SLOT(OnClearGeomAnnotations()));

    analyze_menu->addSeparator();
    new CurrentMIGLWidgetAction("Add &H-bond",
                                "Build an H-bond between the last two picked atoms",
                                analyze_menu,
                                SLOT(OnGeomAddsinglehbond()),
                                SLOT(OnUpdateGeomAddsinglehbond(QAction*)));

    new CurrentMIGLWidgetAction("B&uild H-Bonds",
                                "Build all the H-bonds in the model",
                                analyze_menu,
                                SLOT(OnGeomHbonds()));

    new CurrentMIGLWidgetAction("&Clear H-bonds",
                                "Clear all the H-bonds in the model",
                                analyze_menu,
                                SLOT(OnGeomClearhbonds()));

    analyze_menu->addSeparator();

    new CurrentMIGLWidgetAction("Bui&ld Contacts",
                                "Build all the Contacts in the model",
                                analyze_menu,
                                SLOT(OnGeomNeighbours()),
                                SLOT(OnUpdateGeomNeighbours(QAction*)));

    new CurrentMIGLWidgetAction("Cl&ear Contacts",
                                "Clear all the Contacts in the model",
                                analyze_menu,
                                SLOT(OnGeomClearneighbours()));

    // Refine menu actions which have shortcuts must be immediately updated
    // when refinement state changes. They can't wait for the update which occurs
    // when the menu displays.
    refineResidueAction = new CurrentMIGLWidgetAction("R&efine Residue\tCtrl+R", "Real-space refine the last picked residue", refi_menu, SLOT(OnRefiResidue()), SLOT(OnUpdateRefiResidue(QAction*)));

    connect(MIFitGeomRefiner(), SIGNAL(isRefiningChanged(bool)),
            this, SLOT(updateIsRefining(bool)));
    new CurrentMIGLWidgetAction("Re&fine Local Region", "Real-space refine the last picked residue and its 2 neighbours", refi_menu, SLOT(OnRefiRegion()), SLOT(OnUpdateRefiRegion(QAction*)));

    new CurrentMIGLWidgetAction("Ref&ine Range", "Real-space refine the last 2 picks and the intervening residues", refi_menu, SLOT(OnRefiRange()), SLOT(OnUpdateRefiRange(QAction*)));

    new CurrentMIGLWidgetAction("Refi&ne Molecule", "Real-space refine the entire molecule in all 6 dimensions", refi_menu, SLOT(OnRefiMolecule()), SLOT(OnUpdateRefiMolecule(QAction*)));

    refi_menu->addSeparator();
    new CurrentMIGLWidgetAction("Ri&gid-Body Refine Current Atoms", "Rigid Body Refine the current atoms (cyan color)", refi_menu, SLOT(OnRefiRigidBody()), SLOT(OnUpdateRefiRigidBody(QAction*)));

    new CurrentMIGLWidgetAction("Fin&d Ligand Fit and Conformer", "Search for the best conformer ligand fit and torsion angles", refi_menu, SLOT(OnRefiLigandFit()), SLOT(OnUpdateRefiLigandFit(QAction*)));

    refi_menu->addSeparator();
    acceptRefineAction = new CurrentMIGLWidgetAction("&Accept Refine\tCtrl+Shift+R", "Accept the refinement and finalize the atoms positions", refi_menu, SLOT(OnRefiAccept()), SLOT(OnUpdateRefiAccept(QAction*)));

    new CurrentMIGLWidgetAction("&Reset Refine", "Reset the current real-space refinement by putting atoms back where they were", refi_menu, SLOT(OnRefiReset()), SLOT(OnUpdateRefiReset(QAction*)));

    new CurrentMIGLWidgetAction("&Cancel Refine", "Cancel the current real-space refinement by putting atoms back where they were", refi_menu, SLOT(OnRefiCancel()), SLOT(OnUpdateRefiCancel(QAction*)));

    refi_menu->addSeparator();
    new CurrentMIGLWidgetAction("&Undo Refine", "", refi_menu, SLOT(OnRefiUndo()), SLOT(OnUpdateRefiUndo(QAction*)));

    new CurrentMIGLWidgetAction("Re&do Refine", "", refi_menu, SLOT(OnRefiReDo()), SLOT(OnUpdateRefiRedo(QAction*)));

    refi_menu->addSeparator();
    new CurrentMIGLWidgetAction("Refine &Options", "Real-space refine options", refi_menu, SLOT(OnRefiOptions()));


    miGLWidgetAction = new CurrentMIGLWidgetAction("Fit R&esidue\tF", "Set the last picked residue for fitting", fit_menu, SLOT(OnFitResidue()), SLOT(OnUpdateFitResidue(QAction*)));
    miGLWidgetAction->setCheckable(true);
    miGLWidgetAction->setIcon(QIcon(fit_xpm));
    toolBarActions += miGLWidgetAction;

    miGLWidgetAction = new CurrentMIGLWidgetAction("Fit A&tom", "Set the last picked atom for fitting", fit_menu, SLOT(OnFitSingleatom()), SLOT(OnUpdateFitSingleatom(QAction*)));
    miGLWidgetAction->setCheckable(true);

    miGLWidgetAction = new CurrentMIGLWidgetAction("Fit Re&sidues", "Set all the picked residues for fitting", fit_menu, SLOT(OnFitResidues()), SLOT(OnUpdateFitResidues(QAction*)));
    miGLWidgetAction->setCheckable(true);

    miGLWidgetAction = new CurrentMIGLWidgetAction("Fit Res&idue Range", "Fit the range marked by the top two residues on the stack", fit_menu, SLOT(OnFitRange()), SLOT(OnUpdateFitRange(QAction*)));
    miGLWidgetAction->setCheckable(true);

    miGLWidgetAction = new CurrentMIGLWidgetAction("Fit At&oms", "Set all the picked atoms for fitting", fit_menu, SLOT(OnFitAtoms()), SLOT(OnUpdateFitAtoms(QAction*)));
    miGLWidgetAction->setCheckable(true);

    miGLWidgetAction = new CurrentMIGLWidgetAction("Fit &Molecule", "Set the selected molecule for fitting", fit_menu, SLOT(OnFitFitmolecule()), SLOT(OnUpdateFitFitmolecule(QAction*)));
    miGLWidgetAction->setCheckable(true);

    fit_menu->addSeparator();
    miGLWidgetAction = new CurrentMIGLWidgetAction("&Accept Fit\t;", "Accept fitting, leave atoms at new positions", fit_menu, SLOT(OnFitApply()), SLOT(OnUpdateFitApply(QAction*)));
    miGLWidgetAction->setIcon(QIcon(apply_xpm));
    toolBarActions += miGLWidgetAction;

    new CurrentMIGLWidgetAction("&Reset Fit", "Put atoms back to start and continue fitting", fit_menu, SLOT(OnFitReset()), SLOT(OnUpdateFitReset(QAction*)));

    miGLWidgetAction = new CurrentMIGLWidgetAction("&Cancel Fit", "Stop fitting and put atoms back to starting positions", fit_menu, SLOT(OnFitCancel()), SLOT(OnUpdateFitCancel(QAction*)));
    miGLWidgetAction->setIcon(QIcon(cancel_xpm));
    toolBarActions += miGLWidgetAction;

    fit_menu->addSeparator();
    miGLWidgetAction = new CurrentMIGLWidgetAction("&Rotate", "Set right mouse to rotate mode", fit_menu, SLOT(OnFitRotate()), SLOT(OnUpdateFitRotate(QAction*)));
    miGLWidgetAction->setCheckable(true);
    miGLWidgetAction->setIcon(QIcon(rotate_xpm));
    toolBarActions += miGLWidgetAction;
    miGLWidgetAction = new CurrentMIGLWidgetAction("&Translate", "Set right mouse to translate mode", fit_menu, SLOT(OnFitTranslate()), SLOT(OnUpdateFitTranslate(QAction*)));
    miGLWidgetAction->setCheckable(true);
    miGLWidgetAction->setIcon(QIcon(translate_xpm));
    toolBarActions += miGLWidgetAction;
    miGLWidgetAction = new CurrentMIGLWidgetAction("&Torsion", "Set right mouse to torsion mode", fit_menu, SLOT(OnFitTorsion()), SLOT(OnUpdateFitTorsion(QAction*)));
    miGLWidgetAction->setCheckable(true);
    miGLWidgetAction->setIcon(QIcon(torsion_xpm));
    toolBarActions += miGLWidgetAction;
    miGLWidgetAction = new CurrentMIGLWidgetAction("C&enter", "Set right mouse to move screen center (default mode)", fit_menu, SLOT(OnFitCentermode()), SLOT(OnUpdateFitCentermode(QAction*)));
    miGLWidgetAction->setCheckable(true);
    miGLWidgetAction->setIcon(QIcon(center_xpm));
    toolBarActions += miGLWidgetAction;
    fit_menu->addSeparator();
    new CurrentMIGLWidgetAction("&Set Up Torsion", "Torsion about last two atoms picked, second end moves", fit_menu, SLOT(OnFitSetuptorsion()), SLOT(OnUpdateFitSetuptorsion(QAction*)));
    new CurrentMIGLWidgetAction("C&lear Torsion", "Clear Torsion, leave in fitting mode", fit_menu, SLOT(OnFitCleartorsion()), SLOT(OnUpdateFitCleartorsion(QAction*)));
    fit_menu->addSeparator();
    QMenu *fitsurf_menu = fit_menu->addMenu("Surface &Fit Atoms");
    miGLWidgetAction = new CurrentMIGLWidgetAction("&Van Der Waal Surface", "Surface the current atoms (fit) with a van der Waal surface", fitsurf_menu, SLOT(OnFitSurfaceVdw()), SLOT(OnUpdateFitSurfaceVdw(QAction*)));
    miGLWidgetAction->setCheckable(true);

    miGLWidgetAction = new CurrentMIGLWidgetAction("&Extended Surface", "Surface the current atoms (fit) with an extended (2x) surface", fitsurf_menu, SLOT(OnFitSurfaceExtended()), SLOT(OnUpdateFitSurfaceExtended(QAction*)));
    miGLWidgetAction->setCheckable(true);

    miGLWidgetAction = new CurrentMIGLWidgetAction("&Contact Surface", "Surface the current atoms (fit) wherever there is a contact", fitsurf_menu, SLOT(OnFitSurfaceProbe()), SLOT(OnUpdateFitSurfaceProbe(QAction*)));
    miGLWidgetAction->setCheckable(true);

    miGLWidgetAction = new CurrentMIGLWidgetAction("&No Surface", "Turn off the surface for the current atoms", fitsurf_menu, SLOT(OnFitSurfaceNone()), SLOT(OnUpdateFitSurfaceNone(QAction*)));
    miGLWidgetAction->setCheckable(true);

    QMenu *pentamer_menu = fit_menu->addMenu("Fix &Backbone");
    new CurrentMIGLWidgetAction("&Suggest Backbone Match", "Find a pentamer that best matches the CA and CB atoms (if present)", pentamer_menu, SLOT(OnFindPentamer()), SLOT(OnUpdateFindPentamer(QAction*)));

    new CurrentMIGLWidgetAction("&Replace Middle 3", "Replace the backbone atoms of the middle 3 residues", pentamer_menu, SLOT(OnReplaceMiddle3()), SLOT(OnUpdateReplaceMiddle3(QAction*)));

    new CurrentMIGLWidgetAction("R&eplace First 4", "Replace the backbone atoms of the first 4 residues", pentamer_menu, SLOT(OnReplaceFirst4()), SLOT(OnUpdateReplaceFirst4(QAction*)));

    new CurrentMIGLWidgetAction("Re&place Last 4", "Replace the backbone atoms of the last 4 residues", pentamer_menu, SLOT(OnReplaceLast4()), SLOT(OnUpdateReplaceLast4(QAction*)));

    new CurrentMIGLWidgetAction("Rep&lace All 5", "Replace the backbone atoms of all 5 residues", pentamer_menu, SLOT(OnReplaceAll()), SLOT(OnUpdateReplaceAll(QAction*)));

    //  Removed the Clear backbone match menu entry because it failed to update the tree properly, leading to a crash
    new CurrentMIGLWidgetAction("&Clear Backbone Match", "Clear the backbone pentamer from the screen (and memory)", pentamer_menu, SLOT(OnClearPentamer()), SLOT(OnUpdateClearPentamer(QAction*)));

    new CurrentMIGLWidgetAction("&Flip Peptide", "Flip petide bond to put carbonyl on opposite side", pentamer_menu, SLOT(OnFlipPeptide()), SLOT(OnUpdateFlipPeptide(QAction*)));



    QMenu *disorder_menu = fit_menu->addMenu("&Disorder");
    new CurrentMIGLWidgetAction("Split at &Torsion", "Adds disorder (A and B parts) at a torsion angle", disorder_menu, SLOT(OnFitSplitTorsion()), SLOT(OnUpdateFitSplitTorsion(QAction*)));

    new CurrentMIGLWidgetAction("Split &Fit", "duplicates residues being fit making an A and B part", disorder_menu, SLOT(OnFitSplitFit()), SLOT(OnUpdateFitSplitFit(QAction*)));

    fit_menu->addSeparator();
    new CurrentMIGLWidgetAction("&Undo Fit", "Undo Last fit", fit_menu, SLOT(OnFitUndo()), SLOT(OnUpdateFitUndo(QAction*)));
    new CurrentMIGLWidgetAction("&Redo Fit", "Undo the Undo Last fit", fit_menu, SLOT(OnFitRedo()), SLOT(OnUpdateFitRedo(QAction*)));
    fit_menu->addSeparator();
    new CurrentMIGLWidgetAction("S&uperimpose...", "Superimpose models by LSQ fitting", fit_menu, SLOT(OnFitLsqsuperpose()));

    new CurrentMIGLWidgetAction("&New Model...", "Adds a new blank model", model_menu, SLOT(OnNewModel()));
    model_menu->addSeparator();
    new CurrentMIGLWidgetAction("&Add residue", "Add or insert a new residue after the last pick", model_menu, SLOT(OnFitInsertresidue()), SLOT(OnUpdateFitInsertresidue(QAction*)));
    new CurrentMIGLWidgetAction("&Replace residue", "Replace a residue with another type", model_menu, SLOT(OnFitReplacewith()), SLOT(OnUpdateFitReplacewith(QAction*)));
    new CurrentMIGLWidgetAction("Replace and &Fit\tR", "Replace a residue with another type and search confomers", model_menu, SLOT(OnFitReplaceAndFit()), SLOT(OnUpdateFitReplaceAndFit(QAction*)));
    new CurrentMIGLWidgetAction("Replace with &Sequence", "Replace a residue with the lower sequence", model_menu, SLOT(OnReplaceSequence()), SLOT(OnUpdateReplaceSequence(QAction*)));
    new CurrentMIGLWidgetAction("Next &Conformer", "Replace a residue with another type", model_menu, SLOT(OnFitNextConfomer()), SLOT(OnUpdateFitNextConfomer(QAction*)));
    new CurrentMIGLWidgetAction("&Delete residue\tD", "Delete the last picked residue", model_menu, SLOT(OnFitDeleteresidue()), SLOT(OnUpdateFitDeleteresidue(QAction*)));
    new CurrentMIGLWidgetAction("R&ename residue", "Rename(renumber) the last picked residue", model_menu, SLOT(OnFitRenameresidue()), SLOT(OnUpdateFitRenameresidue(QAction*)));
    new CurrentMIGLWidgetAction("&Delete atom\tDelete", "Delete the last picked atom", model_menu, SLOT(OnDeleteAtom()), SLOT(OnUpdateDeleteAtom(QAction*)));
    new CurrentMIGLWidgetAction("Add &Waters...", "Sprinkle water throughout map", model_menu, SLOT(OnAddWater()), SLOT(OnUpdateAddWater(QAction*)));
    model_menu->addSeparator();
    new CurrentMIGLWidgetAction("&Bond",
                                "Force a bond between the last two picked atoms",
                                model_menu,
                                SLOT(OnGeomBond()),
                                SLOT(OnUpdateGeomBond(QAction*)));
    new CurrentMIGLWidgetAction("B&reak Bond",
                                "Break the bond between the last two picked atoms",
                                model_menu,
                                SLOT(OnGeomUnbond()),
                                SLOT(OnUpdateGeomUnbond(QAction*)));
    model_menu->addSeparator();
    new CurrentMIGLWidgetAction("Add MRK &Before\t<", "Adds a MRK (C-alpha) before the focus residue for tracing the chain", model_menu, SLOT(OnAddMarkBefore()), SLOT(OnUpdateAddMarkBefore(QAction*)));
    new CurrentMIGLWidgetAction("Add MRK &After\t>", "Adds a MRK (C-alpha) after the focus residue for tracing the chain", model_menu, SLOT(OnAddMarkAfter()), SLOT(OnUpdateAddMarkAfter(QAction*)));
    new CurrentMIGLWidgetAction("Poly-Ala &Range", "Turns the model into poly-alanine between the two picks", model_menu, SLOT(OnPolyAla()), SLOT(OnUpdatePolyAla(QAction*)));
    new CurrentMIGLWidgetAction("Poly-Ala &Chain", "Turns the chain containing the last pick into poly-alinine", model_menu, SLOT(OnPolyAlaChain()), SLOT(OnUpdatePolyAlaChain(QAction*)));

    new CurrentMIGLWidgetAction("Go To &N-terminus\t[", "Follows the current chain to its N-terminal end", model_menu, SLOT(OnGotoNter()), SLOT(OnUpdateGotoNter(QAction*)));
    new CurrentMIGLWidgetAction("Go To &C-terminus\t]", "Follows the current chain to its C-terminal end", model_menu, SLOT(OnGotoCter()), SLOT(OnUpdateGotoCter(QAction*)));
    model_menu->addSeparator();
    miGLWidgetAction = new CurrentMIGLWidgetAction("&Checkpoint Model", "Checkpoints the model - use to save before", model_menu, SLOT(OnCheckpointModel()), SLOT(OnUpdateCheckpointModel(QAction*)));
    miGLWidgetAction->setIcon(QIcon(chekpoin_xpm));
    toolBarActions.insert(7, miGLWidgetAction);

    new CurrentMIGLWidgetAction("&Revert Model...", "Reverts the model to an earlier version saved", model_menu, SLOT(OnRevertModel()), SLOT(OnUpdateRevertModel(QAction*)));
    miGLWidgetAction = new CurrentMIGLWidgetAction("Auto-checkpoint Model", "Automatically checkpoints the model every few minutes", model_menu, SLOT(OnAutoCheckpointModel()), SLOT(OnUpdateAutoCheckpointModel(QAction*)));
    miGLWidgetAction->setCheckable(true);


    di_menu->addAction(tr("Load &New Dictionary..."), this, SLOT(OnLoadDictReplace()));
    di_menu->addAction(tr("Load & &Append Dictionary..."), this, SLOT(OnLoadDictAppend()));

    QMenu *ligImportMenu = di_menu->addMenu("&Import Ligand");
    ligImportMenu->addAction(tr("&CIF"), this, SLOT(OnLoadLigCif()));
    ligImportMenu->addAction(tr("&Mol"), this, SLOT(OnLoadLigMol()));
    ligImportMenu->addAction(tr("&Pdb"), this, SLOT(OnLoadLigPdb()));
    ligImportMenu->addAction(tr("&Smiles"), this, SLOT(OnLoadLigSmi()));

    di_menu->addAction(tr("&Save Dictionary..."), this, SLOT(OnSaveDict()));
    di_menu->addAction(tr("&Edit Residue..."), this, SLOT(OnEditDictResidue()));


    JobManager = new BatchJobManager();
    QTimer::singleShot(0, this, SLOT(fillJobMenu()));


    hide_menu = new QMenu(this);
    new CurrentMIGLWidgetAction("Whole Model",
                                "Hide the whole model",
                                hide_menu,
                                SLOT(OnHideModel()),
                                SLOT(OnUpdateHideModel(QAction*)));
    new CurrentMIGLWidgetAction("Hide last picked residue",
                                "Hide the last picked residue",
                                hide_menu,
                                SLOT(OnObjectResidueTurnoff()),
                                SLOT(OnUpdateObjectResidueTurnoff(QAction*)));
    new CurrentMIGLWidgetAction("Hide all picked residues",
                                "Hide all the residues on the stack",
                                hide_menu,
                                SLOT(OnObjectResiduesTurnoff()),
                                SLOT(OnUpdateObjectResiduesTurnoff(QAction*)));
    new CurrentMIGLWidgetAction("Hide residue range",
                                "Hide the residue at the top of stack",
                                hide_menu,
                                SLOT(OnObjectResiduerangeTurnoff()),
                                SLOT(OnUpdateObjectResiduerangeTurnoff(QAction*)));

    new CurrentMIGLWidgetAction("Hide Last Picked Atom",
                                "Hide last the atom on the stack",
                                hide_menu,
                                SLOT(OnShowPickedatomTurnoff()),
                                SLOT(OnUpdateShowPickedatomTurnoff(QAction*)));
    new CurrentMIGLWidgetAction("Hide All Picked Atoms",
                                "Hide all the atoms on the stack",
                                hide_menu,
                                SLOT(OnShowAllpickedatomsTurnoff()),
                                SLOT(OnUpdateShowAllpickedatomsTurnoff(QAction*)));




    windowMenu = menuBar()->addMenu(tr("&Window"));
    updateWindowMenu();
    connect(windowMenu, SIGNAL(aboutToShow()), this, SLOT(updateWindowMenu()));

    viewMenu = menuBar()->addMenu(tr("View"));

    QMenu *help_menu = menuBar()->addMenu(tr("&Help"));
    help_menu->addAction(tr("&Help..."), this, SLOT(OnHelp()));
    help_menu->addAction(tr("&About..."), this, SLOT(OnAbout()));
    help_menu->addAction(tr("&OpenGL format..."), this, SLOT(OnGLFormat()));

    foreach (QAction *action, toolBarActions)
    {
        if (action)
            tool_bar->addAction(action);
        else
            tool_bar->addSeparator();
    }
}


void MIMainWindow::createToolBars()
{
}


void MIMainWindow::createStatusBar()
{
    statusBar()->showMessage(tr("Ready"));
}

MIGLWidget*MIMainWindow::currentMIGLWidget()
{
    QMdiSubWindow *activeSubWindow = mdiArea->currentSubWindow();
    if (activeSubWindow)
    {
        return dynamic_cast<MIGLWidget*>(activeSubWindow->widget());
    }
    return 0;
}

QMdiArea*MIMainWindow::getMdiArea()
{
    return mdiArea;
}

void MIMainWindow::setActiveMIGLWidget(MIGLWidget *w)
{
    foreach (QMdiSubWindow *window, mdiArea->subWindowList())
    {
        if (window->widget() == w)
            mdiArea->setActiveSubWindow(window);
    }
}


QMdiSubWindow*MIMainWindow::findMIGLWidget(const QString &fileName)
{
    QString canonicalFilePath = QFileInfo(fileName).canonicalFilePath();

    foreach (QMdiSubWindow *window, mdiArea->subWindowList())
    {
        MIGLWidget *mdiChild = dynamic_cast<MIGLWidget*>(window->widget());
        if (mdiChild && QString(mdiChild->GetFilename().c_str()) == canonicalFilePath)
            return window;
    }
    return 0;
}

void MIMainWindow::switchLayoutDirection()
{
    if (layoutDirection() == Qt::LeftToRight)
        qApp->setLayoutDirection(Qt::RightToLeft);
    else
        qApp->setLayoutDirection(Qt::LeftToRight);
}

void MIMainWindow::setActiveSubWindow(QWidget *window)
{
    if (!window)
        return;
    mdiArea->setActiveSubWindow(qobject_cast<QMdiSubWindow*>(window));
}

void MIMainWindow::OnNew()
{
    newWindow();
}


void MIMainWindow::openRecentFile()
{
    QAction *action = qobject_cast<QAction*>(sender());
    if (action)
    {
        std::vector<std::string> fnames;
        fnames.push_back(action->data().toString().toStdString());
        OpenFiles(fnames);
    }
}

void MIMainWindow::addRecentFileActions(QMenu *fileMenu)
{
    fileMenu->addSeparator();
    for (int i = 0; i < MaxRecentFiles; ++i)
        fileMenu->addAction(recentFileActs[i]);
    fileMenu->addSeparator();
}

void MIMainWindow::setCurrentFile(const std::string &fname)
{
    QString fileName(fname.c_str());

    QSettings *settings = MIGetQSettings(); // could use MIConfig (it's the same file), but this api is easier here
    QStringList files = settings->value("recentFileList").toStringList();
    files.removeAll(fileName);
    files.prepend(fileName);
    while (files.size() > MaxRecentFiles)
        files.removeLast();
    settings->setValue("recentFileList", files);
    updateRecentFileActions();
}

void MIMainWindow::updateRecentFileActions()
{
    QSettings *settings = MIGetQSettings(); // could use MIConfig (it's the same file), but this api is easier here
    QStringList files = settings->value("recentFileList").toStringList();

    int numRecentFiles = qMin(files.size(), (int)MaxRecentFiles);

    for (int i = 0; i < numRecentFiles; ++i)
    {
        QString text = QString("&%1 %2").arg(i + 1).arg(files[i]);
        recentFileActs[i]->setText(text);
        recentFileActs[i]->setData(files[i]);
        recentFileActs[i]->setVisible(true);
    }
    for (int j = numRecentFiles; j < MaxRecentFiles; ++j)
        recentFileActs[j]->setVisible(false);
}

RamaPlotMgr*MIMainWindow::RamaPlotManager()
{
    return ramaPlotMgr;
}

//some shortcuts for MIMainWindow::instance()->Log, etc
void MIMainWindowLog(const std::string &str)
{
    return MIMainWindow::instance()->Log(str);
}

void MIMainWindowDebug(const std::string &str)
{
    return MIMainWindow::instance()->Debug(str);
}

void MIMainWindowLeftFooter(const std::string &str, int timeout)
{
    MIMainWindow::instance()->LeftFooter(str, timeout);
}

void MIMainWindowMiddleFooter(const std::string &str)
{
    MIMainWindow::instance()->MiddleFooter(str);
}

void MIMainWindowRightFooter(const std::string &str)
{
    MIMainWindow::instance()->RightFooter(str);
}

void MIMainWindow::updateShowMenu()
{
    if (solidSurfMenuAction)
    {
        show_menu->removeAction(solidSurfMenuAction);
        solidSurfMenuAction = NULL;
    }
    MIGLWidget *view = currentMIGLWidget();
    if (view)
    {
        solidSurfMenuAction = view->solidSurfMenuAction();
    }
    else
    {
        solidSurfMenuAction = new QAction(this);
        solidSurfMenuAction->setText(tr("Solid surface"));
    }
    show_menu->insertAction(canvas_menu->menuAction(), solidSurfMenuAction);
}

void MIMainWindow::updateIsRefining(bool isRefining)
{
    refineResidueAction->setEnabled(!isRefining);
    acceptRefineAction->setEnabled(isRefining);
}

void MIMainWindow::setRenderLineThickness(int thickness)
{
    if (currentMIGLWidget())
        currentMIGLWidget()->setViewpointLineThickness(thickness);
}

void MIMainWindow::updateRenderMenu()
{
    bool lineThicknessEnabled = false;
    if (currentMIGLWidget())
    {
        ViewPoint *viewpoint = currentMIGLWidget()->viewpoint;
        int thickness = viewpoint->GetLineThickness();
        if (thickness > 0 && renderLineThicknessMenu_->actions().size() <= thickness)
        renderLineThicknessMenu_->actions().at(thickness-1)->setChecked(true);
        lineThicknessEnabled = viewpoint->GetBallandStick() == ViewPoint::BALLANDSTICK
                  || viewpoint->GetBallandStick() == ViewPoint::STICKS;
    }
    renderLineThicknessMenu_->setEnabled(lineThicknessEnabled);
}

void MIMainWindow::fillJobMenu()
{
    createLocalSocketScript();
    Tools::instance().FillToolsMenu(_jobMenu);
}

QMenu *MIMainWindow::jobMenu() const
{
    return _jobMenu;
}

void MIMainWindow::saveJobMenu()
{
    GetJobManager()->saveJobMenu(_jobMenu);
}

void MIMainWindow::createLocalSocketScript()
{
    _scriptSocket = new LocalSocketScript(this);
}

QString MIMainWindow::scriptPort() const
{
    if (_scriptSocket)
        return _scriptSocket->name();
    return QString::null;
}

void MIMainWindow::addJob(const QString &menuName, const QString &jobName, const QString &executable, const QStringList &arguments, const QString &workingDirectory)
{
    QAction *jobAction = GetJobManager()->customJobAction(menuName, jobName, executable, arguments, workingDirectory);
    jobMenu()->addAction(jobAction);
}

void MIMainWindow::runPythonScript(const std::string &file)
{
    BatchJob *pythonJob = GetJobManager()->CreateJob();
    pythonJob->setProgram(BatchJobManager::pythonExe());
    pythonJob->setArguments(file.c_str());

    pythonJob->StartJob();
}
