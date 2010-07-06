#ifndef MAINWINDOW_H
#define MAINWINDOW_H

#include <QMainWindow>

class QAction;
class QMenu;
class QMdiArea;
class QMdiSubWindow;
class MIGLWidget;
class QSignalMapper;
class QListWidget;
class QTextEdit;
class QDockWidget;
class QLabel;
class QCursor;

class GLOverviewCanvas;
class BatchJobManager;
class RamaPlotMgr;
class DependantActions;
class LocalSocketScript;

namespace chemlib
{
    class MIMolInfo;
}

//cursor constants
const unsigned int imhCross = 0; // Cross cursor used to pick on canvas.
const unsigned int imhTranslate = 1; // Translate cursor used to indicate a translation operation.
const unsigned int imhRotate = 2; // Rotate cursor used to indicate a rotation operation.
const unsigned int imhZRotate = 3; // Z rotation cursor used to indicate a Z-rotation operation.
const unsigned int imhZCursor = 4; // Indicates cursor in the z rotation area which ids the top 30 pixels of the canvas.
const unsigned int imhTorsion = 5; // Torsion operation cusor indicating a bond is being twisted.
const unsigned int imhCenter = 6; // Centering operation cursor.
const unsigned int imhScale = 7; // Scale/zoom operation cursor.
const unsigned int imhSlab = 8; // Slab operation cursor indicating the slab plane is being moved.
const unsigned int imhWait1 = 9; // Animated wait cursor.
const unsigned int imhWait2 = 10; // Animated wait cursor.
const unsigned int imhWait3 = 11; // Animated wait cursor.
const unsigned int imhWait4 = 12; // Animated wait cursor.
const unsigned int imhWait5 = 13; // Animated wait cursor.
const unsigned int imhWait6 = 14; // Animated wait cursor.
const unsigned int imhSlabDrag = 15; // Slab drag operation
const unsigned int imhCount = 16; // Number of cursors


class ViewSyncedPanel;

class MIMainWindow : public QMainWindow
{
    Q_OBJECT

    static MIMainWindow *_instance;
    MIMainWindow();

    QDockWidget *modelsDock;
    QDockWidget *displayDock;
    QDockWidget *jobsDock;
    QDockWidget *logDock;
    QDockWidget *navigatorDock;
    QDockWidget *ramaDock;

    QMenu *canvas_menu;
    DependantActions *refineDependantActions;

public:
    static MIMainWindow *instance();
    ~MIMainWindow();

public slots:
    void Log(const std::string &str);
    void Debug(const std::string &str);
    void RightFooter(const std::string &str);
    void MiddleFooter(const std::string &str);
    void LeftFooter(const std::string &str, int timeout = 0);
    void updateToolBar();
    void addJob(const QString &menuName, const QString &jobName,
                const QString &executable, const QStringList &arguments,
                const QString &workingDirectory);

    void runPythonScript(const std::string &file);

public:
    MIGLWidget *currentMIGLWidget();
    QMdiArea *getMdiArea();
    void setActiveMIGLWidget(MIGLWidget *w);

    QDockWidget *AddAsDockWidget(QWidget *w, const std::string &name, Qt::DockWidgetArea area);

    // if widget is 0, set on the whole app
    // if id < 0, restore default cursor
    void SetCursor(int id, QWidget *widget = 0);

    void fill_surf_menu(QMenu*);

    void OpenFiles(const std::vector<std::string> &files, bool newWindow = false);

    BatchJobManager *GetJobManager();
    bool isJobLimit();

    void updateNavigator();
    RamaPlotMgr *RamaPlotManager();

    ViewSyncedPanel *GetModelsTree()
    {
        return modelsView;
    }

    void addRecentFileActions(QMenu*);
    void setCurrentFile(const std::string &fname);

    QMenu *jobMenu() const;
    void saveJobMenu();
    QString scriptPort() const;

    void showPreferences(int page = 0);

private slots:
    void OnClose();
    void OnExit();
    void OnNew();
    void OnFileOpen();
    void OnFileOpenNew();

    void OnFileSave();
    void OnFileSaveAs();
    void OnExportModel();
    void OnPrint();
    void OnEditCopy();
    void OnExportImage();

    void OnScript();

    //dictionary menu handlers
    void OnSaveDict();
    void OnLoadLigCif();
    void OnLoadLigMol();
    void OnLoadLigPdb();
    void OnLoadLigSmi();
    void OnLoadDictAppend();
    void OnLoadDictReplace();
    void OnEditDictResidue();

    void OnBackgroundColor();
    void OnSideChainTool();
    void OnHideTool();
    void OnColorTool();
    void OnShowTool();

    void OnUpdateStereoToggle();
    void OnUpdateHardwareStereo();
    void OnStereoToggle();
    void OnHardwareStereo();
    void OnAbout();
    void OnGLFormat();
    void OnHelp();
    void OnMenuValidate();
    void OnDefineAtomColors();
    void OnDefineBValueColors();
    void OnManageCrystals();
    void OnPreferences();

    void fileOpen();
    void newWindow();
    void updateMenus();
    void childActivated(QMdiSubWindow*);
    void updateWindowMenu();
    QWidget *createMIGLWidget();
    void switchLayoutDirection();
    void setActiveSubWindow(QWidget *window);
    void AfterInit();

    void openRecentFile();

    void updateShowMenu();
    void updateIsRefining(bool isRefining);

    void updateFileMenu();

    void setRenderLineThickness(int thickness);
    void updateRenderMenu();

    void fillJobMenu();

    void createLocalSocketScript();

private:
    enum { MaxRecentFiles = 6 };
    QAction *recentFileActs[MaxRecentFiles];
    void updateRecentFileActions();

    std::string OnLoadLigand(std::string wildcard, std::string filename, std::string smiles, std::string code);
    void ShowDictEditor(const char *type);
    bool LoadSmiles(std::string &smiles, std::string &tlcode, chemlib::MIMolInfo &mi);

    void closeEvent(QCloseEvent *event);
    void saveLayout();

    void createActions();
    void createMenus();
    void createToolBars();
    void createStatusBar();
    void createDockWindows();

    QMdiSubWindow *findMIGLWidget(const QString &fileName);

    QToolBar *tool_bar;
    QToolBar *displayToolBar;

    QMdiArea *mdiArea;
    QSignalMapper *windowMapper;

    QTextEdit *logWindow;

    QMenu *side_menu;
    QMenu *hide_menu;
    QMenu *color_menu;
    QMenu *view_menu;
    QMenu *show_menu;
    QMenu *render_menu;
    QMenu *model_menu;
    QMenu *fit_menu;
    QMenu *refi_menu;
    QMenu *analyze_menu;
    QAction *solidSurfMenuAction;
    QMenu *windowMenu;
    QMenu *viewMenu;

    // file menu actions;
    QAction *fileNewAction;
    QAction *fileOpenAction;
    QAction *fileSaveAction;
    QAction *fileSaveAsAction;
    QAction *exportModelAction;
    QAction *filePrintAction;
    QAction *copyCanvasAction;
    QAction *exportImageAction;
    QAction *closeAction;

    // window menu actions
    QAction *closeAct;
    QAction *closeAllAct;
    QAction *tileAct;
    QAction *cascadeAct;
    QAction *nextAct;
    QAction *previousAct;

    // edit menu actions
    //QAction *cutAct;
    //QAction *copyAct;
    //QAction *pasteAct;

    QAction *refineResidueAction;
    QAction *acceptRefineAction;

    QLabel *middleFooter;
    QLabel *rightFooter;

    QCursor *cursors[imhCount];


    GLOverviewCanvas *navigator;

    BatchJobManager *JobManager;
    ViewSyncedPanel *modelsView;
    RamaPlotMgr *ramaPlotMgr;

    QMenu *renderLineThicknessMenu_;
    QMenu *_jobMenu;

    QAction *_hardwareStereoAction;
    QAction *_stereoToggleAction;

    LocalSocketScript* _scriptSocket;
};

//some shortcuts for MIMainWindow::instance()->Log, etc
void MIMainWindowLog(const std::string&);
void MIMainWindowDebug(const std::string&);
void MIMainWindowLeftFooter(const std::string&, int timeout = 0);
void MIMainWindowMiddleFooter(const std::string&);
void MIMainWindowRightFooter(const std::string&);

#endif // ifndef MAINWINDOW_H
