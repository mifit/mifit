#ifndef MIFIT_UI_MYAPP_H_
#define MIFIT_UI_MYAPP_H_

/**
 * The top level class from which the application is started.
 */

#include <map/maplib.h>
#include <QtCore/QDir>
#include <QtGui/QApplication>
#include "core/corelib.h"

class MISingleInstanceChecker;
class GeomRefiner;
class QSettings;

class Application : public QApplication
{
    Q_OBJECT
  
    // eventually, we should revoke this friendship, but it's grandfathered in now
    friend class EnvironmentPreferences;
    friend class GeneralPreferences;
    friend class Tools;
    friend class CMapHeader;

    friend class MIGLWidget;
    friend class MIMainWindow;

    GeomRefiner *geomrefiner;

    bool dimNonactiveModels;
    bool hardwareStereoAvailable;
    bool incrementallyColorModels;
    bool onCloseSaveActiveModelToPdb;
    bool silent_mode;
    bool xfitMouseMode;

    std::string CrystalData;
    std::string HTMLBrowser;
    std::string MIbin;
    std::string MolimageHome;
    std::string PentDir;
    std::string ShelxHome;
    std::string SmilesDbCommand;
    std::string XFitDict;
    std::string XFitDictSetting;
    std::string appname;
    std::string checkpointDirectory;
    std::string userSettingsDir;

    QDir jobLogsDir;


    MIPalette *lpal;
    void BuildPalette();
    PaletteColor BackgroundColor;
    int GammaCorrection;
    int concurrentJobLimit;

    chemlib::Residue *ResidueBuffer;

    void WriteProfiles(); // Write the profiles to the registry (Windows) or the ~/.MIFitrc file.
    void WritePalette(); // Write the profiles to the registry (Windows) or the ~/.MIFitrc file.
    void WriteSecStr();  // Write the secondary structure settings to the registry (Windows) or the ~/.MIFitrc file.
    void WriteBValues(); // Write the b-value settings to the registry (Windows) or the ~/.MIFitrc file.
    void WriteAtomTypes(); // Write the atom type coloring info to the registry (Windows) or the ~/.MIFitrc file.

    void MixBackgroundColor(); // Mixes the background color into the palette for proper depthcueing.

    // Gamma is currently defaulted to 1 and not changed, so all this is kind of useless
    int ApplyGammaCorrection(const int color);
    void SetGammaCorrection(const double gamma)
    {
        GammaCorrection = ROUND(gamma*1000.0);
    }
    double GetGammaCorrection()
    {
        return (double)GammaCorrection/1000.0;
    }

    void LoadDictionary();   // reload the current dictionary as specified by  XFitDictSetting
    void saveDict();

public:
    static Application *instance();

    Application(int &argc, char **argv);
    ~Application();
    void Init();
public slots:
    void AfterInit();
public:

    void Write();
    void SetEnv();

    std::string getDictionary();
    bool SetDictionary(const std::string &dict);

    bool GetSilentMode()
    {
        return silent_mode;
    }
    void SetSilentMode(bool t)
    {
        silent_mode = t;
    }

    const std::string&GetBinDirectory();

    const std::string&GetMolimageHome() const
    {
        return MolimageHome;
    }
    void SetMolimageHome(bool reset = false);

    const std::string&GetShelxHome();
    void SetShelxHome(bool reset = false);

    const std::string&GetCrystalData()
    {
        return CrystalData;
    }
    std::string GetCrystalData(const char *crystal, const char *key);
    void SetCrystalData(bool reset = false);

    QString latestFileBrowseDirectory(const QString &path);
    QDir jobLogsDirectory() const
    {
        return jobLogsDir;
    }

    static QString getOpenFileName(QWidget *parent, const QString &caption, const QString &filter = QString());
    static QString getExistingDirectory(QWidget *parent, const QString &caption, const QString &filter);


    bool isHardwareStereoAvailable()
    {
        return hardwareStereoAvailable;
    }
    void setHardwareStereoAvailable(bool t)
    {
        hardwareStereoAvailable = t;
    }

    bool GetCrystalCell(const char *crystal, float &a, float &b, float &c, float &alpha, float &beta, float &gamma);
    int GetCrystalNCSSymmops(const char *crystal, float ncrsymm[MISymmop::MAXNCRSYMMOPS][12]);

    MIPalette *GetLPpal()
    {
        return lpal;
    }                                   // Returns a pointer to the palette.

    void backgroundColor();
    PaletteColor GetBackgroundColor()
    {
        return BackgroundColor;
    }                                                           //Returns the background color.  The colors are defined in colors.h.
    void SetBackgroundColor(PaletteColor c)
    {
        BackgroundColor = c;
    }                                                              // Change the background color.  The colors are defined in colors.h.

    bool GetXfitMouseMode()
    {
        return xfitMouseMode;
    }                                               // Returns the mouse mode - true if it is Xfit emulation mode.
    void SetXfitMouseMode(bool m)
    {
        xfitMouseMode = m;
    }                                                  // Sets the Xfit mouse emulation mode

    bool LabelPicks; // if true then picked atoms are labeled.
    bool LabelToggle; // If true then picking an atom for the second time toggles the label off.

    void ClearResidueBuffer(); // Frees and clears the residue buffer.
    chemlib::Residue *GetResidueBuffer()
    {
        return ResidueBuffer;
    }                                                            // Returns a pointer to the residue buffer.
    bool CopyResidueBuffer(const chemlib::Residue *buffer); // Copy one residue to the to the end of the ResidueBuffer list.

    GeomRefiner *GetGeomRefiner()
    {
        return geomrefiner;
    }                                                    //The geometry refiner is a helper class for refining model geometry and real-space refinement

    void toggleStereo();
    void toggleHardwareStereo();

    long GetPid();
};

GeomRefiner *MIFitGeomRefiner();
chemlib::MIMolDictionary *MIFitDictionary();

class QProgressDialog;

class MIBusyManager
{
public:
    static MIBusyManager *instance();

    bool Busy() const ;
    void SetBusy(bool busy);

    void SetLabel(const char*);

    void ForceAbort();
    bool CheckForAbort();

    void SetWaitCursor();
    void StartWaitCursor();
    void StopWaitCursor();

private:
    MIBusyManager();
    void SetCursor(int);

    unsigned int m_busy;
    bool MustAbortNow;
    unsigned int m_cursor_number;
    QProgressDialog *PROGRESS;
    unsigned int PROG_VALUE;
    int PROG_DIRECTION;

    static MIBusyManager *_instance;
};


#endif /*MIFIT_UI_MYAPP_H_*/
