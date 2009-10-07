#include <sys/types.h>
#include <sys/stat.h>
#include <errno.h>


#include <nongui/nonguilib.h>

#include <chemlib/chemlib.h>
#include <chemlib/RESIDUE_.h>
#include <conflib/conflib.h>
#include <util/utillib.h>
#include "core/corelib.h"

#include "Application.h"
#include "tools.h"
#include "DictEditCanvas.h"
#include "GLFormatEdit.h"

#include "ui/MIDialog.h"
#include <util/FileIo.h>
#include "MIMainWindow.h"
#include "molw.h" // for cursors

#include <QApplication>
#include <QDir>
#include <QFileDialog>
#include <QFileInfo>
#include <QMessageBox>
#include <QProcess>
#include <QProgressDialog>
#include <QSettings>

#ifndef _WIN32
#include <unistd.h>
#include <pwd.h>
#include <sys/param.h>
#include <sys/types.h>
#include <sys/stat.h>
#else
#include <process.h>
#define strncasecmp strnicmp
#endif

#include "core/Version.h"

int PaletteChanged = true;
#define GLOBAL


using namespace std;
using namespace chemlib;

/*
 * tilde expand a file name
 *
 *   returns: the expanded name of the file (in a static buffer)
 *            or NIL(char)
 */
const char* TildeExpand(const char* filename) {
#ifdef _WIN32
  /* tilde's don't make sense in windows */
  return filename;
#else
  static char dummy[1024];
  char username[20];
  const char* loc;
  struct passwd* pw;

  if ((filename == NULL) || !strcmp(filename, "")) {
    return ("");
  }

  if (filename[0] != '~') {
    (void) strcpy(dummy, filename);
    return (dummy);
  }

  /* tilde at the beginning now */
  if ((filename[1] == '\0') || (filename[1] == '/') || (filename[1] == ' ')) {
    /* current user */
    char* home;

    if ((home = getenv("HOME")) == (char*) 0) {
      /* fall back on /etc/passwd */
      if ((pw = getpwuid(getuid())) == (struct passwd*) NULL) {
        return ("");
      }
      strcpy(dummy, pw->pw_dir);
    } else {
      strcpy(dummy, home);
    }
    strcat(dummy, &filename[1]);

  } else {
    loc = strchr(filename, '/');
    if (loc == 0) {
      strcpy(username, &filename[1]);
    } else {
      (void) strncpy(username, &filename[1], loc - &filename[1]);
      username[loc - &filename[1]] = '\0';
    }
    if ((pw = getpwnam(username)) == NULL) {
      return ("");
    }
    strcpy(dummy, pw->pw_dir);
    if (loc != 0) {
      strcat(dummy, loc);
    }
  }
  return (dummy);
#endif
}


// Find the absolute path where this application has been run from.
std::string findMIFitBin() {
  QFileInfo fi(QCoreApplication::applicationDirPath());
  return fi.canonicalFilePath().toStdString();
}

class MIConfigImpl : public MIConfig {
  public:
    MIConfigImpl(const std::string& name, bool read_enabled=true) {
      _instance = this;
      _settings = new QSettings("MIFit", name.c_str());
      QSettings oldSettings("Rigaku", name.c_str());
      foreach (QString key, oldSettings.allKeys())
          _settings->setValue(key, oldSettings.value(key));
      _read_enabled=read_enabled;
    }

    QSettings *GetQSettings() { return _settings; }

    void SetReadEnabled(bool state) {
      _read_enabled=state;
    }


    ~MIConfigImpl() {
      delete _settings;
    }

    bool Read(const std::string& key, std::string &value, const std::string& defaultValue, bool assignDefaultValue = true) {
      if (_read_enabled && _settings->contains(key.c_str())) {
        QByteArray ba=_settings->value(key.c_str(), QByteArray(defaultValue.c_str(),defaultValue.size())).toByteArray();
        value=std::string(ba.constData(),ba.size());
        return true;
      }

      if (assignDefaultValue) {
        value=defaultValue;
      }
      return false;
    }

    bool Read(const std::string& key, long* value, long defaultValue, bool assignDefaultValue = true) {
      if (_read_enabled && _settings->contains(key.c_str())) {
        qlonglong l=(_settings->value(key.c_str(), (qlonglong)defaultValue).toLongLong());
        *value=(long)l;
        return true;
      }
      if (assignDefaultValue) {
        *value = defaultValue;
      }
      return false;
    }

    bool Read(const std::string& key, bool* value, bool defaultValue, bool assignDefaultValue = true) {
      if (_read_enabled && _settings->contains(key.c_str())) {
        *value=_settings->value(key.c_str(), defaultValue).toBool();
        return true;
      }
      if (assignDefaultValue) {
        *value = defaultValue;
      }
      return false;
    }

    bool Read(const std::string& key, double* value, double defaultValue, bool assignDefaultValue = true) {
      if (_read_enabled && _settings->contains(key.c_str())) {
        *value=_settings->value(key.c_str(), defaultValue).toDouble();
        return true;
      }
      if (assignDefaultValue) {
        *value = defaultValue;
      }
      return false;
    }

    bool Read(const std::string& key, std::string& value) {
      return MIConfig::Read(key, value);
    }

    bool Read(const std::string& key, long* value) {
      return MIConfig::Read(key, value);
    }

    bool Read(const std::string& key, bool* value) {
      return MIConfig::Read(key, value);
    }

    bool Read(const std::string& key, double* value) {
      return MIConfig::Read(key, value);
    }

    bool Write(const std::string& key, const std::string& value) {
      _settings->setValue(key.c_str(), QByteArray(value.c_str(),value.size()));
      return true;
    }

    bool Write(const std::string& key, long value) {
      _settings->setValue(key.c_str(), (qlonglong)value);
      return true;
    }

    bool Write(const std::string& key, bool value) {
      _settings->setValue(key.c_str(), value);
      return true;
    }

    bool Write(const std::string& key, double value) {
      _settings->setValue(key.c_str(), value);
      return true;
    }

    bool HasKey(const std::string &key) {
      return _settings->contains(key.c_str());
    }


    void Flush() {
      _settings->sync();
    }


  private:
    QSettings* _settings;
    bool _read_enabled;
};

Application* Application::instance_;

Application *Application::instance() {
  if (!instance_)
    instance_=new Application();
  return instance_;
}

namespace {
    const QString GL_FORMAT_GROUP("glformat");
    const QString GL_FORMAT_DEFAULT("default");
}

Application::Application(void) {

  instance_ = this;

  FileIo::setAsDefaultIo();

  static const char* name = "MIFit";
  config = new MIConfigImpl(name);

  appname = std::string(name);
  MIbin = "";

  silent_mode = false;
  hardwareStereoAvailable = false;
  XFitDict = "";

  CrystalData = "";
  MolimageHome = "";

  XFitDictSetting = "";
  PentDir = "";
  HTMLBrowser = "";
  SmilesDbCommand = "";
  checkpointDirectory = "";
  ShelxHome = "";

  jobLogsDir = QDir::home().absoluteFilePath(".MIFit_jobLogs");
  if (!jobLogsDir.exists()) {
      jobLogsDir.mkpath(".");
  }

  QSettings* settings = MIGetQSettings();
  if (settings->childGroups().contains(GL_FORMAT_GROUP)) {
    settings->beginGroup(GL_FORMAT_GROUP);
    if (settings->childGroups().contains(GL_FORMAT_DEFAULT)) {
      settings->beginGroup(GL_FORMAT_DEFAULT);
      QGLFormat glformat = GLFormatEdit::readSettings(*settings);
      QGLFormat::setDefaultFormat(glformat);
      settings->endGroup();
    }
    settings->endGroup();
  }

  Init();
}

Application::~Application() {
  config->WriteProfileInt("Options", "MouseMode", xfitMouseMode);
  config->WriteProfileInt("Options", "IncrementallyColorModels", incrementallyColorModels);
  config->WriteProfileInt("Options", "DimNonactiveModels", dimNonactiveModels);
  config->WriteProfileInt("Options", "ConcurrentJobLimit", concurrentJobLimit);

  Write();

  QSettings* settings = MIGetQSettings();

  settings->beginGroup(GL_FORMAT_GROUP);
  settings->beginGroup(GL_FORMAT_DEFAULT);
  QGLFormat glformat = QGLFormat::defaultFormat();
  GLFormatEdit::writeSettings(*settings, glformat);
  settings->endGroup();
  settings->endGroup();

  delete config;
  delete lpal;
  ClearResidueBuffer();

  if (geomrefiner) { //Cleanup the GeomRefiner
    if (geomrefiner->dict.GetModified()) {
        if (QMessageBox::question(MIMainWindow::instance(), "Dictionary Modified",
                                  "The Dictionary has been modified\nDo you want to save it?",
                                  QMessageBox::Yes | QMessageBox::No, QMessageBox::Yes) == QMessageBox::Yes) {
        saveDict();
      }
    }
  }
  delete geomrefiner;


}


int Application::ApplyGammaCorrection(const int color) {
  double gamma = GetGammaCorrection();
  if (gamma <= 0.0) {
    gamma = 0.2;
  }
  double fcolor;
  gamma = 1.0/gamma;
  fcolor = (double)color/255.0;
  fcolor = pow(fcolor, gamma);
  fcolor *= 255.0;
  if (fcolor > 255.0) {
    return 255;
  }
  if (fcolor < 0) {
    return 0;
  }
  return ROUND(fcolor);
}

void Application::BuildPalette() {
  lpal = new MIPalette;
  lpal->colors.resize(Colors_NUMBERPALETTE);
  lpal->colors[Colors::PBLACK].red = config->GetProfileInt("Palette", "Black.red", 0);
  lpal->colors[Colors::PBLACK].green = config->GetProfileInt("Palette", "Black.green", 0);
  lpal->colors[Colors::PBLACK].blue = config->GetProfileInt("Palette", "Black.blue", 0);
  lpal->colors[Colors::PRED].red = config->GetProfileInt("Palette", "Red.red", 255);
  lpal->colors[Colors::PRED].green = config->GetProfileInt("Palette", "Red.green", 20);
  lpal->colors[Colors::PRED].blue = config->GetProfileInt("Palette", "Red.blue", 20);
  lpal->colors[Colors::PBLUE].red = config->GetProfileInt("Palette", "Blue.red", 20);
  lpal->colors[Colors::PBLUE].green = config->GetProfileInt("Palette", "Blue.green", 20);
  lpal->colors[Colors::PBLUE].blue = config->GetProfileInt("Palette", "Blue.blue", 255);
  lpal->colors[Colors::PGREEN].red = config->GetProfileInt("Palette", "Green.red", 20);
  lpal->colors[Colors::PGREEN].green = config->GetProfileInt("Palette", "Green.green", 255);
  lpal->colors[Colors::PGREEN].blue = config->GetProfileInt("Palette", "Green.blue", 20);
  lpal->colors[Colors::PCYAN].red = config->GetProfileInt("Palette", "Cyan.red", 20);
  lpal->colors[Colors::PCYAN].green = config->GetProfileInt("Palette", "Cyan.green", 255);
  lpal->colors[Colors::PCYAN].blue = config->GetProfileInt("Palette", "Cyan.blue", 255);
  lpal->colors[Colors::PMAGENTA].red = config->GetProfileInt("Palette", "Magenta.red", 255);
  lpal->colors[Colors::PMAGENTA].green = config->GetProfileInt("Palette", "Magenta.green", 20);
  lpal->colors[Colors::PMAGENTA].blue = config->GetProfileInt("Palette", "Magenta.blue", 255);
  lpal->colors[Colors::PYELLOW].red = config->GetProfileInt("Palette", "Yellow.red", 255);
  lpal->colors[Colors::PYELLOW].green = config->GetProfileInt("Palette", "Yellow.green", 255);
  lpal->colors[Colors::PYELLOW].blue = config->GetProfileInt("Palette", "Yellow.blue", 20);
  lpal->colors[Colors::PWHITE].red = config->GetProfileInt("Palette", "White.red", 253);
  lpal->colors[Colors::PWHITE].green = config->GetProfileInt("Palette", "White.green", 253);
  lpal->colors[Colors::PWHITE].blue = config->GetProfileInt("Palette", "White.blue", 253);
  lpal->colors[Colors::PPINK].red = config->GetProfileInt("Palette", "Pink.red", 255);
  lpal->colors[Colors::PPINK].green = config->GetProfileInt("Palette", "Pink.green", 128);
  lpal->colors[Colors::PPINK].blue = config->GetProfileInt("Palette", "Pink.blue", 128);
  lpal->colors[Colors::PORANGE].red = config->GetProfileInt("Palette", "Orange.red", 255);
  lpal->colors[Colors::PORANGE].green = config->GetProfileInt("Palette", "Orange.green", 128);
  lpal->colors[Colors::PORANGE].blue = config->GetProfileInt("Palette", "Orange.blue", 0);
  lpal->colors[Colors::PBROWN].red = config->GetProfileInt("Palette", "Brown.red", 0);
  lpal->colors[Colors::PBROWN].green = config->GetProfileInt("Palette", "Brown.green", 0);
  lpal->colors[Colors::PBROWN].blue = config->GetProfileInt("Palette", "Brown.blue", 0);

  lpal->colors[Colors::PCUSTOM1].red = config->GetProfileInt("User Colors", "User1.red", 255);
  lpal->colors[Colors::PCUSTOM1].green = config->GetProfileInt("User Colors", "User1.green", 72);
  lpal->colors[Colors::PCUSTOM1].blue = config->GetProfileInt("User Colors", "User1.blue", 96);
  lpal->colors[Colors::PCUSTOM2].red = config->GetProfileInt("User Colors", "User2.red", 90);
  lpal->colors[Colors::PCUSTOM2].green = config->GetProfileInt("User Colors", "User2.green", 255);
  lpal->colors[Colors::PCUSTOM2].blue = config->GetProfileInt("User Colors", "User2.blue", 90);
  lpal->colors[Colors::PCUSTOM3].red = config->GetProfileInt("User Colors", "User3.red", 102);
  lpal->colors[Colors::PCUSTOM3].green = config->GetProfileInt("User Colors", "User3.green", 164);
  lpal->colors[Colors::PCUSTOM3].blue = config->GetProfileInt("User Colors", "User3.blue", 255);
  lpal->colors[Colors::PCUSTOM4].red = config->GetProfileInt("User Colors", "User4.red", 162);
  lpal->colors[Colors::PCUSTOM4].green = config->GetProfileInt("User Colors", "User4.green", 115);
  lpal->colors[Colors::PCUSTOM4].blue = config->GetProfileInt("User Colors", "User4.blue", 20);
  lpal->colors[Colors::PCUSTOM5].red = config->GetProfileInt("User Colors", "User5.red", 255);
  lpal->colors[Colors::PCUSTOM5].green = config->GetProfileInt("User Colors", "User5.green", 143);
  lpal->colors[Colors::PCUSTOM5].blue = config->GetProfileInt("User Colors", "User5.blue", 32);
  lpal->colors[Colors::PCUSTOM6].red = config->GetProfileInt("User Colors", "User6.red", 255);
  lpal->colors[Colors::PCUSTOM6].green = config->GetProfileInt("User Colors", "User6.green", 143);
  lpal->colors[Colors::PCUSTOM6].blue = config->GetProfileInt("User Colors", "User6.blue", 190);
  lpal->colors[Colors::PCUSTOM7].red = config->GetProfileInt("User Colors", "User7.red", 150);
  lpal->colors[Colors::PCUSTOM7].green = config->GetProfileInt("User Colors", "User7.green", 0);
  lpal->colors[Colors::PCUSTOM7].blue = config->GetProfileInt("User Colors", "User7.blue", 0);
  lpal->colors[Colors::PCUSTOM8].red = config->GetProfileInt("User Colors", "User8.red", 0);
  lpal->colors[Colors::PCUSTOM8].green = config->GetProfileInt("User Colors", "User8.green", 150);
  lpal->colors[Colors::PCUSTOM8].blue = config->GetProfileInt("User Colors", "User8.blue", 0);
  lpal->colors[Colors::PCUSTOM9].red = config->GetProfileInt("User Colors", "User9.red", 0);
  lpal->colors[Colors::PCUSTOM9].green = config->GetProfileInt("User Colors", "User9.green", 0);
  lpal->colors[Colors::PCUSTOM9].blue = config->GetProfileInt("User Colors", "User9.blue", 150);
  lpal->colors[Colors::PCUSTOM10].red = config->GetProfileInt("User Colors", "User10.red", 150);
  lpal->colors[Colors::PCUSTOM10].green = config->GetProfileInt("User Colors", "User10.green", 150);

  lpal->colors[Colors::PCUSTOM10].blue = config->GetProfileInt("User Colors", "User10.blue", 0);

  lpal->colors[Colors::PMAP1].red = config->GetProfileInt("Map Colors", "Map1.red", 0);
  lpal->colors[Colors::PMAP1].green = config->GetProfileInt("Map Colors", "Map1.green", 0);
  lpal->colors[Colors::PMAP1].blue = config->GetProfileInt("Map Colors", "Map1.blue", 230);
  lpal->colors[Colors::PMAP2].red = config->GetProfileInt("Map Colors", "Map2.red", 128);
  lpal->colors[Colors::PMAP2].green = config->GetProfileInt("Map Colors", "Map2.green", 0);
  lpal->colors[Colors::PMAP2].blue = config->GetProfileInt("Map Colors", "Map2.blue", 230);
  lpal->colors[Colors::PMAP3].red = config->GetProfileInt("Map Colors", "Map3.red", 230);
  lpal->colors[Colors::PMAP3].green = config->GetProfileInt("Map Colors", "Map3.green", 0);
  lpal->colors[Colors::PMAP3].blue = config->GetProfileInt("Map Colors", "Map3.blue", 230);
  lpal->colors[Colors::PMAP4].red = config->GetProfileInt("Map Colors", "Map4.red", 230);
  lpal->colors[Colors::PMAP4].green = config->GetProfileInt("Map Colors", "Map4.green", 0);
  lpal->colors[Colors::PMAP4].blue = config->GetProfileInt("Map Colors", "Map4.blue", 128);
  lpal->colors[Colors::PMAP5].red = config->GetProfileInt("Map Colors", "Map5.red", 230);
  lpal->colors[Colors::PMAP5].green = config->GetProfileInt("Map Colors", "Map5.green", 0);
  lpal->colors[Colors::PMAP5].blue = config->GetProfileInt("Map Colors", "Map5.blue", 0);
  lpal->colors[Colors::PMAP6].red = config->GetProfileInt("Map Colors", "Map6.red", 0);
  lpal->colors[Colors::PMAP6].green = config->GetProfileInt("Map Colors", "Map6.green", 230);
  lpal->colors[Colors::PMAP6].blue = config->GetProfileInt("Map Colors", "Map6.blue", 0);
  lpal->colors[Colors::PMAP7].red = config->GetProfileInt("Map Colors", "Map7.red", 0);
  lpal->colors[Colors::PMAP7].green = config->GetProfileInt("Map Colors", "Map7.green", 230);
  lpal->colors[Colors::PMAP7].blue = config->GetProfileInt("Map Colors", "Map7.blue", 128);
  lpal->colors[Colors::PMAP8].red = config->GetProfileInt("Map Colors", "Map8.red", 0);
  lpal->colors[Colors::PMAP8].green = config->GetProfileInt("Map Colors", "Map8.green", 230);
  lpal->colors[Colors::PMAP8].blue = config->GetProfileInt("Map Colors", "Map8.blue", 230);
  lpal->colors[Colors::PMAP9].red = config->GetProfileInt("Map Colors", "Map9.red", 0);
  lpal->colors[Colors::PMAP9].green = config->GetProfileInt("Map Colors", "Map9.green", 128);
  lpal->colors[Colors::PMAP9].blue = config->GetProfileInt("Map Colors", "Map9.blue", 230);
  lpal->colors[Colors::PMAP10].red = config->GetProfileInt("Map Colors", "Map10.red", 128);
  lpal->colors[Colors::PMAP10].green = config->GetProfileInt("Map Colors", "Map10.green", 128);
  lpal->colors[Colors::PMAP10].blue = config->GetProfileInt("Map Colors", "Map10.blue", 230);

  lpal->colors[Colors::PCONTOUR1].red = config->GetProfileInt("Contour Colors", "Contour1.red", 218);
  lpal->colors[Colors::PCONTOUR1].green = config->GetProfileInt("Contour Colors", "Contour1.green", 152);
  lpal->colors[Colors::PCONTOUR1].blue = config->GetProfileInt("Contour Colors", "Contour1.blue", 207);
  lpal->colors[Colors::PCONTOUR2].red = config->GetProfileInt("Contour Colors", "Contour2.red", 113);
  lpal->colors[Colors::PCONTOUR2].green = config->GetProfileInt("Contour Colors", "Contour2.green", 87);
  lpal->colors[Colors::PCONTOUR2].blue = config->GetProfileInt("Contour Colors", "Contour2.blue", 185);
  lpal->colors[Colors::PCONTOUR3].red = config->GetProfileInt("Contour Colors", "Contour3.red", 72);
  lpal->colors[Colors::PCONTOUR3].green = config->GetProfileInt("Contour Colors", "Contour3.green", 107);
  lpal->colors[Colors::PCONTOUR3].blue = config->GetProfileInt("Contour Colors", "Contour3.blue", 254);
  lpal->colors[Colors::PCONTOUR4].red = config->GetProfileInt("Contour Colors", "Contour4.red", 75);
  lpal->colors[Colors::PCONTOUR4].green = config->GetProfileInt("Contour Colors", "Contour4.green", 230);
  lpal->colors[Colors::PCONTOUR4].blue = config->GetProfileInt("Contour Colors", "Contour4.blue", 251);
  lpal->colors[Colors::PCONTOUR5].red = config->GetProfileInt("Contour Colors", "Contour5.red", 0);
  lpal->colors[Colors::PCONTOUR5].green = config->GetProfileInt("Contour Colors", "Contour5.green", 255);
  lpal->colors[Colors::PCONTOUR5].blue = config->GetProfileInt("Contour Colors", "Contour5.blue", 0);
  lpal->colors[Colors::PCONTOUR6].red = config->GetProfileInt("Contour Colors", "Contour6.red", 255);
  lpal->colors[Colors::PCONTOUR6].green = config->GetProfileInt("Contour Colors", "Contour6.green", 255);
  lpal->colors[Colors::PCONTOUR6].blue = config->GetProfileInt("Contour Colors", "Contour6.blue", 30);
  lpal->colors[Colors::PCONTOUR7].red = config->GetProfileInt("Contour Colors", "Contour7.red", 255);
  lpal->colors[Colors::PCONTOUR7].green = config->GetProfileInt("Contour Colors", "Contour7.green", 200);
  lpal->colors[Colors::PCONTOUR7].blue = config->GetProfileInt("Contour Colors", "Contour7.blue", 30);
  lpal->colors[Colors::PCONTOUR8].red = config->GetProfileInt("Contour Colors", "Contour8.red", 245);
  lpal->colors[Colors::PCONTOUR8].green = config->GetProfileInt("Contour Colors", "Contour8.green", 141);
  lpal->colors[Colors::PCONTOUR8].blue = config->GetProfileInt("Contour Colors", "Contour8.blue", 3);
  lpal->colors[Colors::PCONTOUR9].red = config->GetProfileInt("Contour Colors", "Contour9.red", 255);
  lpal->colors[Colors::PCONTOUR9].green = config->GetProfileInt("Contour Colors", "Contour9.green", 60);
  lpal->colors[Colors::PCONTOUR9].blue = config->GetProfileInt("Contour Colors", "Contour9.blue", 00);
  lpal->colors[Colors::PCONTOUR10].red = config->GetProfileInt("Contour Colors", "Contour10.red", 220);
  lpal->colors[Colors::PCONTOUR10].green = config->GetProfileInt("Contour Colors", "Contour10.green", 0);
  lpal->colors[Colors::PCONTOUR10].blue = config->GetProfileInt("Contour Colors", "Contour10.blue", 0);

  BackgroundColor=PaletteColor(config->GetProfileInt("View Parameters", "BackgroundColorRed", 0), config->GetProfileInt("View Parameters", "BackgroundColorGreen", 0), config->GetProfileInt("View Parameters", "BackgroundColorBlue", 0));
  MixBackgroundColor();

  Colors::sec_colors[Colors::HELIX] = config->GetProfileInt("SecStr Colors", "Helix", Colors::MAGENTA);
  Colors::sec_colors[Colors::SHEET] = config->GetProfileInt("SecStr Colors", "Sheet", Colors::YELLOW);
  Colors::sec_colors[Colors::COIL] = config->GetProfileInt("SecStr Colors", "Coil", Colors::BLUE);

  Colors::BValueColors[0] = config->GetProfileInt("BValue Colors", "Color1", Colors::CONTOUR1);
  Colors::BValueRanges[0] = config->GetProfileInt("BValue Colors", "Range1", 300);
  Colors::BValueColors[1] = config->GetProfileInt("BValue Colors", "Color2", Colors::CONTOUR2);
  Colors::BValueRanges[1] = config->GetProfileInt("BValue Colors", "Range2", 600);
  Colors::BValueColors[2] = config->GetProfileInt("BValue Colors", "Color3", Colors::CONTOUR3);
  Colors::BValueRanges[2] = config->GetProfileInt("BValue Colors", "Range3", 1000);
  Colors::BValueColors[3] = config->GetProfileInt("BValue Colors", "Color4", Colors::CONTOUR4);
  Colors::BValueRanges[3] = config->GetProfileInt("BValue Colors", "Range4", 1500);
  Colors::BValueColors[4] = config->GetProfileInt("BValue Colors", "Color5", Colors::CONTOUR5);
  Colors::BValueRanges[4] = config->GetProfileInt("BValue Colors", "Range5", 2000);
  Colors::BValueColors[5] = config->GetProfileInt("BValue Colors", "Color6", Colors::CONTOUR6);
  Colors::BValueRanges[5] = config->GetProfileInt("BValue Colors", "Range6", 2500);
  Colors::BValueColors[6] = config->GetProfileInt("BValue Colors", "Color7", Colors::CONTOUR7);
  Colors::BValueRanges[6] = config->GetProfileInt("BValue Colors", "Range7", 3500);
  Colors::BValueColors[7] = config->GetProfileInt("BValue Colors", "Color8", Colors::CONTOUR8);
  Colors::BValueRanges[7] = config->GetProfileInt("BValue Colors", "Range8", 5000);
  Colors::BValueColors[8] = config->GetProfileInt("BValue Colors", "Color9", Colors::CONTOUR9);
  Colors::BValueRanges[8] = config->GetProfileInt("BValue Colors", "Range9", 7500);
  Colors::BValueColors[9] = config->GetProfileInt("BValue Colors", "Color10", Colors::CONTOUR10);
  Colors::BValueRanges[9] = config->GetProfileInt("BValue Colors", "Range10", 10000);

  /* temporary code used to write out the colros to a file for the manual
     8 left in because it seemed it might be useful at a future date - dem 8/9/2005
     FILE * fp = fopen("color_names.txt","w");
     for(int i=0;i<NCOLORS; i++){
     fprintf(fp, "%d\t%4d%4d%4d\t%s\n", i,
     (int)lpal->colors[PaletteIndex(i)].red,
     (int)lpal->colors[PaletteIndex(i)].green,
     (int)lpal->colors[PaletteIndex(i)].blue,
     colornames[i]);
     }
   */
}

void Application::MixBackgroundColor() {
  int i, ir, ib, r, g, b, d, rcolor, gcolor, bcolor;
  if (!lpal) {
    return;
  }
  rcolor = GetBackgroundColor().red;
  gcolor = GetBackgroundColor().green;
  bcolor = GetBackgroundColor().blue;
  for (ir = 0; ir < Colors_NUMBERCOLORS; ir++) {
    Colors::RPallette[PaletteIndex(ir)] = lpal->colors[PaletteIndex(ir)].red;
    Colors::GPallette[PaletteIndex(ir)] = lpal->colors[PaletteIndex(ir)].green;
    Colors::BPallette[PaletteIndex(ir)] = lpal->colors[PaletteIndex(ir)].blue;
    for (d = 1; d < 10; d++) {
      ib = ir*10+ d;
      r = lpal->colors[PaletteIndex(ir)].red*(255-d*20)/255+ rcolor * d*20/ 255;
      g = lpal->colors[PaletteIndex(ir)].green*(255-d*20)/255+ gcolor * d*20/ 255;
      b = lpal->colors[PaletteIndex(ir)].blue*(255-d*20)/255+ bcolor * d*20/ 255;
      lpal->colors[ib].red = (unsigned char)r;
      lpal->colors[ib].green = (unsigned char)g;
      lpal->colors[ib].blue = (unsigned char)b;
      Colors::RPallette[ib] = r;
      Colors::GPallette[ib] = g;
      Colors::BPallette[ib] = b;
    }
  }
  for (i = 0; i < Colors_NUMBERPALETTE; i++) {
    Colors::RPallette[i] = ApplyGammaCorrection(Colors::RPallette[i]);
    Colors::GPallette[i] = ApplyGammaCorrection(Colors::GPallette[i]);
    Colors::BPallette[i] = ApplyGammaCorrection(Colors::BPallette[i]);
  }
}

void Application::backgroundColor() {
  PaletteColor color = MIGetColorFromUser(0, BackgroundColor);
  SetBackgroundColor(color);
}

bool Application::CopyResidueBuffer(const RESIDUE* buffer) {
  RESIDUE* newres = new RESIDUE(*buffer);
  if (ResidueBuffer) {
    RESIDUE* end = ResidueBuffer;
    while (end->next() != NULL) {
      end = end->next();
    }
    end->insertResidue(newres);
  } else {
    ResidueBuffer = newres;
  }
  return true;
}

void Application::ClearResidueBuffer() {
  if (ResidueBuffer) {
    FreeResidueList(ResidueBuffer);
  }
  ResidueBuffer = NULL;
}




class myColorSetter : public MIColorSetter {
public:
  bool operator()(MIAtom* atom) const {
    atom->setColor(color_by_name(atom->name()));
    if (atom->color() == Colors::BLACK) {
      atom->setColor(Colors::WHITE);
    }
    return true;
  }

  bool operator()(MIAtom* atom, char c) const {
    atom->setColor(color_by_name(atom->name(), c));
    if (atom->color() == Colors::BLACK) {
      atom->setColor(Colors::WHITE);
    }
    return true;
  }

};

class myTorsionWritePrompt : public MITorsionWritePrompt {
public:
  bool operator()() {
    return (MIMessageBox("Write torsions to file?", "Write Torsions?",
              MIDIALOG_YES_NO | MIDIALOG_NO_DEFAULT) == MI_YES);
  }

};




class myMolPrefsHandler : public MolPrefsHandler {
public:
  void operator()(bool* breakByDiscontinuity, bool* breakByNonpeptide) {
    MIConfig* config = MIConfig::Instance();
    config->Read("Options/breakByDiscontinuity", breakByDiscontinuity, true);
    config->Read("Options/breakByNonpeptide", breakByNonpeptide, false);
  }

};

// static int StereoAttribs[] = {
//   WX_GL_STEREO,
//   WX_GL_RGBA,
//   WX_GL_MIN_RED, 1,
//   WX_GL_MIN_GREEN, 1,
//   WX_GL_MIN_BLUE, 1,
//   WX_GL_MIN_ALPHA, 1,
//   WX_GL_DEPTH_SIZE,  1,
//   WX_GL_STENCIL_SIZE, 1,
//   WX_GL_DOUBLEBUFFER,
//   0
// };

// static int StereoAttribsLowSpec[] = {
//   WX_GL_STEREO,
//   WX_GL_RGBA,
//   WX_GL_DOUBLEBUFFER,
//   0
// };

// static int MonoAttribs[] = {
//   WX_GL_RGBA,
//   WX_GL_MIN_RED, 1,
//   WX_GL_MIN_GREEN, 1,
//   WX_GL_MIN_BLUE, 1,
//   WX_GL_MIN_ALPHA, 1,
//   WX_GL_DEPTH_SIZE,  1,
//   WX_GL_STENCIL_SIZE, 1,
//   WX_GL_DOUBLEBUFFER,
//   0
// };

// static int MonoAttribsLowSpec[] = {
//   WX_GL_RGBA,
//   WX_GL_DOUBLEBUFFER,
//   0
// };


void Application::Write() {
  if (!CrystalData.empty()) {
    config->Write("CRYSTALDATA", CrystalData);
  }
  config->Write("XFITDICT", XFitDictSetting);
  if (!HTMLBrowser.empty()) {
    config->Write("HTMLBROWSER", HTMLBrowser);
  }
  if (!ShelxHome.empty()) {
    config->Write("SHELXHOME", ShelxHome);
  }
  if (!SmilesDbCommand.empty()) {
    config->Write("SmilesDbCommand", SmilesDbCommand);
  }
  if (!checkpointDirectory.empty()) {
    config->Write("checkpointDirectory", checkpointDirectory);
  }
  config->Write("onCloseSaveActiveModelToPdb", onCloseSaveActiveModelToPdb);

  WriteProfiles();

  config->Flush();
}

void Application::Init() {

  SetMolimageHome();
  SetCrystalData();
  std::string tmp;
  bool ret = config->Read("XFITDICT", tmp);
  if (ret)
    XFitDictSetting=tmp.c_str();
  SetDictionary(XFitDictSetting);
  PentDir = MolimageHome;
#ifndef _WIN32
  PentDir += "/data/pdbvec";
#else
  PentDir += "\\data\\pdbvec";
#endif

  tmp="";
  ret = config->Read("HTMLBROWSER", tmp);
  if (ret)
    HTMLBrowser=tmp.c_str();
  if (!ret) {
#ifndef _WIN32
#ifdef __APPLE__
    HTMLBrowser = "open";
#else
    HTMLBrowser = "mozilla";
#endif
#else
    HTMLBrowser = "C:\\Program Files\\Internet Explorer\\iexplore.exe";
#endif
  }
  tmp="";
  if (config->Read("SmilesDbCommand", tmp))
    SmilesDbCommand=tmp.c_str();
  tmp="";
  if (config->Read("SHELXHOME", tmp))
    ShelxHome=tmp.c_str();
  tmp="";
  if (config->Read("checkpointDirectory", tmp))
    checkpointDirectory=tmp.c_str();
  long boolValue;
  config->Read("onCloseSaveActiveModelToPdb", &boolValue, 0);
  onCloseSaveActiveModelToPdb = boolValue != 0;
  SetEnv();

  xfitMouseMode = config->GetProfileInt("Options", "MouseMode", 0) != 0;
  incrementallyColorModels = config->GetProfileInt("Options", "IncrementallyColorModels", 1) != 0;
  dimNonactiveModels = config->GetProfileInt("Options", "DimNonactiveModels", 1) != 0;
  concurrentJobLimit = config->GetProfileInt("Options", "ConcurrentJobLimit", 1);

  SetGammaCorrection(1.0);
  BuildPalette();

  ResidueBuffer = NULL;

  LabelToggle = config->GetProfileInt("View Parameters", "LabelToggle", 1) != 0;
  LabelPicks = config->GetProfileInt("View Parameters", "LabelPicks", 1);

  ATOMLABEL::defaultStyle(config->GetProfileInt("View Parameters", "LabelStyle", 0));
  ATOMLABEL::defaultColor(config->GetProfileInt("View Parameters", "LabelColorRed", 255),
    config->GetProfileInt("View Parameters", "LabelColorGreen", 255),
    config->GetProfileInt("View Parameters", "LabelColorBlue", 255));
  ATOMLABEL::defaultSize(config->GetProfileInt("View Parameters", "LabelSize", 12));


  int natomnames = config->GetProfileInt("Atom Colors", "NoEntries", 0);
  if (natomnames == 0) {
    init_colornames();
  } else {
    Colors::atomnames.clear();
    Colors::atomcolors.clear();
    for (int i = 0; i < natomnames; i++) {
      std::string b = format("Color%d", i);
      int ci = config->GetProfileInt("Atom Colors", b.c_str(), Colors::WHITE);
      Colors::atomcolors.push_back(Colors::colornames[ci]);
      b = format("Name%d", i);
      b = config->GetProfileString("Atom Colors", b, "");
      Colors::atomnames.push_back(b);
    }
  }

  LoadDictionary();

  // enable setting of colors when reading files
  MIRegisterColorSetter(new myColorSetter());
  MISetTorsionWritePrompt(new myTorsionWritePrompt());
  MISetMolPrefsHandler(new myMolPrefsHandler());
}

//FIXME: previous dictionary is leaked!
void Application::LoadDictionary() {
  geomrefiner = new GeomRefiner();
  geomrefiner->dict.LoadDefaultDictionary(
    Application::instance()->getDictionary().c_str(),
    Application::instance()->MolimageHome.c_str());
  MISetDictionary(&geomrefiner->dict);
}

void Application::saveDict() {
  QFileInfo fi(Application::instance()->XFitDictSetting.c_str());
  const std::string& s = MIFileSelector("Choose a dictionary file",
                                        fi.path().toStdString(), "mydict.pdb", "pdb",
                                        "Dictionary files (*.pdb)|*.pdb|All files (*.*)|*.*",
                                        MI_SAVE_MODE);
  if (s.size()) {
    if (MIFitDictionary()->SaveDictionary(s.c_str()) ) {
      if (s != std::string(Application::instance()->XFitDictSetting).c_str()) {
        if (MIMessageBox("Do you want to make this the default dictionary?", "Change dictionary",
              MIDIALOG_YES_NO) == MI_YES) {
          Application::instance()->XFitDictSetting = s.c_str();
        }
      }
    }
  }
}



void Application::WriteProfiles() {
  WritePalette();
  WriteSecStr();
  WriteBValues();
  WriteAtomTypes();
}

void Application::WritePalette() {
  config->WriteProfileInt("Palette", "Red.red", lpal->colors[Colors::PRED].red);
  config->WriteProfileInt("Palette", "Red.green", lpal->colors[Colors::PRED].green);
  config->WriteProfileInt("Palette", "Red.blue", lpal->colors[Colors::PRED].blue);
  config->WriteProfileInt("Palette", "Blue.red", lpal->colors[Colors::PBLUE].red);
  config->WriteProfileInt("Palette", "Blue.green", lpal->colors[Colors::PBLUE].green);
  config->WriteProfileInt("Palette", "Blue.blue", lpal->colors[Colors::PBLUE].blue);
  config->WriteProfileInt("Palette", "Green.red", lpal->colors[Colors::PGREEN].red);
  config->WriteProfileInt("Palette", "Green.green", lpal->colors[Colors::PGREEN].green);
  config->WriteProfileInt("Palette", "Green.blue", lpal->colors[Colors::PGREEN].blue);
  config->WriteProfileInt("Palette", "Cyan.red", lpal->colors[Colors::PCYAN].red);
  config->WriteProfileInt("Palette", "Cyan.green", lpal->colors[Colors::PCYAN].green);
  config->WriteProfileInt("Palette", "Cyan.blue", lpal->colors[Colors::PCYAN].blue);
  config->WriteProfileInt("Palette", "Magenta.red", lpal->colors[Colors::PMAGENTA].red);
  config->WriteProfileInt("Palette", "Magenta.green", lpal->colors[Colors::PMAGENTA].green);
  config->WriteProfileInt("Palette", "Magenta.blue", lpal->colors[Colors::PMAGENTA].blue);
  config->WriteProfileInt("Palette", "Yellow.red", lpal->colors[Colors::PYELLOW].red);
  config->WriteProfileInt("Palette", "Yellow.green", lpal->colors[Colors::PYELLOW].green);
  config->WriteProfileInt("Palette", "Yellow.blue", lpal->colors[Colors::PYELLOW].blue);
  config->WriteProfileInt("Palette", "White.red", lpal->colors[Colors::PWHITE].red);
  config->WriteProfileInt("Palette", "White.green", lpal->colors[Colors::PWHITE].green);
  config->WriteProfileInt("Palette", "White.blue", lpal->colors[Colors::PWHITE].blue);
  config->WriteProfileInt("Palette", "Black.red", lpal->colors[Colors::PBLACK].red);
  config->WriteProfileInt("Palette", "Black.green", lpal->colors[Colors::PBLACK].green);
  config->WriteProfileInt("Palette", "Black.blue", lpal->colors[Colors::PBLACK].blue);

  config->WriteProfileInt("User Colors", "User1.red", lpal->colors[Colors::PCUSTOM1].red);
  config->WriteProfileInt("User Colors", "User1.green", lpal->colors[Colors::PCUSTOM1].green);
  config->WriteProfileInt("User Colors", "User1.blue", lpal->colors[Colors::PCUSTOM1].blue);
  config->WriteProfileInt("User Colors", "User2.red", lpal->colors[Colors::PCUSTOM2].red);
  config->WriteProfileInt("User Colors", "User2.green", lpal->colors[Colors::PCUSTOM2].green);
  config->WriteProfileInt("User Colors", "User2.blue", lpal->colors[Colors::PCUSTOM2].blue);
  config->WriteProfileInt("User Colors", "User3.red", lpal->colors[Colors::PCUSTOM3].red);
  config->WriteProfileInt("User Colors", "User3.green", lpal->colors[Colors::PCUSTOM3].green);
  config->WriteProfileInt("User Colors", "User3.blue", lpal->colors[Colors::PCUSTOM3].blue);
  config->WriteProfileInt("User Colors", "User4.red", lpal->colors[Colors::PCUSTOM4].red);
  config->WriteProfileInt("User Colors", "User4.green", lpal->colors[Colors::PCUSTOM4].green);
  config->WriteProfileInt("User Colors", "User4.blue", lpal->colors[Colors::PCUSTOM4].blue);
  config->WriteProfileInt("User Colors", "User5.red", lpal->colors[Colors::PCUSTOM5].red);
  config->WriteProfileInt("User Colors", "User5.green", lpal->colors[Colors::PCUSTOM5].green);
  config->WriteProfileInt("User Colors", "User5.blue", lpal->colors[Colors::PCUSTOM5].blue);
  config->WriteProfileInt("User Colors", "User6.red", lpal->colors[Colors::PCUSTOM6].red);
  config->WriteProfileInt("User Colors", "User6.green", lpal->colors[Colors::PCUSTOM6].green);
  config->WriteProfileInt("User Colors", "User6.blue", lpal->colors[Colors::PCUSTOM6].blue);
  config->WriteProfileInt("View Parameters", "BackgroundColorRed", BackgroundColor.red);
  config->WriteProfileInt("View Parameters", "BackgroundColorGreen", BackgroundColor.green);
  config->WriteProfileInt("View Parameters", "BackgroundColorBlue", BackgroundColor.blue);
  config->WriteProfileInt("View Parameters", "LabelPicks", LabelPicks);
  config->WriteProfileInt("View Parameters", "LabelToggle", LabelToggle);
  config->WriteProfileInt("View Parameters", "LabelStyle", ATOMLABEL::defaultStyle());
  config->WriteProfileInt("View Parameters", "LabelSize", ATOMLABEL::defaultSize());
  config->WriteProfileInt("View Parameters", "LabelColorRed", ATOMLABEL::defaultRed());
  config->WriteProfileInt("View Parameters", "LabelColorGreen", ATOMLABEL::defaultGreen());
  config->WriteProfileInt("View Parameters", "LabelColorBlue", ATOMLABEL::defaultBlue());
}

void Application::WriteSecStr() {
  config->WriteProfileInt("SecStr Colors", "Helix", Colors::sec_colors[Colors::HELIX]);
  config->WriteProfileInt("SecStr Colors", "Sheet", Colors::sec_colors[Colors::SHEET]);
  config->WriteProfileInt("SecStr Colors", "Coil", Colors::sec_colors[Colors::COIL]);
}

void Application::WriteAtomTypes() {
  config->WriteProfileInt("Atom Colors", "NoEntries", Colors::atomcolors.size());
  for (unsigned int i = 0; i < Colors::atomcolors.size(); i++) {
    std::string b = format("Color%d", i);
    config->WriteProfileInt("Atom Colors", b.c_str(), Colors::findColorNumber(Colors::atomcolors[i]));
    b = format("Name%d", i);
    config->WriteProfileString("Atom Colors", b.c_str(), Colors::atomnames[i]);
  }
}

void Application::WriteBValues() {
  config->WriteProfileInt("BValue Colors", "Color1", Colors::BValueColors[0]);
  config->WriteProfileInt("BValue Colors", "Range1", Colors::BValueRanges[0]);
  config->WriteProfileInt("BValue Colors", "Color2", Colors::BValueColors[1]);
  config->WriteProfileInt("BValue Colors", "Range2", Colors::BValueRanges[1]);
  config->WriteProfileInt("BValue Colors", "Color3", Colors::BValueColors[2]);
  config->WriteProfileInt("BValue Colors", "Range3", Colors::BValueRanges[2]);
  config->WriteProfileInt("BValue Colors", "Color4", Colors::BValueColors[3]);
  config->WriteProfileInt("BValue Colors", "Range4", Colors::BValueRanges[3]);
  config->WriteProfileInt("BValue Colors", "Color5", Colors::BValueColors[4]);
  config->WriteProfileInt("BValue Colors", "Range5", Colors::BValueRanges[4]);
  config->WriteProfileInt("BValue Colors", "Color6", Colors::BValueColors[5]);
  config->WriteProfileInt("BValue Colors", "Range6", Colors::BValueRanges[5]);
  config->WriteProfileInt("BValue Colors", "Color7", Colors::BValueColors[6]);
  config->WriteProfileInt("BValue Colors", "Range7", Colors::BValueRanges[6]);
  config->WriteProfileInt("BValue Colors", "Color8", Colors::BValueColors[7]);
  config->WriteProfileInt("BValue Colors", "Range8", Colors::BValueRanges[7]);
  config->WriteProfileInt("BValue Colors", "Color9", Colors::BValueColors[8]);
  config->WriteProfileInt("BValue Colors", "Range9", Colors::BValueRanges[8]);
  config->WriteProfileInt("BValue Colors", "Color10", Colors::BValueColors[9]);
  config->WriteProfileInt("BValue Colors", "Range10", Colors::BValueRanges[9]);
}


const std::string& Application::GetBinDirectory() {
  MIbin = findMIFitBin();
  return MIbin;
}

void Application::SetMolimageHome(bool reset) {
  if (reset || MolimageHome.empty()) {
    MolimageHome = findMIFitBin();
    if (strncasecmp(MolimageHome.c_str(),"/bin",4) == 0 ||
        strncasecmp(MolimageHome.c_str(),"\\bin",4) == 0) {
      MolimageHome.erase(MolimageHome.size()-4);
    }
#ifdef __APPLE__
    MolimageHome+="/../Resources";
#endif
  }
}

void Application::SetShelxHome(bool reset) {
  std::string tmp;
  if (config->Read("SHELXHOME", tmp))
    ShelxHome=tmp.c_str();
  if (reset || ShelxHome.empty() || !QFileInfo(ShelxHome.c_str()).exists()) {
    MISelectDirectoryDialog dlg(NULL, "Find the Shelx home directory");
    MIData data;
    if (dlg.GetResults(data)) {
      ShelxHome = data["dir"].str.c_str();
    }
  }
}

void Application::SetEnv() {
  std::string t = format("MOLIMAGEHOME=%s", MolimageHome.c_str());
  t += MolimageHome;

  //NOTE: on linux, at least putenv expects a char*, not a const char* as
  //the pointed-to block becomes part of the environment, and altering the
  //block alters the environment.  A prior version passed t.c_str()
  //directly, and it's a wonder this ever worked, b/c when t goes out of
  //scope, the environment is now pointing to garbage.
  //
  // we leak a slight amount of memory here to avoid this by calling strdup
  putenv(strdup(t.c_str()));
  t = format("CRYSTALDATA=%s", CrystalData.c_str());
  putenv(strdup(t.c_str()));
}

void Application::SetCrystalData(bool reset) {
  if (reset || CrystalData.empty()) {
    std::string tmp;
    if(config->Read("CRYSTALDATA", tmp))
      CrystalData=tmp.c_str();
    if (CrystalData.empty()) {
#ifndef _WIN32
      CrystalData = TildeExpand("~/xtal_info");
#else
      CrystalData = MolimageHome;
      CrystalData += "\\data\\xtal_info";
#endif
    }
  }
  SetEnv();
}


std::string Application::GetCrystalData(const char* crystal, const char* key) {
  std::string data;
  FILE* fp;
  char buff[3200];
  char filename[BUFSIZ];

  if (crystal == NULL) {
    return data;
  }
  if (key == NULL) {
    return data;
  }
  if (crystal[0] == '\0') {
    return data;
  }
  if (key[0] == '\0') {
    return data;
  }

  // first look in current directory for an override
  // make sure it is not a directory
  strcpy(filename, crystal);
  if (filename[0] == '~') {
    strcpy(filename, TildeExpand(filename));
  }
  if ((fp = fopen(crystal, "r")) == NULL) {
    strcpy(filename, CrystalData.c_str());
    strcat(filename, "/");
    strcat(filename, crystal);
    if ((fp = fopen(filename, "r")) == NULL) {
      Logger::message("Error: can't find file for crystal: %s", crystal);
      return data;
    }
  }
  while (fgets(buff, 3200, fp) != NULL) {
    if (strncmp(key, buff, strlen(key)) == 0) {
      /* find position of first space */
      data = buff;
      data = MIAfterFirst(data,' ');
      MIStringTrim(data,false);
      MIStringTrim(data,true);
      break;
    }
  }
  fclose(fp);
  return data;
}

bool Application::GetCrystalCell(const char* crystal, float& a, float& b, float& c,
                                 float& alpha, float& beta, float& gamma) {
  if (crystal == NULL) {
    return false;
  }
  if (crystal[0] == '\0') {
    return false;
  }

  std::string line = GetCrystalData(crystal, "cell");
  if (sscanf(line.c_str(), "%f%f%f%f%f%f", &a, &b, &c, &alpha, &beta, &gamma) == 6) {
    return true;
  }
  return false;
}

int Application::GetCrystalNCSSymmops(const char* crystal, float ncrsymm[MISymmop::MAXNCRSYMMOPS][12]) {
  int nops = 0;
  std::string key;
  std::string buf;

  for (unsigned int i = 1; i <= MISymmop::MAXNCRSYMMOPS; i++) {
    key=format("ncrsymm%d", i);
    buf = GetCrystalData(crystal, key.c_str());
    if (buf.size() > 0) {
      sscanf(buf.c_str(), "%f%f%f%f%f%f%f%f%f%f%f%f", &ncrsymm[nops][0], &ncrsymm[nops][1], &ncrsymm[nops][2],
        &ncrsymm[nops][3], &ncrsymm[nops][4], &ncrsymm[nops][5], &ncrsymm[nops][6], &ncrsymm[nops][7],
        &ncrsymm[nops][8], &ncrsymm[nops][9], &ncrsymm[nops][10], &ncrsymm[nops][11]);
      nops++;
    } else {break;}
  }
  return (nops);
}

const std::string& Application::GetShelxHome() {
  if (ShelxHome.empty()) {
    SetShelxHome();
  }
  return ShelxHome;
}

bool Application::SetDictionary(const std::string& dict) {
  bool valueChanged = false;
  std::string newDict = dict;
  if (newDict.empty()) {
    newDict = MolimageHome;
#ifndef _WIN32
    newDict += "/data/dict.noh.pdb";
#else
    newDict += "\\data\\dict.noh.pdb";
#endif
  }
  if (XFitDict != newDict) {
    XFitDict = newDict;
    valueChanged = true;
  }
  return valueChanged;
}

std::string Application::getDictionary() {
  return XFitDict;
}


void Application::toggleStereo() {
  bool stereo = config->GetProfileInt("View Parameters", "stereo", 0) != 0;
  config->WriteProfileInt("View Parameters", "stereo", !stereo);
}

long Application::GetPid() {
  return getpid();
}


//FIXME: for Qt version, try to get a hardware stereo buffer at startup.
// If we can get one, do there's no saying we have to *use* the stereo
// planes.  This will avoid the silly "re-open to enable stereo"
// requirement.
void Application::toggleHardwareStereo() {
  bool hardwareStereo = config->GetProfileInt("View Parameters", "hardwareStereo", 0) != 0;
  config->WriteProfileInt("View Parameters", "hardwareStereo", !hardwareStereo);
//  MIMessageBox("Changing the hardware stereo setting does not affect\n"
//    "open document or dictionary windows. Re-open documents\n"
//    "to enable the new setting.", "Hardware Stereo Notice", MIDIALOG_ICON_WARNING, this);
}



MIBusyManager *MIBusyManager::_instance=0;
#define QPROGRESS_MAX 100
//#define QPROGRESS_MAX 0

MIBusyManager *MIBusyManager::instance() {
  if (!_instance)
    _instance=new MIBusyManager();
  return _instance;
}


MIBusyManager::MIBusyManager() {
  m_cursor_number=0;
  m_busy = 0;
  MustAbortNow=0;
  PROGRESS=new QProgressDialog("MIFit is busy","Abort",0,QPROGRESS_MAX);
  PROGRESS->setWindowModality(Qt::ApplicationModal);
}

bool MIBusyManager::Busy() const {
  return m_busy > 0;
}

void MIBusyManager::SetBusy(bool busy) {
  if (busy) {
    if (m_busy == 0) {
      MustAbortNow = false;
      PROG_VALUE=0;
      PROG_DIRECTION=1;
      PROGRESS->reset();
      PROGRESS->setValue(0);
      PROGRESS->show();
    }
    m_busy++;
  } else {
    if (m_busy > 0) {
      m_busy--;
    } else {
#ifdef DEBUG
      //FIXME Log("Warning: End of busy without corresponding start busy");
#endif
    }
  }
  if (m_busy == 0) {
    PROGRESS->reset();
    PROGRESS->setLabelText("MIFit is busy");
  }
}

bool MIBusyManager::CheckForAbort() {
  if (PROG_VALUE == 0) {
    PROG_DIRECTION=1;
  } else if (PROG_VALUE + 1 == QPROGRESS_MAX) {
    PROG_DIRECTION=-1;
  }

  PROG_VALUE = (PROG_VALUE + PROG_DIRECTION);
  PROGRESS->setValue(PROG_VALUE);
  if (PROGRESS->wasCanceled() || MustAbortNow)
  {
    return true;
  }
  return false;
}


void MIBusyManager::SetLabel(const char *op)
{
  PROGRESS->setLabelText(op);
}

void MIBusyManager::ForceAbort() {
  MustAbortNow = true;
}

void MIBusyManager::SetWaitCursor() {
  switch (m_cursor_number) {
    case 0: SetCursor(imhWait1); break;
    case 1: SetCursor(imhWait2); break;
    case 2: SetCursor(imhWait3); break;
    case 3: SetCursor(imhWait4); break;
    case 4: SetCursor(imhWait5); break;
    case 5: SetCursor(imhWait6); break;
  }
  m_cursor_number=(m_cursor_number+1)%6;
}

void MIBusyManager::StartWaitCursor() {
  SetWaitCursor();
}

void MIBusyManager::StopWaitCursor() {
  SetCursor(-1);
}

void MIBusyManager::SetCursor(int id) {
  MIMainWindow::instance()->SetCursor(id);
}


GeomRefiner* MIFitGeomRefiner() {
  return Application::instance()->GetGeomRefiner();
}

chemlib::MIMolDictionary* MIFitDictionary() {
  return &MIFitGeomRefiner()->dict;
}

QSettings *MIGetQSettings() {
  return ((MIConfigImpl*)MIConfig::Instance())->GetQSettings();
}

QString Application::latestFileBrowseDirectory(const QString& path)
{
    const QString LATEST_FILE_BROWSE_DIRECTORY("latestFileBrowseDirectory");
    QSettings* settings = MIGetQSettings();
    if (settings) {
        if (path.isEmpty()) {
            QString settingsPath = settings->value(LATEST_FILE_BROWSE_DIRECTORY).toString();
            if (QFile::exists(settingsPath)) {
                return settingsPath;
            }
        } else {
            settings->setValue(LATEST_FILE_BROWSE_DIRECTORY, path);
        }
    }
    return path;
}

QString Application::getOpenFileName(QWidget* parent, const QString& caption, const QString& filter)
{
    QString dir = Application::instance()->latestFileBrowseDirectory("");
    QString file = QFileDialog::getOpenFileName(parent, caption, dir, filter);
    if (!file.isEmpty()) {
        QFileInfo fileInfo(file);
        Application::instance()->latestFileBrowseDirectory(fileInfo.absolutePath());
    }
    return file;
}

QString Application::getExistingDirectory(QWidget* parent, const QString& caption, const QString& filter)
{
    QString dir = Application::instance()->latestFileBrowseDirectory("");
    dir = QFileDialog::getExistingDirectory(parent, caption, dir);
    if (!dir.isEmpty()) {
        Application::instance()->latestFileBrowseDirectory(dir);
    }
    return dir;
}
