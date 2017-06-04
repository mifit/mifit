#include <sys/types.h>
#include <sys/stat.h>
#include <errno.h>


#include <nongui/nonguilib.h>

#include <chemlib/chemlib.h>
#include <chemlib/Monomer.h>
#include <conflib/conflib.h>
#include <util/utillib.h>
#include <core/corelib.h>

#include "Application.h"
#include "tools.h"
#include "DictEditCanvas.h"
#include "GLFormatEdit.h"

#include <util/FileIo.h>
#include "MIMainWindow.h"
#include "molw.h" // for cursors

#include <QApplication>
#include <QColorDialog>
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


using namespace chemlib;
using namespace std;

/*
 * tilde expand a file name
 *
 *   returns: the expanded name of the file (in a static buffer)
 *            or NIL(char)
 */
const char *TildeExpand(const char *filename)
{
#ifdef _WIN32
    /* tilde's don't make sense in windows */
    return filename;
#else
    static char dummy[1024];
    char username[20];
    const char *loc;
    struct passwd *pw;

    if ((filename == NULL) || !strcmp(filename, ""))
    {
        return ("");
    }

    if (filename[0] != '~')
    {
        (void) strcpy(dummy, filename);
        return (dummy);
    }

    /* tilde at the beginning now */
    if ((filename[1] == '\0') || (filename[1] == '/') || (filename[1] == ' '))
    {
        /* current user */
        char *home;

        if ((home = getenv("HOME")) == (char*) 0)
        {
            /* fall back on /etc/passwd */
            if ((pw = getpwuid(getuid())) == (struct passwd*) NULL)
            {
                return ("");
            }
            strcpy(dummy, pw->pw_dir);
        }
        else
        {
            strcpy(dummy, home);
        }
        strcat(dummy, &filename[1]);

    }
    else
    {
        loc = strchr(filename, '/');
        if (loc == 0)
        {
            strcpy(username, &filename[1]);
        }
        else
        {
            (void) strncpy(username, &filename[1], loc - &filename[1]);
            username[loc - &filename[1]] = '\0';
        }
        if ((pw = getpwnam(username)) == NULL)
        {
            return ("");
        }
        strcpy(dummy, pw->pw_dir);
        if (loc != 0)
        {
            strcat(dummy, loc);
        }
    }
    return (dummy);
#endif // ifdef _WIN32
}


// Find the absolute path where this application has been run from.
std::string findMIFitBin()
{
    QFileInfo fi(QCoreApplication::applicationDirPath());
    return fi.canonicalFilePath().toStdString();
}

Application *Application::instance()
{
    return static_cast<Application*>(qApp);
}

namespace
{
const QString GL_FORMAT_GROUP("glformat");
const QString GL_FORMAT_DEFAULT("default");
}

Application::Application(int &argc, char **argv)
    : QApplication(argc, argv)
{

    static const char *name = "MIFit";
    setOrganizationName(name);
    setApplicationName(name);

    FileIo::setAsDefaultIo();

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
    if (!jobLogsDir.exists())
    {
        jobLogsDir.mkpath(".");
    }

    QSettings settings;
    if (settings.childGroups().contains(GL_FORMAT_GROUP))
    {
        settings.beginGroup(GL_FORMAT_GROUP);
        if (settings.childGroups().contains(GL_FORMAT_DEFAULT))
        {
            settings.beginGroup(GL_FORMAT_DEFAULT);
            QGLFormat glformat = GLFormatEdit::readSettings(settings);
            QGLFormat::setDefaultFormat(glformat);
            settings.endGroup();
        }
        settings.endGroup();
    }

    Init();
}

Application::~Application()
{
    QSettings settings;
    settings.setValue("Options/MouseMode", xfitMouseMode);
    settings.setValue("Options/IncrementallyColorModels", incrementallyColorModels);
    settings.setValue("Options/DimNonactiveModels", dimNonactiveModels);
    settings.setValue("Options/ConcurrentJobLimit", concurrentJobLimit);

    Write();

    settings.beginGroup(GL_FORMAT_GROUP);
    settings.beginGroup(GL_FORMAT_DEFAULT);
    QGLFormat glformat = QGLFormat::defaultFormat();
    GLFormatEdit::writeSettings(settings, glformat);
    settings.endGroup();
    settings.endGroup();

    delete lpal;
    ClearResidueBuffer();

    if (geomrefiner) //Cleanup the GeomRefiner
    {
        if (geomrefiner->dict.GetModified())
        {
            if (QMessageBox::question(0, "Dictionary Modified",
                                      "The Dictionary has been modified\nDo you want to save it?",
                                      QMessageBox::Yes | QMessageBox::No, QMessageBox::Yes) == QMessageBox::Yes)
            {
                saveDict();
            }
        }
    }
    delete geomrefiner;

}


int Application::ApplyGammaCorrection(const int color)
{
    double gamma = GetGammaCorrection();
    if (gamma <= 0.0)
    {
        gamma = 0.2;
    }
    double fcolor;
    gamma = 1.0/gamma;
    fcolor = (double)color/255.0;
    fcolor = pow(fcolor, gamma);
    fcolor *= 255.0;
    if (fcolor > 255.0)
    {
        return 255;
    }
    if (fcolor < 0)
    {
        return 0;
    }
    return ROUND(fcolor);
}

void Application::BuildPalette()
{
    QSettings settings;
    lpal = new MIPalette;
    lpal->colors.resize(Colors_NUMBERPALETTE);
    lpal->colors[Colors::PBLACK].red = settings.value("Palette/Black.red", 0).toInt();
    lpal->colors[Colors::PBLACK].green = settings.value("Palette/Black.green", 0).toInt();
    lpal->colors[Colors::PBLACK].blue = settings.value("Palette/Black.blue", 0).toInt();
    lpal->colors[Colors::PRED].red = settings.value("Palette/Red.red", 255).toInt();
    lpal->colors[Colors::PRED].green = settings.value("Palette/Red.green", 20).toInt();
    lpal->colors[Colors::PRED].blue = settings.value("Palette/Red.blue", 20).toInt();
    lpal->colors[Colors::PBLUE].red = settings.value("Palette/Blue.red", 20).toInt();
    lpal->colors[Colors::PBLUE].green = settings.value("Palette/Blue.green", 20).toInt();
    lpal->colors[Colors::PBLUE].blue = settings.value("Palette/Blue.blue", 255).toInt();
    lpal->colors[Colors::PGREEN].red = settings.value("Palette/Green.red", 20).toInt();
    lpal->colors[Colors::PGREEN].green = settings.value("Palette/Green.green", 255).toInt();
    lpal->colors[Colors::PGREEN].blue = settings.value("Palette/Green.blue", 20).toInt();
    lpal->colors[Colors::PCYAN].red = settings.value("Palette/Cyan.red", 20).toInt();
    lpal->colors[Colors::PCYAN].green = settings.value("Palette/Cyan.green", 255).toInt();
    lpal->colors[Colors::PCYAN].blue = settings.value("Palette/Cyan.blue", 255).toInt();
    lpal->colors[Colors::PMAGENTA].red = settings.value("Palette/Magenta.red", 255).toInt();
    lpal->colors[Colors::PMAGENTA].green = settings.value("Palette/Magenta.green", 20).toInt();
    lpal->colors[Colors::PMAGENTA].blue = settings.value("Palette/Magenta.blue", 255).toInt();
    lpal->colors[Colors::PYELLOW].red = settings.value("Palette/Yellow.red", 255).toInt();
    lpal->colors[Colors::PYELLOW].green = settings.value("Palette/Yellow.green", 255).toInt();
    lpal->colors[Colors::PYELLOW].blue = settings.value("Palette/Yellow.blue", 20).toInt();
    lpal->colors[Colors::PWHITE].red = settings.value("Palette/White.red", 253).toInt();
    lpal->colors[Colors::PWHITE].green = settings.value("Palette/White.green", 253).toInt();
    lpal->colors[Colors::PWHITE].blue = settings.value("Palette/White.blue", 253).toInt();
    lpal->colors[Colors::PPINK].red = settings.value("Palette/Pink.red", 255).toInt();
    lpal->colors[Colors::PPINK].green = settings.value("Palette/Pink.green", 128).toInt();
    lpal->colors[Colors::PPINK].blue = settings.value("Palette/Pink.blue", 128).toInt();
    lpal->colors[Colors::PORANGE].red = settings.value("Palette/Orange.red", 255).toInt();
    lpal->colors[Colors::PORANGE].green = settings.value("Palette/Orange.green", 128).toInt();
    lpal->colors[Colors::PORANGE].blue = settings.value("Palette/Orange.blue", 0).toInt();
    lpal->colors[Colors::PBROWN].red = settings.value("Palette/Brown.red", 0).toInt();
    lpal->colors[Colors::PBROWN].green = settings.value("Palette/Brown.green", 0).toInt();
    lpal->colors[Colors::PBROWN].blue = settings.value("Palette/Brown.blue", 0).toInt();

    lpal->colors[Colors::PCUSTOM1].red = settings.value("User Colors/User1.red", 255).toInt();
    lpal->colors[Colors::PCUSTOM1].green = settings.value("User Colors/User1.green", 72).toInt();
    lpal->colors[Colors::PCUSTOM1].blue = settings.value("User Colors/User1.blue", 96).toInt();
    lpal->colors[Colors::PCUSTOM2].red = settings.value("User Colors/User2.red", 90).toInt();
    lpal->colors[Colors::PCUSTOM2].green = settings.value("User Colors/User2.green", 255).toInt();
    lpal->colors[Colors::PCUSTOM2].blue = settings.value("User Colors/User2.blue", 90).toInt();
    lpal->colors[Colors::PCUSTOM3].red = settings.value("User Colors/User3.red", 102).toInt();
    lpal->colors[Colors::PCUSTOM3].green = settings.value("User Colors/User3.green", 164).toInt();
    lpal->colors[Colors::PCUSTOM3].blue = settings.value("User Colors/User3.blue", 255).toInt();
    lpal->colors[Colors::PCUSTOM4].red = settings.value("User Colors/User4.red", 162).toInt();
    lpal->colors[Colors::PCUSTOM4].green = settings.value("User Colors/User4.green", 115).toInt();
    lpal->colors[Colors::PCUSTOM4].blue = settings.value("User Colors/User4.blue", 20).toInt();
    lpal->colors[Colors::PCUSTOM5].red = settings.value("User Colors/User5.red", 255).toInt();
    lpal->colors[Colors::PCUSTOM5].green = settings.value("User Colors/User5.green", 143).toInt();
    lpal->colors[Colors::PCUSTOM5].blue = settings.value("User Colors/User5.blue", 32).toInt();
    lpal->colors[Colors::PCUSTOM6].red = settings.value("User Colors/User6.red", 255).toInt();
    lpal->colors[Colors::PCUSTOM6].green = settings.value("User Colors/User6.green", 143).toInt();
    lpal->colors[Colors::PCUSTOM6].blue = settings.value("User Colors/User6.blue", 190).toInt();
    lpal->colors[Colors::PCUSTOM7].red = settings.value("User Colors/User7.red", 150).toInt();
    lpal->colors[Colors::PCUSTOM7].green = settings.value("User Colors/User7.green", 0).toInt();
    lpal->colors[Colors::PCUSTOM7].blue = settings.value("User Colors/User7.blue", 0).toInt();
    lpal->colors[Colors::PCUSTOM8].red = settings.value("User Colors/User8.red", 0).toInt();
    lpal->colors[Colors::PCUSTOM8].green = settings.value("User Colors/User8.green", 150).toInt();
    lpal->colors[Colors::PCUSTOM8].blue = settings.value("User Colors/User8.blue", 0).toInt();
    lpal->colors[Colors::PCUSTOM9].red = settings.value("User Colors/User9.red", 0).toInt();
    lpal->colors[Colors::PCUSTOM9].green = settings.value("User Colors/User9.green", 0).toInt();
    lpal->colors[Colors::PCUSTOM9].blue = settings.value("User Colors/User9.blue", 150).toInt();
    lpal->colors[Colors::PCUSTOM10].red = settings.value("User Colors/User10.red", 150).toInt();
    lpal->colors[Colors::PCUSTOM10].green = settings.value("User Colors/User10.green", 150).toInt();

    lpal->colors[Colors::PCUSTOM10].blue = settings.value("User Colors/User10.blue", 0).toInt();

    lpal->colors[Colors::PMAP1].red = settings.value("Map Colors/Map1.red", 0).toInt();
    lpal->colors[Colors::PMAP1].green = settings.value("Map Colors/Map1.green", 0).toInt();
    lpal->colors[Colors::PMAP1].blue = settings.value("Map Colors/Map1.blue", 230).toInt();
    lpal->colors[Colors::PMAP2].red = settings.value("Map Colors/Map2.red", 128).toInt();
    lpal->colors[Colors::PMAP2].green = settings.value("Map Colors/Map2.green", 0).toInt();
    lpal->colors[Colors::PMAP2].blue = settings.value("Map Colors/Map2.blue", 230).toInt();
    lpal->colors[Colors::PMAP3].red = settings.value("Map Colors/Map3.red", 230).toInt();
    lpal->colors[Colors::PMAP3].green = settings.value("Map Colors/Map3.green", 0).toInt();
    lpal->colors[Colors::PMAP3].blue = settings.value("Map Colors/Map3.blue", 230).toInt();
    lpal->colors[Colors::PMAP4].red = settings.value("Map Colors/Map4.red", 230).toInt();
    lpal->colors[Colors::PMAP4].green = settings.value("Map Colors/Map4.green", 0).toInt();
    lpal->colors[Colors::PMAP4].blue = settings.value("Map Colors/Map4.blue", 128).toInt();
    lpal->colors[Colors::PMAP5].red = settings.value("Map Colors/Map5.red", 230).toInt();
    lpal->colors[Colors::PMAP5].green = settings.value("Map Colors/Map5.green", 0).toInt();
    lpal->colors[Colors::PMAP5].blue = settings.value("Map Colors/Map5.blue", 0).toInt();
    lpal->colors[Colors::PMAP6].red = settings.value("Map Colors/Map6.red", 0).toInt();
    lpal->colors[Colors::PMAP6].green = settings.value("Map Colors/Map6.green", 230).toInt();
    lpal->colors[Colors::PMAP6].blue = settings.value("Map Colors/Map6.blue", 0).toInt();
    lpal->colors[Colors::PMAP7].red = settings.value("Map Colors/Map7.red", 0).toInt();
    lpal->colors[Colors::PMAP7].green = settings.value("Map Colors/Map7.green", 230).toInt();
    lpal->colors[Colors::PMAP7].blue = settings.value("Map Colors/Map7.blue", 128).toInt();
    lpal->colors[Colors::PMAP8].red = settings.value("Map Colors/Map8.red", 0).toInt();
    lpal->colors[Colors::PMAP8].green = settings.value("Map Colors/Map8.green", 230).toInt();
    lpal->colors[Colors::PMAP8].blue = settings.value("Map Colors/Map8.blue", 230).toInt();
    lpal->colors[Colors::PMAP9].red = settings.value("Map Colors/Map9.red", 0).toInt();
    lpal->colors[Colors::PMAP9].green = settings.value("Map Colors/Map9.green", 128).toInt();
    lpal->colors[Colors::PMAP9].blue = settings.value("Map Colors/Map9.blue", 230).toInt();
    lpal->colors[Colors::PMAP10].red = settings.value("Map Colors/Map10.red", 128).toInt();
    lpal->colors[Colors::PMAP10].green = settings.value("Map Colors/Map10.green", 128).toInt();
    lpal->colors[Colors::PMAP10].blue = settings.value("Map Colors/Map10.blue", 230).toInt();

    lpal->colors[Colors::PCONTOUR1].red = settings.value("Contour Colors/Contour1.red", 218).toInt();
    lpal->colors[Colors::PCONTOUR1].green = settings.value("Contour Colors/Contour1.green", 152).toInt();
    lpal->colors[Colors::PCONTOUR1].blue = settings.value("Contour Colors/Contour1.blue", 207).toInt();
    lpal->colors[Colors::PCONTOUR2].red = settings.value("Contour Colors/Contour2.red", 113).toInt();
    lpal->colors[Colors::PCONTOUR2].green = settings.value("Contour Colors/Contour2.green", 87).toInt();
    lpal->colors[Colors::PCONTOUR2].blue = settings.value("Contour Colors/Contour2.blue", 185).toInt();
    lpal->colors[Colors::PCONTOUR3].red = settings.value("Contour Colors/Contour3.red", 72).toInt();
    lpal->colors[Colors::PCONTOUR3].green = settings.value("Contour Colors/Contour3.green", 107).toInt();
    lpal->colors[Colors::PCONTOUR3].blue = settings.value("Contour Colors/Contour3.blue", 254).toInt();
    lpal->colors[Colors::PCONTOUR4].red = settings.value("Contour Colors/Contour4.red", 75).toInt();
    lpal->colors[Colors::PCONTOUR4].green = settings.value("Contour Colors/Contour4.green", 230).toInt();
    lpal->colors[Colors::PCONTOUR4].blue = settings.value("Contour Colors/Contour4.blue", 251).toInt();
    lpal->colors[Colors::PCONTOUR5].red = settings.value("Contour Colors/Contour5.red", 0).toInt();
    lpal->colors[Colors::PCONTOUR5].green = settings.value("Contour Colors/Contour5.green", 255).toInt();
    lpal->colors[Colors::PCONTOUR5].blue = settings.value("Contour Colors/Contour5.blue", 0).toInt();
    lpal->colors[Colors::PCONTOUR6].red = settings.value("Contour Colors/Contour6.red", 255).toInt();
    lpal->colors[Colors::PCONTOUR6].green = settings.value("Contour Colors/Contour6.green", 255).toInt();
    lpal->colors[Colors::PCONTOUR6].blue = settings.value("Contour Colors/Contour6.blue", 30).toInt();
    lpal->colors[Colors::PCONTOUR7].red = settings.value("Contour Colors/Contour7.red", 255).toInt();
    lpal->colors[Colors::PCONTOUR7].green = settings.value("Contour Colors/Contour7.green", 200).toInt();
    lpal->colors[Colors::PCONTOUR7].blue = settings.value("Contour Colors/Contour7.blue", 30).toInt();
    lpal->colors[Colors::PCONTOUR8].red = settings.value("Contour Colors/Contour8.red", 245).toInt();
    lpal->colors[Colors::PCONTOUR8].green = settings.value("Contour Colors/Contour8.green", 141).toInt();
    lpal->colors[Colors::PCONTOUR8].blue = settings.value("Contour Colors/Contour8.blue", 3).toInt();
    lpal->colors[Colors::PCONTOUR9].red = settings.value("Contour Colors/Contour9.red", 255).toInt();
    lpal->colors[Colors::PCONTOUR9].green = settings.value("Contour Colors/Contour9.green", 60).toInt();
    lpal->colors[Colors::PCONTOUR9].blue = settings.value("Contour Colors/Contour9.blue", 00).toInt();
    lpal->colors[Colors::PCONTOUR10].red = settings.value("Contour Colors/Contour10.red", 220).toInt();
    lpal->colors[Colors::PCONTOUR10].green = settings.value("Contour Colors/Contour10.green", 0).toInt();
    lpal->colors[Colors::PCONTOUR10].blue = settings.value("Contour Colors/Contour10.blue", 0).toInt();

    BackgroundColor = PaletteColor(settings.value("View Parameters/BackgroundColorRed", 0).toInt(), settings.value("View Parameters/BackgroundColorGreen", 0).toInt(), settings.value("View Parameters/BackgroundColorBlue", 0).toInt());
    MixBackgroundColor();

    Colors::sec_colors[Colors::HELIX] = settings.value("SecStr Colors/Helix", Colors::MAGENTA).toInt();
    Colors::sec_colors[Colors::SHEET] = settings.value("SecStr Colors/Sheet", Colors::YELLOW).toInt();
    Colors::sec_colors[Colors::COIL] = settings.value("SecStr Colors/Coil", Colors::BLUE).toInt();

    Colors::BValueColors[0] = settings.value("BValue Colors/Color1", Colors::CONTOUR1).toInt();
    Colors::BValueRanges[0] = settings.value("BValue Colors/Range1", 300).toInt();
    Colors::BValueColors[1] = settings.value("BValue Colors/Color2", Colors::CONTOUR2).toInt();
    Colors::BValueRanges[1] = settings.value("BValue Colors/Range2", 600).toInt();
    Colors::BValueColors[2] = settings.value("BValue Colors/Color3", Colors::CONTOUR3).toInt();
    Colors::BValueRanges[2] = settings.value("BValue Colors/Range3", 1000).toInt();
    Colors::BValueColors[3] = settings.value("BValue Colors/Color4", Colors::CONTOUR4).toInt();
    Colors::BValueRanges[3] = settings.value("BValue Colors/Range4", 1500).toInt();
    Colors::BValueColors[4] = settings.value("BValue Colors/Color5", Colors::CONTOUR5).toInt();
    Colors::BValueRanges[4] = settings.value("BValue Colors/Range5", 2000).toInt();
    Colors::BValueColors[5] = settings.value("BValue Colors/Color6", Colors::CONTOUR6).toInt();
    Colors::BValueRanges[5] = settings.value("BValue Colors/Range6", 2500).toInt();
    Colors::BValueColors[6] = settings.value("BValue Colors/Color7", Colors::CONTOUR7).toInt();
    Colors::BValueRanges[6] = settings.value("BValue Colors/Range7", 3500).toInt();
    Colors::BValueColors[7] = settings.value("BValue Colors/Color8", Colors::CONTOUR8).toInt();
    Colors::BValueRanges[7] = settings.value("BValue Colors/Range8", 5000).toInt();
    Colors::BValueColors[8] = settings.value("BValue Colors/Color9", Colors::CONTOUR9).toInt();
    Colors::BValueRanges[8] = settings.value("BValue Colors/Range9", 7500).toInt();
    Colors::BValueColors[9] = settings.value("BValue Colors/Color10", Colors::CONTOUR10).toInt();
    Colors::BValueRanges[9] = settings.value("BValue Colors/Range10", 10000).toInt();

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

void Application::MixBackgroundColor()
{
    int i, ir, ib, r, g, b, d, rcolor, gcolor, bcolor;
    if (!lpal)
    {
        return;
    }
    rcolor = GetBackgroundColor().red;
    gcolor = GetBackgroundColor().green;
    bcolor = GetBackgroundColor().blue;
    for (ir = 0; ir < Colors_NUMBERCOLORS; ir++)
    {
        Colors::RPallette[PaletteIndex(ir)] = lpal->colors[PaletteIndex(ir)].red;
        Colors::GPallette[PaletteIndex(ir)] = lpal->colors[PaletteIndex(ir)].green;
        Colors::BPallette[PaletteIndex(ir)] = lpal->colors[PaletteIndex(ir)].blue;
        for (d = 1; d < 10; d++)
        {
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
    for (i = 0; i < Colors_NUMBERPALETTE; i++)
    {
        Colors::RPallette[i] = ApplyGammaCorrection(Colors::RPallette[i]);
        Colors::GPallette[i] = ApplyGammaCorrection(Colors::GPallette[i]);
        Colors::BPallette[i] = ApplyGammaCorrection(Colors::BPallette[i]);
    }
}

void Application::backgroundColor()
{
    QColor initColor = QColor(BackgroundColor.red, BackgroundColor.green, BackgroundColor.blue);
    QColor c = QColorDialog::getColor(initColor, 0);
    if (!c.isValid())
    {
        return;
    }
    PaletteColor color;
    color.red = (unsigned char)c.red();
    color.green = (unsigned char)c.green();
    color.blue = (unsigned char)c.blue();

    SetBackgroundColor(color);
}

bool Application::CopyResidueBuffer(const Residue *buffer)
{
    Residue *newres = new Residue(*buffer);
    if (ResidueBuffer)
    {
        Residue *end = ResidueBuffer;
        while (end->next() != NULL)
        {
            end = end->next();
        }
        end->insertResidue(newres);
    }
    else
    {
        ResidueBuffer = newres;
    }
    return true;
}

void Application::ClearResidueBuffer()
{
    if (ResidueBuffer)
    {
        FreeResidueList(ResidueBuffer);
    }
    ResidueBuffer = NULL;
}




class myColorSetter
    : public MIColorSetter
{
public:
    bool operator()(MIAtom *atom) const
    {
        atom->setColor(color_by_name(atom->name()));
        if (atom->color() == Colors::BLACK)
        {
            atom->setColor(Colors::WHITE);
        }
        return true;
    }

    bool operator()(MIAtom *atom, char c) const
    {
        atom->setColor(color_by_name(atom->name(), c));
        if (atom->color() == Colors::BLACK)
        {
            atom->setColor(Colors::WHITE);
        }
        return true;
    }

};

class myTorsionWritePrompt
    : public MITorsionWritePrompt
{
public:
    bool operator()()
    {
        return QMessageBox::question(0, "Write Torsions?", "Write torsions to file?", QMessageBox::Yes | QMessageBox::No, QMessageBox::No) == QMessageBox::Yes;
    }

};




class myMolPrefsHandler
    : public MolPrefsHandler
{
public:
    void operator()(bool *breakByDiscontinuity, bool *breakByNonpeptide)
    {
        QSettings settings;
        *breakByDiscontinuity = settings.value("Options/breakByDiscontinuity", true).toBool();
        *breakByNonpeptide = settings.value("Options/breakByNonpeptide", false).toBool();
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


void Application::Write()
{
    QSettings settings;
    if (!CrystalData.empty())
    {
        settings.setValue("CRYSTALDATA", CrystalData.c_str());
    }
    settings.setValue("XFITDICT", XFitDictSetting.c_str());
    if (!HTMLBrowser.empty())
    {
        settings.setValue("HTMLBROWSER", HTMLBrowser.c_str());
    }
    if (!ShelxHome.empty())
    {
        settings.setValue("SHELXHOME", ShelxHome.c_str());
    }
    if (!SmilesDbCommand.empty())
    {
        settings.setValue("SmilesDbCommand", SmilesDbCommand.c_str());
    }
    if (!checkpointDirectory.empty())
    {
        settings.setValue("checkpointDirectory", checkpointDirectory.c_str());
    }
    settings.setValue("onCloseSaveActiveModelToPdb", onCloseSaveActiveModelToPdb);

    WriteProfiles();

}

void Application::Init()
{

    SetMolimageHome();
    SetCrystalData();
    QSettings settings;
    if (settings.contains("XFITDICT"))
        XFitDictSetting = settings.value("XFITDICT").toString().toStdString();
    SetDictionary(XFitDictSetting);
    PentDir = MolimageHome;
#ifndef _WIN32
    PentDir += "/data/pdbvec";
#else
    PentDir += "\\data\\pdbvec";
#endif

    if (settings.contains("HTMLBROWSER"))
        HTMLBrowser = settings.value("HTMLBROWSER").toString().toStdString();
    else
    {
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
    if (settings.contains("SmilesDbCommand"))
        SmilesDbCommand = settings.value("SmilesDbCommand").toString().toStdString();
    if (settings.contains("SHELXHOME"))
        ShelxHome = settings.value("SHELXHOME").toString().toStdString();
    if (settings.contains("checkpointDirectory"))
        checkpointDirectory = settings.value("checkpointDirectory").toString().toStdString();
    onCloseSaveActiveModelToPdb = settings.value("onCloseSaveActiveModelToPdb", false).toBool();
    SetEnv();

    xfitMouseMode = settings.value("Options/MouseMode", false).toBool();
    incrementallyColorModels = settings.value("Options/IncrementallyColorModels", true).toBool();
    dimNonactiveModels = settings.value("Options/DimNonactiveModels", true).toBool();
    concurrentJobLimit = settings.value("Options/ConcurrentJobLimit", 1).toInt();

    SetGammaCorrection(1.0);
    BuildPalette();

    ResidueBuffer = NULL;

    LabelToggle = settings.value("View Parameters/LabelToggle", true).toBool();
    LabelPicks = settings.value("View Parameters/LabelPicks", true).toBool();

    ATOMLABEL::defaultStyle(settings.value("View Parameters/LabelStyle", 0).toInt());
    ATOMLABEL::defaultColor(settings.value("View Parameters/LabelColorRed", 255).toInt(),
                            settings.value("View Parameters/LabelColorGreen", 255).toInt(),
                            settings.value("View Parameters/LabelColorBlue", 255).toInt());
    ATOMLABEL::defaultSize(settings.value("View Parameters/LabelSize", 12).toInt());


    int natomnames = settings.value("Atom Colors/NoEntries", 0).toInt();
    if (natomnames == 0)
    {
        init_colornames();
    }
    else
    {
        Colors::atomnames.clear();
        Colors::atomcolors.clear();
        for (int i = 0; i < natomnames; i++)
        {
            std::string b = format("Color%d", i);
            int ci = settings.value(("Atom Colors/" + b).c_str(), Colors::WHITE).toInt();
            Colors::atomcolors.push_back(Colors::colornames[ci]);
            b = format("Name%d", i);
            b = settings.value(("Atom Colors/" + b).c_str(), "").toString().toStdString();
            Colors::atomnames.push_back(b);
        }
    }

    geomrefiner = new GeomRefiner();

    // enable setting of colors when reading files
    MIRegisterColorSetter(new myColorSetter());
    MISetTorsionWritePrompt(new myTorsionWritePrompt());
    MISetMolPrefsHandler(new myMolPrefsHandler());
}

void Application::AfterInit()
{
    LoadDictionary();
}

//FIXME: previous dictionary is leaked!
void Application::LoadDictionary()
{
    geomrefiner->dict.LoadDefaultDictionary(
        Application::instance()->getDictionary().c_str(),
        Application::instance()->MolimageHome.c_str());
    MISetDictionary(&geomrefiner->dict);
}

void Application::saveDict()
{
    QFileInfo fi(Application::instance()->XFitDictSetting.c_str());
    QString s = QFileDialog::getSaveFileName(0, "Choose a dictionary file",
                                          fi.filePath(), "Dictionary files (*.pdb);;All files (*.*)");
    if (!s.isEmpty())
    {
        if (MIFitDictionary()->SaveDictionary(s.toUtf8().constData()) )
        {
            if (s.toStdString() != Application::instance()->XFitDictSetting)
            {
                if (QMessageBox::question(0, "Change dictionary", "Do you want to make this the default dictionary?", QMessageBox::Yes | QMessageBox::No) == QMessageBox::Yes)
                {
                    Application::instance()->XFitDictSetting = s.toStdString();
                }
            }
        }
    }
}



void Application::WriteProfiles()
{
    WritePalette();
    WriteSecStr();
    WriteBValues();
    WriteAtomTypes();
}

void Application::WritePalette()
{
    QSettings settings;
    settings.setValue("Palette/Red.red", lpal->colors[Colors::PRED].red);
    settings.setValue("Palette/Red.green", lpal->colors[Colors::PRED].green);
    settings.setValue("Palette/Red.blue", lpal->colors[Colors::PRED].blue);
    settings.setValue("Palette/Blue.red", lpal->colors[Colors::PBLUE].red);
    settings.setValue("Palette/Blue.green", lpal->colors[Colors::PBLUE].green);
    settings.setValue("Palette/Blue.blue", lpal->colors[Colors::PBLUE].blue);
    settings.setValue("Palette/Green.red", lpal->colors[Colors::PGREEN].red);
    settings.setValue("Palette/Green.green", lpal->colors[Colors::PGREEN].green);
    settings.setValue("Palette/Green.blue", lpal->colors[Colors::PGREEN].blue);
    settings.setValue("Palette/Cyan.red", lpal->colors[Colors::PCYAN].red);
    settings.setValue("Palette/Cyan.green", lpal->colors[Colors::PCYAN].green);
    settings.setValue("Palette/Cyan.blue", lpal->colors[Colors::PCYAN].blue);
    settings.setValue("Palette/Magenta.red", lpal->colors[Colors::PMAGENTA].red);
    settings.setValue("Palette/Magenta.green", lpal->colors[Colors::PMAGENTA].green);
    settings.setValue("Palette/Magenta.blue", lpal->colors[Colors::PMAGENTA].blue);
    settings.setValue("Palette/Yellow.red", lpal->colors[Colors::PYELLOW].red);
    settings.setValue("Palette/Yellow.green", lpal->colors[Colors::PYELLOW].green);
    settings.setValue("Palette/Yellow.blue", lpal->colors[Colors::PYELLOW].blue);
    settings.setValue("Palette/White.red", lpal->colors[Colors::PWHITE].red);
    settings.setValue("Palette/White.green", lpal->colors[Colors::PWHITE].green);
    settings.setValue("Palette/White.blue", lpal->colors[Colors::PWHITE].blue);
    settings.setValue("Palette/Black.red", lpal->colors[Colors::PBLACK].red);
    settings.setValue("Palette/Black.green", lpal->colors[Colors::PBLACK].green);
    settings.setValue("Palette/Black.blue", lpal->colors[Colors::PBLACK].blue);

    settings.setValue("User Colors/User1.red", lpal->colors[Colors::PCUSTOM1].red);
    settings.setValue("User Colors/User1.green", lpal->colors[Colors::PCUSTOM1].green);
    settings.setValue("User Colors/User1.blue", lpal->colors[Colors::PCUSTOM1].blue);
    settings.setValue("User Colors/User2.red", lpal->colors[Colors::PCUSTOM2].red);
    settings.setValue("User Colors/User2.green", lpal->colors[Colors::PCUSTOM2].green);
    settings.setValue("User Colors/User2.blue", lpal->colors[Colors::PCUSTOM2].blue);
    settings.setValue("User Colors/User3.red", lpal->colors[Colors::PCUSTOM3].red);
    settings.setValue("User Colors/User3.green", lpal->colors[Colors::PCUSTOM3].green);
    settings.setValue("User Colors/User3.blue", lpal->colors[Colors::PCUSTOM3].blue);
    settings.setValue("User Colors/User4.red", lpal->colors[Colors::PCUSTOM4].red);
    settings.setValue("User Colors/User4.green", lpal->colors[Colors::PCUSTOM4].green);
    settings.setValue("User Colors/User4.blue", lpal->colors[Colors::PCUSTOM4].blue);
    settings.setValue("User Colors/User5.red", lpal->colors[Colors::PCUSTOM5].red);
    settings.setValue("User Colors/User5.green", lpal->colors[Colors::PCUSTOM5].green);
    settings.setValue("User Colors/User5.blue", lpal->colors[Colors::PCUSTOM5].blue);
    settings.setValue("User Colors/User6.red", lpal->colors[Colors::PCUSTOM6].red);
    settings.setValue("User Colors/User6.green", lpal->colors[Colors::PCUSTOM6].green);
    settings.setValue("User Colors/User6.blue", lpal->colors[Colors::PCUSTOM6].blue);
    settings.setValue("View Parameters/BackgroundColorRed", BackgroundColor.red);
    settings.setValue("View Parameters/BackgroundColorGreen", BackgroundColor.green);
    settings.setValue("View Parameters/BackgroundColorBlue", BackgroundColor.blue);
    settings.setValue("View Parameters/LabelPicks", LabelPicks);
    settings.setValue("View Parameters/LabelToggle", LabelToggle);
    settings.setValue("View Parameters/LabelStyle", ATOMLABEL::defaultStyle());
    settings.setValue("View Parameters/LabelSize", ATOMLABEL::defaultSize());
    settings.setValue("View Parameters/LabelColorRed", ATOMLABEL::defaultRed());
    settings.setValue("View Parameters/LabelColorGreen", ATOMLABEL::defaultGreen());
    settings.setValue("View Parameters/LabelColorBlue", ATOMLABEL::defaultBlue());
}

void Application::WriteSecStr()
{
    QSettings settings;
    settings.setValue("SecStr Colors/Helix", Colors::sec_colors[Colors::HELIX]);
    settings.setValue("SecStr Colors/Sheet", Colors::sec_colors[Colors::SHEET]);
    settings.setValue("SecStr Colors/Coil", Colors::sec_colors[Colors::COIL]);
}

void Application::WriteAtomTypes()
{
    QSettings settings;
    settings.setValue("Atom Colors/NoEntries", (unsigned int)Colors::atomcolors.size());
    for (unsigned int i = 0; i < Colors::atomcolors.size(); i++)
    {
        std::string b = format("Color%d", i);
        settings.setValue(("Atom Colors/" + b).c_str(), Colors::findColorNumber(Colors::atomcolors[i]));
        b = format("Name%d", i);
        settings.setValue(("Atom Colors/" + b).c_str(), Colors::atomnames[i].c_str());
    }
}

void Application::WriteBValues()
{
    QSettings settings;
    settings.setValue("BValue Colors/Color1", Colors::BValueColors[0]);
    settings.setValue("BValue Colors/Range1", Colors::BValueRanges[0]);
    settings.setValue("BValue Colors/Color2", Colors::BValueColors[1]);
    settings.setValue("BValue Colors/Range2", Colors::BValueRanges[1]);
    settings.setValue("BValue Colors/Color3", Colors::BValueColors[2]);
    settings.setValue("BValue Colors/Range3", Colors::BValueRanges[2]);
    settings.setValue("BValue Colors/Color4", Colors::BValueColors[3]);
    settings.setValue("BValue Colors/Range4", Colors::BValueRanges[3]);
    settings.setValue("BValue Colors/Color5", Colors::BValueColors[4]);
    settings.setValue("BValue Colors/Range5", Colors::BValueRanges[4]);
    settings.setValue("BValue Colors/Color6", Colors::BValueColors[5]);
    settings.setValue("BValue Colors/Range6", Colors::BValueRanges[5]);
    settings.setValue("BValue Colors/Color7", Colors::BValueColors[6]);
    settings.setValue("BValue Colors/Range7", Colors::BValueRanges[6]);
    settings.setValue("BValue Colors/Color8", Colors::BValueColors[7]);
    settings.setValue("BValue Colors/Range8", Colors::BValueRanges[7]);
    settings.setValue("BValue Colors/Color9", Colors::BValueColors[8]);
    settings.setValue("BValue Colors/Range9", Colors::BValueRanges[8]);
    settings.setValue("BValue Colors/Color10", Colors::BValueColors[9]);
    settings.setValue("BValue Colors/Range10", Colors::BValueRanges[9]);
}


const std::string&Application::GetBinDirectory()
{
    MIbin = findMIFitBin();
    return MIbin;
}

void Application::SetMolimageHome(bool reset)
{
    if (reset || MolimageHome.empty())
    {
        MolimageHome = findMIFitBin();
        if (strncasecmp(MolimageHome.c_str(), "/bin", 4) == 0
            || strncasecmp(MolimageHome.c_str(), "\\bin", 4) == 0)
        {
            MolimageHome.erase(MolimageHome.size()-4);
        }
#ifdef __APPLE__
        MolimageHome += "/../Resources";
#endif
    }
}

void Application::SetShelxHome(bool reset)
{
    QSettings settings;
    ShelxHome = settings.value("SHELXHOME").toString().toStdString();
    if (reset || ShelxHome.empty() || !QFileInfo(ShelxHome.c_str()).exists())
    {
        QString str = QFileDialog::getExistingDirectory(0, "Find the Shelx home directory");
        if (!str.isEmpty())
            ShelxHome = str.toStdString();
    }
}

void Application::SetEnv()
{
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

void Application::SetCrystalData(bool reset)
{
    if (reset || CrystalData.empty())
    {
        CrystalData = QSettings().value("CRYSTALDATA").toString().toStdString();
        if (CrystalData.empty())
        {
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


std::string Application::GetCrystalData(const char *crystal, const char *key)
{
    std::string data;
    FILE *fp;
    char buff[3200];
    char filename[BUFSIZ];

    if (crystal == NULL)
    {
        return data;
    }
    if (key == NULL)
    {
        return data;
    }
    if (crystal[0] == '\0')
    {
        return data;
    }
    if (key[0] == '\0')
    {
        return data;
    }

    // first look in current directory for an override
    // make sure it is not a directory
    strcpy(filename, crystal);
    if (filename[0] == '~')
    {
        strcpy(filename, TildeExpand(filename));
    }
    if ((fp = fopen(crystal, "r")) == NULL)
    {
        strcpy(filename, CrystalData.c_str());
        strcat(filename, "/");
        strcat(filename, crystal);
        if ((fp = fopen(filename, "r")) == NULL)
        {
            Logger::message("Error: can't find file for crystal: %s", crystal);
            return data;
        }
    }
    while (fgets(buff, 3200, fp) != NULL)
    {
        if (strncmp(key, buff, strlen(key)) == 0)
        {
            /* find position of first space */
            data = buff;
            data = MIAfterFirst(data, ' ');
            MIStringTrim(data, false);
            MIStringTrim(data, true);
            break;
        }
    }
    fclose(fp);
    return data;
}

bool Application::GetCrystalCell(const char *crystal, float &a, float &b, float &c,
                                 float &alpha, float &beta, float &gamma)
{
    if (crystal == NULL)
    {
        return false;
    }
    if (crystal[0] == '\0')
    {
        return false;
    }

    std::string line = GetCrystalData(crystal, "cell");
    if (sscanf(line.c_str(), "%f%f%f%f%f%f", &a, &b, &c, &alpha, &beta, &gamma) == 6)
    {
        return true;
    }
    return false;
}

int Application::GetCrystalNCSSymmops(const char *crystal, float ncrsymm[MISymmop::MAXNCRSYMMOPS][12])
{
    int nops = 0;
    std::string key;
    std::string buf;

    for (unsigned int i = 1; i <= MISymmop::MAXNCRSYMMOPS; i++)
    {
        key = format("ncrsymm%d", i);
        buf = GetCrystalData(crystal, key.c_str());
        if (buf.size() > 0)
        {
            sscanf(buf.c_str(), "%f%f%f%f%f%f%f%f%f%f%f%f", &ncrsymm[nops][0], &ncrsymm[nops][1], &ncrsymm[nops][2],
                   &ncrsymm[nops][3], &ncrsymm[nops][4], &ncrsymm[nops][5], &ncrsymm[nops][6], &ncrsymm[nops][7],
                   &ncrsymm[nops][8], &ncrsymm[nops][9], &ncrsymm[nops][10], &ncrsymm[nops][11]);
            nops++;
        }
        else
        {
            break;
        }
    }
    return (nops);
}

const std::string&Application::GetShelxHome()
{
    if (ShelxHome.empty())
    {
        SetShelxHome();
    }
    return ShelxHome;
}

bool Application::SetDictionary(const std::string &dict)
{
    bool valueChanged = false;
    std::string newDict = dict;
    if (newDict.empty())
    {
        newDict = MolimageHome;
#ifndef _WIN32
        newDict += "/data/dict.noh.pdb";
#else
        newDict += "\\data\\dict.noh.pdb";
#endif
    }
    if (XFitDict != newDict)
    {
        XFitDict = newDict;
        valueChanged = true;
    }
    return valueChanged;
}

std::string Application::getDictionary()
{
    return XFitDict;
}


void Application::toggleStereo()
{
    QSettings settings;
    bool stereo = settings.value("View Parameters/stereo", false).toBool();
    settings.setValue("View Parameters/stereo", !stereo);
}

long Application::GetPid()
{
    return getpid();
}


//FIXME: for Qt version, try to get a hardware stereo buffer at startup.
// If we can get one, do there's no saying we have to *use* the stereo
// planes.  This will avoid the silly "re-open to enable stereo"
// requirement.
void Application::toggleHardwareStereo()
{
    QSettings settings;
    bool hardwareStereo = settings.value("View Parameters/hardwareStereo", false).toBool();
    settings.setValue("View Parameters/hardwareStereo", !hardwareStereo);
    //  QMessageBox::warning(this, "Hardware Stereo Notice",
    //    "Changing the hardware stereo setting does not affect\n"
    //    "open document or dictionary windows. Re-open documents\n"
    //    "to enable the new setting.");
}



MIBusyManager*MIBusyManager::_instance = 0;
#define QPROGRESS_MAX 100
//#define QPROGRESS_MAX 0

MIBusyManager*MIBusyManager::instance()
{
    if (!_instance)
        _instance = new MIBusyManager();
    return _instance;
}


MIBusyManager::MIBusyManager()
{
    m_cursor_number = 0;
    m_busy = 0;
    MustAbortNow = 0;
    PROGRESS = new QProgressDialog("MIFit is busy", "Abort", 0, QPROGRESS_MAX);
    PROGRESS->setWindowModality(Qt::ApplicationModal);
}

bool MIBusyManager::Busy() const
{
    return m_busy > 0;
}

void MIBusyManager::SetBusy(bool busy)
{
    if (busy)
    {
        if (m_busy == 0)
        {
            MustAbortNow = false;
            PROG_VALUE = 0;
            PROG_DIRECTION = 1;
            PROGRESS->reset();
            PROGRESS->setValue(0);
            PROGRESS->show();
        }
        m_busy++;
    }
    else
    {
        if (m_busy > 0)
        {
            m_busy--;
        }
        else
        {
#ifdef DEBUG
            //FIXME Log("Warning: End of busy without corresponding start busy");
#endif
        }
    }
    if (m_busy == 0)
    {
        PROGRESS->reset();
        PROGRESS->setLabelText("MIFit is busy");
    }
}

bool MIBusyManager::CheckForAbort()
{
    if (PROG_VALUE == 0)
    {
        PROG_DIRECTION = 1;
    }
    else if (PROG_VALUE + 1 == QPROGRESS_MAX)
    {
        PROG_DIRECTION = -1;
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

void MIBusyManager::ForceAbort()
{
    MustAbortNow = true;
}

void MIBusyManager::SetWaitCursor()
{
    switch (m_cursor_number)
    {
    case 0: SetCursor(imhWait1);
        break;
    case 1: SetCursor(imhWait2);
        break;
    case 2: SetCursor(imhWait3);
        break;
    case 3: SetCursor(imhWait4);
        break;
    case 4: SetCursor(imhWait5);
        break;
    case 5: SetCursor(imhWait6);
        break;
    }
    m_cursor_number = (m_cursor_number+1)%6;
}

void MIBusyManager::StartWaitCursor()
{
    SetWaitCursor();
}

void MIBusyManager::StopWaitCursor()
{
    SetCursor(-1);
}

void MIBusyManager::SetCursor(int id)
{
    MIMainWindow::instance()->SetCursor(id);
}


GeomRefiner *MIFitGeomRefiner()
{
    return Application::instance()->GetGeomRefiner();
}

chemlib::MIMolDictionary *MIFitDictionary()
{
    return &MIFitGeomRefiner()->dict;
}

QString Application::latestFileBrowseDirectory(const QString &path)
{
    const QString LATEST_FILE_BROWSE_DIRECTORY("latestFileBrowseDirectory");
    QSettings settings;
    if (path.isEmpty())
    {
        QString settingsPath = settings.value(LATEST_FILE_BROWSE_DIRECTORY).toString();
        if (QFile::exists(settingsPath))
            return settingsPath;
    }
    else
        settings.setValue(LATEST_FILE_BROWSE_DIRECTORY, path);
    return path;
}

QString Application::getOpenFileName(QWidget *parent, const QString &caption, const QString &filter)
{
    QString dir = Application::instance()->latestFileBrowseDirectory("");
    QString file = QFileDialog::getOpenFileName(parent, caption, dir, filter);
    if (!file.isEmpty())
    {
        QFileInfo fileInfo(file);
        Application::instance()->latestFileBrowseDirectory(fileInfo.absolutePath());
    }
    return file;
}

QString Application::getExistingDirectory(QWidget *parent, const QString &caption, const QString &filter)
{
    QString dir = Application::instance()->latestFileBrowseDirectory("");
    dir = QFileDialog::getExistingDirectory(parent, caption, dir);
    if (!dir.isEmpty())
    {
        Application::instance()->latestFileBrowseDirectory(dir);
    }
    return dir;
}
