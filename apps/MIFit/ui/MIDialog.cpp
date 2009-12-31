#include <QMessageBox>
#include <QInputDialog>
#include <QColorDialog>
#include <QTextEdit>
#include <QDialogButtonBox>
#include <QLabel>
#include <QMdiArea>
#include <QFontDialog>
#include <QFileDialog>

#include <map/maplib.h>

#include "ui/MIMainWindow.h"

#include "MIDialog.h"

// local dialog class definitions
#include "MIColorPickerDlg.h"
#include "SelectCrystal.h"
#include "SmilesDialog.h"
#include "ContourOptions.h"
#include "RefinementOptionsDialog.h"
#include "PhaseFileLoadDialog.h"
#include "BValueColors.h"
#include "LSQFitDialog.h"
#include "CustomJobDialog.h"
#include "GenericDataDialog.h"

MIDialog::MIDialog(QWidget *parent, const std::string &name)
    : _qparent(parent),
      _name(name)
{
    // Note:
    //
    // we can't use MIMainWindow directly as a dialog parent, b/c if it is,
    // the active window changes to the dialog instead of the MIGLWidget,
    // (MIMainWindow gets subWindowActivated(0) signal), and that can bork
    // up handling of dialog results.  It appears to be safe to use the
    // MdiArea as the parent. FMH.

    if (!_qparent || _qparent == MIMainWindow::instance())
        _qparent = MIMainWindow::instance()->currentMIGLWidget();
    if (!_qparent)
        _qparent = MIMainWindow::instance()->getMdiArea();
}

MIDialog::~MIDialog()
{
}

void ValidateData(const MIData &data)
{
    std::map<std::string, MIDatum>::const_iterator i = data.begin();
    for (; i!= data.end(); ++i)
    {
        MIDatum datum = i->second;
        if (datum.radio!=UINT_MAX && datum.radio_count==UINT_MAX)
        {
#ifdef DEBUG
            Logger::message("Programmer error: radio set, but radio_count not set!");
#endif
        }
    }
}

bool MIDialog::GetResults(MIData &data)
{

    bool ret;

    ValidateData(data);
    ret = PromptForResults(data);

    return ret;
}



//
// File
//
MIFileDialog::MIFileDialog(QWidget *parent, const std::string &message,
                           const std::string &deftDir,
                           const std::string &deftFile,
                           const std::string &filter,
                           unsigned int mode)
    : MIDialog(parent, message),
      _deftDir(deftDir),
      _deftFile(deftFile),
      _filter(filter),
      _mode(mode)
{
}

static void stringSplit(std::string str,
                        const std::string &delim,
                        std::vector<std::string> &results)
{
    unsigned int cutAt;
    while ( (cutAt = str.find_first_of(delim)) != str.npos)
    {
        if (cutAt > 0)
        {
            results.push_back(str.substr(0, cutAt));
        }
        str = str.substr(cutAt+1);
    }
    if (str.length() > 0)
    {
        results.push_back(str);
    }
}


bool MIFileDialog::PromptForResults(MIData &data)
{
    std::string path = _deftDir;
    if (_deftFile.size())
    {
        path = path+ "/" + _deftFile;
    }
    if (data["path"].str.size() && data["path"].str != MIDatum::INVALID_STRING)
    {
        path = data["path"].str;
    }

    // build filter string/ vector;
    QString filter;
    std::vector<std::string> filterList;
    MIStringReplace(_filter, ",", " ");
    if (_filter.size())
    {
        std::vector<std::string> results;
        stringSplit(_filter, "|", results);
        for (unsigned int i = 0; i<results.size(); i += 2)
        {
            filterList.push_back(results[i]);
            if (i > 0)
                filter += ";;";
            filter += QString(results[i].c_str());
        }
    }

    QString selectedFilter;
    QString response;
    QStringList responseList;

    std::vector<std::string> &pathlist = data["pathList"].strList;
    pathlist.clear();

    path = Application::instance()->latestFileBrowseDirectory(path.c_str()).toStdString();

    if (_mode == MI_SAVE_MODE)
    {
        response = QFileDialog::getSaveFileName(_qparent, _name.c_str(), path.c_str(), filter, &selectedFilter);
        if (response.isEmpty())
            return false;
        data["path"].str = response.toStdString();
        pathlist.push_back(response.toStdString());
    }
    else if (_mode==MI_OPEN_MODE)
    {
        response = QFileDialog::getOpenFileName(_qparent, _name.c_str(), path.c_str(), filter, &selectedFilter);
        if (response.isEmpty())
            return false;
        data["path"].str = response.toStdString();
        pathlist.push_back(response.toStdString());
    }

    QFileInfo fileInfo(data["path"].str.c_str());
    Application::instance()->latestFileBrowseDirectory(fileInfo.absolutePath());

    // set selected filter index
    std::string selFilter = selectedFilter.toStdString();
    for (size_t i = 0; i<filterList.size(); ++i)
    {
        if (filterList[i] == selFilter)
        {
            data["filterIndex"].radio = i;
            break;
        }
    }

    return true;
}


//
// wraps MIFileDialog
//
std::string MIFileSelector(const std::string &title,
                           const std::string &defaultDirString,
                           const std::string &defaultFileNameString,
                           const std::string &defaultExtension,
                           const std::string &filter,
                           unsigned int flags,
                           QWidget *parent)
{

    std::string dir = Application::instance()->latestFileBrowseDirectory(defaultDirString.c_str()).toStdString();

    std::string filter2 = filter;
    if (defaultExtension.size() && !filter.size() )
    {
        filter2 = std::string("*.") + defaultExtension;
    }

    MIFileDialog fileDialog(parent, title, dir,
                            defaultFileNameString, filter2,
                            flags);

    MIData data;
    data["path"].str = "";
    data["pathList"].strList = std::vector<std::string>();

    std::string filename;
    if (fileDialog.GetResults(data))
    {
        filename = data["path"].str;
    }
    QFileInfo fileInfo(filename.c_str());
    Application::instance()->latestFileBrowseDirectory(fileInfo.absolutePath());
    return filename;
}

