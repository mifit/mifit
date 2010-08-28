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
#include "ContourOptionsWidget.h"
#include "RefinementOptionsDialog.h"
#include "PhaseFileLoadDialog.h"
#include "BValueColors.h"
#include "LSQFitDialog.h"
#include "CustomJobDialog.h"
#include "GenericDataDialog.h"

//
// File
//
MIFileDialog::MIFileDialog(QWidget *parent, const std::string &message,
                           const std::string &deftDir,
                           const std::string &deftFile,
                           const std::string &filter,
                           unsigned int mode)
    : _qparent(parent),
      _name(message),
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


bool MIFileDialog::GetResults(Data &data)
{
    std::string path = _deftDir;
    if (_deftFile.size())
    {
        path = path+ "/" + _deftFile;
    }
    if (!data.path.empty())
    {
        path = data.path;
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

    std::vector<std::string> &pathlist = data.pathlist;
    pathlist.clear();

    path = Application::instance()->latestFileBrowseDirectory(path.c_str()).toStdString();

    if (_mode == MI_SAVE_MODE)
    {
        response = QFileDialog::getSaveFileName(_qparent, _name.c_str(), path.c_str(), filter, &selectedFilter);
        if (response.isEmpty())
            return false;
        data.path = response.toStdString();
        pathlist.push_back(response.toStdString());
    }
    else if (_mode==MI_OPEN_MODE)
    {
        response = QFileDialog::getOpenFileName(_qparent, _name.c_str(), path.c_str(), filter, &selectedFilter);
        if (response.isEmpty())
            return false;
        data.path = response.toStdString();
        pathlist.push_back(response.toStdString());
    }

    QFileInfo fileInfo(data.path.c_str());
    Application::instance()->latestFileBrowseDirectory(fileInfo.absolutePath());

    // set selected filter index
    std::string selFilter = selectedFilter.toStdString();
    for (size_t i = 0; i<filterList.size(); ++i)
    {
        if (filterList[i] == selFilter)
        {
            data.filterIndex = i;
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

    MIFileDialog::Data data;
    std::string filename;
    if (fileDialog.GetResults(data))
    {
        filename = data.path;
    }
    QFileInfo fileInfo(filename.c_str());
    Application::instance()->latestFileBrowseDirectory(fileInfo.absolutePath());
    return filename;
}

