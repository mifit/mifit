#include "OpenJobResults.h"

#include <util/utillib.h>
#include "ui/uilib.h"
#include "ui/MIMainWindow.h"
#include "ui/MIDialog.h"
#include "ui/GenericDataDialog.h"
#include <QFileInfo>
#include <QDateTime>
#include <QDir>
#include <QMessageBox>


#include <vector>
#include <algorithm>


namespace
{
bool later(const QPair<QString, QDateTime> &f1, const QPair<QString, QDateTime> &f2)
{
    return f1 > f2;
}

void sortFileNameByTime(QStringList &names)
{
    typedef QPair<QString, QDateTime> FileTime;
    QList<FileTime> fileTimes;
    foreach (const QString &name, names)
    {
        QFileInfo fn(name);
        QDateTime time = fn.lastModified();
        fileTimes += qMakePair(name, time);
    }
    qSort(fileTimes.begin(), fileTimes.end(), later);

    names.clear();
    foreach (const FileTime &fileTime, fileTimes)
    {
        names += fileTime.first;
    }
}
}


void OpenJobResults::prompt(const std::string &workdir, const std::string &jobName)
{
    QDir dir(workdir.c_str());
    if (!dir.exists())
    {
        return;
    }

    QStringList mlwFiles;
    QStringList pdbFiles;
    QStringList mtzFiles;

    QStringList fileList = dir.entryList(QDir::Files | QDir::Readable | QDir::NoDotAndDotDot);
    foreach (QString file, fileList)
    {
        QFileInfo fi(file);

        if (fi.suffix() == "mlw")
        {
            mlwFiles += file;
        }
        else if (fi.suffix() == "pdb")
        {
            pdbFiles += file;
        }
        else if (fi.suffix() == "mtz")
        {
            mtzFiles += file;
        }
    }

    sortFileNameByTime(mlwFiles);
    sortFileNameByTime(pdbFiles);
    sortFileNameByTime(mtzFiles);


    if (!(mlwFiles.size() > 0 || pdbFiles.size() > 0 || mtzFiles.size() > 0))
    {
        std::string s = format("Sorry, unable to open any results because\nno session, PDB, or MTZ files were found in job working directory:\n%s", workdir.c_str());
        QMessageBox::warning(MIMainWindow::instance(), "Open Results unsuccessful", s.c_str());
        return;
    }


    GenericDataDialog dlg;
    dlg.setWindowTitle(QString("Job %1 Finished").arg(jobName.c_str()));

    int fields = 0;
    int mlwField = -1;
    int mtzField = -1;
    int pdbField = -1;
    if (!mlwFiles.empty())
    {
        mlwFiles += "None";
        dlg.addComboField("Load session:", mlwFiles, 0);
        mlwField = fields;
        ++fields;
    }
    if (!mtzFiles.empty())
    {
        mtzFiles += "None";
        dlg.addComboField("Load data:", mtzFiles, 0);
        mtzField = fields;
        ++fields;
    }
    if (!pdbFiles.empty())
    {
        pdbFiles += "None";
        dlg.addComboField("Load PDB:", pdbFiles, 0);
        pdbField = fields;
        ++fields;
    }
    if (dlg.exec() != QDialog::Accepted)
        return;

    std::vector<std::string> files;
    bool newWin = true;
    if (mlwField != -1)
    {
        QString res = mlwFiles.at(dlg.value(mlwField).toInt());
        if (res != "None")
        {
            files.push_back(dir.absoluteFilePath(res).toStdString());
            newWin = false;
        }
    }
    if (pdbField != -1)
    {
        QString res = pdbFiles.at(dlg.value(pdbField).toInt());
        if (res != "None")
        {
            files.push_back(dir.absoluteFilePath(res).toStdString());
        }
    }
    if (mtzField != -1)
    {
        QString res = mtzFiles.at(dlg.value(mtzField).toInt());
        if (res != "None")
        {
            files.push_back(dir.absoluteFilePath(res).toStdString());
        }
    }
    MIMainWindow::instance()->OpenFiles(files, newWin);
}
