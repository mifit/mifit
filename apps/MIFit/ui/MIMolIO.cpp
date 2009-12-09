#include "MIMolIO.h"
#include "macafxwin.h"
#include <QFileDialog>
#include <QInputDialog>
#include <QStringList>

using namespace chemlib;

MIMolIO::MIMolIO()
    : MIMolIOBase()
{
}

bool MIMolIO::Write(MIMolInfo &mi, const std::string &fname, int writerIndex) const
{
    QString filename = fname.c_str();
    if (filename.isEmpty())
    {
        QString filter;
        WriterList::const_iterator iter = writers.begin();
        for (iter = writers.begin(); iter != writers.end(); ++iter)
        {
            Writer *writer = *iter;
            filter += writer->getDescription().c_str();
            filter += ";;";
        }
        filter += "All files (*.*)";
        filename = QFileDialog::getSaveFileName(0, "Choose a file to save to",
                                  "", filter);
    }
    if (filename.isEmpty())
    {
        return false;
    }

    if (writerIndex == -1)
    {
        writerIndex = getWriterIndex(filename.toStdString());
    }
    if (writerIndex < 0)
    {
        QStringList writerDescriptions;
        for (unsigned int i = 0; i < GetWriterCount(); ++i)
        {
            writerDescriptions += writers[i]->getDescription().c_str();
        }
        QString selectedWriter = QInputDialog::getItem(0, "Output format", "Select Output Format",
                                             writerDescriptions);
        writerIndex = writerDescriptions.indexOf(selectedWriter);
    }
    if (writerIndex < 0)
    {
        return false;
    }

    return MIMolIOBase::Write(mi, filename.toStdString(), writerIndex);
}

bool MIMolIO::Read(MIMolInfo &mi, const std::string &fname, int readerIndex) const
{
    QString filename = fname.c_str();
    if (filename.size() == 0)
    {
        QString filter;
        ReaderList::const_iterator iter = readers.begin();
        for (iter = readers.begin(); iter != readers.end(); ++iter)
        {
            Reader *reader = *iter;
            filter += reader->getDescription().c_str();
            filter += ";;";
        }
        filter += "All files (*.*)";
        filename = QFileDialog::getOpenFileName(0, "Choose a file to read",
                                  "", filter);
    }
    if (filename.isEmpty())
    {
        return false;
    }

    if (readerIndex == -1)
    {
        readerIndex = getReaderIndex(filename.toStdString());
    }
    if (readerIndex < 0)
    {
        QStringList readerDescriptions;
        for (unsigned int i = 0; i < GetReaderCount(); ++i)
        {
            readerDescriptions += readers[i]->getDescription().c_str();
        }
        QString selectedReader = QInputDialog::getItem(0, "Input format", "Select Input Format",
                                             readerDescriptions);
        readerIndex = readerDescriptions.indexOf(selectedReader);
    }
    if (readerIndex < 0)
    {
        return false;
    }

    return MIMolIOBase::Read(mi, filename.toStdString(), readerIndex);
}

