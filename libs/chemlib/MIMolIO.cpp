#include <ui/Logger.h>
#include "MIMolIOBase.h"
#include "mmCIF.h"
#include "PDB.h"
#include "molfile.h"
#include "SMILES.h"
#include <util/system.h>

using namespace std;

namespace chemlib
{

MIMolIOBase::MIMolIOBase()
{
    registerReader(new mmCIF());
    registerWriter(new mmCIF());

    registerReader(new PDB());
    registerWriter(new PDB());

    registerReader(new molfile());
    registerWriter(new molfile());

    registerReader(new SMILES());
}

MIMolIOBase::~MIMolIOBase()
{
    ReaderList::iterator reader = readers.begin();
    while (reader != readers.end())
    {
        delete *reader;
        ++reader;
    }
    WriterList::iterator writer = writers.begin();
    while (writer != writers.end())
    {
        delete *writer;
        ++writer;
    }
}

void MIMolIOBase::registerReader(Reader *reader)
{
    readers.push_back(reader);
}

void MIMolIOBase::registerWriter(Writer *writer)
{
    writers.push_back(writer);
}

int MIMolIOBase::getReaderIndex(const std::string &fname) const
{
    std::string::size_type pos = fname.find_last_of(".");
    if (pos == std::string::npos || pos == fname.size())
    {
        return -1;
    }
    std::string extension = fname.substr(pos+1);

    for (unsigned int i = 0; i < readers.size(); ++i)
    {
        std::string tmp = "*." + extension;
        if (strcasecmp(tmp.c_str(),
                       readers[i]->getExtension().c_str()) == 0)
        {
            return i;
        }
    }
    return -1;
}

int MIMolIOBase::getWriterIndex(const std::string &fname) const
{
    std::string::size_type pos = fname.find_last_of(".");
    if (pos == std::string::npos || pos == fname.size())
    {
        return -1;
    }
    std::string extension = fname.substr(pos+1);

    for (unsigned int i = 0; i < writers.size(); ++i)
    {
        std::string tmp = "*." + extension;
        if (strcasecmp(tmp.c_str(),
                       writers[i]->getExtension().c_str()) == 0)
        {
            return i;
        }
    }
    return -1;
}

std::string MIMolIOBase::GetReaderDescription(unsigned int idx) const
{
    if (idx > readers.size())
    {
        return "";
    }
    return readers[idx]->getDescription();
}

std::string MIMolIOBase::GetWriterDescription(unsigned int idx) const
{
    if (idx > writers.size())
    {
        return "";
    }
    return writers[idx]->getDescription();
}

bool MIMolIOBase::Read(MIMolInfo &mi, const std::string &filename, int readerIndex) const
{
    if (filename.size() == 0 || readerIndex >= (int)readers.size())
    {
        return false;
    }

    if (readerIndex < 0)
    {
        readerIndex = getReaderIndex(filename);
    }
    if (readerIndex < 0)
    {
        return false;
    }

    FILE *input = fopen(filename.c_str(), "r");
    if (input == NULL)
    {
        Logger::message("Error opening file for reading. Permissions problem?");
        return false;
    }

    // Read input
    Reader *fio = readers[readerIndex];
    if (!fio->Read(input, mi))
    {
        Logger::message("Error reading file.");
        fclose(input);
        return false;
    }
    fclose(input);
    return true;
}

bool MIMolIOBase::Write(MIMolInfo &mi, const std::string &filename, int writerIndex) const
{
    if (filename.size() == 0 || writerIndex >= (int)writers.size())
    {
        return false;
    }

    if (writerIndex < 0)
    {
        writerIndex = getWriterIndex(filename);
    }
    if (writerIndex < 0)
    {
        return false;
    }


    // Setup the output file
    FILE *output = fopen(filename.c_str(), "w");
    if (output == NULL)
    {
        Logger::message("Error opening file for writing. Permissions problem?");
        return false;
    }

    // Write output
    Writer *fio = writers[writerIndex];
    if (!fio->Write(output, mi))
    {
        Logger::message("Error writing file. Out of disk space?");
        fclose(output);
        return false;
    }
    fclose(output);
    return true;
}

} // namespace chemlib


