#include "Cfiles.h"

CArchive::CArchive(const char *pathname, int m)
{
    // store or load?
    mode = m;
    if (m == CArchive::store)
    {
        file = fopen(pathname, "w");
    }
    else
    {
        file = fopen(pathname, "r");
    }
}

bool CArchive::IsOpened() const
{
    return file != 0;
}

CArchive::~CArchive()
{
    if (file)
    {
        if (IsStoring())
        {
            fflush(file);
        }
        fclose(file);
    }
}

void CArchive::Flush()
{
    if (file)
        fflush(file);
}

long CArchive::Read(void *buffer, long nread)
{
    if (!file)
        return 0L;
    return fread(buffer, 1, nread, file);
}

void CArchive::Write(const void *InBuf, long nInBuf)
{
    if (!file)
        return;
    fwrite(InBuf, 1, nInBuf, file);
}

void CArchive::Write(std::string &f)
{
    Write(f.c_str(), f.size());
}

