#ifndef CFILES_H_
#define CFILES_H_

#include <string>
#include <stdio.h>

//@{
//  A class to wrap writing to a file.
//  Modeled afer the MFC class of the same name.
//@}
class CArchive
{
private:
    FILE *file;
    int mode;

public:
    int filetype;
    enum Mode { store = 0, load = 1 };
    bool IsLoading() const
    {
        return (mode & CArchive::load) != 0;
    }

    bool IsStoring() const
    {
        return (mode & CArchive::load) == 0;
    }

    CArchive(const char *pathname, int store_or_load);
    ~CArchive();
    long Read(void*, long);
    void Write(const void*, long);
    void Write(std::string&);
    void Flush();
    void Rewind()
    {
        if (file)
            rewind(file);
    }

    bool IsOpened() const;
};

#endif //CFILES_H_
