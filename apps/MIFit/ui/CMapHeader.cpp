
#include "CMapHeader.h"
#include "Application.h"

#include <QFileInfo>


using namespace chemlib;

CMapHeader::CMapHeader()
    : CMapHeaderBase()
{
}

bool CMapHeader::LoadCrystal(const char *crystal)
{
    std::string s;
    int n;
    char buf[1000];
    if (Application::instance()->GetCrystalCell(crystal, a, b, c, alpha, beta, gamma))
    {
        // reset the rest
        spgpno = 0;
        spgpname.clear();
        nsym = 0;
        nNCRSymmops = 0;
    }
    else
    {
        return false;
    }
    orthog(a, b, c, alpha, beta, gamma, ctof);
    uinv(ctof, ftoc);
    s = Application::instance()->GetCrystalData(crystal, "spgp");
    if (s.size() > 0)
    {
        int ret = sscanf(s.c_str(), "%s %d", buf, &n);
        if (ret >= 1)
        {
            spgpname = buf;
        }
        if (ret >= 2)
        {
            spgpno = n;
        }
        SymInfoString = s;
    }
    s = Application::instance()->GetCrystalData(crystal, "name");
    if (s.size() > 0)
    {
        title = "name ";
        title += s;
    }
    else
    {
        title = std::string();
    }
    s = Application::instance()->GetCrystalData(crystal, "symm");
    if (s.size() > 0)
    {
        char opstring[MISymmop::MAXSYMMOPS][MISymmop::MAXSTRING];
        nsym = scan_symmops(s.c_str(), symops, opstring);
        if (nsym > 0)
        {
            SymopsString.clear();
            for (int i = 0; i < nsym; i++)
            {
                SymopsString.push_back(std::string(opstring[i]));
            }
        }
    }
    nNCRSymmops = Application::instance()->GetCrystalNCSSymmops(crystal, NCRSymmops);

    QFileInfo qfile(crystal);
    crystal_name = qfile.fileName().toStdString();
    return true;
}

bool CMapHeader::SaveCrystal(const std::string &crystal)
{
    std::string d, f, s;
    std::vector<std::string> oldfile;
    char buf[1024];
    int i;

    int oldfile_exists = false;
    if (crystal.c_str()[0] != '.' && !Application::instance()->CrystalData.empty())
    {
        d = Application::instance()->CrystalData;
        d += "/";
    }
    else
    {
        d = "";
    }
    f = d;
    f += crystal.c_str();
    FILE *fpold = fopen(f.c_str(), "r");
    /* if an old file already exists, read it in */
    if (fpold)
    {
        oldfile_exists = true;
        while (fgets(buf, sizeof buf, fpold) != NULL)
        {
            oldfile.push_back(std::string(buf));
        }
        fclose(fpold);
    }
    FILE *fpnew = fopen(f.c_str(), "w");
    if (!fpnew)
    {
        Logger::message("Error: couldn't save crystal!  Could be a write permission error, or directory %s may not exist.", d.c_str());
        return false;
    }
    fprintf(fpnew, "%s", title.c_str());
    if (title[title.size()-1] != '\n')
    {
        fprintf(fpnew, "\n");
    }
    fprintf(fpnew, "cell %f %f %f %f %f %f\n", a, b, c, alpha, beta, gamma);
    fprintf(fpnew, "spgp %s %d\n", spgpname.c_str(), (int)spgpno);
    fprintf(fpnew, "symm ");
    if (SymopsString.size() != 0)
    {
        for (i = 0; i < nsym; i++)
        {
            fprintf(fpnew, "%s", MIToLower(SymopsString[i]).c_str());
            if (i < nsym-1)
            {
                fprintf(fpnew, "; ");
            }
        }
    }
    fprintf(fpnew, ".\n");
    // ncrsymmetry
    if (nNCRSymmops > 0)
    {
        for (i = 0; i < nNCRSymmops; i++)
        {
            fprintf(fpnew, "ncrsymm%d %8.4f %8.4f %8.4f %8.4f %8.4f %8.4f %8.4f %8.4f %8.4f %8.4f %8.4f %8.4f\n", i+1,
                    NCRSymmops[i][0], NCRSymmops[i][1], NCRSymmops[i][2], NCRSymmops[i][3],
                    NCRSymmops[i][4], NCRSymmops[i][5], NCRSymmops[i][6], NCRSymmops[i][7],
                    NCRSymmops[i][8], NCRSymmops[i][9], NCRSymmops[i][10], NCRSymmops[i][11]);
        }
    }
    if (oldfile_exists)
    {
        for (i = 0; i < (int)oldfile.size(); i++)
        {
            s = oldfile[i];
            if (strncmp(s.c_str(), "symm", 4)==0
                || strncmp(s.c_str(), "cell", 4)==0
                || strncmp(s.c_str(), "spgp", 4)==0
                || strncmp(s.c_str(), "name", 4)==0
                || strncmp(s.c_str(), "ncrs", 4)==0)
            {
            }
            else
            {
                fprintf(fpnew, "%s", s.c_str());
            }
        }
    }
    fclose(fpnew);
    return true;
}

