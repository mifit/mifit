#define NOMINMAX
#include <cmath>
#include <climits>
#include <cstring>

#include <nongui/nonguilib.h>
#include <math/mathlib.h>
#include <chemlib/chemlib.h>

#include "Molecule.h"
#include "seqtypes.h"
#include "RESIDUE.h"
#include <chemlib/Monomer.h>

using namespace chemlib;


int Molecule::SeqMax()
{
    //int seqmax = 0;
    ResidueListIterator ri = residuesBegin();
    ResidueListIterator res = ri;
    for (; ri != residuesEnd(); ++ri)
        res = ri;
    if (res == residuesEnd())
        return 0;
    return std::max((size_t)res->seqpos(), alt_seq.length());
}

char Molecule::GetSeq(int index)
{
    if (index >= (int)alt_seq.length())
    {
        return ' ';
    }
    else
    {
        return alt_seq.c_str()[index];
    }
}

void Molecule::SetSeq(char /* aa */, int index)
{
    if (index < (int)alt_seq.length())
    {
        alt_seq = MIToUpper(alt_seq);
        //if(islower(aa)) aa = toupper(aa);
        //alt_seq.SetAt(index, aa);
    }
}

std::string Molecule::SeqString()
{
    return alt_seq;
}

void Molecule::SetSequence(std::string s)
{
    int i /*, iseq=0*/;
    char c;
    alt_seq = "";
    for (i = 0; i < (int)s.length(); i++)
    {
        c = s.c_str()[i];
        if (isalpha(c) || c == '.')
        {
            if (islower(c))
            {
                c = toupper(c);
            }
            alt_seq += c;
        }
    }
}

void Molecule::ReadSequence(std::string path, int /* type */, int /* skiplines */)
{
    FILE *fp = fopen(path.c_str(), "r");
    if (!fp)
    {
        Logger::message("Cannot open sequence file");
    }
    std::string file;
    char buf[1024];
    while (fgets(buf, sizeof buf, fp) != NULL)
    {
        if (!(buf[0] == '#' || buf[0] == '<' || buf[0] == ';' || buf[0] == '>'))
        {
            file += buf;
        }
    }
    SetSequence(file);
    fclose(fp);
}

int Molecule::WriteSequence(std::string path, int type)
{
    FILE *fp = fopen(path.c_str(), "w");
    char c;
    if (!fp)
    {
        Logger::message("Cannot open sequence file");
        return 0;
    }
    fprintf(fp, ">%s\n\n", compound.c_str());
    //std::string s = SeqString();
    int l = alt_seq.length();
    int i;
    if (type == SEQ_FORMAT_SIMPLE)
    {
        for (i = 0; i < l; i++)
        {
            c = alt_seq.c_str()[i];
            fprintf(fp, "%c", alt_seq.c_str()[i]);
            if (i > 0 && i%10 == 0 && i%100 != 0)
            {
                fprintf(fp, " ");
            }
            if (i > 0 && i%100 == 0)
            {
                fprintf(fp, "\n");
            }
        }
    }
    else if (type == SEQ_FORMAT_MODEL)
    {
        i = 0;
        ResidueListIterator res = residuesBegin();
        for (; res != residuesEnd(); ++res)
        {
            c = (char)res->name1();
            if (!(c == '?' || c == 'B' || c == 'J' || c == 'O' || c == 'U' || c == 'X' || c == 'Z') && isupper(c))
            {
                fprintf(fp, "%c", c);
                if (i > 0 && i%10 == 0 && i%100 != 0)
                {
                    fprintf(fp, " ");
                }
                if (i > 0 && i%100 == 0)
                {
                    fprintf(fp, "\n");
                }
                i++;
            }
        }
    }
    else if (type == SEQ_FORMAT_MSF)
    {
        int nout = 0;
        for (i = 0; i < l; i++)
        {
            c = toupper(alt_seq.c_str()[i]);
            if (!(c == '?' || c == 'B' || c == 'J' || c == 'O' || c == 'U' || c == 'X' || c == 'Z') && isalpha(c))
            {
                fprintf(fp, "%c", alt_seq.c_str()[i]);
                //if(nout>0 && nout%10==0 && nout%100!=0)fprintf(fp," ");
                if (nout > 0 && nout%80 == 0)
                {
                    fprintf(fp, "\n");
                }
                nout++;
            }
        }
    }
    fprintf(fp, "\n");
    fclose(fp);
    return 1;
}

void Molecule::InsertGap(ResidueListIterator gap_point)
{
    if (gap_point == NULL)
    {
        return;
    }
    ResidueListIterator res = residuesBegin();
    for (; res != residuesEnd(); ++res)
    {
        if (gap_point == res)
        {
            while (res != residuesEnd())
            {
                res->setSeqpos(res->seqpos() + 1);
                ++res;
            }
            return;
        }
    }
}

void Molecule::DeleteGap(ResidueListIterator gap_point)
{
    if (gap_point == NULL)
    {
        return;
    }
    ResidueListIterator prev = residuesEnd();
    ResidueListIterator res = residuesBegin();
    for (; res != residuesEnd(); ++res)
    {
        if (gap_point == res)
        {
            if (prev != residuesEnd())
            {
                if (prev->seqpos() < gap_point->seqpos()-1)
                {
                    while (res != residuesEnd())
                    {
                        res->setSeqpos(res->seqpos()-1);
                        ++res;
                    }
                }
            }
            else if (res == residuesBegin() && res->seqpos() > 0)
            {
                while (res != residuesEnd())
                {
                    res->setSeqpos(res->seqpos()-1);
                    ++res;
                }
            }
            return;
        }
        prev = res;
    }
}

void Molecule::InsertLowerGap(ResidueListIterator gap_point)
{
    if (gap_point == residuesEnd())
    {
        return;
    }
    int p = gap_point->seqpos();
    int l = alt_seq.length();
    if (p >= l)
    {
        return;
    }
    alt_seq.insert(p, ".");
}

void Molecule::DeleteLowerGap(ResidueListIterator gap_point)
{
    if (gap_point == residuesEnd())
    {
        return;
    }
    int p = gap_point->seqpos();
    int l = alt_seq.length();
    if (p >= l)
    {
        return;
    }
    if (alt_seq.c_str()[p] != '.')
    {
        return;
    }
    alt_seq.erase(p, 1);
}

int Molecule::SequenceIdentities()
{
    int n = 0;
    int l = alt_seq.length();
    if (l > 0)
    {
        ResidueListIterator res = residuesBegin();
        for (; res != residuesEnd(); ++res)
        {
            if (res->seqpos() < (unsigned int)l)
            {
                if (alt_seq[res->seqpos()] == res->name1() && isalpha(res->name1()))
                {
                    n++;
                }
            }
        }
    }
    return n;
}

