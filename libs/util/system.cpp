#include "system.h"
#include <cstdio>
#include <cctype>

std::string MIBeforeFirst(const std::string &s, char sep)
{
    std::string::size_type pos = s.find_first_of(sep);
    if (pos == std::string::npos)
    {
        return s;
    }
    if (pos == 0)
    {
        return "";
    }
    return s.substr(0, pos);
}

std::string MIAfterFirst(const std::string &s, char sep)
{
    std::string::size_type pos = s.find_first_of(sep);
    if (pos == std::string::npos)
    {
        return "";
    }
    if (pos+1 == s.size())
    {
        return "";
    }
    return s.substr(pos+1);
}

std::string MIBeforeLast(const std::string &s, char sep)
{
    std::string::size_type pos = s.find_last_of(sep);
    if (pos == std::string::npos || pos == 0)
    {
        return "";
    }
    return s.substr(0, pos);
}

std::string MIAfterLast(const std::string &s, char sep)
{
    std::string::size_type pos = s.find_last_of(sep);
    if (pos == std::string::npos)
    {
        return s;
    }
    if (pos+1 == s.size())
    {
        return "";
    }
    return s.substr(pos+1);
}

std::string MIToUpper(const std::string &s)
{
    std::string str = s;
    for (unsigned int i = 0; i < str.size(); ++i)
    {
        str[i] = toupper(str[i]);
    }
    return str;
}

std::string MIToLower(const std::string &s)
{
    std::string str = s;
    for (unsigned int i = 0; i < str.size(); ++i)
    {
        str[i] = tolower(str[i]);
    }
    return str;
}

// If 'from' matches 'to' or 'from' is empty,
// does not parse 's', returns std::string::npos
// Otherwise returns number of replacements done

std::string::size_type MIStringReplace(std::string &s,
                                       const std::string &from,
                                       const std::string &to)
{
    std::string::size_type cnt(std::string::npos);

    if (from != to && !from.empty())
    {
        std::string::size_type pos1(0);
        std::string::size_type pos2(0);
        const std::string::size_type from_len(from.size());
        const std::string::size_type to_len(to.size());
        cnt = 0;

        while ((pos1 = s.find(from, pos2)) != std::string::npos)
        {
            s.replace(pos1, from_len, to);
            pos2 = pos1 + to_len;
            ++cnt;
        }
    }

    return cnt;
}

void MISplitPath(const std::string &fullname,
                 std::string *path,
                 std::string *name,
                 std::string *ext)
{
    char separator;
#ifdef _WIN32
    separator = '\\';
#else
    separator = '/';
#endif

    if (ext)
    {
        // dot may be (and commonly -- at least under Unix -- is) the first
        // character of the filename, don't treat the entire filename as
        // extension in this case
        if (fullname[0] == '.' && fullname.size() > 1)
        {
            *ext = MIAfterLast(&fullname[1], '.');
        }
        else
        {
            *ext = MIAfterLast(fullname, '.');
        }
    }


    if (name)
    {
        *name = MIBeforeLast(fullname, '.');
        *name = MIAfterLast(*name, separator);
    }

    if (path)
    {
        *path = MIBeforeLast(fullname, separator);
    }
}

bool MIStringToNumber(const std::string &s, long &l)
{
    return (sscanf(s.c_str(), "%ld", &l) == 1);
}

bool MIStringToNumber(const std::string &s, int &i)
{
    return (sscanf(s.c_str(), "%d", &i) == 1);
}

bool MIStringToNumber(const std::string &s, float &f)
{
    return (sscanf(s.c_str(), "%f", &f) == 1);
}

void MIStringTrim(std::string &str, bool fromRight)
{
    if (fromRight)
    {
        std::string::size_type pos = str.find_last_not_of(" \f\n\r\t\v");
        if (pos != std::string::npos)
        {
            str.erase(pos + 1);
        }
        else
        {
            str.erase(str.begin(), str.end());
        }
        return;
    }

    std::string::size_type pos = str.find_first_not_of(" \f\n\r\t\v");
    if (pos != std::string::npos)
    {
        str.erase(0, pos);
    }
    else
    {
        str.erase(str.begin(), str.end());
    }
}

void MIStringSplit(const std::string &s, const std::string &delim, std::vector<std::string> &results)
{
    std::string str = s;

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
