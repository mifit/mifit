#ifndef MI_SYSTEM_H
#define MI_SYSTEM_H


// a bunch of system-level commands (string utilities, etc), not really
// dependent on chemlib, could be broken out into a separate library
#include <string>
#include <vector>

std::string MIBeforeFirst(const std::string &s, char sep);
std::string MIAfterFirst(const std::string &s, char sep);
std::string MIBeforeLast(const std::string &s, char sep);
std::string MIAfterLast(const std::string &s, char sep);

std::string MIToUpper(const std::string &s);
std::string MIToLower(const std::string &s);

bool MIStringToNumber(const std::string &s, long &l);
bool MIStringToNumber(const std::string &s, int &i);
bool MIStringToNumber(const std::string &s, float &f);

// If 'from' matches 'to' or 'from' is empty,
// does not parse 's', returns std::string::npos
// Otherwise returns number of replacements done
std::string::size_type MIStringReplace(std::string &s,
                                       const std::string &from,
                                       const std::string &to);

void MISplitPath(const std::string &fullname,
                 std::string *path,
                 std::string *name,
                 std::string *ext);

void MIStringTrim(std::string &str, bool fromRight = true);

void MIStringSplit(const std::string &str,
                   const std::string &delim,
                   std::vector<std::string> &results);


#endif // ifndef MI_SYSTEM_H
