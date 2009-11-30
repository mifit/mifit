#include "MIConfig.h"

MIConfig::MIConfig()
{
}
MIConfig*MIConfig::_instance = NULL;

bool MIConfig::Read(const std::string &key, std::string &value)
{
    return Read(key, value, "", false);
}

bool MIConfig::Read(const std::string &key, long *value)
{
    return Read(key, value, long(), false);
}

bool MIConfig::Read(const std::string &key, bool *value)
{
    return Read(key, value, bool(), false);
}

bool MIConfig::Read(const std::string &key, double *value)
{
    return Read(key, value, double(), false);
}

int MIConfig::GetProfileInt(const std::string &s, const std::string &i, int d)
{
    std::string key = s;
    key += "/";
    key += i;
    long l;
    Read(key, &l, (long)d);
    return (int)l;
}

std::string MIConfig::GetProfileString(const std::string &s, const std::string &i, const std::string &d)
{
    std::string key = s;
    key += "/";
    key += i;
    std::string str;
    Read(key, str, d);
    return str;
}

void MIConfig::WriteProfileInt(const std::string &s, const std::string &i, int d)
{
    std::string key = s;
    key += "/";
    key += i;
    Write(key, (long)d);
}

void MIConfig::WriteProfileString(const std::string &s, const std::string &i, const std::string &d)
{
    std::string key = s;
    key += "/";
    key += i;
    Write(key, d);
}
