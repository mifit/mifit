#include "paths.h"

#include <boost/algorithm/string.hpp>
#include <algorithm>
#include <cassert>
#include <iostream>
#include <iterator>
#include <sstream>
#include <stdexcept>

using namespace std;

namespace nb
{

paths::paths() {}
paths::paths(initializer_list<fs::path> s): _paths(s) {}
paths::paths(vector<fs::path> s): _paths(s) {}
paths::paths(fs::path s): _paths() { _paths.push_back(s); }

paths& paths::operator+=(const fs::path& s)
{
    _paths.push_back(s);
    return *this;
}

paths& paths::operator+=(const initializer_list<fs::path>& s)
{
    _paths.insert(_paths.end(), s.begin(), s.end());
    return *this;
}

paths& paths::operator+=(const paths& s)
{
    _paths.insert(_paths.end(), s.cbegin(), s.cend());
    return *this;
}

paths& paths::operator-=(const fs::path& s)
{
    _paths.erase(remove(_paths.begin(), _paths.end(), s),
                   _paths.end());
    return *this;
}

paths& paths::operator-=(const initializer_list<fs::path>& s)
{
    auto new_end = _paths.end();
    for (const auto& i : s) {
        new_end = remove(_paths.begin(), new_end, i);
    }
    _paths.erase(new_end, _paths.end());
    return *this;
}

paths& paths::operator-=(const paths& s)
{
    auto new_end = _paths.end();
    for (const auto& i : s._paths) {
        new_end = remove(_paths.begin(), new_end, i);
    }
    _paths.erase(new_end, _paths.end());
    return *this;
}

bool paths::empty() const { return _paths.empty(); }
size_t paths::size() const { return _paths.size(); }

bool paths::contains(const fs::path& s) const
{
    return find(_paths.cbegin(), _paths.cend(), s) != _paths.cend();
}

const fs::path& paths::operator[](int i) const
{
    return _paths[i];
}

paths::iterator paths::begin() { return _paths.begin(); }
paths::iterator paths::end() { return _paths.end(); }
paths::const_iterator paths::begin() const { return _paths.begin(); }
paths::const_iterator paths::end() const { return _paths.end(); }
paths::const_iterator paths::cbegin() const { return _paths.cbegin(); }
paths::const_iterator paths::cend() const { return _paths.cend(); }

paths paths::prefix_copy(const fs::path& prefix) const
{
    paths result(_paths);
    for (auto& s : result)
        s = prefix / s;
    return result;
}

paths paths::suffix_copy(const fs::path& suffix) const
{
    paths result(_paths);
    for (auto& s : result)
        s /= suffix;
    return result;
}

void paths::swap(paths &other)
{
    _paths.swap(other._paths);
}

paths operator+(const paths& a, const paths& b)
{
    paths result(a);
    result += b;
    return result;
}

paths operator+(const std::string& a, const paths& b)
{
    paths result(b._paths);
    for (auto& s : result)
        s = a + s.generic_string();
    return result;
}

paths operator+(const paths& a, const std::string& b)
{
    paths result(a._paths);
    for (auto& s : result)
        s = s.generic_string() + b;
    return result;
}

paths operator/(const fs::path& a, const paths& b)
{
    return b.prefix_copy(a);
}

paths operator/(const paths& a, const fs::path& b)
{
    return a.suffix_copy(b);
}

ostream& operator<<(ostream& out, const paths& s)
{
    copy(s._paths.begin(), s._paths.end(),
         ostream_iterator<fs::path>(out, " "));
    return out;
}

paths::operator strings() const
{
    strings results;
    for (const auto& p : _paths) {
        results += p.generic_string();
    }
    return results;
}

} // namespace nb
