#include "strings.h"

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

strings::strings() {}
strings::strings(initializer_list<string> s): _strings(s) {}
strings::strings(vector<string> s): _strings(s) {}
strings::strings(string s): _strings() { _strings.push_back(s); }

strings& strings::operator+=(const string& s)
{
    _strings.push_back(s);
    return *this;
}

strings& strings::operator+=(const initializer_list<string>& s)
{
    _strings.insert(_strings.end(), s.begin(), s.end());
    return *this;
}

strings& strings::operator+=(const strings& s)
{
    _strings.insert(_strings.end(), s.cbegin(), s.cend());
    return *this;
}

strings& strings::operator-=(const string& s)
{
    _strings.erase(remove(_strings.begin(), _strings.end(), s),
                   _strings.end());
    return *this;
}

strings& strings::operator-=(const initializer_list<string>& s)
{
    auto new_end = _strings.end();
    for (const auto& i : s) {
        new_end = remove(_strings.begin(), new_end, i);
    }
    _strings.erase(new_end, _strings.end());
    return *this;
}

strings& strings::operator-=(const strings& s)
{
    auto new_end = _strings.end();
    for (const auto& i : s._strings) {
        new_end = remove(_strings.begin(), new_end, i);
    }
    _strings.erase(new_end, _strings.end());
    return *this;
}

bool strings::empty() const { return _strings.empty(); }
size_t strings::size() const { return _strings.size(); }

bool strings::contains(const string& s) const
{
    return find(_strings.cbegin(), _strings.cend(), s) != _strings.cend();
}

const string& strings::operator[](int i) const
{
    return _strings[i];
}

std::string strings::join(const string& separator) const
{
    return boost::algorithm::join(_strings, separator);
}

strings::iterator strings::begin() { return _strings.begin(); }
strings::iterator strings::end() { return _strings.end(); }
strings::const_iterator strings::cbegin() const { return _strings.cbegin(); }
strings::const_iterator strings::cend() const { return _strings.cend(); }

strings strings::prefix_copy(const string& prefix) const
{
    strings result(_strings);
    for (auto& s : result)
        s.insert(0, prefix);
    return result;
}

strings strings::suffix_copy(const string& suffix) const
{
    strings result(_strings);
    for (auto& s : result)
        s.append(suffix);
    return result;
}

void strings::swap(strings &other)
{
    _strings.swap(other._strings);
}

strings strings::split(const string &s)
{
    vector<string> results;
    if (!s.empty()) {
        boost::split(results, s, boost::is_any_of(" \t\n"), boost::token_compress_on);
    }
    return strings(results);
}

strings operator+(const strings& a, const strings& b)
{
    strings result(a);
    result += b;
    return result;
}

strings operator+(const string& a, const strings& b)
{
    return b.prefix_copy(a);
}

strings operator+(const strings& a, const string& b)
{
    return a.suffix_copy(b);
}

ostream& operator<<(ostream& out, const strings& s)
{
    copy(s._strings.begin(), s._strings.end(),
         ostream_iterator<string>(out, " "));
    return out;
}

} // namespace nb
