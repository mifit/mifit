#pragma once

#include <iosfwd>
#include <string>
#include <vector>

namespace nb
{

class strings
{
public:
    strings();
    strings(std::initializer_list<std::string> s);
    strings(std::vector<std::string> s);
    strings(std::string s);
    strings& operator+=(const std::string& s);
    strings& operator+=(const std::initializer_list<std::string>& s);
    strings& operator+=(const strings& s);

    strings& operator-=(const std::string& s);
    strings& operator-=(const std::initializer_list<std::string>& s);
    strings& operator-=(const strings& s);

    bool empty() const;
    size_t size() const;
    bool contains(const std::string& s) const;

    const std::string& operator[](int i) const;

    std::string join(const std::string& separator) const;

    typedef std::vector<std::string>::iterator iterator;
    typedef std::vector<std::string>::const_iterator const_iterator;
    iterator begin();
    iterator end();
    const_iterator cbegin() const;
    const_iterator cend() const;

    strings prefix_copy(const std::string& prefix) const;
    strings suffix_copy(const std::string& suffix) const;

    //! Exchanges the contents of the container with those of other.
    void swap(strings& other);

    static strings split(const std::string& s);

private:
    std::vector<std::string> _strings;
    friend strings operator+(const strings& a, const strings& b);
    friend std::ostream& operator<<(std::ostream& a, const strings& s);
};

strings operator+(const strings& a, const strings& b);
strings operator+(const std::string& a, const strings& b);
strings operator+(const strings& a, const std::string& b);
std::ostream& operator<<(std::ostream& out, const strings& s);

} // namespace nb
