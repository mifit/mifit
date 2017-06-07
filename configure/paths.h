#pragma once

#include "strings.h"
#include <boost/filesystem.hpp>
#include <iosfwd>
#include <string>
#include <vector>

namespace nb
{
namespace fs = boost::filesystem;

class paths
{
public:
    paths();
    paths(std::initializer_list<fs::path> s);
    paths(std::vector<fs::path> s);
    paths(fs::path s);
    paths& operator+=(const fs::path& s);
    paths& operator+=(const std::initializer_list<fs::path>& s);
    paths& operator+=(const paths& s);

    paths& operator-=(const fs::path& s);
    paths& operator-=(const std::initializer_list<fs::path>& s);
    paths& operator-=(const paths& s);

    bool empty() const;
    size_t size() const;
    bool contains(const fs::path& s) const;

    const fs::path& operator[](int i) const;

    typedef std::vector<fs::path>::iterator iterator;
    typedef std::vector<fs::path>::const_iterator const_iterator;
    iterator begin();
    iterator end();
    const_iterator begin() const;
    const_iterator end() const;
    const_iterator cbegin() const;
    const_iterator cend() const;

    paths prefix_copy(const fs::path& prefix) const;
    paths suffix_copy(const fs::path& suffix) const;

    //! Exchanges the contents of the container with those of other.
    void swap(paths& other);

    explicit operator strings() const;

private:
    std::vector<fs::path> _paths;
    friend paths operator+(const paths& a, const paths& b);
    friend std::ostream& operator<<(std::ostream& a, const paths& s);
    friend paths operator+(const std::string& a, const paths& b);
    friend paths operator+(const paths& a, const std::string& b);
};

paths operator+(const paths& a, const paths& b);
paths operator+(const std::string& a, const paths& b);
paths operator+(const paths& a, const std::string& b);
paths operator/(const fs::path& a, const paths& b);
paths operator/(const paths& a, const fs::path& b);
std::ostream& operator<<(std::ostream& out, const paths& s);

} // namespace nb
