#pragma once

#include "paths.h"
#include "strings.h"
#include <boost/variant.hpp>
#include <fstream>
#include <iosfwd>
#include <string>

namespace nb
{
namespace fs = boost::filesystem;

class rule;
class build;
class ninjas;

typedef std::vector<std::pair<std::string, std::string>> variables;

class ninja_file
{
public:
    ninja_file(fs::path file, size_t width = 80);
    ~ninja_file();

    fs::path file() const { return _file; }
    void newline();
    void comment(std::string text);
    void variable(std::string key, std::string value, int indent = 0);
    void include(fs::path file);
    void subninja(fs::path file);
    void defaultTarget(fs::path target);
    void defaultTarget(paths targets);

    ninja_file& operator<<(const rule& r);
    ninja_file& operator<<(const build& b);
    ninja_file& operator<<(const ninjas& n);

    template<typename F>
    ninja_file& operator<<(F f);

    operator paths() { return _lastTargetOutputs; }

private:
    fs::path _file;
    std::size_t _width;
    std::ofstream _out;
    paths _lastTargetOutputs;

    void _line(std::string text, int indent = 0);
};

template<typename F>
ninja_file& ninja_file::operator<<(F f)
{
  return f(*this);
}

typedef std::function<ninja_file& (ninja_file&)> ninja_manip;

ninja_file& newline(ninja_file& n);
ninja_manip comment(std::string comment);
ninja_manip variable(std::string key, std::string value, int indent = 0);
ninja_manip variable(variables v);
ninja_manip include(fs::path file);
ninja_manip include(paths files);
ninja_manip subninja(fs::path file);
ninja_manip subninja(paths files);
ninja_file& default_target(ninja_file& n);
ninja_manip set_default_target(paths targets);


class rule
{
public:
    rule(std::string name, std::string command);
    rule& description(std::string text);
    rule& depfile(std::string depfile);
    rule& generator(bool value);
    rule& restat(bool value);

    ninja_file& operator<<(ninja_file& n);

private:
    friend class ninja_file;
    std::string _name;
    std::string _command;
    std::string _description;
    std::string _depfile;
    bool _generator;
    bool _restat;
};

class build
{
public:
    build(fs::path output, std::string rule);
    build(paths outputs, std::string rule);
    build(fs::path output, std::string rule, fs::path input);
    build(fs::path output, std::string rule, paths inputs);
    build(paths output, std::string rule, fs::path input);
    build(paths outputs, std::string rule, paths inputs);

    build& inputs(fs::path input);
    build& inputs(paths inputs);
    build& implicit(paths targets);
    build& order_only(paths targets);
    build& variables(nb::variables targets);

private:
    friend class ninja_file;
    friend class builds;
    paths _outputs;
    std::string _rule;
    paths _inputs;
    paths _implicit;
    paths _order_only;
    std::vector<std::pair<std::string, std::string>> _variables;
};

typedef boost::variant<rule, build> ninja;

class ninjas
{
public:
    ninjas() {}
    ninjas(ninja n) { _ninjas.push_back(n); }
    ninjas& operator+=(const ninja& n);
    ninjas& operator+=(const ninjas& n);

    ninjas& description(std::string text);
    ninjas& depfile(std::string depfile);
    ninjas& generator(bool value);
    ninjas& restat(bool value);

    ninjas& implicit(paths targets);
    ninjas& order_only(paths targets);
    ninjas& variables(std::vector<std::pair<std::string, std::string>> targets);

private:
    friend class ninja_file;
    std::vector<ninja> _ninjas;
};

} // namespace nb
