#include "ninja_file.h"

#include "lambda_visitor.h"
#include <boost/algorithm/string.hpp>

using namespace std;

namespace nb
{

ninja_file::ninja_file(fs::path file, size_t width):
    _file(move(file)), _width(width), _out(_file.generic_string())
{
    assert(_out.is_open());
}

ninja_file::~ninja_file()
{
    _out << flush;
    _out.close();
}

void ninja_file::newline()
{
    _out << '\n';
}


void ninja_file::comment(string text)
{
    boost::algorithm::replace_all(text, "\n", " ");
    string leading = "# ";
    while (text.length() > _width) {
        // The text is too wide; wrap if possible.

        // Find the rightmost space that would obey our width constraint.
        int available_space = _width - leading.length();
        auto space = text.rfind(' ', available_space);
        if (space == string::npos) {
            // No such space; just use the first space we can find.
            space = text.find(' ', available_space);
        }
        if (space == string::npos) {
            // Give up on breaking.
            break;
        }

        _out << leading << text.substr(0, space) << '\n';
        text.erase(0, space+1);
    }
    _out << leading << text << '\n';
}

void ninja_file::variable(string key, string value, int indent)
{
    if (value.empty()) return;
    _line(key + " = " + value, indent);
}

void ninja_file::include(fs::path file)
{
    _line("include " + file.generic_string());
}

void ninja_file::subninja(fs::path file)
{
    _line("subninja " + file.generic_string());
}

void ninja_file::defaultTarget(fs::path target)
{
    _line("default " +  target.generic_string());
}

void ninja_file::defaultTarget(paths targets)
{
    _line("default " + strings(targets).join(" "));
}

// Write 'text' word-wrapped at _width characters.
void ninja_file::_line(std::string text, int indent)
{
    string leading_space(2*indent, ' ');
    boost::algorithm::replace_all(text, "\n", " ");
    while (text.length() > _width) {
        // The text is too wide; wrap if possible.

        // Find the rightmost space that would obey our width constraint.
        int available_space = _width - leading_space.length() - 2;
        auto space = text.rfind(' ', available_space);
        if (space == string::npos) {
            // No such space; just use the first space we can find.
            space = text.find(' ', available_space);
        }
        if (space == string::npos) {
            // Give up on breaking.
            break;
        }

        _out << leading_space << text.substr(0, space) <<  " $\n";
        text.erase(0, space+1);

        // Subsequent lines are continuations, so indent them.
        leading_space = string(2*(indent+2), ' ');
    }
    _out << leading_space << text << '\n';
}

ninja_file& newline(ninja_file& n)
{
    n.newline();
    return n;
}

ninja_manip comment(std::string comment)
{
    return [=](ninja_file& n) -> ninja_file& {
        n.comment(comment);
        return n;
    };
}

ninja_manip variable(std::string key, std::string value, int indent)
{
    return [=](ninja_file& n) -> ninja_file& {
        n.variable(key, value, indent);
        return n;
    };
}

ninja_manip variable(variables v)
{
    return [=](ninja_file& n) -> ninja_file& {
        for (const auto& pair : v)
            n.variable(pair.first, pair.second);
        return n;
    };
}

ninja_manip include(fs::path file)
{
    return [=](ninja_file& n) -> ninja_file& {
        n.include(file);
        return n;
    };
}

ninja_manip include(paths files)
{
    return [=](ninja_file& n) -> ninja_file& {
        for (auto file : files)
            n.include(file);
        return n;
    };
}

ninja_manip subninja(fs::path file)
{
    return [=](ninja_file& n) -> ninja_file& {
        n.subninja(file);
        return n;
    };
}

ninja_manip subninja(paths files)
{
    return [=](ninja_file& n) -> ninja_file& {
        for (auto file : files)
            n.subninja(file);
        return n;
    };
}

ninja_file& default_target(ninja_file& n)
{
    n.defaultTarget(paths(n));
    return n;
}

ninja_manip set_default_target(paths targets)
{
    return [=](ninja_file& n) -> ninja_file& {
        n.defaultTarget(targets);
        return n;
    };
}


ninja_file& ninja_file::operator<<(const rule& r)
{
    _line("rule " + r._name);
    variable("command", r._command, 1);
    if (!r._description.empty())
        variable("description", r._description, 1);
    if (!r._depfile.empty())
        variable("depfile", r._depfile, 1);
    if (r._generator)
        variable("generator", "true", 1);
    if (r._restat)
        variable("restat", "true", 1);
    return *this;
}

rule::rule(std::string name, std::string command):
    _name(name), _command(command), _generator(false), _restat(false)
{}

rule& rule::description(std::string text) { _description = text; return *this; }
rule& rule::depfile(std::string depfile) { _depfile = depfile; return *this; }
rule& rule::generator(bool value) { _generator = value; return *this; }
rule& rule::restat(bool value) { _restat = value; return *this; }

ninja_file& ninja_file::operator<<(const build& b)
{
    if (b._rule.empty()) {
        return *this;
    }
    auto all_inputs = b._inputs;

    if (!b._implicit.empty()) {
        all_inputs += "|";
        all_inputs += b._implicit;
    }
    if (!b._order_only.empty()) {
        all_inputs += "||";
        all_inputs += b._order_only;
    }

    _line("build " + strings(b._outputs).join(" ") + ": " + b._rule + " "
          + strings(all_inputs).join(" "));

    for (auto pair : b._variables) {
        variable(pair.first, pair.second, 1);
    }

    _lastTargetOutputs = b._outputs;
    return *this;
}

build::build(fs::path output, std::string rule):
    _outputs( { output }), _rule(rule)
{}

build::build(paths outputs, std::string rule):
    _outputs(outputs), _rule(rule)
{}

build::build(fs::path output, std::string rule, fs::path input):
    _outputs( { output }), _rule(rule), _inputs({ input })
{}

build::build(fs::path output, std::string rule, paths inputs):
    _outputs( { output }), _rule(rule), _inputs(inputs)
{}

build::build(paths outputs, std::string rule, fs::path input):
    _outputs(outputs), _rule(rule), _inputs({ input })
{}

build::build(paths outputs, std::string rule, paths inputs):
    _outputs(outputs), _rule(rule), _inputs(inputs)
{}

build& build::inputs(fs::path input) { _inputs += input; return *this; }
build& build::inputs(paths inputs) { _inputs += inputs; return *this; }
build& build::implicit(paths targets) { _implicit += targets; return *this; }
build& build::order_only(paths targets) { _order_only += targets; return *this; }
build& build::variables(nb::variables targets)
{
    _variables.insert(_variables.end(), targets.begin(), targets.end());
    return *this;
}

ninjas& ninjas::operator+=(const ninja& n)
{
    _ninjas.push_back(n);
    return *this;
}

ninjas& ninjas::operator+=(const ninjas& n)
{
    _ninjas.insert(_ninjas.end(), n._ninjas.begin(), n._ninjas.end());
    return *this;
}

ninjas& ninjas::implicit(paths targets)
{
    auto v = make_lambda_visitor<void>(
            [](rule&) {},
            [&](build& b) { b.implicit(targets); });
    for (auto& n : _ninjas) boost::apply_visitor(v, n);
    return *this;
}

ninjas& ninjas::order_only(paths targets)
{
    auto v = make_lambda_visitor<void>(
            [](rule&) {},
            [&](build& b) { b.order_only(targets); });
    for (auto& n : _ninjas) boost::apply_visitor(v, n);
    return *this;
}

ninjas& ninjas::variables(std::vector<std::pair<std::string, std::string>> targets)
{
    auto v = make_lambda_visitor<void>(
            [](rule&) {},
            [&](build& b) { b.variables(targets); });
    for (auto& n : _ninjas) boost::apply_visitor(v, n);
    return *this;
}

ninja_file& ninja_file::operator<<(const ninjas& n)
{
    paths results;
    auto v = make_lambda_visitor<void>(
            [&](const rule& r) { (*this) << r; },
            [&](const build& b) {  (*this) << b; results += _lastTargetOutputs;  });
    for (const auto& i : n._ninjas) {
        boost::apply_visitor(v, i);
    }
    _lastTargetOutputs = results;
    return *this;
}

} // namespace nb
