#include <fstream>
#include <iostream>
#include <boost/filesystem.hpp>
#include <set>
#include <sstream>
#include <string>
#include <vector>
#include "make_relative.h"
#include "ninja_file.h"

using namespace std;
namespace fs = boost::filesystem;

int grep(const fs::path &filename, const string &keyword)
{
    int counter = 0;
    string line;
    ifstream in(filename.string());
    if (in.is_open())
    {
        while (getline(in, line))
        {
            if (line.find(keyword) != string::npos)
                counter++;
        }
    }
    return counter;
}

int main(int argc, char *argv[])
{
    try {
        vector<string> args { argv+1, argv+argc };

        if (args.empty())
        {
            cerr << "Usage: configure <platform>\n"
                 << "\n"
                 << "  platform   mingw or trusty\n"
                 << endl;
            return EXIT_FAILURE;
        }
        auto platform = args.front();

        static set<string> ignoredDirs {
            ".git", "o", "data", "docs", "examples", "packaging", "configure",
            "apps/MIFit/python", "apps/MIFit/script", "apps/MIFlex", "libs/nongui"
        };
        static set<string> ignoredFiles {
            "libs/umtz/mmtztest.c",
            "libs/umtz/umtzdiff.c",
            "libs/umtz/umtztest.c",
        };

        static set<string> extensions {
            ".cpp", ".c", ".h", ".ui",
        };

        static set<fs::path> results;

        static auto outputFile = [&](fs::directory_entry entry) {
            auto relative_path = make_relative(fs::current_path(), entry.path()).generic_string();
            if (fs::is_regular_file(entry))
            {
                if (extensions.find(entry.path().extension().string()) != end(extensions))
                    results.insert(relative_path);
            }
        };
        static auto isIgnoredDir = [&](fs::directory_entry entry) {
            if (!fs::is_directory(entry))
                return false;
            auto relative_path = make_relative(fs::current_path(), entry.path()).generic_string();
            if (ignoredDirs.find(relative_path) != ignoredDirs.end())
                return true;
            return false;
        };

        fs::path base_path(fs::current_path());
        for_each(fs::directory_iterator(base_path), fs::directory_iterator(),
                 [](fs::directory_entry entry) {
            if (isIgnoredDir(entry)) return;
            if (fs::is_directory(entry)) {
                auto dirIter = fs::recursive_directory_iterator(entry);
                for_each(dirIter, fs::recursive_directory_iterator(),
                         [&](fs::directory_entry entry) {
                    if (isIgnoredDir(entry))
                    {
                        dirIter.no_push();
                        return;
                    }
                    outputFile(entry);
                });
            }
        });

        nb::ninja_file ninjas { "build.ninja", 160 };
        ninjas.include(platform + ".ninja");
        ninjas.include("common.ninja");
        ninjas.variable("o_dir", "o");
        ninjas.variable("base_dir", base_path.generic_string());
        fs::path object_dir { "$o_dir" };
        fs::path base_dir { "$base_dir" };

        nb::paths uic_targets;
        nb::paths object_files;
        for (auto s : results) {
            auto relative_path = make_relative(fs::current_path(), s).generic_string();
            if (ignoredFiles.find(relative_path) != ignoredFiles.end())
                continue;
            if (s.extension() == ".h")
            {
                if (grep(s, "Q_OBJECT"))
                {
                    fs::path moc_file = s;
                    moc_file.replace_extension("moc.cpp");
                    ninjas << nb::build(object_dir/moc_file, "moc", base_dir/s);

                    auto object_file = moc_file;
                    object_file.replace_extension(".o");

                    ninjas << nb::build(object_dir/object_file, "cxx", object_dir/moc_file);
                    object_files += object_dir/object_file;
                }
            }
            else if (s.extension() == ".cpp")
            {
                auto object_file = s;
                object_file.replace_extension(".o");

                nb::paths implicit_paths;
                auto ui_file = s;
                ui_file.replace_extension(".ui");
                if (fs::exists(ui_file))
                {
                    auto ui_header_file = s.parent_path()/("ui_" + s.filename().string());
                    ui_header_file.replace_extension(".h");
                    implicit_paths += base_dir/ui_header_file;
                }

                if (grep(s, "Q_OBJECT"))
                {
                    fs::path moc_file = s;
                    moc_file.replace_extension("moc");
                    ninjas << nb::build(base_dir/moc_file, "moc", base_dir/s);
                    implicit_paths += base_dir/moc_file;
                }

                ninjas << nb::build(object_dir/object_file, "cxx", base_dir/s)
                          .implicit(implicit_paths);
                object_files += object_dir/object_file;

            }
            else if (s.extension() == ".c")
            {
                auto object_file = s;
                object_file.replace_extension(".o");
                ninjas << nb::build(object_dir/object_file, "cc", base_dir/s);
                object_files += object_dir/object_file;
            }
            else if (s.extension() == ".ui")
            {
                auto ui_header_file = s.parent_path()/("ui_" + s.filename().string());
                ui_header_file.replace_extension(".h");
                ninjas << nb::build(base_dir/ui_header_file, "uic", base_dir/s);
                uic_targets += base_dir/ui_header_file;
            }
        }

        ninjas << nb::build("uics", "phony", uic_targets);

        ninjas << nb::build("MIFit.exe", "link", object_files);

        return EXIT_SUCCESS;
    }
    catch (const string &error)
    {
        cerr << "fatal: " << error << endl;
    }
    catch (...)
    {
        cerr << "unknown exception\n";
    }
    return EXIT_FAILURE;
}
