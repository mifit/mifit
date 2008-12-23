#include "mi_python.h"

#include <boost/python.hpp>
#include <boost/python/reference_existing_object.hpp>
#include <boost/python/return_value_policy.hpp>
#include <boost/python/suite/indexing/vector_indexing_suite.hpp>

#include <vector>

#include "uilib.h"
#include "container_conversions.h"

using namespace boost::python;
using namespace container_conversions;

BOOST_PYTHON_MEMBER_FUNCTION_OVERLOADS(MIMainWindow_openFiles_overload, OpenFiles, 1, 2)

BOOST_PYTHON_MODULE(mi) {
  from_python_sequence< std::vector<std::string>, variable_capacity_policy>();

  class_< std::vector<std::string> >("vec_string")
    .def(vector_indexing_suite< std::vector<std::string>, true>());

  class_<Application>("Application", no_init)
    .def("instance", &Application::instance, return_value_policy<reference_existing_object>())
    .staticmethod("instance")
    .add_property("binDir", make_function(&Application::GetBinDirectory,
        return_value_policy<copy_const_reference>()))
    .add_property("mifitHome", make_function(&Application::GetMolimageHome,
        return_value_policy<copy_const_reference>()))
    .def("getBinDir", &Application::GetBinDirectory,
        return_value_policy<copy_const_reference>())
    .def("getMIFitHome", &Application::GetMolimageHome,
        return_value_policy<copy_const_reference>());

  def("app", &Application::instance, return_value_policy<reference_existing_object>());

  class_<MIMainWindow, boost::noncopyable>("MIMainWindow", no_init)
    .def("instance", &MIMainWindow::instance, return_value_policy<reference_existing_object>())
    .staticmethod("instance")
    .def("openFiles", &MIMainWindow::OpenFiles,
        MIMainWindow_openFiles_overload(args("files", "newWin"),
            "Open files (optionally into new windows)"))
    .def("log", &MIMainWindow::Log)
    .def("debug", &MIMainWindow::Debug);

  def("app", &Application::instance, return_value_policy<reference_existing_object>());
  def("window", &MIMainWindow::instance, return_value_policy<reference_existing_object>());
}
