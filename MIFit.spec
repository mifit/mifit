%define name MIFit
%define is_mandrake %(test -e /etc/mandrake-release && echo 1 || echo 0)
%define is_suse %(test -e /etc/SuSE-release && echo 1 || echo 0)
%define is_fedora %(test -e /etc/fedora-release && echo 1 || echo 0)
%define qmake qmake

%if %is_fedora
%define distr %(cat /etc/fedora-release)
%define qmake qmake-qt4
%endif
%if %is_suse
%define distr %(head -1 /etc/SuSE-release)
%endif
%if %is_mandrake
%define distr %(cat /etc/mandrake-release)
%endif

Summary: MIFit
Name: %{name}
License: GPLv2
URL: http://code.google.com/p/mifit
Version: 2009.09
Release: 1
Group: Science/Crystallography
%if %{?suse_version}
Requires:  libqt4 >= 4.4.3
%endif
%if %{?fedora_version} || %{?rhel_version} || %{?centos_version}
Requires: qt
%endif

Source: %{name}-%{version}.tar.gz
BuildRoot:      %{_tmppath}/%{name}-%{version}-build
%if %{?suse_version}
BuildRequires:  libqt4-devel >= 4.4.3
%endif
BuildRequires:  boost-devel
%if %{?fedora_version} || %{?rhel_version} || %{?centos_version}
BuildRequires:  qt-devel
%endif

%description
MIFit is a cross-platform interactive graphics application
for molecular modeling, fitting, and refinement of protein
structures from x-ray crystallography.
%prep
%setup
%build
./configure --prefix=$RPM_BUILD_ROOT
make
%install
make install
%clean
rm -rf $RPM_BUILD_ROOT

%files
%defattr(-,root,root)
/MIExpert
/MIFit
/README.txt
/data
/examples

%changelog
* Thu Sep 10 2009 MIFit 2009.09
- MIFit 2009.09
