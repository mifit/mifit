Name:           MIFit
Version:        2009.09
Release:        1
Summary:        A protein crystallography fitting application
License:        GPLv2
URL:            http://code.google.com/p/mifit
Group:          Science/Crystallography

Source:         %{name}-%{version}-2.tar.gz
BuildRoot:      %{_tmppath}/%{name}-%{version}-build
BuildRequires:  boost-devel

%if 0%{?suse_version}
Requires:       libqt4
BuildRequires:  libqt4-devel
%endif

%if 0%{?fedora_version} || 0%{?rhel_version} || 0%{?centos_version}
Requires:       qt4 >= 4.3
BuildRequires:  qt4-devel >= 4.3
BuildRequires:  gcc-c++
%endif

BuildRequires:  -post-build-checks

%description
MIFit is a cross-platform interactive graphics application
for molecular modeling, fitting, and refinement of protein
structures from x-ray crystallography.

%define mifitdir /opt/%{name}

%if 0%{?fedora_version}
%define qmake qmake-qt4
%else
%if 0%{?rhel_version} || 0%{?centos_version}
%define qmake /usr/lib/qt4/bin/qmake
%else
%define qmake qmake
%endif
%endif

%prep
%setup

%build
%qmake -r MI.pro PREFIX=$RPM_BUILD_ROOT/%{mifitdir}
%__make

%install
%__make install

%clean
rm -rf $RPM_BUILD_ROOT

%files
%defattr(-,root,root)
%{mifitdir}/MIExpert
%{mifitdir}/MIFit
%{mifitdir}/README.txt
%{mifitdir}/license.txt
%{mifitdir}/data
%{mifitdir}/examples

%changelog
