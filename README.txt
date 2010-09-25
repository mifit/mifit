
MIFit is a cross-platform interactive graphics application
for molecular modeling, fitting, and refinement of protein
structures from x-ray crystallography.

MIFit project web site: http://code.google.com/p/mifit

MIFit is licensed under the GNU GPL. See the license.txt
file for details.

Installation from binaries

  Binary distributions may be available for your operating
  system. Check the MIFit project web site for more
  information.

Installation from source code

  MIFit may be built from source code using the Qt SDK
  (version 2010.05 or later) and Boost (version 1.39.0 or
  later)

  Windows

    Start the Qt command prompt
    Change to the MIFit source directory
    qmake -r MI.pro PREFIX=C:\MIFit BOOST=c:\boost_1_39_0
    mingw32-make.exe release
    mingw32-make.exe release install

  Linux

    Change to the MIFit source directory
    qmake -r MI.pro PREFIX=/opt/MIFit
    make
    make install

  Mac OS

    The Linux procedure may work, but it has not been tested.
