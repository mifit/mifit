
MIFit is a cross-platform interactive graphics application
for molecular modeling, fitting, and refinement of protein
structures from x-ray crystallography.

MIFit project web site: http://code.google.com/p/mifit


Installation from binaries

  Binary distributions may be available for your operating
  system. Check the MIFit project web site for more
  information.

Installation from source code

  MIFit may be built from source code using the Qt SDK
  (version 2009.03 or later) and Boost (versino 1.39.0 or
  later)

  Windows

    Start the Qt command prompt
    Change to the MIFit source directory
    .\configure.exe --boost=C:\boost_1_39_0
    mingw32-make.exe
    mingw32-make.exe install

    Run "configure.exe --help" to see more options.

  Linux

    Change to the MIFit source directory
    ./configure
    make
    make install

    Run "configure --help" to see more options.

  Mac OS

    The Linux procedure may work, but it has not been tested.
