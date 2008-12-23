TOP_SRCDIR = $$[mifitDir]
ExternalDir = $$TOP_SRCDIR/external
BoostDir = $$ExternalDir/boost_1_36_0
PythonDir = $$ExternalDir/Python-2.5.2

CONFIG += $$[config]

DEFINES += \
    USE_ASPLOT \
    USE_NAV_WINDOW \
    USE_QT_RAMAPLOT \
    MI_USE_JOBS \
    USE_DICT_EDITOR \
    MI_USE_TREE

CONFIG(debug, debug|release) {
  DEFINES += DEBUG
}

CONFIG += \
    USE_ASPLOT \
    USE_NAV_WINDOW \
    USE_QT_RAMAPLOT \
    MI_USE_JOBS \
    USE_DICT_EDITOR \
    MI_USE_TREE

QT += opengl xml
CONFIG += no_keywords USE_ASPLOT MI_USE_JOBS opengl ordered

!USE_SHARED_LIBS {
  CONFIG += staticlib 
}


USE_SHARED_LIBS {
  mac {
    LIB_EXTENSION=.dylib
  }
  linux {
    LIB_EXTENSION=.so
  }
  win32 {
    LIB_EXTENSION=.dll
  }
} else {
  unix {
    LIB_EXTENSION=.a
  }
  win32 {
    LIB_EXTENSION=.lib
  }
}

INCLUDEPATH += . \
   $$TOP_SRCDIR/libs \
   $$TOP_SRCDIR/libs/chemlib \
   $$TOP_SRCDIR/libs/conflib \
   $$TOP_SRCDIR/libs/ligand \
   $$TOP_SRCDIR/libs/map \
   $$TOP_SRCDIR/libs/math \
   $$TOP_SRCDIR/libs/molopt \
   $$TOP_SRCDIR/libs/nongui \
   $$TOP_SRCDIR/libs/util \
   $$TOP_SRCDIR/libs/umtz \
   $$TOP_SRCDIR/libs/jacgrid \
   $$TOP_SRCDIR/apps/MIFit \
   $$TOP_SRCDIR/apps/MIFit/core \
   $$TOP_SRCDIR/apps/MIFit/images \
   $$TOP_SRCDIR/apps/MIFit/wxdr \
   $$TOP_SRCDIR/apps/MIFit/ui \
   $$TOP_SRCDIR/apps/MIFit/python \
   $$PythonDir/Include \
   $$PythonDir \
   $$BoostDir

USE_ASPLOT {
  INCLUDEPATH += $$TOP_SRCDIR/apps/MIFit/figurelib
}
MI_USE_JOBS {
  INCLUDEPATH += $$TOP_SRCDIR/apps/MIFit/jobs
}

DEFINES += BOOST_PYTHON_NO_LIB
CONFIG(debug, debug|release) {
  DEFINES += BOOST_DEBUG_PYTHON
}

win32 {
  DEFINES -= UNICODE
  DEFINES += _CRT_SECURE_NO_DEPRECATE _CRT_NONSTDC_NO_DEPRECATE
  LIBS += -L$$BoostDir/stage/lib
  CONFIG(debug, debug|release) {
    DEFINES += DEBUG
    LIBS += $$PythonDir/PCbuild8/win32debug/python25_d.lib $$BoostDir/stage/lib/boost_python-vc80-mt-gyd-1_36.lib
  }
  CONFIG(release, debug|release) {
    LIBS += $$PythonDir/PCbuild8/win32release/python25.lib $$BoostDir/stage/lib/boost_python-vc80-mt-1_36.lib
  }
} else {
  LIBS += $$BoostDir/libboost_signals.a $$BoostDir/libboost_python.a $$PythonDir/libpython2.5.a -lutil
}


unix {
  QMAKE_CXXFLAGS_WARN_ON += -Wno-unknown-pragmas
  QMAKE_CFLAGS_WARN_ON += -Wno-unknown-pragmas
  INCLUDEPATH += /usr/include/X11
}

mac {
  QMAKE_LFLAGS_SHLIB += -undefined dynamic_lookup
  QMAKE_CXXFLAGS_WARN_ON += -Wno-unknown-pragmas
  QMAKE_CFLAGS_WARN_ON += -Wno-unknown-pragmas
}
  
