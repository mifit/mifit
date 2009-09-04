TOP_SRCDIR = $$PWD
win32 {
  BoostDir = $$TOP_SRCDIR/../external/boost_1_39_0
}

libsDir = $$TOP_SRCDIR/libs

CONFIG -= debug release
CONFIG += release

CONFIG(debug, debug|release) {
  DEFINES += DEBUG
}

QT += opengl xml
CONFIG += no_keywords opengl ordered

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
    LIB_EXTENSION=.a
  }
}

INCLUDEPATH += $$BoostDir

foo {
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
   $$TOP_SRCDIR/apps/MIFit/figurelib \
   $$TOP_SRCDIR/apps/MIFit/jobs \
   $$BoostDir
}

win32 {
  DEFINES -= UNICODE
  DEFINES += _CRT_SECURE_NO_DEPRECATE _CRT_NONSTDC_NO_DEPRECATE
  LIBS += -lglu32 -lopengl32 -L$$BoostDir/stage/lib
  CONFIG(debug, debug|release) {
    DEFINES += DEBUG
    LIBS += -L$$BoostDir/stage/lib -llibboost_signals-mgw34-mt-d
  }
  CONFIG(release, debug|release) {
    LIBS += -L$$BoostDir/stage/lib -llibboost_signals-mgw34-mt
  }
} else {
  LIBS += -lboost_signals
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
  
