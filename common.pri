TOP_SRCDIR = $$PWD
win32 {
  BoostDir = $$TOP_SRCDIR/../external/boost_1_39_0
}

libsDir = $$TOP_SRCDIR/libs

CONFIG -= debug release
CONFIG += release
#CONFIG += create_prl link_prl depend_prl

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

win32 {
  DEFINES -= UNICODE
  DEFINES += _CRT_SECURE_NO_DEPRECATE _CRT_NONSTDC_NO_DEPRECATE
  LIBS += -lglu32 -lopengl32
  CONFIG(debug, debug|release) {
    DEFINES += DEBUG
  }
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
  
