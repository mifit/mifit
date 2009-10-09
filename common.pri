include(conf.pri)

!isEmpty(PREFIX) {
  linux {
    PREFIX = /opt/MIFit
  }
  win32 {
    PREFIX = /MIFit
  }
}

TOP_SRCDIR = $$PWD

libsDir = $$TOP_SRCDIR/libs

INCLUDEPATH *= $$libsDir

MOC_DIR = moc

CONFIG(debug, debug|release) {
  DEFINES += DEBUG
}

QT += opengl xml network script
CONFIG += opengl ordered

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

!isEmpty(BOOST) {
  INCLUDEPATH += $$BOOST
}

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
  
