include(../MIFit.pri)

TEMPLATE = app
TARGET = MIFit
DESTDIR = ../../..

target.path = $$PREFIX
INSTALLS += target

win32 {
  RC_FILE=mifit.rc
}

mac {
  ICON=../images/mifit.icns
}

INCLUDEPATH *= ..
LIBS *= -L..
appLibs = ui figurelib jobs core
for(l, appLibs) {
  LIBS *= -l$$qtLibraryTarget($${l})
  PRE_TARGETDEPS += ../lib$$qtLibraryTarget($${l})$$LIB_EXTENSION
}

INCLUDEPATH *= ../../../libs
LIBS *= -L../../../libs

otherLibs = nongui ligand map molopt conflib chemlib miopengl mimath miutil jacgrid umtz
for(l, otherLibs) {
  LIBS *= -l$$qtLibraryTarget($${l})
  PRE_TARGETDEPS += ../../../libs/lib$$qtLibraryTarget($${l})$$LIB_EXTENSION
}

win32 {
  LIBS -= -lglu32 -lopengl32
  LIBS += -lglu32 -lopengl32
}

include(../../../rpath.pri)

SOURCES += *.cpp
