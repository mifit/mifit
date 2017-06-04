include(../MIFit.pri)

TEMPLATE = app
TARGET = MIFit
DESTDIR = $$PWD/../../..

target.path = $$PREFIX
INSTALLS += target

win32 {
  RC_FILE=mifit.rc
}

mac {
  ICON=$$PWD/../images/mifit.icns
}

INCLUDEPATH *= $$OUT_PWD/..
LIBS *= -L$$OUT_PWD/..
appLibs = ui jobs script core
for(l, appLibs) {
  LIBS *= -l$$qtLibraryTarget($${l})
  PRE_TARGETDEPS += $$OUT_PWD/../lib$$qtLibraryTarget($${l})$$LIB_EXTENSION
}

INCLUDEPATH *= $$OUT_PWD/../../../libs
LIBS *= -L$$OUT_PWD/../../../libs

otherLibs = nongui ligand map molopt conflib chemlib miopengl mimath miutil jacgrid umtz
for(l, otherLibs) {
  LIBS *= -l$$qtLibraryTarget($${l})
  PRE_TARGETDEPS += $$OUT_PWD/../../../libs/lib$$qtLibraryTarget($${l})$$LIB_EXTENSION
}

win32 {
  LIBS -= -lglu32 -lopengl32
  LIBS += -lglu32 -lopengl32
}

unix {
  LIBS += -lGLU
}

include(../../../rpath.pri)

SOURCES += $$files(*.cpp)
