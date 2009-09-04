
include(../../../common.pri)
include(../MIFit.pri)

TEMPLATE = app
TARGET = MIFit
INCLUDEPATH += .
DESTDIR = ../../..

# to copy executable from build dir to source distribution dir
target.path = $DESTDIR
INSTALLS += target

CONFIG += no_lflags_merge
CONFIG -= staticlib

win32 {
  RC_FILE=mifit.rc
}

mac {
  ICON=../images/mifit.icns
}

include(../ui/ui.pri)
include(../figurelib/figurelib.pri)
include(../jobs/jobs.pri)
include(../wxdr/wxdr.pri)
LIBS += -lui
include(../core/core.pri)
include($${libsDir}/nongui/nongui.pri)
include($${libsDir}/ligand/ligand.pri)
include($${libsDir}/map/map.pri)
include($${libsDir}/molopt/molopt.pri)
include($${libsDir}/conflib/conflib.pri)
include($${libsDir}/chemlib/chemlib.pri)
include($${libsDir}/opengl/opengl.pri)
include($${libsDir}/math/math.pri)
include($${libsDir}/util/util.pri)
include($${libsDir}/jacgrid/jacgrid.pri)
include($${libsDir}/umtz/umtz.pri)

win32 {
  LIBS -= -lglu32 -lopengl32 -llibboost_signals-mgw34-mt
  LIBS += -lglu32 -lopengl32 -llibboost_signals-mgw34-mt
}

include(../../../rpath.pri)

SOURCES += *.cpp
