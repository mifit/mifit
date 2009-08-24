
include(../../../common.pri)

TEMPLATE = app
TARGET = MIFit
INCLUDEPATH += .
DESTDIR = ../../..

# to copy executable from build dir to source distribution dir
target.path = $DESTDIR
INSTALLS += target

CONFIG += no_lflags_merge
CONFIG -= staticlib

unix { 
  LIBPREFIX = /lib
  LIBSPATH = ../../../libs
}
win32 {
  LIBSPATH = ../../../libs
  LIBPREFIX = /lib
  RC_FILE=mifit.rc
}

mac {
ICON=../images/mifit.icns
}

# specify library dependencies
PRE_TARGETDEPS += \
        ..$${LIBPREFIX}core$${LIB_EXTENSION} \
        ..$${LIBPREFIX}ui$${LIB_EXTENSION} \
        ..$${LIBPREFIX}wxdr$${LIB_EXTENSION} \
        ..$${LIBPREFIX}figurelib$${LIB_EXTENSION} \
        ..$${LIBPREFIX}jobs$${LIB_EXTENSION} \
        $${LIBSPATH}$${LIBPREFIX}ligand$${LIB_EXTENSION} \
        $${LIBSPATH}$${LIBPREFIX}map$${LIB_EXTENSION} \
        $${LIBSPATH}$${LIBPREFIX}molopt$${LIB_EXTENSION} \
        $${LIBSPATH}$${LIBPREFIX}conflib$${LIB_EXTENSION} \
        $${LIBSPATH}$${LIBPREFIX}chemlib$${LIB_EXTENSION} \
        $${LIBSPATH}$${LIBPREFIX}miopengl$${LIB_EXTENSION} \
        $${LIBSPATH}$${LIBPREFIX}mimath$${LIB_EXTENSION} \
        $${LIBSPATH}$${LIBPREFIX}miutil$${LIB_EXTENSION} \
        $${LIBSPATH}$${LIBPREFIX}jacgrid$${LIB_EXTENSION} \
        $${LIBSPATH}$${LIBPREFIX}umtz$${LIB_EXTENSION}

MYLIBS=-L../ -lui -lfigurelib -ljobs -lwxdr -lcore -L../../../libs -lligand -lmap -lmolopt -lconflib -lchemlib -lmiopengl -lmimath -lmiutil -ljacgrid -lumtz
LIBS =  $$MYLIBS $$MYLIBS $$LIBS

unix:!mac{
  QMAKE_LFLAGS += -Wl,--rpath=\\\$\$ORIGIN
  QMAKE_LFLAGS += -Wl,--rpath=\\\$\$ORIGIN/libs
  QMAKE_RPATH=
}

mac {
  QMAKE_LFLAGS+=-Wl,-search_paths_first
}

SOURCES += *.cpp
