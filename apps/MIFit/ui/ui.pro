
include(../../../common.pri)

TEMPLATE = lib
SOURCES = *.cpp
HEADERS = *.h
FORMS = *.ui
DESTDIR = ..
INCLUDEPATH += .
CONFIG += no_keywords

!win32 {
  # If using shared libraries, set the run-time library path.
  !staticlib{ 
    # The $${DOLLAR} represents a $ after being interpreted by qmake.
    DOLLAR=$
    # The command passes $$ORIGIN/source and $$ORIGIN../../libs on to the linker.
    QMAKE_RPATH =-Wl,-rpath,'$${DOLLAR}$${DOLLAR}ORIGIN/source',-rpath,'$${DOLLAR}$${DOLLAR}ORIGIN/../../libs',-rpath,
  }

}
