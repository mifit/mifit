
include(../../common.pri)

TEMPLATE = lib
CONFIG += staticlib
CONFIG -= qt moc
DESTDIR=..

# Input
HEADERS += jacgrid.h atom.h cell_table.h surface.h grid.h plane.h
SOURCES += grid.cpp atom.cpp isosurface.cpp surface.cpp plane.cpp
