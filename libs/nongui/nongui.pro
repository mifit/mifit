include(../libs.pri)

TARGET = $$qtLibraryTarget(nongui)
TEMPLATE = lib
SOURCES = $$files(*.cpp)
HEADERS = $$files(*.h)
