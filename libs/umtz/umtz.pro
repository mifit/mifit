
include(../../common.pri)

TEMPLATE = lib
CONFIG += staticlib
CONFIG -= qt moc
DESTDIR=..

HEADERS = *.h
SOURCES = mmtzlib.c umtzlib.c library.c
win32{
  DEFINES += _MVS i386 WINVER=0x400
}
