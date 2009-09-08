include(../libs.pri)

TEMPLATE = lib
CONFIG -= qt moc

HEADERS = *.h
SOURCES = mmtzlib.c umtzlib.c library.c
win32{
  DEFINES += _MVS i386 WINVER=0x400
}
