
unix:!mac{
  QMAKE_LFLAGS += -Wl,--rpath=\\\$\$ORIGIN
  QMAKE_LFLAGS += -Wl,--rpath=\\\$\$ORIGIN/lib
  QMAKE_LFLAGS += -Wl,--rpath=\\\$\$ORIGIN/libs
  QMAKE_RPATH=
}
