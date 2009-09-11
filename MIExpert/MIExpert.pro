include(../common.pri)

# There's nothing to do for this build, we only need to copy files
pyfiles.path = $$PREFIX/MIExpert
pyfiles.files = *.py
INSTALLS += pyfiles
PRE_TARGETDEPS += install_pyfiles

# trick qmake in to not doing any "normal" processing
TEMPLATE=lib
QMAKE_LINK=echo
QMAKE_AR=echo
QMAKE_RANLIB=echo

