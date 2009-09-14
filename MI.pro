include(conf.pri)
include(common.pri)

TEMPLATE = subdirs
SUBDIRS += libs apps MIExpert

data.path = $$PREFIX/data
data.files = data/*

examples.path = $$PREFIX/examples
examples.files = examples/*

OTHER_FILES += README.txt
other.path = $$PREFIX/README.txt
other.files = $$OTHER_FILES

INSTALLS += data examples other
