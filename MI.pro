exists(conf.pri) {
    include(conf.pri)
}
include(common.pri)

TEMPLATE = subdirs
SUBDIRS += libs apps

data.path = $$PREFIX/data
data.files = data/*

examples.path = $$PREFIX/examples
examples.files = examples/*

OTHER_FILES += README.txt license.txt
other.path = $$PREFIX
other.files = $$OTHER_FILES

INSTALLS += data examples other
