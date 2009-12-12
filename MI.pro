include(conf.pri)
include(common.pri)

TEMPLATE = subdirs
SUBDIRS += libs apps MIExpert

data.path = $$PREFIX/data
data.files = data/*

examples.path = $$PREFIX/examples
examples.files = examples/*

OTHER_FILES += README.txt license.txt jobs.js
other.path = $$PREFIX
other.files = $$OTHER_FILES

INSTALLS += data examples other
