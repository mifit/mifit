include(conf.pri)

TEMPLATE = subdirs
SUBDIRS += libs apps MIExpert

data.path = $$PREFIX/data
data.files = data/*

examples.path = $$PREFIX/examples
examples.files = examples/*

INSTALLS += data examples
