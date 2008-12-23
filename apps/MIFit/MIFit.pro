
include(../../common.pri)

TEMPLATE = subdirs

USE_ASPLOT {
  FIGURELIB = figurelib
}
MI_USE_JOBS{
  JOBSLIB = jobs
}

SUBDIRS += core $$JOBSLIB $$FIGURELIB wxdr ui python main
