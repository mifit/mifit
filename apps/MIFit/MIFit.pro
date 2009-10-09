include(MIFit.pri)

TEMPLATE = subdirs

SUBDIRS += part_core part_script part_jobs part_figurelib part_ui part_main

part_core.subdir = core

part_script.subdir = script
part_script.depend = core

part_jobs.subdir = jobs
part_jobs.depends = core

part_figurelib.subdir = figurelib
part_figurelib.depends = core

part_ui.subdir = ui
part_ui.depends = core

part_main.subdir = main
part_main.depends += part_core part_jobs part_figurelib part_ui
