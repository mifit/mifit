import sys, os, os.path, getopt, subprocess, pickle
from PyQt4 import QtCore, QtGui, uic


def Usage():
    print "Usage: %s [options]" % sys.argv[0]
    print "Options are:"
    print "  -w, --workdir=DIR      working dir to use"
    print "  -r, --refprogram=PROG  shelx (default), cns, or cnx"
    print ""
    print "  -?, --help              for this help message"
    print ""


if __name__ == '__main__':

    app = QtGui.QApplication(sys.argv)

    config = {}
    settings = QtCore.QSettings("MIFit", "MIExpert")
    appSettings = settings.value("convertlib").toString()
    if not appSettings.isEmpty():
        config = pickle.loads(str(appSettings))

    if config.has_key('workdir'):
        workingdir = config['workdir']
    else:
        workingdir = os.getcwd()
    refprogram = None
    optlist, args = getopt.getopt(sys.argv[1:], 'r:?', [ 'refprogram=', 'help' ])
    for opt in optlist:
        arg = opt[0]
        if arg == '-?' or arg == '--help':
            Usage()
            exit(-1)
        elif arg == '-r' or arg =='--refprogram':
            refprogram = opt[1]
    
    if refprogram is None:
        print "--refprogram must be specified"
        Usage()
        exit(-1)
        
    filename = QtGui.QFileDialog.getOpenFileName(None, "Choose a CIF file",
        workingdir, "CIF Files (*.cif);;All files (*.*)");

    dir = QtCore.QDir.toNativeSeparators(QtCore.QFileInfo(filename).absolutePath())
    filename = QtCore.QDir.toNativeSeparators(filename)

    config['workdir'] = str(dir)
    settings.setValue("convertlib", pickle.dumps(config))
    
    miexpert = os.path.join(os.path.dirname(sys.argv[0]), "MIExpert.py")
    convertArgs = [ sys.executable, miexpert, "convertlib" ]
    convertArgs += [ "--workdir=" + str(dir) ]
    convertArgs += [ "--refprogram=" + refprogram ]
    convertArgs += [ "--cif=" + str(filename) ]
    result = subprocess.call(convertArgs)
    exit(result)
