import sys, os, os.path, getopt, subprocess, pickle
from PyQt4 import QtCore, QtGui, uic
import mifit

def Usage():
    print "Usage: %s [options]" % sys.argv[0]
    print "Options are:"
    print "  -w, --workdir=DIR      working dir to use"
    print ""
    print "  -?, --help              for this help message"
    print ""


if __name__ == '__main__':

    app = QtGui.QApplication(sys.argv)

    config = {}
    settings = QtCore.QSettings("MIFit", "MIExpert")
    appSettings = settings.value("restraints").toString()
    if not appSettings.isEmpty():
        config = pickle.loads(str(appSettings))

    if config.has_key('workdir'):
        workingdir = config['workdir']
    else:
        workingdir = os.getcwd()
    refprogram = None
    optlist, args = getopt.getopt(sys.argv[1:], '?', [ 'help' ])
    for opt in optlist:
        arg = opt[0]
        if arg == '-?' or arg == '--help':
            Usage()
            exit(-1)
    
        
    filename = QtGui.QFileDialog.getOpenFileName(None, "Choose a PDB file",
        workingdir, "PDB files (*.pdb);;All files (*.*)");

    dir = QtCore.QDir.toNativeSeparators(QtCore.QFileInfo(filename).absolutePath())
    filename = QtCore.QDir.toNativeSeparators(filename)

    config['workdir'] = str(dir)
    settings.setValue("restraints", pickle.dumps(config))
    
    miexpert = os.path.join(os.path.dirname(sys.argv[0]), "MIExpert.py")
    restraintArgs = [ sys.executable, miexpert, "restraints" ]
    restraintArgs += [ "--pdbfile", str(filename) ]
    restraintArgs += [ "--workdir", str(dir) ]
    result = subprocess.call(restraintArgs)
    exit(result)
