import os
from PyQt4 import QtCore, QtGui, QtNetwork

__doc__ = """MIFit interaction module

Defines functions for interfacing with a MIFit session.
The environment must have MIFIT_SOCKET_ID for an
associated MIFit session. This will be the case for
scripts run from the MIFit Jobs menu.
"""

socketId = None
mifitDir = None
shelxDir = None

if 'MIFIT_SOCKET_ID' in os.environ.keys():
    socketId = os.environ['MIFIT_SOCKET_ID']

if 'MIFIT_DIR' in os.environ.keys():
    mifitDir = os.environ['MIFIT_DIR']

if 'SHELX_DIR' in os.environ.keys():
    shelxDir = os.environ['SHELX_DIR']


def exec_script(script):
    """Executes the given script in the associated MIFit session"""
    global socketId
    result = None
    sock = QtNetwork.QLocalSocket()
    sock.connectToServer(socketId)
    if sock.waitForConnected():

        stream = QtCore.QDataStream(sock)
        stream.setVersion(QtCore.QDataStream.Qt_4_0)

        scriptString = QtCore.QString(script + "\b")
        stream << scriptString
        sock.flush()

        dataSize = 0
        waitForData = True
        while waitForData:
            if not sock.waitForReadyRead(5000):
                break
            if dataSize == 0:
                if sock.bytesAvailable() < 8:
                  continue
                dataSize = stream.readInt64()
            if stream.atEnd():
                break
            if sock.bytesAvailable() < dataSize:
                continue
            waitForData = False
            result = QtCore.QString()
            stream >> result

        sock.close()
    else:
        print 'error connecting to local socket', socketId
    sock = None
    return result

def version():
    """Returns the version of MIFit"""
    return str(exec_script("mifit.version"))

def directory():
    """Returns the installation directory for MIFit"""
    return str(exec_script("mifit.directory()"))

def scriptPort():
    """Returns the identifier for the local sockect port used for interaction with MIFit"""
    global socketId
    return socketId

def setScriptPort(id):
    """Sets the identifier for the local sockect port used for interaction with MIFit.
Usually this will be obtained automatically from the MIFIT_SOCKET_ID environment variable."""
    global socketId
    socketId = id

def setJobWorkDir(dir):
    """Sets the working directory for the current job (when the script is run from MIFit)"""
    script = "mifit.setJobWorkDir('" + os.environ['MIFIT_JOB_ID'] + "', '" + dir + "')"
    return exec_script(script)

def writeCurrentModel(file):
    """Causes MIFit to output the current model to the given file"""
    return exec_script("mifit.writeCurrentModel('" + file + "')")

def dictionaryResidueList():
    """Returns the dictionary residue list from MIFit"""
    result = exec_script("mifit.dictionaryResidueList()")
    if result:
        return str(result).split(',')
    return None

def spacegroupList():
    """Returns the spacegroup list from MIFit"""
    result = exec_script("mifit.spacegroupList()")
    if result:
        return str(result).split(',')
    return None

def addJob(menuName, jobName, executable, arguments, workingDirectory):
    """Adds a new custom job to the MIFit Job menu"""
    args = "[ '" + "', '".join(arguments) + "' ]"
    script = "mifit.addJob('%s', '%s', '%s', %s, '%s')" % (menuName, jobName, executable, args, workingDirectory)
    exec_script(script)
