import os
from PyQt4 import QtCore, QtGui, QtNetwork

__doc__ = """MIFit interaction module

Defines functions for interfacing with a MIFit session.
The environment must have MIFIT_SOCKET_ID for an
associated MIFit session. This will be the case for
scripts run from the MIFit Jobs menu.
"""

mifitDir = None
shelxDir = None

if 'MIFIT_DIR' in os.environ.keys():
    mifitDir = os.environ['MIFIT_DIR']

if 'SHELX_DIR' in os.environ.keys():
    shelxDir = os.environ['SHELX_DIR']


def exec_script(script):
    """Executes the given script in the associated MIFit session"""
    socketId = os.environ['MIFIT_SOCKET_ID']
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
        print 'error connecting to local socket', sys.argv[1]
    sock = None
    return result

def version():
    """Returns the version of MIFit"""
    return exec_script("mifit.version()")

def directory():
    """Returns the installation directory for MIFit"""
    return exec_script("mifit.directory()")

def scriptPort():
    """Returns the identifier for the local sockect port used for interaction with MIFit"""
    return os.environ['MIFIT_SOCKET_ID']

def setJobWorkDir(dir):
    """Sets the working directory for the current job (when the script is run from MIFit)"""
    script = "mifit.setJobWorkDir('" + os.environ['MIFIT_JOB_ID'] + "', '" + dir + "')"
    return exec_script(script)

def writeCurrentModel(file):
    """Causes MIFit to output the current model to the given file"""
    return exec_script("mifit.writeCurrentModel('" + file + "')")

def dictionaryResidueList():
    """Returns the dictionary residue list from MIFit"""
    return exec_script("mifit.dictionaryResidueList()")

def spacegroupList():
    """Returns the spacegroup list from MIFit"""
    return exec_script("mifit.spacegroupList()")

