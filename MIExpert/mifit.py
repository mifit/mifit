import os
from PyQt4 import QtCore, QtGui, QtNetwork

def exec_script(script):
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

def setJobWorkDir(dir):
    script = "mifit.setJobWorkDir('" + os.environ['MIFIT_JOB_ID'] + "', '" + dir + "')"
    exec_script(script)
