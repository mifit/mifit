#!/usr/bin/python

from util import *
import os
import sys

version = "9.b1"

qtDir = os.getenv('QTDIR')

def NameToolCommand(libname,fname):
    return "install_name_tool -change " + libname + ".framework/Versions/4/" + libname + " @executable_path/../Frameworks/" + libname + ".framework/Versions/4/" + libname + " " + fname

def makeDist():
  distName = "MIFit-" + version

  if Windows:
    osDistDir = 'dist/Windows'
    mkdir("dist")
    mkdir("dist/Windows")
    
    requiredDlls = [
      "QtCore4.dll", "QtGui4.dll", "QtOpenGL4.dll", "QtXml4.dll",
#      "QtNetwork4.dll",
    ]
    missingDlls = []
    for dll in requiredDlls:
      if not os.path.exists(dll):
        missingDlls += [dll]

    if len(missingDlls) > 0:
      print "Unable to build distribution; missing dlls: " + ", ".join(missingDlls)
      return
    mkdir("dist")
    mkdir(osDistDir)
    distDir = os.path.join(osDistDir, 'MIFit')
    mkdir(distDir)

    fileList = [ 'MIFitManual.pdf', 'license.txt' ]
    dirs = [ 'data', 'examples' ]
    fileList += [ '*.dll', "MIFit.exe" ]
    dirs += [r'C:\Program Files\Microsoft Visual Studio 8\VC\redist\x86\Microsoft.VC80.CRT']
    dirs += [r'MIExpert']
    for f in fileList:
      copy(f, distDir)
    for d in dirs:
      copy(d, distDir)

    # qt plugins, too
    # create qt.conf file to tell qt to in app dir for plugins
    conf=open(distDir + "/qt.conf","w")
    conf.write("[Paths]\n")
    conf.write("Plugins=plugins\n")
    conf.close()

    # do plugins, too
    qtDir = os.getenv('QTDIR')
    if qtDir==None:
        qtDir="c:\Qt\4.3.4"
    pluginDir = qtDir + "\plugins"
    try:
        execute("rmdir /s /q dist\Windows\MIFit\plugins")
    except:
        pass
    execute(r"mkdir dist\Windows\MIFit\plugins\accessible")
    execute(r"mkdir dist\Windows\MIFit\plugins\imageformats")
    execute(r"copy " + pluginDir + r"\accessible\*.dll dist\Windows\MIFit\plugins\accessible")
    execute(r"copy " + pluginDir + r"\imageformats\*.dll dist\Windows\MIFit\plugins\imageformats")
    execute(r"del dist\Windows\MIFit\Qt*d4.dll") # remove debug versions
    execute(r"del dist\Windows\MIFit\plugins\imageformats\qsvg*") # not used, introduces QtSql and QtXML dependency
    execute(r"del dist\Windows\MIFit\plugins\imageformats\*d4.dll") # remove debug versions
    execute(r"del dist\Windows\MIFit\plugins\accessible\*d4.dll") # remove debug versions
    execute(r"del dist\Windows\MIFit\boost*d-1_36.dll") # remove debug versions

    dir = os.getcwd()
    chdir(osDistDir)
    zip("MIFit", distName + ".zip")
    chdir(dir)
    
    # TODO fix creation of Windows installer
    copy('MIFit.nsi', osDistDir)
    copy('license.txt', osDistDir)
    chdir(osDistDir)
    nsisVersion = version
    cmd = r'"c:\Program Files\NSIS\makensis.exe" /DVERSION=' + nsisVersion + ' MIFit.nsi'
    execute(cmd)
    chdir(dir)

  elif Linux:
    osDistDir = 'dist/Linux'

    requiredSos = [
      "libQtCore.so.4", "libQtGui.so.4", "libQtOpenGL.so.4", "libQtXml.so.4",
      "libpython2.6.so.1.0",
#      "libQtNetwork.so.4",
    ]
    missingSos = []
    for so in requiredSos:
      if not os.path.exists(so):
        missingSos += [so]

    if len(missingSos) > 0:
      print "Unable to build distribution; missing shared libraries: " + ", ".join(missingSos)
      return
    mkdir("dist")
    mkdir(osDistDir)
    distDir = os.path.join(osDistDir, 'MIFit')
    mkdir(distDir)

    fileList = [ 'MIFitManual.pdf', 'license.txt' ]
    dirs = [ 'data', 'examples' ]
    fileList += [ "MIFit", "MIFit_log" ]
    dirs += [r'MIExpert']
    for f in fileList:
      copy(f, distDir)
    for d in dirs:
      copy(d, distDir)

    #libraries
    libsDir=distDir + "/libs"
    mkdir(libsDir)
    fileList = [ '*.so.*' ]
    for f in fileList:
      copy(f, libsDir)
    
    # qt plugins, too
    # create qt.conf file to tell qt to in app dir for plugins
    conf=open(distDir + "/qt.conf","w")
    conf.write("[Paths]\n")
    conf.write("Plugins=plugins\n")
    conf.close()

    # do plugins, too
    qtDir = os.getenv('QTDIR')
    if qtDir==None:
        qtDir="/usr/lib/qt4"
    pluginDir = qtDir + "/plugins"
    execute("rm -rf dist/Linux/MIFit/plugins")
    execute("mkdir dist/Linux/MIFit/plugins")
    execute("cp -R " + pluginDir + "/imageformats dist/Linux/MIFit/plugins")
    execute("cp -R " + pluginDir + "/inputmethods dist/Linux/MIFit/plugins")
    execute("rm -rf dist/Linux/MIFit/plugins/imageformats/libqsvg*") # not used, introduces QtSql and QtXML dependency

    cmd = [ 'tar', 'zcvf', "../" + distName + ".tar.gz", "MIFit" ]
    dir = os.getcwd()
    chdir(osDistDir)
    execute(cmd)
    chdir(dir)

  elif Darwin:
    osDistDir = 'dist/darwin'
    mkdir("dist")
    mkdir(osDistDir)
    distDir = os.path.join(osDistDir, 'MIFit')
    mkdir(distDir)

    fileList = ['license.txt' ]
    dirs = []
    fileList += [ "MIFit.app", "MIFit_log" ]

    for f in fileList:
      copy(f, distDir)
    for d in dirs:
      copy(d, distDir)

    # install Qt libraries, etc., into Frameworks in app bundle
    execute("rm -rf dist/darwin/MIFit/MIFit.app/Contents/Frameworks")
    execute("mkdir dist/darwin/MIFit/MIFit.app/Contents/Frameworks")
    execute("cp -R /Library/Frameworks/QtCore.framework dist/darwin/MIFit/MIFit.app/Contents/Frameworks")
    execute("cp -R /Library/Frameworks/QtGui.framework dist/darwin/MIFit/MIFit.app/Contents/Frameworks")
    execute("cp -R /Library/Frameworks/QtOpenGL.framework dist/darwin/MIFit/MIFit.app/Contents/Frameworks")

    # tell libraries their new names
    execute("install_name_tool -id @executable_path/../Frameworks/QtCore.framework/Versions/4/QtCore dist/darwin/MIFit/MIFit.app/Contents/Frameworks/QtCore.framework/Versions/4/QtCore")
    execute("install_name_tool -id @executable_path/../Frameworks/QtGui.framework/Versions/4/QtGui dist/darwin/MIFit/MIFit.app/Contents/Frameworks/QtGui.framework/Versions/4/QtGui")
    execute("install_name_tool -id @executable_path/../Frameworks/QtOpenGL.framework/Versions/4/QtOpenGL dist/darwin/MIFit/MIFit.app/Contents/Frameworks/QtOpenGL.framework/Versions/4/QtOpenGL")
    
    # tell MIFit the new library names
    execute(NameToolCommand("QtCore","dist/darwin/MIFit/MIFit.app/Contents/MacOs/MIFit"))
    execute(NameToolCommand("QtGui","dist/darwin/MIFit/MIFit.app/Contents/MacOs/MIFit"))
    execute(NameToolCommand("QtOpenGL","dist/darwin/MIFit/MIFit.app/Contents/MacOs/MIFit"))

    # tell libraries new names for their dependent libraries
    execute(NameToolCommand("QtCore","dist/darwin/MIFit/MIFit.app/Contents/Frameworks/QtGui.framework/Versions/4/QtGui"))
    execute(NameToolCommand("QtCore","dist/darwin/MIFit/MIFit.app/Contents/Frameworks/QtOpenGL.framework/Versions/4/QtOpenGL"))
    execute(NameToolCommand("QtGui","dist/darwin/MIFit/MIFit.app/Contents/Frameworks/QtOpenGL.framework/Versions/4/QtOpenGL"))

    # do plugins, too
    execute("rm -rf dist/darwin/MIFit/MIFit.app/Contents/plugins")
    execute("mkdir dist/darwin/MIFit/MIFit.app/Contents/plugins")
    execute("cp -R /Developer/Applications/Qt/plugins/imageformats dist/darwin/MIFit/MIFit.app/Contents/plugins")
    execute("cp -R /Developer/Applications/Qt/plugins/accessible dist/darwin/MIFit/MIFit.app/Contents/plugins")
    execute("rm -rf dist/darwin/MIFit/MIFit.app/Contents/plugins/imageformats/libqsvg.dylib") # not used, introduces QtSql and QtXML dependency
    execute("rm -rf dist/darwin/MIFit/MIFit.app/Contents/plugins/accessible/libqtaccessiblecompatwidgets.dylib") # not used
    
    imgplugins=["libqgif.dylib","libqico.dylib","libqjpeg.dylib","libqmng.dylib","libqtiff.dylib"]
    accplugins=["libqtaccessiblewidgets.dylib"]
    for p in imgplugins:
      pi="dist/darwin/MIFit/MIFit.app/Contents/plugins/imageformats/" + p
      execute(NameToolCommand("QtCore",pi))
      execute(NameToolCommand("QtGui",pi))
    for p in accplugins:
      pi="dist/darwin/MIFit/MIFit.app/Contents/plugins/accessible/" + p
      execute(NameToolCommand("QtCore",pi))
      execute(NameToolCommand("QtGui",pi))

    # create qt.conf file to tell qt to look in bundle for plugins
    conf=open("dist/darwin/MIFit/MIFit.app/Contents/Resources/qt.conf","w")
    conf.write("[Paths]\n")
    conf.write("Plugins=plugins\n")
    conf.close()

    # make disk image of bundle
    execute("cd dist/darwin/MIFit/; rm *.dmg*; hdiutil create -fs HFS+ -format UDZO -volname " + distName + " -srcfolder . " + distName)


def mifitVersion():
  mifitVersionFile = "apps/MIFit/core/Version.cpp"
  if not os.path.exists(mifitVersionFile):
    f = open(mifitVersionFile, "w")
    assert f, "unable to open file '" + mifitVersionFile + "'"
    f.write('char MIFit_version[] = "' + version + '";\n')
    f.close()

def main(argv=None):
  if argv is None:
    argv = sys.argv

  release = False
  help = False
  qmake = False
  mifit = False
  clean = False
  dist = False
  external = None
  
  for arg in argv[1:]:
    value = None
    if arg.find('=') >= 0:
      arg, value = arg.split('=')
    if arg == 'release':
      release = True
    elif arg == 'qmake':
      qmake = True
    elif arg == 'external':
      external = os.path.abspath(value)
    elif arg == 'clean':
      clean = True
    elif arg == 'mifit':
      mifit = True
    elif arg == 'dist':
      dist = True
    elif arg == "-h" or arg == "/?" or arg == "help" or arg == "-help" or arg == "--help":
      help = True
    else:
      print "Unrecognized argument: " + arg
      help = True
  
  if not mifit and not clean:
    mifit = True
    
  if not qmake:
    if Windows and not (os.path.exists('mi.sln') or os.path.exists('Makefile')):
      qmake = True
    elif (Linux or Darwin) and not os.path.exists('Makefile'):
      qmake = True

  if help:
    print "build [options]"
    print "  Options: (status)"
    print "    release      Build release rather than debug (" + str(release) + ")"
    print "    qmake        Regenerate platform build files (" + str(qmake) + ")"
    print "    external=dir Location of external libraries (" + str(external) + ")"
    print "    mifit        Build MIFit (" + str(mifit) + ")"
    print "    clean        Clean build files (" + str(clean) + ")"
    print "    dist         Create distributions (" + str(dist) + ")"
    return

  mifitVersion()

  config = release and 'release' or 'debug'

  execute('qmake -set mifitDir "' + os.getcwd() + '"')
  execute('qmake -set config ' + config)
  if external:
    execute('qmake -set externalDir "' + external + '"')
  else:
    external = "external"
  
  if not os.path.exists('apps/MIFit/python/sipAPIPythonEngine.h'):
    import PyQt4.pyqtconfig
    config = PyQt4.pyqtconfig.Configuration()
    execute([config.sip_bin] + config.pyqt_sip_flags.split() + [
        "-I", config.pyqt_sip_dir, "-c", "apps/MIFit/python",
        "apps/MIFit/python/*.sip"])

  if Windows:
    if qmake:
      cmd = "qmake.exe -r MI.pro"
      execute(cmd)
    if clean:
      execute("qmake.exe -r MI.pro")
      execute("nmake clean")
    if mifit:

      qtDir = os.getenv('QTDIR')
      if qtDir==None:
        qtDir="c:\Qt\4.4.1"
      if release:
          dlls = {
              'QtCore4.dll': qtDir + "/bin/QtCore4.dll",
              'QtGui4.dll': qtDir + "/bin/QtGui4.dll",
              'QtOpenGL4.dll': qtDir + "/bin/QtOpenGL4.dll",
              'QtXml.dll': qtDir + "/bin/QtXml4.dll",
              'boost_python-vc90-mt-1_36.dll': "external/boost_1_36_0/stage/lib/boost_python-vc90-mt-1_36.dll",
#              'QtNetwork4.dll': qtDir + "/bin/QtNetwork4.dll",
              }
      else:
          dlls = {
              'QtCored4.dll': qtDir + "/bin/QtCored4.dll",
              'QtGuid4.dll': qtDir + "/bin/QtGuid4.dll",
              'QtOpenGLd4.dll': qtDir + "/bin/QtOpenGLd4.dll",
              'QtXmld4.dll': qtDir + "/bin/QtXmld4.dll",
#              'QtNetworkd4.dll': qtDir + "/bin/QtNetworkd4.dll",
              }
      for k, v in dlls.iteritems():
        if not os.path.exists(k):
          copy(v, ".")
      execute('nmake')
  
  elif Linux or Darwin:
    if qmake:
      cmd = "qmake -r"
      execute(cmd)
    if clean:
      execute('make clean')
    if mifit:
      copyIfNeeded(".so", external, ".")
      if Linux:
        qtDir = os.getenv('QTDIR')
        pythonDir = "python"
        if qtDir==None:
            qtDir="/usr"
        dlls = {
          'libQtCore.so.4': qtDir + "/lib/libQtCore.so.4",
          'libQtGui.so.4': qtDir + "/lib/libQtGui.so.4",
          'libQtOpenGL.so.4': qtDir + "/lib/libQtOpenGL.so.4",
          'libQtXml.so.4': qtDir + "/lib/libQtXml.so.4",
#           'libQtNetwork.so.4': qtDir + "/lib/libQtNetwork.so.4",
          'libpython2.6.so.1.0': pythonDir + "/lib/libpython2.6.so.1.0",
          }
      elif Darwin:
        dlls = {}
      for k,v in dlls.iteritems():
        if not os.path.exists(k):
          copy(v, ".")
      execute('make')
      if Darwin:
        # create Resources bundle and links to resources
        # links will be followed when building distribution
        if not os.path.exists("MIFit.app/Contents/Resources"):
          execute("mkdir MIFit.app/Contents/Resources")
        if not os.path.exists("MIFit.app/Contents/Resources/MIExpert"):
          execute("cd MIFit.app/Contents/Resources; ln -s ../../../MIExpert MIExpert")
        if not os.path.exists("MIFit.app/Contents/Resources/data"):
          execute("cd MIFit.app/Contents/Resources; ln -s ../../../data data")
        if not os.path.exists("MIFit.app/Contents/Resources/examples"):
          execute("cd MIFit.app/Contents/Resources; ln -s ../../../examples examples")
        if not os.path.exists("MIFit.app/Contents/Resources/MIFitManual.pdf"):
          execute("cd MIFit.app/Contents/Resources; ln -s ../../../MIFitManual.pdf MIFitManual.pdf")
  
  if dist:
    makeDist()

if __name__ == "__main__":
  main()
