import os
import sys
import types
import re

# Flags identifying platform type
Linux = sys.platform.find('linux') > -1
Darwin = sys.platform ==  'darwin'
Windows = not Darwin and sys.platform.find('win') > -1
Cygwin = sys.platform == 'cygwin'

def mkdir(d):
  if not os.path.exists(d):
    os.mkdir(d)

def chdir(d):
  os.chdir(d)

def execute(*args):
  cmd = []
  for a in args:
    if isinstance(a, types.ListType):
      cmd += a
    else:
      cmd += [a]
  cmd = " ".join(cmd)
  print(cmd)
  if not globals().has_key('testing') or not globals()['testing']:
    result = os.system(cmd)
    if result != 0:
      assert False, "command failed: " + cmd

def copy(src, dest):
  if Windows:
    src = src.replace("/", "\\")
    dest = dest.replace("/", "\\")
    srcre = re.compile(".*\\\\")
    srcname = srcre.sub('', src)
    if os.path.isdir(src):
      dest = os.path.join(dest, srcname)
      mkdir(dest)
      command = 'xcopy "' + src + '" "' + dest + '" /S /Y'
    else:
      command = 'copy "' + src + '" "' + dest + '"'
    
  elif Linux:
    command = "cp -rH " + src + " " + dest
  else:
    command = "cp -r " + src + " " + dest
  execute(command)

def copyIfNeeded(extension, fromDir, toDir):
  existingFiles = {}
  extIndex = -len(extension)
  for name in os.listdir(toDir):
    if name[extIndex:] == extension:
      existingFiles[name] = True
    
  for name in os.listdir(fromDir):
    if name[extIndex:] == extension:
      if not existingFiles.has_key(name):
        copy(fromDir + "/" + name, toDir)

def zip(dir, file):
  import zipfile
  
  z = zipfile.ZipFile(file, mode="w", compression=zipfile.ZIP_DEFLATED)
  
  cwd = os.getcwd()
  os.chdir(dir)
  try:
    for dirpath,dirs,files in os.walk(''):
      for a_file in files:
        a_path = os.path.join(dirpath, a_file)
        z.write(a_path, os.path.join(dir, a_path))
    z.close()
  finally:
    os.chdir(cwd)
