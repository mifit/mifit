import os

def ccp4check():
  # set path and check the CCP4 installation
  class Ccp4Settings: pass

  settings = Ccp4Settings()
  
  env_vars = os.environ.keys()
  required_env_vars = ['CCP4', 'CCP4_SCR', 'CLIBD_MON', 'CLIBD', 'CBIN']
  for v in required_env_vars:
    if env_vars.count(v) == 0:
      return False, 'The environment variable ' + v + ' was not found - is the CCP4 environment established ?'

  settings.ccp4 = os.environ['CCP4']
  settings.scr = os.environ['CCP4_SCR']
  settings.mon = os.environ['CLIBD_MON']
  settings.libd = os.environ['CLIBD']
  settings.bin = os.environ['CBIN']

  settings.symfile = 'symop.lib'
  settings.names = 'full_names.list'

  settings.entitylist = os.path.join(settings.mon, settings.names)
  settings.symmetrylib = os.path.join(settings.libd, settings.symfile)

  if not os.path.exists(settings.entitylist):
    return False, 'The CCP4 entity list was not found'

  if not os.path.exists(settings.symmetrylib):
    return False, 'The CCP4 symmetry library was not found'

  # Check CCP4 version

  originalDir = os.getcwd()
  os.chdir(settings.scr)

  settings.version = 'none'

  file = open('mi_mtzdump.inp','w')
  file.write('HEADER\n')
  file.write('END\n')
  file.close()

  os.system('mtzdump < mi_mtzdump.inp > mi_mtzdump.log 2> mi_mtzdump_err.log')
  os.remove('mi_mtzdump.inp')

  file = open('mi_mtzdump.log','r')
  allLines = file.readlines()
  file.close()

  fileexists = os.path.exists('mi_mtzdump.log')
  if fileexists != 0:

      file = open('mi_mtzdump.log','r')
      allLines = file.readlines()
      file.close()

      os.remove('mi_mtzdump.log')
      os.remove('mi_mtzdump_err.log')   

      for eachLine in allLines:
          if eachLine.find('CCP4 6') > -1:
              settings.version = '6'

  if settings.version != '6':
      return False, 'This script requires CCP4 v6'

  os.chdir(originalDir)
  return settings,None

