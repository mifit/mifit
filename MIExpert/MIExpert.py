#!/usr/bin/env python
#####################################################################
#                                                                   #
# Master script for all mi_ python scripts                          #
# Copyright: Molecular Images   2007                                #
#                                                                   #
# This script is distributed under the same conditions as MIFit     #
#                                                                   #
#####################################################################

import sys, os, time, ConfigParser, StringIO, traceback

def Usage():
    print "Usage: %s [command] [args]" % sys.argv[0]
    print "Command may be:"
    print "  bng [options]"                # long
    print "  convertlib [options]"         # long
    print "  dataprep [options]"           # long
    print "  deposit3d"                    # long
    print "  integrate [options]"          # long
    print "  integrate_mosflm [options]"   # long
    print "  ligandoverlap [options]"      # long
    print "  molrep [options]"             # long
    print "  ncsmodeler [options]"         # long
    print "  refine [options]"             # long
    print "  restraints [options]"         # long
    print "  sadphase [inputfile]"         # long
    print "  --version"
    print ""
    print "Run \"%s command --help\" for help on that command" % sys.argv[0]

def Process():
    cmd=sys.argv[1]
    sys.argv=sys.argv[1:] # shift arguments

    if cmd=="bng":
        import mi_bng
        return mi_bng.Run(sys.argv)
    elif cmd=="convertlib":
        import mi_convertlib
        return mi_convertlib.Run(sys.argv)
    elif cmd=="dataprep":
        import mi_dataprep
        return mi_dataprep.Run(sys.argv)
    elif cmd=="deposit3d":
        import mi_deposit3d
        return mi_deposit3d.Run(sys.argv)
    elif cmd=="integrate":
        import mi_integrate
        return mi_integrate.Run(sys.argv)
    elif cmd=="integrate_mosflm":
        import mi_integrate_mosflm
        return mi_integrate_mosflm.Run(sys.argv)
    elif cmd=="ligandfit":
        import mi_ligandfit
        return mi_ligandfit.Run(sys.argv)
    elif cmd=="ligandoverlap":
        import mi_ligandoverlap
        return mi_ligandoverlap.Run(sys.argv)
    elif cmd=="molrep":
        import mi_molrep
        return mi_molrep.Run(sys.argv)
    elif cmd=="ncsmodeler":
        import mi_ncsmodeler
        return mi_ncsmodeler.Run(sys.argv)
    elif cmd=="refine":
        import mi_refine
        return mi_refine.Run(sys.argv)
    elif cmd=="restraints":
        import mi_restraints
        return mi_restraints.Run(sys.argv)
    elif cmd=="sadphase":
        import mi_sadphase
        return mi_sadphase.Run(sys.argv)
    else:
        print "Unknown command"
        return 1
    return 1

def printVersion():
  versionString = "MIExpert %s" % version
  print versionString

def Run():
   retval=1
   if len(sys.argv) < 3:
       Usage()
       return
     
   for i in range(0, len(sys.argv)):
     if sys.argv[i] == "--version":
       printVersion()
       sys.exit(retval)
     elif sys.argv[i] == "--output":
       sys.stdout = open(sys.argv[i+1], 'w', 0)
       sys.stderr = sys.stdout
       del sys.argv[i:i+2]
       break
   try:
       retval=Process()
   except:
       print "An error occured processing the request"
       traceback.print_exc(file=sys.stdout)
       if Windows:
           if os.path.exists("C:\MIEXPERT.DBG"): raise
       else:
           if os.path.exists("/tmp/miexpert.dbg"): raise
           
   sys.exit(retval)

version = "undefined"
Windows = sys.platform.find("win") > -1
Cygwin = sys.platform == "cygwin"
Linux = sys.platform == "linux2"
molimageHome = os.path.dirname(sys.executable)
if __name__ == "__main__":
   sys.exit(Run())
