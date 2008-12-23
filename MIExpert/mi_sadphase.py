#####################################################################
# Script: mi_sadphase.py                                            #
# Release: All                                                      #
#                                                                   #
# SAD phasing script using SHELX and CCP4                           # 
#                                                                   #
# Copyright: Molecular Images   2006                                #
#                                                                   #
# This script is distributed under the same conditions as MIFit     #
#                                                                   #
# Compatible with CCP4 6.0.1                                        # 
#                                                                   #
#####################################################################

import sys
import os
import time
import string
import math
import getopt
import ccp4check

def Usage():
    print "Usage: %s [options]" % sys.argv[0]
    print "Options are:"
    print "  -d,--workdir=DIR            The working directory"
    print "     --saddatafile=FILE       Reflection data file"
    print "     --sitefile=FILE          The site file. Default: none, determine with shelxd "
    print "     --sitenumber=NUM         Number of expected anomalous scattering sites. Default: 0"
    print "     --scatterer=ELEM         The scattering element. Default: S"
    print "     --solventfraction=NUM    The solvent fraction.  No default"
    print "     --separation=NUM         The minimum separation distance. Default: 3.5"
    print "     --bothhands=yes or no    Try both hands? Default=yes."
    print "     --ssnumber=NUM           Number SS bridges. Default=0"
    print "     --spacegroup_no=NUM      Specify a SG number. Default: deduce from data file"
    print "     --siterefinemethod=type  One of bp3 or mlphare. Default: mlphare."
    print "  -s,--shelx_dir=DIR          Path to shelx executables. Default $SHELXBIN."
    print "  -h,--molimagehome           Path to MIFit."
    print "     --writemap=yes or no     Write ccp4 maps. Default: no"
    print "  -?,--help                   This help file"

def Run(argv=None):
    if argv is None:
        argv=sys.argv

    # Backdoor to hardcoded path to MIFit

    mifit_root = 'none'

    # Initialize

    quote = """'"""

    inputfile = 'mi_runsadphase.txt'
    shelx_basename = 'mi_sad_'
    shelxd = 'none'
    do_mlw_files = 'true'

    workingdir = 'none'
    datafile = 'none'
    sitefile = 'none'
    number_sites = 'none'
    solvent_fraction = 'none'
    scatterer_in = 'S'
    separation = 'none'
    run_invert = 'yes'
    number_ss = 'none'
    shelx_directory = 'none'
    write_map = 'no'

    wavelength = '2.2909'
    auto_site_find = 'no'
    datatype = 'none'
    acell = 'none'
    bcell = 'none'
    ccell = 'none'
    alpha = 'none'
    beta = 'none'
    gamma = 'none'
    spacegroup = 'none'
    spacegroup_input = 'none'
    spacegroup_no = 'none'
    spacegroup_lattice = 'none'

    cc = 'none'
    cc_weak = 'none'

    refine_method = 'mlphare'

    filename_local_in = 'mi_sad.hkl'
    filename_mtz = 'mi_sad.mtz'
    filename_out_mtz = 'mi_sad_sort.mtz'
    filename_out_mtz_f = 'mi_sad_sort_f.mtz'

    html_site_summary = 'mi_site_summary.html'
    html_phase_summary = 'mi_phase_summary.html'

    aLine = []
    symList = []
    atomList = []
    atomList_i = []
    atomList_elements = []

    aList_res_files = []
    aList_cc_weak = []
    aList_res = []
    aList_cc = []

    aList_resolution_sites = []
    aList_numberrefs_sites = []
    aList_fom_shell = []

    aList_atom_out_element = []
    aList_atom_out_x = []
    aList_atom_out_y = []
    aList_atom_out_z = []
    aList_atom_out_occ = []
    aList_atom_out_b = []

    ################################################
    # Check for 3rd party software and environment #
    ################################################

    test_platform = sys.platform

    ############################
    # Set input information    #
    ############################
    args = argv[1:]
    optlist, args = getopt.getopt(
        args,'d:h:s:?',
        ["workdir=","saddatafile=","sitefile=","sitenumber=",
         "scatterer=","solventfraction=","separation=",
         "bothhands=","ssnumber=","spacegroup_no=",
         "siterefinemethod=","shelx_dir=","molimagehome=","writemap=","help"])
    number_of_inputs = len(optlist)
    if number_of_inputs==0:
        Usage()
        return
    count = 0
    while count < number_of_inputs:
        aList = optlist[count]
        number_of_list_inputs = len(aList)
        arg_value=""
        if number_of_list_inputs >=1:
            arg_value = aList[0]
            if arg_value == '-?' or arg_value=='--help':
                Usage()
                return
        if number_of_list_inputs >=2:
            param_value = aList[1]
            if arg_value=='--workdir' or arg_value=='-d':
                workingdir = param_value
            elif arg_value=='--saddatafile':
                datafile = param_value
            elif arg_value=='--sitefile':
                sitefile = param_value
            elif arg_value=='--sitenumber':
                number_sites = param_value
            elif arg_value=='--scatterer':
                scatterer_in = param_value
            elif arg_value=='--solventfraction':
                solvent_fraction = param_value
            elif arg_value=='--separation':
                separation = param_value
            elif arg_value=='--bothhands':
                run_invert = param_value
            elif arg_value=='--ssnumber':
                number_ss = param_value
            elif arg_value=='--spacegroup_no':
                spacegroup_input = param_value
            elif arg_value=='--siterefinemethod':
                refine_method = param_value
            elif arg_value=='--shelx_dir' or arg_value=='-s':
                shelx_directory = param_value
            elif arg_value=='--molimagehome' or arg_value=='-h':
                mifit_root = param_value
            elif arg_value=='--writemap':
                write_map = param_value
        count=count+1

    ccp4,error = ccp4check.ccp4check()
    if not ccp4:
      print '\n' + error + '\n'
      time.sleep(4)
      return 1
    
    # Check the inputs
    if sitefile == 'None' or sitefile == 'NONE':
        sitefile = 'none'

    if solvent_fraction == 'None' or solvent_fraction == 'NONE':
        solvent_fraction = 'none'

    if number_sites == 'None' or number_sites == 'NONE' or number_sites == '0':
        number_sites = 'none'

    if number_ss == 'None' or number_ss == 'NONE' or number_ss == '0':
        number_ss = 'none'

    if refine_method == 'None' or refine_method == 'NONE' or refine_method == 'none' or refine_method == 'BP3':
        refine_method = 'bp3'

    if refine_method == 'MLPHARE':
        refine_method = 'mlphare'

    if refine_method != 'bp3' and refine_method != 'mlphare':
        print 'The site refinement method must be bp3 or mlphare'
        time.sleep(4)
        return 1

    fileexists = os.path.exists(datafile)
    if fileexists == 0:
        print 'The SAD reflection data file was not found ',datafile
        time.sleep(4)
        return 1

    fileexists = os.path.exists(workingdir)
    if fileexists == 0:
        print 'The working directory was not found ',workingdir
        time.sleep(4)
        return 1

    fileexists = os.path.exists(sitefile)
    if fileexists == 0 and sitefile != 'none':
        print 'The scatterer site file was not found ',sitefile
        time.sleep(4)
        return 1

    if solvent_fraction == 'none':
        print 'The solvent faction was not given'
        time.sleep(4)
        return 1
    else:
        float_solvent_fraction = float(solvent_fraction)
        if float_solvent_fraction < 0.1 or float_solvent_fraction > 0.90:
            print 'The solvent fraction should be a reasonable FRACTION'
            time.sleep(4)
            return 1

    if number_sites == 'none':
        print 'The number of expected anomalous scattering sites must be given'
        time.sleep(4)
        return 1

    #
    # If needed, get determine SHELX bin directory path
    #

    if sitefile == 'none':

        # Determine from environment variable (LINUX)

        find_shelx_bin = os.environ.keys().count('SHELXBIN')
        if find_shelx_bin != 0:
            shelx_dir_linux = os.environ['SHELXBIN']
            shelxd = os.path.join(shelx_dir_linux,'shelxd') 

        # Determine from input parameter

        if shelx_directory != 'none':

            fileexists = os.path.exists(shelx_directory)
            if fileexists != 0:
                if test_platform.find('win') > -1:
                    shelxd = os.path.join(shelx_directory,'shelxd.exe')
                else:
                    shelxd = os.path.join(shelx_directory,'shelxd')

        # Confirm we have SHELXD

        fileexists = os.path.exists(shelxd)
        if fileexists == 0:
            print 'The SHELXD executable was not found ',shelxd
            time.sleep(4)
            return 1


    ###################################################################################
    # Copy reflection data to local working area and collect information on file type #
    ###################################################################################

    os.chdir(workingdir)

    if datafile.find('.mtz') > -1:
        datatype = 'scala'

        file = open(datafile,'rb')
        allLines = file.readlines()
        file.close()

        file = open(filename_local_in,'wb')
        file.writelines(allLines)
        file.close()

    else:

        file = open(datafile,'r')
        allLines = file.readlines()
        file.close()

        file = open(filename_local_in,'w')

        for eachLine in allLines:

            if eachLine.find('SOURCE_WAVELENGTH') > -1:
                datatype = 'dstartrek'

            file.write(eachLine)

        file.close()

    if datatype == 'none':
        datatype = 'scalepack'

    ##############################################
    # Convert to MTZ as standard working format  #
    ##############################################

    print '\nPREPARING DATA\n'

    #
    # 1. D*TREK case
    #

    if datatype == 'dstartrek':

        # Obtain cell and spacegroup

        for eachLine in allLines:

            if eachLine.find('SOURCE_WAVELENGTH') > -1:
                aList = eachLine.split()
                wavelength = aList[2]
                wavelength = wavelength.replace(';','')

            if eachLine.find('CRYSTAL_UNIT_CELL') > -1:
                aList = eachLine.split()
                acell = aList[1]
                bcell = aList[2]
                ccell = aList[3]
                alpha = aList[4]
                beta = aList[5]
                gamma = aList[6]
                gamma = gamma.replace(';','')

            if eachLine.find('CRYSTAL_SPACEGROUP') > -1:
                aList = eachLine.split('=')
                spacegroup_no = aList[1]
                spacegroup_no = spacegroup_no.replace(';','')

        if acell == 'none' or bcell == 'none' or ccell == 'none' \
           or alpha == 'none' or 'beta' == 'none' or gamma == 'none' or spacegroup_no == 'none':
            print 'Cell or spacegroup were not found in the reflection data file'
            time.sleep(4)
            return 1

        # Run DTREK2MTZ conversion

        print 'Input file appears to be from D*TREK'

        filename_inp = 'mi_rundstartrek.inp'
        filename_log = 'mi_rundstartrek.log'

        file = open(filename_inp,'w')
        file.write('title Converted by MIFit')
        file.write('\nsymm ')
        file.write(spacegroup_no)
        file.write('\n')    
        file.write('WAVE ')
        file.write(wavelength)
        file.write('\n')
        file.write('end\n')
        file.close()

        rundtrek2mtz = 'dtrek2mtz hklin ' + filename_local_in + ' hklout ' + filename_mtz + ' < ' + filename_inp + ' > ' + filename_log

        os.system(rundtrek2mtz)

        fileexists = os.path.exists(filename_mtz)
        if fileexists != 0:
            os.remove(filename_local_in)
            os.remove(filename_inp)
            os.remove(filename_log)
        else:
            print 'Output mtz file from DTREK2MTZ was not found'
            time.sleep(4)
            return 1

    #
    # 2. SCALEPACK case 
    #

    if datatype == 'scalepack':

        # Obtain cell and spacegroup

        count = 0
        for eachLine in allLines:

            count = count + 1
            aLine = eachLine.split()

            if count == 3:
                acell = aLine[0]
                bcell = aLine[1]
                ccell = aLine[2]
                alpha = aLine[3]
                beta = aLine[4]
                gamma = aLine[5]
                spacegroup = aLine[6]

                spacegroup = spacegroup.upper()

        if acell == 'none' or bcell == 'none' or ccell == 'none' \
           or alpha == 'none' or 'beta' == 'none' or gamma == 'none' or spacegroup == 'none':
            print 'Cell or spacegroup were not found in the reflection data file'
            time.sleep(4)
            return 1

        # Run SCALEPACK2MTZ conversion

        print 'Input file appears to be from SCALEPACK'

        filename_inp = 'mi_runscalepack2mtz.inp'
        filename_log = 'mi_scalepack2mtz.log'

        file = open(filename_inp,'w')
        file.write('name proj ')
        file.write(shelx_basename)
        file.write('\nsymm ')
        file.write(spacegroup)
        file.write('\nend\n')
        file.close()

        runscalepack2mtz = 'scalepack2mtz hklin ' + filename_local_in + ' hklout ' + filename_mtz + ' < ' + filename_inp + ' > ' + filename_log

        os.system(runscalepack2mtz)

        fileexists = os.path.exists(filename_mtz)
        if fileexists != 0:
            os.remove(filename_local_in)
            os.remove(filename_inp)

            # Get spacegroup number

            file = open(filename_log,'r')
            allLines = file.readlines()
            file.close()

            for eachLine in allLines:
                if eachLine.find('* Space group =') > -1:
                    aLine = eachLine.split('number')
                    spacegroup_no = aLine[1]
                    spacegroup_no = spacegroup_no.replace(')','')
                    spacegroup_no = spacegroup_no.strip()

            os.remove(filename_log)
        else:
            print 'Output mtz file from SCALEPACK2MTZ was not found'
            time.sleep(4)
            return 1    

    #
    # 3. SCALA case (not tested)
    #

    if datatype == 'scala':

        file = open('mi_mtzdump.inp','w')
        file.write('HEADER\n')
        file.write('END\n')
        file.close()

        runmtz = 'mtzdump HKLIN ' + filename_local_in + ' < mi_mtzdump.inp > mi_mtzdump.log'
        os.system(runmtz)

        file = open('mi_mtzdump.log','r')
        allLines = file.readlines()
        file.close()

        os.remove('mi_mtzdump.log')
        os.remove('mi_mtzdump.inp')

        read_columns = 'no'
        read_labels = 'no'
        read_cell = 'no'
        imean = 'no'
        sigimean = 'no'
        iplus = 'no'
        sigiplus = 'no'
        iminus = 'no'
        sigiminus = 'no'

        # Parse for cell, spacegroup and check for standard data types/names

        for eachLine in allLines:

            line_length = len(eachLine)

            if read_columns == 'yes' and line_length > 1:
                colList = eachLine.split()
                read_columns = 'no'

            if read_labels == 'yes' and line_length > 1:
                labelList = eachLine.split()
                read_labels = 'no'

            if eachLine.find('* Column Labels :') > -1:
                read_columns = 'yes'

            if eachLine.find('* Column Types :') > -1:
                read_labels = 'yes'

            if eachLine.find('* Space group') > -1:
                aLine = eachLine.split('number')
                spacegroup_no = aLine[1]
                spacegroup_no = spacegroup_no.replace(')','')
                spacegroup_no = spacegroup_no.strip()

            if eachLine.find('* Cell dimensions ') > -1:
                read_cell = 'yes'

            if readcell == 'yes' and line_length > 1:
                aLine = eachLine.split()
                number_items = len(aLine)

                if number_items == 6:
                    acell = aLine[0]
                    bcell = aLine[1]
                    ccell = aLine[2]
                    alpha = aLine[3]
                    beta = aLine[4]
                    gamma = aLine[5]

                    readcell = 'no'

        # check standard data

        list_length = len(labelList)
        count = 0

        while count < list_length:
            if colList[count] == 'IMEAN':
                imean = 'yes'

            if colList[count] == 'SIGIMEAN':
                sigimean = 'yes'

            if colList[count] == 'I(+)':
                iplus = 'yes'

            if colList[count] == 'SIGI(+)':
                sigiplus = 'yes'

            if colList[count] == 'I(-)':
                iminus = 'yes'

            if colList[count] == 'SIGI(-)':
                sigiminus = 'yes'           

            if labelList[count] == 'J' and ilabel == 'none':
                ilabel = colList[count]

            if labelList[count] == 'Q' and sigilabel == 'none':
                sigilabel = colList[count]

            count = count + 1

        if ilabel == 'none' or sigilabel == 'none':
            print 'MTZ types for intensity data could not be established'
            time.sleep(4)
            return 1

        if imean == 'no' or sigimean == 'no' or iplus == 'no' or \
           sigiplus == 'no' or iminus == 'no' or sigiminus == 'no':
            print 'Standard MTZ labels IMEAN SIGIMEAN I(+) SIGI(+) I(-) SIGI(-) were not found'
            time.sleep(4)
            return 1

        if acell == 'none' or bcell == 'none' or ccell == 'none' \
           or alpha == 'none' or 'beta' == 'none' or gamma == 'none' or spacegroup_no == 'none':
            print 'Cell or spacegroup were not found in the reflection data file'
            time.sleep(4)
            return 1

    # Setup the cell parameter

    crystal_cell = acell + ' ' + bcell + ' ' + ccell + ' ' + alpha + ' ' + beta + ' ' + gamma


    # Correct the spacegroup number if it given by user

    if spacegroup_input != 'none':
        spacegroup_no = spacegroup_input

        print 'Changing spacegroup to ',spacegroup_no

    #########################################
    # Fix sort/asymmetric unit to standard  #
    #########################################

    print 'Running CCP4/CAD to set standard sort order'

    filename_inp = 'mi_cad.inp'
    filename_log = 'mi_cad.log'

    file = open(filename_inp,'w')
    file.write('LABIN FILE_NUMBER 1 ALL\n')
    file.write('SORT H K L \n')

    if spacegroup_input != 'none':
        file.write('SYMMETRY ')
        file.write(spacegroup_input)
        file.write('\n')

    file.write('END\n')
    file.close()

    runcad = 'cad HKLIN1 ' + filename_mtz + ' HKLOUT ' + filename_out_mtz + ' < ' + filename_inp + ' > ' + filename_log

    os.system(runcad)

    fileexists = os.path.exists(filename_out_mtz)
    if fileexists != 0:
        os.remove(filename_inp)
        os.remove(filename_log)
    else:
        print 'Output mtz file from CAD sort was not found'
        time.sleep(4)
        return 1

    ##########################################################
    # Run TRUNCATE to get F's from I's and CCP4 ANOM columns #
    ##########################################################

    print 'Running CCP4/TRUNCATE to reduce data from I to F'

    filename_inp = 'mi_truncate.inp'

    file = open(filename_inp,'w')
    file.write('title ')
    file.write(shelx_basename)
    file.write('\ntruncate yes\n')
    file.write('NRESIDUE 500\n')
    file.write('NOHARVEST\n')
    file.write('END\n')
    file.close()

    runtruncate = 'truncate hklin ' + filename_out_mtz + ' hklout ' + filename_out_mtz_f + ' < ' + filename_inp + ' > mi_truncate.log'

    os.system(runtruncate)

    fileexists = os.path.exists(filename_out_mtz_f)
    if fileexists != 0:
        os.remove(filename_inp)
    else:
        print 'Output mtz file from TRUNCATE was not found'
        time.sleep(4)
        return 1

    # Check the TRUNCATE log 

    fileexists = os.path.exists('mi_truncate.log')
    if fileexists != 0:
        file = open('mi_truncate.log','r')
        allLines = file.readlines()
        file.close()

        read_resolution = 'no'
        resolution_mtz = 'none'

        for eachLine in allLines:

            line_length = len(eachLine)

            if eachLine.find('Beware-serious ANISOTROPY') > -1:
                print '\nWarning - data is seriously anisotropic\n'
                time.sleep(4)

            if eachLine.find('* Space group =') > -1:
                aLine = eachLine.split('number')
                spacegroup_no = aLine[1]
                spacegroup_no = spacegroup_no.replace(')','')
                spacegroup_no = spacegroup_no.strip()

                aLine = eachLine.split(quote)
                spacegroup_lattice = aLine[1]
                spacegroup_lattice = spacegroup_lattice[0:1]
                spacegroup_lattice = spacegroup_lattice.strip()

            if read_resolution == 'yes' and line_length > 1:
                aLine = eachLine.split()
                resolution_mtz = aLine[5]
                read_resolution = 'no'

            if eachLine.find('*  Resolution Range :') > -1:
                read_resolution = 'yes'


    if spacegroup_lattice == 'none':
        print 'Space group lattice was not parsed'
        time.sleep(4)
        return 1

    if spacegroup_no == 'none':
        print 'Space group number was not parsed'
        time.sleep(4)
        return 1

    if resolution_mtz == 'none':
        print 'Resolution was not parsed'
        time.sleep(4)
        return 1

    print 'Resolution of input data is',resolution_mtz,'A'
    print 'Space group',spacegroup_no
    print 'Log file from CCP4/TRUNCATE:  mi_truncate.log'
    print 'Data file from CCP4/TRUNCATE:',filename_out_mtz_f

    #################################################
    # Find sites with SHELXD if not given as input  #
    #################################################

    if sitefile == 'none':

        print '\nFINDING SITES WITH SHELXD\n'

        auto_site_find = 'yes'

        ###################################  
        # Prepare Fanom data for SHELXD   #
        ###################################

        print 'Running CCP4/MTZ2VARIOUS to write anomalous differences\n'

        file_anom_shelx = shelx_basename + '.hkl'

        file = open('mi_mtz2various.inp','w')
        file.write('LABIN DP=DANO SIGDP=SIGDANO\n')
        file.write('OUTPUT SHELXDIff\n')
        file.write('END\n')
        file.close()

        runmtz = 'mtz2various HKLIN ' + filename_out_mtz_f + ' HKLOUT ' + file_anom_shelx + ' < mi_mtz2various.inp > mi_mtz2various.log'
        os.system(runmtz)

        fileexists = os.path.exists('mi_mtz2various.inp')
        if fileexists != 0:
            os.remove('mi_mtz2various.inp')

        fileexists = os.path.exists('mi_mtz2various.log')
        if fileexists != 0:
            os.remove('mi_mtz2various.log')

        # Write reflection list without header and obtain symmetry

        file = open(file_anom_shelx,'r')
        allLines = file.readlines()
        file.close()

        os.remove(file_anom_shelx)

        file = open(file_anom_shelx,'w')

        for eachLine in allLines:

            tag = eachLine[0:5]
            tag = tag.strip()

            if tag != 'TITLE' and tag != 'CELL' and tag != 'ZERR' and tag != 'LATT' and tag != 'SYMM' and tag != 'HKLF':
                file.write(eachLine)

            if tag == 'SYMM':
                symList.append(eachLine)

        file.close()

        num_sym = len(symList)

        #################################
        # Execute SHELXD to find sites  #
        #################################

        fileexists = os.path.exists('mi_runshelxd.log')
        if fileexists != 0:
            os.remove('mi_runshelxd.log')

        resolution_mtz = float(resolution_mtz)

        # Start search at reduced resolution

        if resolution_mtz < 2.4:
            resolution_search = resolution_mtz + 0.5
        else:
            resolution_search = resolution_mtz

        while resolution_search < 3.5:

            resolution_use = str(resolution_search)

            shelxd_ins = shelx_basename + '.ins'
            file = open(shelxd_ins,'w')

            file.write('TITL SAD phasing\n')
            file.write('CELL ')
            file.write(wavelength)
            file.write(' ')
            file.write(crystal_cell)
            file.write('\n')
            file.write('LATT -')

            if spacegroup_lattice == 'P':
                file.write('1')

            if spacegroup_lattice == 'I':
                file.write('2')

            if spacegroup_lattice == 'R' or spacegroup_lattice == 'H':
                file.write('3')

            if spacegroup_lattice == 'F':
                file.write('4')

            if spacegroup_lattice == 'C':
                file.write('7')

            file.write('\n')

            count = 0
            while count < num_sym:
                sym_op = symList[count]
                file.write(sym_op)

                count = count + 1

            file.write('SFAC S\n')
            file.write('UNIT 288\n')
            file.write('SHELL 15 ')
            file.write(resolution_use)
            file.write('\n')
            file.write('PATS\n')
            file.write('FIND ')
            file.write(number_sites)
            file.write('\n')

            if separation != 'none':
                file.write('MIND -')
                file.write(separation)
                file.write('\n')
            else:
                file.write('MIND -3.5\n')

            if number_ss != 'none':
                file.write('DSUL ')
                file.write(number_ss)
                file.write('\n')

            file.write('NTRY 50\n')
            file.write('SEED 1\n')
            file.write('HKLF 3\n')
            file.write('END\n')
            file.close()

            runshelxd = '"' + shelxd + '" ' + shelx_basename + ' > mi_runshelxd.log'
            os.system(runshelxd)

            # Check sites file output for CC value

            shelxd_res = shelx_basename + '.res'

            fileexists = os.path.exists(shelxd_res)
            if fileexists != 0:

                file = open(shelxd_res,'r')
                allLines = file.readlines()
                file.close()

                for eachLine in allLines:

                    if eachLine.find('CC(weak)') > -1:
                        aLine = eachLine.split()

                        number_items = len(aLine)
                        if number_items > 6:
                            cc = aLine[4]
                            cc_weak = aLine[6]

                os.remove(shelxd_res)

                # Rewrite results to file encoded with resolution

                shelxd_res_out = shelx_basename + resolution_use + '.res'

                file = open(shelxd_res_out,'w')
                for eachLine in allLines:
                    file.write(eachLine)
                file.close()

                aList_res_files.append(shelxd_res_out)
                aList_cc_weak.append(cc_weak)
                aList_cc.append(cc)
                aList_res.append(resolution_use)

                print 'Solution file:',shelxd_res_out
                print 'Best CC/CC(weak):',cc,'/',cc_weak,'at resolution',resolution_use,'A\n'

            else:

                print 'The SHELXD results file was not located:',shelxd_res
                time.sleep(4)
                return 1

            os.remove('mi_runshelxd.log')

            # decrease resolution for next test

            resolution_search = resolution_search + 0.2

        ##############################
        # Assess and Log SHELXD runs #
        ##############################

        # Collect the best cc_weak solution from the list

        number_files = len(aList_res_files)
        cc_weak_best = 0.0

        count = 0
        while count < number_files:
            cc_weak = aList_cc_weak[count]
            cc_weak = float(cc_weak)

            if cc_weak > cc_weak_best:
                use_sitefile = aList_res_files[count]
                cc_weak_best = cc_weak

            count = count + 1

        # Record to HTML file

        print 'HTML summary and index:',html_site_summary

        runtime = time.ctime(time.time())

        file = open(html_site_summary,'w')
        file.write('<html>\n')
        file.write('<head><title>SHELXD site finder summary</title></head>\n')
        file.write('<body bgcolor = "white">\n')
        file.write('<h1><center>SHELXD site finder summary</center></h1>\n')
        file.write('<p>\n')

        file.write('<table border = 1>\n')
        file.write('<tr><td><b>Job time</b></td><td>')
        file.write(runtime)
        file.write('</td></tr>\n')
        file.write('<tr><td><b>Working directory</b></td><td>')
        file.write('<a href = "file:///')
        file.write(workingdir)
        file.write('">')
        file.write(workingdir)
        file.write('</td></tr>\n')
        file.write('<tr><td><b>Data file</b></td><td>')
        file.write(filename_out_mtz_f)
        file.write('</td></tr>\n')
        file.write('</table>\n')

        file.write('<p>\n')
        file.write('<table border=1>\n')
        file.write('<tr>\n')
        file.write('<tr bgcolor = "yellow">\n')
        file.write('<td>File</td>')
        file.write('<td>Res. (&#197)</td>')
        file.write('<td>CC</td>')
        file.write('<td>CC(weak)</td>')
        file.write('</tr>\n')
        file.write('<tr>\n')

        # Insert table data

        count = 0
        while count < number_files:

            resfile = aList_res_files[count]
            cc = aList_cc[count]
            ccweak = aList_cc_weak[count]
            resolution = aList_res[count]

            resfile_full = os.path.join(workingdir,resfile)

            file.write('<td><a href = "file:///')
            file.write(resfile_full)
            file.write('">')
            file.write(resfile)
            file.write('</a></td>\n')

            file.write('<td>')
            file.write(resolution)
            file.write('</td>\n')

            file.write('<td>')
            file.write(cc)
            file.write('</td>\n')

            file.write('<td>')
            file.write(ccweak)
            file.write('</td>\n')

            count = count + 1

            file.write('</tr>\n')

        # Write HTML tail

        file.write('</table>')
        file.write('</body>\n')
        file.write('</html>\n')

        file.close()

    else:

        # Bring the imported site file local to avoid path issues on Windows

        print 'Using input scatter sites:',sitefile

        file = open(sitefile,'r')
        allLines = file.readlines()
        file.close()

        if sitefile.find('.pdb') > -1:
            use_sitefile = 'mi_sitefile.pdb'
        else:
            use_sitefile = 'mi_sitefile.res'

        file = open(use_sitefile,'w')

        for eachLine in allLines:
            file.write(eachLine)

        file.close()

    # Clean-up temporary file debris

    fileexists = os.path.exists(filename_local_in)
    if fileexists != 0:
        os.remove(filename_local_in)

    fileexists = os.path.exists(filename_mtz)
    if fileexists != 0:
        os.remove(filename_mtz)

    fileexists = os.path.exists(filename_out_mtz)
    if fileexists != 0:
        os.remove(filename_out_mtz)

    filename = shelx_basename + '.hkl'
    fileexists = os.path.exists(filename)
    if fileexists != 0:
        os.remove(filename)

    filename = shelx_basename + '.ins'
    fileexists = os.path.exists(filename)
    if fileexists != 0:
        os.remove(filename)

    filename = shelx_basename + '.lst'
    fileexists = os.path.exists(filename)
    if fileexists != 0:
        os.remove(filename)

    filename = shelx_basename + '.pdb'
    fileexists = os.path.exists(filename)
    if fileexists != 0:
        os.remove(filename)

    ##################################################################
    # Obtain fractions coordinate list from input site coordinates   #
    ##################################################################

    print '\nSITE REFINEMENT AND PHASING\n'

    file = open(use_sitefile)
    allLines = file.readlines()
    file.close()

    if use_sitefile.find('.pdb') == -1 and use_sitefile.find('.res') == -1:
        print 'Input site files should be in PDB (.pdb) or SHELX (.res) formats'
        time.sleep(4)
        return 1

    # Sites in PDB format

    if use_sitefile.find('.pdb') > -1:

        # Check integrity

        pdb_cryst1 = 'no'

        for eachLine in allLines:
            if eachLine.find('CRYST1') > -1:
                pdb_cryst1 = 'yes'

        if pdb_cryst1 == 'no':
            print 'The input PDB file must contain a CRYST1 record'
            time.sleep(4)
            return 1

        # Get element code

        for eachLine in allLines:
            tag = eachLine[0:6]
            tag = tag.strip()

            if tag == 'ATOM' or tag == 'HETATM':
                element = eachLine[76:78]
                element = element.strip()

                length_element = len(element)
                if length_element > 0:
                    scatterer = element
                else:
                    scatterer = scatterer_in

                atomList_elements.append(scatterer)

        # Convert PDB input format to fractional

        file = open('mi_coordconv.inp','w')
        file.write('CELL ')
        file.write(crystal_cell)
        file.write('\nINPUT PDB\n')
        file.write('OUTPUT FRAC\n')
        file.write('END\n')
        file.close()

        run_coordconv = 'coordconv xyzin ' + use_sitefile + ' xyzout mi_sites.frc < mi_coordconv.inp > mi_coordconv.log'
        os.system(run_coordconv)

        fileexists = os.path.exists('mi_sites.frc')
        if fileexists != 0:

            os.remove('mi_coordconv.inp')
            os.remove('mi_coordconv.log')

            file = open('mi_sites.frc')
            allLines = file.readlines()
            file.close()

            os.remove('mi_sites.frc')

            for eachLine in allLines:

                aLine = eachLine.split()
                x = aLine[1]
                y = aLine[2]
                z = aLine[3]

                atom_out = x + ' ' + y + ' ' + z
                atomList.append(atom_out)

                # Precompute inversion

                x = float(x)
                y = float(y)
                z = float(z)
                x = -x
                y = -y
                z = -z
                x = str(x)
                y = str(y)
                z = str(z)

                atom_out = x + ' ' + y + ' ' + z
                atomList_i.append(atom_out)

        else:

            print 'Conversion of PDB coordinates to fractional coordinates failed'
            time.sleep(4)
            return 1

    if use_sitefile.find('.res') > -1:

        # Otherwise assume SHELX res file format

        for eachLine in allLines:
            tag = eachLine[0:4]
            tag = tag.strip()

            if tag != 'REM' and tag != 'TITL' and tag != 'CELL' and tag !='LATT' and tag != 'SYMM' and \
               tag != 'SFAC' and tag != 'UNIT' and tag !='HKLF' and tag !='END':

                aLine = eachLine.split()
                line_length = len(aLine)

                if line_length == 7:

                    # Obtain element

                    if scatterer_in == 'S':

                        atom_name = aLine[0]
                        element1 = atom_name[0]
                        element2 = atom_name[1]

                        if element2 == '0':
                            scatterer = element1
                        else:
                            scatterer = element1 + element2

                    else:

                        scatterer = scatterer_in

                    atomList_elements.append(scatterer)

                    # Coordinates

                    x = aLine[2]
                    y = aLine[3]
                    z = aLine[4]

                    atom_out = x + ' ' + y + ' ' + z
                    atomList.append(atom_out)

                    # Inverted coordinates

                    x = float(x)
                    y = float(y)
                    z = float(z)
                    x = -x
                    y = -y
                    z = -z
                    x = str(x)
                    y = str(y)
                    z = str(z)

                    atom_out = x + ' ' + y + ' ' + z
                    atomList_i.append(atom_out)

    # Check for mismatch in site number inputs

    number_sites = int(number_sites)
    if number_ss != 'none' and use_sitefile != 'none':
        number_sites = number_sites + int(number_ss)

    number_fracsites = len(atomList)
    if number_fracsites == 0:
        print 'No sites are available in fractional coordinates'
        time.sleep(4)
        return 1

    if number_sites > number_fracsites:
        number_sites = number_fracsites

    # Delete temporary copy files when input site file was given

    fileexists = os.path.exists('mi_sitefile.pdb')
    if fileexists != 0:
        os.remove('mi_sitefile.pdb')

    fileexists = os.path.exists('mi_sitefile.res')
    if fileexists != 0:
        os.remove('mi_sitefile.res')    

    # Set output file names for site refinement

    phased_mtz = shelx_basename + 'phased.mtz'
    phased_mtz_i = shelx_basename + 'phased_i.mtz'

    fileexists = os.path.exists(phased_mtz)
    if fileexists != 0:
        os.remove(phased_mtz)

    fileexists = os.path.exists(phased_mtz_i)
    if fileexists != 0:
        os.remove(phased_mtz_i)

    ########################################
    # Refine the sites/phase with BP3      #
    ########################################

    if refine_method == 'bp3':

        print 'Running CCP4/BP3 to refine sites and phase on hand 1'

        phased_pdb = shelx_basename + 'phased'
        phased_pdb_i = shelx_basename + 'phased_i'
        phased_pdb_actual = phased_pdb + '1.pdb'
        phased_pdb_i_actual = phased_pdb_i + '1.pdb'
        phased_pdb_full = os.path.join(workingdir,phased_pdb_actual)
        phased_pdb_full_i = os.path.join(workingdir,phased_pdb_i_actual)

        fileexists = os.path.exists(phased_pdb)
        if fileexists != 0:
            os.remove(phased_pdb)

        fileexists = os.path.exists(phased_pdb_i)
        if fileexists != 0:
            os.remove(phased_pdb_i)

        phased_log = shelx_basename + 'bp3.log'
        phased_log_i = shelx_basename + 'bp3_i.log'

        file = open('mi_bp3.inp','w')
        file.write('Xtal DER1\n')

        # Refine the expected number of sites

        count = 0
        while count < number_sites:

            atom_out = atomList[count]
            scatter = atomList_elements[count]

            file.write('  ATOM ')
            file.write(scatterer)
            file.write('\n')

            file.write('    XYZ ')
            file.write(atom_out)
            file.write('\n')

            file.write('    OCCU 1.0\n')
            file.write('    BISO 25.0\n')

            count = count + 1

        # Column assignments

        file.write(' DNAME PEAK\n')
        file.write('   COLUmn F+=F(+) SF+=SIGF(+) F-=F(-) SF-=SIGF(-)\n')

        # Set f' and f'' for S and Cr edge (should be OK for Cu)

        if scatterer == 'S':
            file.write('   FORM ')
            file.write(scatterer)
            file.write(' FP=0.36 FPP=0.69\n')

        file.write('ALLIn\n')

        file.write('OUTPut ')
        file.write(phased_pdb)
        file.write('\n')

        file.close()

        # Execute

        run_bp3 = 'bp3 HKLIN ' + filename_out_mtz_f  + ' HKLOUT ' + phased_mtz + ' < mi_bp3.inp > ' + phased_log
        os.system(run_bp3)

        fileexists = os.path.exists(phased_mtz)
        if fileexists != 0:
            os.remove('mi_bp3.inp')

            print 'Log file from CCP4/BP3: ',phased_log
            print 'Sites from CCP4/BP3:    ',phased_pdb_actual
            print 'Data file from CCP4/BP3:',phased_mtz
            print 'MTZ labels for Fobs, FOM, Phase are F,FOM,PHIB'
        else:
            print 'BP3 site refinement failed'
            time.sleep(4)
            return 1


        # Run with inverted constellation

        if run_invert == 'yes':

            print '\nRunning CCP4/BP3 to refine sites and phase on hand 2'

            file = open('mi_bp3.inp','w')
            file.write('Xtal DER1\n')

            # Refine the expected number of sites

            count = 0
            while count < number_sites:

                atom_out = atomList_i[count]
                scatter = atomList_elements[count]

                file.write('  ATOM ')
                file.write(scatterer)
                file.write('\n')

                file.write('    XYZ ')
                file.write(atom_out)
                file.write('\n')

                file.write('    OCCU 1.0\n')
                file.write('    BISO 25.0\n')

                count = count + 1

            # Column assignments

            file.write(' DNAME PEAK\n')
            file.write('   COLUmn F+=F(+) SF+=SIGF(+) F-=F(-) SF-=SIGF(-)\n')

            # Set f' and f'' for S and Cr edge (should be OK for Cu)

            if scatterer == 'S':
                file.write('   FORM ')
                file.write(scatterer)
                file.write(' FP=0.36 FPP=0.69\n')

            file.write('ALLIn\n')

            file.write('OUTPut ')
            file.write(phased_pdb_i)
            file.write('\n')

            file.close()

            # Execute

            run_bp3 = 'bp3 HKLIN ' + filename_out_mtz_f  + ' HKLOUT ' + phased_mtz_i + ' < mi_bp3.inp > ' + phased_log_i
            os.system(run_bp3)

            fileexists = os.path.exists(phased_mtz_i)
            if fileexists != 0:
                os.remove('mi_bp3.inp')

                print 'Log file from CCP4/BP3: ',phased_log_i
                print 'Sites from CCP4/BP3:    ',phased_pdb_i_actual
                print 'Data file from CCP4/BP3:',phased_mtz_i        
            else:
                print 'BP3 inverted site refinement failed'
                time.sleep(4)
                return 1

        # Extract diagnostic information from BP3 log

        file = open(phased_log,'r')
        allLines = file.readlines()
        file.close()

        read_fom = 'no'
        read_atoms = 'no'
        read_atom_qb = 'no'

        for eachLine in allLines:

            # Capture FOM tables

            if eachLine.find('TOTAL') > -1 and read_fom == 'yes':
                read_fom = 'no'
                read_atoms = 'yes'

            if eachLine.find('Bin   HiRes  LoRes  1/Res^2   Refls  Centr    Refls  Acentr   Refls   All') > -1:
                read_fom = 'yes'

            if read_fom == 'yes':
                aLine = eachLine.split()
                num_parameters = len(aLine)

                if num_parameters == 10:
                    res = aLine[2]
                    fom = aLine[7]
                    number_refs = aLine[8]

                    aList_resolution_sites.append(res)
                    aList_numberrefs_sites.append(number_refs)
                    aList_fom_shell.append(fom)

            # Capture atom data

            if read_atoms == 'yes':

                if eachLine.find('Atom') > -1:
                    aLine = eachLine.split()
                    element = aLine[1]

                if read_atom_qb == 'yes':
                    aLine = eachLine.split()        
                    occ = aLine[0]
                    b = aLine[1]

                    aList_atom_out_element.append(element)
                    aList_atom_out_occ.append(occ)
                    aList_atom_out_b.append(b)

                    read_atom_qb = 'no'

                if eachLine.find('O        B') > -1:
                    read_atom_qb = 'yes'

        # Clean-up

        fileexists = os.path.exists('mi_sad_phased.xml')
        if fileexists != 0:
            os.remove('mi_sad_phased.xml')

        fileexists = os.path.exists('mi_sad_phased.sh')
        if fileexists != 0:
            os.remove('mi_sad_phased.sh')   

        fileexists = os.path.exists('mi_sad_phased_i.xml')
        if fileexists != 0:
            os.remove('mi_sad_phased_i.xml')

        fileexists = os.path.exists('mi_sad_phased_i.sh')
        if fileexists != 0:
            os.remove('mi_sad_phased_i.sh')

    ########################################
    # Refine the sites/phase with MLPHARE  #
    ########################################

    if refine_method == 'mlphare':

        print 'Running CCP4/MLPHARE to refine sites and phase on hand 1'

        phased_log = shelx_basename + 'mlphare.log'
        phased_log_i = shelx_basename + 'mlphare_i.log'
        phased_pdb = shelx_basename + 'phased1.pdb'
        phased_pdb_i = shelx_basename + 'phased_i1.pdb'
        phased_pdb_full = os.path.join(workingdir,phased_pdb)
        phased_pdb_full_i = os.path.join(workingdir,phased_pdb_i)
        phased_pdb_actual = phased_pdb
        phased_pdb_i_actual = phased_pdb_i

        fileexists = os.path.exists(phased_pdb)
        if fileexists != 0:
            os.remove(phased_pdb)

        fileexists = os.path.exists(phased_pdb_i)
        if fileexists != 0:
            os.remove(phased_pdb_i)

        file = open('mi_mlphare.inp','w')
        file.write('TITLE SAD site refinement\n')
        file.write('CYCLE 12\n')
        file.write('HLOUT\n')
        file.write('COORDS\n')
        file.write('LABIN FP=F SIGFP=SIGF FPH1=F SIGFPH1=SIGF DPH1=DANO SIGDPH1=SIGDANO\n')
        file.write('LABOUT FP=F SIGFP=SIGF\n')
        file.write('NOHARVEST\n')
        file.write(' DERIV SAD\n')
        file.write(' DCYCLE PHASE ALL REFCYC ALL KBOV ALL\n')

        # Refine the expected number of sites

        count = 0
        while count < number_sites:

            atom_out = atomList[count]
            scatter = atomList_elements[count]

            file.write(' ATOM  ')
            file.write(scatterer)
            file.write('   ')
            file.write(atom_out)
            file.write(' 0.0 1.0 BFAC 25.00\n')
            file.write(' ATREF AX ALL AY ALL AZ ALL AOCC ALL AB 9 10 11 12 \n')

            count = count + 1

        file.write('END\n')
        file.close()

        run_mlphare = 'mlphare HKLIN ' + filename_out_mtz_f  + ' HKLOUT ' + phased_mtz + \
                      ' XYZOUT ' + phased_pdb + ' < mi_mlphare.inp > ' + phased_log
        os.system(run_mlphare)

        fileexists = os.path.exists(phased_mtz)
        if fileexists != 0:
            os.remove('mi_mlphare.inp')

            print 'Log file from CCP4/MLPHARE: ',phased_log
            print 'Sites from CCP4/MLPHARE:    ',phased_pdb
            print 'Data file from CCP4/MLPHARE:',phased_mtz
            print 'MTZ labels for Fobs, FOM, Phase are F,FOM,PHIB'
        else:
            print 'MLPHARE site refinement failed'
            time.sleep(4)
            return 1

        # Fix MLPHARE pseudo PDB file to a usable PDB file

        fileexists = os.path.exists(phased_pdb)
        if fileexists != 0:

            file = open(phased_pdb,'r')
            allLines = file.readlines()
            file.close()

            os.remove(phased_pdb)

        else:

            print 'MLPHARE output site file not found'
            time.sleep(4)
            return 1

        file = open(phased_pdb,'w')

        atom_count = 0
        for eachLine in allLines:

            eachLine = eachLine.strip()

            tag = eachLine[0:6]
            tag = tag.strip()

            if tag == 'ATOM':

                start_Line = eachLine[0:17]
                end_Line = eachLine[26:80]

                atom_count = atom_count + 1
                res_number = str(atom_count)
                res_number = res_number.rjust(4)
                middle_Line = 'DUM A' + res_number

                eachLine = start_Line + middle_Line + end_Line

            file.write(eachLine)
            file.write('\n')

        file.close()

        # Run with inverted constellation

        if run_invert == 'yes':

            print '\nRunning CCP4/MLPHARE to refine sites and phase on hand 2'

            file = open('mi_mlphare.inp','w')
            file.write('TITLE SAD site refinement\n')
            file.write('CYCLE 12\n')
            file.write('HLOUT\n')
            file.write('COORDS\n')
            file.write('LABIN FP=F SIGFP=SIGF FPH1=F SIGFPH1=SIGF DPH1=DANO SIGDPH1=SIGDANO\n')
            file.write('LABOUT FP=F SIGFP=SIGF\n')
            file.write('NOHARVEST\n')
            file.write(' DERIV SAD\n')
            file.write(' DCYCLE PHASE ALL REFCYC ALL KBOV ALL\n')

            # Refine the expected number of sites

            count = 0
            while count < number_sites:

                atom_out = atomList_i[count]
                scatter = atomList_elements[count]

                file.write(' ATOM  ')
                file.write(scatterer)
                file.write('   ')
                file.write(atom_out)
                file.write(' 0.0 1.0 BFAC 25.00\n')
                file.write(' ATREF AX ALL AY ALL AZ ALL AOCC ALL AB 9 10 11 12 \n')

                count = count + 1

            file.write('END\n')
            file.close()

            run_mlphare = 'mlphare HKLIN ' + filename_out_mtz_f  + ' HKLOUT ' + phased_mtz_i + \
                          ' XYZOUT ' + phased_pdb_i + ' < mi_mlphare.inp > ' + phased_log_i
            os.system(run_mlphare)

            fileexists = os.path.exists(phased_mtz_i)
            if fileexists != 0:
                os.remove('mi_mlphare.inp')

                print 'Log file from CCP4/MLPHARE: ',phased_log_i
                print 'Sites from CCP4/MLPHARE:    ',phased_pdb_i
                print 'Data file from CCP4/MLPHARE:',phased_mtz_i        
            else:
                print 'MLPHARE inverted site refinement failed'
                time.sleep(4)
                return 1

            # Fix MLPHARE pseudo PDB file to a usable PDB file

            fileexists = os.path.exists(phased_pdb_i)
            if fileexists != 0:

                file = open(phased_pdb_i,'r')
                allLines = file.readlines()
                file.close()

                os.remove(phased_pdb_i)

            else:

                print 'MLPHARE inverted output site file not found'
                time.sleep(4)
                return 1

            file = open(phased_pdb_i,'w')

            atom_count = 0
            for eachLine in allLines:

                eachLine = eachLine.strip()

                tag = eachLine[0:6]
                tag = tag.strip()

                if tag == 'ATOM':

                    start_Line = eachLine[0:17]
                    end_Line = eachLine[26:80]

                    atom_count = atom_count + 1
                    res_number = str(atom_count)
                    res_number = res_number.rjust(4)
                    middle_Line = 'DUM A' + res_number

                    eachLine = start_Line + middle_Line + end_Line

                file.write(eachLine)
                file.write('\n')

            file.close()

        # Extract diagnostic information from MLPHARE log

        file = open(phased_log,'r')
        allLines = file.readlines()
        file.close()

        read_phase_data = 'no'
        read_res_shell = 'no'
        read_refs_shell = 'no'
        read_fom_shell = 'no'
        fom_done = 'no'

        for eachLine in allLines:

            # Capture FOM tables

            if read_phase_data == 'yes':

                if read_res_shell == 'yes':
                    aList_resolution_sites = eachLine.split()
                    read_res_shell = 'no'

                if read_refs_shell == 'yes':
                    aList_numberrefs_sites = eachLine.split()
                    read_refs_shell = 'no'

                if read_fom_shell == 'yes' and fom_done == 'no':
                    aList_fom_shell = eachLine.split()
                    read_fom_shell = 'no'
                    fom_done = 'yes'

                if eachLine.find('Resolution of each shell in angstroms:') > -1:
                    read_res_shell = 'yes'

                if eachLine.find('Number of Measurements phased -ACENTRIC') > -1:
                    read_refs_shell = 'yes'

                if eachLine.find('Mean Figure of Merit') > -1:
                    read_fom_shell = 'yes'

            if eachLine.find('This is the last cycle - PHASING ONLY') > -1:
                read_phase_data = 'yes'

            # Capture atom data

            if read_phase_data == 'yes' and eachLine.find('ATOM') > -1:

                aLine = eachLine.split()

                element = aLine[1]
                occ = aLine[6]
                b = aLine[8]

                aList_atom_out_element.append(element)
                aList_atom_out_occ.append(occ)
                aList_atom_out_b.append(b)

    #######################################
    # Report phasing results to terminal  #
    #######################################

    # Report FOM data

    print '\nFOM as a function of resolution shell\n'
    print ' Resolution   No. data    FOM'

    number_args = len(aList_resolution_sites)

    count = 0
    while count < number_args:

        res = aList_resolution_sites[count]
        number_refs = aList_numberrefs_sites[count]
        fom = aList_fom_shell[count]

        res = res.rjust(6)
        number_refs = number_refs.rjust(6)
        fom = fom.rjust(6)

        aLine = ' ' + res + '      ' + number_refs + '      ' + fom
        print aLine

        count = count + 1

    # Report atom data

    number_args = len(aList_atom_out_element)

    print '\nSite occupancies and B-factors\n'
    print ' Site Atom   Occ.     B'

    count = 0
    while count < number_args:

        site_number = count + 1
        site_number = str(site_number)
        element = aList_atom_out_element[count]
        occ = aList_atom_out_occ[count]
        b = aList_atom_out_b[count]

        site_number = site_number.rjust(4)
        element = element.rjust(3)
        occ = occ.rjust(7)
        b = b.rjust(8)

        aLine = site_number + ' ' + element + '  ' + occ + ' ' + b
        print aLine

        count = count + 1

    ######################################
    # Density modification with CCP4/DM  #
    ######################################

    # Setup input

    file = open('mi_dm.inp','w')
    file.write('SOLC ')
    file.write(solvent_fraction)
    file.write('\n')
    file.write('MODE SOLV HIST MULT\n')
    file.write('COMBINE PERT\n')
    file.write('NCYCLE AUTO\n')
    file.write('SCHEME ALL\n')
    file.write('LABIN FP=F SIGFP=SIGF PHIO=PHIB FOMO=FOM HLA=HLA HLB=HLB HLC=HLC HLD=HLD\n')
    file.write('LABOUT FDM=FDM PHIDM=PHIDM FOMDM=FOMDM\n')
    file.write('NOHARVEST\n')
    file.close()

    # Run for hand 1

    print '\nRunning CCP4/DM to refine phases on hand 1\n'

    protein_1 = '0.0'
    residual_1 = 'none'

    dm_mtz = shelx_basename + 'dm.mtz'
    dm_mtz_log = shelx_basename + 'dm.log'

    rundm = 'dm HKLIN ' + phased_mtz + ' HKLOUT ' + dm_mtz + ' < mi_dm.inp > ' + dm_mtz_log
    os.system(rundm)

    fileexists = os.path.exists(dm_mtz)
    if fileexists != 0:

        file = open(dm_mtz_log)
        allLines = file.readlines()
        file.close()

        for eachLine in allLines:
            if eachLine.find('Protein Resid') > -1:
                residual_1 = eachLine.strip()
    else:

        print 'DM calculation on hand 1 failed'
        time.sleep(4)
        return 1

    if residual_1 != 'none':

        aLine = residual_1.split()
        protein_1 = aLine[4]
        free_protein_1 = aLine[8]

        print 'Protein Residual:',protein_1

    print 'Refined mtz phase data from CCP4/DM from hand-1:',dm_mtz

    # Run for hand 2

    if run_invert == 'yes':

        print '\nRunning CCP4/DM to refine phases on hand 2\n'

        residual_2 = 'none'
        protein_2 = '0.0'

        dm_mtz_i = shelx_basename + 'dm_i.mtz'
        dm_mtz_log_i = shelx_basename + 'dm_i.log'

        rundm = 'dm HKLIN ' + phased_mtz_i + ' HKLOUT ' + dm_mtz_i + ' < mi_dm.inp > ' + dm_mtz_log_i
        os.system(rundm)

        fileexists = os.path.exists(dm_mtz_i)
        if fileexists != 0:

            file = open(dm_mtz_log_i)
            allLines = file.readlines()
            file.close()

            for eachLine in allLines:
                if eachLine.find('Protein Resid') > -1:
                    residual_2 = eachLine.strip()            
        else:

            print 'DM calculation on hand 2 failed'
            time.sleep(4)
            return 1

        if residual_2 != 'none':
            aLine = residual_2.split()
            protein_2 = aLine[4]
            free_protein_2 = aLine[8]

            print 'Protein Residual:',protein_2

        print 'Refined mtz phase data from CCP4/DM for hand-2:',dm_mtz_i

    os.remove('mi_dm.inp')

    # Select solution

    if run_invert == 'yes':
        test_protein_1 = float(protein_1)
        test_protein_2 = float(protein_2)

        if test_protein_1 < test_protein_2:
            print '\nHand 1 appears to be the correct solution'
        else:
            print '\nHand 2 appears to be the correct solution'

    #############################################
    # Create anomalous difference map-1 file    #
    #############################################

    anom_mtz = shelx_basename + 'anom_dm.mtz'
    anom_map = shelx_basename + 'anom.map'

    # Safeguard against generation of huge maps

    resolution_map = 'none'
    if resolution_mtz < 2.0:
        resolution_map = '2.0'

    file = open('mi_cad.inp','w')
    file.write('LABIN FILE 1 E1=F E2=SIGF E3=DANO E4=SIGDANO\n')
    file.write('LABIN FILE 2 E1=FOMDM E2=PHIDM\n')
    file.write('END\n')
    file.close()

    runcad = 'cad hklin1 ' + filename_out_mtz_f + ' hklin2 ' + dm_mtz + ' hklout ' + anom_mtz + ' < mi_cad.inp > mi_cad.log'
    os.system(runcad)

    fileexists = os.path.exists(anom_mtz)
    if fileexists == 0:
        print 'Output mtz file from CAD was not found'
        time.sleep(4)
        return 1

    # FFT to anomalous difference map-1. Note this automatically rotates phases.

    file = open('mi_fft.inp','w')
    file.write('LABIN F1=DANO W=FOMDM PHI=PHIDM\n')

    if resolution_map != 'none':
        file.write('RESOLUTION 100.0 ')
        file.write(resolution_map)
        file.write('\n')
    
    file.write('END\n')
    file.close()

    runfft = 'fft HKLIN ' + anom_mtz + ' MAPOUT ' + anom_map + ' < mi_fft.inp > mi_fft.log'
    os.system(runfft)

    fileexists = os.path.exists(anom_map)
    if fileexists == 0:
        print 'FFT for map display failed'
        time.sleep(4)
        return 1
    else:
        os.remove(anom_mtz)
        print 'Anomalous difference map from hand-1:',anom_map

    #############################################
    # Create anomalous difference map-2 file    #
    #############################################

    if run_invert == 'yes':

        anom_mtz_i = shelx_basename + 'anom_dm_i.mtz'
        anom_map_i = shelx_basename + 'anom_i.map'

        runcad = 'cad hklin1 ' + filename_out_mtz_f + ' hklin2 ' + dm_mtz_i + ' hklout ' + anom_mtz_i + ' < mi_cad.inp > mi_cad.log'
        os.system(runcad)

        fileexists = os.path.exists(anom_mtz_i)
        if fileexists == 0:
            print 'Output mtz file from CAD was not found'
            time.sleep(4)
            return 1

        # FFT to anomalous difference map-2 (note this automatically rotates phases)

        runfft = 'fft HKLIN ' + anom_mtz_i + ' MAPOUT ' + anom_map_i + ' < mi_fft.inp > mi_fft.log'
        os.system(runfft)

        fileexists = os.path.exists(anom_map_i)
        if fileexists == 0:
            print 'FFT for map display failed'
            time.sleep(4)
            return 1
        else:
            os.remove(anom_mtz_i)
            print 'Anomalous difference map for hand-2:',anom_map_i

    # Tidy

    print 'MTZ labels for F, [F(+)-F(-)], FOM, Phase are F,DANO,FOMDM,PHIDM'

    os.remove('mi_fft.inp')
    os.remove('mi_fft.log')
    os.remove('mi_cad.inp')
    os.remove('mi_cad.log')

    ###########################################
    # Option to write CCP4 maps for SS viewer #
    ###########################################

    if write_map == 'yes':

        # Map names

        dm_map = shelx_basename + 'dm.map'
        dm_map_i = shelx_basename + 'dm_i.map'

        # FFT

        file = open('mi_fft.inp','w')
        file.write('LABIN F1=F W=FOMDM PHI=PHIDM\n')
        file.write('END\n')
        file.close()

        runfft = 'fft HKLIN ' + anom_mtz + ' MAPOUT ' + dm_map + ' < mi_fft.inp > mi_fft.log'
        os.system(runfft)

        fileexists = os.path.exists(dm_map)
        if fileexists == 0:
            print 'FFT for map display failed'
            time.sleep(4)
            #return 1
        else:
            os.remove('mi_fft.inp')
            os.remove('mi_fft.log')

        # Option to write map with inverted atom constellation

        if run_invert == 'yes':

            file = open('mi_fft.inp','w')
            file.write('LABIN F1=F W=FOMDM PHI=PHIDM\n')
            file.write('END\n')
            file.close()

            runfft = 'fft HKLIN ' + anom_mtz_i + ' MAPOUT ' + dm_map_i + ' < mi_fft.inp > mi_fft.log'
            os.system(runfft)

            fileexists = os.path.exists(dm_map_i)
            if fileexists == 0:
                print 'FFT for map display failed'
                time.sleep(4)
                #return 1
            else:
                os.remove('mi_fft.inp')
                os.remove('mi_fft.log')

    ####################################################
    # Establish full names for session files and HTML  #
    ####################################################

    phased_log_full = os.path.join(workingdir,phased_log)
    dm_mtz_log_full = os.path.join(workingdir,dm_mtz_log)

    if run_invert == 'yes':
        phased_log_full_i = os.path.join(workingdir,phased_log_i)
        dm_mtz_log_full_i = os.path.join(workingdir,dm_mtz_log_i)

    sitefile_print = os.path.basename(use_sitefile)

    if auto_site_find == 'yes':
        sitefile_print_full = os.path.join(workingdir,sitefile_print)
    else:
        sitefile_print_full = sitefile

    html_site_summary_full = os.path.join(workingdir,html_site_summary)

    #############################################
    # Setup crystal and session files for MIFit #
    #############################################

    if do_mlw_files != 'none':

        mlw_file = 'mi_sad.mlw'       

        if run_invert == 'yes':
            print 'Session file contains maps for for hand-1 and hand-2:',mlw_file
        else:
            print 'Session file contains map for for hand-1:',mlw_file            

        file = open(mlw_file,'w')

        file.write('MapColumns FO=F FOM=FOMDM PHI=PHIDM\n')
        file.write('LoadMapPhase 1 ')
        file.write(dm_mtz)
        file.write('\n')
        file.write('silentmode\n')
        file.write('coefficients Fo*fom\n')
        file.write('fftapply\n')
        file.write('maptocont 1\n')
        file.write('maplinewidth 1.000000\n')
        file.write('contourlevels 21\n')
        file.write('contourleveldefault 50.000000 100.000000 150.000000 200.000000 250.000000\n')
        file.write('contourradius 14.000000\n')
        file.write('color 21\n')
        file.write('contourcolor 1\n')
        file.write('color 22\n')
        file.write('contourcolor 2\n')
        file.write('color 23\n')
        file.write('contourcolor 3\n')
        file.write('color 24\n')
        file.write('contourcolor 4\n')
        file.write('color 25\n')
        file.write('contourcolor 5\n')
        file.write('contourmap 1\n')

        # Load hand-2 map into same mlw file

        if run_invert == 'yes':

            file.write('MapColumns FO=F FOM=FOMDM PHI=PHIDM\n')
            file.write('LoadMapPhase 2 ')
            file.write(dm_mtz_i)
            file.write('\n')
            file.write('silentmode\n')
            file.write('coefficients Fo*fom\n')
            file.write('fftapply\n')
            file.write('maptocont 2\n')
            file.write('maplinewidth 1.000000\n')
            file.write('contourlevels 21\n')
            file.write('contourleveldefault 50.000000 100.000000 150.000000 200.000000 250.000000\n')
            file.write('contourradius 14.000000\n')
            file.write('color 26\n')
            file.write('contourcolor 1\n')
            file.write('color 27\n')
            file.write('contourcolor 2\n')
            file.write('color 28\n')
            file.write('contourcolor 3\n')
            file.write('color 29\n')
            file.write('contourcolor 4\n')
            file.write('color 30\n')
            file.write('contourcolor 5\n')
            file.write('contourmap 2\n')

        #
        
        file.write('translation 30.0 30.0 30.0\n')    
        file.write('rotation 1.0000 0.0000 0.0000 0.0000 1.0000 0.0000 0.0000 0.0000 1.0000\n')
        file.write('zoom 21\n')
        file.write('perspective 0.000\n')
        file.write('frontclip 4.0\n')
        file.write('backclip -4.0\n')
        file.write('transform\n')
        file.write('stereo off\n')
        file.close()


    ####################
    # Log phasing runs #
    ####################

    print '\nHTML summary and index:',html_phase_summary

    runtime = time.ctime(time.time())

    # Header

    file = open(html_phase_summary,'w')
    file.write('<html>\n')
    file.write('<head><title>Site refinement and phasing summary</title></head>\n')
    file.write('<body bgcolor = "white">\n')
    file.write('<h1><center>Site refinement and phasing summary</center></h1>\n')

    # I/O information

    file.write('<table border=1>\n')

    file.write('<tr><td>')
    file.write('<b>Job time</b></td><td>')
    file.write(runtime)
    file.write('</td></tr>\n')
    file.write('<tr><td>')
    file.write('<b>Working directory</b></td><td>')
    file.write('<a href = "file:///')
    file.write(workingdir)
    file.write('">')
    file.write(workingdir)
    file.write('</a></td></tr>\n')
    file.write('<tr><td>')
    file.write('<b>Input mtz data file</b></td><td>')
    file.write(filename_out_mtz_f)
    file.write('</a></td></tr>\n')
    file.write('<tr><td>')

    if sitefile != 'none':
        file.write('<b>Input site data file</b></td><td>')
        file.write('<a href = "file:///')
        file.write(sitefile_print_full)
        file.write('">')
        file.write(sitefile_print_full)
    else:
        file.write('<b>Automated site selection</b></td><td>')
        file.write('<a href = "file:///')
        file.write(html_site_summary_full)
        file.write('">')
        file.write(html_site_summary)

    file.write('</a></td></tr>\n')
    file.write('</table>\n')

    #

    file.write('<p><h2>Site and phase refinement data</h2>\n')
    file.write('<table border=1>')

    if do_mlw_files != 'none':
        file.write('<tr><td>Session file for MIFit</td>')
        file.write('<td>')
        file.write(mlw_file)
        file.write('</td></tr>\n')
    
    file.write('<tr><td>Hand 1 log file from site refinement</td>')
    file.write('<td>')
    file.write('<a href = "file:///')
    file.write(phased_log_full)
    file.write('">')
    file.write(phased_log)
    file.write('</a>')
    file.write('</td></tr>\n')
    file.write('<tr><td>Hand 1 site file</td>')
    file.write('<td>')
    file.write('<a href = "file:///')
    file.write(phased_pdb_full)
    file.write('">')
    file.write(phased_pdb_actual)
    file.write('</a>')
    file.write('</td></tr>')
    file.write('</td></tr>\n')
    file.write('<tr><td>Hand 1 anomalous difference map</td><td>')
    file.write(anom_map)    
    file.write('<tr><td>Hand 1 mtz data file from site refinement</td>')
    file.write('<td>')
    file.write(phased_mtz)
    file.write('</td></tr>\n')
    file.write('<tr><td>Hand 1 log file from CCP4/DM</td>')
    file.write('<td>')
    file.write('<a href = "file:///')
    file.write(dm_mtz_log_full)
    file.write('">')
    file.write(dm_mtz_log)
    file.write('</a>')
    file.write('</td></tr>\n')
    file.write('<tr><td>Hand 1 mtz data file from CCP4/DM</td>')
    file.write('<td>')
    file.write(dm_mtz)
    file.write('</td></tr>\n')

    if run_invert == 'yes':

        file.write('<tr><td>Hand 2 log file from site refinement</td>')
        file.write('<td>')
        file.write('<a href = "file:///')
        file.write(phased_log_full_i)
        file.write('">')
        file.write(phased_log_i)
        file.write('</a>')
        file.write('</td></tr>\n')
        file.write('<tr><td>Hand 2 sites from site refinement</td>')
        file.write('<td>')
        file.write('<a href = "file:///')
        file.write(phased_pdb_full_i)
        file.write('">')
        file.write(phased_pdb_i_actual)
        file.write('</a>')
        file.write('</td></tr>')
        file.write('</td></tr>\n')
        file.write('<tr><td>Hand 2 anomalous difference map</td><td>')
        file.write(anom_map_i)
        file.write('<tr><td>Hand 2 mtz data file from site refinement</td>')
        file.write('<td>')
        file.write(phased_mtz_i)
        file.write('</td></tr>\n')
        file.write('<tr><td>Hand 2 log file from CCP4/DM</td>')
        file.write('<td>')
        file.write('<a href = "file:///')
        file.write(dm_mtz_log_full_i)
        file.write('">')
        file.write(dm_mtz_log_i)
        file.write('</a>')
        file.write('</td></tr>\n')
        file.write('<tr><td>Hand 2 mtz data file from CCP4/DM</td>')
        file.write('<td>')
        file.write(dm_mtz_i)
        file.write('</td></tr>\n')

    file.write('<tr><td>MTZ label for F<sub>obs</sub></td><td>F</td></tr>')
    file.write('<tr><td>MTZ label for anomalous difference</td><td>DANO</td></tr>')
    file.write('<tr><td>MTZ label for initial FOM</td><td>FOM</td></tr>')
    file.write('<tr><td>MTZ label for initial phase</td><td>PHIB</td></tr>')
    file.write('<tr><td>MTZ label for modified FOM</td><td>FOMDM</td></tr>')
    file.write('<tr><td>MTZ label for modified phase</td><td>PHIDM</td></tr>')

    file.write('<tr><td>Hand 1 protein residual from CCP4/DM</td>')
    file.write('<td>')
    file.write(protein_1)
    file.write('</td></tr>')

    if run_invert == 'yes':
        file.write('<tr><td>Hand 2 protein residual from CCP4/DM</td>')
        file.write('<td>')
        file.write(protein_2)
        file.write('</td></tr>')

    file.write('</table>')

    # Insert FOM table data

    file.write('<h2>FOM as a function of resolution after site refinement</h2>')
    file.write('<table border=1>\n')
    file.write('<tr>\n')
    file.write('<tr bgcolor = "yellow">\n')
    file.write('<td>Res. (&#197)</td>')
    file.write('<td>No. phased data</td>')
    file.write('<td>FOM</td>')
    file.write('</tr>\n')

    number_args = len(aList_resolution_sites)

    count = 0
    while count < number_args:

        res = aList_resolution_sites[count]
        number_refs = aList_numberrefs_sites[count]
        fom = aList_fom_shell[count]

        file.write('<tr>')
        file.write('<td>')
        file.write(res)
        file.write('</td>')
        file.write('<td>')
        file.write(number_refs)
        file.write('</td>')
        file.write('<td>')
        file.write(fom)
        file.write('</td>')
        file.write('</tr>\n')

        count = count + 1

    file.write('</table>')
    file.write('<p>\n')

    # Insert atom parameter table

    number_args = len(aList_atom_out_element)

    file.write('<h2>Occupancies and B-factors after site refinement</h2>')
    file.write('<table border=1>\n')
    file.write('<tr>\n')
    file.write('<tr bgcolor = "yellow">\n')
    file.write('<td>Site</td>')
    file.write('<td>Atom</td>')
    file.write('<td>Occ</td>')
    file.write('<td>B</td>')
    file.write('</tr>\n')
    file.write('<tr>\n')

    count = 0
    while count < number_args:

        site_number = count + 1
        site_number = str(site_number)
        element = aList_atom_out_element[count]
        occ = aList_atom_out_occ[count]
        b = aList_atom_out_b[count]

        file.write('<tr>')
        file.write('<td>')
        file.write(site_number)
        file.write('</td>')
        file.write('<td>')
        file.write(element)
        file.write('</td>')
        file.write('<td>')
        file.write(occ)
        file.write('</td>')
        file.write('<td>')
        file.write(b)
        file.write('</td>')
        file.write('</tr>\n')

        count = count + 1

    file.write('</table>')

    # Write HTML tail

    file.write('</body>\n')
    file.write('</html>\n')

    file.close()

    # End
    print '\nNow go look at your maps !'
    time.sleep(4)

    # End 
    return 0

if __name__ == "__main__":
   sys.exit(Run())






