#####################################################################
#                                                                   #
# SAD phasing script using CCP4                                     # 
#                                                                   #
# Automated site determination within this script                   #
# possibly using SHELXD may be a future project                     #
#                                                                   #
# Copyright: Molecular Images   2006                                #
#                                                                   #
# This script is distributed under the same conditions as MIFit     #
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
    print "     --sitefile=FILE          The anomalous scattere site file. "
    print "     --sitenumber=NUM         Number of top anomalous scatterer sites read from sitefile. "
    print "     --scatterer=ELEM         The scattering element. Default: S"
    print "     --solventfraction=NUM    The solvent fraction.  No default"
    print "     --bothhands=yes or no    Flip site hand? Default=no."
    print "     --spacegroup_no=NUM      Specify a SG number. Default: deduce from data file"
    print "  -?,--help                   This help file"

def Run(argv=None):
    if argv is None:
        argv=sys.argv

    # Initialize

    inputfile = 'mi_runsadphase.txt'
    workingdir = 'none'
    datafile = 'none'
    sitefile = 'none'
    number_sites = 'none'
    solvent_fraction = 'none'
    scatterer_in = 'S'
    run_invert = 'no'
    
    #wavelength = '2.2909'
    wavelength = '0.0'
    acell = 'none'
    bcell = 'none'
    ccell = 'none'
    alpha = 'none'
    beta = 'none'
    gamma = 'none'
    spacegroup = 'none'
    spacegroup_input = 'none'
    spacegroup_no = 'none'
    filename_local_in = 'mi_sad.hkl'
    filename_mtz = 'mi_sad.mtz'
    filename_out_mtz = 'mi_sad_sort.mtz'

    refine_method = 'phaser'
    aList_x = []
    aList_y = []
    aList_z = []
    aLine = []

    ############################
    # Set input information    #
    ############################
    
    args = argv[1:]
    optlist, args = getopt.getopt(
        args,'d:h:s:?',
        ["workdir=","saddatafile=","sitefile=","sitenumber=",
         "scatterer=","solventfraction=",
         "bothhands=","spacegroup_no=","help"])
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
            elif arg_value=='--bothhands':
                run_invert = param_value
            elif arg_value=='--spacegroup_no':
                spacegroup_input = param_value

        count=count+1

    # Checks

    ccp4,error = ccp4check.ccp4check()
    if not ccp4:
      print '\n' + error + '\n'
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
    if fileexists == 0:
        print 'The anomalous scatterer site file was not found ',sitefile
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

    if number_sites == 'none' or number_sites == '0':
        print 'The number of expected anomalous scattering sites should be given'
        time.sleep(4)
        return 1

    if sitefile.find('.pdb') == -1 and sitefile.find('.res') == -1:
        print 'Input site files should be in either the PDB (.pdb) or SHELX (.res) format'
        time.sleep(4)
        return 1

    # Establish all filenames distinguishing for site inversion

    datafile_basename = os.path.basename(datafile)
    aLine = datafile_basename.split('.')
    datafile_basename = aLine[0]

    filename_out_mtz_f = datafile_basename + '_truncated.mtz'
    filename_out_mtz_f_log = datafile_basename + '_truncated.log'

    if refine_method == 'phaser':
        phasing_engine = '_phaser'
    else:
        phasing_engine = '_mlphare'

    if run_invert == 'yes':
        phasing_root = datafile_basename + '_inv' + phasing_engine
        dm_root = datafile_basename + '_inv_dm'
        runsummaryfile = '00run_summary_inv.txt'
    else:
        phasing_root = datafile_basename + phasing_engine
        dm_root = datafile_basename + '_dm'
        runsummaryfile = '00run_summary.txt'        
 
    phased_pdb = phasing_root + '.pdb'
    phased_mtz = phasing_root + '.mtz'
    phased_log = phasing_root + '.sol'
    dm_mtz = dm_root + '.mtz'
    dm_mtz_log = dm_root + '.log'

    # Remove previous versions of refined sites and density modification

    fileexists = os.path.exists(phased_pdb)
    if fileexists != 0:
        os.remove(phased_pdb)

    fileexists = os.path.exists(phased_log)
    if fileexists != 0:
        os.remove(phased_log)

    fileexists = os.path.exists(phased_mtz)
    if fileexists != 0:
        os.remove(phased_mtz)

    fileexists = os.path.exists(dm_mtz)
    if fileexists != 0:
        os.remove(dm_mtz)

    fileexists = os.path.exists(dm_mtz_log)
    if fileexists != 0:
        os.remove(dm_mtz_log)

    # Do all calculations in working directory

    os.chdir(workingdir)
    runtime = time.ctime(time.time())

    # Banner
    
    file = open(runsummaryfile,'w')
    aLine = '* SAD PHASING *'    
    file.write('\n')
    file.write(aLine)
    file.write('\n')
    print '\n',aLine
    aLine = 'Start time: ' + runtime
    file.write('\n')
    file.write(aLine)
    file.write('\n')
    print '\n',aLine   
    file.close()

    ##############################################
    # Convert to MTZ as standard working format  #
    ##############################################

    file = open(runsummaryfile,'a')
    aLine = 'PREPARING DATA'    
    file.write('\n')
    file.write(aLine)
    file.write('\n')
    print '\n',aLine
    file.close()

    #
    # 1. D*TREK case
    #

    if datafile.find('.ref') > -1:

        file = open(datafile,'r')
        allLines = file.readlines()
        file.close()

        file = open(filename_local_in,'w')
        file.writelines(allLines)
        file.close()        

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

        file = open(runsummaryfile,'a')
        aLine = 'Input file appears to be from D*TREK'
        file.write(aLine)
        file.write('\n')
        print '\n',aLine
        file.close()

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

    if datafile.find('.sca') > -1:

        file = open(datafile,'r')
        allLines = file.readlines()
        file.close()

        file = open(filename_local_in,'w')
        file.writelines(allLines)
        file.close()      
        
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

        file = open(runsummaryfile,'a')
        aLine = 'Input file appears to be from SCALEPACK'
        file.write(aLine)
        file.write('\n')
        print '\n',aLine
        file.close()

        filename_inp = 'mi_runscalepack2mtz.inp'
        filename_log = 'mi_scalepack2mtz.log'

        file = open(filename_inp,'w')
        file.write('name proj ')
        file.write(datafile_basename)
        file.write(' crystal 1 dataset 1\n')
        file.write('wave ')
        file.write(wavelength)
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
    # 3. SCALA case
    #

    if datafile.find('.mtz') > -1:
        datatype = 'scala'

        file = open(runsummaryfile,'a')
        aLine = 'Input file appears to be from SCALA (MTZ)'
        file.write(aLine)
        file.write('\n')
        print '\n',aLine
        file.close()

        file = open(datafile,'rb')
        allLines = file.readlines()
        file.close()

        file = open(filename_local_in,'wb')
        file.writelines(allLines)
        file.close()

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

    # Correct the spacegroup number if overridden by user

    if spacegroup_input != 'none':
        spacegroup_no = spacegroup_input

        file = open(runsummaryfile,'a')
        aLine =  'Changing spacegroup to ' + str(spacegroup_no)        
        file.write(aLine)
        file.write('\n')
        print '\n',aLine
        file.close()

    #########################################
    # Fix sort/asymmetric unit to standard  #
    #########################################

    file = open(runsummaryfile,'a')
    aLine = 'Running CCP4/CAD to set standard sort order'    
    file.write(aLine)
    file.write('\n')
    print '\n',aLine
    file.close()

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

    file = open(runsummaryfile,'a')
    aLine = 'Running CCP4/TRUNCATE to reduce data from I to F'
    file.write(aLine)
    file.write('\n')
    print '\n',aLine
    file.close()

    filename_inp = 'mi_truncate.inp'

    file = open(filename_inp,'w')
    file.write('title ')
    file.write(datafile_basename)
    file.write('\ntruncate yes\n')
    file.write('NRESIDUE 500\n')
    file.write('NOHARVEST\n')
    file.write('END\n')
    file.close()

    runtruncate = 'truncate hklin ' + filename_out_mtz + ' hklout ' + filename_out_mtz_f + ' < ' + filename_inp + ' > ' \
                  + filename_out_mtz_f_log

    os.system(runtruncate)

    fileexists = os.path.exists(filename_out_mtz_f)
    if fileexists != 0:
        os.remove(filename_inp)
        
        file = open(runsummaryfile,'a')
        aLine =  'Log file from CCP4/TRUNCATE: ' + filename_out_mtz_f_log
        file.write(aLine)
        file.write('\n')
        print '\n',aLine
        aLine = 'Data file from CCP4/TRUNCATE:' + filename_out_mtz_f
        file.write(aLine)
        file.write('\n')
        print '\n',aLine
        file.close()
        
    else:
        print 'Output mtz file from TRUNCATE was not found'
        time.sleep(4)
        return 1

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

    ###################################
    # Obtain input site coordinates   #
    ###################################

    file = open(runsummaryfile,'a')
    aLine = 'SITE REFINEMENT'
    file.write('\n')
    file.write(aLine)
    file.write('\n')
    print '\n',aLine
    aLine = 'Number of sites set by user: ' + str(number_sites)
    file.write(aLine)
    file.write('\n')
    print '\n',aLine
    file.close()

    file = open(sitefile,'r')
    allLines = file.readlines()
    file.close()

    # Parse sites in PDB format

    if sitefile.find('.pdb') > -1:

        xyz_type = 'ORTH'

        for eachLine in allLines:
            tag = eachLine[0:6]
            tag = tag.strip()

            if tag == 'ATOM' or tag == 'HETATM':
                x = eachLine[30:38]
                y = eachLine[38:46]
                z = eachLine[46:54]
                aList_x.append(x)
                aList_y.append(y)
                aList_z.append(z)

    # or parse sites in SHELXD RES format

    elif sitefile.find('.res') > -1:

        xyz_type = 'FRAC'

        for eachLine in allLines:
            tag = eachLine[0:4]
            tag = tag.strip()

            if tag != 'REM' and tag != 'TITL' and tag != 'CELL' and tag !='LATT' and tag != 'SYMM' and \
               tag != 'SFAC' and tag != 'UNIT' and tag !='HKLF' and tag !='END':

                aLine = eachLine.split()
                line_length = len(aLine)

                if line_length == 7:
                    x = aLine[2]
                    y = aLine[3]
                    z = aLine[4]
                    aList_x.append(x)
                    aList_y.append(y)
                    aList_z.append(z)
        
    # Check that number of sites entered by user does not exceed number in file

    number_sites = int(number_sites)
    number_sites_file = len(aList_x)
    
    if number_sites_file < number_sites or number_sites == 0:
        number_sites = number_sites_file

    ############################
    # Refine sites with PHASER # 
    ############################

    if refine_method == 'phaser':

        file = open(runsummaryfile,'a')
        aLine = 'Running CCP4/PHASER to refine sites'
        file.write(aLine)
        file.write('\n')
        print '\n',aLine
        file.close()

        file = open('mi_phaser.inp','w')
        file.write('MODE EP_AUTO\n')
        file.write('TITLe sad phasing input file written by MIFit\n')
        file.write('HKLIn ')
        file.write(filename_out_mtz_f)
        file.write('\n')
        file.write('COMPosition SOLVent ')
        file.write(solvent_fraction)
        file.write('\n')
        file.write('CRYStal myprotein DATAset sad LABIn F+=F(+) SIG+=SIGF(+) F-=F(-) SIG-=SIGF(-)\n')
        file.write('LLGComplete CRYStal myprotein COMPLETE ON SCATtering ELEMent ')
        file.write(scatterer_in)
        file.write('\n')

        # Write coordinates
        
        count = 0
        while count < number_sites:

            x = aList_x[count]
            y = aList_y[count]
            z = aList_z[count]
            aLine = 'ATOM CRYstal myprotein ELEMent ' + scatterer_in + ' ' + xyz_type + ' ' + x + ' ' + y + ' ' + z + ' OCCUpancy 1.0'
            file.write(aLine)
            file.write('\n')

            count = count + 1
        
        file.write('ROOT ')
        file.write(phasing_root)
        file.write('\n')

        # Option to invert hand of site constellation
        
        if run_invert == 'yes':
            file.write('HAND ON\n')
        
        file.close()

        run_phaser = 'phaser <  mi_phaser.inp > mi_phaser.log'
        os.system(run_phaser)        

        fileexists = os.path.exists(phased_mtz)
        if fileexists != 0:

            file = open(runsummaryfile,'a')
            aLine = 'Log file from CCP4/PHASER: ' + str(phased_log)
            file.write(aLine)
            file.write('\n')
            print '\n',aLine
            aLine = 'Sites from CCP4/PHASER:    ' + str(phased_pdb)
            file.write(aLine)
            file.write('\n')
            print '\n',aLine
            aLine = 'Data file from CCP4/PHASER:' + str(phased_mtz)
            file.write(aLine)
            file.write('\n')
            aLine = 'MTZ labels for Fobs, FOM, Phase are F,FOM,PHIB'
            file.write(aLine)
            file.write('\n')
            print '\n',aLine
            file.close()
        
        else:
            print 'PHASER site refinement failed'
            time.sleep(4)
            return 1

    #########################################################
    # Refine the sites/phase with MLPHARE (legacy code)     #
    # Note coordinates must be fractional in this case      #
    #########################################################

    if refine_method == 'mlphare':

        file = open(runsummaryfile,'a') 
        aLine =  'Running CCP4/MLPHARE to refine sites'
        file.write(aLine)
        file.write('\n')
        print '\n',aLine
        file.close()

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

        # Refine the given number of sites (must be fractional)

        count = 0
        while count < number_sites:

            x = aList_x[count]
            y = aList_y[count]
            z = aList_z[count]

            if run_invert == 'yes':
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
            file.write(' ATOM  ')
            file.write(scatterer_in)
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
            
            file = open(runsummaryfile,'a')
            aLine = 'Log file from CCP4/MLPHARE: ' + str(phased_log)
            file.write(aLine)
            file.write('\n')
            print '\n',aLine
            aLine = 'Sites from CCP4/MLPHARE:    ' + str(phased_pdb)
            file.write(aLine)
            file.write('\n')
            print '\n',aLine
            aLine = 'Data file from CCP4/MLPHARE:' + str(phased_mtz)
            file.write(aLine)
            file.write('\n')            
            aLine = 'MTZ labels for Fobs, FOM, Phase are F,FOM,PHIB'
            file.write(aLine)
            file.write('\n')
            print '\n',aLine
            file.close()
        
        else:
            print 'MLPHARE site refinement failed'
            time.sleep(4)
            return 1

    ######################################
    # Density modification with CCP4/DM  #
    ######################################

    file = open(runsummaryfile,'a')
    aLine = 'PHASE REFINEMENT'
    file.write('\n')
    file.write(aLine)
    file.write('\n')
    print '\n',aLine
    aLine = 'Running CCP4/DM to refine phases'
    file.write(aLine)
    file.write('\n')
    print '\n',aLine
    file.close()

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

    rundm = 'dm HKLIN ' + phased_mtz + ' HKLOUT ' + dm_mtz + ' < mi_dm.inp > ' + dm_mtz_log
    os.system(rundm)

    fileexists = os.path.exists(dm_mtz)
    if fileexists != 0:

        file = open(runsummaryfile,'a')
        aLine =  'Refined mtz phase data from CCP4/DM:' + str(dm_mtz)
        file.write(aLine)
        file.write('\n')
        print '\n',aLine
        aLine =  'MTZ labels for maps using F, FOM, Phase are FDM,FOMDM,PHIDM'
        file.write(aLine)
        file.write('\n')
        print '\n',aLine
        file.close()
        
    else:
        print 'DM calculation failed'
        time.sleep(4)
        return 1

    ########
    # End  #
    ########

    runtime = time.ctime(time.time())

    file = open(runsummaryfile,'a')
    aLine = 'Now go look at your map !'
    file.write('\n')
    file.write(aLine)
    file.write('\n')
    print '\n',aLine    
    aLine = 'Note: It may be necessary to try an inverted site constellation'
    file.write('\n')
    file.write(aLine)
    file.write('\n')
    print '\n',aLine
    aLine = 'End time: ' + runtime
    file.write('\n')
    file.write(aLine)
    file.write('\n')
    print '\n',aLine 
    file.close()
        
    return 0

if __name__ == "__main__":
   sys.exit(Run())






