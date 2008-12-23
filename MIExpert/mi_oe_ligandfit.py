#####################################################################
# Script: mi_oe_ligandfit.py                                        #
# Release: Consortium                                               #
#                                                                   #
# Manages OpenEye process for ligand density docking                #
#                                                                   #
# Copyright: Molecular Images   2006                                #
#                                                                   #
# This script is distributed under the same conditions as MIFit     #
#                                                                   #
#####################################################################

import sys
import os
import getopt
import time
import string
import math
import ccp4check

def Usage():
    print "Usage: %s [options]" % sys.argv[0]
    print "Options are:"
    print "  -p,--pdbfile=FILE          pdb file to use"
    print "  -m,--mtzfile=FILE          mtz file to use"
    print "  -l,--ligandfile=FILE       ligand file to use"
    print "  -x,--xcenter=FLOAT         X coordinate of center of ligand region"
    print "  -y,--ycenter=FLOAT         Y coordinate of center of ligand region"
    print "  -z,--zcenter=FLOAT         Z coordinate of center of ligand region"
    print "  -s,--oe_path=DIR           Path to OpenEye executable"
    print "  -n,--max_confs=NUM         Maximum number of confs to write per blob. Default: 1, max: 5"
    print "  -b,--input_border=FLOAT    half-length of size of region to search. Default: 10, max: 15"
    print "  -a,--f_map=LABEL           label of F1 coefficient. Default: deduce from mtz file"
    print "  -f,--phase_map=LABEL       label of PHI coefficient. Default: PHIC"
    print "  -?,--help                  This help"

def Run(argv=None):
    if argv is None:
        argv=sys.argv

    # Initialize

    quote = """'"""

    oe_program = 'none'

    workingdir = 'none'
    pdbfile = 'none'
    mtzfile = 'none'
    ligandfile = 'none'
    x_center = 'none'
    y_center = 'none'
    z_center = 'none'

    # May hardwire a path here
    #oe_path = 'C:/OpenEye/AFITT'
    oe_path = 'none'

    input_max_confs = 'none'
    input_border = 'none'

    fmap = 'none'
    phase_map = 'PHIC'

    a_cell = 'none'
    b_cell = 'none'
    c_cell = 'none'
    alpha_cell = 'none'
    beta_cell = 'none'
    gamma_cell = 'none'
    a11 = 'none'
    a12 = 'none'
    a13 = 'none'
    a21 = 'none'
    a22 = 'none'
    a23 = 'none'
    a31 = 'none'
    a32 = 'none'
    a33 = 'none'
    spacegroup_name = 'none'
    crystal_pdb = 'none'
    blob_number = '1'
    prev_blob_number = '1'

    grid_x = 'none'
    grid_y = 'none'
    grid_z = 'none'
    axis_order = 'none'

    aList_chain = ['Y','H','I','J','K','L','M','N','O','P','Q','R','S','T','U','V','W','X','Z','a','b','c','d','e','f','g']

    # These define the volume of map to search (10A half-side) and number of solutions to display.
    # Note that default is to produce 5 solutions per blob

    border = 10.0
    max_confs = 1

    # Parse args
    args = argv[1:]
    optlist, args = getopt.getopt(
        args,'m:p:l:h:x:y:z:s:n:b:a:f:?',
        ['help','mtzfile=','pdbfile=','ligandfile=',
         'xcenter=','ycenter=','zcenter=','oe_path=','max_confs=','input_border=',
         'f_map=','phase_map='])
    number_of_inputs = len(optlist)
    if number_of_inputs==0:
        Usage()
        return
    count = 0
    while count < number_of_inputs:
        aList = optlist[count]
        number_of_list_inputs = len(aList)
        if number_of_list_inputs >=1:
            arg_value = aList[0]
            if arg_value== '-?' or arg_value=='--help':
                Usage()
                return
        if number_of_list_inputs >=2:
            param_value = aList[1]
            if arg_value == '-m' or arg_value=='--mtzfile':
                mtzfile = param_value
            elif arg_value == '-p' or arg_value=='--pdbfile':
                pdbfile = param_value
            elif arg_value == '-l' or arg_value=='--ligandfile':
                ligandfile = param_value
            elif arg_value == '-x' or arg_value=='--xcenter':
                x_center = param_value
            elif arg_value == '-y' or arg_value=='--ycenter':
                y_center = param_value
            elif arg_value == '-z' or arg_value=='--zcenter':
                z_center = param_value
            elif (arg_value == '-s' or arg_value=='--oe_path') and oe_path=='none':
                oe_path = param_value
            elif arg_value == '-n' or arg_value=='--max_confs':
                input_max_confs = param_value
            elif arg_value == '-b' or arg_value=='--input_border':
                input_border = param_value
            elif arg_value == '-a' or arg_value=='--f_map':
                fmap = param_value
            elif arg_value == '-f' or arg_value=='--phase_map':
                phase_map = param_value
        else:
            Usage()
            return

        count = count + 1

    ccp4,error = ccp4check.ccp4check()
    if not ccp4:
      print '\n' + error + '\n'
      time.sleep(4)
      return 1

    # Verify input 

    fileexists = os.path.exists(mtzfile)
    if fileexists == 0:
        print 'The input mtz file was not found: ',mtzfile
        time.sleep(4)
        return 1

    fileexists = os.path.exists(pdbfile)
    if fileexists == 0:    
        print 'The input pdb file was not found: ',pdbfile
        time.sleep(4)
        return 1
    else:
        pdbfile_oe = pdbfile

    fileexists = os.path.exists(ligandfile)
    if fileexists == 0:
        print 'The input ligand file was not found: ',ligandfile
        time.sleep(4)
        return 1

    fileexists = os.path.exists(oe_path)
    if fileexists == 0:
        print 'The OpenEye directory was not found: ',oe_path
        time.sleep(4)
        return 1

    if x_center == 'none':
        print 'The center x position was not found'
        time.sleep(4)
        return 1    

    if y_center == 'none':
        print 'The center y position was not found'
        time.sleep(4)
        return 1

    if z_center == 'none':
        print 'The center z position was not found'
        time.sleep(4)
        return 1

    if input_max_confs != 'none':
        max_confs = int(input_max_confs)

    if input_border != 'none':
        border = float(input_border)

    if max_confs > 5:
        max_confs = 5

    if border > 15.0:
        border = 15.0

    # Go to working area if full path was given

    if os.path.basename(pdbfile) != pdbfile:
        workingdir = os.path.dirname(pdbfile)
        os.chdir(workingdir)

    # Copy MTZ file local

    file = open(mtzfile,'rb')
    allLines = file.readlines()
    file.close()

    file = open('mi_oe.mtz','wb')
    file.writelines(allLines)
    file.close()

    ########################################################
    # For pseudo-interactive runs deduce what kind of sf   #
    # coefficients we have in mtz file written by MIFit.   #
    # Note that mtz coefficients from MIFit are hardcoded  #
    # FP SIGFP FC FMAP PHIC FOM FreeR_Flag                 #
    ########################################################

    rfactor = -1.0
    use_2fp_fc = 'no'

    if fmap == 'none':

        # Run CCP4/SCALEIT to check FP and FMAP

        fileexists = os.path.exists('mi_oe_scaled.mtz')
        if fileexists != 0:
            os.remove('mi_oe_scaled.mtz')

        file = open('mi_scaleit.inp','w')
        file.write('LABIN FP=FP SIGFP=SIGFP FPH1=FMAP SIGFPH1=SIGFP\n')
        file.write('END\n')
        file.close()

        runscaleit = 'scaleit HKLIN mi_oe.mtz HKLOUT mi_oe_scaled.mtz < mi_scaleit.inp > mi_scaleit.log'
        os.system(runscaleit)

        # Parse SCALEIT log

        fileexists = os.path.exists('mi_oe_scaled.mtz')
        if fileexists != 0:

            file = open('mi_scaleit.log','r')
            allLines = file.readlines()
            file.close()

            os.remove('mi_scaleit.inp')
            os.remove('mi_scaleit.log')
            os.remove('mi_oe_scaled.mtz')

            for eachLine in allLines:

                if eachLine.find('THE TOTALS') > -1:
                    aLine = eachLine.split()

                    num_params = len(aLine)

                    if num_params > 7:
                        rfactor = aLine[7]
                        rfactor = float(rfactor)
        else:

            print 'SCALEIT run to check data coefficients failed'
            return 1
            time.sleep(4)

        if rfactor < 0.1 and rfactor > -1.0:

            # FP=FMAP so the user entered precomputed map data (did not use any MIFit computation)                

            print '\nUsing input FP as precomputed map amplitude'
            fmap = 'FP'

        else:

            # If FP != FMAP disagree so get scaling for a 2Fo-Fc map calculation

            print '\nUsing input FP and FC without further scaling to form 2Fo-Fc map amplitude'
            use_2fp_fc = 'yes'

    ############################################
    # Create map around ligand target density  #
    ############################################

    fileexists = os.path.exists('mi.map')
    if fileexists != 0:
        os.remove('mi.map')

    file = open('mi_fft.inp','w')

    # Map is from user-defined 'FP' coefficient or as supplied FP,FC

    if use_2fp_fc  == 'no':
        file.write('LABIN F1=')
        file.write(fmap)
        file.write(' PHI=')
        file.write(phase_map)
    else:
        file.write('LABIN F1=FP F2=FC PHI=PHIC\n')
        file.write('SCALE F1 2.0 0.0 F2 1.0 0.0\n')

    file.write('\nEND\n')
    file.close()

    runfft = 'fft HKLIN mi_oe.mtz MAPOUT mi.map < mi_fft.inp > mi_fft.log'
    os.system(runfft)

    fileexists = os.path.exists('mi.map')
    if fileexists == 0:
        print 'FFT for ligand_fitting failed'
        time.sleep(4)
        return 1
    else:
        os.remove('mi_fft.inp')
        os.remove('mi_oe.mtz')

    # Parse FFT data

    file = open('mi_fft.log','r')
    allLines = file.readlines()
    file.close()

    os.remove('mi_fft.log')

    read_cell = 'no'

    for eachLine in allLines:

        if eachLine.find('* Space group') > -1:
            aLine = eachLine.split(quote)
            number_args = len(aLine)

            if number_args == 3:
                spacegroup_name = aLine[1]

        if read_cell == 'yes':
            aLine = eachLine.split()
            number_args = len(aLine)

            if number_args == 6:
                a_cell = aLine[0]
                b_cell = aLine[1]
                c_cell = aLine[2]
                alpha_cell = aLine[3]
                beta_cell = aLine[4]
                gamma_cell = aLine[5]

                read_cell = 'no'

        if eachLine.find('* Cell Dimensions') > -1:
            read_cell = 'yes'

        if eachLine.find('Grid sampling on x, y, z') > -1:
            cutLine = eachLine[61:100]
            aLine = cutLine.split()
            number_args = len(aLine)

            if number_args == 3:
                grid_x = aLine[0]
                grid_y = aLine[1]
                grid_z = aLine[2]

        if eachLine.find('Fast, medium, slow axes') > -1:
            axis_order = eachLine[61:100]

    # Verify

    if a_cell == 'none' or b_cell == 'none' or c_cell == 'none' or \
       alpha_cell == 'none' or beta_cell == 'none' or gamma_cell =='none':
        print 'Cell dimensions were not found'
        time.sleep(4)
        return 1
    else:
        cell = a_cell + ' ' + b_cell + ' ' + c_cell + ' ' + alpha_cell + ' ' + beta_cell + ' ' + gamma_cell

    if spacegroup_name == 'none':
        print 'Space group number was not found'
        time.sleep(4)
        return 1

    if grid_x == 'none' or grid_y == 'none' or grid_z == 'none':
        print 'FFT grid spacings were not found'
        time.sleep(4)
        return 1
    else:
        grid = grid_x + ' ' + grid_y + ' ' + grid_z

    if axis_order == 'none':
        print 'Axis order was not found'
        time.sleep(4)
        return 1

    print 'Map calculation done'

    ############################################################################
    # Get deorthogonalization matrix by forcing CCP4/PDBSET to recalculate it  #
    ############################################################################

    file = open(pdbfile,'r')
    allLines = file.readlines()
    file.close()

    file = open('mi_temp.pdb','w')

    for eachLine in allLines:

        tag = eachLine[0:6]
        tag = tag.strip()

        if tag == 'ATOM' or tag == 'HETATM' or tag == 'CRYST1':
            file.write(eachLine)

    file.close()

    file = open('mi_pdbset.inp','w')
    file.write('CELL ')
    file.write(cell)
    file.write('\n')
    file.write('SPACEGROUP ')
    file.write(quote)
    file.write(spacegroup_name)
    file.write(quote)
    file.write('\n')
    file.close()

    runpdbset = 'pdbset XYZIN mi_temp.pdb XYZOUT mi_temp_out.pdb < mi_pdbset.inp > mi_pdbset.log'
    os.system(runpdbset)

    fileexists = os.path.exists('mi_temp_out.pdb')
    if fileexists == 0:
        print 'PDBSET to establish crystal parameters failed'
        time.sleep(4)
        return 1
    else:
        os.remove('mi_pdbset.inp')
        os.remove('mi_pdbset.log')

    # Parse out CRYST1 record and SCALE records

    file = open('mi_temp_out.pdb','r')
    allLines = file.readlines()
    file.close()

    os.remove('mi_temp_out.pdb')

    for eachLine in allLines:

        tag = eachLine[0:6]
        tag = tag.strip()

        eachLine = eachLine.strip()

        if tag == 'CRYST1':
            crystal_pdb = eachLine.strip()

        if tag == 'SCALE1':
            aList = eachLine.split()
            a11 = aList[1]
            a12 = aList[2]
            a13 = aList[3]
            a11 = float(a11)
            a12 = float(a12)
            a13 = float(a13)

        if tag == 'SCALE2':
            aList = eachLine.split()
            a21 = aList[1]
            a22 = aList[2]
            a23 = aList[3]
            a21 = float(a21)
            a22 = float(a22)
            a23 = float(a23)

        if tag == 'SCALE3':
            aList = eachLine.split()
            a31 = aList[1]
            a32 = aList[2]
            a33 = aList[3]
            a31 = float(a31)
            a32 = float(a32)
            a33 = float(a33)

    if crystal_pdb == 'none':
        print 'CRYST1 record from PDBSET output was not found'
        time.sleep(4)
        return 1    

    if a11 == 'none' or a22 == 'none' or a33 == 'none':
        print 'SCALE records from PDBSET output were not found'
        time.sleep(4)
        return 1

    ##########################################################
    # Check if input PDB file has a CRYST1 record and insert #
    ##########################################################

    file = open(pdbfile,'r')
    allLines_protein = file.readlines()
    file.close()

    found_cryst1 = 'no'

    for eachLine in allLines:
        if eachLine.find('CRYST1') > -1:
            found_cryst1 = 'yes'

    if found_cryst1 == 'no':

        print 'Adding CRYST1 record to coordinates for OpenEye software\n'
        print crystal_pdb

        os.remove('mi_temp.pdb')

        file = open('mi_temp.pdb','w')
        file.write(crystal_pdb)
        file.write('\n')

        for eachLine in allLines:

            tag = eachLine[0:6]
            tag = tag.strip()

            eachLine = eachLine.strip()

            if tag == 'ATOM' or tag == 'HETATM' or tag == 'TER':

                file.write(eachLine)
                file.write('\n')

        file.write('END\n')
        file.close()

        pdbfile_oe = 'mi_temp.pdb'

    ####################################################################
    # Mask map to zero protein regions (OE prototype not doing it yet) # 
    ####################################################################

    # Run SFALL to create atom model for mask (1.5A radius)

    file = open('mi_sfall.inp','w')
    file.write('MODE ATMMAP\n')
    file.write('GRID ')
    file.write(grid)
    file.write('\n')
    file.write('SYMMETRY ')
    file.write(quote)
    file.write(spacegroup_name)
    file.write(quote)
    file.write('\n')
    file.write('SFSG 1\n')
    file.write('VDWR 1.5\n')
    file.write('END\n')
    file.close()

    runsfall = 'sfall XYZIN mi_temp.pdb MAPOUT mi_atoms.map < mi_sfall.inp > mi_sfall.log'
    os.system(runsfall)

    fileexists = os.path.exists('mi_atoms.map')
    if fileexists == 0:
        print 'CCP4/SFALL for masking failed'
        time.sleep(4)
        return 1
    else:
        os.remove('mi_sfall.inp')
        os.remove('mi_sfall.log')

    # Run MAPMASK to convert to mask with 1 protein, 0 solvent and fix axis order

    file = open('mi_mapmask.inp','w')
    file.write('AXIS ')
    file.write(axis_order)
    file.write('\n')
    file.write('MASK CUT 0.000001\n')
    file.close()

    runmapmask = 'mapmask MAPIN mi_atoms.map MSKOUT mi_atoms1.mask < mi_mapmask.inp > mi_mapmask.log'
    os.system(runmapmask)

    fileexists = os.path.exists('mi_atoms1.mask')
    if fileexists == 0:
        print 'CCP4/MAPMASK for map to mask conversion failed'
        time.sleep(4)
        return 1
    else:
        os.remove('mi_mapmask.inp')
        os.remove('mi_mapmask.log')
        os.remove('mi_atoms.map')

    # Run MAPMASK to invert mask to 0 protein, 1 solvent

    file = open('mi_mapmask.inp','w')
    file.write('SCALE FACTOR -1.0 1.0')
    file.close()

    runmapmask = 'mapmask MSKIN mi_atoms1.mask MSKOUT mi_atoms2.mask < mi_mapmask.inp > mi_mapmask.log'
    os.system(runmapmask)

    fileexists = os.path.exists('mi_atoms2.mask')
    if fileexists == 0:
        print 'CCP4/MAPMASK to invert mask failed'
        time.sleep(4)
        return 1
    else:
        os.remove('mi_mapmask.inp')
        os.remove('mi_mapmask.log')
        os.remove('mi_atoms1.mask')

    # Run MAPMASK to zero protein density

    file = open('mi_mapmask.inp','w')
    file.write('MAPS MULT')
    file.close()

    runmapmask = 'mapmask MAPIN1 mi.map MSKIN2 mi_atoms2.mask MAPOUT mi_masked.map < mi_mapmask.inp > mi_mapmask.log'
    os.system(runmapmask)

    fileexists = os.path.exists('mi_masked.map')
    if fileexists == 0:
        print 'CCP4/MAPMASK for protein masking failed'
        time.sleep(4)
        return 1
    else:
        os.remove('mi_mapmask.inp')
        os.remove('mi_mapmask.log')
        os.remove('mi_atoms2.mask')

        os.remove('mi.map')
        os.rename('mi_masked.map','mi.map')

    ################################################
    # Use CCP4/MAPMASK to box out required region  #
    ################################################

    # Obtain box coordinates

    x_center = float(x_center)
    y_center = float(y_center)
    z_center = float(z_center)

    xmin_x = x_center - border
    xmin_y = y_center 
    xmin_z = z_center

    ymin_x = x_center
    ymin_y = y_center - border
    ymin_z = z_center

    zmin_x = x_center
    zmin_y = y_center
    zmin_z = z_center - border

    xmax_x = x_center + border
    xmax_y = y_center
    xmax_z = z_center

    ymax_x = x_center
    ymax_y = y_center + border
    ymax_z = z_center

    zmax_x = x_center
    zmax_y = y_center
    zmax_z = z_center + border

    xmin = a11*xmin_x + a12*xmin_y + a13*xmin_z
    ymin = a21*ymin_x + a22*ymin_y + a23*ymin_z
    zmin = a31*zmin_x + a32*zmin_y + a33*zmin_z
    xmax = a11*xmax_x + a12*xmax_y + a13*xmax_z
    ymax = a21*ymax_x + a22*ymax_y + a23*ymax_z
    zmax = a31*zmax_x + a32*zmax_y + a33*zmax_z

    xmin = round(xmin,4)
    ymin = round(ymin,4)
    zmin = round(zmin,4)
    xmax = round(xmax,4)
    ymax = round(ymax,4)
    zmax = round(zmax,4)
    xmin = str(xmin)
    ymin = str(ymin)
    zmin = str(zmin)
    xmax = str(xmax)
    ymax = str(ymax)
    zmax = str(zmax)

    xyz_limits = xmin + ' ' + xmax + ' ' + ymin + ' ' + ymax + ' ' + zmin + ' ' + zmax

    fileexists = os.path.exists('mi_oe_map.ccp4')
    if fileexists != 0:
        os.remove('mi_oe_map.ccp4')

    file = open('mi_mapmask.inp','w')
    file.write('XYZLIM ')
    file.write(xyz_limits)
    file.write('\n')
    file.write('EXTEND XTAL\n')
    file.write('END\n')
    file.close()

    runmapmask = 'mapmask MAPIN mi.map MAPOUT mi_oe_map.ccp4 < mi_mapmask.inp > mi_mapmask.log'
    os.system(runmapmask)

    fileexists = os.path.exists('mi_oe_map.ccp4')
    if fileexists == 0:
        print 'CCP4/MAPMASK for target site failed'
        time.sleep(4)
        return 1
    else:
        f = os.popen('mapdump mapin mi_oe_map.ccp4')
        allLines = f.readlines()
        f.close()
        print 'Information for input map:'
        print '    Path:', os.path.abspath('mi_oe_map.ccp4')
        for eachLine in allLines:
            if eachLine.find('.....') > -1:
                print '    ' + eachLine.strip()
        print
        os.remove('mi_mapmask.inp')
        os.remove('mi_mapmask.log')
        os.remove('mi.map')

    print 'Calculations to reduced and mask map done'

    ####################
    # Run OE software  #
    ####################

    test_platform = sys.platform
    if test_platform.find('win') > -1:
        flynn = os.path.join(oe_path,'flynn.exe')
        afitt = os.path.join(oe_path,'afitt.bat')
    else:
        flynn = os.path.join(oe_path,'flynn')
        afitt = os.path.join(oe_path,'afitt')    

    fraglib = os.path.join(oe_path,'fraglib.oeb.gz')

    # Run

    runtime = time.ctime(time.time())

    #############################################
    # Standalone (flynn) or full afitt option   #
    #############################################

    print '\nLaunching OpenEye ligand fitting process'

    # Generate a unique file name for OpenEye output for repeated runs

    idcode = '000000'

    gmt = time.gmtime(time.time())
    fmt = '%H%M%S'
    idcode = time.strftime(fmt,gmt)

    oe_ligandfits = 'mi_oe_ligandfit_' + idcode + '.pdb'

    # Run Process

    fileexists = os.path.exists(flynn)
    if fileexists != 0:    
        runflynn = '"' + flynn + '" -fraglib ' + '"' + fraglib + '" -grid "' + os.path.abspath('mi_oe_map.ccp4') + '" ' + '-prot "' \
                   + pdbfile_oe + '" -verbose -in "' + ligandfile + '" -out "' + os.path.abspath(oe_ligandfits) + '"'

        oe_program = 'flynn'

    fileexists = os.path.exists(afitt)
    if fileexists != 0:
        runflynn = '"' + afitt + '" -e mifit_automate -grid "' + os.path.abspath('mi_oe_map.ccp4') + '" -ligand ' + '"' + ligandfile + '" -protein "' + pdbfile_oe + '" -exit'

        oe_program = 'afitt'

    if oe_program == 'none':
        print 'Neither flynn nor afitt were found'
        time.sleep(4)
        return 1

    if oe_program == 'flynn':
        print 'Using OpenEye/flynn'

    if oe_program == 'afitt':
        print 'Using OpenEye/afitt'

    print 'Will write',max_confs,'conformers for each density blob'
    print 'OE Command line',runflynn
    print 'Ligand results wil be stored in file',oe_ligandfits
    print 'Start time:',runtime

    file = open('mi_flynn.bat','w')
    file.write(runflynn)
    file.write('\n')
    file.close()

    test_platform = sys.platform
    if test_platform.find('win') == -1:
        os.system('chmod +x mi_flynn.bat')

    os.system('mi_flynn.bat')

    # Check run.
    # Note that the afitt automation.py script needs to be configured to write results.pdb by
    # changing the line that contains results%s.sdf to results%s.pdb 

    if oe_program == 'flynn':

        fileexists = os.path.exists(oe_ligandfits)
        if fileexists != 0:
            os.remove('mi_flynn.bat')
        else:
            print 'No output coordinate file from OE ligand docking run'
            time.sleep(4)
            return 1

    if oe_program == 'afitt':

        fileexists = os.path.exists('results.pdb')
        if fileexists != 0:
            os.remove('mi_flynn.bat')

            file = open('results.pdb','r')
            allLines = file.readlines()
            file.close()   

            os.remove('results.pdb')

            file = open(oe_ligandfits,'w')
            file.writelines(allLines)
            file.close()

        else:
            print 'No output coordinate file - results.pdb - from OE ligand docking run'
            time.sleep(4)
            return 1

    runtime = time.ctime(time.time())
    print 'End time:',runtime

    ###################################################
    # Splice fitted ligand conformers back into file  #
    ###################################################

    # Read protein pdb

    file = open(pdbfile,'r')
    allLines_protein = file.readlines()
    file.close()

    # Read ligand fits

    file = open(oe_ligandfits,'r')
    allLines_ligand = file.readlines()
    file.close()

    ligand_write = 'no'
    conformer_count = 0
    conformer_count_chain = 0
    chain_id = 'Y'

    # Write back into standard filename for MIFit

    file = open('mi_oe_ligandfit.pdb','w')

    for eachLine_protein in allLines_protein:

        tag = eachLine_protein[0:6]
        tag = tag.strip()

        # Re write selected PDB file items

        if tag.find('HEADER') > -1 or tag.find('CRYST1') > -1 or tag.find('SCALE') > -1 or tag.find('TER') > -1 \
            or tag.find('HETATM') > -1 or tag.find('ATOM') > -1:

            # Put ligand 

            if eachLine_protein.find('HOH') > -1 and ligand_write != 'done':

                for eachLine_ligand in allLines_ligand:

                    tag_ligand = eachLine_ligand[0:6]
                    tag_ligand = tag_ligand.strip()

                    if tag_ligand == 'CMPND':
                        aLine = tag_ligand.split(':')
                        blob_number = aLine[2]
                        blob_number = blob_number.strip()

                        if blob_number != prev_blob_number:
                            conformer_count = 0
                            prev_blob_number = blob_number

                    if tag_ligand == 'TER':
                        conformer_count = conformer_count + 1
                        conformer_count_chain = conformer_count_chain + 1
                        chain_id = aList_chain[conformer_count_chain]

                    if tag_ligand == 'ATOM' or tag_ligand == 'HETATM':

                        element = eachLine_ligand[76:80]
                        element = element.strip()

                        if element != 'H':

                            eachLine_ligand = eachLine_ligand.strip()

                            if conformer_count < max_confs:
                                eachLine_ligand_out = eachLine_ligand[0:21] + chain_id + eachLine_ligand[22:80]
                                file.write(eachLine_ligand_out)
                                file.write('\n')

                ligand_write = 'done'

            file.write(eachLine_protein)

    # Write ligand for no-water case

    if ligand_write != 'done':

        for eachLine_ligand in allLines_ligand:

            tag_ligand = eachLine_ligand[0:6]
            tag_ligand = tag_ligand.strip()

            if tag_ligand == 'CMPND':
                aLine = tag_ligand.split(':')
                blob_number = aLine[2]
                blob_number = blob_number.strip()

                if blob_number != prev_blob_number:
                    conformer_count = 0
                    prev_blob_number = blob_number

            if tag_ligand == 'TER':
                conformer_count = conformer_count + 1
                conformer_count_chain = conformer_count_chain + 1
                chain_id = aList_chain[conformer_count_chain]

            if tag_ligand == 'ATOM' or tag_ligand == 'HETATM':

                element = eachLine_ligand[76:80]
                element = element.strip()

                if element != 'H':

                    eachLine_ligand = eachLine_ligand.strip()

                    if conformer_count < max_confs:
                        eachLine_ligand_out = eachLine_ligand[0:21] + chain_id + eachLine_ligand[22:80]
                        file.write(eachLine_ligand_out)
                        file.write('\n')

    # Close out

    file.write('END\n')
    file.close()

    print '\nLigand and protein coordinates are combined in file: mi_oe_ligandfit.pdb' 

    #############
    # Clean-up  #
    #############

    fileexists = os.path.exists('mi_temp.pdb')
    if fileexists != 0:
        os.remove('mi_temp.pdb')

    fileexists = os.path.exists('mi_oe_map.ccp4')
    if fileexists != 0:
        os.remove('mi_oe_map.ccp4')            

    return 0

if __name__ == "__main__":
    sys.exit(Run())
