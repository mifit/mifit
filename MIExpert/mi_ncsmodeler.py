#####################################################################
#                                                                   #
# Remodel NCS copies of protein to a template model                 #
# and average map                                                   #
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
    print "  -p,--pdbfile=FILE          the pdb target file"
    print "  -c,--targetchain=ID        the target chain id"
    print "  -m,--mtzfile=FILE          the mtz file to use. Default: no file"
    print "  -k,--preserve_model=FILE   file containing fixed model segments. Default: no file"
    print "  -n,--maskextras=FILE       file defining extra entities for NCS masking. Default: no file"
    print "  -h,--phase_prob=yes or no  Use experimental phase probability. Default: no"
    print "  -s,--phase_comb=yes or no  Do phase combination.  Default: no"
    print "  -?,--help                  This help"
    print ""
    print "These options are only used when phase_prob and phase_comb are both no:"
    print "  -l,--fom_label=LABEL       the figure of merit column label."
    print "  -f,--phase_label=LABEL     the phase column label."


def Run(argv=None):
    if argv is None:
        argv=sys.argv

    # Initialize
    workingdir = '.'
    pdbtargetfile = 'none'
    targetchain = 'none'
    preservemodel_file = 'none'
    mtzfile = 'none'
    maskextras_file = 'none'

    crystal_record = 'none'
    pdb_spacegroup_name = 'P1'
    solvent_percent = '0.5'

    aList_chain = ['?']
    aList_ca_count = [0]
    aList_chain_start = ['?']
    aList_chain_end = ['?']
    aList_rmsd = []
    aList_output_atoms = []
    aList_input_atoms = []
    aList_deviations = []

    aList_fix_chain = []
    aList_fix_res_start = []
    aList_fix_res_end = []

    aList_ncs_chain = []
    aList_ncs_res_number = []

    aList_preserve_atom = []
    aList_preserve_chain = []
    aList_preserve_res_number = []
    aList_preserve_atom_name = []

    aList_rotation_matrix = []
    aList_translation_matrix = []

    flabel = 'none'
    sigflabel = 'none'
    rfreelabel = 'none'
    fomlabel = 'none'
    phaselabel = 'none'
    fomlabel_exp = 'none'
    phaselabel_exp = 'none'
    phase_prob = 'no'
    phase_comb = 'no'
    calc_sfs = 'yes'

    aList_hl = []

    aList_protein = ['ALA','ARG','ASN','ASP','CYS','GLN','GLU','GLY','HIS','ILE','LEU',\
                     'LYS','MET','PHE','PRO','SER','THR','TRP','TYR','VAL','MSE','PTR','SEP','TPO']

    number_protein_list = len(aList_protein)

    # Get a time string to almost uniquely encode the temp file names

    gmt = time.gmtime(time.time())
    fmt = '%H%M%S'
    idcode = time.strftime(fmt,gmt)

    # Temporary files

    mi_compar_pdb = 'mi_temp_compar_' + idcode + '.pdb'
    mi_refmac_pdb_use = 'mi_temp_use_' + idcode + '.pdb'
    mi_refmac_mtz_out = 'mi_temp_ref_' + idcode + '.mtz'
    mi_refmac_pdb_out = 'mi_temp_ref_' + idcode + '.pdb'
    mi_ncsmask_msk = 'mi_temp_ncsmask_' + idcode + '.msk'
    mi_ncsmask_pdb = 'mi_temp_ncsmask_' + idcode + '.pdb'
    mi_dm_mtz_use = 'mi_temp_dm_' + idcode + '.mtz'
    mi_dm_mtz_out = 'mi_temp_dm_out_' + idcode + '.mtz'

    # Parse args
    args = argv[1:]
    optlist, args = getopt.getopt(
        args,'p:c:k:m:f:l:n:h:s:?',
        ['pdbfile=','targetchain=','mtzfile=','preserve_model=',
         'maskextras=','phase_prob=','phase_comb=','help',
         'fom_label=','phase_label='])
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
            if arg_value == '-?' or arg_value=='--help':
                Usage()
                return
        if number_of_list_inputs >=2:
            param_value = aList[1]
            if arg_value == '-p' or arg_value=='--pdbfile':
                pdbtargetfile = param_value
            elif arg_value == '-c' or arg_value=='--targetchain':
                targetchain = param_value
            elif arg_value == '-k' or arg_value=='--preserve_model':
                preservemodel_file = param_value
            elif arg_value == '-m' or arg_value=='--mtzfile':
                mtzfile = param_value
            elif arg_value == '-l' or arg_value=='--fom_label':
                fomlabel_exp = param_value
            elif arg_value == '-f' or arg_value=='--phase_label':
                phaselabel_exp = param_value
            elif arg_value == '-n' or arg_value=='--maskextras':
                maskextras_file = param_value
            elif arg_value == '-h' or arg_value=='--phase_prob':
                phase_prob = param_value
            elif arg_value == '-s' or arg_value=='--phase_comb':
                phase_comb = param_value
        count = count + 1

    ccp4,error = ccp4check.ccp4check()
    if not ccp4:
      print '\n' + error + '\n'
      time.sleep(4)
      return 1
    
    # Setup full CCP4 program paths for ncsmask because an Xfit program has same name

    test_platform = sys.platform
    if test_platform.find('win') > -1:
        ncsmask_exe = 'ncsmask.exe'
    else:
        ncsmask_exe = 'ncsmask'

    ncsmask = os.path.join(ccp4.bin,ncsmask_exe)

    # Verify
    fileexists = os.path.exists(pdbtargetfile)
    if fileexists == 0:
        print 'The target PDB file was not found ',pdbtargetfile
        time.sleep(4)
        return 1

    if targetchain == 'none':
        print 'The target chain id must be given'
        time.sleep(4)
        return 1

    fileexists = os.path.exists(preservemodel_file)
    if fileexists == 0 and preservemodel_file != 'none':
        print 'The file containing fixed model segments was not found ',preservemodel_file
        time.sleep(4)
        return 1

    fileexists = os.path.exists(mtzfile)
    if fileexists == 0 and mtzfile != 'none':
        print 'The file containing mtz data was not found ',mtzfile
        time.sleep(4)
        return 1

    # Experimental phase inputs

    if phase_prob == 'yes' and phase_comb == 'yes':
        print 'Experimental phasing and phase combination are incompatible options'
        time.sleep(4)
        return 1

    if phase_prob == 'yes':
        calc_sfs = 'no'

    if phase_comb == 'yes':
        fomlabel = 'FOM'
        phaselabel = 'PHCOMB'

    if phase_prob == 'no' and phase_comb == 'no':
        if fomlabel_exp != 'none' and phaselabel_exp != 'none':
            fomlabel = fomlabel_exp
            phaselabel = phaselabel_exp
            calc_sfs = 'no'

    # Set output file names

    filename_ncs = os.path.basename(pdbtargetfile)
    filename_ncs = 'ncs_' + filename_ncs
    filename_ncs_mtz = filename_ncs.replace('.pdb','.mtz')
    filename_ncs_dmlog = filename_ncs.replace('.pdb','_dm.log')

    # Go to working area if full path was given

    if os.path.basename(pdbtargetfile) != pdbtargetfile:
        workingdir = os.path.dirname(pdbtargetfile)
        os.chdir(workingdir)

    #####################################################
    # Parse the file that defines fixed model segments  #
    #####################################################

    if preservemodel_file != 'none':

        file = open(preservemodel_file,'r')
        allLines = file.readlines()
        file.close()

        for eachLine in allLines:
            aLine = eachLine.split()

            number_args = len(aLine)
            if number_args == 3:

                chain = aLine[0]
                res_start = aLine[1]
                res_end = aLine[2]
                res_start = int(res_start)
                res_end = int(res_end)

                aList_fix_chain.append(chain)
                aList_fix_res_start.append(res_start)
                aList_fix_res_end.append(res_end)

    number_preserve_records = len(aList_fix_chain)

    ##############################################################
    # Parse the file that defines extra entities for NCS masking #
    ##############################################################

    if maskextras_file != 'none':

        file = open(maskextras_file,'r')
        allLines = file.readlines()
        file.close()

        for eachLine in allLines:
            aLine = eachLine.split()

            number_args = len(aLine)
            if number_args == 2:

                chain = aLine[0]
                res_number = aLine[1]

                aList_ncs_chain.append(chain)
                aList_ncs_res_number.append(res_number)

    number_ncs_extras = len(aList_ncs_chain)

    ##############################
    # Obtain chain list data     #
    ##############################

    file = open(pdbtargetfile,'r')
    allLines = file.readlines()
    file.close()

    for eachLine in allLines:

        eachLine = eachLine.strip()

        tag = eachLine[0:6]
        tag = tag.strip()

        if tag == 'CRYST1':
            crystal_record = eachLine.strip()
            pdb_spacegroup_name = crystal_record[55:66]
            pdb_spacegroup_name = pdb_spacegroup_name.replace(' ','')

        if tag == 'ATOM' or tag == 'HETATM':

            aList_input_atoms.append(eachLine)

            # Compile list of chains

            chain_id = eachLine[21:22]

            new_chain = 'yes'
            number_chains = len(aList_chain)

            count = 0
            while count < number_chains:
                chain_id_old = aList_chain[count]

                if chain_id == chain_id_old:
                    new_chain = 'no'

                count = count + 1

            if new_chain == 'yes':
                aList_chain.append(chain_id)
                aList_ca_count.append(0)
                aList_chain_end.append('?')
                aList_chain_start.append('?')

    # Obtain CA count for each chain

    for eachLine in allLines:

        tag = eachLine[0:6]
        tag = tag.strip()

        if tag == 'ATOM' or tag == 'HETATM':

            chain_id = eachLine[21:22]
            res_number = eachLine[22:26]
            res_number = res_number.strip()
            atom_name = eachLine[13:16]
            atom_name = atom_name.strip()

            count = 0
            while count < number_chains:
                chain_id_old = aList_chain[count]

                if chain_id == chain_id_old:

                    if atom_name == 'CA':
                        atom_count = aList_ca_count[count]
                        atom_count = atom_count + 1
                        aList_ca_count[count] = atom_count

                        aList_chain_end[count] = res_number

                    if aList_chain_start[count] == '?':
                        aList_chain_start[count] = res_number

                count = count + 1

    # Record chains contents for user (first entry is dummy)

    print '\n CHAIN DIAGNOSTICS\n'
    print ' Chain-id  No. CA atoms  Start-no   End-no'

    number_chains = len(aList_chain)
    count = 1
    while count < number_chains:

        chain_id = aList_chain[count]
        atom_count = aList_ca_count[count]
        res_start = aList_chain_start[count]
        res_end = aList_chain_end[count]

        print ' ',chain_id,'        ',atom_count,'        ',res_start,'     ',res_end

        count = count + 1

    print '\n LSQ DIAGNOSTICS\n'

    #################################################################
    # Split input file into files with independent chain fragments  #
    #################################################################

    count = 1
    while count < number_chains:

        chain_id = aList_chain[count]

        filename = 'mi_split_' + chain_id + '.pdb'

        fileexists = os.path.exists(filename)
        if fileexists != 0:
            os.remove(filename)

        file = open(filename,'w')

        for eachLine in allLines:

            eachLine = eachLine.strip()

            tag = eachLine[0:6]
            tag = tag.strip()

            if tag == 'ATOM' or tag == 'HETATM':

                chain_id_atom = eachLine[21:22]
                res_number = eachLine[22:26]
                res_number = res_number.strip()
                res_number = int(res_number)

                if chain_id_atom == chain_id:

                    # Exclude preserved zones from use as matching template in superposition calculations

                    write_flag = 'yes'

                    count_preserve = 0
                    while count_preserve < number_preserve_records:

                        chain_id_fix = aList_fix_chain[count_preserve]
                        res_start = aList_fix_res_start[count_preserve]
                        res_end = aList_fix_res_end[count_preserve]

                        if chain_id_fix == chain_id:
                            if res_number >= res_start and res_number <= res_end:
                                write_flag = 'no'

                        count_preserve = count_preserve + 1

                    if write_flag == 'yes':
                        file.write(eachLine)
                        file.write('\n')

        file.close()

        count = count + 1

    ################################################################################
    # Superimpose the targetchain onto all the independent chains with CCP4/LSQKAB #
    ################################################################################

    filename_target = 'mi_split_' + targetchain + '.pdb'

    count = 1
    while count < number_chains:

        chain_id = aList_chain[count]
        resnumber1 = aList_chain_start[count]
        resnumber2 = aList_chain_end[count]
        chain_extent = aList_ca_count[count]

        filename = 'mi_split_' + chain_id + '.pdb'
        output_file = 'mi_ncs_' + chain_id + '_from_' + targetchain + '.pdb'
        output_log = 'mi_ncs_' + chain_id + '.log'
        output_rms = 'mi_rms_analysis_' + chain_id + '.log'

        fileexists = os.path.exists(output_file)
        if fileexists != 0:
            os.remove(output_file)

        if chain_id != targetchain and chain_extent > 10:

            file = open('mi_runlsqkab.inp','w')
            file.write('FIT RESIDUE CA ')
            file.write(resnumber1)
            file.write(' TO ')
            file.write(resnumber2)

            if chain_id != ' ':
                file.write(' CHAIN ')
                file.write(targetchain)

            file.write('\nMATCH RESIDUE  ')
            file.write(resnumber1)
            file.write(' TO ')
            file.write(resnumber2)
            file.write(' CHAIN ')
            file.write(chain_id)
            file.write('\n')
            file.write('OUTPUT XYZ\n')
            file.write('END\n')
            file.close()

            # Execute

            runlsqkab = 'lsqkab xyzinm ' + filename_target + ' xyzinf ' + filename +\
                        ' xyzout ' + output_file + ' < mi_runlsqkab.inp > ' + output_log

            os.system(runlsqkab)

            # Verify and reset chain id

            fileexists = os.path.exists(output_file)
            if fileexists == 0:
                print 'Superposition with CCP4/LSQKAB failed'
                time.sleep(4)
                return 1
            else:

                # Parse rmsd and matrix from log

                file = open(output_log,'r')
                allLines = file.readlines()
                file.close()

                read_matrix = 'no'
                read_translation = 'yes'

                for eachLine in allLines:

                    # RMS

                    if eachLine.find('RMS     XYZ DISPLACEMENT') > -1:
                        aLine = eachLine.split('=')
                        rmsd = aLine[1]
                        rmsd = rmsd.strip()
                        aList_rmsd.append(rmsd)

                    # Matrix

                    if read_matrix == 'yes' and eachLine.find('TRANSLATION VECTOR') > -1 and read_translation == 'yes':
                        read_matrix = 'no'
                        read_translation = 'no'
                        aLine = eachLine.split()
                        translation = aLine[4] + ' ' + aLine[5] + ' ' + aLine[6]
                        aList_translation_matrix.append(translation)

                    if read_matrix == 'yes':
                        aLine = eachLine.strip()
                        aList_rotation_matrix.append(aLine)

                    if eachLine.find('ROTATION MATRIX:') > -1:
                        read_matrix = 'yes'

                os.remove(output_log)

                print '  Did LSQ for chain-id',chain_id,'RMSD',rmsd,'Angstroms'

                # Collect new coordinates for this chain_id

                file = open(output_file,'r')
                allLines = file.readlines()
                file.close()

                os.remove(output_file)

                for eachLine in allLines:

                    eachLine = eachLine.strip()

                    tag = eachLine[0:6]
                    tag = tag.strip()

                    if tag == 'ATOM' or tag == 'HETATM':
                        record1 = eachLine[0:21]
                        record2 = eachLine[22:80]
                        outLine = record1 + chain_id + record2

                        aList_output_atoms.append(outLine)

                file.close()

                # Compare with CCP4/COMPAR

                fileexists = os.path.exists('mi_rmstab.log')
                if fileexists != 0:
                    os.remove('mi_rmstab.log')

                fileexists = os.path.exists('mi_temp_rms.txt')
                if fileexists != 0:
                    os.remove('mi_temp_rms.txt')

                file = open(mi_compar_pdb,'w')

                number_atoms = len(aList_output_atoms)

                count_atoms = 0
                while count_atoms < number_atoms:

                    aLine = aList_output_atoms[count_atoms]

                    file.write(aLine)
                    file.write('\n')

                    count_atoms = count_atoms + 1

                file.close()

                file = open('mi_runcompar.inp','w')
                file.write('MI chain comparison\n')
                file.write('2\n')
                file.write('1.0 100\n')
                file.close()

                runcompar = 'compar XYZIN1 ' + mi_compar_pdb + ' XYZIN2 ' + filename +\
                        ' RMSTAB mi_temp_rms.txt < mi_runcompar.inp > mi_rmstab.log'

                os.system(runcompar)

                fileexists = os.path.exists('mi_rmstab.log')
                if fileexists != 0:

                    # Extract large main chain deviations

                    file = open('mi_rmstab.log','r')
                    allLines = file.readlines()
                    file.close()

                    for eachLine in allLines:
                        if eachLine.find('XYZ SHIFT >') > -1:

                            aLine = eachLine[22:38]
                            if aLine.find('CA') > -1:
                                aList_deviations.append(aLine)

                    os.remove('mi_rmstab.log')

                    fileexists = os.path.exists('mi_temp_rms.txt')
                    if fileexists != 0:  
                        os.remove('mi_temp_rms.txt')

                else:

                    print 'Analysis with CCP4/COMPAR failed'
                    time.sleep(4)
                    return 1

                # Clean-up

                fileexists = os.path.exists(mi_compar_pdb)
                if fileexists != 0:  
                    os.remove(mi_compar_pdb)

                fileexists = os.path.exists(filename)
                if fileexists != 0:   
                    os.remove(filename)

        else:

            print '  No  LSQ for chain-id',chain_id

            # Maintain original coordinates for non-LSQ chains

            file = open(filename,'r')
            allLines = file.readlines()
            file.close()

            if chain_id != targetchain:
                os.remove(filename)

            for eachLine in allLines:
                eachLine = eachLine.strip()
                aList_output_atoms.append(eachLine)

        count = count + 1

    # Clean-up

    fileexists = os.path.exists('mi_runcompar.inp')
    if fileexists != 0:   
        os.remove('mi_runcompar.inp')

    fileexists = os.path.exists('mi_runlsqkab.inp')
    if fileexists != 0:   
        os.remove('mi_runlsqkab.inp')

    fileexists = os.path.exists(filename_target)
    if fileexists != 0:   
        os.remove(filename_target)

    # Write deviations

    number_deviations = len(aList_deviations)

    if number_deviations > 0:

        print '\n LIST OF CA DEVIATIONS GREATER THAN 1.0 ANGSTROM\n'

        filename_rms = os.path.basename(pdbtargetfile)
        filename_rms = filename_rms.replace('.pdb','_deviations.txt')
        file = open(filename_rms,'w')

        count = 0
        while count < number_deviations:
            aLine = aList_deviations[count]

            file.write(aLine)
            file.write('\n')
            print aLine

            count  = count + 1

        file.close()

    else:
        print '\n NO CA ATOMS DEVIATED BY MORE THAN 1.0 ANSTROMS\n'

    number_output_atoms = len(aList_output_atoms)

    #######################################################
    # Apply any exclusions to NCS to output atom records  #
    #######################################################

    if preservemodel_file != 'none':

        # 1. For efficiency, pull records for atoms we will want to keep into a list

        number_input_atoms = len(aList_input_atoms)


        count = 0
        while count < number_input_atoms:

            aLine = aList_input_atoms[count]

            chain_id = aLine[21:22]
            res_number = aLine[22:26]
            atom_name = aLine[12:16]
            res_number = res_number.strip()
            res_number = int(res_number)
            atom_name = atom_name.strip()

            # Test if this atom record is in the preservation zone and collect it into a shorter list

            count_preserve = 0
            while count_preserve < number_preserve_records:

                chain_id_fix = aList_fix_chain[count_preserve]
                res_start = aList_fix_res_start[count_preserve]
                res_end = aList_fix_res_end[count_preserve]

                if chain_id_fix == chain_id:
                    if res_number >= res_start and res_number <= res_end:
                        aList_preserve_atom.append(aLine)
                        aList_preserve_chain.append(chain_id)
                        aList_preserve_res_number.append(res_number)
                        aList_preserve_atom_name.append(atom_name)

                count_preserve = count_preserve + 1

            count = count + 1

        number_atom_preserve_records = len(aList_preserve_atom)

        # 2. Search through NCS averaged model and put back preserved atoms

        count = 0
        while count < number_output_atoms:

            # Get parameters from this output atom record

            atom_record = aList_output_atoms[count]

            chain_id = atom_record[21:22]
            res_number = atom_record[22:26]
            atom_name = atom_record[12:16]
            res_number = res_number.strip()
            res_number = int(res_number)
            atom_name = atom_name.strip()

            # Test if this atom record is in the preservation zone and replace it

            count_preserve = 0
            while count_preserve < number_preserve_records:

                chain_id_fix = aList_fix_chain[count_preserve]
                res_start_fix = aList_fix_res_start[count_preserve]
                res_end_fix = aList_fix_res_end[count_preserve]

                if chain_id_fix == chain_id:
                    if res_number >= res_start_fix and res_number <= res_end_fix:


                        # Remove all atom records for this amino acid

                        aList_output_atoms[count] = ' '

                        # For a replaceable atom record loop through and find record to replace it with

                        count_atom_preserve_records = 0
                        while count_atom_preserve_records < number_atom_preserve_records:

                            aLine = aList_preserve_atom[count_atom_preserve_records]
                            chain_preserve = aList_preserve_chain[count_atom_preserve_records]
                            res_preserve = aList_preserve_res_number[count_atom_preserve_records]
                            atom_preserve = aList_preserve_atom_name[count_atom_preserve_records]

                            if chain_preserve == chain_id:
                                if res_preserve == res_number:
                                    if atom_preserve == atom_name:
                                        aList_output_atoms[count] = aLine

                            count_atom_preserve_records = count_atom_preserve_records + 1

                        # end loop

                count_preserve = count_preserve + 1

            #

            count = count + 1

    ############################################################
    # Check if any water molecules clash with rebuilt protein  #
    ############################################################

    aList_water_x = []
    aList_water_y = []
    aList_water_z = []
    aList_output_delete_flag = []

    number_output_atoms = len(aList_output_atoms)

    # Obtain list of waters

    count = 0
    while count < number_output_atoms:

        aLine = aList_output_atoms[count]
        tag = aLine[0:6]
        tag = tag.strip()

        if tag == 'ATOM' or tag == 'HETATM':
            res_name = aLine[30:33]
            res_name = res_name.strip()

            if res_name == 'HOH':
                water_x = aLine[30:38]
                water_y = aLine[38:46]
                water_z = aLine[46:54]
                water_x = float(water_x)
                water_y = float(water_y)
                water_z = float(water_z)
                aList_water_x.append(water_x)
                aList_water_y.append(water_y)
                aList_water_z.append(water_z)

        count = count + 1

    # Establish deletions

    number_waters = len(aList_water_x)

    count = 0
    while count < number_output_atoms:

        aList_output_delete_flag.append('no')

        aLine = aList_output_atoms[count]
        tag = aLine[0:6]
        tag = tag.strip()

        if tag == 'ATOM' or tag == 'HETATM':

            protein_x = aLine[30:38]
            protein_y = aLine[38:46]
            protein_z = aLine[46:54]
            protein_x = float(protein_x)
            protein_y = float(protein_y)
            protein_z = float(protein_z)

            count_waters = 0
            while count_waters < number_waters:

                water_x = aList_water_x[count_waters]
                water_y = aList_water_y[count_waters]
                water_z = aList_water_z[count_waters]

                dx = protein_x - water_x
                dy = protein_y - water_y
                dz = protein_z - water_z

                contact = dx*dx + dy*dy + dz*dz

                if contact < 6.25:
                    aList_output_delete_flag[count] = 'yes'

                count_waters = count_waters + 1

        count = count + 1

    ##########################################
    # Write the full model with NCS applied  #
    ##########################################

    print '\n Writing NCS rebuilt model:',filename_ncs,'\n'

    fileexists = os.path.exists(filename_ncs)
    if fileexists != 0:
        os.remove(filename_ncs)

    number_output_atoms = len(aList_output_atoms)

    file = open(filename_ncs,'w')

    if crystal_record != 'none':
        file.write(crystal_record)
        file.write('\n')

    count = 0
    while count < number_output_atoms:

        aLine = aList_output_atoms[count]
        delete_flag = aList_output_delete_flag[count]

        if aLine.find('ATOM') > -1 or aLine.find('HETATM') > -1:
            if delete_flag == 'no':
                file.write(aLine)
                file.write('\n')

        count = count + 1

    file.write('END\n')
    file.close()

    #########################
    # Map averaging option  #
    #########################

    if mtzfile == 'none':

        print '\n No MTZ file so no map averaging\n'

    else:

        print ' RUNNING MAP AVERAGING CALCULATION\n'

        if phase_comb == 'yes':
            print ' Input phases will be from phase combination\n'

        if phase_prob == 'yes':
            print ' Input phases will be from experimental phase probability\n'

        number_av_ops = len(aList_translation_matrix)

        # Copy data file and find labels

        file = open(mtzfile,'rb')
        allLines = file.readlines()
        file.close()

        file = open(mi_dm_mtz_use,'wb')
        file.writelines(allLines)
        file.close() 

        file = open('mi_mtzdump.inp','w')
        file.write('HEADER\n')
        file.write('END\n')
        file.close()

        runmtz = 'mtzdump HKLIN ' + mi_dm_mtz_use + ' < mi_mtzdump.inp > mi_mtzdump.log'
        os.system(runmtz)

        file = open('mi_mtzdump.log','r')
        allLines = file.readlines()
        file.close()

        fileexists = os.path.exists('mi_mtzdump.log')
        if fileexists != 0:
            os.remove('mi_mtzdump.log')

        fileexists = os.path.exists('mi_mtzdump.inp')
        if fileexists != 0:
            os.remove('mi_mtzdump.inp')

        read_columns = 'no'
        read_labels = 'no'
        read_cell = 'no'
        read_spacegroup = 'no'

        for eachLine in allLines:

            line_length = len(eachLine)

            # Crystal cell

            if eachLine.find('* Space group') > -1:
                aLine = eachLine.split('number')
                space_group_out = aLine[1]
                space_group = space_group_out.replace(')','')
                space_group = space_group.strip()

            if read_cell == 'yes' and line_length > 1:
                aLine = eachLine.split()
                acell = aLine[0]
                bcell = aLine[1]
                ccell = aLine[2]
                alpha = aLine[3]
                beta = aLine[4]
                gamma = aLine[5]
                read_cell = 'no'

            if eachLine.find('* Cell Dimensions :') > -1:
                read_cell = 'yes'

            # Reflection data

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

        list_length = len(labelList)

        count = 0
        while count < list_length:
            if labelList[count] == 'F' and flabel == 'none':
                if colList[count] == 'FP' or colList[count] == 'F' or colList[count] == 'Fav':
                    flabel = colList[count]

            if labelList[count] == 'Q' and sigflabel == 'none':
                if colList[count] == 'SIGFP' or colList[count] == 'SIGF' or colList[count] == 'SIGFav':
                    sigflabel = colList[count]

            if labelList[count] == 'I' and rfreelabel == 'none':
                rfreelabel = colList[count]

            if labelList[count] == 'A':
                hl = colList[count]
                aList_hl.append(hl)

            count = count + 1

        if flabel == 'none' or sigflabel == 'none':
            print '\nMTZ labels for all of F,sd(F) could not be established\n'
            time.sleep(4)
            return 1
        else:
            print ' F label    :',flabel
            print ' SD(F) label:',sigflabel

        # Check phase probabilities exist 

        if phase_comb == 'yes' or phase_prob == 'yes':

            number_hl = len(aList_hl)

            if number_hl < 4:
                print '\nHL phase probability coefficients were not found\n'
                time.sleep(4)
                return 1
            else:
                hla = aList_hl[0]
                hlb = aList_hl[1]
                hlc = aList_hl[2]
                hld = aList_hl[3]
                print ' HLA:',hla
                print ' HLB:',hlb
                print ' HLC:',hlc
                print ' HLD:',hld

        ################################################################################################
        # Do a point structure factor calculation (fixed Babinet for stability) to obtain FOM and PHIC #
        ################################################################################################

        if calc_sfs == 'yes':

            print '\n Computing structure factors from input protein entities'

            fileexists = os.path.exists(mi_refmac_pdb_use)
            if fileexists != 0:
                os.remove(mi_refmac_pdb_use)

            fileexists = os.path.exists(mi_refmac_pdb_out)
            if fileexists != 0:
                os.remove(mi_refmac_pdb_out)

            fileexists = os.path.exists(mi_refmac_mtz_out)
            if fileexists != 0:
                os.remove(mi_refmac_mtz_out)

            fileexists = os.path.exists('output.refmac')
            if fileexists != 0:
                os.remove('output.refmac')

            # Adjust coordinate file to preserve protein entities and liquify others

            count_water = 1000

            file = open(pdbtargetfile,'r')
            allLines = file.readlines()
            file.close()

            file = open(mi_refmac_pdb_use,'w')

            for eachLine in allLines:

                eachLine = eachLine.strip()

                tag = eachLine[0:6]
                tag = tag.strip()

                if tag == 'CRYST1':
                    file.write(eachLine)
                    file.write('\n')

                # Check it is protein

                if tag == 'ATOM' or tag == 'HETATM':

                    res_name = eachLine[17:20]
                    res_name = res_name.strip()

                    protein_res = 'no'

                    count = 0
                    while count < number_protein_list:
                        protein_residue = aList_protein[count]

                        if res_name == protein_residue:

                            protein_res = 'yes'

                            file.write(eachLine)
                            file.write('\n')

                        count = count + 1

                    # If not protein retain scattering as pseudo-water dummy atoms

                    if protein_res == 'no':

                        count_water = count_water + 1

                        aLine1 = eachLine[0:12]
                        aLine2 = eachLine[27:80]
                        str_count_water = str(count_water)

                        eachLine = aLine1 + ' O   HOH z' + str_count_water + aLine2

                        file.write(eachLine)
                        file.write('\n')

            # Structure factor calculation

            filename = 'mi_runrefmac5.inp'
            file = open(filename,'w')

            file.write('LABIN FP=')
            file.write(flabel)
            file.write(' SIGFP=')
            file.write(sigflabel)

            if rfreelabel != 'none':
                file.write(' FREE=')
                file.write(rfreelabel)

            if phase_comb == 'yes':
                file.write(' HLA=')
                file.write(hla)
                file.write(' HLB=')
                file.write(hlb)
                file.write(' HLC=')
                file.write(hlc)
                file.write(' HLD=')
                file.write(hld)

            file.write('\nREFI TYPE RESTrained\n')
            file.write('REFI RESI MLKF\n')
            file.write('REFI BREF ISOT METH CGMAT\n')
            file.write('SCAL TYPE BULK LSSC ANIS FIXBulk SCBULk 0.78 BBULk 180.0 \n')
            file.write('SOLVENT NO\n')
            file.write('MAKE_RESTRAINTS HYDR N\n')
            file.write('MAKE_RESTRAINTS NEWLigand Noexit\n')
            file.write('MAKE_RESTRAINTS SS Y\n')
            file.write('MAKE_RESTRAINTS CISP Y\n')
            file.write('NCYC 0\n')
            file.write('MONI DIST 6.0\n')
            file.write('MONI ANGL 8.0\n')
            file.write('MONI TORSION 10.0\n')
            file.write('MONI PLANE 10.0\n')
            file.write('MONI VANderwaals 4.25\n')
            file.write('MONI CHIRAL 8.0\n')
            file.write('MONI BFACTOR 99.0\n')
            file.write('USECWD\n')
            file.write('PNAME NOID\n')
            file.write('DNAME output_ncs\n')
            file.write('END\n')
            file.close()

            runrefmac5 = 'refmac5 HKLIN ' + mi_dm_mtz_use + ' XYZIN ' + mi_refmac_pdb_use + ' HKLOUT ' + mi_refmac_mtz_out + \
                         ' XYZOUT ' + mi_refmac_pdb_out + ' < mi_runrefmac5.inp > mi_refmac.out'

            os.system(runrefmac5)

            # Clean and check

            fileexists = os.path.exists('output_ncs.refmac')
            if fileexists != 0:            
                os.remove('output_ncs.refmac')

            fileexists = os.path.exists(mi_refmac_pdb_use)
            if fileexists != 0:             
                os.remove(mi_refmac_pdb_use)        

            fileexists = os.path.exists(mi_refmac_pdb_out)
            if fileexists == 0:
                print '\nREFMAC5 calculation failed - check mi_refmac.out\n'
                print 'The usual problem is atom names inconsistent with the PDB residue code\n'
                time.sleep(4)
                return 1
            else:
                os.remove(mi_refmac_pdb_out)            

            fileexists = os.path.exists(mi_refmac_mtz_out)
            if fileexists == 0:
                print '\nREFMAC5 calculation failed - check mi_refmac.out\n'
                print 'The usual problem is atom names inconsistent with the PDB residue code\n'
                time.sleep(4)
                return 1
            else:
                fileexists = os.path.exists(mi_dm_mtz_use)
                if fileexists != 0:           
                    os.remove(mi_dm_mtz_use)

                os.rename(mi_refmac_mtz_out,mi_dm_mtz_use)            

            fileexists = os.path.exists('mi_refmac.out')
            if fileexists != 0:
                os.remove('mi_refmac.out')

            fileexists = os.path.exists('mi_runrefmac5.inp')
            if fileexists != 0: 
                os.remove('mi_runrefmac5.inp')

            fomlabel = 'FOM'
            phaselabel = 'PHIC'

        #########################################
        # Build mask from selected coordinates  #
        #########################################

        print '\n Creating mask from target coordinates'

        found_cryst = 'no'
        number_res_in_target = 0

        file = open(pdbtargetfile,'r')
        allLines = file.readlines()
        file.close()

        file = open(mi_ncsmask_pdb,'w')

        for eachLine in allLines:

            eachLine = eachLine.strip()

            tag = eachLine[0:6]
            tag = tag.strip()

            if tag == 'CRYST1':
                file.write(eachLine)
                file.write('\n')

                found_cryst = 'yes'

            if tag == 'ATOM' or tag == 'HETATM':

                chain_id = eachLine[21:22]
                atom_name = eachLine[13:16]
                atom_name = atom_name.strip()
                res_number = eachLine[22:26]
                res_number = res_number.strip()

                # Use if in target chain

                if chain_id == targetchain:

                    if atom_name == 'CA':
                        number_res_in_target = number_res_in_target + 1

                    file.write(eachLine)
                    file.write('\n')

                # Or use if user-specified

                else:

                    count = 0
                    while count < number_ncs_extras:

                        chain_id_extra = aList_ncs_chain[count]
                        res_number_extra = aList_ncs_res_number[count]

                        if chain_id == chain_id_extra and res_number == res_number_extra:
                            file.write(eachLine)
                            file.write('\n')

                        count = count + 1

        file.close()

        if found_cryst == 'no':
            print '\nNo CRYST1 record in coordinates for map averaging\n'
            time.sleep(4)
            return 1

        # Run CCP4/NCSMASK to setup mask format using 3.5A radii about selected atoms

        fileexists = os.path.exists(mi_ncsmask_msk)
        if fileexists != 0:
            os.remove(mi_ncsmask_msk)

        file = open('mi_ncsmask.inp','w')
        file.write('RADIUS 3.5\n')
        file.write('END\n')
        file.close()

        run_ncsmask = ncsmask + ' XYZIN ' + mi_ncsmask_pdb + ' MSKOUT ' + mi_ncsmask_msk + ' < mi_ncsmask.inp > mi_ncsmask.log'
        os.system(run_ncsmask)

        fileexists = os.path.exists(mi_ncsmask_msk)
        if fileexists == 0:
            print '\nCCP4/NCSMASK calculation failed\n'
            time.sleep(4)
            return 1

        fileexists = os.path.exists(mi_ncsmask_pdb)
        if fileexists != 0:
            os.remove(mi_ncsmask_pdb)

        fileexists = os.path.exists('mi_ncsmask.inp')
        if fileexists != 0:
            os.remove('mi_ncsmask.inp')

        fileexists = os.path.exists('mi_ncsmask.log')
        if fileexists != 0:
            os.remove('mi_ncsmask.log')             

        # Obtain estimate of solvent content

        num_res_in_crystal_au = number_res_in_target * (number_av_ops + 1)
        str_num_res_in_crystal_au = str(num_res_in_crystal_au)

        filename = 'mi_runmatthews_coef.inp'
        file = open(filename, 'w')
        file.write('CELL ')
        file.write(acell)
        file.write(' ')
        file.write(bcell)
        file.write(' ')
        file.write(ccell)
        file.write(' ')
        file.write(alpha)
        file.write(' ')
        file.write(beta)
        file.write(' ')
        file.write(gamma)
        file.write(' ')
        file.write('\nSYMM ')
        file.write(space_group)
        file.write('\nNRES ')
        file.write(str_num_res_in_crystal_au)
        file.write('\n')
        file.close()

        runmatthews_coeff = 'matthews_coef < mi_runmatthews_coef.inp > mi_matthews_coef.out'
        os.system(runmatthews_coeff)

        fileexists = os.path.exists('mi_matthews_coef.out')
        if fileexists == 0:
            print '\nMATTHEWS_COEF calculation failed to execute\n'
            time.sleep(4)
            return 1

        file = open('mi_matthews_coef.out','r')
        allLines = file.readlines()
        file.close()

        for eachLine in allLines:

            if eachLine.find('Assuming protein density') > -1:

                aLine = eachLine.split(':')
                solvent_percent = aLine[1]
                solvent_percent = solvent_percent.strip()
                solvent_percent = 0.01 * float(solvent_percent)
                str_solvent_percent = str(solvent_percent)


        fileexists = os.path.exists('mi_runmatthews_coef.inp')
        if fileexists != 0:
            os.remove('mi_runmatthews_coef.inp')

        fileexists = os.path.exists('mi_matthews_coef.out')
        if fileexists != 0:
            os.remove('mi_matthews_coef.out')

        ############################## 
        # Run CCP4/DM for averaging  #
        ##############################

        print ' Running map averaging over',number_av_ops + 1,' copies'

        fileexists = os.path.exists(mi_dm_mtz_out)
        if fileexists != 0:
            os.remove(mi_dm_mtz_out)

        fileexists = os.path.exists(filename_ncs_dmlog)
        if fileexists != 0:
            os.remove(filename_ncs_dmlog)    

        file = open('mi_dm.inp','w')
        file.write('SOLC ')
        file.write(str_solvent_percent)
        file.write('\n')
        file.write('NCYC 10\n')
        file.write('MODE SOLV HIST MULT AVER\n')
        file.write('COMBINE PERT\n')

        file.write('AVER REFI\n')
        file.write('TRAN 0 0 0\n')
        file.write('ROTA MATRIX 1 0 0 0 1 0 0 0 1\n')

        count = 0
        while count < number_av_ops:

            count1 = count + 1
            count2 = count + 2
            rotation1 = aList_rotation_matrix[count]
            rotation2 = aList_rotation_matrix[count1]
            rotation3 = aList_rotation_matrix[count2]

            rotation = rotation1 + ' ' + rotation2 + ' ' + rotation3

            translation = aList_translation_matrix[count]

            file.write('AVER REFI\n')
            file.write('TRAN ')
            file.write(translation)
            file.write('\n')
            file.write('ROTA MATRIX ')
            file.write(rotation)
            file.write('\n')

            count = count + 1

        file.write('LABI FP=')
        file.write(flabel)
        file.write(' SIGFP=')
        file.write(sigflabel)

        # Add experimental phase probability information if requested

        if phase_prob == 'yes':
            file.write(' HLA=')
            file.write(hla)
            file.write(' HLB=')
            file.write(hlb)
            file.write(' HLC=')
            file.write(hlc)
            file.write(' HLD=')
            file.write(hld)
        else:
            file.write(' PHIO=')
            file.write(phaselabel)
            file.write(' FOMO=')
            file.write(fomlabel)

        file.write('\n')
        file.write('LABO FOMDM=FOMDM PHIDM=PHIDM\n')
        file.write('NOHARVEST\n')
        file.write('END\n')
        file.close()

        run_dm = 'dm HKLIN ' + mi_dm_mtz_use + ' HKLOUT ' + mi_dm_mtz_out +\
                 ' NCSIN1 ' + mi_ncsmask_msk + ' < mi_dm.inp > ' + filename_ncs_dmlog

        os.system(run_dm)

        fileexists = os.path.exists(mi_dm_mtz_use)
        if fileexists != 0:        
            os.remove(mi_dm_mtz_use)

        fileexists = os.path.exists(mi_ncsmask_msk)
        if fileexists != 0:          
            os.remove(mi_ncsmask_msk)    

        fileexists = os.path.exists(mi_dm_mtz_out)
        if fileexists == 0:
            print '\nCCP4/DM calculation failed\n'
            time.sleep(4)
            return 1

        fileexists = os.path.exists('mi_dm.inp')
        if fileexists != 0:          
            os.remove('mi_dm.inp')        

        # Filter and set final output file

        fileexists = os.path.exists(filename_ncs_mtz)
        if fileexists != 0:
            os.remove(filename_ncs_mtz)

        file = open('mi_cad.inp','w')
        file.write('LABIN FILE_NUMBER 1 E1=')
        file.write(flabel)
        file.write(' E2=')
        file.write(sigflabel)
        file.write(' E3=FOMDM E4=PHIDM\n')
        file.write('END\n')
        file.close()

        runcad = 'cad HKLIN1 ' + mi_dm_mtz_out + ' HKLOUT ' + filename_ncs_mtz + ' < mi_cad.inp > mi_cad.log'
        os.system(runcad)

        fileexists = os.path.exists(mi_dm_mtz_out)
        if fileexists != 0:
            os.remove(mi_dm_mtz_out)

        fileexists = os.path.exists('mi_cad.log')
        if fileexists != 0:
            os.remove('mi_cad.log')

        fileexists = os.path.exists('mi_cad.inp')
        if fileexists != 0:
            os.remove('mi_cad.inp')    

        fileexists = os.path.exists(filename_ncs_mtz)
        if fileexists != 0:

            # Check file size (bytes) to be sure it really worked

            filesize = os.path.getsize(filename_ncs_mtz)

            if filesize < 1000:
                print 'Output file does not seem populated - removing'
                os.remove(filename_ncs_mtz)
                time.sleep(4)
                return 1

        else:
            print 'The CAD run to filter the data seems to have failed'
            time.sleep(4)
            return 1

        #

        print '\n Writing NCS averaged data:',filename_ncs_mtz,'\n'
        print ' FOMs from averaging  : FOMDM'
        print ' Phases from averaging: PHIDM\n'
        print ' CCP4/DM log file:', filename_ncs_dmlog

    # View

    time.sleep(2)
    return 0

if __name__ == "__main__":
    sys.exit(Run())
