#####################################################################
#                                                                   #
# Overlap and compare ligands in related structures                 # 
#                                                                   #
# Copyright: Molecular Images   2005                                #
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
    print "  -p,--targetpdb=FILE         target pdb file" 
    print "  -d,--workdir=DIR            the working directory to use"
    print "  -t,--targetsite=\"x y z\"     center of site of interest"
    print "     --pdbdir=DIR             directory containing the pdb files"
    print "  -?,--help                   this help file"

def Run(argv=None):
    if argv is None:
        argv=sys.argv
    # Initialize

    projectlog = 'project_history.txt'
    workingdir = 'none'
    pdbtargetfile = 'none'
    targetsite = 'none'
    pdbList = []
    aList_chain = []
    aList_resnumber = []
    aList = []

    site_range = 100.0

    # Standard non-ligand

    aList_protein = ['ALA','ARG','ASN','ASP','CYS','GLN','GLU','GLY','HIS','ILE','LEU',\
                     'LYS','MET','PHE','PRO','SER','THR','TRP','TYR','VAL','MSE','PTR','SEP','TPO','HOH']

    number_protein_list = len(aList_protein)

    aList_new_chain = ['a','b','c','d','e','f','g','h','i','j','k','l','m','n','o','p','q','r','s','t','u'\
                      ,'v','w','x','y','z','A','B','C','D','E','F','G','H','I','J','K','L','M','N','O'\
                      ,'P','Q','R','S','T','U','V','W','X','Y','Z']

    number_new_chain = len(aList_new_chain)

    # parse args
    number_of_args = len(argv)
    args = argv[1:]
    optlist, args = getopt.getopt(
        args,'p:d:t:?',
        ['targetpdb=','workdir=','pdbdir=','targetsite=','help'])
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
            if arg_value == '-p' or arg_value=='--targetpdb':
                pdbtargetfile = param_value
            elif arg_value == '-d' or arg_value=='--workdir':
                workingdir = param_value
            elif arg_value == '--pdbdir':
                pdblist = param_value
            elif arg_value == '--targetsite' or arg_value=='-t':
                targetsite = param_value
        count=count+1

    ccp4,error = ccp4check.ccp4check()
    if not ccp4:
      print '\n' + error + '\n'
      time.sleep(4)
      return 1
    
    fileexists = os.path.exists(pdbtargetfile)
    if fileexists == 0:
        print 'The target PDB file was not found ',pdbtargetfile
        time.sleep(4)
        return 1

    fileexists = os.path.exists(workingdir)
    if fileexists == 0:
        print 'The working directory was not found ',workingdir
        time.sleep(4)
        return 1

    fileexists = os.path.exists(pdblist)
    if fileexists == 0:
        print 'The directory containing pdb files was not found',pdblist
        time.sleep(4)
        return 1

    if targetsite == 'none':
        print 'The target XYZ coordinates must be given'
        time.sleep(4)
        return 1
    else:

        aString = targetsite.split()
        num_markers = len(aString)

        if num_markers != 3:
            print 'The target site parameter must contain 3 numbers',targetsite
            time.sleep(4)
            return 1
        else:
            xmarker = aString[0]
            ymarker = aString[1]
            zmarker = aString[2]

            xmarker = float(xmarker)
            ymarker = float(ymarker)
            zmarker = float(zmarker)

    # Go to working area

    os.chdir(workingdir)

    ##############################################
    # Create the list of target coordinate files #
    ##############################################

    aList = os.listdir(pdblist)
    num_files = len(aList)

    # Extract coordinate files

    count = 0
    while count < num_files:

        pdb_file = aList[count]
        if pdb_file.find('.pdb') > -1:
            pdbList.append(pdb_file)

        count = count + 1

    num_files = len(pdbList)

    if num_files < 2:
        print 'There must be more than one .pdb file in the working directory'
        time.sleep(4)
        return 1

    ######################################################
    # Obtain the protein CA atoms around the target site #
    ######################################################

    file = open(pdbtargetfile,'r')
    allLines = file.readlines()
    file.close()

    for eachLine in allLines:

        tag = eachLine[0:6]
        tag = tag.strip()

        if tag == 'ATOM' or tag == 'HETATM':

            chain_id = eachLine[21:22]
            res_number = eachLine[22:26]
            res_name = eachLine[17:20]
            res_name = res_name.strip()
            atom_name = eachLine[13:16]
            atom_name = atom_name.strip()

            if chain_id == ' ':
                chain_id = '.'

            if atom_name == 'CA':

                protein = 'no'
                count = 0
                while count < number_protein_list:
                    res_name_standard = aList_protein[count]

                    if res_name_standard == res_name:
                        protein = 'yes'

                    count = count + 1

                if protein == 'yes':

                    xcoord = eachLine[30:38]
                    ycoord = eachLine[39:46]
                    zcoord = eachLine[47:54]
                    xcoord = float(xcoord)
                    ycoord = float(ycoord)
                    zcoord = float(zcoord)

                    dx = xcoord - xmarker
                    dy = ycoord - ymarker
                    dz = zcoord - zmarker

                    dxsq = dx*dx + dy*dy + dz*dz

                    if dxsq < site_range:
                        aList_chain.append(chain_id)
                        aList_resnumber.append(res_number)

    num_cas = len(aList_resnumber)

    if num_cas < 4:
        print 'Fewer than four CA atoms surround the target site'
        time.sleep(4)
        return 1    

    # Summarize

    print '\nMultiple structure superposition based on target site CA atoms\n'
    print 'Number of structures to superimpose: ',num_files
    print 'Number of CA atoms used for overlap: ',num_cas

    ######################################################
    # Do multiple superposition, writing a file for each #
    ######################################################

    # Copy target file to working directory as mi_rot_0.pdb

    file = open(pdbtargetfile,'r')
    allLines = file.readlines()
    file.close()

    file = open('mi_rot_0.pdb','w')

    for eachLine in allLines:

        tag = eachLine[0:6]
        tag = tag.strip()

        if tag == 'ATOM' or tag == 'HETATM' or tag == 'TER' or tag == 'END':
            file.write(eachLine)

    file.close()

    # Create CCP4/LSQKAB input

    file = open('mi_runlsqkab.inp','w')

    count = 0
    while count < num_cas:

        resnumber = aList_resnumber[count]
        chainid = aList_chain[count]

        file.write('FIT RESIDUE CA ')
        file.write(resnumber)
        file.write(' TO ')
        file.write(resnumber)

        if chainid != '.':
            file.write(' CHAIN ')
            file.write(chainid)

        file.write('\nMATCH RESIDUE  ')
        file.write(resnumber)
        file.write(' TO ')
        file.write(resnumber)
        file.write(' CHAIN ')
        file.write(chainid)
        file.write('\n')

        count = count + 1

    file.write('OUTPUT XYZ\n')
    file.write('END\n')
    file.close()

    # Loop over pdb file targets

    count = 0
    while count < num_files:

        input_file = pdbList[count]
        input_file_full = os.path.join(pdblist,input_file)

        # Make local copy

        file = open(input_file_full,'r')
        allLines = file.readlines()
        file.close()

        file = open('mi_move.pdb','w')

        for eachLine in allLines:

            tag = eachLine[0:6]
            tag = tag.strip()

            if tag == 'ATOM' or tag == 'HETATM' or tag == 'TER' or tag == 'END':
                file.write(eachLine)

        file.close()

        # Name output file

        count_id = count + 1
        count_str = str(count_id)
        output_file = 'mi_rot_' + count_str + '.pdb'

        fileexists = os.path.exists(output_file)
        if fileexists != 0:
            os.remove(output_file)

        # Execute

        runlsqkab = 'lsqkab xyzinm mi_move.pdb xyzinf mi_rot_0.pdb xyzout ' + output_file + \
                    ' < mi_runlsqkab.inp 1> mi_runlsqkab.log 2> mi_runlsqkab_errors.log'

        os.system(runlsqkab)

        # Verify

        fileexists = os.path.exists(output_file)
        if fileexists == 0:
            print 'Superposition with LSQKAB failed for',input_file_full
            time.sleep(4)

        fileexists = os.path.exists('mi_runlsqkab_errors.log')
        if fileexists != 0:
            os.remove('mi_runlsqkab_errors.log')

        os.remove('mi_move.pdb')

        count = count + 1

    os.remove('mi_runlsqkab.log')
    os.remove('mi_runlsqkab.inp')
    os.remove('mi_rot_0.pdb')

    ##############################################
    # Capture all the ligands into a single file #
    ##############################################

    print 'All superimposed ligands are loaded into file: allligands.pdb\n'

    # Ligand master file

    filemaster = open('allligands.pdb','w')

    count = 0
    while count < num_files:

        # Obtain contents of superimposed file

        count_id = count + 1
        count_str = str(count_id)
        output_file = 'mi_rot_' + count_str + '.pdb'

        fileexists = os.path.exists(output_file)
        if fileexists != 0:
            file = open(output_file,'r')
            allLines = file.readlines()
            file.close()

            # Obtain only ligands

            for eachLine in allLines:

                tag = eachLine[0:6]
                tag = tag.strip()

                if tag == 'ATOM' or tag == 'HETATM':

                    res_name = eachLine[17:20]
                    res_name = res_name.strip()

                    protein_flag = 'no'

                    count1 = 0
                    while count1 < number_protein_list:

                        protein_entity = aList_protein[count1]

                        if protein_entity == res_name:
                            protein_flag = 'yes'

                        count1 = count1 + 1

                    # For ligands distinguish by chain-ids

                    if protein_flag == 'no':

                        if count_id < number_new_chain:
                            new_chain = aList_new_chain[count]
                        else:
                            new_chain = '1'

                        # Write coordinates to master file

                        atom_out = eachLine[0:21] + new_chain + eachLine[22:80]
                        filemaster.write(atom_out)
                        filemaster.write('\n')

            os.remove(output_file)

        # End of loop over this file

        count = count + 1

    filemaster.write('END')
    filemaster.close()

    #
    return 0           

            
            
if __name__ == "__main__":
    sys.exit(Run())
