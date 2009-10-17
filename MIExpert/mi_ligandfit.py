#####################################################################
#                                                                   #
# Example of single conformer ligand density docking with           #
# CCP4/FFFEAR. Intended to be called from mi_bng.py                 #
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
import dircache
import getopt
import ccp4check

def Usage():
    print "Usage:"
    print " %s -p pdbfile -m mtzfile -f fragfitfile -o fragfitfile_out \\" % sys.argv[0]
    print "    -c frag_center -d workingdir \\"
    print "    [-? or --help] for this help"

def Run(argv=None):
    if argv is None:
        argv=sys.argv
    # Initialize

    pdbfile = 'none'
    mtzfile = 'none'
    workingdir = 'none'
    fragfitfile = 'none'
    fragfitfile_out = 'none'
    frag_center = 'none'

    autoligand_res_number = 0

    # parse args
    number_of_args = len(argv)
    args = argv[1:]
    optlist, args = getopt.getopt(args,'p:m:d:f:o:c:?',['help'])
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
        if number_of_list_inputs >=2:
            param_value = aList[1]
            
        if arg_value == '-p':
            pdbfile = param_value
        if arg_value == '-m':
            mtzfile = param_value
        if arg_value == '-d':
            workingdir = param_value
        if arg_value == '-f':
            fragfitfile = param_value
        if arg_value == '-o':
            fragfitfile_out = param_value
        if arg_value == '-c':
            frag_center = param_value
        if arg_value == '-?' or arg_value=="--help":
            Usage()
            return
        count=count+1

    ccp4,error = ccp4check.ccp4check()
    if not ccp4:
      print '\n' + error + '\n'
      time.sleep(4)
      return 1
    
    # verify parameters
    fileexists = os.path.exists(workingdir)
    if fileexists == 0:
        print 'The working directory was not found ',workingdir
        time.sleep(4)
        return 1
    else:
        # Go to working area - PDB and MTZ should be local to this directory, ligand file may not be
        os.chdir(workingdir)

    fileexists = os.path.exists(pdbfile)
    if fileexists == 0:
        print 'The PDB file was not found ',pdbfile
        time.sleep(4)
        return 1

    fileexists = os.path.exists(mtzfile)
    if fileexists == 0:
        print 'The MTZ file was not found ',mtzfile
        time.sleep(4)
        return 1

    fileexists = os.path.exists(fragfitfile)
    if fileexists == 0:
        print 'The ligand file was not found ',fragfitfile
        time.sleep(4)
        return 1                 

    if frag_center == 'none':
        print 'The fragment target point was not given'
        time.sleep(4)
        return 1

    if fragfitfile_out == 'none':
        print 'The name for the output coordinate file was not given'
        time.sleep(4)
        return 1

    #################
    # Ligand search #
    #################

    # Obtain the local name for the input data file

    mtzfile_local_ligand = os.path.basename(mtzfile)

    # Check through protein file to determine a suitable chain-resno for the ligand
    # Set new ligand to Z 1, add additional ligands to the Z-chain

    file = open(pdbfile,'r')
    allLines = file.readlines()
    file.close()

    for eachLine in allLines:
        tag = eachLine[0:6]
        tag = tag.strip()

        if tag == 'ATOM' or tag == 'HETATM':

            chain = eachLine[21:22]
            chain = chain.strip()
            res_number = eachLine[22:26]
            res_number = res_number.strip()

            if chain == 'Z':
                autoligand_res_number = int(res_number)

    # Check through ligand fragment file to set B-factor and get ligand name

    file = open(fragfitfile,'r')
    allLines = file.readlines()
    file.close()

    num_lig_atoms = 0

    file = open('mi_frag.pdb','w')    

    for eachLine in allLines:

        tag = eachLine[0:6]
        tag = tag.strip()

        if tag == 'ATOM' or tag == 'HETATM':

            lig_res_name = eachLine[17:20]
            lig_res_name = lig_res_name.strip()

            atom = eachLine[0:61] + '20.00' + eachLine[66:80]
            atom = atom.strip()
            file.write(atom)
            file.write('\n')

            num_lig_atoms = num_lig_atoms + 1        

    file.write('END\n')
    file.close()  

    # Check it is a ligand

    if num_lig_atoms > 100:
        print 'This ligand contains more than 100 atoms - seems to be too many for a ligand'
        time.sleep(4)
        return 1

    # Setup atom record

    autoligand_res_number = autoligand_res_number + 1

    if autoligand_res_number < 10:
        autoligand_label = lig_res_name + ' Z   ' + str(autoligand_res_number)
    else:
        autoligand_label = lig_res_name + ' Z  ' + str(autoligand_res_number)

    # Resolution needed for effective search depends on ligand size

    if num_lig_atoms < 16:
        search_angle = '15'
        search_res = '2.4'
    else:
        search_angle = '20'
        search_res = '2.8'

    # Coarse 6D convolution search with CCP4/FFFEAR against difference map coefficients

    fileexists = os.path.exists('mi_frag_out_0.pdb')
    if fileexists != 0:
        os.remove('mi_frag_out_0.pdb')

    file = open('mi_fffear.inp','w')
    file.write('SOLC 0.95\n')
    file.write('SEARCH STEP ')
    file.write(search_angle)
    file.write(' PEAKS 1\n')
    file.write('SCALE NATURAL\n')
    file.write('MODEL RADIUS 2.5\n')
    file.write('RESO 1000.0 ')
    file.write(search_res)
    file.write('\nCENTRE ORTH ')
    file.write(frag_center)
    file.write('\nLABIN FP=DELFWT SIGFP=SIGFP PHIO=PHDELFWT\n')
    file.write('END\n')
    file.close()

    runfffear = 'fffear HKLIN ' + mtzfile_local_ligand + ' XYZIN mi_frag.pdb XYZOUT mi_frag_out_0.pdb \
    MAPOUT mi_fffear.map < mi_fffear.inp > mi_fffear_0.log'

    os.system(runfffear)

    fileexists = os.path.exists('mi_frag_out_0.pdb')
    if fileexists == 0:        
        print 'Ligand fitting process failed'
        time.sleep(4)
        return 1

    # Reset ligand name and B-factors to initial value

    file = open('mi_frag_out_0.pdb','r')
    allLines = file.readlines()
    file.close()

    os.remove('mi_frag_out_0.pdb')

    file = open('mi_frag_out_0.pdb','w')    

    for eachLine in allLines:

        tag = eachLine[0:6]
        tag = tag.strip()

        if tag == 'ATOM' or tag == 'HETATM':
            atom_out = eachLine.replace('ALA A   1',autoligand_label)
            atom = atom_out[0:61] + '20.00' + atom_out[66:80]
            atom = atom.strip()
            file.write(atom)
            file.write('\n')

    file.write('END\n')
    file.close()  

    # Improve the first solution with a fine 6D convolution search

    print 'Optimizing ligand fit'

    fileexists = os.path.exists('mi_frag_out_1.pdb')
    if fileexists != 0:
        os.remove('mi_frag_out_1.pdb')

    file = open('mi_fffear.inp','w')
    file.write('SOLC 0.95\n')
    file.write('SEARCH STEP 3 RANGE 12 PEAKS 1\n')
    file.write('SCALE NATURAL\n')
    file.write('MODEL RADIUS 2.5\n')
    file.write('RESO 1000.0 ')
    file.write(search_res)
    file.write('\nCENTRE ORTH ')
    file.write(frag_center)
    file.write('\nLABIN FP=DELFWT SIGFP=SIGFP PHIO=PHDELFWT\n')
    file.write('END\n')
    file.close()

    runfffear = 'fffear HKLIN ' + mtzfile_local_ligand + ' XYZIN mi_frag_out_0.pdb XYZOUT mi_frag_out_1.pdb \
    MAPOUT mi_fffear.map < mi_fffear.inp > mi_fffear_1.log'               

    os.system(runfffear)

    fileexists = os.path.exists('mi_frag_out_1.pdb')
    if fileexists == 0:
        print 'Ligand fitting optimization process failed'
        time.sleep(4)
        return 1 

    ################################################
    # Splice solution into protein coordinate file #
    ################################################

    file = open(pdbfile,'r')
    allLines_protein = file.readlines()
    file.close()

    file = open('mi_frag_out_1.pdb','r')
    allLines_ligand = file.readlines()
    file.close()

    ligand_write = 'no'

    file = open(fragfitfile_out,'w')

    for eachLine_protein in allLines_protein:

        tag = eachLine_protein[0:6]
        tag = tag.strip()

        # Re write selected PDB file items

        if tag.find('HEADER') > -1 or tag.find('CRYST1') > -1 or tag.find('SCALE') > -1 or tag.find('TER') > -1 \
            or tag.find('HETATM') > -1 or tag.find('ATOM') > -1:

            # Put ligand 

            if eachLine_protein.find('HOH') > -1 and ligand_write != 'done':

                for eachLine_ligand in allLines_ligand:
                    if eachLine_ligand.find('ATOM') > -1:                       
                        atom = str(eachLine_ligand)

                        # Fix residue name and chain-id/res-number and reset B-factor

                        atom_out = atom.replace('ALA A   1',autoligand_label)
                        atom = atom_out[0:61] + '20.00' + atom_out[66:80]
                        atom = atom.strip()                            
                        file.write(atom)
                        file.write('\n')

                ligand_write = 'done'

            file.write(eachLine_protein)

    # Write ligand for no-water case

    if ligand_write != 'done':

        for eachLine_ligand in allLines_ligand:
            if eachLine_ligand.find('ATOM') > -1:                       
                atom = str(eachLine_ligand)

                # Fix residue name and chain-id/res-number and reset B-factor

                atom_out = atom.replace('ALA A   1',autoligand_label)
                atom = atom_out[0:61] + '20.00' + atom_out[66:80]
                atom = atom.strip()
                file.write(atom)
                file.write('\n')

    # Close out

    file.write('END\n')
    file.close()

    ############
    # Clean up #
    ############

    fileexists = os.path.exists('mi_frag_out_0.pdb')
    if fileexists != 0:    
        os.remove('mi_frag_out_0.pdb')

    fileexists = os.path.exists('mi_frag_out_1.pdb')
    if fileexists != 0:    
        os.remove('mi_frag_out_1.pdb')

    fileexists = os.path.exists('mi_fffear.inp')
    if fileexists != 0:    
        os.remove('mi_fffear.inp')

    fileexists = os.path.exists('mi_frag.pdb')
    if fileexists != 0:
        os.remove('mi_frag.pdb')

    fileexists = os.path.exists('mi_fffear_0.log')
    if fileexists != 0:
        os.remove('mi_fffear_0.log')

    fileexists = os.path.exists('mi_fffear_1.log')
    if fileexists != 0:
        os.remove('mi_fffear_1.log')

    fileexists = os.path.exists('mi_fffear.map')
    if fileexists != 0:
        os.remove('mi_fffear.map')

    fileexists = os.path.exists('mi_fffear_out.mtz')
    if fileexists != 0:
        os.remove('mi_fffear_out.mtz')  
    return 0


if __name__ == "__main__":
    sys.exit(Run())






