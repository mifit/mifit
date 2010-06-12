#####################################################################
#                                                                   #
# Refinement Script with SHELX and CCP4/REFMAC5                     #
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
import dircache
import getopt
import ccp4check

def Usage():
    print "Usage: %s [options]" % sys.argv[0]
    print "Options are:"
    print "  -p,--pdbfile=FILE           the pdb file"
    print "  -m,--mtzfile=FILE           the mtz file"
    print "  -l,--libfile=FILE           the library file. Default: no file"
    print "  -d,--workdir=DIR            The working directory"
    print "  -e,--engine=ENGINE          One of refmac5 (default), shelx, or rigid"
    print "  -w,--weight=NUM             The weighting factor. Default: 0.1"
    print "  -c,--cycles=NUM             Number of refinement cycles to run"
    print "     --water_cycles=NUM       Number of water cycles to run. Default: 0"
    print "  -t,--tls_file=FILE          TLS specification file. Default: no file"
    print "  -s,--shelx_dir=DIR          Path to shelx executables. Default: $SHELXBIN"
    print "  -h,--mifithome=DIR          Path to MIFit. Default: no path"
    print "     --bref_type=TYPE         B-factor refinement type: anisotropic or none. Default: none" 
    print "     --max_res=NUM            Maximum resolution. Default: no value"
    print "  -?,--help                   this help"

def Run(argv=None):
    if argv is None:
        argv=sys.argv

    # Path to MIFit installation to find phi-psi data

    mifit_root = 'none'

    # Initialize

    quote = """'"""

    job_prefix = 'refine_'
    pdbfile = 'none'
    mtzfile = 'none'
    libfile = 'none'
    workingdir = 'none'
    runid = '1'
    projectlog = 'project_history.txt'
    ref_engine = 'refmac5'
    flabel = 'none'
    sigflabel = 'none'
    rfreelabel = 'none'
    anomlabel = 'none'
    siganomlabel = 'none'
    resolution_mtz = 'none'
    rwork = 'none'
    rfree = 'none'
    rmsd_bonds = 'none'
    percent_phi_psi = 'none'
    percent_rotamer = 'none'
    shelxpro = 'none'
    shelxh = 'none'

    max_res = 'none'
    weight = 'none'
    bref_type = 'none'
    freeflag = '0'
    number_molecules = '1'
    res_number_X_high = 0
    int_res_number = 0
    number_sym_ops = 1
    cycles = '5'
    water_pick = 'no'
    water_cycles = 0
    big_cycles = 1
    water_count = 0
    missing_protein_chain = 'no'
    validate = 'yes'
    max_conformers = 10
    disorder = 'no'
    extra_links = 'no'
    resolution_output = '?'
    tlsfile = 'none'
    shelx_directory = 'none'
 
    seq_chain_prev = '?'

    filename_log_full = 'none'
    filename_pdb_full = 'none'
    filename_mtz_full = 'none'
    filename_refmac_full = 'none'
    filename_anom_full = 'none'

    number_chain_list = 0
    aList_chains = []
    aList_nterm = []
    aList_cterm = []
    parseLine = []
    aLine = []
    labelList = []
    colList = []

    aList_chain_store = []
    aList_res_number_store = []
    aList_res_name_store = []

    aList_sequence_chain_id = []
    pdb_annotate = []
    aList_sequence_chain = []
    aList_sequence_number = []
    aList_sequence_resname = []
    aList_missing_residues = []
    aList_current_residue_atoms = []
    aList_allatoms_chain = []
    aList_allatoms_res_number = []
    aList_allatoms_res_name = []
    aList_allatoms_atom_name = []
    aList_SEQRES = []

    aList_rotamer_chain = []
    aList_rotamer_resno = []
    aList_rotamer_resname = []
    aList_bonds_chain = []
    aList_bonds_resno = []
    aList_bonds_resname = []
    aList_angles_chain = []
    aList_angles_resno = []
    aList_angles_resname = []
    aList_contacts_chain = []
    aList_contacts_resno = []
    aList_contacts_resname = []
    aList_chiral_chain = []
    aList_chiral_resno = []
    aList_chiral_resname = []
    aList_cis_chain = []
    aList_cis_resno = []
    aList_cis_resname = []
    aList_rama_chain = []
    aList_rama_resno = []
    aList_rama_resname = []
    aList_omega_chain = []
    aList_omega_resno = []
    aList_omega_resname = []
    aList_density_chain = []
    aList_density_resno = []
    aList_density_resname = []
    aList_disorder_chain = []
    aList_disorder_resno = []
    aList_disorder_resname = []
    aList_errors = []
    aList_phi_all = []
    aList_psi_all = []
    aList_phipsi_prob_all = []
    aList_phi_gly = []
    aList_psi_gly = []
    aList_phipsi_prob_gly = []
    aList_phi_pro = []
    aList_psi_pro = []
    aList_phipsi_prob_pro = []
    aList_peak_x = []
    aList_peak_y = []
    aList_peak_z = []
    aList_x = []
    aList_y = []
    aList_z = []

    bond_list = 'no'
    angle_list = 'no'
    contact_list = 'no'
    chiral_list = 'no'
    iteration_final = 'no'
    read_error_log = 'no'
    chain_id_prev = '?'
    res_number_prev = '?'
    num_residues = 0.0
    amino_acid_count = 0.0
    count_phipsi = 0.0
    count_rotamer = 0.0
    phipsi_gen_datafile = 'none'
    phipsi_gly_datafile = 'none'
    phipsi_pro_datafile = 'none'
    number_phipsi_gen_table = 0
    number_phipsi_gly_table = 0
    number_phipsi_pro_table = 0

    # Omega threshhold is 4 x true sd from peak at 178.9
    # allowed phi-psi: gen 99.95% data, 41.5% area, GLY 99.8% data, 63% area, PRO 99.8% data, 18.1% area

    omega_peak = 178.9
    omega_thresh = 4.0 * 5.6
    phipsi_thresh_gen = 0.00847
    phipsi_thresh_gly = 0.00384
    phipsi_thresh_pro = 0.0015

    # Atom lists

    aList_GLY_atoms = ['N','CA','C','O']
    aList_ALA_atoms = ['N','CA','C','O','CB']
    aList_VAL_atoms = ['N','CA','C','O','CB','CG1','CG2']
    aList_ILE_atoms = ['N','CA','C','O','CB','CG1','CG2','CD1']
    aList_LEU_atoms = ['N','CA','C','O','CB','CG','CD1','CD2']
    aList_PHE_atoms = ['N','CA','C','O','CB','CG','CD1','CE1','CZ','CE2','CD2']
    aList_PRO_atoms = ['N','CA','C','O','CB','CG','CD']
    aList_MET_atoms = ['N','CA','C','O','CB','CG','SD','CE']
    aList_TRP_atoms = ['N','CA','C','O','CB','CG','CD1','NE1','CE2','CD2','CE3','CZ3','CH2','CZ2']
    aList_CYS_atoms = ['N','CA','C','O','CB','SG']
    aList_SER_atoms = ['N','CA','C','O','CB','OG']
    aList_THR_atoms = ['N','CA','C','O','CB','OG1','CG2']
    aList_ASN_atoms = ['N','CA','C','O','CB','CG','OD1','ND2']
    aList_GLN_atoms = ['N','CA','C','O','CB','CD','CG','OE1','NE2']
    aList_TYR_atoms = ['N','CA','C','O','CB','CG','CD1','CE1','CZ','OH','CE2','CD2']
    aList_HIS_atoms = ['N','CA','C','O','CB','CG','ND1','CE1','NE2','CD2']
    aList_ASP_atoms = ['N','CA','C','O','CB','CG','OD1','OD2']
    aList_GLU_atoms = ['N','CA','C','O','CB','CD','CG','OE1','OE2']
    aList_LYS_atoms = ['N','CA','C','O','CB','CG','CD','CE','NZ']
    aList_ARG_atoms = ['N','CA','C','O','CB','CD','CG','NE','CZ','NH1','NH2']

    # Platform

    test_platform = sys.platform

    # Read args
    number_of_args = len(argv)
    args = argv[1:]
    optlist, args = getopt.getopt(
        args,'p:m:l:d:e:w:c:t:s:h:?',
        ["pdbfile=","mtzfile=","libfile=","workdir=","engine=",
         "weight=","cycles=","water_cycles=",
         "tls_file=","shelx_dir=","mifithome=",
         "bref_type=","max_res=","help"])
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
            if arg_value == '-p' or arg_value=="--pdbfile":
                pdbfile = param_value
            elif arg_value == '-m' or arg_value=="--mtzfile":
                mtzfile = param_value
            elif arg_value == '-l' or arg_value=="--libfile":
                libfile = param_value
            elif arg_value == '-d' or arg_value=="--workdir":
                workingdir = param_value
            elif arg_value == '-e' or arg_value=="--engine":
                ref_engine = param_value
            elif arg_value == '-w' or arg_value=="--weight":
                weight = param_value
            elif arg_value == '--max_res':
                max_res = param_value
            elif arg_value == '-bref_type':
                bref_type = param_value
            elif arg_value == '-c' or arg_value=="--cycles":
                cycles = param_value
            elif arg_value == '--water_cycles':
                water_cycles = int(param_value)
            elif arg_value == '-t' or arg_value=="--tls_file":
                tlsfile = param_value
            elif arg_value == '-s' or arg_value=="--shelx_dir":
                shelx_directory = param_value
            elif arg_value == '-h' or arg_value=="--mifithome":
                mifit_root = param_value
        count = count + 1

    ccp4,error = ccp4check.ccp4check()
    if not ccp4:
      print '\n' + error + '\n'
      time.sleep(4)
      return 1

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

    fileexists = os.path.exists(workingdir)
    if fileexists == 0:
        print 'The working directory was not found ',workingdir
        time.sleep(4)
        return 1

    fileexists = os.path.exists(libfile)
    if fileexists == 0 and os.path.basename(libfile) != 'none':
        print 'The library file was not found ',libfile
        time.sleep(4)
        return 1                    

    fileexists = os.path.exists(tlsfile)
    if fileexists == 0 and os.path.basename(tlsfile) != 'none':
        print 'The TLS specification file was not found ',tlsfile
        time.sleep(4)
        return 1

    if ref_engine != 'shelx' and ref_engine != 'rigid' and ref_engine != 'refmac5':
        print 'The refinement type must be one of shelx/rigid/refmac5'
        time.sleep(4)
        return 1

    if weight == 'none':
        weight = '0.1'

    if water_cycles > 0:
        water_pick = 'yes'
        big_cycles = water_cycles + 1 

    if bref_type == 'none':
        bref_type = 'isotropic'

    # Check MIFit installation to access phi-psi data files (try environment variable, default and direct path)

    if mifit_root != 'none':
        mifit_root_data = os.path.join(mifit_root,'data')
        phipsi_gen_datafile = os.path.join(mifit_root_data,'rama500-general.data')
        phipsi_gly_datafile = os.path.join(mifit_root_data,'rama500-gly-sym-nosec.data')
        phipsi_pro_datafile = os.path.join(mifit_root_data,'rama500-pro.data')

    # If needed, get determine SHELX bin directory path

    if ref_engine == 'shelx':

        # Determine installation from a possible environment variable

        find_shelx_bin = os.environ.keys().count('SHELXBIN')
        if find_shelx_bin != 0:
            shelx_directory = os.environ['SHELXBIN']

        # Determine from input parameter

        if shelx_directory != 'none':

            fileexists = os.path.exists(shelx_directory)
            if fileexists != 0:

                if test_platform.find('win') > -1:
                    shelxpro = os.path.join(shelx_directory,'shelxpro.exe')
                    shelxh = os.path.join(shelx_directory,'shelxh.exe')
                else:
                    shelxpro = os.path.join(shelx_directory,'shelxpro')
                    shelxh = os.path.join(shelx_directory,'shelxh')

            # Confirm we have SHELXH and SHELXPRO

            fileexists = os.path.exists(shelxpro)
            if fileexists == 0:
                print 'The SHELXPRO executable was not found ',shelxpro
                time.sleep(4)
                return 1

            fileexists = os.path.exists(shelxh)
            if fileexists == 0:
                print 'The SHELXH executable was not found ',shelxh
                time.sleep(4)
                return 1

    # Create a CCP4 temp space for temporary files

    idcode = '000000'

    gmt = time.gmtime(time.time())
    fmt = '%H%M%S'
    idcode = time.strftime(fmt,gmt)

    path_scratch = 'temp_' + idcode

    working_ccp4_scratch = os.path.join(ccp4.scr,path_scratch)

    fileexists = os.path.exists(working_ccp4_scratch)
    if fileexists == 0:
        os.mkdir(working_ccp4_scratch)

    os.environ['CCP4_SCR'] = working_ccp4_scratch

    # Go to working area

    os.chdir(workingdir)

    # Copy MTZ file to working area

    file = open(mtzfile,'rb')
    allLines = file.readlines()
    file.close()

    file = open('mi_refine_unsorted.mtz','wb')
    file.writelines(allLines)
    file.close()

    # Use CCP4/CAD to ensure that data is sorted properly for subsequent CCP4/FFT process

    fileexists = os.path.exists('mi_refine.mtz')
    if fileexists != 0:
        os.remove('mi_refine.mtz')

    file = open('mi_cad.inp','w')
    file.write('LABIN FILE_NUMBER 1 ALL\n')
    file.write('END\n')
    file.close()

    runcad = 'cad HKLIN1 mi_refine_unsorted.mtz HKLOUT mi_refine.mtz < mi_cad.inp > mi_cad.log'
    os.system(runcad)

    fileexists = os.path.exists('mi_refine.mtz')
    if fileexists != 0:
        os.remove('mi_refine_unsorted.mtz')
        os.remove('mi_cad.log')
        os.remove('mi_cad.inp')
    else:
        print 'The CAD run to resort the data seems to have failed'
        time.sleep(4)
        return 1

    # Extract MTZ labels for F and SD(F) and Rfree set

    file = open('mi_mtzdump.inp','w')
    file.write('HEADER\n')
    file.write('SYMMETRY\n')
    file.write('END\n')
    file.close()

    runmtz = 'mtzdump HKLIN mi_refine.mtz < mi_mtzdump.inp > mi_mtzdump.log'
    os.system(runmtz)

    file = open('mi_mtzdump.log','r')
    allLines = file.readlines()
    file.close()

    os.remove('mi_mtzdump.log')
    os.remove('mi_mtzdump.inp')

    read_columns = 'no'
    read_labels = 'no'
    read_resolution = 'no'
    read_cell = 'no'

    for eachLine in allLines:

        line_length = len(eachLine)

        if read_columns == 'yes' and line_length > 1:
            colList = eachLine.split()
            read_columns = 'no'

        if read_labels == 'yes' and line_length > 1:
            labelList = eachLine.split()
            read_labels = 'no'

        if read_resolution == 'yes' and line_length > 1:
            parseLine = eachLine.split()
            resolution_mtz = parseLine[5]
            read_resolution = 'no'

        if read_cell == 'yes' and line_length > 1:
            parseLine = eachLine.split()
            acell_mtz = parseLine[0]
            bcell_mtz = parseLine[1]
            ccell_mtz = parseLine[2]
            alpha_mtz = parseLine[3]
            beta_mtz = parseLine[4]
            gamma_mtz = parseLine[5]

            acell_mtz = float(acell_mtz)
            bcell_mtz = float(bcell_mtz)
            ccell_mtz = float(ccell_mtz)
            alpha_mtz = float(alpha_mtz)
            beta_mtz = float(beta_mtz)
            gamma_mtz = float(gamma_mtz)

            acell_mtz = round(acell_mtz,3)
            bcell_mtz = round(bcell_mtz,3)
            ccell_mtz = round(ccell_mtz,3)
            alpha_mtz = round(alpha_mtz,2)
            beta_mtz = round(beta_mtz,2)
            gamma_mtz = round(gamma_mtz,2)

            read_cell = 'no'

        if eachLine.find('* Number of Symmetry Operations') > -1:
            parseLine = eachLine.split('=')
            number_sym_ops = parseLine[1]
            number_sym_ops = int(number_sym_ops)

        if eachLine.find('* Column Labels :') > -1:
            read_columns = 'yes'

        if eachLine.find('* Column Types :') > -1:
            read_labels = 'yes'

        if eachLine.find('Resolution Range') > -1:
            read_resolution = 'yes'

        if eachLine.find('* Cell Dimensions :') > -1:
            read_cell = 'yes'

        if eachLine.find('* Space Group =') > -1:

            # SG number and name

            eachLine = eachLine.strip()
            parseLine = eachLine.split('=')
            space_group_out = parseLine[1]
            parseLine = space_group_out.split(quote)
            space_group = parseLine[0]
            space_group = space_group.strip()           
            space_group_mtz = parseLine[1]
            space_group_mtz = space_group_mtz.strip()

    list_length = len(labelList)
    count = 0

    while count < list_length:
        if labelList[count] == 'F' and flabel == 'none':
            flabel = colList[count]

        if labelList[count] == 'Q' and sigflabel == 'none':
            sigflabel = colList[count]

        if labelList[count] == 'I' and rfreelabel == 'none':
            rfreelabel = colList[count]

        if labelList[count] == 'D' and anomlabel == 'none':
            anomlabel = colList[count]

        if labelList[count] == 'Q' and siganomlabel == 'none':
            if anomlabel != 'none' and colList[count] != sigflabel:
                siganomlabel = colList[count]

        count = count + 1

    if flabel == 'none' or sigflabel == 'none' or rfreelabel == 'none':
        print 'MTZ labels for F,sd(F) or Rfree-data could not be established'
        time.sleep(4)
        return 1

    # Use CCP4/CAD to capture any anomalous difference data

    if anomlabel != 'none' and siganomlabel != 'none':
        
        fileexists = os.path.exists('mi_anommap.mtz')
        if fileexists != 0:
            os.remove('mi_anommap.mtz')

        file = open('mi_cad.inp','w')
        file.write('LABIN FILE_NUMBER 1 E1=')
        file.write(anomlabel)
        file.write(' E2=')
        file.write(siganomlabel)
        file.write('\n')

        # Avoid issues of the anomalous data exceeding the refinement resolution
        if max_res != 'none':
            file.write('RESOLUTION FILE_NUMBER 1 1000.0 ')
            file.write(max_res)
            file.write('/n')
                       
        file.write('END\n')
        file.close()

        runcad = 'cad HKLIN1 mi_refine.mtz HKLOUT mi_anommap.mtz < mi_cad.inp > mi_cad.log'
        os.system(runcad)

        fileexists = os.path.exists('mi_anommap.mtz')
        if fileexists != 0:
            os.remove('mi_cad.log')
            os.remove('mi_cad.inp')
        else:
            print 'The CAD run to extract anomalous difference data seems to have failed'
            time.sleep(4)
            return 1       

    # Set resolution

    if max_res == 'none':
        resolution_output = resolution_mtz
    else:
        resolution_output = max_res

    # Copy coordinates to local working area and collect information on content

    offset = 0
    res_count = 0

    cryst_found = 'no'

    file = open(pdbfile,'r')
    allLines = file.readlines()
    file.close()

    # Precheck to look for CRYST1 record

    for eachLine in allLines:

        tag = eachLine[0:6]
        tag = tag.strip()

        if tag == 'CRYST1':
            cryst_found = 'yes'

    file = open('mi_refine.pdb','w')

    # Add CRYST1 record if not present

    if cryst_found == 'no':

        acell_mtz = '%.3f'%(acell_mtz)                
        acell_mtz = str(acell_mtz)
        acell_mtz = acell_mtz.rjust(9)

        bcell_mtz = '%.3f'%(bcell_mtz)
        bcell_mtz = str(bcell_mtz)
        bcell_mtz = bcell_mtz.rjust(9)

        ccell_mtz = '%.3f'%(ccell_mtz)
        ccell_mtz = str(ccell_mtz)
        ccell_mtz = ccell_mtz.rjust(9)

        alpha_mtz = '%.2f'%(alpha_mtz)
        alpha_mtz = str(alpha_mtz)
        alpha_mtz = alpha_mtz.rjust(7)

        beta_mtz = '%.2f'%(beta_mtz)
        beta_mtz = str(beta_mtz)
        beta_mtz = beta_mtz.rjust(7)

        gamma_mtz = '%.2f'%(gamma_mtz)
        gamma_mtz = str(gamma_mtz)
        gamma_mtz = gamma_mtz.rjust(7)

        aLine = 'CRYST1' + acell_mtz + bcell_mtz + ccell_mtz \
                   + alpha_mtz + beta_mtz + gamma_mtz + ' ' + space_group_mtz

        file.write(aLine)
        file.write('\n')

    # Read/write file contents

    seq_chain_prev = '?'

    for eachLine in allLines:

        tag = eachLine[0:6]
        tag = tag.strip()

        # Parse SEQRES records into chain,name,dummy-number lists

        if tag == 'SEQRES':

            SEQRES = eachLine.strip()
            aList_SEQRES.append(SEQRES)

            seq_chain = eachLine[11:12]
            seq_chain = seq_chain.strip()
            seqLine = eachLine[19:70]
            parseLine = seqLine.split()

            if seq_chain != seq_chain_prev:
                aList_sequence_chain_id.append(seq_chain)

            seq_chain_prev = seq_chain

            length = len(parseLine)

            count = 0
            while count < length:

                res_name = parseLine[count]
            
                aList_sequence_chain.append(seq_chain)
                aList_sequence_number.append('?')
                aList_sequence_resname.append(res_name)

                count = count + 1

        # Ensure pdb/mtz cell compatibility and check CRYST1 integrity

        if tag == 'CRYST1':

            parseLine = eachLine.split()
            length = len(parseLine)

            if length > 7:
                acell_pdb = parseLine[1]
                bcell_pdb = parseLine[2]
                ccell_pdb = parseLine[3]
                alpha_pdb = parseLine[4]
                beta_pdb = parseLine[5]
                gamma_pdb = parseLine[6]

                acell_pdb = float(acell_pdb)
                bcell_pdb = float(bcell_pdb)
                ccell_pdb = float(ccell_pdb)
                alpha_pdb = float(alpha_pdb)
                beta_pdb = float(beta_pdb)
                gamma_pdb = float(gamma_pdb)

                a_dif = acell_pdb - acell_mtz
                b_dif = bcell_pdb - bcell_mtz
                c_dif = ccell_pdb - ccell_mtz
                alpha_dif = alpha_pdb - alpha_mtz
                beta_dif = beta_pdb - beta_mtz
                gamma_dif = gamma_pdb - gamma_mtz

                a_dif = abs(a_dif)
                b_dif = abs(b_dif)
                c_dif = abs(c_dif)
                alpha_dif = abs(alpha_dif)
                beta_dif = abs(beta_dif)
                gamma_dif = abs(gamma_dif)

                # Fix to mtz values if disagreeing and write a new CRYST1 record

                if a_dif > 0.1 or b_dif > 0.1 or c_dif > 0.1 \
                   or alpha_dif > 0.1 or beta_dif > 0.1 or gamma_dif > 0.1:

                    acell_mtz = '%.3f'%(acell_mtz)                
                    acell_mtz = str(acell_mtz)
                    acell_mtz = acell_mtz.rjust(9)

                    bcell_mtz = '%.3f'%(bcell_mtz)
                    bcell_mtz = str(bcell_mtz)
                    bcell_mtz = bcell_mtz.rjust(9)

                    ccell_mtz = '%.3f'%(ccell_mtz)
                    ccell_mtz = str(ccell_mtz)
                    ccell_mtz = ccell_mtz.rjust(9)

                    alpha_mtz = '%.2f'%(alpha_mtz)
                    alpha_mtz = str(alpha_mtz)
                    alpha_mtz = alpha_mtz.rjust(7)

                    beta_mtz = '%.2f'%(beta_mtz)
                    beta_mtz = str(beta_mtz)
                    beta_mtz = beta_mtz.rjust(7)

                    gamma_mtz = '%.2f'%(gamma_mtz)
                    gamma_mtz = str(gamma_mtz)
                    gamma_mtz = gamma_mtz.rjust(7)

                    eachLine = 'CRYST1' + acell_mtz + bcell_mtz + ccell_mtz\
                               + alpha_mtz + beta_mtz + gamma_mtz + ' ' + space_group_mtz

        # Check atom information

        if tag == 'ATOM' or tag == 'HETATM':

            chain_id = eachLine[21:22]
            chain_id = chain_id.strip()
            res_number = eachLine[22:26]
            res_number = res_number.strip()
            res_name = eachLine[17:20]
            res_name = res_name.strip()
            atom_name = eachLine[12:16]

            # Check and fix potential atom justification issues for common ions

            atom_justify = 'OK'

            if atom_name == ' NA ':
                atom_name = 'NA  '
                atom_justify = 'notOK'

            if atom_name == ' MG ':
                atom_name = 'MG  '
                atom_justify = 'notOK'

            if atom_name == ' CL ':
                atom_name = 'CL  '
                atom_justify = 'notOK'

            if atom_name == ' CR ':
                atom_name = 'CR  '
                atom_justify = 'notOK'

            if atom_name == ' MN ':
                atom_name = 'MN  '
                atom_justify = 'notOK'

            if atom_name == ' FE ':
                atom_name = 'FE  '
                atom_justify = 'notOK'

            if atom_name == ' CO ':
                atom_name = 'CO  '
                atom_justify = 'notOK'

            if atom_name == ' NI ':
                atom_name = 'NI  '
                atom_justify = 'notOK'

            if atom_name == ' CU ':
                atom_name = 'CU  '
                atom_justify = 'notOK'

            if atom_name == ' ZN ':
                atom_name = 'ZN  '
                atom_justify = 'notOK'

            if atom_name == ' SE ':
                atom_name = 'SE  '
                atom_justify = 'notOK'

            if atom_name == ' BR ':
                atom_name = 'BR  '
                atom_justify = 'notOK'

            if atom_name == ' CS ':
                atom_name = 'CS  '
                atom_justify = 'notOK'

            if atom_justify == 'notOK':
                eachLine_fix = eachLine[0:12] + atom_name + eachLine[16:80]
                eachLine = eachLine_fix

            # Obtain residue (CA) count

            if eachLine.find(' CA ') > -1:
                res_count = res_count + 1

            # Check for waters

            if res_name == 'HOH':            
                water_count = water_count + 1

            # Get highest residue number in the chain X we will assign for waters

            if chain_id == 'X':
                int_res_number = int(res_number)
                if int_res_number > res_number_X_high:
                    res_number_X_high = int_res_number

            # Obtain protein chain names and terminii

            count = 0
            found = 'no'

            if res_name != 'HOH':

                if chain_id == ' ':
                    missing_protein_chain = 'yes'

                while count < number_chain_list:
                    if chain_id == aList_chains[count]:
                        found = 'yes'

                    count = count + 1

                if found == 'no':
                    aList_nterm.append(res_number)
                    aList_chains.append(chain_id)
                    number_chain_list = len(aList_chains)

                    if number_chain_list > 1:
                        aList_cterm.append(res_number_prev)

                res_number_prev = res_number

        # Write record but eliminate SCALE records to avoid issues with changed cells and CISPEP
        # records since it it better if they are recomputed following refitting

        if tag != 'SCALE1' and tag != 'SCALE2' and tag != 'SCALE3' and tag != 'CISPEP':
            eachLine = eachLine.strip()
            file.write(eachLine)

        file.write('\n')

    file.close()

    aList_cterm.append(res_number_prev)        
    number_chain_list = len(aList_chains)

    if cryst_found == 'no':
        print 'There was no CRYST1 record in the coordinate file - stopping\n'
        time.sleep(4)
        return 1

    # Set water picking defaults depending on molecule size

    water_add = 0.25*res_count
    water_add = int(water_add)

    # Copy library file (if any) to temp area. REFMAC5 read requires a full path.

    if os.path.basename(libfile) != 'none':

        file = open(libfile,'r')
        allLines = file.readlines()
        file.close()

        # Check for extra protein-ligand covalent links

        for eachLine in allLines:
            if eachLine.find('data_link_list') > -1:
                extra_links = 'yes'

        temp_lib = os.path.join(ccp4.scr,'mi_templib.lib')

        file = open(temp_lib,'w')
        file.writelines(allLines)
        file.close()

    else:

        temp_lib = 'none'

    # Set general file names 

    fileexists = os.path.exists(projectlog)
    if fileexists != 0:

        file = open(projectlog,'r')
        allLines = file.readlines()
        file.close()

        for eachLine in allLines:
            if eachLine.find('Job ID') > -1 and eachLine.find('refine') > -1:
                aList = eachLine.split('_')
                runid = aList[1]
                runid_int = int(runid)            
                runid_int = runid_int + 1
                runid = str(runid_int)

    job_id = job_prefix + runid

    # Fix case where there are unrecorded refinement files - get highest serial number

    filename_pdb = job_id + '.pdb'
    filename_pdb_full = os.path.join(workingdir,filename_pdb)

    fileexists = os.path.exists(filename_pdb_full)
    if fileexists != 0:

        runid_prev = 0

        aList_dir = dircache.listdir(workingdir)
        number_files = len(aList_dir)

        count = 0
        while count < number_files:

            afile = aList_dir[count]

            if afile.find('refine_') > -1 and afile.find('.pdb') > -1:

                afile_tag = afile.replace('refine_','')
                afile_tag = afile_tag.replace('.pdb','')

                if afile_tag.isdigit() == 1:               
                    runid_int = int(afile_tag)

                    if runid_int > runid_prev:
                        runid_prev = runid_int
                        runid_int = runid_int + 1
                        runid = str(runid_int)

                        job_id = job_prefix + runid

            count = count + 1
    #

    filename_log = job_id + '.log'
    filename_pdb = job_id + '.pdb'
    filename_mtz = job_id + '.mtz'
    filename_tls = job_id + '.tls'
    errorfile = job_id + '_errors.txt'
    filename_log_full = os.path.join(workingdir,filename_log)
    filename_pdb_full = os.path.join(workingdir,filename_pdb)
    filename_mtz_full = os.path.join(workingdir,filename_mtz)
    filename_errors_full = os.path.join(workingdir,errorfile)
    filename_tls_full = os.path.join(workingdir,filename_tls)

    fileexists = os.path.exists(filename_log)
    if fileexists != 0:
        os.remove(filename_log)

    fileexists = os.path.exists(filename_pdb)
    if fileexists != 0:
        os.remove(filename_pdb)

    fileexists = os.path.exists(filename_mtz)
    if fileexists != 0:
        os.remove(filename_mtz)

    fileexists = os.path.exists(filename_tls)
    if fileexists != 0:
        os.remove(filename_tls)

    #################################################
    # Start rigid-body refinement (REFMAC5) section #
    #################################################

    if ref_engine == 'rigid':

        print '\nStarting REFMAC5 rigid-body refinement process'
        print 'CCP4 scratch space:',working_ccp4_scratch
        print 'Job-ID:',job_id
        print 'Using mtz data:',flabel,',',sigflabel,',',rfreelabel

        # REFMAC specific file names

        filename_in = job_id + '.inp'
        filename_refmac_temp = job_id + '.refmac'
        filename_refmac = job_id + '_cif.txt'
        filename_refmac_full = os.path.join(workingdir,filename_refmac)

        # Setup

        file = open(filename_in,'w')
        file.write('LABIN FP=')
        file.write(flabel)
        file.write(' SIGFP=')
        file.write(sigflabel)
        file.write(' FREE=')
        file.write(rfreelabel)
        file.write('\nLABOUT FC=FC PHIC=PHIC DELFWT=DELFWT PHDELWT=PHDELFWT FWT=FWT FOM=FOM\n')
        
        if max_res != 'none':
            file.write('RESOLUTION 100.0 ')
            file.write(max_res)
            file.write('\n')

        file.write('FREE ')
        file.write(freeflag)
        file.write('\nREFI TYPE RIGID\n')
        file.write('REFI RESI MLKF\n')
        file.write('REFI BREF OVERall METH CGMAT\n')
        file.write('SCAL TYPE BULK LSSC ANIS FIXBulk SCBUlk 0.78 BBULK 180.0\n')
        file.write('SOLVENT NO\n')
        file.write('MAKE_RESTRAINTS HYDR N\n')
        file.write('RIGIDbody NCYC 12\n')

        # Set group definitions by chain-id

        count = 0
        while count < number_chain_list:

            group_number = count + 1
            group_number = str(group_number)
            chain_id = aList_chains[count]
            nterm = aList_nterm[count]
            cterm = aList_cterm[count]

            file.write('RIGIDbody GROUP ')
            file.write(group_number)
            file.write(' FROM ')
            file.write(nterm)
            file.write(' ')
            file.write(chain_id)
            file.write(' TO ')
            file.write(cterm)
            file.write(' ')
            file.write(chain_id)
            file.write('\n')

            count = count + 1
        #

        file.write('MONI FEW\n')
        file.write('BINS 10\n')
        file.write('USECWD\n')
        file.write('PNAME noid\n')
        file.write('DNAME ')
        file.write(job_id)
        file.write('\n')
        file.write('END\n')
        file.close()

        # Run

        runrefine = 'refmac5 HKLIN mi_refine.mtz XYZIN mi_refine.pdb XYZOUT mi_refine_out.pdb HKLOUT '\
                    + filename_mtz + ' < ' + filename_in + ' > ' + filename_log

        os.system(runrefine)

        # Clean-up and rename

        os.remove(filename_in)

        fileexists = os.path.exists('mi_refine_out.pdb')
        if fileexists != 0:
            os.rename('mi_refine_out.pdb',filename_pdb)
            print 'Output PDB file:',filename_pdb
        else:
            'REFMAC5 rigid-body o/p coordinate file was not found'
            time.sleep(4)
            return 1

        fileexists = os.path.exists(filename_mtz)
        if fileexists != 0:
            print 'Output MTZ file:',filename_mtz
        else:
            'REFMAC rigid-body o/p phased data file was not found'
            time.sleep(4)
            return 1

        fileexists = os.path.exists(filename_refmac)
        if fileexists != 0:
            os.remove(filename_refmac)

        fileexists = os.path.exists('mi_refine.mtz')
        if fileexists != 0:
            os.remove('mi_refine.mtz')        

        fileexists = os.path.exists(filename_refmac_temp)
        if fileexists != 0:
            os.rename(filename_refmac_temp,filename_refmac)
            print 'Output CIF log file:',filename_refmac

            # Parse global summary

            file = open(filename_refmac,'r')
            allLines = file.readlines()
            file.close()

            for eachLine in allLines:

                if eachLine.find('_refine.ls_R_factor_R_work') > -1:
                    parseLine = eachLine.split()
                    rwork = parseLine[1]

                if eachLine.find('_refine.ls_R_factor_R_free') > -1:
                    parseLine = eachLine.split()
                    rfree = parseLine[1]

        else:
            print 'REFMAC o/p CIF log file was not found'
            return 1

        fileexists = os.path.exists(filename_log)
        if fileexists != 0:
            print 'Output REFMAC5 log:',filename_log
        else:
            print 'The REFMAC5 log file was not found'
            time.sleep(4)
            return 1

        print 'Rwork=',rwork,' Rfree=',rfree

    ###########################
    # Start REFMAC5 section   #
    ###########################

    if ref_engine == 'refmac5':

        print '\nStarting REFMAC5 process'
        print 'CCP4 scratch space:',working_ccp4_scratch
        print 'Job-ID:',job_id
        print 'Using mtz data:',flabel,',',sigflabel,',',rfreelabel

        if os.path.basename(libfile) != 'none':
            print 'Using library file:',libfile

        # REFMAC specific file names

        filename_in = job_id + '.inp'
        filename_refmac_temp = job_id + '.refmac'
        filename_refmac = job_id + '_cif.txt'
        filename_refmac_full = os.path.join(workingdir,filename_refmac)

        # Establish TLS file

        if os.path.basename(tlsfile) != 'none':

            file = open(tlsfile,'r')
            allLines = file.readlines()
            file.close()

            file = open('mi_temp.tls','w')
            file.writelines(allLines)
            file.close()

            tls_files = ' TLSIN mi_temp.tls TLSOUT ' + filename_tls + ' '

        else:

            tls_files = ' '

        # Setup

        file = open(filename_in,'w')
        file.write('LABIN FP=')
        file.write(flabel)
        file.write(' SIGFP=')
        file.write(sigflabel)
        file.write(' FREE=')
        file.write(rfreelabel)
        file.write('\nLABOUT FC=FC PHIC=PHIC DELFWT=DELFWT PHDELWT=PHDELFWT FWT=FWT FOM=FOM\n')

        file.write('FREE ')
        file.write(freeflag)
        file.write('\n')

        # Options/defaults

        if max_res != 'none':
            file.write('RESOLUTION 100.0 ')
            file.write(max_res)
            file.write('\n')

        if bref_type == 'anisotropic':
            file.write('REFI BREF ANISotropic METH CGMAT\n')
        else:
            file.write('REFI BREF ISOT METH CGMAT\n')

        # Standard setup

        file.write('SCAL TYPE SIMPLE LSSC ANIS\n')
        file.write('SOLVENT YES\n')           
        file.write('REFI TYPE RESTtrained\n')
        file.write('REFI RESI MLKF\n')

        # TLS option - set uniform B, establish TLS then refine residual B-factors

        if os.path.basename(tlsfile) != 'none':
            file.write('REFI TLSC 20\n')
            file.write('BFAC SET 30.0\n')

        file.write('WEIGH MATRIX ')
        file.write(weight)
        file.write('\n')

        if extra_links == 'yes':
            file.write('MAKE_RESTRAINTS LINK Y\n')

        file.write('MAKE_RESTRAINTS CISP Y\n')
        file.write('MAKE_RESTRAINTS SS Y\n')
        file.write('MAKE_RESTRAINTS HYDR N\n')
        file.write('BFAC 1 2.0 2.5 3.0 4.5\n')
        file.write('NCYC ')
        file.write(cycles)
        file.write('\n')

        # validation monitor

        file.write('MONI DIST 6.0\n')
        file.write('MONI ANGL 8.0\n')   
        file.write('MONI TORSION 10.0\n')
        file.write('MONI PLANE 10.0\n')
        file.write('MONI VANderwaals 4.25\n')
        file.write('MONI CHIRAL 8.0\n')
        file.write('MONI BFACTOR 4.0\n')
        file.write('BINS 20\n')

        file.write('USECWD\n')
        file.write('PNAME noid\n')
        file.write('DNAME ')
        file.write(job_id)
        file.write('\n')
        file.write('END\n')
        file.close()

        # Run process over number of big cycles

        count = 0
        while count < big_cycles:

            print 'Refining'

            if os.path.basename(libfile) == 'none':
                runrefine = 'refmac5 HKLIN mi_refine.mtz XYZIN mi_refine.pdb XYZOUT mi_refine_out.pdb HKLOUT '\
                            + filename_mtz + tls_files + ' < ' + filename_in + ' > ' + filename_log
            else:
                runrefine = 'refmac5 HKLIN mi_refine.mtz XYZIN mi_refine.pdb LIBIN ' + temp_lib + \
                            ' XYZOUT mi_refine_out.pdb HKLOUT ' + filename_mtz + tls_files + ' < ' + filename_in + ' > ' + filename_log

            os.system(runrefine)

            # Check run completed

            fileexists = os.path.exists(filename_mtz)
            if fileexists == 0:
                print 'REFMAC o/p phased data file was not found'
                time.sleep(4)
                return 1

            fileexists = os.path.exists('mi_refine_out.pdb')
            if fileexists == 0:
                print 'REFMAC o/p coordinate file was not found'
                time.sleep(4)
                return 1

            #
            # Apply water picking option
            #

            if water_pick == 'yes' and count < water_cycles:

                print 'Water picking'

                # 1FF map

                file = open('mi_fft.inp','w')
                file.write('LABIN F1=DELFWT PHI=PHDELFWT\n')
                file.write('END\n')
                file.close()

                runfft = 'fft HKLIN ' + filename_mtz + ' MAPOUT mi_1ff.map < mi_fft.inp > mi_fft.log'
                os.system(runfft)

                fileexists = os.path.exists('mi_1ff.map')
                if fileexists == 0:
                    print 'FFT for water picking failed'
                    time.sleep(4)
                    return 1
                else:
                    os.remove('mi_fft.inp')
                    os.remove('mi_fft.log')

                # Setup crystal araound the protein

                file = open('mi_mapmask.inp','w')
                file.write('BORDER 5.0\n')
                file.write('EXTEND XTAL\n')
                file.write('END\n')
                file.close()

                runmapmask = 'mapmask XYZIN mi_refine_out.pdb MAPIN mi_1ff.map MAPOUT mi_1ff_masked.map < mi_mapmask.inp > mi_mapmask.log'
                os.system(runmapmask)

                fileexists = os.path.exists('mi_1ff_masked.map')
                if fileexists == 0:
                    print 'MAPMASK for water picking failed'
                    time.sleep(4)
                    return 1
                else:
                    os.remove('mi_mapmask.inp')
                    os.remove('mi_mapmask.log')
                    os.remove('mi_1ff.map')

                # Water peak picking

                file = open('mi_peakmax.inp','w')
                file.write('THRESHOLD RMS 4.0\n')
                file.write('OUTPUT PDB\n')
                file.write('BFACTOR 30.0 1.0\n')
                file.write('RESIDUE HOH\n')
                file.write('ATNAME O\n')
                file.write('CHAIN X\n')
                file.write('NUMPEAKS 500\n')
                file.write('EXCLUDE EDGE\n')
                file.write('END\n')
                file.close()
                
                runpeakmax = 'peakmax MAPIN mi_1ff_masked.map XYZOUT mi_refine_peaks.pdb < mi_peakmax.inp > mi_peakmax_wat.log'

                os.system(runpeakmax)

                fileexists = os.path.exists('mi_refine_peaks.pdb')
                if fileexists == 0:
                    print 'PEAKMAX run failed'
                    time.sleep(4)
                    return 1                

                # Water peak reduction by symmetry and protein proximity 
                
                file = open('mi_watpeak.inp','w')
                file.write('DISTANCE 2.3 3.5\n')
                file.write('CHAIN X\n')
                file.write('SYMMETRY ')
                file.write(space_group)
                file.write('\nEND\n')
                file.close()
                
                runwatpeak = 'watpeak XYZIN mi_refine_out.pdb PEAKS mi_refine_peaks.pdb XYZOUT mi_refine_wat.pdb < mi_watpeak.inp > mi_watpeak.log'

                os.system(runwatpeak)

                fileexists = os.path.exists('mi_refine_wat.pdb')
                if fileexists == 0:
                    print 'WATPEAK run failed'
                    time.sleep(4)
                    return 1

                # Capture water atom records up to limit
                # Adjust for ascending numbering within the water chain

                file = open('mi_refine_wat.pdb','r')
                allLines = file.readlines()
                file.close()

                aList_waters = []
                water_pick_counter = 1

                for eachLine in allLines:

                    tag = eachLine[0:6]
                    tag = tag.strip()

                    if tag == 'ATOM' or tag == 'HETATM':

                        water_pick_counter = water_pick_counter + 1
                        if water_pick_counter < water_add:                            
                            res_number_X_high = res_number_X_high + 1
                            str_res_number = str(res_number_X_high)
                            str_res_number = str_res_number.rjust(4)
                            atom_water_record = eachLine[0:22] + str_res_number + eachLine[26:80]
                            atom_water_record = atom_water_record.strip()
                            aList_waters.append(atom_water_record)

                # Clean-up coordinate file debris from water picking

                fileexists = os.path.exists('mi_refine_peaks.pdb')
                if fileexists != 0:
                    os.remove('mi_refine_peaks.pdb')

                fileexists = os.path.exists('mi_refine_wat.pdb')
                if fileexists != 0:                    
                    os.remove('mi_refine_wat.pdb')
                

                # Rewrite current PDB file ready to append new waters

                file = open('mi_refine_out.pdb','r')
                allLines = file.readlines()
                file.close()

                os.remove('mi_refine_out.pdb')

                file = open('mi_refine_out.pdb','w')

                for eachLine in allLines:

                    tag = eachLine[0:6]
                    tag = tag.strip()

                    if tag != 'END' and tag != 'CONECT':
                        file.write(eachLine)

                # Add new waters

                number_water_list = len(aList_waters)

                count_rec = 0
                while count_rec < number_water_list:
                    aLine = aList_waters[count_rec]
                    file.write(aLine)
                    file.write('\n')

                    count_rec = count_rec + 1

                file.write('END\n')
                file.close()

                print 'Number of waters added:  ',number_water_list
                        
                # Clean-up

                os.remove('mi_1ff_masked.map')
                os.remove('mi_peakmax.inp')
                os.remove('mi_peakmax_wat.log')
                os.remove('mi_watpeak.inp')
                os.remove('mi_watpeak.log')               

            #        
            # end of water picking option
            #

            os.remove('mi_refine.pdb')
            os.rename('mi_refine_out.pdb','mi_refine.pdb')

            count = count + 1

        # Clean-up and rename

        os.remove(filename_in)

        print 'Output MTZ file:',filename_mtz

        fileexists = os.path.exists('mi_refine.pdb')
        if fileexists != 0:
            os.rename('mi_refine.pdb',filename_pdb)
            print 'Output PDB file:',filename_pdb
        else:
            'REFMAC5 o/p coordinate file was not found'
            time.sleep(4)
            return 1

        fileexists = os.path.exists(filename_refmac)
        if fileexists != 0:
            os.remove(filename_refmac)

        fileexists = os.path.exists('mi_refine.mtz')
        if fileexists != 0:
            os.remove('mi_refine.mtz')

        fileexists = os.path.exists('mi_temp.tls')
        if fileexists != 0:
            os.remove('mi_temp.tls')

        fileexists = os.path.exists(filename_refmac_temp)
        if fileexists != 0:
            os.rename(filename_refmac_temp,filename_refmac)
            print 'Output CIF log file:',filename_refmac

            # Parse global summary

            file = open(filename_refmac,'r')
            allLines = file.readlines()
            file.close()

            for eachLine in allLines:

                if eachLine.find('_refine.ls_R_factor_R_work') > -1:
                    parseLine = eachLine.split()
                    rwork = parseLine[1]

                if eachLine.find('_refine.ls_R_factor_R_free') > -1:
                    parseLine = eachLine.split()
                    rfree = parseLine[1]

                if eachLine.find('r_bond_refined_d') > -1:
                    parseLine = eachLine.split()
                    rmsd_bonds = parseLine[2]

        else:
            'REFMAC o/p CIF log file was not found'
            return 1

        fileexists = os.path.exists(filename_log)
        if fileexists != 0:
            print 'Output REFMAC5 log:',filename_log
        else:
            print 'The REFMAC5 log file was not found'
            time.sleep(4)
            return 1

        print 'Rwork=',rwork,' Rfree=',rfree,' RMSD(bonds)=',rmsd_bonds

    ##########################
    # End of REFMAC5 section #
    ##########################

    #######################
    # Start SHELX section #
    #######################

    if ref_engine == 'shelx':

        print '\nStarting SHELX refinement process'
        print 'CCP4 scratch space:',working_ccp4_scratch
        print 'Job-ID:',job_id
        print 'Using mtz data:',flabel,',',sigflabel

        if missing_protein_chain == 'yes':
            print 'Chain identifiers must be assigned for all protein atoms'
            time.sleep(4)
            return 1

        if os.path.basename(libfile) != 'none':
            print 'Using library file:',libfile

            filelib = open(libfile,'r')
            allLiblines = filelib.readlines()
            filelib.close()

        # Base .hkl and .ins root files name on pdb file name

        ins_file = job_id + '.ins'
        hkl_file = job_id + '.hkl'
        ins_file_full = os.path.join(workingdir,ins_file)
        hkl_file_full = os.path.join(workingdir,hkl_file)

        filename_lst = job_id + '.lst'
        filename_res = job_id + '.res'
        filename_fcf = job_id + '.fcf'
        filename_mtz = job_id + '.mtz'
        filename_lst_full = os.path.join(workingdir,filename_lst)
        filename_res_full = os.path.join(workingdir,filename_res)
        filename_fcf_full = os.path.join(workingdir,filename_fcf)
        filename_mtz_full = os.path.join(workingdir,filename_mtz)

        fileexists = os.path.exists(ins_file)
        if fileexists != 0:
            os.remove(ins_file)

        fileexists = os.path.exists(hkl_file)
        if fileexists != 0:
            os.remove(hkl_file)

        # Setup reflection file in SHELX format

        file = open('mi_mtz2various.inp','w')
        file.write('LABIN FP=')
        file.write(flabel)
        file.write(' SIGFP=')
        file.write(sigflabel)
        file.write(' FREE=')
        file.write(rfreelabel)
        file.write('\n')

        if max_res != 'none':
            file.write('RESOLUTION 1000.0 ')
            file.write(max_res)
            file.write('\n')

        file.write('OUTPUT SHELX\n')
        file.write('EXCLUDE FREER 0\n')
        file.write('END\n')
        file.close()

        runmtz = 'mtz2various HKLIN mi_refine.mtz HKLOUT mi_refine.hkl < mi_mtz2various.inp > mi_mtz2various.log'
        os.system(runmtz)

        fileexists = os.path.exists('mi_mtz2various.inp')
        if fileexists != 0:
            os.remove('mi_mtz2various.inp')

        fileexists = os.path.exists('mi_mtz2various.log')
        if fileexists != 0:
            os.remove('mi_mtz2various.log')

        fileexists = os.path.exists('mi_refine.hkl')
        if fileexists != 0:

            file = open('mi_refine.hkl','r')
            allLines = file.readlines()
            file.close()

            os.remove('mi_refine.hkl')

            file = open(hkl_file,'w')

            for eachLine in allLines:

                write_reflection = 'yes'

                if eachLine.find('TITLE') > -1 or eachLine.find('CELL') > -1 or eachLine.find('ZERR') > -1\
                   or eachLine.find('LATT') > -1 or eachLine.find('SYMM') > -1 or eachLine.find('HKLF') > -1:
                    write_reflection = 'no'

                if write_reflection == 'yes':
                    file.write(eachLine)

            file.close()

        else:

            print 'File format conversion for SHELXH seems to have failed'
            time.sleep(4)
            return 1

        # Setup input coordinates/restraints file (.ins) file with SHELXPRO

        print 'Running SHELXPRO'

        fileexists = os.path.exists('mi_shelxpro.inp')
        if fileexists != 0:
            os.remove('mi_shelxpro.inp')

        file = open('mi_shelxpro.inp','w')

        if test_platform.find('win') > -1:
            file.write('mi_shelxpro\n')

        file.write('I\n')
        file.write('\n')
        file.write(ins_file)
        file.write('\nmi_refine.pdb\n')
        file.write('Written by MIFit\n')
        file.write('\n')
        file.write('\n')
        file.write('\n')
        file.write('\n')
        file.write('C\n')

        # Chain offsets (+1000 etc) since SHELX does not support chains

        if number_chain_list > 0:

            count = 0
            while count < number_chain_list:
                file.write('\n')

                count = count + 1

            # List N-terminii

            final_chain = number_chain_list - 1

            count = 0
            while count < number_chain_list:
                nterm = aList_nterm[count]
                file.write(nterm)

                if count < final_chain:
                    file.write('=\n')
                else:
                    file.write('\n')

                count = count + 1

            # List C-terminii 

            count = 0
            while count < number_chain_list:
                cterm = aList_cterm[count]
                file.write(cterm)

                if count < final_chain:
                    file.write('=\n')
                else:
                    file.write('\n')

                count = count + 1
        file.write('\n')

        file.write('\n')
        file.write('N\n')
        file.write('3\n')
        file.write('\n')
        file.write('Q\n')
        file.write('\n')
        file.close()

        # Execute SHELXPRO to obtain mi_refine.ins

        runshelxpro = '"' + shelxpro + '"' + ' < mi_shelxpro.inp > mi_shelxpro.log'
        os.system(runshelxpro)

        os.remove('mi_shelxpro.inp')
        os.remove('mi_shelxpro.log')

        # Adjust run parameters, insert restraints and rename 

        fileexists = os.path.exists(ins_file)
        if fileexists != 0:

            file = open(ins_file,'r')
            allLines = file.readlines()
            file.close()

            os.remove(ins_file)

            # Read/write to adjust ins file

            file = open(ins_file,'w')

            for eachLine in allLines:

                write_flag = 'no'

                # Capture number of molecules for PDB write later

                if eachLine.find('ZERR') > -1:
                    parseLine = eachLine.split()
                    number_molecules = parseLine[1]

                if eachLine.find('WGHT') > -1:

                    # Insert any extra restraints

                    if os.path.basename(libfile) != 'none':
                        for eachLibline in allLiblines:
                            file.write(eachLibline)

                        file.write('\n')

                    # Insert weight

                    file.write('WGHT ')
                    file.write(weight)

                    if bref_type == 'anisotropic':
                        file.write('\nANIS\n')

                    write_flag = 'yes'

                # Number of refinement cycles

                if eachLine.find('CGLS') > -1:
                    file.write('CGLS ')
                    file.write(cycles)
                    file.write('\n')
                    write_flag = 'yes'

                if write_flag == 'no':
                    file.write(eachLine)

            file.close()

        else:

            print 'SHELXPRO run failed to generate .ins file'
            time.sleep(4)
            return 1

        fileexists = os.path.exists('mi_shelxpro.pro')
        if fileexists != 0:
            os.remove('mi_shelxpro.pro')

        fileexists = os.path.exists('mi_shelxpro.ps')
        if fileexists != 0:        
            os.remove('mi_shelxpro.ps')

        # Execute SHELXH refinement job

        print 'Running SHELXH'

        runshelx = '"' + shelxh + '"' + ' ' + job_id + ' > mi_shelxh.log'
        os.system(runshelx)

        fileexists = os.path.exists('mi_shelxh.log')
        if fileexists != 0:

            file = open('mi_shelxh.log','r')
            allLines = file.readlines()
            file.close()

            for eachLine in allLines:
                if eachLine.find('R1') > -1:
                    parseLine = eachLine.split()
                    rwork = parseLine[2]

            os.rename('mi_shelxh.log',filename_log)

        else:
            print 'SHELXH job did not run'
            return 1

        fileexists = os.path.exists(filename_lst)
        if fileexists == 0:
            print 'SHELXH lst file was not created'
            time.sleep(4)
            return 1

        fileexists = os.path.exists(filename_res)
        if fileexists == 0:
            print 'SHELXH res file was not created'
            time.sleep(4)
            return 1        

        fileexists = os.path.exists(filename_fcf)
        if fileexists == 0:
            print 'SHELXH fcf file was not created'
            time.sleep(4)
            return 1   

        print 'Rwork=',rwork

        # Append water information (note that waters are renumbered 1,2,3..)

        if water_count > 0:
            pr_water_count = str(water_count)
            aList_chains.append(' ')
            aList_nterm.append('1')
            aList_cterm.append(pr_water_count)

            number_chain_list = number_chain_list + 1

        # Back convert to PDB format with SHELXPRO

        print 'Running SHELXPRO'

        file = open('mi_shelxpro.inp','w')

        if test_platform.find('win') > -1:
            file.write('mi_shelxpro\n')

        file.write('G\n')
        file.write('\n')
        file.write('S\n')
        file.write(filename_res)
        file.write('\n')
        file.write('N\n')
        file.write('Y\n')
        file.write('Y\n')
        file.write('N\n')
        file.write('K\n')
        file.write(number_molecules)
        file.write('\n')
        file.write('\n')
        file.write(filename_pdb)
        file.write('\n')
        file.write('Written by a MIFit application\n')

        # Loop over chains to put back correct chain-number pairs

        if number_chain_list > 0:

            count = 0
            while count < number_chain_list:

                file.write('$\n')

                chain_id = aList_chains[count]
                nterm = aList_nterm[count]
                cterm = aList_cterm[count]
                nterm_current = int(nterm) + (count + 1) * 1000
                cterm_current = int(cterm) + (count + 1) * 1000
                nterm_current = str(nterm_current)
                cterm_current = str(cterm_current)

                file.write(chain_id)
                file.write('\n')
                file.write('\n')
                file.write(nterm_current)
                file.write(' ')
                file.write(cterm_current)
                file.write('\n')
                file.write(nterm)
                file.write('\n')

                count = count + 1

        file.write('\n')
        file.write('Q\n')
        file.write('\n')

        file.close()

        # Execute SHELXPRO 

        runshelxpro = '"' + shelxpro + '"' + ' < mi_shelxpro.inp > mi_shelxpro.log'
        os.system(runshelxpro)

        fileexists = os.path.exists(filename_pdb)
        if fileexists == 0:
            print 'SHELXPRO job failed to generate PDB file'
            time.sleep(4)
            return 1

        # Back convert the fcf file phased data information into mtz format for easy MIFit load

        print 'Converting output data to mtz format for MIFit input'
        print 'Note:FWT contains pre computed 2Fo-Fc map coefficients'
        print '     and DELFWT contains precomputed Fo-Fc map coefficients'

        file = open(filename_fcf,'r')
        allLines = file.readlines()
        file.close()

        file = open(hkl_file,'w')

        for eachLine in allLines:

            tag = eachLine[0:1]
            tag = tag.strip()
            
            if tag != ' ' and tag != '_' and tag != '#':
                parseLine = eachLine.split()

                num_args = len(parseLine)
                if num_args == 7:

                    h = parseLine[0]
                    k = parseLine[1]
                    l = parseLine[2]
                    fobs_sq = parseLine[3]
                    fcalc = parseLine[5]
                    phase = parseLine[6]
                    
                    h = int(h)
                    k = int(k)
                    l = int(l)
                    fobs_sq = float(fobs_sq)
                    fcalc = float(fcalc)
                    phase = float(phase)

                    fobs = math.sqrt(fobs_sq)
                    twofofc = 2.0*fobs - fcalc
                    twofofc = round(twofofc,3)                    
                    fofc = fobs - fcalc
                    fofc = round(fofc,3)

                    aLine = str(h) + ' ' + str(k) + ' ' + str(l) + \
                            ' ' + str(twofofc) + ' ' + str(fofc) + ' ' + str(phase) + ' ' + str(phase)

                    file.write(aLine)
                    file.write('\n')

        file.close()

        # Step 2, convert ascii to mtz

        aLine = str(acell_mtz) + ' ' + str(bcell_mtz) + ' ' + str(ccell_mtz)\
                + ' ' + str(alpha_mtz) + ' ' + str(beta_mtz) + ' ' + str(gamma_mtz)

        file = open('mi_f2mtz.inp','w')

        file.write('NAME PROJECT Shelx_map_coeffs CRYSTAL 1 DATASET 1\n')
        file.write('CELL ')
        file.write(aLine)
        file.write('\n')
        file.write('SYMMETRY ')
        file.write(space_group)
        file.write('\n')
        file.write('LABOUT H K L FWT DELFWT PHWT PHDELFWT\n')
        file.write('CTYPOUT H H H F F P P\n')
        file.write('END\n')

        file.close()

        runf2mtz = 'f2mtz HKLIN ' + hkl_file + ' HKLOUT ' + filename_mtz + ' < mi_f2mtz.inp > mi_f2mtz.log'
        os.system(runf2mtz)

        fileexists = os.path.exists(filename_mtz)
        if fileexists == 0:
            print 'F2MTZ failed to convert to output mtz file'
            time.sleep(4)
            return 1
        else:
            os.remove('mi_f2mtz.inp')
            os.remove('mi_f2mtz.log')

        # Clean-up various intermediate files

        fileexists = os.path.exists(ins_file_full)
        if fileexists != 0:
            os.remove(ins_file_full)

        fileexists = os.path.exists(hkl_file_full)
        if fileexists != 0:
            os.remove(hkl_file_full)           

        fileexists = os.path.exists(filename_fcf_full)
        if fileexists != 0:
            os.remove(filename_fcf_full)

        fileexists = os.path.exists(filename_res_full)
        if fileexists != 0:
            os.remove(filename_res_full)           

        fileexists = os.path.exists('mi_refine.pdb')
        if fileexists != 0:
            os.remove('mi_refine.pdb')

        fileexists = os.path.exists('mi_refine.mtz')
        if fileexists != 0:
            os.remove('mi_refine.mtz')           

        fileexists = os.path.exists('mi_shelxpro.pro')
        if fileexists != 0:
            os.remove('mi_shelxpro.pro')

        fileexists = os.path.exists('mi_shelxpro.ps')
        if fileexists != 0:        
            os.remove('mi_shelxpro.ps')

        fileexists = os.path.exists('shelxpro.pro')
        if fileexists != 0:
            os.remove('shelxpro.pro')

        fileexists = os.path.exists('shelxpro.ps')
        if fileexists != 0:        
            os.remove('shelxpro.ps')
            
        fileexists = os.path.exists('mi_shelxpro.inp')
        if fileexists != 0:  
            os.remove('mi_shelxpro.inp')

        fileexists = os.path.exists('mi_shelxpro.log')
        if fileexists != 0:          
            os.remove('mi_shelxpro.log')    

    ########################
    # End of SHELX section #
    ########################

    ########################
    # Structure validation #
    ########################

    if validate == 'yes' and ref_engine == 'refmac5':

        print 'Checking structure'

        # Get entity list of current model

        file = open(filename_pdb,'r')
        allLines = file.readlines()
        file.close()

        chain_id_prev = '?'
        res_number_prev = '?'

        for eachLine in allLines:

            tag = eachLine[0:6]
            tag = tag.strip()

            if tag == 'ATOM' or tag == 'HETATM':

                chain_id = eachLine[21:22]
                res_number = eachLine[22:26]
                res_number = res_number.strip()
                res_name = eachLine[17:20]
                res_name = res_name.strip()
                atom_name = eachLine[12:16]
                atom_name = atom_name.strip()
                
                disorder_id = eachLine[16:17]
                disorder_id = disorder_id.strip()
                
                if disorder_id != '':
                    aList_disorder_chain.append(chain_id)
                    aList_disorder_resno.append(res_number)
                    aList_disorder_resname.append(res_name)

                # Form all atom list

                aList_allatoms_chain.append(chain_id)
                aList_allatoms_res_number.append(res_number)
                aList_allatoms_res_name.append(res_name)
                aList_allatoms_atom_name.append(atom_name)

                # Form residue list

                if res_name != 'HOH':
                    
                    if chain_id != chain_id_prev or res_number != res_number_prev:
                        aList_chain_store.append(chain_id)
                        aList_res_number_store.append(res_number)
                        aList_res_name_store.append(res_name)        

                chain_id_prev = chain_id
                res_number_prev = res_number

            # Identify any non-PRO cis peptide links

            if tag == 'CISPEP':
                chain = eachLine[29:30]
                resnumber = eachLine[32:35]
                resname = eachLine[25:28]
                resnumber = resnumber.strip()
                resname = resname.strip()

                if resname != 'PRO':
                    aList_cis_chain.append(chain)
                    aList_cis_resno.append(resnumber)
                    aList_cis_resname.append(resname)

        #######################################
        # Parse refmac stereochemical scores  #
        #######################################

        file = open(filename_log,'r')
        allLines = file.readlines()
        file.close()

        for eachLine in allLines:

            # Parse section limits

            if eachLine.find('****') > -1 or eachLine.find('----') > -1:
                bond_list = 'no'
                angle_list = 'no'
                chiral_list = 'no'
                contact_list = 'no'
                chiral_list = 'no'

            # Start logging on finding final iteration number

            if eachLine.find('CGMAT cycle number') > -1 and eachLine.find(cycles) > -1:
                iteration_final = 'yes'

            if iteration_final == 'yes' and eachLine.find('Restraint type') > -1:
                read_error_log = 'yes'

            # get abnormal bond list

            if bond_list == 'yes' and read_error_log == 'yes':

                chain = eachLine[0:1]
                chain = chain.strip()
                if chain != '':
                    resnumber = eachLine[1:5]
                    resname = eachLine[6:9]
                    resnumber = resnumber.strip()
                    resname = resname.strip()
                    aList_bonds_chain.append(chain)
                    aList_bonds_resno.append(resnumber)
                    aList_bonds_resname.append(resname)

            if eachLine.find('Bond distance deviations ') > -1:
                bond_list = 'yes'

            # get abnormal bond angle list

            if angle_list == 'yes' and read_error_log == 'yes':

                chain = eachLine[0:1]
                chain = chain.strip()
                if chain != '':
                    resnumber = eachLine[1:5]
                    resname = eachLine[6:9]
                    resnumber = resnumber.strip()
                    resname = resname.strip()
                    aList_angles_chain.append(chain)
                    aList_angles_resno.append(resnumber)
                    aList_angles_resname.append(resname)

            if eachLine.find('Bond angle deviations ') > -1:
                angle_list = 'yes'

            # get abnormal contacts list

            if contact_list == 'yes' and read_error_log == 'yes':

                chain = eachLine[0:1]
                chain = chain.strip()
                if chain != '':
                    resnumber = eachLine[1:5]
                    resnumber = resnumber.strip()
                    resname = eachLine[6:9]
                    resname = resname.strip()
                    chain2 = eachLine[18:19]
                    chain2 = chain2.strip()
                    resnumber2 = eachLine[19:24]
                    resnumber2 = resnumber2.strip()
                    resname2 = eachLine[24:27]
                    resname2 = resname2.strip()
                    disorder1 = eachLine[14:15]
                    disorder2 = eachLine[32:33]
                    disorder1 = disorder1.strip()
                    disorder2 = disorder2.strip()

                    # Skip intra-residue interactions

                    if chain != chain2 or resnumber != resnumber2:
                        if disorder1 == '.' and disorder2 == '.':
                            aList_contacts_chain.append(chain)
                            aList_contacts_resno.append(resnumber)
                            aList_contacts_resname.append(resname)
                            aList_contacts_chain.append(chain2)
                            aList_contacts_resno.append(resnumber2)
                            aList_contacts_resname.append(resname2)

            if eachLine.find('VDW deviations ') > -1:
                contact_list = 'yes'

            # get severe chiral center violations

            if chiral_list == 'yes' and read_error_log == 'yes':

                chain = eachLine[0:1]
                chain = chain.strip()
                if chain != '':
                    resnumber = eachLine[1:5]
                    resname = eachLine[6:9]
                    resnumber = resnumber.strip()
                    resname = resname.strip()
                    aList_chiral_chain.append(chain)
                    aList_chiral_resno.append(resnumber)
                    aList_chiral_resname.append(resname)

            if eachLine.find('Chiral volume deviations') > -1:
                chiral_list = 'yes'

        ###########################################################
        # Run omega check and phi-psi check using Richardson data #
        ###########################################################

        # Read Richardson data 

        # General (non-GLY, non-PRO) data

        fileexists = os.path.exists(phipsi_gen_datafile)
        if fileexists == 0:        
            print 'WARNING - Unable to locate general phi-psi validation data'

        else:

            file = open(phipsi_gen_datafile,'r')
            allLines = file.readlines()
            file.close

            for eachLine in allLines:

                tag = eachLine[0:1]

                if tag != '#':
                    aLine = eachLine.split()
                    phi = aLine[0]
                    psi = aLine[1]
                    phipsi_prob = aLine[2]
                    phi = float(phi)
                    psi = float(psi)
                    phipsi_prob = float(phipsi_prob)

                    aList_phi_all.append(phi)
                    aList_psi_all.append(psi)
                    aList_phipsi_prob_all.append(phipsi_prob)

            number_phipsi_gen_table = len(aList_phi_all)

        # GLY data

        fileexists = os.path.exists(phipsi_gly_datafile)
        if fileexists == 0:        
            print 'WARNING - Unable to locate GLY phi-psi validation data'

        else:

            file = open(phipsi_gly_datafile,'r')
            allLines = file.readlines()
            file.close

            for eachLine in allLines:

                tag = eachLine[0:1]

                if tag != '#':
                    aLine = eachLine.split()
                    phi = aLine[0]
                    psi = aLine[1]
                    phipsi_prob = aLine[2]
                    phi = float(phi)
                    psi = float(psi)
                    phipsi_prob = float(phipsi_prob)

                    aList_phi_gly.append(phi)
                    aList_psi_gly.append(psi)
                    aList_phipsi_prob_gly.append(phipsi_prob)

            number_phipsi_gly_table = len(aList_phi_gly)

        # PRO data

        fileexists = os.path.exists(phipsi_pro_datafile)
        if fileexists == 0:        
            print 'WARNING - Unable to locate PRO phi-psi validation data'

        else:

            file = open(phipsi_pro_datafile,'r')
            allLines = file.readlines()
            file.close

            for eachLine in allLines:

                tag = eachLine[0:1]

                if tag != '#':
                    aLine = eachLine.split()
                    phi = aLine[0]
                    psi = aLine[1]
                    phipsi_prob = aLine[2]
                    phi = float(phi)
                    psi = float(psi)
                    phipsi_prob = float(phipsi_prob)

                    aList_phi_pro.append(phi)
                    aList_psi_pro.append(psi)
                    aList_phipsi_prob_pro.append(phipsi_prob)

            number_phipsi_pro_table = len(aList_phi_pro)

        # run SECSTR to compute phi,psi,omega

        file = open(filename_pdb)
        allLines = file.readlines()
        file.close()

        file = open('mi_secstr.new','w')
        for eachLine in allLines:

            tag = eachLine[0:6]
            tag = tag.strip()
            if tag == 'ATOM' or tag == 'HETATM':
                file.write(eachLine)

        file.close()

        file = open('mi_secstr.inp','w')
        file.write('mi_secstr.new\n')
        file.close()

        runsecstr = 'secstr < mi_secstr.inp > mi_secstr.log'
        os.system(runsecstr)

        fileexists = os.path.exists('mi_secstr.rin')
        if fileexists == 0:

            print 'Phi-psi calculation failed'
            time.sleep(4)
            return 1

        else:

            file = open('mi_secstr.rin','r')
            allLines = file.readlines()
            file.close()

            os.remove('mi_secstr.new')
            os.remove('mi_secstr.inp')
            os.remove('mi_secstr.log')
            os.remove('mi_secstr.rin')

            for eachLine in allLines:

                res_name = eachLine[4:7]
                res_number = eachLine[9:13]
                chain_id = eachLine[8:9]
                phi = eachLine[15:22]
                psi = eachLine[22:29]
                omega = eachLine[29:36]

                res_name = res_name.strip()
                res_number = res_number.strip()
                chain_id = chain_id.strip()
                phi = phi.strip()
                psi = psi.strip()
                omega = omega.strip()

                omega = float(omega)
                phi = float(phi)
                psi = float(psi)

                amino_acid_count = amino_acid_count + 1.0

                #######################################################
                # Search for outliers versus phi-psi probability data #
                #######################################################

                if phi < 180 and psi < 180: 
                    lookup = 'yes'
                else:
                    lookup = 'no'

                # Jump to useful region of table (minus safety margin)

                phi_point = phi + 180
                phi_point = phi_point * 9.0
                phi_point = math.floor(phi_point) - 19
                phi_point = int(phi_point)

                if phi_point > 0:
                    count = phi_point
                else:
                    count = 0  

                if res_name != 'GLY' and res_name != 'PRO':

                    while count < number_phipsi_gen_table and lookup == 'yes':

                        phi_table = aList_phi_all[count]
                        phi_diff = phi - phi_table
                        phi_diff = abs(phi_diff)

                        if phi_diff < 2.0:

                            psi_table = aList_psi_all[count]
                            psi_diff = psi - psi_table
                            psi_diff = abs(psi_diff)

                            if psi_diff < 2.0:
                                phipsi_prob = aList_phipsi_prob_all[count]

                                if phipsi_prob < phipsi_thresh_gen:
                                    aList_rama_chain.append(chain_id)
                                    aList_rama_resno.append(res_number)
                                    aList_rama_resname.append(res_name)

                                lookup = 'no'

                        count = count + 1

                if res_name == 'GLY':

                    while count < number_phipsi_gly_table and lookup == 'yes':

                        phi_table = aList_phi_gly[count]
                        phi_diff = phi - phi_table
                        phi_diff = abs(phi_diff)

                        if phi_diff < 2.0:

                            psi_table = aList_psi_gly[count]
                            psi_diff = psi - psi_table
                            psi_diff = abs(psi_diff)

                            if psi_diff < 2.0:
                                phipsi_prob = aList_phipsi_prob_gly[count]

                                if phipsi_prob < phipsi_thresh_gly:
                                    aList_rama_chain.append(chain_id)
                                    aList_rama_resno.append(res_number)
                                    aList_rama_resname.append(res_name)

                                lookup = 'no'

                        count = count + 1

                if res_name == 'PRO':

                    while count < number_phipsi_pro_table and lookup == 'yes':

                        phi_table = aList_phi_pro[count]
                        phi_diff = phi - phi_table
                        phi_diff = abs(phi_diff)

                        if phi_diff < 2.0:

                            psi_table = aList_psi_pro[count]
                            psi_diff = psi - psi_table
                            psi_diff = abs(psi_diff)

                            if psi_diff < 2.0:
                                phipsi_prob = aList_phipsi_prob_pro[count]

                                if phipsi_prob < phipsi_thresh_pro:
                                    aList_rama_chain.append(chain_id)
                                    aList_rama_resno.append(res_number)
                                    aList_rama_resname.append(res_name)

                                lookup = 'no'

                        count = count + 1

                #############################        
                # Search for omega outliers #
                #############################

                omega = float(omega)
                if omega < 180.0:
                    if omega < 0.0:
                        omega = -omega

                    omega_deviation = omega_peak - omega
                    omega_deviation = abs(omega_deviation)
                    if omega_deviation > omega_thresh:
                        aList_omega_chain.append(chain_id)
                        aList_omega_resno.append(res_name)
                        aList_omega_resname.append(res_number)

        ####################################
        # Run sidechain check with ROTAMER #
        ####################################

        file = open('mi_rotamer.inp','w')
        file.write('DELT 45\n')
        file.write('END\n')
        file.close()

        runrotamer = 'rotamer XYZIN ' + filename_pdb + ' < mi_rotamer.inp > mi_rotamer.log'
        os.system(runrotamer)

        fileexists = os.path.exists('mi_rotamer.log')
        if fileexists == 0:

            print 'Rotamer validation check failed'
            time.sleep(4)
            return 1

        else:

            # Parse for chi-1 deviations (greater than 45 degrees)

            file = open('mi_rotamer.log','r')
            allLines = file.readlines()
            file.close()

            for eachLine in allLines:

                if eachLine.find(')') > -1 and eachLine.find('(') > -1:

                    if eachLine[11:12] == '*':
                        chain_id = eachLine[0:1]
                        residue_id = eachLine[1:5]
                        residue_name = eachLine[6:9]                    
                        chain_id = chain_id.strip()
                        residue_id = residue_id.strip()
                        residue_name = residue_name.strip()

                        aList_rotamer_chain.append(chain_id)
                        aList_rotamer_resno.append(residue_id)
                        aList_rotamer_resname.append(residue_name)

            os.remove('mi_rotamer.log')
            os.remove('mi_rotamer.inp')

        #################################################
        # Locate difference density features on protein #
        #################################################

        # Calculate 1FF map

        file = open('mi_fft.inp','w')
        file.write('LABIN F1=DELFWT PHI=PHDELFWT\n')
        file.write('END\n')
        file.close()

        runfft = 'fft HKLIN ' + filename_mtz + ' MAPOUT mi_1ff.map < mi_fft.inp > mi_fft.log'
        os.system(runfft)

        fileexists = os.path.exists('mi_1ff.map')
        if fileexists == 0:
            print 'FFT for density test failed'
            time.sleep(4)
            return 1
        else:
            os.remove('mi_fft.inp')
            os.remove('mi_fft.log')

        # Build density around protein

        file = open('mi_mapmask.inp','w')
        file.write('EXTEND XTAL\n')
        file.write('BORDER 2.0\n')
        file.write('END\n')
        file.close()

        runmapmask = 'mapmask MAPIN mi_1ff.map XYZIN ' + filename_pdb + ' MAPOUT mi_1ff_masked.map < mi_mapmask.inp > mi_mapmask.log'
        os.system(runmapmask)

        fileexists = os.path.exists('mi_1ff_masked.map')
        if fileexists == 0:
            print 'MAPMASK for density test failed'
            return 1
        else:
            os.remove('mi_mapmask.inp')
            os.remove('mi_mapmask.log')
            os.remove('mi_1ff.map')

        # Peak/hole pick near protein

        file = open('mi_peakmax.inp','w')
        file.write('THRESHOLD RMS 4.0 NEGATIVE\n')
        file.write('END\n')
        file.close()

        runpeakmax = 'peakmax MAPIN mi_1ff_masked.map XYZOUT mi_peakmax.pdb < mi_peakmax.inp > mi_peakmax.log 2> mi_peakmax_err.log'
        os.system(runpeakmax)

        fileexists = os.path.exists('mi_peakmax_err.log')
        if fileexists != 0:
            os.remove('mi_peakmax_err.log')

        fileexists = os.path.exists('mi_peakmax.log')
        if fileexists == 0:
            print 'PEAKMAX for density test failed'
            time.sleep(4)
            return 1
        else:
            os.remove('mi_peakmax.inp')
            os.remove('mi_1ff_masked.map')
            os.remove('mi_peakmax.log')

        # Identify amino acids within 2.0A of any 4 sigma peaks/holes

        fileexists = os.path.exists('mi_peakmax.pdb')
        if fileexists != 0:

            file = open('mi_peakmax.pdb','r')
            allLines = file.readlines()
            file.close()

            os.remove('mi_peakmax.pdb')

            for eachLine in allLines:
                tag = eachLine[0:6]
                tag = tag.strip()

                if tag == 'ATOM' or tag == 'HETATM':
                    x = eachLine[30:38]
                    y = eachLine[38:46]
                    z = eachLine[46:54]
                    x = float(x)
                    y = float(y)
                    z = float(z)
                    aList_peak_x.append(x)
                    aList_peak_y.append(y)
                    aList_peak_z.append(z)

            number_peaks = len(aList_peak_x)

            file = open(filename_pdb,'r')
            allLines = file.readlines()
            file.close()

            count = 0
            while count < number_peaks:

                xp = aList_peak_x[count]
                yp = aList_peak_y[count]
                zp = aList_peak_z[count]

                for eachLine in allLines:
                    tag = eachLine[0:6]
                    tag = tag.strip()

                    if tag == 'ATOM' or tag == 'HETATM':
                        x = eachLine[30:38]
                        y = eachLine[38:46]
                        z = eachLine[46:54]
                        x = float(x)
                        y = float(y)
                        z = float(z)

                        dist = (xp - x) ** 2 + (yp - y)**2 + (zp - z)** 2

                        if dist < 4.0:

                            chain = eachLine[21:22]
                            res_number = eachLine[22:26]
                            res_number = res_number.strip()
                            res_name = eachLine[17:20]
                            res_name = res_name.strip()

                            aList_density_chain.append(chain)
                            aList_density_resno.append(res_number)
                            aList_density_resname.append(res_name)

                count = count + 1

        ###################################################################
        # Build tidy error lists by combining error types for each entity #
        ###################################################################

        entity_count = len(aList_chain_store)
        bond_count = len(aList_bonds_chain)
        angles_count = len(aList_angles_chain)
        chiral_count = len(aList_chiral_chain)
        contacts_count = len(aList_contacts_chain)
        cis_count = len(aList_cis_chain)
        rotamer_count = len(aList_rotamer_chain)
        omega_count = len(aList_omega_chain)
        rama_count = len(aList_rama_chain)
        density_count = len(aList_density_chain)

        count = 0
        while count < entity_count:

            geom_error = '.'
            contacts_error = '.'
            omega_error = '.'
            phipsi_error = '.'
            rotamer_error = '.'
            cis_error = '.'
            density_error = '.'
            error_flag = 'no'

            chain_store = aList_chain_store[count]
            res_number_store = aList_res_number_store[count]
            res_name_store = aList_res_name_store[count]

            # bonds

            count1 = 0
            while count1 < bond_count:

                chain = aList_bonds_chain[count1]
                res_number = aList_bonds_resno[count1]
                res_name = aList_bonds_resname[count1]

                if chain == chain_store and res_number == res_number_store:
                    geom_error = 'G'
                    error_flag = 'yes'

                count1 = count1 + 1

            # angles

            count1 = 0
            while count1 < angles_count:

                chain = aList_angles_chain[count1]
                res_number = aList_angles_resno[count1]
                res_name = aList_angles_resname[count1]

                if chain == chain_store and res_number == res_number_store:
                    geom_error = 'G'
                    error_flag = 'yes'

                count1 = count1 + 1

            # contacts

            count1 = 0
            while count1 < contacts_count:

                chain = aList_contacts_chain[count1]
                res_number = aList_contacts_resno[count1]
                res_name = aList_contacts_resname[count1]

                if chain == chain_store and res_number == res_number_store:
                    contacts_error = 'V'
                    error_flag = 'yes'

                count1 = count1 + 1

            # chiral

            count1 = 0
            while count1 < chiral_count:

                chain = aList_chiral_chain[count1]
                res_number = aList_chiral_resno[count1]
                res_name = aList_chiral_resname[count1]

                if chain == chain_store and res_number == res_number_store:
                    geom_error = 'G'
                    error_flag = 'yes'

                count1 = count1 + 1

            # cis peptide

            count1 = 0
            while count1 < cis_count:

                chain = aList_cis_chain[count1]
                res_number = aList_cis_resno[count1]
                res_name = aList_cis_resname[count1]

                if chain == chain_store and res_number == res_number_store:
                    cis_error = 'C'
                    error_flag = 'yes'

                count1 = count1 + 1        

            # rotamer

            count1 = 0
            while count1 < rotamer_count:

                chain = aList_rotamer_chain[count1]
                res_number = aList_rotamer_resno[count1]
                res_name = aList_rotamer_resname[count1]

                if chain == chain_store and res_number == res_number_store:
                    rotamer_error = 'R'
                    error_flag = 'yes'

                count1 = count1 + 1 

            # omega angles

            count1 = 0
            while count1 < omega_count:

                chain = aList_omega_chain[count1]
                res_number = aList_omega_resno[count1]
                res_name = aList_omega_resname[count1]

                if chain == chain_store and res_number == res_number_store:
                    omega_error = 'O'
                    error_flag = 'yes'

                count1 = count1 + 1           

            # phi-psi

            count1 = 0
            while count1 < rama_count:

                chain = aList_rama_chain[count1]
                res_number = aList_rama_resno[count1]
                res_name = aList_rama_resname[count1]

                if chain == chain_store and res_number == res_number_store:
                    phipsi_error = 'P'
                    error_flag = 'yes'

                count1 = count1 + 1

            # Density

            count1 = 0
            while count1 < density_count:

                chain = aList_density_chain[count1]
                res_number = aList_density_resno[count1]
                res_name = aList_density_resname[count1]

                if chain == chain_store and res_number == res_number_store:
                    density_error = 'D'
                    error_flag = 'yes'

                count1 = count1 + 1 

            # Write all error types  for this residue

            if error_flag == 'yes':

                # count phi-psi and sidechain errors

                if phipsi_error == 'P':
                    count_phipsi = count_phipsi + 1.0

                if rotamer_error == 'R':
                    count_rotamer = count_rotamer + 1.0

                # Tidy output

                res_number_field = len(res_number_store)
                print_res_number = res_number_store    
                if res_number_field == 1:
                    print_res_number = '   ' + res_number_store 
                if res_number_field == 2:
                    print_res_number = '  ' + res_number_store
                if res_number_field == 3:
                    print_res_number = ' ' + res_number_store

                aLine = ' ' + chain_store + ' ' + print_res_number + ' ' + res_name_store + ' ' + geom_error + ' ' \
                        + contacts_error + ' ' + omega_error + ' ' + phipsi_error + ' ' + cis_error + ' ' \
                        + rotamer_error + ' ' + density_error

                aList_errors.append(aLine)

            count = count + 1

        # Local error counts

        percent_phi_psi = 100.0 * count_phipsi / amino_acid_count
        percent_phi_psi = round(percent_phi_psi,2)
        percent_phi_psi = str(percent_phi_psi)

        percent_rotamer = 100.0 * count_rotamer / amino_acid_count
        percent_rotamer = round(percent_rotamer,2)
        percent_rotamer = str(percent_rotamer)

        # Write error list

        number_errors = len(aList_errors)
        percent_errors = 100.0 * float(number_errors)/entity_count
        percent_errors = round(percent_errors,1)
        percent_errors = str(percent_errors)

        print 'Output putative error list:',errorfile
        print 'Percentage of residues in error list:', percent_errors

        file = open(errorfile,'w')
        file.write('#\n')
        file.write('# Working directory: ')
        file.write(workingdir)
        file.write('\n# Coordinates: ')
        file.write(filename_pdb)
        file.write('\n# Data: ')
        file.write(filename_mtz)
        file.write('\n#\n')
        file.write('# Rwork: ')
        file.write(rwork)
        file.write('\n# Rfree: ')
        file.write(rfree)
        file.write('\n# Percentage of residues outside Richardson phi-psi core: ')
        file.write(percent_phi_psi)
        file.write('\n# Percentage of residues with abnormal rotamers: ')
        file.write(percent_rotamer)
        file.write('\n# Percentage of residues flagged: ')
        file.write(percent_errors)
        file.write('\n#\n')
        file.write('# Residue list codes for severe abnormality types:\n')
        file.write('# (G)eometry, (V)an der Waals, (O)mega, (P)hi-psi, (C)is peptide,\n')
        file.write('# (R)otamer chi-1, (D)ensity\n')
        file.write('#\n')

        count = 0
        while count < number_errors:
            aLine = aList_errors[count]
            file.write(aLine)
            file.write('\n')

            count = count + 1

        file.write('#\n')
        file.close()

        #####################################################################
        # Establish records diagnostics in PDB REMARK 465,470,500, format   #
        #####################################################################

        # Determination of missing amino acids from SEQRES records if there were any

        number_sequence = len(aList_sequence_resname)
        number_chains = len(aList_sequence_chain_id)

        print '\nNumber amino acids in SEQRES:',number_sequence,'over',number_chains,'chains\n'

        if number_sequence > 0:

            # Load data for a particular chain

            count_chains = 0
            while count_chains < number_chains:

                # initialize lists that will be used for this sequence/structure comparison

                aList_sequence_resname_temp = []
                aList_sequence_resnumber_temp = []
                aList_structure_resname_temp = []
                aList_structure_resnumber_temp = []
            
                working_seq_chain = aList_sequence_chain_id[count_chains]

                sequence_match = 'no'

                # Load sequence data for the current chain

                count1 = 0
                while count1 < number_sequence:

                    seq_chain = aList_sequence_chain[count1]

                    if working_seq_chain == seq_chain:

                        resname = aList_sequence_resname[count1]
 
                        aList_sequence_resname_temp.append(resname)
                        aList_sequence_resnumber_temp.append('?')

                    count1 = count1 + 1

                # Load structure data for the current chain

                count1 = 0
                while count1 < entity_count:

                    structure_chain = aList_chain_store[count1]

                    if working_seq_chain == structure_chain:

                        resname = aList_res_name_store[count1]
                        resnumber = aList_res_number_store[count1]
 
                        aList_structure_resname_temp.append(resname)
                        aList_structure_resnumber_temp.append(resnumber)

                    count1 = count1 + 1

                # Algorithm for establishing numbering in sequence-structure comparison
            
                # match to leading pentamer in this structure along sequence

                number_structure_temp = len(aList_structure_resname_temp)
                number_sequence_temp = len(aList_sequence_resname_temp)
                number_sequence_search = number_sequence_temp - 6

                test_structure_resname_1 = aList_structure_resname_temp[0]
                test_structure_resname_2 = aList_structure_resname_temp[1]
                test_structure_resname_3 = aList_structure_resname_temp[2]
                test_structure_resname_4 = aList_structure_resname_temp[3]
                test_structure_resname_5 = aList_structure_resname_temp[4]

                structure_resnumber_start = aList_structure_resnumber_temp[0]
                structure_resnumber_start = int(structure_resnumber_start)

                count1 = 0
                while count1 < number_sequence_search:

                    count2 = count1 + 1
                    count3 = count1 + 2
                    count4 = count1 + 3
                    count5 = count1 + 4

                    test_sequence_resname_1 = aList_sequence_resname_temp[count1]
                    test_sequence_resname_2 = aList_sequence_resname_temp[count2]
                    test_sequence_resname_3 = aList_sequence_resname_temp[count3]
                    test_sequence_resname_4 = aList_sequence_resname_temp[count4]
                    test_sequence_resname_5 = aList_sequence_resname_temp[count5]

                    if test_structure_resname_1 == test_sequence_resname_1:
                        if test_structure_resname_2 == test_sequence_resname_2:
                            if test_structure_resname_3 == test_sequence_resname_3:
                                if test_structure_resname_4 == test_sequence_resname_4:
                                    if test_structure_resname_5 == test_sequence_resname_5:
                                         
                                        sequence_resnumber_start = structure_resnumber_start - count1
                                        count1 = number_sequence_search
                                        sequence_match = 'yes'

                    count1 = count1 + 1

                # Now setup sequence numbering List

                if sequence_match == 'yes':

                    count1=0
                    while count1 < number_sequence_temp:

                        sequence_resnumber_put = sequence_resnumber_start + count1
                        sequence_resnumber_put = str(sequence_resnumber_put)
                        aList_sequence_resnumber_temp[count1] = sequence_resnumber_put

                        count1 = count1 + 1

                    # Now analyse and catch missing residues

                    count1=0
                    while count1 < number_sequence_temp:

                        sequence_resname_put = aList_sequence_resname_temp[count1]
                        sequence_resnumber_put = aList_sequence_resnumber_temp[count1]

                        count2=0
                        while count2 < number_structure_temp:

                            structure_resnumber_put = aList_structure_resnumber_temp[count2]

                            find_error = 'yes'
                            if sequence_resnumber_put == structure_resnumber_put:
                                find_error = 'no'
                                count2 = number_structure_temp

                            count2= count2 + 1

                        if find_error == 'yes':
                            out_line = 'REMARK 465   1 ' + sequence_resname_put + ' ' + working_seq_chain + ' ' + sequence_resnumber_put
                            aList_missing_residues.append(out_line)

                        count1 = count1 + 1

                # End of loop over chains

                count_chains = count_chains + 1

            # Write missing residues

            number_missing_residues = len(aList_missing_residues)

            if number_missing_residues > 0:

                pdb_annotate.append('REMARK 465')
                pdb_annotate.append('REMARK 465 MISSING RESIDUES')
                pdb_annotate.append('REMARK 465 THE FOLLOWING RESIDUES WERE NOT LOCATED IN THE')
                pdb_annotate.append('REMARK 465 EXPERIMENT. (M=MODEL NUMBER; RES=RESIDUE NAME; C=CHAIN')
                pdb_annotate.append('REMARK 465 IDENTIFIER; SSSEQ=SEQUENCE NUMBER; I=INSERTION CODE.)')
                pdb_annotate.append('REMARK 465')
                pdb_annotate.append('REMARK 465   M RES C SSSEQI')

                count1 = 0
                while count1 < number_missing_residues:

                    out_line = aList_missing_residues[count1]
                    pdb_annotate.append(out_line)
                
                    count1 = count1 + 1

        ################################################################
        # Obtain a list of missing atoms in each residue  (REMARK 470) #
        ################################################################

        pdb_annotate.append('REMARK 470')
        pdb_annotate.append('REMARK 470 MISSING ATOM')
        pdb_annotate.append('REMARK 470 THE FOLLOWING RESIDUES HAVE MISSING ATOMS (M=MODEL NUMBER;')
        pdb_annotate.append('REMARK 470 RES=RESIDUE NAME; C=CHAIN IDENTIFIER; SSEQ=SEQUENCE NUMBER;')
        pdb_annotate.append('REMARK 470 I=INSERTION CODE):')
        pdb_annotate.append('REMARK 470   M RES CSSEQI  ATOMS') 

        number_allatoms = len(aList_allatoms_chain)
          
        count1 = 0
        while count1 < entity_count:
            chain_store = aList_chain_store[count1]
            res_number_store = aList_res_number_store[count1]
            res_name_store = aList_res_name_store[count1]

            # Get all atoms for this residue

            aList_current_residue_atoms = []
            count = 0
            while count < number_allatoms:
                chain_all = aList_allatoms_chain[count]
                res_number_all = aList_allatoms_res_number[count]
                res_name_all = aList_allatoms_res_name[count]
                atom_name_all = aList_allatoms_atom_name[count]

                if chain_store == chain_all and res_number_store == res_number_all:
                    aList_current_residue_atoms.append(atom_name_all)

                count = count + 1

            # Process to find missing atoms in this residue

            number_residue_atoms = len(aList_current_residue_atoms)
            aList_atoms_expected = []

            if res_name_store == 'GLY':
                aList_atoms_expected = aList_GLY_atoms                   

            if res_name_store == 'ALA':
                aList_atoms_expected = aList_ALA_atoms

            if res_name_store == 'VAL':
                aList_atoms_expected = aList_VAL_atoms

            if res_name_store == 'ILE':
                aList_atoms_expected = aList_ILE_atoms

            if res_name_store == 'LEU':
                aList_atoms_expected = aList_LEU_atoms

            if res_name_store == 'PHE':
                aList_atoms_expected = aList_PHE_atoms

            if res_name_store == 'PRO':
                aList_atoms_expected = aList_PRO_atoms

            if res_name_store == 'MET':
                aList_atoms_expected = aList_MET_atoms

            if res_name_store == 'TRP':
                aList_atoms_expected = aList_TRP_atoms

            if res_name_store == 'CYS':
                aList_atoms_expected = aList_CYS_atoms

            if res_name_store == 'SER':
                aList_atoms_expected = aList_SER_atoms

            if res_name_store == 'THR':
                aList_atoms_expected = aList_THR_atoms

            if res_name_store == 'ASN':
                aList_atoms_expected = aList_ASN_atoms

            if res_name_store == 'GLN':
                aList_atoms_expected = aList_GLN_atoms

            if res_name_store == 'TYR':
                aList_atoms_expected = aList_TYR_atoms

            if res_name_store == 'HIS':
                aList_atoms_expected = aList_HIS_atoms

            if res_name_store == 'ASP':
                aList_atoms_expected = aList_ASP_atoms

            if res_name_store == 'GLU':
                aList_atoms_expected = aList_GLU_atoms

            if res_name_store == 'LYS':
                aList_atoms_expected = aList_LYS_atoms

            if res_name_store == 'ARG':
                aList_atoms_expected = aList_ARG_atoms                      

            number_atoms_expected = len(aList_atoms_expected)

            if number_atoms_expected > 0:

                # Check each expected atomname to see if it is found
            
                aList_atoms_expected_flag = []
            
                count2 = 0
                while count2 < number_atoms_expected:
                    atom_name_expected = aList_atoms_expected[count2]

                    found = 'no'
                    count_current_atoms = 0
                    while count_current_atoms < number_residue_atoms:
                        atom_name = aList_current_residue_atoms[count_current_atoms]

                        if atom_name_expected == atom_name:
                            found = 'yes'

                        count_current_atoms = count_current_atoms + 1

                    if found == 'yes':
                        aList_atoms_expected_flag.append('yes')
                    else:
                        aList_atoms_expected_flag.append('no')
                    
                    count2 = count2 + 1
        
                # Collect missing atoms for this residue into a list and create formatted REMARK 470

                write_flag = 'no'
                out_list = ''
                count2 = 0
                while count2 < number_atoms_expected:
                    found = aList_atoms_expected_flag[count2]

                    if found == 'no':
                        atom_name = aList_atoms_expected[count2]

                        number_chars = len(atom_name)
                        if number_chars == 1:
                            atom_name = atom_name + '  '
                        if number_chars == 2:
                            atom_name = atom_name + ' '

                        out_list = out_list + '  ' + atom_name
                        write_flag = 'yes'

                    count2 = count2 + 1
            
                if write_flag == 'yes':

                    str_res_number_store = str(res_number_store)
                    number_chars = len(str_res_number_store)
                    if number_chars == 1:
                        str_res_number_store = str_res_number_store + '   '
                    if number_chars == 2:
                        str_res_number_store = str_res_number_store + '  '
                    if number_chars == 3:
                        str_res_number_store = str_res_number_store + ' '                   

                
                    out_line = 'REMARK 470   1 ' + res_name_store + ' ' + chain_store + ' ' + str_res_number_store + ' ' + out_list
                    pdb_annotate.append(out_line)
                
            # End of loop over entities
            
            count1 = count1 + 1

        ##########################################################################
        # MI stereochemistry subtopic to identify discrete and errors REMARK 500 #
        ##########################################################################

        pdb_annotate.append('REMARK 500')
        pdb_annotate.append('REMARK 500 GEOMETRY AND STEREOCHEMISTRY')
        pdb_annotate.append('REMARK 500 SUBTOPIC: DISCRETE DISORDER')
        pdb_annotate.append('REMARK 500')
        pdb_annotate.append('REMARK 500 RESIDUES IN MULTIPLE CONFORMATIONS')
        pdb_annotate.append('REMARK 500 THE FOLLOWING RESIDUES WERE DESCRIBED BY MULTIPLE')
        pdb_annotate.append('REMARK 500 CONFORMATIONS.(M=MODEL NUMBER; RES=RESIDUE NAME; C=CHAIN')
        pdb_annotate.append('REMARK 500 IDENTIFIER; SSSEQ=SEQUENCE NUMBER; I=INSERTION CODE.)')
        pdb_annotate.append('REMARK 500')
        pdb_annotate.append('REMARK 500   M RES C SSSEQI')
        pdb_annotate.append('REMARK 500')

        # Reduce disorder atom list to residue list

        disorder_count = len(aList_disorder_chain)

        write_error = 'no'
        count = 0
        while count < entity_count:
            chain_store = aList_chain_store[count]
            res_number_store = aList_res_number_store[count]
        
            count1 = 0
            found_error = 'no'
            while count1 < disorder_count:                  
                chain = aList_disorder_chain[count1]
                res_number = aList_disorder_resno[count1]
                res_name = aList_disorder_resname[count1]

                if chain == chain_store and res_number == res_number_store and found_error == 'no':
                    out_line = 'REMARK 500   1 ' + res_name + ' ' + chain + ' ' + res_number
                    pdb_annotate.append(out_line)

                    found_error = 'yes'
                    write_error = 'yes'

                count1 = count1 + 1

            count = count + 1

        if write_error == 'no':
            pdb_annotate.append('REMARK 500   THERE WERE NO RESIDUES IN MULTIPLE CONFORMATIONS')        

        # Covalent bond lengths

        pdb_annotate.append('REMARK 500')
        pdb_annotate.append('REMARK 500 GEOMETRY AND STEREOCHEMISTRY')
        pdb_annotate.append('REMARK 500 SUBTOPIC: COVALENT BOND LENGTHS (MI ERROR LIST)')
        pdb_annotate.append('REMARK 500')
        pdb_annotate.append('REMARK 500 THE FOLLOWING RESIDUES CONTAINED ABNORMAL BOND LENGTHS')
        pdb_annotate.append('REMARK 500 (M=MODEL NUMBER; RES=RESIDUE NAME; C=CHAIN')
        pdb_annotate.append('REMARK 500 IDENTIFIER; SSSEQ=SEQUENCE NUMBER; I=INSERTION CODE.)')
        pdb_annotate.append('REMARK 500')
        pdb_annotate.append('REMARK 500   M RES C SSSEQI')
        pdb_annotate.append('REMARK 500')       

        write_error = 'no'
        count = 0
        while count < entity_count:
            chain_store = aList_chain_store[count]
            res_number_store = aList_res_number_store[count]
        
            count1 = 0
            found_error = 'no'
            while count1 < bond_count:                  
                chain = aList_bonds_chain[count1]
                res_number = aList_bonds_resno[count1]
                res_name = aList_bonds_resname[count1]

                if chain == chain_store and res_number == res_number_store and found_error == 'no':
                    out_line = 'REMARK 500   1 ' + res_name + ' ' + chain + ' ' + res_number
                    pdb_annotate.append(out_line)

                    found_error = 'yes'
                    write_error = 'yes'

                count1 = count1 + 1

            count = count + 1

        if write_error == 'no':
            pdb_annotate.append('REMARK 500   THERE WERE NO SEVERE ABNORMALITIES IN THIS CATEGORY')
        
        # Covalent bond angles

        pdb_annotate.append('REMARK 500')
        pdb_annotate.append('REMARK 500 GEOMETRY AND STEREOCHEMISTRY')
        pdb_annotate.append('REMARK 500 SUBTOPIC: COVALENT BOND ANGLES (MI ERROR LIST)')
        pdb_annotate.append('REMARK 500')
        pdb_annotate.append('REMARK 500 THE FOLLOWING RESIDUES CONTAINED ABNORMAL BOND ANGLES')
        pdb_annotate.append('REMARK 500 (M=MODEL NUMBER; RES=RESIDUE NAME; C=CHAIN')
        pdb_annotate.append('REMARK 500 IDENTIFIER; SSSEQ=SEQUENCE NUMBER; I=INSERTION CODE.)')
        pdb_annotate.append('REMARK 500')
        pdb_annotate.append('REMARK 500   M RES C SSSEQI')    
        pdb_annotate.append('REMARK 500')

        write_error = 'no'
        count = 0
        while count < entity_count:
            chain_store = aList_chain_store[count]
            res_number_store = aList_res_number_store[count]
        
            count1 = 0
            found_error = 'no'
            while count1 < angles_count:                  
                chain = aList_angles_chain[count1]
                res_number = aList_angles_resno[count1]
                res_name = aList_angles_resname[count1]

                if chain == chain_store and res_number == res_number_store and found_error == 'no':
                    out_line = 'REMARK 500   1 ' + res_name + ' ' + chain + ' ' + res_number
                    pdb_annotate.append(out_line)

                    found_error = 'yes'
                    write_error = 'yes'

                count1 = count1 + 1

            count = count + 1

        if write_error == 'no':
            pdb_annotate.append('REMARK 500   THERE WERE NO SEVERE ABNORMALITIES IN THIS CATEGORY')
        
        # Chiral centers

        pdb_annotate.append('REMARK 500')
        pdb_annotate.append('REMARK 500 GEOMETRY AND STEREOCHEMISTRY')
        pdb_annotate.append('REMARK 500 SUBTOPIC: CHIRAL CENTERS (MI ERROR LIST)')
        pdb_annotate.append('REMARK 500')
        pdb_annotate.append('REMARK 500 THE FOLLOWING RESIDUES CONTAINED ABNORMAL CHIRAL CENTERS')
        pdb_annotate.append('REMARK 500 (M=MODEL NUMBER; RES=RESIDUE NAME; C=CHAIN')
        pdb_annotate.append('REMARK 500 IDENTIFIER; SSSEQ=SEQUENCE NUMBER; I=INSERTION CODE.)')
        pdb_annotate.append('REMARK 500')
        pdb_annotate.append('REMARK 500   M RES C SSSEQI')
        pdb_annotate.append('REMARK 500')

        write_error = 'no'
        count = 0
        while count < entity_count:
            chain_store = aList_chain_store[count]
            res_number_store = aList_res_number_store[count]
        
            count1 = 0
            found_error = 'no'
            while count1 < chiral_count:                  
                chain = aList_chiral_chain[count1]
                res_number = aList_chiral_resno[count1]
                res_name = aList_chiral_resname[count1]

                if chain == chain_store and res_number == res_number_store and found_error == 'no':
                    out_line = 'REMARK 500   1 ' + res_name + ' ' + chain + ' ' + res_number
                    pdb_annotate.append(out_line)

                    found_error = 'yes'
                    write_error = 'yes'

                count1 = count1 + 1

            count = count + 1

        if write_error == 'no':
            pdb_annotate.append('REMARK 500   THERE WERE NO SEVERE ABNORMALITIES IN THIS CATEGORY')

        # Abnormal omega 

        pdb_annotate.append('REMARK 500')
        pdb_annotate.append('REMARK 500 GEOMETRY AND STEREOCHEMISTRY')
        pdb_annotate.append('REMARK 500 SUBTOPIC: NON-CIS, NON-TRANS (MI ERROR LIST)')
        pdb_annotate.append('REMARK 500')
        pdb_annotate.append('REMARK 500 THE FOLLOWING RESIDUES CONTAINED ABNORMAL OMEGA ANGLES')
        pdb_annotate.append('REMARK 500 (M=MODEL NUMBER; RES=RESIDUE NAME; C=CHAIN')
        pdb_annotate.append('REMARK 500 IDENTIFIER; SSSEQ=SEQUENCE NUMBER; I=INSERTION CODE.)')
        pdb_annotate.append('REMARK 500')
        pdb_annotate.append('REMARK 500   M RES C SSSEQI')
        pdb_annotate.append('REMARK 500')

        write_error = 'no'
        count = 0
        while count < entity_count:
            chain_store = aList_chain_store[count]
            res_number_store = aList_res_number_store[count]
        
            count1 = 0
            found_error = 'no'
            while count1 < omega_count:                  
                chain = aList_omega_chain[count1]
                res_number = aList_omega_resno[count1]
                res_name = aList_omega_resname[count1]

                if chain == chain_store and res_number == res_number_store and found_error == 'no':
                    out_line = 'REMARK 500   1 ' + res_name + ' ' + chain + ' ' + res_number
                    pdb_annotate.append(out_line)

                    found_error = 'yes'
                    write_error = 'yes'

                count1 = count1 + 1

            count = count + 1

        if write_error == 'no':
            pdb_annotate.append('REMARK 500   THERE WERE NO SEVERE ABNORMALITIES IN THIS CATEGORY')

        # Close contacts (note - all)

        pdb_annotate.append('REMARK 500')
        pdb_annotate.append('REMARK 500 GEOMETRY AND STEREOCHEMISTRY')
        pdb_annotate.append('REMARK 500 SUBTOPIC: CLOSE CONTACTS (MI ERROR LIST)')
        pdb_annotate.append('REMARK 500')
        pdb_annotate.append('REMARK 500 THE FOLLOWING RESIDUES CONTAINED ABNORMAL CONTACT DISTANCES')
        pdb_annotate.append('REMARK 500 (M=MODEL NUMBER; RES=RESIDUE NAME; C=CHAIN')
        pdb_annotate.append('REMARK 500 IDENTIFIER; SSSEQ=SEQUENCE NUMBER; I=INSERTION CODE.)')
        pdb_annotate.append('REMARK 500')
        pdb_annotate.append('REMARK 500   M RES C SSSEQI')
        pdb_annotate.append('REMARK 500')

        write_error = 'no'
        count = 0
        while count < entity_count:
            chain_store = aList_chain_store[count]
            res_number_store = aList_res_number_store[count]
        
            count1 = 0
            found_error = 'no'
            while count1 < contacts_count:                  
                chain = aList_contacts_chain[count1]
                res_number = aList_contacts_resno[count1]
                res_name = aList_contacts_resname[count1]

                if chain == chain_store and res_number == res_number_store and found_error == 'no':
                    out_line = 'REMARK 500   1 ' + res_name + ' ' + chain + ' ' + res_number
                    pdb_annotate.append(out_line)

                    found_error = 'yes'
                    write_error = 'yes'

                count1 = count1 + 1

            count = count + 1
        
        if write_error == 'no':
            pdb_annotate.append('REMARK 500   THERE WERE NO SEVERE ABNORMALITIES IN THIS CATEGORY')

        # PHI-PSI data
    
        pdb_annotate.append('REMARK 500')
        pdb_annotate.append('REMARK 500 GEOMETRY AND STEREOCHEMISTRY')
        pdb_annotate.append('REMARK 500 SUBTOPIC: TORSION ANGLES (MI ERROR LIST)')
        pdb_annotate.append('REMARK 500')
        pdb_annotate.append('REMARK 500 THE FOLLOWING RESIDUES CONTAINED ABNORMAL PHI-PSI ANGLES')
        pdb_annotate.append('REMARK 500 (M=MODEL NUMBER; RES=RESIDUE NAME; C=CHAIN')
        pdb_annotate.append('REMARK 500 IDENTIFIER; SSSEQ=SEQUENCE NUMBER; I=INSERTION CODE.)')
        pdb_annotate.append('REMARK 500')
        pdb_annotate.append('REMARK 500   M RES C SSSEQI')
        pdb_annotate.append('REMARK 500') 

        write_error = 'no'
        count = 0
        while count < entity_count:
            chain_store = aList_chain_store[count]
            res_number_store = aList_res_number_store[count]
        
            count1 = 0
            found_error = 'no'
            while count1 < rama_count:                  
                chain = aList_rama_chain[count1]
                res_number = aList_rama_resno[count1]
                res_name = aList_rama_resname[count1]

                if chain == chain_store and res_number == res_number_store and found_error == 'no':
                    out_line = 'REMARK 500   1 ' + res_name + ' ' + chain + ' ' + res_number
                    pdb_annotate.append(out_line)

                    found_error = 'yes'
                    write_error = 'yes'

                count1 = count1 + 1

            count = count + 1

        if write_error == 'no':
            pdb_annotate.append('REMARK 500   THERE WERE NO SEVERE ABNORMALITIES IN THIS CATEGORY')

        # invented MI remark to identify density issues

        pdb_annotate.append('REMARK 501')
        pdb_annotate.append('REMARK 501 OTHER VALIDATION')   
        pdb_annotate.append('REMARK 501 SUBTOPIC: ELECTRON DENSITY (MI ERROR LIST)')
        pdb_annotate.append('REMARK 501')    
        pdb_annotate.append('REMARK 501 THE FOLLOWING RESIDUES ARE NEAR DENSITY DIFFERENCE FEATURES')
        pdb_annotate.append('REMARK 501 (M=MODEL NUMBER; RES=RESIDUE NAME; C=CHAIN')
        pdb_annotate.append('REMARK 501 IDENTIFIER; SSSEQ=SEQUENCE NUMBER; I=INSERTION CODE.)')
        pdb_annotate.append('REMARK 501')
        pdb_annotate.append('REMARK 501   M RES C SSSEQI')
        pdb_annotate.append('REMARK 501')   

        write_error = 'no'    
        count = 0
        while count < entity_count:
            chain_store = aList_chain_store[count]
            res_number_store = aList_res_number_store[count]
        
            count1 = 0
            found_error = 'no'
            while count1 < density_count:                  
                chain = aList_density_chain[count1]
                res_number = aList_density_resno[count1]
                res_name = aList_density_resname[count1]

                if chain == chain_store and res_number == res_number_store and found_error == 'no':
                    out_line = 'REMARK 501   1 ' + res_name + ' ' + chain + ' ' + res_number
                    pdb_annotate.append(out_line)

                    found_error = 'yes'
                    write_error = 'yes'

                count1 = count1 + 1

            count = count + 1

        if write_error == 'no':
            pdb_annotate.append('REMARK 501   THERE WERE NO SEVERE ABNORMALITIES IN THIS CATEGORY')


        # invented MI remark to identify cis-pep

        pdb_annotate.append('REMARK 501')
        pdb_annotate.append('REMARK 501 OTHER VALIDATION')   
        pdb_annotate.append('REMARK 501 SUBTOPIC: CIS PEPTIDE (MI ERROR LIST)')
        pdb_annotate.append('REMARK 501')    
        pdb_annotate.append('REMARK 501 THE FOLLOWING RESIDUES HAVE CIS PEPTIDE BONDS')
        pdb_annotate.append('REMARK 501 (M=MODEL NUMBER; RES=RESIDUE NAME; C=CHAIN')
        pdb_annotate.append('REMARK 501 IDENTIFIER; SSSEQ=SEQUENCE NUMBER; I=INSERTION CODE.)')
        pdb_annotate.append('REMARK 501')
        pdb_annotate.append('REMARK 501   M RES C SSSEQI')
        pdb_annotate.append('REMARK 501')   

        write_error = 'no'    
        count = 0
        while count < entity_count:
            chain_store = aList_chain_store[count]
            res_number_store = aList_res_number_store[count]
        
            count1 = 0
            found_error = 'no'
            while count1 < cis_count:                  
                chain = aList_cis_chain[count1]
                res_number = aList_cis_resno[count1]
                res_name = aList_cis_resname[count1]

                if chain == chain_store and res_number == res_number_store and found_error == 'no':
                    out_line = 'REMARK 501   1 ' + res_name + ' ' + chain + ' ' + res_number
                    pdb_annotate.append(out_line)

                    found_error = 'yes'
                    write_error = 'yes'

                count1 = count1 + 1

            count = count + 1

        if write_error == 'no':
            pdb_annotate.append('REMARK 501   THERE WERE NO SEVERE ABNORMALITIES IN THIS CATEGORY')

        # invented MI remark to identify rotamer (chi-1)

        pdb_annotate.append('REMARK 501')
        pdb_annotate.append('REMARK 501 OTHER VALIDATION')   
        pdb_annotate.append('REMARK 501 SUBTOPIC: ROTAMER (MI ERROR LIST)')
        pdb_annotate.append('REMARK 501')    
        pdb_annotate.append('REMARK 501 THE FOLLOWING RESIDUES HAVE CHI-1 ANGLES WHICH DEVIATE MORE')
        pdb_annotate.append('REMARK 501 THAN 45 DEGREES FROM A KNOWN ROTAMER')
        pdb_annotate.append('REMARK 501 (M=MODEL NUMBER; RES=RESIDUE NAME; C=CHAIN')
        pdb_annotate.append('REMARK 501 IDENTIFIER; SSSEQ=SEQUENCE NUMBER; I=INSERTION CODE.)')
        pdb_annotate.append('REMARK 501')
        pdb_annotate.append('REMARK 501   M RES C SSSEQI')
        pdb_annotate.append('REMARK 501')   

        write_error = 'no'    
        count = 0
        while count < entity_count:
            chain_store = aList_chain_store[count]
            res_number_store = aList_res_number_store[count]
        
            count1 = 0
            found_error = 'no'
            while count1 < rotamer_count:                  
                chain = aList_rotamer_chain[count1]
                res_number = aList_rotamer_resno[count1]
                res_name = aList_rotamer_resname[count1]

                if chain == chain_store and res_number == res_number_store and found_error == 'no':
                    out_line = 'REMARK 501   1 ' + res_name + ' ' + chain + ' ' + res_number
                    pdb_annotate.append(out_line)

                    found_error = 'yes'
                    write_error = 'yes'

                count1 = count1 + 1

            count = count + 1

        if write_error == 'no':
            pdb_annotate.append('REMARK 501   THERE WERE NO SEVERE ABNORMALITIES IN THIS CATEGORY')

        ################################
        # Insert annotation into PDB   #
        ################################

        num_lines = len(pdb_annotate)
        num_SEQRES = len(aList_SEQRES)
    
        file = open(filename_pdb_full,'r')
        allLines = file.readlines()
        file.close()

        os.remove(filename_pdb_full)

        file = open(filename_pdb_full,'w')

        tag_prev = '?'

        for eachLine in allLines:
        
            tag = eachLine[0:6]
            tag = tag.strip()

            if tag != 'REMARK' and tag_prev == 'REMARK':
   
                count = 0
                while count < num_lines:

                    out_line = pdb_annotate[count]
                    file.write(out_line)
                    file.write('\n')

                    count = count + 1

                count = 0
                while count < num_SEQRES:

                    out_line = aList_SEQRES[count]
                    file.write(out_line)
                    file.write('\n')

                    count = count + 1

            file.write(eachLine)

            tag_prev = tag
 
        file.close()

    else:

        print 'Structure checking only enabled for REFMAC5 refinement'

        file = open(errorfile,'w')
        file.write('#\n')
        file.write('# Generation of error lists requires REFMAC5 refinement\n')
        file.write('#\n')
        file.close()        

    ############################################################
    # Standard file for nomalous difference map when available #
    ############################################################

    if anomlabel != 'none' and siganomlabel != 'none' and ref_engine == 'refmac5':

        fileexists = os.path.exists('anom_diffmap')
        if fileexists != 0:
            os.remove('anom_diffmap.map')

        fileexists = os.path.exists('mi_anommap_out.mtz')
        if fileexists != 0:
            os.remove('mi_anommap_out.mtz')    

        # Combine anomalous coefficients with refined phases back onto the refined mtz

        file = open('mi_cad.inp','w')
        file.write('LABIN FILE_NUMBER 1 ALL\n')
        file.write('LABIN FILE_NUMBER 2 ALL\n')                       
        file.write('END\n')
        file.close()

        runcad = 'cad HKLIN1 ' + filename_mtz + ' HKLIN2 mi_anommap.mtz HKLOUT mi_anommap_out.mtz < mi_cad.inp > mi_cad.log'
        os.system(runcad)

        fileexists = os.path.exists('mi_anommap_out.mtz')
        if fileexists != 0:
            os.remove('mi_cad.log')
            os.remove('mi_cad.inp')
            os.remove('mi_anommap.mtz')
        else:
            print 'The CAD run to reattach anomalous difference data seems to have failed'
            time.sleep(4)
            return 1       

        # Use special CCP4/FFT condition that rotates phases for anomalous difference maps

        file = open('mi_fft.inp','w')
        file.write('LABIN DANO=')
        file.write(anomlabel)
        file.write(' PHI=PHIC\n')
        file.write('END\n')
        file.close()

        runfft = 'fft HKLIN mi_anommap_out.mtz MAPOUT mi_1ff.map < mi_fft.inp 1> mi_fft.log 2>mi_fft_err.log'
        os.system(runfft)

        os.remove('mi_fft.inp')
        
        fileexists = os.path.exists('mi_fft.log')
        if fileexists != 0:
            os.remove('mi_fft.log')
            
        fileexists = os.path.exists('mi_fft_err.log')
        if fileexists != 0:
            os.remove('mi_fft_err.log')          

        # Note that FFT may fail if anom columns are present but unfilled so not a stop

        fileexists = os.path.exists('mi_1ff.map')
        if fileexists != 0:

            print 'Creating anomalous difference map file: anom_diffmap.map'

            # Build cell around the protein in the anomalous difference map with CCP4/MAPMASK

            file = open('mi_mapmask.inp','w')
            file.write('BORDER 5.0\n')
            file.write('EXTEND XTAL\n')
            file.write('END\n')
            file.close()

            runmapmask = 'mapmask XYZIN ' + filename_pdb + ' MAPIN mi_1ff.map MAPOUT anom_diffmap.map < mi_mapmask.inp > mi_mapmask.log'
            os.system(runmapmask)

            fileexists = os.path.exists('anom_diffmap.map')
            if fileexists == 0:
                print 'MAPMASK for anomalous difference map failed'
                time.sleep(4)
                return 1
            else:
                os.remove('mi_mapmask.inp')
                os.remove('mi_mapmask.log')
                os.remove('mi_1ff.map')

            # Rename mtz carrying anomalous data to standard refinement output

            os.remove(filename_mtz)
            os.rename('mi_anommap_out.mtz',filename_mtz)

            filename_anom_full = os.path.join(workingdir,'anom_diffmap.map')

    ######################
    # Append project log #
    ######################

    print 'Writing project log'

    runtime = time.ctime(time.time())

    file = open(projectlog,'a')
    file.seek(0,2)
    file.write('Job ID: ')
    file.write(job_id)
    file.write('\nDate: ')
    file.write(runtime)
    file.write('\nInput atoms: ')
    file.write(pdbfile)
    file.write('\nInput data: ')
    file.write(mtzfile)
    file.write('\nInput library: ')
    file.write(libfile)  

    if ref_engine == 'rigid':
        file.write('\nOutput atoms: ')
        file.write(filename_pdb_full)    
        file.write('\nOutput phased data: ')
        file.write(filename_mtz_full)
        file.write('\nOutput log: ')
        file.write(filename_log_full)
        file.write('\nOutput CIF log: ')
        file.write(filename_refmac_full)
        file.write('\nOptions: none\n')
        file.write('Summary: REFMAC5 rigid-body Rwork=')
        file.write(rwork)
        file.write(' Rfree=')
        file.write(rfree)
        file.write(' Resolution=')
        file.write(resolution_output)

    if ref_engine == 'refmac5':
        file.write('\nOutput atoms: ')
        file.write(filename_pdb_full)    
        file.write('\nOutput phased data: ')
        file.write(filename_mtz_full)
        file.write('\nOutput log: ')
        file.write(filename_log_full)
        file.write('\nOutput CIF log: ')
        file.write(filename_refmac_full)
        file.write('\nOutput error list: ')
        file.write(filename_errors_full)
        file.write('\nOutput anomalous difference map: ')
        file.write(filename_anom_full)

        if water_pick == 'yes':
            file.write('\nOptions: water-pick\n')
        else:
            file.write('\nOptions: none\n')            

        file.write('Summary: REFMAC5 Rwork=')
        file.write(rwork)
        file.write(' Rfree=')
        file.write(rfree)
        file.write(' RMSD(bonds)=')
        file.write(rmsd_bonds)
        file.write(' Resolution=')
        file.write(resolution_output)

    if ref_engine == 'shelx':
        file.write('\nOutput pdb file: ')        
        file.write(filename_pdb_full)
        file.write('\nOutput precomputed map data file: ')
        file.write(filename_mtz_full)    
        file.write('\nOutput log file: ')
        file.write(filename_log_full)
        file.write('\nOutput lst file: ')
        file.write(filename_lst_full)
        file.write('\nOptions: none\n')
        file.write('Summary: SHELXH Rwork=')
        file.write(rwork)
        file.write(' Resolution=')
        file.write(resolution_output)

    if ref_engine == 'refmac5' or ref_engine == 'shelx':
        file.write('\nParameters: Weight=')
        file.write(weight)
        file.write(' Cycles=')
        file.write(cycles)
        file.write(' Bfactor=')
        file.write(bref_type)
        file.write(' TLS_input_file=')
        file.write(tlsfile)

    file.write('\n---------------\n')
    file.close()

    #  Clean-up job-specific temporary CCP4_SCR space

    fileexists = os.path.exists(temp_lib)
    if fileexists != 0:
        os.remove(temp_lib)

    fileexists = os.path.exists(working_ccp4_scratch)
    if fileexists != 0:

        dir_list = os.listdir(working_ccp4_scratch)
        number_files = len(dir_list)

        count = 0
        while count < number_files:
            target_file = dir_list[count]
            target_file_full_path = os.path.join(working_ccp4_scratch,target_file)
            os.remove(target_file_full_path)

            count = count + 1

        os.rmdir(working_ccp4_scratch)

    time.sleep(4)

    #
    return 0

if __name__ == "__main__":
    sys.exit(Run())
