#################################################################
#                                                               #
# Prepares an mmCIF annotation file that may be used for        #
# report generation or PDB structure deposition                 #
#                                                               #
# Based on 'deposit3d.py'   Badger et al,                       #
#                           Acta Cryst F61 818-820 (2005)       #
#                                                               #
# Copyright: Molecular Images   2005                            #
#                                                               #
# This script is distributed under the same conditions as MIFit #
#                                                               #
#################################################################

import sys
import os
import string
import time
import math
import getopt
import ccp4check

def Usage():
    print "Usage: %s [options]" % sys.argv[0]
    print "Options are:"
    print "  -p,--pdbfile=FILE                the pdb file"
    print "  -m,--mtzfile=FILE                the mtz file"
    print "  -l,--libfile=FILE                the library file. Default: no file"
    print "  -d,--workdir=DIR                 The working directory"
    print "     --seqfile=FILE                The FASTA sequence file. Default: no file"
    print "     --templatefile=FILE           The annotation template file. Default no file"
    print "     --datalogfile=FILE            Data merging log file. Default: no file"
    print "     --cif_write=yes or no         Write the cif file. Default: no"
    print "     --text_write=yes or no        Write an annotation text file. Default: no"
    print "     --html_write=yes or no        Write an annotation HTML file. Default: no"
    print "     --hkl_write=yes or no         Write cif reflection file. Default: no"
    print "     --map_write=yes or no         Write maps of the unit cell and around the ligand(s). Default: no"
    print "     --map_border=NUM              Half-side dimension (A) for map around extent of ligand. Default: 6.0" 
    print "     --image1=FILE                 Image file 1. Default: no file"
    print "     --image2=FILE                 Image file 2. Default: no file"
    print "     --image3=FILE                 Image file 3. Default: no file"
    print "     --title=\"STRING\"              Title for the structure. Default: '?'"                
    print "     --rootname=name               Root for output files: Default: pdbdeposit"
    print "     --nfree_exclude=NUM           RFree data flag. Default: 0"
    print "  -h,--molimagehome=DIR            Path to MIFit."
    print "  -?,--help                        This help file"
    print ""
    print '\nNotes on input file formats: \n'
    print '1. Each protein in the coordinate file should be identified by a separate'
    print '   chain-id (usually A,B,C,...). Non-protein entities should be contained in'
    print '   chains. The file should contain a CRYST1 record.\n'
    print '2. The reflection file should contain only one column of type F (structure'
    print '   factor amplitude), one column of type Q (standard deviation on structure '
    print '   factor amplitude) and one colume of type I (flags for working and test'
    print '   data for validation). The CCP4 convention with Rfree data flagged by "0"'
    print '   is applied. This may be changed via script default parameter "nfree_exclude".'
    print '   All data are used in the calculation, without cutoff on sd(F). i.e. Rall and'
    print '   Robs are synonymous.\n'
    print '3. The FASTA sequence file may contain blank lines and a title line identified'
    print '   by a ">" symbol.\n'
    print '4. The parsing of SCALEPACK log files has only been lightly tested with'
    print '   HKL2000.\n'

def Run(argv=None):
    if argv is None:
        argv=sys.argv

    # Hardcoded backdoor path to scripts (only needed to find phi-psi data)
    mifit_root = 'none'

    #####################################
    # Annotation default initialization #
    #####################################

    # Set CCP4 default flag for cross-validation set '0','mask' bulk solvent correction and TRUNCATE defaults for reduction

    nfree_exclude = '0'
    ref_bulksolvent = 'mask'

    truncate_default_i = '-4.0'
    truncate_default_f = '0.0'

    # Null initializations

    audit_author_name = '?'
    audit_contact_author_name = '?'
    audit_contact_author_email = '?'
    audit_contact_author_address = '?'
    audit_contact_author_phone = '?'
    audit_contact_author_fax = '?'

    citation_title = '?'
    citation_journal_abbrev = '?'
    citation_journal_volume = '?'
    citation_page_first = '?'
    citation_page_last = '?'
    citation_year = '?'
    citation_author_name = '?'

    data_collection_temp_K = '?'
    data_collection_date = '?'
    wavelengths = '?'
    beamline = '?'
    detector_type = '?'
    detector_maker = '?'
    monochromator_type = '?'
    xray_method = '?'

    computing_data_collection = '?'
    computing_data_reduction = '?'
    computing_structure_solution = '?'
    computing_molecular_graphics = '?'
    computing_structure_refinement = '?'

    protein_name = '?'
    protein_ec_number = '?'

    structure_title = '?'
    structure_class = '?'
    structure_keywords = '?'

    biological_unit = '?'

    sequence_databasename = '?'
    sequence_databasecode = '?'

    source_common_name = '?'
    source_scientific_name = '?'
    source_gene_name = '?'
    source_host_common_name = '?'
    source_host_scientific_name = '?'

    acell = '?'
    bcell = '?'
    ccell = '?'
    alpha = '?'
    beta = '?'
    gamma = '?'
    spgno = '?'
    spgname = '?'

    matthews_coef = '?'
    solvent_percent = '?'

    exptl_crystal_grow_method = '?'
    exptl_crystal_grow_pH = '?'
    exptl_crystal_grow_temp = '?'
    exptl_crystal_grow_components = '?'

    data_num_unmerged = '?'
    data_num = '?'
    data_rlow = '?'
    data_rhigh = '?'
    data_percentobs = '?'
    data_redund = '?'
    data_rmerge = '?'
    data_ioversig = '?'
    datas_num = '?'
    datas_num_unmerged = '?'
    datas_rlow = '?'
    datas_rhigh = '?'
    datas_percentobs = '?'
    datas_redund = '?'
    datas_rmerge = '?'
    datas_ioversig = '?'

    ref_b11 = '?'
    ref_b12 = '?'
    ref_b13 = '?'
    ref_b23 = '?'
    ref_b22 = '?'
    ref_b33 = '?'
    ref_dlow = '?'
    ref_dhigh = '?'
    ref_bmean = '?'
    ref_numobs = '?'
    ref_numall = '?'
    ref_numfree = '?'
    ref_numwork = '?'
    ref_percent = '?'
    ref_rall = '?'
    ref_robs = '?'
    ref_rwork = '?'
    ref_rfree = '?'
    ref_dbond = '?'
    ref_dangle = '?'
    ref_dtorsion = '?'
    ref_dchiral = '?'
    ref_dplane = '?'
    ref_bmbond = '?'
    ref_bmangle = '?'
    ref_bsbond = '?'
    ref_bsangle = '?'
    ref_solvent_vdw_probe_radii = '?'
    ref_solvent_ion_probe_radii = '?'
    ref_solvent_shrinkage_radii = '?'
    ref_ksolv = '?'
    ref_bsolv = '?'
    ref_natom = '?'
    ref_nsolvent = '?'

    ref_number_phi_psi_errors = '?'

    #####################################################
    # Initialize operating defaults and structure data  #
    #####################################################

    quote = """'"""
    inputfile = 'mi_rundeposit3d.txt'

    workingdir = 'none'
    pdbfile = 'none'
    mtzfile = 'none'
    seqfile = 'none'
    templatefile = 'none'
    libfile = 'none'
    datalogfile = 'none'
    mergeprog = 'none'
    local_pdbfile = 'mi_temp.pdb'
    local_mtzfile = 'mi_temp.mtz'
    liganddir = 'none'
    write_cif = 'no'
    write_text = 'no'
    write_html = 'no'
    write_hkl = 'no'
    write_map = 'no'
    write_image1 = 'none'
    write_image2 = 'none'
    write_image3 = 'none'
    write_title = 'none'
    rootname = 'pdbdeposit'

    job_prefix = 'deposit3d_'
    runid = '1'
    projectlog = 'project_history.txt'

    reflineList = []
    dataList = []
    dataList_prev = []
    mtzList = []
    aList_chains = []
    seqList = []
    symList = []
    aList_connect = []
    readdata = 'no'
    read_mtzlabels = 'no'
    read_cell = 'no'
    acell_mtz = '?'
    bcell_mtz = '?'
    ccell_mtz = '?'
    alpha_mtz = '?'
    beta_mtz = '?'
    gamma_mtz = '?'
    found_c = 'no'
    found_n = 'no'
    cryst_flag = 'no'
    ref_anisoflag = 'no'
    water_flag = 'no'
    read_project_count = 0
    atom_count = 0
    water_count = 0
    residue_count = 0
    solvent_count = 0
    ligand_count = 0
    num_connect = 0
    map_border = 6.0
    famp = '?'
    sd = '?'
    freer = '?'
    footnote = ' . '
    spgname_standard = '?'

    # Het group lists

    aList_hets = []
    aList_hets_names = []
    aList_hets_nonPDB = []
    aList_hets_number = []
    aList_hets_asym = []

    # Standard non-ligand

    aList_protein = ['ALA','ARG','ASN','ASP','CYS','GLN','GLU','GLY','HIS','ILE','LEU',\
                     'LYS','MET','PHE','PRO','SER','THR','TRP','TYR','VAL','HOH','MSE','PTR','SEP','TPO']

    number_protein_list = len(aList_protein)

    #########################
    # Read input parameters #
    #########################
    number_of_args = len(argv)
    args = argv[1:]
    optlist, args = getopt.getopt(
        args,'p:m:l:d:h:?',
        ["pdbfile=","mtzfile=","libfile=","workdir=","seqfile=",
         "templatefile=","datalogfile=","cif_write=","text_write=",
         "html_write=","hkl_write=","map_write=","map_border=","image1=","image2=",
         "image3=","title=","rootname=","nfree_exclude=","molimagehome=","help"])
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
            
            if arg_value=='--mtzfile' or arg_value=='-m':
                mtzfile = param_value
            elif arg_value=='--workdir' or arg_value=='-d':
                workingdir = param_value
            elif arg_value=='--pdbfile' or arg_value=='-p':
                pdbfile = param_value
            elif arg_value=='--libfile' or arg_value=='-l':
                libfile = param_value
            elif arg_value=='--seqfile':
                seqfile= param_value
            elif arg_value=='--templatefile':
                templatefile = param_value
            elif arg_value=='--datalogfile':
                datalogfile = param_value
            elif arg_value=='--cif_write':
                write_cif = param_value
            elif arg_value=='--text_write':
                write_text = param_value
            elif arg_value=='--html_write':
                write_html = param_value
            elif arg_value=='--hkl_write':
                write_hkl = param_value
            elif arg_value=='--map_write':
                write_map = param_value
            elif arg_value=='--map_border':
                map_border = param_value
            elif arg_value=='--image1':
                write_image1 = param_value
            elif arg_value=='--image2':
                write_image2 = param_value
            elif arg_value=='--image3':
                write_image3 = param_value
            elif arg_value=='--title':
                write_title = param_value
                structure_title = write_title
            elif arg_value=='--rootname':
                rootname = param_value
            elif arg_value=='--nfree_exclude':
                nfree_exclude= param_value
            elif arg_value=='--molimagehome':
                mifit_root = param_value
        count=count + 1

    # Trap and reset null paths to auxillary files

    if os.path.basename(seqfile) == 'none':
        seqfile = 'none'

    if os.path.basename(libfile) == 'none':
        libfile = 'none'
        
    if os.path.basename(templatefile) == 'none':
        templatefile = 'none'

    if os.path.basename(datalogfile) == 'none':
        datalogfile = 'none'       

    if os.path.basename(write_image1) == 'none':
        write_image1 = 'none'

    if os.path.basename(write_image2) == 'none':
        write_image2 = 'none'

    if os.path.basename(write_image3) == 'none':
        write_image3 = 'none'       
        
    #

    ccp4,error = ccp4check.ccp4check()
    if not ccp4:
      print '\n' + error + '\n'
      time.sleep(4)
      return 1
    
    fileexists = os.path.exists(mtzfile)
    if fileexists == 0:
        print 'The refinement data file was not found ',mtzfile
        time.sleep(4)
        return 1

    fileexists = os.path.exists(workingdir)
    if fileexists == 0:
        print 'The working directory was not found ',workingdir
        time.sleep(4)
        return 1

    fileexists = os.path.exists(pdbfile)
    if fileexists == 0:
        print 'The coordinate file was not found ',pdbfile
        time.sleep(4)
        return 1

    fileexists = os.path.exists(seqfile)
    if fileexists == 0 and seqfile != 'none':
        print 'The FASTA sequence file was not found ',seqfile
        time.sleep(4)
        return 1

    fileexists = os.path.exists(libfile)
    if fileexists == 0 and libfile != 'none':
        print 'The REFMAC5 dictionary file was not found ',libfile
        time.sleep(4)
        return 1

    fileexists = os.path.exists(templatefile)
    if fileexists == 0 and templatefile != 'none':
        print 'The annotation template file was not found ',templatefile
        time.sleep(4)
        return 1

    fileexists = os.path.exists(datalogfile)
    if fileexists == 0 and datalogfile != 'none':
        print 'The data merging log file was not found ',datalogfile
        time.sleep(4)
        return 1

    fileexists = os.path.exists(write_image1)
    if fileexists == 0 and write_image1 != 'none':
        print 'Image file 1 was not found ',write_image1
        time.sleep(4)
        return 1

    fileexists = os.path.exists(write_image2)
    if fileexists == 0 and write_image2 != 'none':
        print 'Image file 2 was not found ',write_image2
        time.sleep(4)
        return 1

    fileexists = os.path.exists(write_image3)
    if fileexists == 0 and write_image3 != 'none':
        print 'Image file 3 was not found ',write_image3
        time.sleep(4)
        return 1

    write_cif.lower()
    write_hkl.lower()
    write_text.lower()
    write_html.lower()
    write_map.lower()

    rootname = os.path.basename(rootname)
    if rootname == 'none':
        rootname = 'pdbdeposit'

    map_border = float(map_border)

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

    # Obtain run information if any previous logged runs

    fileexists = os.path.exists(projectlog)
    if fileexists != 0:

        file = open(projectlog,'r')
        allLines = file.readlines()
        file.close()

        for eachLine in allLines:
            if eachLine.find('Job ID') > -1 and eachLine.find(job_prefix) > -1:
                aList = eachLine.split('_')
                runid = aList[1]
                runid_int = int(runid)            
                runid_int = runid_int + 1
                runid = str(runid_int)

    job_id = job_prefix + runid

    # Copy PDB files to local for robust CCP4 operation

    file = open(pdbfile,'r')
    allLines = file.readlines()
    file.close()

    file = open(local_pdbfile,'w')
    file.writelines(allLines)
    file.close()

    # Copy MTZ files to local for robust CCP4 operation

    file = open(mtzfile,'rb')
    allLines = file.readlines()
    file.close()

    file = open(local_mtzfile,'wb')
    file.writelines(allLines)
    file.close()

    # Copy library file to ccp4_scr to establish valid full path

    if libfile != 'none':

        liganddir = os.path.join(ccp4.scr,'mi_temp.lib')

        file = open(libfile,'r')
        allLines = file.readlines()
        file.close()

        file = open(liganddir,'w')
        file.writelines(allLines)
        file.close()

    # Banner

    print '\n___________________________________________________________________\n'
    print ' ** MI_Deposit3D **\n'
    print 'Using CCP4 version:',ccp4.version
    print '___________________________________________________________________\n'

    ######################################################
    # Inspect coordinate data to obtain entity lists etc #
    ######################################################

    file = open(pdbfile,'r')
    allLines = file.readlines()
    file.close()

    ###############################
    # Get the list of HET groups  #
    ###############################

    for eachLine in allLines:

        tag = eachLine[0:6]
        tag = tag.strip()
        res_name = eachLine[17:20]
        res_name = res_name.strip()

        if tag == 'ATOM' or tag == 'HETATM':

            atom_count = atom_count + 1

            # Check if entity is a HET

            protein_type = 'no'
            count = 0
            while count < number_protein_list:
                protein_residue = aList_protein[count]

                if res_name == protein_residue:
                    protein_type = 'yes'

                count = count + 1

            # Add to HET list unless we already have it

            if protein_type == 'no':
                repeat = 'no'
                count = 0
                count_hets = len(aList_hets)

                while count < count_hets:
                    if res_name == aList_hets[count]:
                        repeat = 'yes'

                    count = count + 1

                if repeat == 'no':
                    aList_hets.append(res_name)

    count_hets = len(aList_hets)

    ###################################################
    # Get the list of protein chains and check input  #
    ###################################################

    for eachLine in allLines:

        tag = eachLine[0:6]
        tag = tag.strip()

        res_name = eachLine[17:20]
        res_name = res_name.strip()
        chain_id = eachLine[21:22]
        chain_id = chain_id.strip()
        atom_name = eachLine[13:16]
        atom_name = atom_name.strip()

        # Check for cell dimensions

        if tag == 'CRYST1':
            cryst_flag = 'yes'

        # Check for anisotropic refinement

        if tag == 'ANISOU':
            ref_anisoflag = 'yes'

        # Chains

        if tag == 'ATOM' or tag == 'HETATM':

            protein = 'yes'

            # Check if chain-id contains a protein atom

            count = 0
            while count < count_hets:
                het_name = aList_hets[count]

                if res_name == het_name:
                    protein = 'no'

                count = count + 1

            if res_name == 'HOH':
                protein = 'no'

            # Collect the chain-id if it is new

            if protein == 'yes':

                # Trap protein without chain id

                length_chain_id = len(chain_id)
                if length_chain_id == 0:
                    print '\nThere is a protein ATOM/HETATM record(s) without a chain-ids.'
                    print 'Each protein in the coordinate file should be identified by a separate'
                    print 'chain-id (usually A,B,C,...).\n'
                    time.sleep(4)
                    return 1

                # Add chain

                repeat_chain = 'no'
                count = 0
                count_chains = len(aList_chains)

                while count < count_chains:
                    if chain_id == aList_chains[count]:
                        repeat_chain = 'yes'

                    count = count + 1

                if repeat_chain == 'no':
                    aList_chains.append(chain_id)

                # Count the number of amino acids for the text report

                if atom_name == 'CA':
                    residue_count = residue_count + 1

            else:

                # Count the number of water atoms and the number of solvent atoms

                if res_name == 'HOH':
                    water_flag = 'yes'
                    water_count = water_count + 1

                solvent_count = solvent_count + 1

    # Integrity checks

    number_chains = len(aList_chains)
    if number_chains == 0:
        print '\nEach protein molecule must be identified by chain-id (A,B,..)\n'
        time.sleep(4)
        return 1

    if cryst_flag == 'no':
        print '\nThe coordinate file must contain a CRYST1 record\n'
        time.sleep(4)
        return 1

    #################
    # Count ligands #
    #################

    chain_id_prev = '?'
    res_number_prev = '?'

    for eachLine in allLines:

        tag = eachLine[0:6]
        tag = tag.strip()

        tag = eachLine[0:6]
        tag = tag.strip()

        if tag == 'ATOM' or tag == 'HETATM':

            chain_id = eachLine[21:22]
            chain_id = chain_id.strip()
            res_name = eachLine[17:20]
            res_name = res_name.strip()
            res_number = eachLine[22:26]
            res_number = res_number.strip()

            if chain_id != chain_id_prev or res_number != res_number_prev:

                count = 0
                while count < count_hets:

                    if res_name == aList_hets[count]:
                        ligand_count = ligand_count + 1

                    count = count + 1

            chain_id_prev = chain_id
            res_number_prev = res_number

    #########################################
    # Add links across peptide chain breaks #
    #########################################

    chain_id_prev = '?'
    res_number_prev = '-9999'

    for eachLine in allLines:

        tag = eachLine[0:6]
        tag = tag.strip()

        x = eachLine[30:38]
        y = eachLine[38:46]
        z = eachLine[46:54]

        if tag == 'ATOM' or tag == 'HETATM':

            # Check for peptide links across non-consecutive residue nos

            if res_name != 'HOH':

                if atom_name == 'C':
                    store_chain_c = chain_id
                    store_number_c = res_number
                    store_name_c = atom_name
                    store_res_c = res_name
                    xc = float(x)
                    yc = float(y)
                    zc = float(z)
                    found_c = 'yes'

                if atom_name == 'N':
                    store_chain_n = chain_id
                    store_number_n = res_number
                    store_name_n = atom_name
                    store_res_n = res_name
                    xn = float(x)
                    yn = float(y)
                    zn = float(z)
                    found_n = 'yes'

                if found_c == 'yes' and found_n == 'yes':
                    store_number_c_int = int(store_number_c)
                    store_number_n_int = int(store_number_n)
                    seqno_diff = store_number_c_int - store_number_n_int
                    seqno_diff = abs(seqno_diff)

                    if seqno_diff > 1:
                        dx = xn - xc
                        dy = yn - yc
                        dz = zn - zc
                        dist = dx*dx + dy*dy + dz*dz

                        # Add special TRANS link record if it seems to be one

                        if dist < 3.0:
                            num_connect = num_connect + 1

                            link_record = 'LINK         C   ' + store_res_c + ' ' + store_chain_c + store_number_c \
                                          + '                 N   ' + store_res_n + ' ' + store_chain_n + store_number_n \
                                          + '                TRANS'
                            aList_connect.append(link_record)
                            found_c = 'no'
                            found_n = 'no'

        chain_id_prev = chain_id
        res_number_prev = res_number

    ############################
    # Obtain the entity lists  #
    ############################

    fileexists = os.path.exists(ccp4.entitylist)
    if fileexists != 0:
        file = open(ccp4.entitylist,'r')
        allLines = file.readlines()
        file.close()
    else:
        print '\nList of PDB entities was not found\n'
        time.sleep(4)
        return 1

    count = 0
    count_hets = len(aList_hets)

    while count < count_hets:

        found = 'no'
        pdbentity = aList_hets[count]
        pdbentity = pdbentity.strip()

        for eachLine in allLines:
            tag = eachLine[0:4]
            tag = tag.strip()

            if tag == 'code':
                entitycode = eachLine[5:8]
                entitycode = entitycode.strip()

            if tag == 'name':
                entityname = eachLine[5:80]
                entityname = entityname.strip()

                if pdbentity == entitycode:
                    found = 'yes'
                    aList_hets_names.append(entityname)

        # Obtain names for novel ligands from dictionary and store

        if found == 'no':

            novelligand = pdbentity

            fileexists = os.path.exists(libfile)
            if fileexists != 0:
                file = open(libfile,'r')
                allLines = file.readlines()
                file.close()

                for eachLine in allLines:
                    if eachLine.find(pdbentity) > -1 and eachLine.find('non-polymer') > -1:
                        aLine = eachLine.split(quote)
                        novelligand = aLine[1]

            aList_hets_names.append(novelligand)
            aList_hets_nonPDB.append(pdbentity)

        count = count + 1

    # List PDB name assignments from the REFMAC5 list

    if count_hets > 0:
        print '\nList of PDB HETNAM assignments'
        print '=============================='

        count = 0
        while count < count_hets:

            prentitycode = aList_hets[count]
            prentityname = aList_hets_names[count]

            print prentitycode,prentityname

            count = count + 1

    # Establish lists of ligand entity codes and pointers

    number_entities = len(aList_hets)
    entity_list_number = 1
    count = 0

    if number_entities > 0:

        while count < number_entities:        
            entity_list_number = entity_list_number + 1
            pr_entity_list_number = str(entity_list_number)
            aList_hets_number.append(pr_entity_list_number)        

            entity_code = aList_hets[count]
            entity_code_asym = entity_code + '_W '
            aList_hets_asym.append(entity_code_asym)

            count = count + 1

    if water_flag == 'yes':
        entity_list_number = entity_list_number + 1
        water_entity_list_number = str(entity_list_number)

    # Obtain atom counts

    if ref_natom == 0:
        print '\nNo atoms were found in this PDB file\n'
        time.sleep(4)
        return 1

    ref_natom = str(atom_count)
    ref_nsolvent = str(solvent_count)

    ###################################################
    # Rewrite clean PDB now everything seems in order #
    ###################################################

    file = open(pdbfile,'r')
    allLines = file.readlines()
    file.close()

    file = open('temp_use.pdb','w')

    # Write any peptide LINK records

    if num_connect > 0:

        count = 0
        while count < num_connect:
            link_record = aList_connect[count]
            file.write(link_record)
            file.write('\n')

            count = count + 1

    # Write CRYST1, ATOM/HETATM, TER, END records

    for eachLine in allLines:
        tag = eachLine[0:6]
        tag = tag.strip()

        if tag == 'CRYST1':
            file.write(eachLine)

        if tag == 'ATOM' or tag == 'HETATM':
            file.write(eachLine)        

        if tag == 'TER':
            file.write(eachLine)

        if tag == 'END':
            file.write(eachLine)

    file.close()

    ##############################################
    # Begin reporting main calculation processes #
    ##############################################

    print '\nProcess Summary'
    print '==============='

    #########################################################
    # Analyse mtz file for label and spacegroup information #
    #########################################################                

    # Execution of MTZDUMP

    filename = 'mi_mtzdump.inp'
    file = open(filename,'w')
    file.write('HEADER\n')
    file.write('END\n')
    file.close()

    runmtzdump = 'mtzdump HKLIN ' + local_mtzfile + ' < ' + filename + ' > mi_mtzdump.out'
    os.system(runmtzdump)

    file = open('mi_mtzdump.out', 'r')
    allLines = file.readlines()
    file.close()

    os.remove('mi_mtzdump.out')
    os.remove('mi_mtzdump.inp')

    read_columns = 'no'
    read_labels = 'no'
    read_resolution = 'no'
    read_spacegroup = 'no'

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

        if eachLine.find('* Column Labels :') > -1:
            read_columns = 'yes'

        if eachLine.find('* Column Types :') > -1:
            read_labels = 'yes'

        if eachLine.find('Resolution Range') > -1:
            read_resolution = 'yes'

        if eachLine.find('* Space group') > -1:
            parseLine = eachLine.split('number')
            space_group_out = parseLine[1]
            spgno = space_group_out.replace(')','')
            spgno = spgno.strip()

        # Determine cell

        if read_cell == 'yes':
            mtzList = eachLine.split()
            mtzList_length = len(mtzList)

            if mtzList_length == 6:
                acell_mtz = mtzList[0]
                bcell_mtz = mtzList[1]
                ccell_mtz = mtzList[2]
                alpha_mtz = mtzList[3]
                beta_mtz = mtzList[4]
                gamma_mtz = mtzList[5]

                read_cell = 'no'

        if eachLine.find('Cell Dimensions') > -1:
            read_cell = 'yes'

        # Try to parse wavelength (available from CCP4 5 mtzdumps)

        if read_project_count > 0:
            read_project_count = read_project_count + 1      

        if eachLine.find('Dataset ID, project/crystal/dataset names, cell dimensions, wavelength') > -1:
            read_project_count = read_project_count + 1

        if read_project_count == 7:
            mtzList = eachLine.split()
            mtzList_length = len(mtzList)
            read_project_count = 0

            if mtzList_length == 1:
                mtzwavelength = mtzList[0]

                if mtzwavelength != '0.00000' and wavelengths == '?':
                    wavelengths = mtzwavelength

    list_length = len(labelList)
    count = 0

    while count < list_length:
        if labelList[count] == 'F' and famp == '?':
            famp = colList[count]

        if labelList[count] == 'Q' and sd == '?':
            sd = colList[count]

        if labelList[count] == 'I' and freer == '?':
            freer = colList[count]

        count = count + 1

    # Check all items were found

    if famp == '?':
        print 'The MTZ file label for structure factor amplitude was not determined (type F)'
        time.sleep(4)
        return 1
    else:
        print 'Using structure factor amplitude:',famp

    if sd == '?':
        print '\nThe MTZ file label for structure factor standard deviation was not determined (type Q)'
        time.sleep(4)
        return 1
    else:
        print 'Using standard deviation on structure factor amplitude:',sd

    if freer == '?':
        print '\nWarning: The MTZ file label for the freeR flag was not determined (type I)\n'
        time.sleep(4)
        #return 1
    else:
        print 'Using Free-R data set defined by flag:',freer

    # Set H-M space group name from space group number

    file = open(ccp4.symmetrylib, 'r')
    allLines = file.readlines()
    file.close()

    for eachLine in allLines:

        symList = eachLine.split()
        length_symList = len(symList)

        if length_symList > 1:
            if symList[0] == spgno:
                symList = eachLine.split(quote)
                spgname = symList[1]

    ################################
    # Parse FASTA sequence file    #
    ################################

    if seqfile != 'none':

        print 'Parsing FASTA file'

        # Load single character series into a large list

        title_count = 0

        file = open(seqfile,'r')
        allLines = file.readlines()
        file.close()

        for eachLine in allLines:

            parse_line = 'yes'

            # skip title line

            tag = eachLine[0:1]
            if tag == '>':
                parse_line = 'no'
                title_count = title_count + 1

            # Trap out multiple sequences deliminated by title

            if title_count > 1:
                parse_line = 'no'

            # skip blank lines

            line_length = len(eachLine)
            if line_length == 0:
                parse_line = 'no'

            # load

            if parse_line == 'yes':
                sequence_line = eachLine.strip()

                line_length = len(sequence_line)
                line_length = line_length + 1

                count = 1
                while count < line_length:
                    j = count 
                    i = count - 1
                    aacode = sequence_line[i:j]
                    seqList.append(aacode)

                    count = count + 1

        # Check that we have the sequence

        sequence_length = len(seqList)
        if sequence_length == 0:
            print '\nFASTA sequence extraction failed\n'
            time.sleep(4)
            return 1

    #####################################################
    # Determine merging log type for subsequent parsing #
    #####################################################

    if datalogfile != 'none':

        file = open(datalogfile,'r')
        allLines = file.readlines()
        file.close()

        for eachLine in allLines:
            if eachLine.find('CCP4') > -1:
                mergeprog = 'scala'

            if eachLine.find('Scalepack') > -1:
                mergeprog = 'scalepack'

            if eachLine.find('reading from a file:') > -1 and eachLine.find('.x') > -1:
                mergeprog = 'scalepack'

            if eachLine.find('d*TREK') > -1:
                mergeprog = 'dstartrek'

        if mergeprog != 'scala' and mergeprog != 'scalepack' and mergeprog != 'dstartrek':
            print '\nThe type of data merging log file was not identified'
            print 'Omit the data log file to complete the Deposit3D run\n'
            time.sleep(4)
            return 1

    if mergeprog == 'scala':
        computing_data_reduction = 'SCALA, TRUNCATE'

    if mergeprog == 'scalepack':
        computing_data_reduction = 'SCALEPACK, TRUNCATE'

    if mergeprog == 'dstartrek':
        computing_data_reduction = 'DSTARTREK, TRUNCATE'    

    ################ ########################
    # Parse D*TREK log file (if available)  #
    ############### #########################

    # This function tested on version 9.2

    if mergeprog == 'dstartrek' and datalogfile != 'none':

        print 'Parsing D*TREK log file'

        file = open(datalogfile,'r')
        allLines = file.readlines()
        file.close()

        for eachLine in allLines:

            if eachLine.find('Summary of data collection statistics'):
                readdata = 'yes'

            if readdata == 'yes':

                eachLine = eachLine.replace('(',' ')
                eachLine = eachLine.replace(')',' ')

                dataList = eachLine.split()

                if eachLine.find('Resolution range') > -1:
                    data_rlow = dataList[2]
                    datas_rlow = dataList[5]                 
                    data_rhigh = dataList[4]
                    datas_rhigh = dataList[7]

                if eachLine.find('Rmerge') > -1:
                    data_rmerge = dataList[1]
                    datas_rmerge = dataList[2]

                if eachLine.find('Total number of reflections') > -1:
                    data_num_unmerged = dataList[4]
                    datas_num_unmerged = '?'

                if eachLine.find('Number of unique reflections') > -1:
                    data_num = dataList[4]
                    datas_num = '?'

                if eachLine.find('Output <I/sigI>') > -1:
                    data_ioversig = dataList[2]
                    datas_ioversig = dataList[3]

                if eachLine.find('% completeness') > -1:
                    data_percentobs = dataList[2]
                    datas_percentobs = dataList[3]

                if eachLine.find('Average redundancy') > -1:
                    data_redund = dataList[2]
                    datas_redund = dataList[3]

    ########################################
    # Parse SCALA log file (if available)  #
    ########################################

    if mergeprog == 'scala' and datalogfile != 'none':

        print 'Parsing SCALA log file'

        file = open(datalogfile,'r')
        allLines = file.readlines()
        file.close()

        for eachLine in allLines:

            dataList = eachLine.split()
            num_items = len(dataList)

            if num_items == 3:
                
                if dataList[0] == 'Overall' and dataList[1] == 'InnerShell' and dataList[2] == 'OuterShell':
                    readdata = 'yes'
                
            if readdata == 'yes':

                if eachLine.find('Low resolution limit') > -1 and num_items == 6:
                    data_rlow = dataList[3]
                    datas_rlow = dataList[5]                 

                if eachLine.find('High resolution limit') > -1 and num_items == 6:
                    data_rhigh = dataList[3]
                    datas_rhigh = dataList[5]

                if eachLine.find('Rmerge') > -1 and num_items == 4:
                    data_rmerge = dataList[1]
                    datas_rmerge = dataList[3]

                if eachLine.find('Total number of observations') > -1 and num_items == 7:
                    data_num_unmerged = dataList[4]
                    datas_num_unmerged = dataList[6]

                if eachLine.find('Total number unique') > -1 and num_items == 6:
                    data_num = dataList[3]
                    datas_num = dataList[5]

                if eachLine.find('Mean((I)/sd(I))') > -1 and num_items == 4:
                    data_ioversig = dataList[1]
                    datas_ioversig = dataList[3]

                if eachLine.find('Completeness') > -1 and num_items == 4:
                    data_percentobs = dataList[1]
                    datas_percentobs = dataList[3]

                if eachLine.find('Multiplicity') > -1 and num_items == 4:
                    data_redund = dataList[1]
                    datas_redund = dataList[3]

    ############################################
    # Parse SCALEPACK log file (if available)  #
    ############################################

    # Tested on HKL2000/SCALEPACK log from a beamline 

    float_num_refs_theoretical = '?'

    if mergeprog == 'scalepack' and datalogfile != 'none':

        print 'Parsing SCALEPACK log file'

        file = open(datalogfile,'r')
        allLines = file.readlines()
        file.close()

        readredundancy = 'no'
        restable = 'no'
        readcompleteness = 'no'
        readnumberrefs = 'no'

        for eachLine in allLines:

            # Resolution, Rmerge, Av-I/error

            if restable == 'yes':
                dataList = eachLine.split()
                num_args = len(dataList)

                if num_args == 8 and dataList[0] == 'All' and dataList[1] == 'reflections':
                    data_rmerge = dataList[6]

                    mean_i = dataList[2]
                    mean_sigi = dataList[3]
                    mean_i = float(mean_i)
                    mean_sigi = float(mean_sigi)
                    netIoveravsigmaI = mean_i/mean_sigi
                    netIoveravsigmaI_out = round(netIoveravsigmaI,1)
                    data_ioversig = str(netIoveravsigmaI_out)

                    restable = 'no'

                if num_args == 8 and dataList[0] != 'All':
                    datas_rmerge = dataList[6]

                    mean_i = dataList[2]
                    mean_sigi = dataList[3]
                    mean_i = float(mean_i)
                    mean_sigi = float(mean_sigi)
                    netIoveravsigmaI = mean_i/mean_sigi
                    netIoveravsigmaI_out = round(netIoveravsigmaI,1)
                    datas_ioversig = str(netIoveravsigmaI_out)

                    datas_rlow = dataList[0]
                    datas_rhigh = dataList[1]
                    data_rhigh = dataList[1]

                if data_rlow == '?':
                    data_rlow = dataList[0]

            if eachLine.find('limit    Angstrom       I   error   stat. Chi**2  R-fac  R-fac') > -1:
                restable = 'yes'

            # Redundancy

            if readredundancy == 'yes':
                dataList = eachLine.split()
                num_args = len(dataList)

                if dataList[0] == 'All' and dataList[1] == 'hkl' and num_args == 3:
                    data_redund = dataList[2]
                    readredundancy = 'no'

                if dataList[0] != 'All' and num_args == 3:
                    datas_redund = dataList[2]

            if eachLine.find('Average Redundancy Per Shell') > -1:
                readredundancy = 'yes'

            # Completeness

            if readcompleteness == 'yes':
                dataList = eachLine.split()
                num_args = len(dataList)

                if dataList[0] == 'All' and dataList[1] == 'hkl' and num_args == 13:
                    data_percentobs = dataList[12]
                    readcompleteness = 'no'

                if dataList[0] != 'All' and num_args == 13:
                    datas_percentobs = dataList[12]

            if eachLine.find('% of reflections with given No. of observations') > -1:
                readcompleteness = 'yes'

            # Total reflection count

            if eachLine.find('All films') > -1:
                dataList = eachLine.split()
                num_args = len(dataList)

                if num_args > 1:
                    data_num_unmerged = dataList[2]

            # Number of unique reflections

            if readnumberrefs == 'yes':
                dataList = eachLine.split()
                num_args = len(dataList)

                if num_args == 11:
                    if dataList[0] == 'All':
                        data_num = dataList[10]
                        readnumberrefs = 'no'

            if eachLine.find('No. of reflections with given No. of observations') > -1:
                readnumberrefs = 'yes'

    ###################################################
    # Parse template information file (if available)  #
    ###################################################

    if templatefile != 'none':

        print 'Parsing template file'

        file = open(templatefile,'r')
        allLines = file.readlines()
        file.close()

        for eachLine in allLines:

            dataList = eachLine.split(':')
            length = len(dataList)

            # Parse each item from the template file

            if length == 2:

                # Section 1

                if eachLine.find('audit_author_name') > -1:
                    audit_author_name = dataList[1]
                    audit_author_name = audit_author_name.strip()

                if eachLine.find('audit_contact_author_name') > -1:
                    audit_contact_author_name = dataList[1]
                    audit_contact_author_name =audit_contact_author_name.strip()

                if eachLine.find('audit_contact_author_email') > -1:
                    audit_contact_author_email = dataList[1]
                    audit_contact_author_email = audit_contact_author_email.strip()

                if eachLine.find('audit_contact_author_address') > -1:
                    audit_contact_author_address = dataList[1]
                    audit_contact_author_address = audit_contact_author_address.strip()

                if eachLine.find('audit_contact_author_phone') > -1:
                    audit_contact_author_phone = dataList[1]
                    audit_contact_author_phone = audit_contact_author_phone.strip()

                if eachLine.find('audit_contact_author_fax') > -1:
                    audit_contact_author_fax = dataList[1]
                    audit_contact_author_fax = audit_contact_author_fax.strip()

                if eachLine.find('citation_title') > -1:
                    citation_title = dataList[1]
                    citation_title = citation_title.strip()

                if eachLine.find('citation_journal_abbrev') > -1:
                    citation_journal_abbrev = dataList[1]
                    citation_journal_abbrev = citation_journal_abbrev.strip()

                if eachLine.find('citation_journal_volume') > -1:
                    citation_journal_volume = dataList[1]
                    citation_journal_volume = citation_journal_volume.strip()

                if eachLine.find('citation_page_first') > -1:
                    citation_page_first = dataList[1]
                    citation_page_first = citation_page_first.strip()

                if eachLine.find('citation_page_last') > -1:
                    citation_page_last = dataList[1]
                    citation_page_last = citation_page_last.strip()

                if eachLine.find('citation_year') > -1:
                    citation_year = dataList[1]
                    citation_year = citation_year.strip()

                if eachLine.find('citation_author_name') > -1:
                    citation_author_name = dataList[1]
                    citation_author_name = citation_author_name.strip()

                # Section 2

                if eachLine.find('data_collection_temp_K') > -1:
                    data_collection_temp_K = dataList[1]
                    data_collection_temp_K = data_collection_temp_K.strip()

                if eachLine.find('wavelengths') > -1:
                    wavelengths = dataList[1]
                    wavelengths = wavelengths.strip()

                if eachLine.find('data_collection_date') > -1:
                    data_collection_date = dataList[1]
                    data_collection_date = data_collection_date.strip()

                if eachLine.find('beamline') > -1:
                    beamline = dataList[1]
                    beamline = beamline.strip()

                if eachLine.find('detector_type') > -1:
                    detector_type = dataList[1]
                    detector_type = detector_type.strip()

                if eachLine.find('detector_maker') > -1:
                    detector_maker = dataList[1]
                    detector_maker = detector_maker.strip()

                if eachLine.find('monochromator_type') > -1:
                    monochromator_type = dataList[1]
                    monochromator_type = monochromator_type.strip()

                if eachLine.find('xray_method') > -1:
                    xray_method = dataList[1]
                    xray_method = xray_method.strip()

                if eachLine.find('computing_data_collection') > -1:
                    computing_data_collection = dataList[1]
                    computing_data_collection = computing_data_collection.strip()

                if eachLine.find('computing_data_reduction') > -1:
                    computing_data_reduction = dataList[1]
                    computing_data_reduction = computing_data_reduction.strip()

                if eachLine.find('computing_structure_solution') > -1:
                    computing_structure_solution = dataList[1]
                    computing_structure_solution = computing_structure_solution.strip()

                if eachLine.find('computing_molecular_graphics') > -1:
                    computing_molecular_graphics = dataList[1]
                    computing_molecular_graphics = computing_molecular_graphics.strip()

                if eachLine.find('computing_structure_refinement') > -1:
                    computing_structure_refinement = dataList[1]
                    computing_structure_refinement = computing_structure_refinement.strip()

                # Section 3

                if eachLine.find('protein_name') > -1:
                    protein_name = dataList[1]
                    protein_name = protein_name.strip()

                if eachLine.find('protein_ec_number') > -1:
                    protein_ec_number = dataList[1]
                    protein_ec_number = protein_ec_number.strip()

                if eachLine.find('structure_title') > -1:
                    structure_title = dataList[1]
                    structure_title = structure_title.strip()

                if eachLine.find('structure_class') > -1:
                    structure_class = dataList[1]
                    structure_class = structure_class.strip()

                if eachLine.find('structure_keywords') > -1:
                    structure_keywords = dataList[1]
                    structure_keywords = structure_keywords.strip()

                if eachLine.find('biological_unit') > -1:
                    biological_unit = dataList[1]
                    biological_unit = biological_unit.strip()

                if eachLine.find('sequence_databasename') > -1:
                    sequence_databasename = dataList[1]
                    sequence_databasename = sequence_databasename.strip()

                if eachLine.find('sequence_databasecode') > -1:
                    sequence_databasecode = dataList[1]
                    sequence_databasecode = sequence_databasecode.strip()

                if eachLine.find('source_common_name') > -1:
                    source_common_name = dataList[1]
                    source_common_name = source_common_name.strip()

                if eachLine.find('source_scientific_name') > -1:
                    source_scientific_name = dataList[1]
                    source_scientific_name = source_scientific_name.strip()

                if eachLine.find('source_gene_name') > -1:
                    source_gene_name = dataList[1]
                    source_gene_name = source_gene_name.strip()

                if eachLine.find('source_host_common_name') > -1:
                    source_host_common_name = dataList[1]
                    source_host_common_name = source_host_common_name.strip()

                if eachLine.find('source_host_scientific_name') > -1:
                    source_host_scientific_name = dataList[1]
                    source_host_scientific_name = source_host_scientific_name.strip()

                # Section 4

                if eachLine.find('exptl_crystal_grow_method') > -1:
                    exptl_crystal_grow_method = dataList[1]
                    exptl_crystal_grow_method = exptl_crystal_grow_method.strip()            

                if eachLine.find('exptl_crystal_grow_pH') > -1:
                    exptl_crystal_grow_pH = dataList[1]
                    exptl_crystal_grow_pH = exptl_crystal_grow_pH.strip()

                if eachLine.find('exptl_crystal_grow_temp') > -1:
                    exptl_crystal_grow_temp = dataList[1]
                    exptl_crystal_grow_temp = exptl_crystal_grow_temp.strip()

                if eachLine.find('exptl_crystal_grow_components') > -1:
                    exptl_crystal_grow_components = dataList[1]
                    exptl_crystal_grow_components = exptl_crystal_grow_components.strip()


    #############################################
    # R-factor and stereochemistry calculation  #
    #############################################

    print 'Running R-factor and stereochemistry calculations'

    filename = 'mi_runrefmac5.inp'
    file = open(filename,'w')
    file.write('LABIN FP=')
    file.write(famp)
    file.write(' SIGFP=')
    file.write(sd)

    if freer != '?':
        file.write(' FREE=')
        file.write(freer)

    file.write('\nLABOUT FC=FC PHIC=PHIC DELFWT=DELFWT PHDELWT=PHDELWT FWT=FWT PHWT=PHWT FOM=FOM\n')
    file.write('FREE ')
    file.write(nfree_exclude)
    file.write('\nREFI TYPE RESTrained\n')
    file.write('REFI RESI MLKF\n')

    if ref_anisoflag == 'no':
        file.write('REFI BREF ISOT METH CGMAT\n')
    else:
        file.write('REFI BREF ANISotropic METH CGMAT\n')  

    if ref_bulksolvent == 'babinet':
        file.write('SCAL TYPE BULK LSSC ANIS\n')
        file.write('SOLVENT NO\n')

    if ref_bulksolvent == 'mask':
        file.write('SCAL TYPE SIMPLE LSSC ANIS\n')
        file.write('SOLVENT YES\n')

    if ref_bulksolvent == 'fixedbabinet':
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
    file.write('DNAME output\n')
    file.write('END\n')
    file.close()

    fileexists = os.path.exists(liganddir)
    if fileexists == 0:
        runrefmac5 = 'refmac5 HKLIN ' + local_mtzfile +\
                     ' XYZIN temp_use.pdb HKLOUT temp_ref.mtz XYZOUT temp_ref.pdb < mi_runrefmac5.inp > mi_refmac.out'
    else:
        runrefmac5 = 'refmac5 HKLIN ' + local_mtzfile + ' LIBIN ' + liganddir +\
                     ' XYZIN temp_use.pdb HKLOUT temp_ref.mtz XYZOUT temp_ref.pdb < mi_runrefmac5.inp > mi_refmac.out'    

    os.system(runrefmac5)

    fileexists = os.path.exists('output.refmac')
    if fileexists == 0:
        print '\nREFMAC5 calculation failed - check mi_refmac.out\n'
        print 'The usual problem is atom names inconsistent with the PDB residue code\n'
        time.sleep(4)
        return 1

    fileexists = os.path.exists('temp_ref.pdb')
    if fileexists == 0:    
        print '\nREFMAC5 calculation failed - check mi_refmac.out\n'
        print 'The usual problem is atom names inconsistent with the PDB residue code\n'
        time.sleep(4)
        return 1

    fileexists = os.path.exists('mi_runrefmac5.inp')
    if fileexists != 0:
        os.remove('mi_runrefmac5.inp')

    fileexists = os.path.exists('temp.lib')
    if fileexists != 0:     
        os.remove('temp.lib')    

    # Parse the required deposition information from output.refmac
    # Note: this file already conveniently uses mmCIF tags..

    file = open('output.refmac','r')
    allLines = file.readlines()
    file.close()

    for eachLine in allLines:

        reflineList = eachLine.split()

        if eachLine.find('_refine.ls_R_factor_R_all') > -1:
            ref_rall = reflineList[1]

        if eachLine.find('_refine.ls_R_factor_R_free') > -1:
            ref_rfree = reflineList[1]

        if eachLine.find('_refine.ls_R_factor_R_work') > -1:
            ref_rwork = reflineList[1]

        if eachLine.find('_refine.ls_d_res_low') > -1:
            ref_dlow = reflineList[1]

        if eachLine.find('_refine.ls_d_res_high') > -1:
            ref_dhigh = reflineList[1]

        if eachLine.find('_refine.ls_R_factor_R_all') > -1:
            ref_rall = reflineList[1]
            ref_robs = reflineList[1]

        if eachLine.find('_refine.ls_number_reflns_R_free') > -1:
            ref_numfree = reflineList[1]

        if eachLine.find('_refine.ls_number_reflns_R_work') > -1:
            ref_numwork = reflineList[1]

        if eachLine.find('_refine.ls_number_reflns_obs') > -1:
            ref_numobs = reflineList[1]
            ref_numall = reflineList[1]

        if eachLine.find('_refine.B_iso_mean') > -1:
            ref_bmean = reflineList[1]

        if eachLine.find('_refine.aniso_B[1][1]') > -1:
            ref_b11 = reflineList[1]

        if eachLine.find('_refine.aniso_B[2][2]') > -1:
            ref_b22 = reflineList[1]

        if eachLine.find('_refine.aniso_B[3][3]') > -1:
            ref_b33 = reflineList[1]

        if eachLine.find('_refine.aniso_B[1][2]') > -1:
            ref_b12 = reflineList[1]

        if eachLine.find('_refine.aniso_B[1][3]') > -1:
            ref_b13 = reflineList[1]

        if eachLine.find('_refine.aniso_B[2][3]') > -1:
            ref_b23 = reflineList[1]       

        if eachLine.find('_cell.length_a') > -1:
            acell = reflineList[1]  

        if eachLine.find('_cell.length_b') > -1:
            bcell = reflineList[1]

        if eachLine.find('_cell.length_c') > -1:
            ccell = reflineList[1]

        if eachLine.find('_cell.angle_alpha') > -1:
            alpha = reflineList[1]

        if eachLine.find('_cell.angle_beta') > -1:
            beta = reflineList[1]

        if eachLine.find('_cell.angle_gamma') > -1:
            gamma = reflineList[1]

        if eachLine.find('_refine.solvent_vdw_probe_radii') > -1:
            ref_solvent_vdw_probe_radii = reflineList[1]

        if eachLine.find('_refine.solvent_ion_probe_radii') > -1:
            ref_solvent_ion_probe_radii = reflineList[1]

        if eachLine.find('_refine.solvent_shrinkage_radii') > -1:
            ref_solvent_shrinkage_radii = reflineList[1]    

        if eachLine.find('_refine.solvent_model_param_ksol') > -1:
            ref_ksolv = reflineList[1]

        if eachLine.find('_refine.solvent_model_param_bsol') > -1:
            ref_bsolv = reflineList[1]    

        if eachLine.find('r_bond_refined_d') > -1:
            ref_dbond = reflineList[2]

        if eachLine.find('r_angle_refined_deg') > -1:
            ref_dangle = reflineList[2]

        if eachLine.find('r_dihedral_angle_1_deg') > -1:
            ref_dtorsion = reflineList[2]

        if eachLine.find('r_chiral_restr') > -1:
            ref_dchiral = reflineList[2]

        if eachLine.find('r_gen_planes_refined') > -1:
            ref_dplane = reflineList[2]

        if eachLine.find('r_mcbond_it') > -1:
            ref_bmbond = reflineList[2]

        if eachLine.find('r_mcangle_it') > -1:
            ref_bmangle = reflineList[2]

        if eachLine.find('r_scbond_it') > -1:
            ref_bsbond = reflineList[2]

        if eachLine.find('r_scangle_it') > -1:
            ref_bsangle = reflineList[2]

    # Fix for cubic space groups

    if ref_b11 == '?':
        ref_b11 = '0.00'

    if ref_b22 == '?':
        ref_b22 = '0.00'

    if ref_b33 == '?':
        ref_b33 = '0.00'

    if ref_b12 == '?':
        ref_b12 = '0.00'

    if ref_b13 == '?':
        ref_b13 = '0.00'

    if ref_b23 == '?':
        ref_b23 = '0.00'

    # Check for incorrect use of cross-validation flag (Rwork differs from Rfree by less than 1%)

    if ref_rfree != '?' and ref_rwork != '?':
        float_ref_rfree = float(ref_rfree)
        float_ref_rwork = float(ref_rwork)
        float_rdif = float_ref_rfree - float_ref_rwork

        if float_rdif < 0.01:
            print '\nRwork and Rfree are extremely close !'
            print 'The default assignment of Rfree flags (0) by script parameter nfree_exclude'
            print 'needs to be changed\n'
            time.sleep(4)

    # Check for calculations where no Rfree flags were assigned

    if ref_numall == '?' and ref_numfree == '?':
        ref_numall = ref_numwork
        ref_rfree = '0'

    ###########################################################
    # Check mi_refmac.out for Babinet bulk solvent parameters #
    # and get  *correct* resolution limits                    #
    ###########################################################

    file = open('mi_refmac.out','r')
    allLines = file.readlines()
    file.close()

    read_data = 'no'

    for eachLine in allLines:

        if eachLine.find('Babinet"s bulk solvent:') > -1:
            ref_ksolv = eachLine[33:40]
            ref_ksolv = ref_ksolv.strip()
            ref_bsolv = eachLine[47:54]
            ref_bsolv = ref_bsolv.strip()

        if read_data == 'yes':

            reflineList = eachLine.split()
            number_items = len(reflineList)

            if number_items == 8:
                ref_dlow = reflineList[3]
                ref_dhigh = reflineList[5]

                # Force o/p format to print 2dp

                ref_dlow = float(ref_dlow)
                ref_dhigh = float(ref_dhigh)
                ref_dlow = round(ref_dlow,2)
                ref_dhigh = round(ref_dhigh,2)
                ref_dlow = str(ref_dlow)
                ref_dhigh = str(ref_dhigh)

                aLine = ref_dlow.split('.')
                number_places = len(aLine[1])
                if number_places == 1:
                    ref_dlow_out = ref_dlow + '0'
                    ref_dlow = ref_dlow_out
                if number_places == 0:
                    ref_dlow_out = ref_dlow + '00'
                    ref_dlow = ref_dlow_out

                aLine = ref_dhigh.split('.')
                number_places = len(aLine[1])
                if number_places == 1:
                    ref_dhigh_out = ref_dhigh + '0'
                    ref_dhigh = ref_dhigh_out
                if number_places == 0:
                    ref_dhigh_out = ref_dhigh + '00'
                    ref_dhigh = ref_dhigh_out

                read_data = 'no'

        if eachLine.find('*  Resolution Range') > -1:
            read_data = 'yes'

    ###########################################################
    # Compute refinement data completeness versus theoretical #
    ###########################################################

    float_num_refs_theoretical = '?'

    filename = 'mi_rununique.inp'
    file = open(filename, 'w')
    file.write('TITLE  unique_data\n')
    file.write('LABOUT  F=XFUNI SIGF=XSIGFUNI\n')
    file.write('RESOLUTION ')
    file.write(ref_dhigh)
    file.write('\nSYMM ')
    file.write(spgno)
    file.write('\nCELL ')
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
    file.write('\nEND\n')
    file.close()

    rununique = 'unique HKLOUT temp_unique.mtz < mi_rununique.inp > mi_unique.out'
    os.system(rununique)

    fileexists = os.path.exists('mi_unique.out')
    if fileexists == 0:
        print '\nrefinement data completeness calculation failed'
        time.sleep(4)
        return 1
    else:
        file = open('mi_unique.out','r')
        allLines = file.readlines()
        file.close()

        for eachLine in allLines:

            if eachLine.find('reflections within resolution limits') > -1:
                dataList = eachLine.split()
                num_refs_theoretical = dataList[0]
                float_num_refs_theoretical = float(num_refs_theoretical)

        if float_num_refs_theoretical != '?':
            float_data_refinement = float(ref_numall)
            float_data_percentobs = 100.0 * float_data_refinement/float_num_refs_theoretical
            ref_percent = round(float_data_percentobs,1)
            ref_percent = str(ref_percent)

    os.remove('mi_unique.out')
    os.remove('mi_rununique.inp')

    fileexists = os.path.exists('temp_unique.mtz')
    if fileexists != 0:
        os.remove('temp_unique.mtz')

    ####################################################
    # Compute phi-psi statistics using Richardson data #
    ####################################################

    # Check MIFit installation to access phi-psi data files

    if mifit_root != 'none':
        mifit_root_data = os.path.join(mifit_root,'data')
        phipsi_gen_datafile = os.path.join(mifit_root_data,'rama500-general.data')
        phipsi_gly_datafile = os.path.join(mifit_root_data,'rama500-gly-sym-nosec.data')
        phipsi_pro_datafile = os.path.join(mifit_root_data,'rama500-pro.data')
    else:
        phipsi_gen_datafile = 'none'
        phipsi_gly_datafile = 'none'
        phipsi_pro_datafile = 'none'

    # If we have MIFit environment read Richardson data

    fileexists = os.path.exists(phipsi_gen_datafile)
    if fileexists != 0:

        print 'Checking phi-psi angles'

        aList_phi_all = []
        aList_psi_all = []
        aList_phipsi_prob_all = []  
        aList_phi_gly = []
        aList_psi_gly = []
        aList_phipsi_prob_gly = []
        aList_phi_pro = []
        aList_psi_pro = []
        aList_phipsi_prob_pro = []
        aList_rama_chain = []
        aList_rama_resno = []
        aList_rama_resname = []

        amino_acid_count = 0.0

        # allowed phi-psi: gen 99.95% data, 41.5% area, GLY 99.8% data, 63% area, PRO 99.8% data, 18.1% area

        phipsi_thresh_gen = 0.00847
        phipsi_thresh_gly = 0.00384
        phipsi_thresh_pro = 0.0015

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

        file = open('temp_use.pdb')
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

        # Count phi-psi errors

        number_phi_psi_errors = len(aList_rama_chain)
        ref_number_phi_psi_errors = str(number_phi_psi_errors)

    ####################################
    # Execution of CCP4/MATTHEWS_COEF  #             
    ####################################

    if seqfile != 'none':

        print 'Running solvent volume calculations'

        # Note that this uses the FASTA sequence, which should be the protein in the crystal sample

        num_res_in_crystal_au = sequence_length * number_chains
        num_res_in_crystal_au = str(num_res_in_crystal_au)

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
        file.write(spgno)
        file.write('\nNRES ')
        file.write(num_res_in_crystal_au)
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

            if eachLine.find('The Matthews Coefficient is') > -1:

                dataList = eachLine.split(':')
                matthews_coef = dataList[1]
                matthews_coef = matthews_coef.strip()

            if eachLine.find('Assuming protein density') > -1:

                dataList = eachLine.split(':')
                solvent_percent = dataList[1]
                solvent_percent = solvent_percent.strip() 

        os.remove('mi_runmatthews_coef.inp')
        os.remove('mi_matthews_coef.out')

    ###################################################################
    # Obtain space group name in standard (not H-M) notation for HTML #
    ###################################################################

    if spgno == '1': spgname_standard = 'P1'
    if spgno == '3': spgname_standard = 'P2'
    if spgno == '4': spgname_standard = 'P2<sub>1</sub>'
    if spgno == '5': spgname_standard = 'C2'
    if spgno == '16': spgname_standard = 'P222'
    if spgno == '17': spgname_standard = 'P222<sub>1</sub>'
    if spgno == '18': spgname_standard = 'P2<sub>1</sub>2<sub>1</sub>2'
    if spgno == '19': spgname_standard = 'P2<sub>1</sub>2<sub>1</sub>2<sub>1</sub>'
    if spgno == '20': spgname_standard = 'C222<sub>1</sub>'
    if spgno == '21': spgname_standard = 'C222'
    if spgno == '22': spgname_standard = 'F222'
    if spgno == '23': spgname_standard = 'I222'
    if spgno == '24': spgname_standard = 'I2<sub>1</sub>2<sub>1</sub>2<sub>1</sub>'
    if spgno == '75': spgname_standard = 'P4'
    if spgno == '76': spgname_standard = 'P4<sub>1</sub>'
    if spgno == '77': spgname_standard = 'P4<sub>2</sub>'
    if spgno == '78': spgname_standard = 'P4<sub>3</sub>'
    if spgno == '79': spgname_standard = 'I4'
    if spgno == '80': spgname_standard = 'I4<sub>1</sub>'
    if spgno == '89': spgname_standard = 'P422'
    if spgno == '90': spgname_standard = 'P42<sub>1</sub>2'
    if spgno == '91': spgname_standard = 'P4<sub>1</sub>22'
    if spgno == '92': spgname_standard = 'P4<sub>1</sub>2<sub>1</sub>2'
    if spgno == '93': spgname_standard = 'P4<sub>2</sub>22'
    if spgno == '94': spgname_standard = 'P4<sub>2</sub>2<sub>1</sub>2'
    if spgno == '95': spgname_standard = 'P4<sub>3</sub>22'
    if spgno == '96': spgname_standard = 'P4<sub>3</sub>2<sub>1</sub>2'
    if spgno == '97': spgname_standard = 'I422'
    if spgno == '98': spgname_standard = 'I4<sub>1</sub>22'
    if spgno == '143': spgname_standard = 'P3'
    if spgno == '144': spgname_standard = 'P3<sub>1</sub>'
    if spgno == '145': spgname_standard = 'P3<sub>2</sub>'
    if spgno == '146': spgname_standard = 'R3'
    if spgno == '149': spgname_standard = 'P3<sub>1</sub>2'
    if spgno == '150': spgname_standard = 'P32<sub>1</sub>'
    if spgno == '151': spgname_standard = 'P3<sub>1</sub>12'
    if spgno == '152': spgname_standard = 'P3<sub>1</sub>21'
    if spgno == '153': spgname_standard = 'P3<sub>2</sub>12'
    if spgno == '154': spgname_standard = 'P3<sub>2</sub>21'
    if spgno == '155': spgname_standard = 'R32'
    if spgno == '168': spgname_standard = 'P6'
    if spgno == '169': spgname_standard = 'P6<sub>1</sub>'
    if spgno == '170': spgname_standard = 'P6<sub>5</sub>'
    if spgno == '171': spgname_standard = 'P6<sub>2</sub>'
    if spgno == '172': spgname_standard = 'P6<sub>4</sub>'
    if spgno == '173': spgname_standard = 'P6<sub>3</sub>'
    if spgno == '177': spgname_standard = 'P622'
    if spgno == '178': spgname_standard = 'P6<sub>1</sub>22'
    if spgno == '179': spgname_standard = 'P6<sub>5</sub>22'
    if spgno == '180': spgname_standard = 'P6<sub>2</sub>22'
    if spgno == '181': spgname_standard = 'P6<sub>4</sub>22'
    if spgno == '182': spgname_standard = 'P6<sub>3</sub>22'
    if spgno == '195': spgname_standard = 'P23'
    if spgno == '196': spgname_standard = 'F23'
    if spgno == '197': spgname_standard = 'I23'
    if spgno == '198': spgname_standard = 'P2<sub>1</sub>3'
    if spgno == '199': spgname_standard = 'I2<sub>1</sub>3'
    if spgno == '207': spgname_standard = 'P432'
    if spgno == '208': spgname_standard = 'P4<sub>2</sub>32'
    if spgno == '209': spgname_standard = 'F432'
    if spgno == '210': spgname_standard = 'F4<sub>1</sub>32'
    if spgno == '211': spgname_standard = 'I432'
    if spgno == '212': spgname_standard = 'P4<sub>3</sub>32'
    if spgno == '213': spgname_standard = 'P4<sub>1</sub>32'
    if spgno == '214': spgname_standard = 'I4<sub>1</sub>32'

    #####################
    # Write everything  #
    #####################

    if write_cif != 'no':
        print 'Writing deposition file'

    cifdepositfile = rootname + '.cif'

    file = open(cifdepositfile, 'w')
    file.write('data_structure_1\n')

    file.write('#\n')
    file.write('##############################\n')
    file.write('# Release information etc    #\n')
    file.write('##############################\n')
    file.write('#\n')
    file.write('_audit_author.name ')
    file.write(quote)
    file.write(audit_author_name)
    file.write(quote)
    file.write('\n#\n')
    file.write('_audit_contact_author.name ')
    file.write(quote)
    file.write(audit_contact_author_name)
    file.write(quote)
    file.write('\n_audit_contact_author.email ')
    file.write(quote)
    file.write(audit_contact_author_email)
    file.write(quote)
    file.write('\n_audit_contact_author.address\n')
    file.write('; ')
    file.write(audit_contact_author_address)
    file.write('\n;\n')           
    file.write('_audit_contact_author.phone ')
    file.write(quote)
    file.write(audit_contact_author_phone)
    file.write(quote)
    file.write('\n_audit_contact_author.fax ')
    file.write(quote) 
    file.write(audit_contact_author_fax)
    file.write(quote)
    file.write('\n#\n')
    file.write('_pdbx_database_status.dep_release_code_coordinates ')
    file.write(quote)
    file.write('HPUB')
    file.write(quote)
    file.write('\n_pdbx_database_status.dep_release_code_struct_fact ')
    file.write(quote)
    file.write('HPUB')
    file.write(quote)
    file.write('\n_pdbx_database_status.dep_release_code_sequence    ')
    file.write(quote)
    file.write('REL')
    file.write(quote)
    file.write('\n#\n')

    if citation_title != '?':

        file.write('################################\n')
        file.write('# Citation and author          #\n')
        file.write('################################\n')
        file.write('#\n')
        file.write('_citation.id primary\n')
        file.write('_citation.title ')
        file.write(quote)
        file.write(citation_title)
        file.write(quote)
        file.write('\n_citation.journal_abbrev ')
        file.write(quote)
        file.write(citation_journal_abbrev)
        file.write(quote)
        file.write('\n_citation.journal_volume ')
        file.write(quote)
        file.write(citation_journal_volume)
        file.write(quote)
        file.write('\n_citation.page_first ')
        file.write(citation_page_first)
        file.write('\n_citation.page_last ')
        file.write(citation_page_last)
        file.write('\n_citation.year ')
        file.write(citation_year)
        file.write('\n#\n')
        file.write('_citation_author.citation_id primary\n')
        file.write('_citation_author.name ')
        file.write(quote)
        file.write(citation_author_name)
        file.write(quote)
        file.write('\n#\n')

    if citation_journal_abbrev == 'unpublished' or citation_journal_abbrev == 'UNPUBLISHED':

        file.write('################################\n')
        file.write('# Citation and author          #\n')
        file.write('################################\n')
        file.write('#\n')
        file.write('_citation.id primary\n')
        file.write('_citation.journal_abbrev unpublished')

    file.write('##############################################################\n')
    file.write('#                                                            #\n')
    file.write('# Extra data collection and processing information for PDB   #\n')
    file.write('#                                                            #\n')
    file.write('##############################################################\n')
    file.write('#\n')
    file.write('_exptl.entry_id 1\n')
    file.write('_exptl.method ')
    file.write(quote)
    file.write('x-ray diffraction')
    file.write(quote)
    file.write('\n_exptl.crystals_number 1\n')
    file.write('#\n')
    file.write('_diffrn.id  1\n')
    file.write('_diffrn.ambient_temp ')
    file.write(data_collection_temp_K)
    file.write('\n')

    file.write('_diffrn_source.diffrn_id 1\n')

    if wavelengths != '1.54':

        file.write('_diffrn_source.source ')
        file.write(quote)
        file.write('Synchrotron')
        file.write(quote)
        file.write('\n_diffrn_source.type ')
        file.write(quote)
        file.write(beamline)
        file.write(quote)
        file.write('\n_diffrn_detector.diffrn_id 1\n')
        file.write('_diffrn_detector.pdbx_collection_date ')
        file.write(quote)
        file.write(data_collection_date)
        file.write(quote)
        file.write('\n_diffrn_detector.detector ')
        file.write(quote)
        file.write(detector_type)
        file.write(quote)    
        file.write('\n_diffrn_detector.type ')
        file.write(quote)
        file.write(detector_maker)
        file.write(quote)
        file.write('\n#\n')    
    else:
        file.write('_diffrn_source.type ')
        file.write(quote)
        file.write('Rotating anode')
        file.write(quote)
        file.write('\n_diffrn_detector.diffrn_id 1\n')
        file.write('_diffrn_detector.pdbx_collection_date ')
        file.write(quote)
        file.write(data_collection_date)
        file.write(quote)
        file.write('\n_diffrn_detector.detector ')
        file.write(quote)
        file.write(detector_type)
        file.write(quote)    
        file.write('\n_diffrn_detector.type ')
        file.write(quote)
        file.write(detector_maker)
        file.write(quote)
        file.write('\n#\n')

    file.write('_diffrn_radiation.diffrn_id  1\n')
    file.write('_diffrn_radiation.wavelength_id  1\n')

    if xray_method == 'MAD':
        file.write('_diffrn_radiation.pdbx_diffrn_protocol MAD\n')
    else:
        file.write('_diffrn_radiation.pdbx_diffrn_protocol ')
        file.write(quote)
        file.write('Single wavelength')
        file.write(quote)

    file.write('\n_diffrn_radiation.monochromator ')
    file.write(quote)
    file.write(monochromator_type)
    file.write(quote)
    file.write('\n_diffrn_radiation.pdbx_wavelength_list ')
    file.write(quote)
    file.write(wavelengths)
    file.write(quote)
    file.write('\n_diffrn_reflns.diffrn_id   1\n')
    file.write('_diffrn_reflns.number ')
    file.write(data_num_unmerged)
    file.write('\n#\n')
    file.write('_computing.entry_id 1\n')
    file.write('_computing.data_collection ')
    file.write(quote)
    file.write(computing_data_collection)
    file.write(quote)
    file.write('\n_computing.data_reduction ')
    file.write(quote)
    file.write(computing_data_reduction)
    file.write(quote)
    file.write('\n_computing.structure_solution ')
    file.write(quote)
    file.write(computing_structure_solution)
    file.write(quote)
    file.write('\n_computing.molecular_graphics ')
    file.write(quote)
    file.write(computing_molecular_graphics)
    file.write(quote)
    file.write('\n_computing.structure_refinement ')
    file.write(quote)
    file.write(computing_structure_refinement)
    file.write(quote)

    file.write('\n#\n')
    file.write('################################\n')
    file.write('# Sequence information         #\n')
    file.write('################################\n')
    file.write('#\n')

    # Write protein sequence (entity 1)

    file.write('_entity_poly.entity_id     1\n')
    file.write('_entity_poly.pdbx_seq_one_letter_code\n')
    file.write(';')

    if seqfile != 'none':

        # Pad out sequence with spaces

        remainder = sequence_length%60
        pad = 60 - remainder

        count = 0
        while count < pad:
            seqList.append(' ')
            count  = count + 1

        sequence_length = len(seqList)

        # Print out in blocks of 60 amino acids

        count = 0
        while count < sequence_length:

            aa = seqList[count]
            file.write(aa)

            count = count + 1
            remainder = count%60

            if remainder == 0:
                file.write('\n')

        file.write(';\n')

    else:

        file.write('?\n')
        file.write(';\n')

    # Write the PDB chain-ids to which the sequence corresponds

    file.write('_entity_poly.pdbx_strand_id ')
    file.write(quote)

    count = 0
    while count < number_chains:
        chain_out = aList_chains[count]
        pr_chain_out = ' ' + chain_out + ' '
        file.write(pr_chain_out)

        if count != number_chains - 1:
            file.write(',')

        count = count + 1

    file.write(quote)
    file.write('\n')

    # Entity annotation

    file.write('#\n')
    file.write('################################\n')
    file.write('# Entity information           #\n')
    file.write('################################\n')
    file.write('#\n')
    file.write('loop_\n')
    file.write('_entity.id\n')
    file.write('_entity.pdbx_description\n')
    file.write('_entity.type\n')
    file.write('_entity.pdbx_ec\n')

    # Protein

    file.write('1')
    file.write('\n; ')
    file.write(protein_name)
    file.write('\n;\n')
    file.write('                  polymer     ')
    file.write(protein_ec_number)
    file.write('\n')

    # Ligand list

    number_entities = len(aList_hets)
    entity_list_number = 1
    count = 0

    if number_entities > 0:

        while count < number_entities:

            entity_number = aList_hets_number[count]
            entity_name = aList_hets_names[count]

            file.write(entity_number)
            file.write('\n; ')
            file.write(entity_name)
            file.write('\n;\n')
            file.write('                  non-polymer ?\n')

            count = count + 1

    # Waters

    if water_flag == 'yes':

        file.write(water_entity_list_number)
        file.write(' water           water       ?\n') 


    # Add the entity keyword ids as separate table for ADIT 

    file.write('#\n')
    file.write('loop_\n')
    file.write('_entity_keywords.entity_id\n')
    file.write('1\n')

    count = 0
    if number_entities > 0:

        while count < number_entities:        
            entity_number = aList_hets_number[count]

            file.write(entity_number)
            file.write('\n')

            count = count + 1

    if water_flag == 'yes':    
        file.write(water_entity_list_number)
        file.write('\n')

    file.write('#\n')
    file.write('##############################################################\n')
    file.write('# Structure annotation                                       #\n')
    file.write('# note: _struct_keywords.pdbx_keywords maps to HEADER        #\n')
    file.write('# and should report function. Use broad enzyme               #\n')
    file.write('# classification or "Structural genomics, unknown function"  #\n')
    file.write('#       _struct_keywords.text can be any words               #\n') 
    file.write('##############################################################\n')
    file.write('#\n')
    file.write('_struct.entry_id 1\n')
    file.write('_struct.title ')
    file.write(quote)
    file.write(structure_title)
    file.write(quote)
    file.write('\n#\n')
    file.write('_struct_keywords.entry_id 1\n')
    file.write('_struct_keywords.pdbx_keywords ')
    file.write(quote)
    file.write(structure_class)
    file.write(quote)
    file.write('\n_struct_keywords.text ')
    file.write(quote)
    file.write(structure_keywords)
    file.write(quote)
    file.write('\n#\n')
    file.write('##############################################\n')
    file.write('# Unique entity map                          #\n') 
    file.write('# note: points to _atom_site.label_asym_id   #\n')
    file.write('##############################################\n')
    file.write('#\n')
    file.write('loop_\n')
    file.write('_struct_asym.id\n')
    file.write('_struct_asym.entity_id\n')

    # Protein

    count = 0
    if number_chains > 0:

        while count < number_chains:

            entity_list_number = '1'
            entity_list_asym = aList_chains[count]

            file.write(entity_list_asym)
            file.write('     ')
            file.write(entity_list_number)
            file.write('\n')

            count = count + 1

    # Ligand list

    count = 0
    if number_entities > 0:

        while count < number_entities:

            entity_list_number = aList_hets_number[count]
            entity_list_asym =  aList_hets_asym[count]

            file.write(entity_list_asym)
            file.write(entity_list_number)
            file.write('\n')

            count = count + 1

    # Waters

    if water_flag == 'yes':

        file.write('W     ')
        file.write(water_entity_list_number)
        file.write('\n')

    file.write('#\n')
    file.write('##############################################\n')
    file.write('# Asymmetric unit description                #\n')
    file.write('##############################################\n')
    file.write('#\n')
    file.write('_struct_biol.id 1\n')
    file.write('_struct_biol.details ')
    file.write(quote)
    file.write(biological_unit)
    file.write(quote)
    file.write('\n#\n')
    file.write('##############################################\n')
    file.write('# Database sequence reference                #\n')
    file.write('##############################################\n')
    file.write('#\n')
    file.write('loop_\n')
    file.write('_struct_ref.id\n')
    file.write('_struct_ref.entity_id\n')
    file.write('_struct_ref.db_name\n')
    file.write('_struct_ref.db_code\n')
    file.write('1 1 ')
    file.write(quote)
    file.write(sequence_databasename)
    file.write(quote)
    file.write(' ')
    file.write(quote)
    file.write(sequence_databasecode)
    file.write(quote)
    file.write('\n#\n')
    file.write('###########################\n')
    file.write('# Source information      #\n')
    file.write('###########################\n')
    file.write('#\n')
    file.write('_entity_src_gen.entity_id                1\n')
    file.write('_entity_src_gen.gene_src_common_name     ')
    file.write(quote)
    file.write(source_common_name)
    file.write(quote)
    file.write('\n_entity_src_gen.pdbx_gene_src_scientific_name  ')
    file.write(quote)
    file.write(source_scientific_name)
    file.write(quote)
    file.write('\n_entity_src_gen.pdbx_gene_src_gene     ')
    file.write(quote)
    file.write(source_gene_name)
    file.write(quote)
    file.write('\n_entity_src_gen.host_org_common_name   ')
    file.write(quote)
    file.write(source_host_common_name)
    file.write(quote)
    file.write('\n_entity_src_gen.pdbx_host_org_scientific_name ')
    file.write(quote)
    file.write(source_host_scientific_name)
    file.write(quote)

    # Write crystal data

    file.write('\n#\n')
    file.write('##########################\n')
    file.write('# Crystal information    #\n')
    file.write('##########################\n')
    file.write('#\n')
    file.write('_cell.entry_id 1\n')
    file.write('_cell.length_a ')
    file.write(acell)
    file.write('\n_cell.length_b ')
    file.write(bcell)
    file.write('\n_cell.length_c ')
    file.write(ccell)
    file.write('\n_cell.angle_alpha ')
    file.write(alpha)
    file.write('\n_cell.angle_beta ')
    file.write(beta)
    file.write('\n_cell.angle_gamma ')
    file.write(gamma)
    file.write('\n#\n')
    file.write('_symmetry.entry_id 1\n')
    file.write('_symmetry.Int_Tables_number ')
    file.write(spgno)
    file.write('\n_symmetry.space_group_name_H-M ')
    file.write(quote)
    file.write(spgname)
    file.write(quote)
    file.write('\n')

    file.write('#\n')
    file.write('_exptl_crystal.id 1\n')
    file.write('_exptl_crystal.density_percent_sol ')
    file.write(solvent_percent)
    file.write('\n_exptl_crystal.density_Matthews ')
    file.write(matthews_coef )
    file.write('\n#\n')

    file.write('_exptl_crystal_grow.crystal_id 1\n')
    file.write('_exptl_crystal_grow.method ')
    file.write(quote)
    file.write(exptl_crystal_grow_method)
    file.write(quote)
    file.write('\n_exptl_crystal_grow.pH ')
    file.write(exptl_crystal_grow_pH)
    file.write('\n_exptl_crystal_grow.temp ')
    file.write(exptl_crystal_grow_temp)
    file.write('\n_exptl_crystal_grow.pdbx_details ')
    file.write(quote)
    file.write(exptl_crystal_grow_components)
    file.write(quote)
    file.write('\n')

    # Write data collection data

    file.write('#\n')
    file.write('###########################\n')
    file.write('# Data collection         #\n')
    file.write('###########################\n')
    file.write('#\n')
    file.write('# Overall processing statistics\n')
    file.write('#\n')
    file.write('_reflns.entry_id 1\n')
    file.write('_reflns.number_all ')
    file.write(data_num)                  
    file.write('\n_reflns.number_obs ')
    file.write(data_num)

    file.write('\n_reflns.observed_criterion_sigma_F ')
    file.write(truncate_default_f)
    file.write('\n_reflns.observed_criterion_sigma_I ')
    file.write(truncate_default_i)

    file.write('\n_reflns.d_resolution_low ')
    file.write(data_rlow)         
    file.write('\n_reflns.d_resolution_high ')
    file.write(data_rhigh)        
    file.write('\n_reflns.percent_possible_obs ')
    file.write(data_percentobs)
    file.write('\n_reflns.pdbx_redundancy ')
    file.write(data_redund)    
    file.write('\n_reflns.pdbx_Rmerge_I_obs ')
    file.write(data_rmerge)    
    file.write('\n_reflns.pdbx_netI_over_av_sigmaI ')
    file.write(data_ioversig)  
    file.write('\n#\n')
    file.write('# Outer shell processing statistics\n')
    file.write('#\n')
    file.write('_reflns_shell.number_measured_all ')
    file.write(datas_num)                  
    file.write('\n_reflns_shell.number_measured_obs ')
    file.write(datas_num) 
    file.write('\n_reflns_shell.d_res_low ')
    file.write(datas_rlow)
    file.write('\n_reflns_shell.d_res_high ')
    file.write(datas_rhigh)
    file.write('\n_reflns_shell.meanI_over_sigI_obs ')
    file.write(datas_ioversig)
    file.write('\n_reflns_shell.Rmerge_I_obs ')
    file.write(datas_rmerge)
    file.write('\n_reflns_shell.percent_possible_all ')
    file.write(datas_percentobs)
    file.write('\n_reflns_shell.pdbx_redundancy ')
    file.write(datas_redund)

    # Write refinement information

    file.write('\n#\n')
    file.write('###########################\n')
    file.write('# Refinement information  #\n')
    file.write('###########################\n')
    file.write('#\n')
    file.write('_refine.entry_id 1\n')

    if xray_method == 'MAD':
        file.write('_refine.pdbx_method_to_determine_struct  ')
        file.write(quote)
        file.write('MAD phasing')
        file.write(quote)

    if xray_method == 'SAD':
        file.write('_refine.pdbx_method_to_determine_struct  ')
        file.write(quote)
        file.write('SAD phasing')
        file.write(quote)

    if xray_method == 'IR':
        file.write('_refine.pdbx_method_to_determine_struct  ')
        file.write(quote)
        file.write('Molecular Replacement')
        file.write(quote)

    if xray_method != 'MAD' and xray_method != 'SAD' and xray_method != 'IR':
        file.write('_refine.pdbx_method_to_determine_struct  ')
        file.write(quote)
        file.write(xray_method)
        file.write(quote)

    file.write('\n#\n')
    file.write('# Data selection\n')
    file.write('#\n')
    file.write('_refine.ls_d_res_low ')
    file.write(ref_dlow)
    file.write('\n_refine.ls_d_res_high ')
    file.write(ref_dhigh)
    file.write('\n#\n')
    file.write('# Bulk solvent scattering model correction\n')
    file.write('#\n')
    file.write('_refine.solvent_model_details\n')

    if ref_bulksolvent == 'babinet' or ref_bulksolvent == 'fixedbabinet':
        file.write('; Babinet bulk solvent correction\n')
        file.write(';\n')
        file.write('_refine.solvent_model_param_ksol ')
        file.write(ref_ksolv)
        file.write('\n_refine.solvent_model_param_bsol ')
        file.write(ref_bsolv)
        file.write('\n')

    if ref_bulksolvent == 'mask':
        file.write('; Mask bulk solvent correction\n')
        file.write(';\n')
        file.write('_refine.pdbx_solvent_vdw_probe_radii ')
        file.write(ref_solvent_vdw_probe_radii)
        file.write('\n_refine.pdbx_solvent_ion_probe_radii ')
        file.write(ref_solvent_ion_probe_radii)
        file.write('\n_refine.pdbx_solvent_shrinkage_radii ')
        file.write(ref_solvent_shrinkage_radii)
        file.write('\n')

    file.write('#\n')
    file.write('# Refinement scaling\n')
    file.write('#\n')
    file.write('_refine.aniso_B[1][1] ')
    file.write(ref_b11)
    file.write('\n_refine.aniso_B[1][2] ')
    file.write(ref_b12)
    file.write('\n_refine.aniso_B[1][3] ')
    file.write(ref_b13)
    file.write('\n_refine.aniso_B[2][2] ')
    file.write(ref_b22)
    file.write('\n_refine.aniso_B[2][3] ')
    file.write(ref_b23)
    file.write('\n_refine.aniso_B[3][3] ')
    file.write(ref_b33)
    file.write('\n#\n')
    file.write('# Mean B-factor\n')
    file.write('#\n')
    file.write('_refine.B_iso_mean ')
    file.write(ref_bmean)

    file.write('\n#\n')
    file.write('# B-factor refinement method\n')
    file.write('#\n')

    if ref_anisoflag == 'yes':
        file.write('_refine.pdbx_isotropic_thermal_model ')
        file.write(quote)
        file.write('anisotropic')
        file.write(quote)
    else:
        file.write('_refine.pdbx_isotropic_thermal_model ')
        file.write(quote)
        file.write('isotropic')
        file.write(quote)

    file.write('\n#\n')
    file.write('# Overall R-factors\n')
    file.write('#\n')
    file.write('_refine.ls_number_reflns_all ')
    file.write(ref_numall)
    file.write('\n_refine.ls_number_reflns_obs ')
    file.write(ref_numobs)
    file.write('\n_refine.ls_number_reflns_R_free ')
    file.write(ref_numfree)
    file.write('\n_refine.ls_percent_reflns_obs ')
    file.write(ref_percent)
    file.write('\n_refine.ls_R_factor_all ')
    file.write(ref_rall)
    file.write('\n_refine.ls_R_factor_obs ')
    file.write(ref_robs)
    file.write('\n_refine.ls_R_factor_R_work ')
    file.write(ref_rwork)
    file.write('\n_refine.ls_R_factor_R_free ')
    file.write(ref_rfree)
    file.write('\n_refine.pdbx_ls_sigma_F 0.0\n')
    file.write('_refine.pdbx_ls_cross_valid_method ')
    file.write(quote)
    file.write('Free R-value')
    file.write(quote)
    file.write('\n_refine.pdbx_R_Free_selection_details ')
    file.write(quote)
    file.write('random')
    file.write(quote)
    file.write('\n_refine.pdbx_stereochemistry_target_values ')
    file.write(quote)
    file.write('Engh-Huber')
    file.write(quote)
    file.write('\n#\n')
    file.write('# Stereochemical agreement\n')
    file.write('#\n')
    file.write('loop_\n')
    file.write('_refine_ls_restr.type\n')
    file.write('_refine_ls_restr.dev_ideal\n')
    file.write('r_bond_d         ')
    file.write(ref_dbond)
    file.write('\nr_angle_d        ')
    file.write(ref_dangle)
    file.write('\nr_planar_tor     ')
    file.write(ref_dtorsion)
    file.write('\nr_chiral_restr   ')
    file.write(ref_dchiral)
    file.write('\nr_plane_restr    ')
    file.write(ref_dplane)
    file.write('\nr_mcbond_it      ')
    file.write(ref_bmbond)
    file.write('\nr_mcangle_it     ')
    file.write(ref_bmangle)
    file.write('\nr_scbond_it      ')
    file.write(ref_bsbond)
    file.write('\nr_scangle_it     ')
    file.write(ref_bsangle)            

    file.write('\n#\n')
    file.write('# Atom counts\n')
    file.write('#\n')
    file.write('_refine_hist.cycle_id 1\n')
    file.write('_refine_hist.d_res_high ')
    file.write(ref_dhigh)
    file.write('\n_refine_hist.d_res_low ')
    file.write(ref_dlow)
    file.write('\n_refine_hist.number_atoms_total ')
    file.write(ref_natom)
    file.write('\n_refine_hist.number_atoms_solvent ')
    file.write(ref_nsolvent)
    file.write('\n#\n')

    # Add coordinates. Note that this uses the REFMAC o/p for better atom-typing

    xyzfile = open('temp_ref.pdb', 'r')
    allLines = xyzfile.readlines()
    xyzfile.close()

    os.remove('temp_use.pdb')
    os.remove('temp_ref.pdb')

    file.write('##############################\n')
    file.write('# Coordinates                #\n')
    file.write('##############################\n')
    file.write('#                             \n')
    file.write('loop_                         \n')
    file.write('_atom_site.type_symbol        \n')
    file.write('_atom_site.label_atom_id      \n')
    file.write('_atom_site.label_comp_id      \n')
    file.write('_atom_site.auth_asym_id       \n')
    file.write('_atom_site.auth_seq_id        \n')
    file.write('_atom_site.label_seq_id       \n')
    file.write('_atom_site.label_alt_id       \n')
    file.write('_atom_site.Cartn_x            \n')
    file.write('_atom_site.Cartn_y            \n')
    file.write('_atom_site.Cartn_z            \n')
    file.write('_atom_site.occupancy          \n')
    file.write('_atom_site.B_iso_or_equiv     \n')
    file.write('_atom_site.footnote_id        \n')
    file.write('_atom_site.label_entity_id    \n')
    file.write('_atom_site.id                 \n')
    file.write('_atom_site.label_asym_id      \n')

    for eachLine in allLines:

        tag = eachLine[0:6]
        tag = tag.strip()

        # Get deorthogonalization matrix

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

        # Get coordinates from ATOM/HETATM records

        if tag == 'ATOM' or tag == 'HETATM':

            atom_serial = eachLine[6:11]
            atom_name = eachLine[12:16]
            atom_alt = eachLine[16:17]
            res_name = eachLine[17:20]
            chain_id = eachLine[21:22]
            res_number = eachLine[22:26]
            insert_code = eachLine[26:27]
            x_coord = eachLine[30:38]
            y_coord = eachLine[38:46]
            z_coord = eachLine[46:54]
            occ_value = eachLine[54:60]
            b_value = eachLine[60:66]
            element = eachLine[76:78]

            # strip or pad some records

            pr_res_name = res_name.strip()

            wr_atom_serial = ' ' + atom_serial + ' '
            wr_atom_name = ' ' + atom_name + ' '
            wr_atom_alt = ' ' + atom_alt + ' '
            wr_res_name = ' ' + res_name + ' '
            wr_chain_id = ' ' + chain_id + ' '
            wr_res_number = ' ' + res_number + ' '
            wr_insert_code = ' ' + insert_code + ' '
            wr_x_coord = ' ' + x_coord + ' '
            wr_y_coord = ' ' + y_coord + ' '
            wr_z_coord = ' ' + z_coord + ' '
            wr_occ_value = ' ' + occ_value + ' '
            wr_b_value = ' ' + b_value + ' '
            wr_element = ' ' + element + ' '

            # Patch null alternate records

            if atom_alt == ' ':
                wr_atom_alt = ' . '

            # Establish label_asym records

            if res_name != 'HOH':

                # Initially as default protein record

                label_asym = '   ' + chain_id
                label_entity = '1'

                # Reset as a ligand record

                count = 0
                while count < number_entities:
                    hetname = aList_hets[count]
                    if hetname == res_name:
                        label_asym = aList_hets_asym[count]
                        label_entity = aList_hets_number[count]

                    count = count + 1

            else:

                label_asym =  '   W'
                label_entity = water_entity_list_number

            #

            wr_label_asym = ' ' + label_asym
            wr_label_entity = ' ' + label_entity + ' '

            # Write CIF atom record

            file.write('  ')
            file.write(wr_element)
            file.write(wr_atom_name)
            file.write(wr_res_name)
            file.write(wr_chain_id)
            file.write(wr_res_number)
            file.write(wr_res_number)
            file.write(wr_atom_alt)
            file.write(wr_x_coord)
            file.write(wr_y_coord)
            file.write(wr_z_coord)
            file.write(wr_occ_value)
            file.write(wr_b_value)
            file.write(footnote)
            file.write(wr_label_entity)
            file.write(wr_atom_serial)
            file.write(label_asym)
            file.write('\n')

    # Add list of anisotropic records of same type and order as PDB file

    if ref_anisoflag == 'yes':

        file.write('################################\n')
        file.write('# Anisotropic B-factor records #\n')
        file.write('################################\n')
        file.write('loop_\n')
        file.write('_atom_site_anisotrop.id\n')
        file.write('_atom_site_anisotrop.type_symbol\n')
        file.write('_atom_site_anisotrop.U[1][1]\n')
        file.write('_atom_site_anisotrop.U[2][2]\n')
        file.write('_atom_site_anisotrop.U[3][3]\n')
        file.write('_atom_site_anisotrop.U[1][2]\n')
        file.write('_atom_site_anisotrop.U[1][3]\n')
        file.write('_atom_site_anisotrop.U[2][3]\n')

        for eachLine in allLines:

            tag = eachLine[0:6]
            tag = tag.strip()

            if tag == 'ANISOU':

                atom_serial = eachLine[6:11]           
                element = eachLine[76:78]
                u11 = eachLine[28:35]
                u22 = eachLine[35:42]
                u33 = eachLine[42:49]
                u12 = eachLine[49:56]
                u13 = eachLine[56:63]
                u23 = eachLine[63:70]

                wr_atom_serial = ' ' + atom_serial + ' '
                wr_element = ' ' + element + ' '
                wr_u11 = ' ' + u11 + ' '
                wr_u22 = ' ' + u22 + ' '
                wr_u33 = ' ' + u33 + ' '
                wr_u12 = ' ' + u12 + ' '
                wr_u13 = ' ' + u13 + ' '
                wr_u23 = ' ' + u23 + ' '

                file.write(wr_atom_serial)
                file.write(wr_element)
                file.write(wr_u11)
                file.write(wr_u22)
                file.write(wr_u33)
                file.write(wr_u12)
                file.write(wr_u13)
                file.write(wr_u23)
                file.write('\n')
    #

    file.write('#\n')
    file.write('###########\n')
    file.write('# The End #\n')
    file.write('###########\n')
    file.write('#\n')

    file.close()

    ###################################################################
    # Write simple CIF reflection file list from the refinement data  #
    ###################################################################

    if write_hkl != 'no':

        hklfile = rootname + '_hkl.cif'

        print 'Creating CIF reflection list'

        filename = 'mi_runmtz2various.inp'
        file = open(filename, 'w')
        file.write('LABIN FP=')
        file.write(famp)
        file.write(' SIGFP=')
        file.write(sd)
        file.write(' FREE=')
        file.write(freer)
        file.write('\nOUTPUT CIF data_1\n')
        file.write('FREEVAL ')
        file.write(nfree_exclude)
        file.write('\nMISS\n')
        file.write('SCALE 10\n')
        file.write('MONITOR 1000\n')
        file.write('END\n')
        file.close()

        runmtz2various = 'mtz2various HKLIN ' + local_mtzfile + ' HKLOUT ' + hklfile + ' < mi_runmtz2various.inp > mi_mtz2various.out'
        os.system(runmtz2various)

        fileexists = os.path.exists(hklfile)
        if fileexists == 0:
            print 'MTZ2VARIOUS run to generate CIF reflection file appears to have failed\n'
            time.sleep(4)
            return 1
        
        fileexists = os.path.exists('mi_mtz2various.out')
        if fileexists != 0:
            os.remove('mi_mtz2various.out')

        fileexists = os.path.exists('mi_runmtz2various.inp')
        if fileexists != 0:            
            os.remove('mi_runmtz2various.inp')

    # Clean-up

    fileexists = os.path.exists(local_mtzfile)
    if fileexists != 0:
        os.remove(local_mtzfile)

    fileexists = os.path.exists(liganddir)
    if fileexists != 0:
        os.remove(liganddir)

    fileexists = os.path.exists('mi_refmac.out')
    if fileexists != 0:
        os.remove('mi_refmac.out')

    fileexists = os.path.exists('output.refmac')
    if fileexists != 0:
        os.remove('output.refmac')

    #########################
    # Report statistics     #
    #########################

    print '\nRefinement Summary'
    print '==================='

    print 'Resolution            :',ref_dhigh
    print 'R(working)            :',ref_rwork
    print 'R(free)               :',ref_rfree
    print 'RMSD(bonds)           :',ref_dbond
    print 'RMSD(angles)          :',ref_dangle
    print 'Outlier phi-psi angles:',ref_number_phi_psi_errors

    ##########################################
    # Write PDB (cif file) deposition notes  #
    ##########################################

    if write_cif != 'no':

        print '\nPDB Deposition Notes'
        print '====================='
        print 'Structure deposition file : ',cifdepositfile

        if write_hkl != 'no':
            print 'X-ray data deposition file: ',hklfile

        print '\nThe mmCIF sequence/structure/entity mappings will need adjustment'
        print 'if there is more than one type of protein in the crystal.'

        print '\nThe mmCIF structure deposition file contains annotation data and coordinates.'
        print 'This file may be deposited to the PDB through the RCSB/PDB ADIT interface'
        print '(http://rcsb-deposit.rutgers.edu/adit/). From the ADIT session ' 
        print 'select file type "mmCIF" and upload the file. After deposition, use the'
        print '"PREVIEW ENTRY" option to check the information uploaded to the RCSB/PDB.'
        print 'Any alterations, missing or additional information may be entered through'
        print 'this interface.\n'

    else:

        os.remove(cifdepositfile)

    ##############################
    # Write annotation text file #
    ##############################

    if write_text != 'no':

        textfile = rootname + '.rtf'

        print 'Creating text file: ',textfile

        celllengthLine = acell + ' ' + bcell + ' ' + ccell
        cellangleLine = alpha + ' ' + beta + ' ' + gamma
        residue_count = str(residue_count)
        ref_dhigh = float(ref_dhigh)
        ref_dhigh = round(ref_dhigh,2)
        ref_dhigh = str(ref_dhigh)
        ref_rwork = float(ref_rwork)
        ref_rwork = round(ref_rwork,3)
        ref_rwork = str(ref_rwork)
        ref_rfree = float(ref_rfree)
        ref_rfree = round(ref_rfree,3)
        ref_rfree = str(ref_rfree)

        file = open(textfile,'w')
        file.write('Resolution: ')
        file.write(ref_dhigh)
        file.write('\nUnit Cell: ')
        file.write(celllengthLine)
        file.write('\n           ')
        file.write(cellangleLine)
        file.write('\nSpace Group: ')
        file.write(spgname)
        file.write('\nR(work): ')
        file.write(ref_rwork)
        file.write('\nR(free): ')
        file.write(ref_rfree)
        file.write('\nPolymer Chain IDs: ')

        count = 0
        while count < number_chains:
            chain_out = aList_chains[count]
            pr_chain_out = ' ' + chain_out + ' '
            file.write(pr_chain_out)

            if count != number_chains - 1:
                file.write(',')

            count = count + 1

        file.write('\nResidues: ')
        file.write(residue_count)
        file.write('\n')

        file.write('Components: ')

        count = 0
        if number_entities > 0:

            while count < number_entities:

                entity_list_asym = aList_hets_asym[count]
                entity_list_asym = entity_list_asym.replace('_W','')
                pr_entity_list_asym = ' ( ' + entity_list_asym + ')'

                entity_list_name = aList_hets_names[count]
                entity_list_name = entity_list_name.lower()

                file.write(entity_list_name)
                file.write(pr_entity_list_asym)
                file.write('            \n')

                count = count + 1

        file.close()

    ##############################
    # Write annotation html file #
    ##############################

    if write_html != 'no':

        htmlfile = rootname + '.htm'

        # Bring any images local for better browser/acrobat compatibility

        if write_image1 != 'none':

            local_image1_base = os.path.basename(write_image1)
            local_image1 = rootname + '_' + local_image1_base

            file = open(write_image1,'rb')
            allLines = file.readlines()
            file.close()

            file = open(local_image1,'wb')
            file.writelines(allLines)
            file.close()

        if write_image2 != 'none':

            local_image2_base = os.path.basename(write_image2)
            local_image2 = rootname + '_' + local_image2_base

            file = open(write_image2,'rb')
            allLines = file.readlines()
            file.close()

            file = open(local_image2,'wb')
            file.writelines(allLines)
            file.close()

        if write_image3 != 'none':

            local_image3_base = os.path.basename(write_image3)
            local_image3 = rootname + '_' + local_image3_base

            file = open(write_image3,'rb')
            allLines = file.readlines()
            file.close()

            file = open(local_image3,'wb')
            file.writelines(allLines)
            file.close()

        # Write html

        print 'Creating HTML file: ',htmlfile

        celllengthline = acell + ' ' + bcell + ' ' + ccell + ' ' + alpha + ' ' + beta + ' ' + gamma
        water_count = str(water_count)
        ligand_count = str(ligand_count)
        pr_number_chains = str(number_chains)
        residue_count = str(residue_count)

        file = open(htmlfile,'w')
        file.write('<html>\n')
        file.write('<head><title>structure report</title></head>\n')
        file.write('<body bgcolor = ')
        file.write(quote)
        file.write('white')
        file.write(quote)
        file.write('>\n')

        file.write('<h2><center>')
        if write_title == 'none':
            file.write('Structure determination and model refinement statistics')
        else:
            file.write(write_title)
        file.write('</center></h2>\n')

        # Crystal information

        if datalogfile != 'none':
            file.write('<p>\n')
            file.write('<b>Crystal characteristics and data collection statistics (outer shell statistics in parenthesis)</b>\n')
            file.write('<p>\n')
        else:
            file.write('<p>\n')
            file.write('<b>Crystal characteristics</b>\n')
            file.write('<p>\n')

        file.write('<table border=0>\n')
        file.write('<tr><td><i>Unit cell (&#197, &#186)</i></td><td>')
        file.write(celllengthline)
        file.write('</td></tr>\n')
        file.write('<tr><td><i>Space group</i></td><td>')
        file.write(spgname_standard)
        file.write('</td></tr>\n')

        # Data statistics

        if datalogfile != 'none':

            file.write('<tr><td><i>Resolution range (&#197)</i></td><td>')
            file.write(data_rlow)
            file.write(' - ')
            file.write(data_rhigh)
            file.write(' (')
            file.write(datas_rlow)
            file.write(' - ')
            file.write(datas_rhigh)
            file.write(')')
            file.write('</td></tr>\n')
            file.write('<tr><td><i>No. of observations</i></td><td>')
            file.write(data_num_unmerged)
            file.write('</td></tr>\n')
            file.write('<tr><td><i>No. of unique reflections</i></td><td>')
            file.write(data_num)
            file.write('</td></tr>\n')
            file.write('<tr><td><i>Redundancy</i></td><td>')
            file.write(data_redund)
            file.write(' (')
            file.write(datas_redund)
            file.write(')')
            file.write('</td></tr>\n')   
            file.write('<tr><td><i>Completeness (%)</i></td><td>')
            file.write(data_percentobs)
            file.write(' (')
            file.write(datas_percentobs)
            file.write(')')
            file.write('</td></tr>\n')
            file.write('<tr><td><i>Mean I/sigma(I)</i></td><td>')
            file.write(data_ioversig)
            file.write(' (')
            file.write(datas_ioversig)
            file.write(')')
            file.write('</td></tr>\n')
            file.write('<tr><td><i>R<sub>merge</sub></i></td><td>')
            file.write(data_rmerge)
            file.write(' (')
            file.write(datas_rmerge)
            file.write(')')
            file.write('</td></tr>\n')

        file.write('</table>\n')
        file.write('<p>\n')

        # Model related

        file.write('<b>Crystallographic data and refinement statistics</b>\n')
        file.write('<p>\n')
        file.write('<table border=0>\n')

        file.write('<tr><td><i>Resolution range (&#197)</i></td><td>')
        file.write(ref_dlow)
        file.write(' - ')
        file.write(ref_dhigh)
        file.write('</td></tr>\n')
        file.write('<tr><td><i>No. of reflections</i></td><td>')
        file.write(ref_numall)
        file.write(' (')
        file.write(ref_numwork)
        file.write(' working set, \n')
        file.write(ref_numfree)
        file.write(' test set)</td></tr>\n')

        # Chains

        file.write('<tr><td><i>No. of protein chains</i></td><td>')
        file.write(pr_number_chains)
        file.write(' (')

        count = 0
        while count < number_chains:
            chain_out = aList_chains[count]
            pr_chain_out = ' ' + chain_out + ' '
            file.write(pr_chain_out)

            if count != number_chains - 1:
                file.write(',')

            count = count + 1    

        file.write(')')
        file.write('</td></tr>\n')

        # Ligand ids

        file.write('<tr><td><i>Ligand id codes</i></td><td>')

        count = 0
        if number_entities > 0:

            while count < number_entities:

                entity_list_asym = aList_hets_asym[count]
                entity_list_asym = entity_list_asym.replace('_W','')
                file.write(entity_list_asym)

                if count != number_entities - 1:
                    file.write(',')

                count = count + 1

        else:

            file.write('-\n')

        file.write('</td></tr>\n')

        # Summary counts and statistics

        file.write('<tr><td><i>No. of protein residues</i></td><td>')
        file.write(residue_count)
        file.write('</td></tr>\n')
        file.write('<tr><td><i>No. of ligands</i></td><td>')
        file.write(ligand_count)
        file.write('</td></tr>\n')
        file.write('<tr><td><i>No. of waters</i></td><td>')
        file.write(water_count)
        file.write('</td></tr>\n')
        file.write('<tr><td><i>R<sub>work</sub></i></td><td>')
        file.write(ref_rwork)
        file.write('</td></tr>\n')
        file.write('<tr><td><i>R<sub>free</sub></i></td><td>')
        file.write(ref_rfree)
        file.write('</td></tr>\n')
        file.write('<tr><td><i>Rmsd bond lengths (&#197)</i></td><td>')
        file.write(ref_dbond)
        file.write('</td></tr>\n')
        file.write('<tr><td><i>Rmsd bond angles (&#186)</i></td><td>')
        file.write(ref_dangle)
        file.write('</td></tr>\n')
        file.write('<tr><td><i>Number of disallowed &#966&#968 angles</i></td><td>')
        file.write(ref_number_phi_psi_errors)
        file.write('</td></tr>\n')
        file.write('</table>\n')
        file.write('<p>\n')

        # Optional addition of molecular images

        if write_image3 != 'none':
            image_size = 'width="210" height="210" border="1"'
        else:
            image_size = 'width="250" height="250" border="1"'

        if write_image1 != 'none':

            file.write('<b>Molecular images</b>\n')
            file.write('<p>\n')
            file.write('<img src="')
            file.write(local_image1)
            file.write('" align="left" ')
            file.write(image_size)
            file.write('>\n')

            if write_image2 != 'none':

                file.write('<img src="')
                file.write(local_image2)
                file.write('" align="center" ')
                file.write(image_size)
                file.write('>\n')            

                if write_image3 != 'none':

                    file.write('<img src="')
                    file.write(local_image3)
                    file.write('" align="center" ')
                    file.write(image_size)
                    file.write('>\n')

        file.write('</body>\n')
        file.write('</html>\n')

        file.close()

    ###################
    # Create CCP4 map #
    ###################

    if write_map == 'yes':

        print '\nCalculating CCP4 map\n'

        # Compute ML-weighted map

        file = open('mi_fft.inp','w')
        file.write('LABIN F1=FWT PHI=PHWT\n')
        file.write('END\n')
        file.close()

        runfft = 'fft HKLIN temp_ref.mtz MAPOUT mi_2ff.map < mi_fft.inp > mi_fft.log'
        os.system(runfft)

        fileexists = os.path.exists('mi_2ff.map')
        if fileexists == 0:
            print 'FFT for map display failed'
            time.sleep(4)
            return 1
        else:
            os.remove('mi_fft.inp')
            os.remove('mi_fft.log')

        # Expand to P1 cell

        mapfilename = rootname + '_fullcell.map'

        file = open('mi_mapmask.inp','w')
        file.write('XYZLIM 0 1 0 1 0 1\n')
        file.write('END\n')
        file.close()

        runmapmask = 'mapmask MAPIN mi_2ff.map MAPOUT ' + mapfilename + ' < mi_mapmask.inp > mi_mapmask.log'
        os.system(runmapmask)

        fileexists = os.path.exists(mapfilename)
        if fileexists == 0:
            print 'MAPMASK for expansion to full cell'
            time.sleep(4)
            return 1
        else:
            os.remove('mi_mapmask.inp')
            os.remove('mi_mapmask.log')

        print 'Writing CCP4 for unit cell:',mapfilename,'\n'

    #############################################
    # Write CCP4 map pieces around the ligands  #
    #############################################

    if write_map == 'yes' and number_entities > 0:

        print '\nWriting CCP4 maps around ligands\n'

        aList_ligand_chain = []
        aList_ligand_res_number = []
        aList_ligand_atom_count = []
        aList_ligand_res_name = []

        aList_ligand_xmin_x = []
        aList_ligand_xmin_y = []
        aList_ligand_xmin_z = []    
        aList_ligand_ymin_x = []
        aList_ligand_ymin_y = []
        aList_ligand_ymin_z = []    
        aList_ligand_zmin_x = []
        aList_ligand_zmin_y = []
        aList_ligand_zmin_z = []    
        aList_ligand_xmax_x = []
        aList_ligand_xmax_y = []
        aList_ligand_xmax_z = []    
        aList_ligand_ymax_x = []
        aList_ligand_ymax_y = []
        aList_ligand_ymax_z = []   
        aList_ligand_zmax_x = []
        aList_ligand_zmax_y = []
        aList_ligand_zmax_z = []
        aList_ligand_limits = []

        chain_id_prev = '?'
        res_number_prev = '?'


        file = open(local_pdbfile,'r')
        allLines = file.readlines()
        file.close()

        # Get ID information for all individual ligands

        for eachLine in allLines:

            tag = eachLine[0:6]
            tag = tag.strip()            

            if tag == 'ATOM' or tag == 'HETATM':

                chain_id = eachLine[21:22]
                chain_id = chain_id.strip()
                res_number = eachLine[22:26]
                res_number = res_number.strip()
                res_name = eachLine[17:20]
                res_name = res_name.strip()

                count = 0
                while count < number_entities:
                    entity = aList_hets[count]

                    if res_name == entity:    

                        if chain_id != chain_id_prev or res_number != res_number_prev:

                            aList_ligand_chain.append(chain_id)
                            aList_ligand_res_number.append(res_number)
                            aList_ligand_res_name.append(res_name)
                            aList_ligand_atom_count.append(0)

                            aList_ligand_xmin_x.append(999.)
                            aList_ligand_xmin_y.append(999.)
                            aList_ligand_xmin_z.append(999.)
                            aList_ligand_ymin_x.append(999.)
                            aList_ligand_ymin_y.append(999.)
                            aList_ligand_ymin_z.append(999.)
                            aList_ligand_zmin_x.append(999.)
                            aList_ligand_zmin_y.append(999.)
                            aList_ligand_zmin_z.append(999.)
                            aList_ligand_xmax_x.append(-999.)
                            aList_ligand_xmax_y.append(-999.)
                            aList_ligand_xmax_z.append(-999.)
                            aList_ligand_ymax_x.append(-999.)
                            aList_ligand_ymax_y.append(-999.)
                            aList_ligand_ymax_z.append(-999.)
                            aList_ligand_zmax_x.append(-999.)
                            aList_ligand_zmax_y.append(-999.)
                            aList_ligand_zmax_z.append(-999.)

                    chain_id_prev = chain_id
                    res_number_prev = res_number

                    count = count + 1

        # Get coordinate limit data for each ligand

        number_pdb_ligands = len(aList_ligand_chain)

        for eachLine in allLines:

            tag = eachLine[0:6]
            tag = tag.strip()

            if tag == 'ATOM' or tag == 'HETATM':

                chain_id = eachLine[21:22]
                chain_id = chain_id.strip()
                res_number = eachLine[22:26]
                res_number = res_number.strip()

                count = 0
                while count < number_pdb_ligands:

                    chain_entity = aList_ligand_chain[count]
                    res_number_entity = aList_ligand_res_number[count]

                    if chain_id == chain_entity and res_number == res_number_entity:

                        atom_count = aList_ligand_atom_count[count] + 1
                        aList_ligand_atom_count[count] = atom_count

                        x = eachLine[30:38]
                        y = eachLine[38:46]
                        z = eachLine[46:54]
                        x = float(x)
                        y = float(y)
                        z = float(z)
                        xmin = aList_ligand_xmin_x[count]
                        ymin = aList_ligand_ymin_y[count]
                        zmin = aList_ligand_zmin_z[count]
                        xmax = aList_ligand_xmax_x[count]
                        ymax = aList_ligand_ymax_y[count]
                        zmax = aList_ligand_zmax_z[count]

                        if x < xmin:
                            aList_ligand_xmin_x[count] = x
                            aList_ligand_xmin_y[count] = y
                            aList_ligand_xmin_z[count] = z

                        if y < ymin:
                            aList_ligand_ymin_x[count] = x
                            aList_ligand_ymin_y[count] = y
                            aList_ligand_ymin_z[count] = z

                        if z < zmin:
                            aList_ligand_zmin_x[count] = x
                            aList_ligand_zmin_y[count] = y
                            aList_ligand_zmin_z[count] = z

                        if x > xmax:
                            aList_ligand_xmax_x[count] = x
                            aList_ligand_xmax_y[count] = y
                            aList_ligand_xmax_z[count] = z

                        if y > ymax:
                            aList_ligand_ymax_x[count] = x
                            aList_ligand_ymax_y[count] = y
                            aList_ligand_ymax_z[count] = z

                        if z > zmax:
                            aList_ligand_zmax_x[count] = x
                            aList_ligand_zmax_y[count] = y
                            aList_ligand_zmax_z[count] = z

                    count = count + 1

        # Obtain coordinate limits with border

        count = 0
        while count < number_pdb_ligands:

            xmin_x = aList_ligand_xmin_x[count] - map_border
            xmin_y = aList_ligand_xmin_y[count]
            xmin_z = aList_ligand_xmin_z[count]
            ymin_x = aList_ligand_ymin_x[count]
            ymin_y = aList_ligand_ymin_y[count] - map_border
            ymin_z = aList_ligand_ymin_z[count]
            zmin_x = aList_ligand_zmin_x[count]
            zmin_y = aList_ligand_zmin_y[count] 
            zmin_z = aList_ligand_zmin_z[count] - map_border
            xmax_x = aList_ligand_xmax_x[count] + map_border
            xmax_y = aList_ligand_xmax_y[count]
            xmax_z = aList_ligand_xmax_z[count]
            ymax_x = aList_ligand_ymax_x[count]
            ymax_y = aList_ligand_ymax_y[count] + map_border
            ymax_z = aList_ligand_ymax_z[count]
            zmax_x = aList_ligand_zmax_x[count]
            zmax_y = aList_ligand_zmax_y[count]
            zmax_z = aList_ligand_zmax_z[count] + map_border

            # Convert key points to fractional coordinates

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
            aList_ligand_limits.append(xyz_limits)

            count = count + 1

        # Build density maps with CCP4/MAPMASK around each ligand

        count = 0
        while count < number_pdb_ligands:

            atom_count = aList_ligand_atom_count[count]
            chain_id = aList_ligand_chain[count]
            res_number = aList_ligand_res_number[count]
            res_name = aList_ligand_res_name[count]
            xyz_limits = aList_ligand_limits[count]

            if atom_count > 6:

                mapfilename = rootname + '_' + chain_id + res_number + '_' + res_name + '.map'
                print 'Writing',mapfilename,' with limits in fractional coords',xyz_limits

                file = open('mi_mapmask.inp','w')
                file.write('XYZLIM ')
                file.write(xyz_limits)
                file.write('\n')
                file.write('EXTEND XTAL\n')
                file.write('END\n')
                file.close()

                runmapmask = 'mapmask MAPIN mi_2ff.map MAPOUT ' + mapfilename + ' < mi_mapmask.inp > mi_mapmask.log'
                os.system(runmapmask)

                fileexists = os.path.exists(mapfilename)
                if fileexists == 0:
                    print 'MAPMASK for ligand maps failed'
                    time.sleep(4)
                    return 1
                else:
                    os.remove('mi_mapmask.inp')
                    os.remove('mi_mapmask.log')

            count = count + 1

        file.close()        


    fileexists = os.path.exists('mi_2ff.map')
    if fileexists != 0:
        os.remove('mi_2ff.map')

    ######################
    # Append project log #
    ######################

    print '\nWriting project log'

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
    file.write('\nInput sequence: ')
    file.write(seqfile)
    file.write('\nInput data log: ')
    file.write(datalogfile)
    file.write('\nInput template: ')
    file.write(templatefile)
    file.write('\n')

    if write_cif != 'no':
        file.write('Output mmCIF deposition file: ')
        file.write(cifdepositfile)
        file.write('\n')
    else:
        file.write('Output mmCIF deposition file: none\n')

    if write_hkl != 'no':
        file.write('Output data deposition file: ')
        file.write(hklfile)
        file.write('\n')
    else:
        file.write('Output data deposition file: none\n')

    if write_text != 'no':
        file.write('Output report text summary: ')
        file.write(textfile)
        file.write('\n')
    else:
        file.write('Output report text summary: none\n')

    if write_html != 'no':
        file.write('Output report html summary: ')
        file.write(htmlfile)
        file.write('\n')
    else:
        file.write('Output report html summary: none\n')

    if write_map != 'no':
        file.write('Output ligand density maps: yes\n')
    else:
        file.write('Output ligand density maps: none\n')

    file.write('---------------\n')
    file.close()

    # Clean-up job-specific temporary CCP4_SCR space and put back original CCP4_SCR

    fileexists = os.path.exists(local_pdbfile)
    if fileexists != 0:  
        os.remove(local_pdbfile)

    fileexists = os.path.exists('temp_ref.mtz')
    if fileexists != 0:     
        os.remove('temp_ref.mtz')

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

    os.environ['CCP4_SCR'] = ccp4.scr

    time.sleep(4)

    #
    return 0

if __name__ == "__main__":
    sys.exit(Run())
