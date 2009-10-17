#####################################################################
#                                                                   #
# Molecular Replacement Script using CCP4/MOLREP or CCP4/PHASER     #
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
import getopt
import ccp4check

def Usage():
    print "Usage: %s [options]" % sys.argv[0]
    print "Options are:"
    print "  -p,--pdbfile=FILE             PDB file to use"
    print "  -m,--mtzfile=FILE             The diffraction data file to use (MTZ, D*TREK or SCALEPACK formats)"
    print "  -f,--fixed_pdb=FILE           The fixed PDB model. default: no fixed file"
    print "  -d,--workdir=DIR              Working dir"
    print "  -a,--multi_search=yes or no   Search over multiple models? Default: no"
    print "  -s,--match_pdbin=yes or no    Symmetry match coordinates to input position: Default no"
    print "  -e,--engine=phaser or molrep  Which MR engine to use. Default: molrep"
    print "  -t,--sg_search=yes or no      Run a spacegroup search. Default no"
    print "  -n,--spacegroup=NUM           Set the spacegroup number. Default: from mtz file"
    print "  -c,--copies=NUM               Set the number of protein copies to search for. Default: auto"
    print "  -h,--help                     This help"

# Initialize
def Run(argv=None):
    if argv is None:
        argv=sys.argv

    quote = """'"""

    pdbfile = 'none'
    mtzfile = 'none'
    pdbfile2 = 'none'
    workingdir = 'none'
    runid = '1'
    runid_int = 0
    projectlog = 'project_history.txt'
    multi_search = 'no'
    match_pdbin = 'no'
    mr_r = 'none'
    ilabel = 'none'
    flabel = 'none'
    sigflabel = 'none'
    seq_insertion_flag = 'no' 
    number_models = 1
    number_molecules = '1'
    mr_method = 'molrep'
    sg_search = 'no'
    mr_spacegroup_no = 'none'
    spacegroup_user_input = 'none'
    spacegroup_name = 'none'
    number_molecules_input = 'auto'
    
    datatype = 'none'
    wavelength = 'none'
    spacegroup_no = 'none'
    filename_local_in = 'mi_data.hkl'
    filename_mtz = 'mi_data.mtz'
    filename_out_mtz = 'mi_data_sorted.mtz'
    filename_out_mtz_f_only = 'mi_data_sorted_f_only.mtz'

    acell_mtz = 'none'
    bcell_mtz = 'none'
    ccell_mtz = 'none'
    alpha_mtz = 'none'
    beta_mtz = 'none'
    gamma_mtz = 'none'
    space_group_mtz = 'none'

    aList_targetPDB = []
    aList_multi_search = []
    aList_spacegroups = []
    aList_true_seq_num = []
    parseLine = []
    labelList = []
    colList = []

    job_prefix = 'molrep_'

    # Space group point group lists

    sg_class_1 = ['1']
    sg_class_2 = ['3','4','5']
    sg_class_222 = ['16','17','18','19','20','21','22','23','24']
    sg_class_4 = ['75','76','77','78','79','80']
    sg_class_422 = ['89','90','91','95','92','96','93','94','97','98']
    sg_class_3 = ['143','144','145','146']
    sg_class_32 = ['149','150','152','154','151','153','145']
    sg_class_6 = ['168','169','171','173','172','170']
    sg_class_622 = ['177','178','179','180','181','182']
    sg_class_23 = ['195','196','197','198','199']
    sg_class_432 = ['207','213','212','208','209','210','211','214']

    ######################################################
    # Read command-line #
    ######################################################

    # Option for short option command-line input
    number_of_args = len(argv)
    if number_of_args > 4:
        args = argv[1:]
        optlist, args = getopt.getopt(
            args,'p:m:f:d:a:s:e:t:n:c:?',
            ['pdbfile=','mtzfile=','fixed_pdb=','workdir=','multi_search=',
             'match_pdbin=','engine=','sg_search=','spacegroup=','copies=','help'])
        number_of_inputs = len(optlist)
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
                    pdbfile = param_value
                elif arg_value == '-m' or arg_value=='--mtzfile':
                    mtzfile = param_value
                elif arg_value == '-f' or arg_value=='--fixed_pdb':
                    pdbfile2 = param_value
                elif arg_value == '-d' or arg_value=='--workdir':
                    workingdir = param_value
                elif arg_value == '-a' or arg_value=='--multi_search':
                    multi_search = param_value
                elif arg_value == '-s' or arg_value=='--match_pdbin':
                    match_pdbin = param_value
                elif arg_value == '-e' or arg_value=='--engine':
                    mr_method = param_value
                elif arg_value == '-t' or arg_value=='--sg_search':
                    sg_search = param_value
                elif arg_value == '-n' or arg_value=='--spacegroup':
                    spacegroup_user_input = param_value
                elif arg_value == '-c' or arg_value=='--copies':
                    number_molecules_input = param_value
            count = count + 1
    else:
        Usage()
        return

    ccp4,error = ccp4check.ccp4check()
    if not ccp4:
      print '\n' + error + '\n'
      time.sleep(4)
      return 1

    # Check input
    fileexists = os.path.exists(pdbfile)
    if fileexists == 0:
        print 'The PDB file was not found ',pdbfile
        time.sleep(4)
        return 1

    fileexists = os.path.exists(mtzfile)
    if fileexists == 0:
        print 'The diffraction data file was not found ',mtzfile
        time.sleep(4)
        return 1

    fileexists = os.path.exists(workingdir)
    if fileexists == 0:
        print 'The working directory was not found ',workingdir
        time.sleep(4)
        return 1

    fileexists = os.path.exists(pdbfile2)
    if fileexists == 0 and pdbfile2 != 'none':
        print 'The fixed PDB model was not found ',pdbfile2
        time.sleep(4)
        return 1               

    if mr_method != 'phaser' and mr_method != 'PHASER':
        mr_method = 'molrep'
    else:
        mr_method = 'phaser'

    # Space group search 

    if sg_search != 'yes' and sg_search != 'YES':
        sg_search = 'no'
    else:
        sg_search = 'yes'
        spacegroup_user_input = 'none'

    # Get file list for multi-search option

    if multi_search == 'yes' or multi_search == 'YES':

        multi_search_dir = os.path.dirname(pdbfile)
        aList_multi_search = os.listdir(multi_search_dir)
        length = len(aList_multi_search)

        print 'Searching over multiple models from:',multi_search_dir

        count = 0
        while count < length:

            target_file = aList_multi_search[count]

            if target_file.find('.pdb') > -1:
                target_file_out = os.path.join(multi_search_dir,target_file)

                # Check it has enough atoms to be a protein model

                file = open(target_file_out,'r')
                allLines = file.readlines()
                file.close()

                count_atoms = 0

                for eachLine in allLines:
                    tag = eachLine[0:6]
                    tag = tag.strip()

                    if tag == 'ATOM' or tag == 'HETATM':
                        count_atoms = count_atoms + 1

                if count_atoms > 150:
                    aList_targetPDB.append(target_file_out)

            count = count + 1

        number_models = len(aList_targetPDB)
        print 'Number of models to search:',number_models,'\n'

    else:

        aList_targetPDB.append(pdbfile)

    ###########################################
    # Setup directories with required files   #
    ###########################################

    # Go to working directory

    os.chdir(workingdir)

    # Establish prefix to force local scratch files using a time string to almost uniquely identify them

    idcode = '000000'

    gmt = time.gmtime(time.time())
    fmt = '%H%M%S'
    idcode = time.strftime(fmt,gmt)

    path_scratch = 'temp_' + idcode + '_'

    #######################################################
    # Evaluate data file type (MTZ, D*TREK or SCALEPACK)  #
    #######################################################

    # Devise file name for case of integrated data

    base_mtzfile = os.path.basename(mtzfile)
    parseLine = base_mtzfile.split('.')
    root_mtzfile = parseLine[0]
    filename_out_mtz_f = root_mtzfile + '_f.mtz'

    fileexists = os.path.exists(filename_mtz)
    if fileexists != 0:
        os.remove(filename_mtz)

    if mtzfile.find('.mtz') > -1:
        datatype = 'ccp4'

        print 'Input file appears to be in CCP4 format'

    else:

        file = open(mtzfile,'r')
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

    # Case: Convert D*TREK to MTZ

    if datatype == 'dstartrek':

        print 'Input file appears to be from D*TREK'        

        for eachLine in allLines:

            if eachLine.find('SOURCE_WAVELENGTH') > -1:
                aList = eachLine.split()
                wavelength = aList[2]
                wavelength = wavelength.replace(';','')

            if eachLine.find('CRYSTAL_SPACEGROUP') > -1:
                aList = eachLine.split('=')
                spacegroup_no = aList[1]
                spacegroup_no = spacegroup_no.replace(';','')

        if spacegroup_no == 'none' or wavelength == 'none':
            print 'Spacegroup or wavelength were not found in the D*TREK reflection data file'
            time.sleep(4)
            return 1

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

    # Case: Convert SCALEPACK to MTZ 

    if datatype == 'scalepack':

        print 'Input file appears to be from SCALEPACK'

        # Obtain initial space group

        count = 0
        for eachLine in allLines:

            count = count + 1

            if count == 3:
                aLine = eachLine.split()
                line_length = len(aLine)
                if line_length == 7:
                    spacegroup_name = aLine[6]
                    spacegroup_name = spacegroup_name.upper()

        if spacegroup_name == 'none':
            print 'Spacegroup was not found in the SCALEPACK reflection data file'
            time.sleep(4)
            return 1

        #        

        filename_inp = 'mi_runscalepack2mtz.inp'
        filename_log = 'mi_scalepack2mtz.log'

        file = open(filename_inp,'w')
        file.write('name project noname crystal 1 dataset 1\n')
        file.write('symm ')
        file.write(spacegroup_name)
        file.write('\n')
        file.write('wave 0.0\n')
        file.write('end\n')
        file.close()

        runscalepack2mtz = 'scalepack2mtz hklin ' + filename_local_in + ' hklout ' + filename_mtz + ' < ' + filename_inp + ' > ' + filename_log

        os.system(runscalepack2mtz)

    # Intensity to F conversion

    if datatype == 'scalepack' or datatype == 'dstartrek':

        fileexists = os.path.exists(filename_mtz)
        if fileexists != 0:
            os.remove(filename_local_in)
            os.remove(filename_inp)
            os.remove(filename_log)
        else:
            print 'Output mtz file from SCALEPACK or D*TREK format conversion was not found'
            time.sleep(4)
            return 1    

        # Fix sort/asymmetric unit to standard

        print 'Running CCP4/CAD to set standard sort order'

        filename_inp = 'mi_cad.inp'
        filename_log = 'mi_cad.log'

        file = open(filename_inp,'w')
        file.write('LABIN FILE_NUMBER 1 ALL\n')
        file.write('SORT H K L \n')
        file.write('END\n')
        file.close()

        runcad = 'cad HKLIN1 ' + filename_mtz + ' HKLOUT ' + filename_out_mtz + ' < ' + filename_inp + ' > ' + filename_log

        os.system(runcad)

        fileexists = os.path.exists(filename_out_mtz)
        if fileexists != 0:
            os.remove(filename_inp)
            os.remove(filename_log)
            os.remove(filename_mtz)
        else:
            print 'Output mtz file from CAD sort was not found'
            time.sleep(4)
            return 1

        # Run TRUNCATE to get F's from I's

        print 'Running CCP4/TRUNCATE to reduce data from I to F'

        filename_inp = 'mi_truncate.inp'
        filename_log = 'mi_truncate.log'

        file = open(filename_inp,'w')
        file.write('title MR\n')
        file.write('\ntruncate yes\n')
        file.write('NRESIDUE 500\n')
        file.write('LABIN IMEAN=IMEAN SIGIMEAN=SIGIMEAN\n')
        file.write('LABOUT F=FP SIGF=SIGFP\n')
        file.write('NOHARVEST\n')
        file.write('END\n')
        file.close()

        runtruncate = 'truncate hklin ' + filename_out_mtz + ' hklout ' + filename_out_mtz_f + ' < ' + filename_inp + ' > ' + filename_log

        os.system(runtruncate)

        fileexists = os.path.exists(filename_out_mtz_f)
        if fileexists != 0:
            os.remove(filename_inp)
            os.remove(filename_out_mtz)
        else:
            print 'Output mtz file from CCP4/TRUNCATE was not found'
            time.sleep(4)
            return 1

        # Rewrite file with F, sd(F) only

        filename_inp = 'mi_cad.inp'
        filename_log = 'mi_cad.log'

        file = open(filename_inp,'w')
        file.write('labi file 1 E1=FP E2=SIGFP\n')
        file.write('\nEND\n')
        file.close()

        runcad = 'cad hklin1 ' + filename_out_mtz_f + ' hklout ' + filename_out_mtz_f_only + ' < ' + filename_inp + ' > ' + filename_log

        os.system(runcad)  

        fileexists = os.path.exists(filename_out_mtz_f_only)
        if fileexists != 0:
            os.remove(filename_inp)
            os.remove(filename_log)
            os.remove(filename_out_mtz_f)
        else:
            print 'Output mtz file from CCP4/CAD was not found'
            time.sleep(4)
            return 1

        # Add Rfree flags

        filename_inp = 'mi_rfree.inp'
        filename_log = 'mi_rfree.log'

        file = open(filename_inp,'w')
        file.write('FREERFRAC 0.05\n')
        file.write('END\n')
        file.close()

        runfreer = 'freerflag hklin ' + filename_out_mtz_f_only + ' hklout ' + filename_out_mtz_f + ' < ' + filename_inp + ' > ' + filename_log

        os.system(runfreer)

        fileexists = os.path.exists(filename_out_mtz_f)
        if fileexists != 0:
            os.remove(filename_out_mtz_f_only)
            os.remove(filename_inp)
            os.remove(filename_log)
        else:
            print 'Output mtz file from CCP4/FREERFLAG was not found'
            time.sleep(4)
            return 1

        # Final mtz file

        mtzfile = os.path.join(workingdir,filename_out_mtz_f)
            
    ###############################
    # Setup fixed coordinate file #
    ###############################

    mw_fixed = 0

    if pdbfile2 != 'none':

        file = open(pdbfile,'r')
        allLines = file.readlines()
        file.close()

        file = open('mi_molrep_fixed.pdb','w')

        for eachLine in allLines:

            eachLine = eachLine.strip()

            tag = eachLine[0:6]
            tag = tag.strip()

            if tag.find('ATOM') > -1 or tag.find('HETATM') > -1:

                if eachLine.find('HOH ') == -1 and eachLine.find('WAT ') == -1:
                    file.write(eachLine)
                    file.write('\n')

                    # Get approximate MW for PHASER

                    protein_element = eachLine[12:14]
                    protein_element = protein_element.strip()

                    if protein_element == 'C':
                        mw_fixed = mw_fixed + 12
                    if protein_element == 'N':
                        mw_fixed = mw_fixed + 14
                    if protein_element == 'O':
                        mw_fixed = mw_fixed + 16
                    if protein_element == 'SD' or protein_element == 'SG':
                        mw_fixed = mw_fixed + 32
            else:

                file.write(eachLine)
                file.write('\n')

        file.close()

    ###################
    # Setup MTZ file  #
    ###################

    file = open(mtzfile,'rb')
    allLines = file.readlines()
    file.close()

    file = open('mi_molrep.mtz','wb')
    file.writelines(allLines)
    file.close()

    # Extraction of cell data and MTZ labels for F and SD(F)

    file = open('mi_mtzdump.inp','w')
    file.write('HEADER\n')
    file.write('END\n')
    file.close()

    runmtz = 'mtzdump HKLIN mi_molrep.mtz < mi_mtzdump.inp > mi_mtzdump.log'
    os.system(runmtz)

    file = open('mi_mtzdump.log','r')
    allLines = file.readlines()
    file.close()

    os.remove('mi_mtzdump.log')
    os.remove('mi_mtzdump.inp')

    read_columns = 'no'
    read_labels = 'no'
    read_cell = 'no'

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

        if eachLine.find('* Space group =') > -1:
            parseLine = eachLine.split('number')
            mr_spacegroup_no_in = parseLine[1]
            mr_spacegroup_no_in = mr_spacegroup_no_in.replace(')','')
            mr_spacegroup_no_in = mr_spacegroup_no_in.strip()

            parseLine = eachLine.split(quote)
            space_group_mtz = parseLine[1]

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

        if eachLine.find('* Cell Dimensions :') > -1:
            read_cell = 'yes'

    # Check and setup cell/spacegroup for a CRYST1 record

    if acell_mtz == 'none' or bcell_mtz == 'none' or ccell_mtz == 'none' or alpha_mtz == 'none'\
       or beta_mtz == 'none' or gamma_mtz == 'none' or space_group_mtz == 'none':
        print 'Unable to evaluate cell or spacegroup from MTZ file'
        time.sleep(4)
        return 1
    else:

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

        cryst1 = 'CRYST1' + acell_mtz + bcell_mtz + ccell_mtz \
                   + alpha_mtz + beta_mtz + gamma_mtz + ' ' + space_group_mtz

    list_length = len(labelList)
    count = 0

    while count < list_length:

        if labelList[count] =='J' and ilabel == 'none':
            ilabel = colList[count]
        
        if labelList[count] == 'F' and flabel == 'none':
            flabel = colList[count]

        if labelList[count] == 'Q' and sigflabel == 'none':
            sigflabel = colList[count]

        count = count + 1

    # Handle case where input was intensity data by reducing with CCP4/TRUNCATE

    if ilabel != 'none' and flabel == 'none' and sigflabel != 'none':

        filename_inp = 'mi_truncate.inp'
        filename_log = 'mi_truncate.log'

        print 'Running CCP4/TRUNCATE to reduce data from I to F'

        file = open(filename_inp,'w')
        file.write('title MR\n')
        file.write('\ntruncate yes\n')
        file.write('NRESIDUE 500\n')
        file.write('LABIN IMEAN=')
        file.write(ilabel)
        file.write(' SIGIMEAN=')
        file.write(sigflabel)
        file.write('\n')
        file.write('LABOUT F=FP SIGF=SIGFP\n')
        file.write('NOHARVEST\n')
        file.write('END\n')
        file.close()

        runtruncate = 'truncate hklin mi_molrep.mtz hklout ' + filename_out_mtz_f + ' < ' + filename_inp + ' > ' + filename_log

        os.system(runtruncate)

        fileexists = os.path.exists(filename_out_mtz_f)
        if fileexists != 0:
            os.remove(filename_inp)
            os.remove('mi_molrep.mtz')
        else:
            print 'Output mtz file from CCP4/TRUNCATE was not found'
            time.sleep(4)
            return 1

        # Rewrite file with F, sd(F) only

        filename_inp = 'mi_cad.inp'
        filename_log = 'mi_cad.log'

        file = open(filename_inp,'w')
        file.write('labi file 1 E1=FP E2=SIGFP\n')
        file.write('END\n')
        file.close()

        runcad = 'cad hklin1 ' + filename_out_mtz_f + ' hklout ' + filename_out_mtz_f_only + ' < ' + filename_inp + ' > ' + filename_log

        os.system(runcad)  

        fileexists = os.path.exists(filename_out_mtz_f_only)
        if fileexists != 0:
            os.remove(filename_inp)
            os.remove(filename_log)
            os.remove(filename_out_mtz_f)
        else:
            print 'Output mtz file from CCP4/CAD was not found'
            time.sleep(4)
            return 1

        # Add Rfree flags

        filename_inp = 'mi_rfree.inp'
        filename_log = 'mi_rfree.log'

        file = open(filename_inp,'w')
        file.write('FREERFRAC 0.05\n')
        file.write('END\n')
        file.close()

        runfreer = 'freerflag hklin ' + filename_out_mtz_f_only + ' hklout ' + filename_out_mtz_f + ' < ' + filename_inp + ' > ' + filename_log

        os.system(runfreer)

        fileexists = os.path.exists(filename_out_mtz_f)
        if fileexists != 0:
            os.remove(filename_out_mtz_f_only)
            os.remove(filename_inp)
            os.remove(filename_log)
        else:
            print 'Output mtz file from CCP4/FREERFLAG was not found'
            time.sleep(4)
            return 1

        # Final mtz file

        file = open(filename_out_mtz_f,'rb')
        allLines = file.readlines()
        file.close()

        file = open('mi_molrep.mtz','wb')
        file.writelines(allLines)
        file.close()

        mtzfile = os.path.join(workingdir,filename_out_mtz_f)

        flabel = 'FP'
        sigflabel = 'SIGFP'

    # Standard case where input was amplitude data

    if flabel == 'none' or sigflabel == 'none':
        print 'MTZ labels for F and sd(F) could not be established'
        time.sleep(4)
        return 1

    ###########################################################
    # Check the point group for a space group search option   #
    ###########################################################

    if sg_search == 'yes':
        if mr_spacegroup_no_in != '1' and mr_spacegroup_no_in != '3' and mr_spacegroup_no_in != '16'\
           and mr_spacegroup_no_in != '75' and mr_spacegroup_no_in != '89' and mr_spacegroup_no_in != '143'\
           and mr_spacegroup_no_in != '149' and mr_spacegroup_no_in != '168' and mr_spacegroup_no_in != '177'\
           and mr_spacegroup_no_in != '195' and mr_spacegroup_no_in != '207':

            print 'Space group search option requires the input data to be in the symmetry point group'
            print '(one of P1, P2, P222, P4, P422, P3, P312, P6, P622, P23, P432)'
            time.sleep(4)
            return 1

    if sg_search == 'yes':

        if mr_spacegroup_no_in == '1':
            aList_spacegroups = sg_class_1
        if mr_spacegroup_no_in == '3':
            aList_spacegroups = sg_class_2    
        if mr_spacegroup_no_in == '16':
            aList_spacegroups = sg_class_222
        if mr_spacegroup_no_in == '75':
            aList_spacegroups = sg_class_4
        if mr_spacegroup_no_in == '89':
            aList_spacegroups = sg_class_422
        if mr_spacegroup_no_in == '143':
            aList_spacegroups = sg_class_3  
        if mr_spacegroup_no_in == '149':
            aList_spacegroups = sg_class_32
        if mr_spacegroup_no_in == '168':
            aList_spacegroups = sg_class_6
        if mr_spacegroup_no_in == '177':
            aList_spacegroups = sg_class_622
        if mr_spacegroup_no_in == '195':
            aList_spacegroups = sg_class_23  
        if mr_spacegroup_no_in == '207':
            aList_spacegroups = sg_class_432

    else:

        if spacegroup_user_input != 'none':
            aList_spacegroups.append(spacegroup_user_input)
        else:
            aList_spacegroups.append(mr_spacegroup_no_in)    

    number_space_groups = len(aList_spacegroups)

    #############################
    # Start of loop over models #
    #############################

    count_models = 0
    while count_models < number_models:

        pdbfile = aList_targetPDB[count_models]

        file = open(pdbfile,'r')
        allLines = file.readlines()
        file.close()

        file = open('mi_molrep.pdb','w')
        file.write(cryst1)
        file.write('\n')

        seq_insertion_flag = 'no'
        aList_true_seq_num = []

        mw = 0

        for eachLine in allLines:

            eachLine = eachLine.strip()

            tag = eachLine[0:6]
            tag = tag.strip()

            if tag == 'ATOM' or tag == 'HETATM':

                # Strip waters

                if eachLine.find('HOH ') == -1 and eachLine.find('WAT ') == -1:
                    file.write(eachLine)
                    file.write('\n')

                    # Get approximate MW for PHASER option

                    protein_element = eachLine[12:14]
                    protein_element = protein_element.strip()

                    if protein_element == 'C':
                        mw = mw + 12
                    if protein_element == 'N':
                        mw = mw + 14
                    if protein_element == 'O':
                        mw = mw + 16
                    if protein_element == 'SD' or protein_element == 'SG':
                        mw = mw + 32

                    # Check for sequence insertion flag

                    insert_code = eachLine[26:27]
                    if insert_code != ' ':
                        seq_insertion_flag = 'yes'

                    # Capture correct sequence id

                    true_seq_num = eachLine[22:27]
                    aList_true_seq_num.append(true_seq_num)

            else:

                file.write(eachLine)
                file.write('\n')

        file.close()

        #########################################
        # Start loop over possible spacegroups  #
        #########################################

        count_space_groups = 0
        while count_space_groups < number_space_groups:

            mr_spacegroup_no = aList_spacegroups[count_space_groups]

            # Set file names

            fileexists = os.path.exists(projectlog)
            if fileexists != 0:

                file = open(projectlog,'r')
                allLines = file.readlines()
                file.close()

                for eachLine in allLines:
                    if eachLine.find('Job ID') > -1 and eachLine.find('molrep') > -1:
                        aList = eachLine.split('_')
                        runid = aList[1]
                        runid_int = int(runid)

                runid_int = runid_int + 1
                runid = str(runid_int)

            job_id = job_prefix + runid

            filename_in = job_id + '.inp'
            filename_log = job_id + '.log'
            filename_pdb = job_id + '.pdb'
            filename_pdb_full = os.path.join(workingdir,filename_pdb)

            filename_pdb_phaser = job_id + '.1.pdb'
            filename_mtz_phaser = job_id + '.1.mtz'
            filename_sol_phaser = job_id + '.sol'
            filename_sum_phaser = job_id + '.sum'

            if mr_method == 'molrep':
                filename_log_full = os.path.join(workingdir,filename_log)
            else:
                filename_log_full = os.path.join(workingdir,filename_sum_phaser)

            ###################################################
            # Write application instructions - MOLREP option  #
            ###################################################

            if mr_method == 'molrep':

                # Clean up any debris

                fileexists = os.path.exists('molrep_mtz.cif')
                if fileexists != 0:
                    os.remove('molrep_mtz.cif')

                fileexists = os.path.exists('molrep.bat')
                if fileexists != 0:
                    os.remove('molrep.bat')

                fileexists = os.path.exists('molrep.xml')
                if fileexists != 0:
                    os.remove('molrep.xml')

                fileexists = os.path.exists('molrep.doc')
                if fileexists != 0:
                    os.remove('molrep.doc')

                # Write and run MOLREP job

                file = open(filename_in,'w')
                file.write('LABIN F=')
                file.write(flabel)
                file.write(' SIGF=')
                file.write(sigflabel)
                file.write('\n')
                file.write('NP 8\n')
                file.write('SURF O\n')
                file.write('RESMAX 3.0\n')
                file.write('SCORE N\n')
                file.write('STICK Y\n')
                file.write('PACK N\n')

                if number_molecules_input != 'auto':
                    file.write('NMON ')
                    file.write(number_molecules_input)
                    file.write('\n')

                if sg_search == 'yes' or spacegroup_user_input != 'none':
                    file.write('NOSG ')
                    file.write(mr_spacegroup_no)
                    file.write('\n')

                file.close()

                # Run process

                if pdbfile2 == 'none':
                    runmr = 'molrep HKLIN mi_molrep.mtz MODEL mi_molrep.pdb PATH_SCR ' + path_scratch + \
                            ' < ' + filename_in + ' > ' + filename_log
                else:
                    runmr = 'molrep HKLIN mi_molrep.mtz MODEL mi_molrep.pdb MODEL2 mi_molrep_fixed.pdb PATH_SCR ' + path_scratch + \
                            ' < ' + filename_in + ' > ' + filename_log

                print '\nStarting MOLREP process'
                print 'Job-ID:',job_id
                print 'Using mtz data:',flabel,',',sigflabel
                print 'Space group:',mr_spacegroup_no

                os.system(runmr)

            ###################################################
            # Write application instructions - PHASER option  #
            ###################################################

            if mr_method == 'phaser':

                print '\nStarting PHASER process'
                print 'Job-ID:',job_id
                print 'Using mtz data:',flabel,',',sigflabel
                print 'Space group:',mr_spacegroup_no

                mw = str(mw)
                mw_fixed = str(mw_fixed)

                # Run to estimate number of molecules

                if number_molecules_input == 'auto':

                    file = open(filename_in,'w')

                    file.write('MODE MR_CCA\n')
                    file.write('HKLIn mi_molrep.mtz\n')
                    file.write('LABIn F=')
                    file.write(flabel)
                    file.write(' SIGF=')
                    file.write(sigflabel)
                    file.write('\n')
                    file.write('COMPosition PROTein MW ')
                    file.write(mw)
                    file.write(' NUM 1\n')

                    if pdbfile2 != 'none':
                        file.write('COMPosition PROTein MW ')
                        file.write(mw_fixed)
                        file.write(' NUM 1\n')

                    file.write('ROOT ')
                    file.write(job_id)
                    file.write('\n')

                    if sg_search == 'yes' or spacegroup_user_input != 'none':
                        file.write('SPACegroup ')
                        file.write(mr_spacegroup_no)
                        file.write('\n')

                    file.close()

                    runmr = 'phaser ' + ' < ' + filename_in + ' > ' + filename_log

                    os.system(runmr)

                    fileexists = os.path.exists(filename_log)
                    if fileexists != 0:

                        file = open(filename_log,'r')
                        allLines = file.readlines()
                        file.close()

                        os.remove(filename_log)

                        for eachLine in allLines:
                            if eachLine.find('<== most probable') > -1:
                                parseLine = eachLine.split()
                                number_molecules = parseLine[0]

                    fileexists = os.path.exists(filename_sol_phaser)
                    if fileexists != 0:
                        os.remove(filename_sol_phaser)

                    os.remove(filename_in)

                else:

                    # User defined number of molecules

                    number_molecules = number_molecules_input

                # Check number of molecules is reasonable and allow a few clashes

                int_number_molecules = int(number_molecules)

                if int_number_molecules > 6:
                    int_number_molecules = 6

                    print '\nWARNING: number of molecules to find in PHASER run limited to 6\n'

                pack = 5*int_number_molecules
                pack = str(pack)

                # Run MR

                file = open(filename_in,'w')

                file.write('MODE MR_AUTO\n')
                file.write('HKLIn mi_molrep.mtz\n')
                file.write('LABIn F=')
                file.write(flabel)
                file.write(' SIGF=')
                file.write(sigflabel)
                file.write('\n')
                file.write('RESOLUTION 3.0\n')
                file.write('ENSEmble mi_molrep PDB mi_molrep.pdb IDENtity 100\n')

                if pdbfile2 != 'none':
                    file.write('ENSEmble mi_molrep_fixed PDB mi_molrep_fixed.pdb IDENtity 100\n')

                file.write('COMPosition PROTein MW ')
                file.write(mw)
                file.write(' NUM ')
                file.write(number_molecules)
                file.write('\n')

                if pdbfile2 != 'none':
                    file.write('COMPosition PROTein MW ')
                    file.write(mw_fixed)
                    file.write(' NUM 1\n')

                file.write('SEARch ENSEmble mi_molrep NUM ')
                file.write(number_molecules)
                file.write('\n')       
                file.write('ROOT ')
                file.write(job_id)
                file.write('\n')
                file.write('VERBose ON')
                file.write('\nPACK ')
                file.write(pack)
                file.write('\n')

                if sg_search == 'yes' or spacegroup_user_input != 'none':
                    file.write('SPACegroup ')
                    file.write(mr_spacegroup_no)
                    file.write('\n')

                file.close()

                runmr = 'phaser ' + ' < ' + filename_in + ' > ' + filename_log

                os.system(runmr)

                # Rename Phaser pdb file name to standard name for postprocessing 

                fileexists = os.path.exists(filename_pdb_phaser)
                if fileexists != 0:
                    os.rename(filename_pdb_phaser,'molrep.pdb')

            # Remove all debris

            fileexists = os.path.exists('molrep_mtz.cif')
            if fileexists != 0:
                os.remove('molrep_mtz.cif')

            fileexists = os.path.exists('molrep.bat')
            if fileexists != 0:
                os.remove('molrep.bat')

            fileexists = os.path.exists('molrep.xml')
            if fileexists != 0:
                os.remove('molrep.xml')

            fileexists = os.path.exists('molrep.doc')
            if fileexists != 0:
                os.remove('molrep.doc')

            fileexists = os.path.exists(filename_mtz_phaser)
            if fileexists != 0:
                os.remove(filename_mtz_phaser)

            fileexists = os.path.exists(filename_sol_phaser)
            if fileexists != 0:
                os.remove(filename_sol_phaser)

            os.remove(filename_in)

            # Check solution coordinate file exists

            fileexists = os.path.exists('molrep.pdb')
            if fileexists == 0:
                file = open('molrep.pdb','w')
                file.write('REMARK  No coordinates were created from this run\n')
                file.write('REMARK  Packing function may have excluded all solutions\n')
                file.write('REMARK  Check MR log file\n')
                file.close()

            ########################
            # PDB File corrections #
            ########################

            # Fix-up PDB format per atom per justifications

            file = open('molrep.pdb','r')
            allLines = file.readlines()
            file.close()

            os.remove('molrep.pdb')

            mr_number_atoms_out = 0

            file = open('molrep_fixed.pdb','w')

            for eachLine in allLines:

                eachLine = eachLine.strip()
                tag = eachLine[0:6]
                tag = tag.strip()

                if tag == 'ATOM' or tag == 'HETATM':

                    mr_number_atoms_out = mr_number_atoms_out + 1

                    # MOLREP and PHASER do not maintain proper justification
                    # Fix for common ions

                    atom_name = eachLine[12:16]

                    if atom_name == ' NA ':
                        atom_name = 'NA  '

                    if atom_name == ' MG ':
                        atom_name = 'MG  '

                    if atom_name == ' CL ':
                        atom_name = 'CL  '

                    if atom_name == ' CR ':
                        atom_name = 'CR  '

                    if atom_name == ' MN ':
                        atom_name = 'MN  '

                    if atom_name == ' FE ':
                        atom_name = 'FE  '

                    if atom_name == ' CO ':
                        atom_name = 'CO  '

                    if atom_name == ' NI ':
                        atom_name = 'NI  '

                    if atom_name == ' CU ':
                        atom_name = 'CU  '

                    if atom_name == ' ZN ':
                        atom_name = 'ZN  '

                    if atom_name == ' SE ':
                        atom_name = 'SE  '

                    if atom_name == ' BR ':
                        atom_name = 'BR  '

                    if atom_name == ' CS ':
                        atom_name = 'CS  '

                    atom_out = eachLine[0:12] + atom_name + eachLine[16:80]

                    file.write(atom_out)
                    file.write('\n')

                else:

                    file.write(eachLine)
                    file.write('\n')

            file.close()

            os.rename('molrep_fixed.pdb','molrep.pdb')

            # If engine was CCP4/MOLREP need to fix sequence numbers for case of sequence inserts

            if mr_method == 'molrep' and seq_insertion_flag == 'yes' and pdbfile2 == 'none' and mr_number_atoms_out > 0:

                print 'Correcting MOLREP output PDB file to match input sequence numbers'

                file = open('molrep.pdb','r')
                allLines = file.readlines()
                file.close()

                os.remove('molrep.pdb')

                mr_number_atoms_in = len(aList_true_seq_num)
                mr_atom_count = 0

                print 'Input and output number of atoms:',mr_number_atoms_in,mr_number_atoms_out

                file = open('molrep_fixed.pdb','w')

                for eachLine in allLines:

                    eachLine = eachLine.strip()
                    tag = eachLine[0:6]
                    tag = tag.strip()

                    if tag == 'ATOM' or tag == 'HETATM':

                        true_seq_num = aList_true_seq_num[mr_atom_count]
                        mr_atom_count = mr_atom_count + 1

                        if mr_atom_count == mr_number_atoms_in:
                            mr_atom_count = 0
                        
                        mr_atom_out = eachLine[0:22] + true_seq_num + eachLine[27:80]
                        
                        file.write(mr_atom_out)
                        file.write('\n')

                    else:

                        file.write(eachLine)
                        file.write('\n')

                file.close()

                os.rename('molrep_fixed.pdb','molrep.pdb')
            
            ################################################################
            # Option to match output to input crystal (lattice) position.  #
            # Sensitive to matching exact ATOM records                     #
            ################################################################

            if match_pdbin == 'yes':

                print 'Applying symmetry to match coordinates to the input position'

                runreforigin = 'reforigin XYZIN molrep.pdb XYZREF mi_molrep.pdb XYZOUT molrep_origin.pdb DMAX 7.0 > reforigin.log 2> reforigin_err.log'
                os.system(runreforigin)

                fileexists = os.path.exists('molrep_origin.pdb')
                if fileexists != 0:
                    os.remove('molrep.pdb')
                    os.rename('molrep_origin.pdb','molrep.pdb')
                else:
                    print 'Warning - REFORIGIN run seems to have failed'
                    time.sleep(4)

                fileexists = os.path.exists('reforigin.log')
                if fileexists != 0:
                    os.remove('reforigin.log')     

                fileexists = os.path.exists('reforigin_err.log')
                if fileexists != 0:
                    os.remove('reforigin_err.log')                   

            ##############################################            
            # Parse diagnostics from MR log and clean-up #
            ##############################################

            fileexists = os.path.exists(filename_pdb)
            if fileexists != 0:
                os.remove(filename_pdb)

            os.rename('molrep.pdb',filename_pdb)
            print 'Output PDB file:',filename_pdb

            fileexists = os.path.exists(filename_log)
            if fileexists != 0:
                file = open(filename_log,'r')
                allLines = file.readlines()
                file.close()

                # Parser for MOLREP

                if mr_method == 'molrep':

                    print 'Output MOLREP log:',filename_log

                    # Parse R per format in CCP4 6.0.1

                    for eachLine in allLines:
                        if eachLine.find('S__ ') > -1:
                            aString = eachLine.strip()
                            parseLine = aString.split()
                            mr_r = parseLine[11]
                            print 'R=',mr_r

                    # Parse R per format in CCP4 5.0.2 

                    if mr_r == 'none':

                        for eachLine in allLines:
                            if eachLine.find('Sol_Mon') > -1:
                                aString = eachLine.strip()
                                parseLine = aString.split()
                                mr_r = parseLine[7]
                                print 'R=',mr_r 

                    if mr_r == 'none':
                        mr_r = '1.000'

                if mr_method == 'phaser':

                    print 'Output PHASER log:',filename_log

                    for eachLine in allLines:
                        if eachLine.find('The R-factor') > -1:
                            aString = eachLine.strip()
                            parseLine = aString.split()
                            mr_r = parseLine[3]
                            mr_r = float(mr_r)
                            mr_r = 0.01 * mr_r
                            mr_r = '%.3f'%(mr_r)
                            mr_r = str(mr_r)

                    if mr_r == 'none':
                        mr_r = '1.000'

            else:

                print 'The MR log file was not found'
                time.sleep(4)
                return 1

            #########################
            # Append project log    #
            #########################

            if projectlog != 'none':

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
                file.write('\nInput fixed atoms: ')
                file.write(pdbfile2)
                file.write('\nOutput atoms: ')
                file.write(filename_pdb_full)
                file.write('\nOutput log: ')
                file.write(filename_log_full)
                file.write('\nMR method: ')
                file.write(mr_method)
                file.write('\nMR space group: ')
                file.write(mr_spacegroup_no)
                file.write('\n')
                file.write('Summary: R=')
                file.write(mr_r)
                file.write('\n')    
                file.write('---------------\n')
                file.close()

            count_space_groups = count_space_groups + 1   

        ##################################
        # End of loop over space groups  #
        ##################################

        count_models = count_models + 1

        fileexists = os.path.exists('mi_molrep.pdb')
        if fileexists != 0:
            os.remove('mi_molrep.pdb')

    ###########################
    # End of loop over models #
    ###########################

    # Clean-up local junk and scratch files

    fileexists = os.path.exists('mi_molrep.mtz')
    if fileexists != 0:
        os.remove('mi_molrep.mtz')

    fileexists = os.path.exists('mi_molrep_fixed.pdb')
    if fileexists != 0:
        os.remove('mi_molrep_fixed.pdb')

    fileexists = os.path.exists('mi_molrep_fixed.pdb')
    if fileexists != 0:
        os.remove('mi_molrep_fixed.pdb')

    fileexists = os.path.exists('molrep.btc')
    if fileexists != 0:
        os.remove('molrep.btc')   

    dir_list = os.listdir(workingdir)
    number_files = len(dir_list)

    count = 0
    while count < number_files:
        target_file = dir_list[count]

        if target_file.find(path_scratch) > -1:
            os.remove(target_file)

        count = count + 1
    #

    time.sleep(4)

    #
    return 0

if __name__ == "__main__":
    sys.exit(Run())
