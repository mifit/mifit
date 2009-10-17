#####################################################################
#                                                                   #
# Script to automatically builds complete restraint libraries for   #
# for CCP4/REFMAC5:                                                 #
# (1) Combines novel ligand dictionaries                            #
# (2) Converts 'minimal' to full dictionaries                       #
# (3) Adds special covalent links given by pdb LINK records         #
#                                                                   #
# Copyright: Molecular Images   2005                                #
#                                                                   #
# This script is distributed under the same conditions as MIFit     #
#                                                                   #
#####################################################################

import sys
import os
import string
import time
import getopt
import ccp4check

def Usage():
    print "Usage: %s [options]" % sys.argv[0]
    print "Options are:"
    print "  -p,--pdbfile=FILE        the pdbfile"
    print "  -d,--workdir=DIR         the working dir"
    print "  -?,--help                this help"

def Run(argv=None):
    if argv is None:
        argv=sys.argv

    # Initialize

    cryst_flag = 'no'
    pdbfile = 'none'
    projectlog = 'none'
    workingdir = 'none'
    parseLine = []
    quote = """'"""

    runid = '1'
    runid_int = 0
    projectlog = 'project_history.txt'
    job_prefix = 'restraints_'

    aList_hets = []
    aList_hets_names = []
    aList_hets_nonPDB = []
    aList_chain_id = []
    aList_res_name = []
    aList_atom_name = []
    aList_res_number = []
    aList_entity_file = []
    aList_x = []
    aList_y = []
    aList_z = []

    aList_connect = []
    aList_restype1 = []
    aList_restype2 = []
    aList_resnumber1 = []
    aList_resnumber2 = []
    aList_atom1 = []
    aList_atom2 = []
    aList_link_id = []

    aList_protein = ['ALA','ARG','ASN','ASP','CYS','GLN','GLU','GLY','HIS','ILE','LEU',\
                     'LYS','MET','PHE','PRO','SER','THR','TRP','TYR','VAL','MSE','HOH']

    number_protein_list = len(aList_protein)

    # Read args
    number_of_args = len(argv)
    args = argv[1:]
    optlist, args = getopt.getopt(args,'p:d:?:',
                                  ['pdbfile=','workdir=','help'])
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
            elif arg_value == '-d' or arg_value=="--workdir":
                workingdir = param_value
        count=count+1

    ccp4,error = ccp4check.ccp4check()
    if not ccp4:
      print '\n' + error + '\n'
      time.sleep(4)
      return 1
    
    # Check inputs
    fileexists = os.path.exists(pdbfile)
    if fileexists == 0:
        print 'The PDB file was not found ',pdbfile
        time.sleep(4)
        return 1

    fileexists = os.path.exists(workingdir)
    if fileexists == 0:
        print 'The working directory was not found ',workingdir
        time.sleep(4)
        return 1
    else:
        os.chdir(workingdir)
        pwd = os.getcwd()

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

    # Temporary files 

    temp_pdb = os.path.join(working_ccp4_scratch,'mi_temp.pdb')
    temp_lib = os.path.join(working_ccp4_scratch,'mi_temp.lib')
    output_lib = os.path.join(working_ccp4_scratch,'mi_temp_all.lib')
    output_pdb = os.path.join(working_ccp4_scratch,'mi_temp_out.pdb')

    # Establish name for final complete output library

    libbase_pdb = os.path.basename(pdbfile)
    libbase_cif = libbase_pdb.replace('.pdb','.lib')
    libbase_full = 'restraints_' + libbase_cif
    final_lib = os.path.join(pwd,libbase_full)

    # Check project history

    fileexists = os.path.exists(projectlog)
    if fileexists != 0:

        file = open(projectlog,'r')
        allLines = file.readlines()
        file.close()

        for eachLine in allLines:
            if eachLine.find('Job ID') > -1 and eachLine.find(job_prefix) > -1:
                parseLine = eachLine.split('_')
                runid = parseLine[1]
                runid_int = int(runid)

        runid_int = runid_int + 1
        runid = str(runid_int)

    job_id = job_prefix + runid

    print '\nChecking and preparing a complete REFMAC5 library'
    print 'Job-ID:',job_id

    ######################################################
    # Inspect coordinate data to obtain entity lists etc #
    ######################################################

    file = open(pdbfile,'r')
    allLines = file.readlines()
    file.close()

    file = open(temp_pdb,'w')

    for eachLine in allLines:

        # Write input PDB to local file

        file.write(eachLine)

        # Check record

        tag = eachLine[0:6]
        tag = tag.strip()

        # Check for CRYST1 record

        if tag == 'CRYST1':
            cryst_flag = 'yes'

        # Capture non-gap link data

        if tag == 'LINK':

            link_record = str(eachLine)

            if link_record.find('gap') == -1:
                aList_connect.append(link_record)

        # Capture atom data for potential use later

        if tag == 'ATOM' or tag == 'HETATM':

            chain_id = eachLine[21:22]
            res_name = eachLine[17:20]
            atom_name = eachLine[13:16]
            res_number = eachLine[22:26]
            x = eachLine[30:38]
            y = eachLine[38:46]
            z = eachLine[46:54]

            aList_chain_id.append(chain_id)
            aList_res_name.append(res_name)
            aList_atom_name.append(atom_name)
            aList_res_number.append(res_number)
            aList_x.append(x)
            aList_y.append(y)
            aList_z.append(z)

            res_name = res_name.strip()

            # Check if the current entity is standard

            res_name = res_name.strip()

            count = 0
            found = 'no'
            while count < number_protein_list:
                if res_name == aList_protein[count]:
                    entity_type = 'P'
                    found = 'yes'

                count = count + 1

            if found == 'no':
                entity_type = 'L'

            # Check if we already have this entity        

            if entity_type == 'L':

                repeat = 'no'
                count = 0
                count_hets = len(aList_hets)
                while count < count_hets:
                    if res_name == aList_hets[count]:
                        repeat = 'yes'

                    count = count + 1

                # If a new entity, add to the list

                if repeat == 'no':
                    aList_hets.append(res_name)

    file.close()

    number_atoms = len(aList_atom_name)
    number_links = len(aList_connect)

    if cryst_flag == 'no':
        print '\nThe coordinate file must contain a CRYST1 record\n'
        time.sleep(4)
        return 1     

    #######################################
    # Obtain the known entity name lists  #
    #######################################

    fileexists = os.path.exists(ccp4.entitylist)
    if fileexists != 0:
        file = open(ccp4.entitylist,'r')
        allLines = file.readlines()
        file.close()

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

        if found == 'no':
            aList_hets_names.append('novel-entity')
            aList_hets_nonPDB.append(pdbentity)

        count = count + 1

    # List PDB name assignments from the REFMAC5 list

    if count_hets > 0:
        print '\nStandard PDB (REFMAC) entity assignments'
        print '========================================'

        count = 0
        while count < count_hets:

            prentitycode = aList_hets[count]
            prentityname = aList_hets_names[count]
            print prentitycode,prentityname

            count = count + 1

    # Record links

    if number_links > 0:
        print '\nNon-gap LINK records found:',number_links
        print 'These distances must be edited from ? in the dictionary file data_link_ loop'

    ##############################################################
    # Check we have dictionary files to deal with novel entities #
    ##############################################################

    number_nonPDB = len(aList_hets_nonPDB)

    if number_nonPDB > 0:

        count = 0
        while count < number_nonPDB:

            entity = aList_hets_nonPDB[count]
            entity_file = entity + '.cif'

            # Check file exists

            fileexists = os.path.exists(entity_file)

            if fileexists == 0:

                print '\nDid not find local dictionary file for entity: ',entity_file,'\n'
                time.sleep(4)
                return 1

            else:

                # Check that this really file contains records for this ligand

                ligand_records = 'no'

                file = open(entity_file, 'r')
                allLines = file.readlines()
                file.close()

                aList_entity_file.append(entity_file)

                for eachLine in allLines:
                    if eachLine.find(entity) > -1:
                        ligand_records = 'yes'

                if ligand_records == 'no':
                    print '\nFile did not contain ',entity,' records\n'
                    time.sleep(4)
                    return 1

            count = count + 1

    ##############################################################################
    # Combine novel ligand dictionaries with LIBCHECK if there is more than one  #
    ##############################################################################

    fileexists = os.path.exists(temp_lib)
    if fileexists != 0:
        os.remove(temp_lib)

    # Copy first input library to initialize copy in CCP4 temporary space

    if number_nonPDB > 0:

        entitylib = aList_entity_file[0]

        file = open(entitylib, 'r')
        allLines = file.readlines()
        file.close()

        file = open(temp_lib, 'w')
        for eachLine in allLines:
            file.write(eachLine)
        file.close()

    # Build-up multiple ligand dictionaries

    if number_nonPDB > 1:

        count = 1
        while count < number_nonPDB:

            # Get next ligand

            entity = aList_hets_nonPDB[count]
            entitylib_add = aList_entity_file[count]

            filename = 'runlibcheck.inp'
            file = open(filename, 'w')
            file.write('N\n')
            file.write('_FILE_L ')
            file.write(entitylib_add)
            file.write('\n_FILE_L2 ')
            file.write(temp_lib)
            file.write('\n')
            file.write('\n')
            file.close()

            runlibcheck = 'libcheck < runlibcheck.inp > libcheck.log'

            os.system(runlibcheck)

            # Reset complete (base) dictionary file

            fileexists = os.path.exists('libcheck.lib')
            if fileexists == 0:
                print 'The output library file from LIBCHECK was not created'
                print 'Entity ',entity,' may already be in the standard dictionaries'
                time.sleep(4)
                return 1
            else:
                os.remove(temp_lib)
                os.rename('libcheck.lib',temp_lib)

            # Clean-up

            fileexists = os.path.exists('libcheck.out')
            if fileexists != 0:
                os.remove('libcheck.out')

            entity_libcheck_pdb = 'libcheck_' + entity + '.pdb'
            entity_libcheck_ps = 'libcheck_' + entity + '.ps'
            entity_libcheck_cif = 'libcheck_' + entity + '.cif'

            fileexists = os.path.exists(entity_libcheck_pdb)
            if fileexists != 0:        
                os.remove(entity_libcheck_pdb)

            fileexists = os.path.exists(entity_libcheck_ps)
            if fileexists != 0:            
                os.remove(entity_libcheck_ps)

            fileexists = os.path.exists(entity_libcheck_cif)
            if fileexists != 0:        
                os.remove(entity_libcheck_cif)           

            count = count + 1

    #############################################
    # Run CCP4/REFMAC run to generate minimal   #
    # description ligands and handle link data  #
    #############################################

    fileexists = os.path.exists(output_lib)
    if fileexists != 0:
        os.remove(output_lib)

    filename = 'runrefmac5.inp'
    file = open(filename, 'w')
    file.write('MODE NEWEntry\n')
    file.write('MAKE_RESTRAINTS CHECK None\n')
    file.write('MAKE_RESTRAINTS EXIT Yes\n')
    file.write('END\n')
    file.close()

    fileexists = os.path.exists(temp_lib)
    if fileexists == 0:        
        runrefmac = 'refmac5 XYZIN ' + temp_pdb + ' XYZOUT ' + output_pdb + ' LIB_OUT ' + output_lib + ' < runrefmac5.inp > mi_refmac.log'
    else:
        runrefmac = 'refmac5 XYZIN ' + temp_pdb + ' XYZOUT ' + output_pdb + ' LIB_IN ' + temp_lib + ' LIB_OUT ' + output_lib + ' < runrefmac5.inp > mi_refmac.log'

    os.system(runrefmac)

    fileexists = os.path.exists(output_pdb)
    if fileexists == 0:
        print 'REFMAC5 dictionary calculation for PDB entities failed - check "mi_refmac.log"'
        print 'A common problem is a standard PDB ligand code but atom names inconsistent with a standard dictionary\n'
        time.sleep(4)
        return 1
    else:
        os.remove(temp_pdb)
        os.remove(output_pdb)
        os.remove('mi_refmac.log')
        os.remove('runrefmac5.inp')

    # If there was no output library (i.e. input library was already sufficient) rename input to output

    fileexists = os.path.exists(output_lib)

    if fileexists == 0:

        fileexists = os.path.exists(temp_lib)
        if fileexists != 0:
            os.rename(temp_lib,output_lib)

    fileexists = os.path.exists(temp_lib)
    if fileexists != 0:
        os.remove(temp_lib)   

    # Splice any additional LINKS into the final dictionary

    if number_links > 0:

        fileexists = os.path.exists(output_lib)
        if fileexists != 0:

            file = open(output_lib, 'r')
            allLines = file.readlines()
            file.close()

            # Check that links are not already inserted

            set_links = 'no'

            for eachLine in allLines: 
                if eachLine.find('data_link_list') > -1:
                    print 'Dictionary file already contains LINK records so none can be added'
                    set_links = 'yes'

            if set_links == 'no':

                print 'Adding covalent LINK records to dictionary file'
                if fileexists != 0:
                    os.remove(output_lib)

                file = open(output_lib, 'w')

                for eachLine in allLines:

                    # Insert list of links

                    if eachLine.find('DESCRIPTION OF MONOMERS') > -1:

                        file.write('#\n')
                        file.write('# ---   LIST OF LINKS ---\n')
                        file.write('#\n')
                        file.write('data_link_list\n')
                        file.write('loop_\n')
                        file.write('_chem_link.id\n')
                        file.write('_chem_link.comp_id_1\n')
                        file.write('_chem_link.mod_id_1\n')
                        file.write('_chem_link.group_comp_1\n')
                        file.write('_chem_link.comp_id_2\n')
                        file.write('_chem_link.mod_id_2\n')
                        file.write('_chem_link.group_comp_2\n')
                        file.write('_chem_link.name\n')

                        count = 0
                        while count < number_links:

                            aLine = aList_connect[count]

                            # parse link information

                            atom1 = aLine[12:16]
                            restype1 = aLine[17:20]
                            chain1 = aLine[21:22]
                            resnumber1 = aLine[22:26]
                            atom2 = aLine[42:46]
                            restype2 = aLine[47:50]
                            chain2 = aLine[51:52]
                            resnumber2 = aLine[52:56]
                            link_id = aLine[72:80]

                            atom1 = atom1.strip()
                            restype1 = restype1.strip()
                            chain1 = chain1.strip()
                            resnumber1 = resnumber1.strip()
                            atom2 = atom2.strip()
                            restype2 = restype2.strip()
                            chain2 = chain2.strip()
                            resnumber2 = resnumber2.strip()
                            link_id = link_id.strip()

                            if atom1.find(quote) > -1:
                                atom1_out = quote + atom1 + quote
                            else:
                                atom1_out = atom1

                            if atom2.find(quote) > -1:
                                atom2_out = quote + atom2 + quote
                            else:
                                atom2_out = atom2

                            aList_restype1.append(restype1)
                            aList_restype2.append(restype2)
                            aList_resnumber1.append(resnumber1)
                            aList_resnumber2.append(resnumber2)
                            aList_atom1.append(atom1_out)
                            aList_atom2.append(atom2_out)
                            aList_link_id.append(link_id)

                            label_out = link_id + '  ' + restype1 + ' . . ' + restype2 + ' . . bond_' + restype1 + '_=_' + restype2

                            file.write(label_out)
                            file.write('\n')

                            count = count + 1

                        file.write('#\n')

                    # Copy out input CIF line

                    file.write(eachLine)

            # Write LINK specifications

            file.write('#\n')
            file.write('# --- DESCRIPTION OF LINKS ---\n')

            count = 0
            while count < number_links:

                restype1 = aList_restype1[count]
                restype2 = aList_restype2[count]
                atom1 = aList_atom1[count]
                atom2 = aList_atom2[count]
                link_id = aList_link_id[count]
                dist = '?'

                label1 = count*2 + 1
                label2 = label1 + 1
                label1 = str(label1)
                label2 = str(label2)

                file.write('#\n')
                link_data = 'data_link_' + link_id
                file.write(link_data)
                file.write('\n#\n')
                file.write('loop_\n')
                file.write('_chem_link_bond.link_id\n')
                file.write('_chem_link_bond.atom_1_comp_id\n')
                file.write('_chem_link_bond.atom_id_1\n')
                file.write('_chem_link_bond.atom_2_comp_id\n')
                file.write('_chem_link_bond.atom_id_2\n')
                file.write('_chem_link_bond.type\n')
                file.write('_chem_link_bond.value_dist\n')
                file.write('_chem_link_bond.value_dist_esd\n')

                label_out = ' ' + link_id + ' ' + label1 + ' ' + atom1 + ' ' + label2 + ' ' + atom2 + ' . ' + dist + ' 0.020'
                file.write(label_out)
                file.write('\n')

                count = count + 1

            file.write('#-------------------------------------------------------\n')
            file.close()

    # Copy the ultimate library back to the working directory as refmac5.lib

    fileexists = os.path.exists(output_lib)
    if fileexists != 0:

        print '\nBuilding a complete REFMAC5 library:',final_lib,'\n'

        fileexists = os.path.exists(final_lib)
        if fileexists !=0:
            os.remove(final_lib)

        os.rename(output_lib,final_lib)
        library_needed = final_lib

    else:

        print '\nNo output library needed - REFMAC5 already has libraries for all entities\n'
        library_needed = 'none needed'

    ########################
    # Append project log   #
    ########################

    print 'Writing project log'
    time.sleep(4)

    runtime = time.ctime(time.time())

    file = open(projectlog,'a')
    file.seek(0,2)
    file.write('Job ID: ')
    file.write(job_id)
    file.write('\nDate: ')
    file.write(runtime)
    file.write('\nInput atoms: ')
    file.write(pdbfile)
    file.write('\nOutput REFMAC5 library: ')
    file.write(library_needed)
    file.write('\n')

    count = 0
    while count < count_hets:
        prentitycode = aList_hets[count]
        prentityname = aList_hets_names[count]
        entity_list_number = count + 1
        prnumber = str(entity_list_number)
        file.write('Entity_')
        file.write(prnumber)
        file.write(': ')
        file.write(prentitycode)
        file.write(' ')
        file.write(prentityname)
        file.write('\n')

        count = count + 1

    file.write('---------------\n')
    file.close()

    #  Clean-up job-specific temporary CCP4_SCR space and put back original CCP4_SCR

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

    #
    return 0

if __name__ == "__main__":
    sys.exit(Run())
