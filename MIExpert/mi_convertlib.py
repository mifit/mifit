#####################################################################
#                                                                   #
# Converts a mmCIF restraint dictionary to SHELX or CNS/CNX format  #
#                                                                   #
# Copyright: Molecular Images    2005                               #
#                                                                   #
# This script is distributed under the same conditions as MIFit     #
#                                                                   #
#####################################################################

import sys
import os
import math
import time
import getopt

def Usage():
    print "Usage: %s [Parameters] " % sys.argv[0]
    print "Parameters are:"
    print "  -c, --cif=FILE         cif file to use"
    print "  -w, --workdir=DIR      working dir to use"
    print "  -r, --refprogram=PROG  shelx (default), cns, or cnx"
    print ""
    print "  -?, --help              for this help message"
    print ""

def Run(argv=None):
    if argv is None:
        argv=sys.argv
    # Initialize

    ciffile = 'none'
    workingdir = 'none'
    refprogram = 'shelx'
    write_mode = 'new'

    runid = '1'
    runid_int = 0
    projectlog = 'project_history.txt'
    job_prefix = 'convertlib_'

    twopi = 6.28
    parseLine = []

    # Defaults for CNX

    bond_weight = '1000.0'
    angle_weight = '500.0'
    improper_weight = '750.0'
    chiral_sign = 'none'
    improper_chiral_angle = '35.264'

    # For atom lists

    aList_atom_id = []
    aList_atom_symbol = []
    aList_atom_type = []
    aList_atom_mass = []
    aList_atom_vdv = []
    aList_atom_x = []
    aList_atom_y = []
    aList_atom_z = []
    chem_comp_atom = []
    num_chem_comp_atom = 0
    chem_comp_atom_head = 'no'
    extract_atom_header = 'no'
    read_atom = 'no'
    atom_x = 'none'
    atom_y = 'none'
    atom_z = 'none'

    # For bonds

    aList_atom1_bond = []
    aList_atom2_bond = []
    aList_atom12_dist_bond = []
    chem_comp_bond = []
    num_chem_comp_bond = 0
    chem_comp_bond_head = 'no'
    extract_bond_header = 'no'
    read_bond = 'no'

    # For angles

    aList_atom1_angle = []
    aList_atom2_angle = []
    aList_atom3_angle = []
    aList_atom13_dist_angle = []
    aList_angle13 = []
    chem_comp_angle = []
    num_chem_comp_angle = 0
    chem_comp_angle_head = 'no'
    extract_angle_header = 'no'
    read_angle = 'no'

    # For chirals

    aList_chiral_id = []
    aList_chiral_atom = []
    aList_chiral_atom1 = []
    aList_chiral_atom2 = []
    aList_chiral_atom3 = []
    aList_chiral_sign = []
    chem_comp_chiral = []
    num_chem_comp_chiral = 0
    chem_comp_chiral_head = 'no'
    extract_chiral_header = 'no'
    read_chiral = 'no'

    # For planes

    aList_plane_id = []
    aList_plane_atom = []
    chem_comp_plane = []
    aList_quartet = []
    aList_improper_planes = []
    num_chem_comp_plane = 0
    chem_comp_plane_head = 'no'
    extract_plane_header = 'no'
    read_plane = 'no'

    # Read args
    args = argv[1:]
    optlist, args = getopt.getopt(
        args,'c:w:r:?',['cif=','workdir=','refprogram=','help'])
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
            if arg_value == '-c' or arg_value=='--cif':
                ciffile = param_value
            elif arg_value == '-w' or arg_value=='--workdir':
                workingdir = param_value
            elif arg_value == '-r' or arg_value=='--refprogram':
                refprogram = param_value
        count=count+1

    # Check inputs exist
    fileexists = os.path.exists(ciffile)
    if fileexists == 0:
        print 'The CIF file was not found ',ciffile
        time.sleep(4)
        return 1

    fileexists = os.path.exists(workingdir)
    if fileexists == 0:
        print 'The working directory was not found ',workingdir
        time.sleep(4)
        return 1
    else:
        os.chdir(workingdir)

    if refprogram != 'shelx' and refprogram != 'cns' and refprogram != 'cnx':
        print 'Refinement program not recognized - options are shelx/cnx/cns'
        time.sleep(4)
        return 1

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

    #############################################
    # Parse CIF library file for restraint data #
    #############################################

    file = open(ciffile,'r')
    allLines = file.readlines()
    file.close()

    for eachLine in allLines:

        # Establish entity name

        if eachLine.find('data_comp_') > -1:
            entityname = eachLine[10:13]

        # Identify CIF tags from types found in dictionary files

        ciftag = 'no'

        if eachLine.find('loop_') > -1 or eachLine.find('_chem') > -1 \
           or eachLine.find('data_') > -1 or eachLine.find('#') > -1:
            ciftag = 'yes'

        #################
        # 0. ATOM LISTS #
        #################

        if ciftag == 'yes' and extract_atom_header == 'done':
            read_atom = 'done'

        if ciftag == 'no':

            # Turn on reading once header mapping is extracted

            if chem_comp_atom_head == 'yes' and extract_atom_header != 'done':
                extract_atom_header = 'yes'
                read_atom = 'yes'

            # Extract fields for atom list mapping

            if extract_atom_header == 'yes':

                num_chem_comp_atom = len(chem_comp_atom)

                count = 0
                while count < num_chem_comp_atom:

                    tag = chem_comp_atom[count]

                    if tag.find('_chem_comp_atom.atom_id') > -1:
                        atom_id_point = count

                    if tag.find('_chem_comp_atom.type_symbol') > -1:
                        atom_symbol_point = count

                    if tag.find('_chem_comp_atom.x') > -1:
                        atom_x_point = count

                    if tag.find('_chem_comp_atom.y') > -1:
                        atom_y_point = count

                    if tag.find('_chem_comp_atom.z') > -1:
                        atom_z_point = count

                    count = count + 1

                    extract_atom_header = 'done'

            # Extract data into lists using pointers

            if read_atom == 'yes':

                parseLine = eachLine.split()
                atom_id = parseLine[atom_id_point]
                atom_symbol = parseLine[atom_symbol_point]

                aList_atom_id.append(atom_id)
                aList_atom_symbol.append(atom_symbol)

                if atom_x != 'none':
                    atom_x = parseLine[atom_x_point]
                    aList_atom_x.append(atom_x)

                if atom_y != 'none':  
                    atom_y = parseLine[atom_y_point]
                    aList_atom_y.append(atom_y)

                if atom_z != 'none': 
                    atom_z = parseLine[atom_z_point]
                    aList_atom_z.append(atom_z)

        # Identify and collect _chem_comp_atom tags

        if eachLine.find('_chem_comp_atom') > -1:
            eachLine.strip()
            chem_comp_atom.append(eachLine)
            chem_comp_atom_head = 'yes'    

        ##############
        # 1. BONDS   #
        ##############

        if ciftag == 'yes' and extract_bond_header == 'done':
            read_bond = 'done'

        if ciftag == 'no':

            # Turn on reading once header mapping is extracted

            if chem_comp_bond_head == 'yes' and extract_bond_header != 'done':
                extract_bond_header = 'yes'
                read_bond = 'yes'

            # Extract fields for bond list mapping

            if extract_bond_header == 'yes':

                num_chem_comp_bond = len(chem_comp_bond)

                count = 0
                while count < num_chem_comp_bond:

                    tag = chem_comp_bond[count]

                    if tag.find('_chem_comp_bond.atom_id_1') > -1:
                        atom1_id_point = count

                    if tag.find('_chem_comp_bond.atom_id_2') > -1:
                        atom2_id_point = count                    

                    if tag.find('_chem_comp_bond.value_dist') > -1 \
                       and tag.find('value_dist_esd') == -1:
                        atom12_dist_point = count

                    count = count + 1

                    extract_bond_header = 'done'

            # Extract data into lists using pointers

            if read_bond == 'yes':

                parseLine = eachLine.split()
                atom1 = parseLine[atom1_id_point]
                atom2 = parseLine[atom2_id_point]
                atom12_dist = parseLine[atom12_dist_point]

                aList_atom1_bond.append(atom1)
                aList_atom2_bond.append(atom2)
                aList_atom12_dist_bond.append(atom12_dist)

        # Identify and collect _chem_comp_bond tags

        if eachLine.find('_chem_comp_bond') > -1:
            eachLine.strip()
            chem_comp_bond.append(eachLine)
            chem_comp_bond_head = 'yes'

        #############
        # 2. ANGLES #
        #############

        if ciftag == 'yes' and extract_angle_header == 'done':
            read_angle = 'done'

        if ciftag == 'no':

            # Turn on reading once header mapping is extracted

            if chem_comp_angle_head == 'yes' and extract_angle_header != 'done':
                extract_angle_header = 'yes'
                read_angle = 'yes'

            # Extract fields for bond list mapping

            if extract_angle_header == 'yes':

                num_chem_comp_angle = len(chem_comp_angle)

                count = 0
                while count < num_chem_comp_angle:

                    tag = chem_comp_angle[count]

                    if tag.find('_chem_comp_angle.atom_id_1') > -1:
                        atom1_id_point = count

                    if tag.find('_chem_comp_angle.atom_id_2') > -1:
                        atom2_id_point = count

                    if tag.find('_chem_comp_angle.atom_id_3') > -1:
                        atom3_id_point = count                    

                    if tag.find('_chem_comp_angle.value_angle') > -1 \
                       and tag.find('value_angle_esd') == -1:
                        atom13_dist_point = count

                    count = count + 1

                    extract_angle_header = 'done'

            # Extract data into lists using pointers

            if read_angle == 'yes':

                parseLine = eachLine.split()
                atom1 = parseLine[atom1_id_point]
                atom2 = parseLine[atom2_id_point]
                atom3 = parseLine[atom3_id_point]
                atom13_dist = parseLine[atom13_dist_point]

                aList_atom1_angle.append(atom1)
                aList_atom2_angle.append(atom2)
                aList_atom3_angle.append(atom3)
                aList_atom13_dist_angle.append(atom13_dist)

        # Identify and collect _chem_comp_angle tags

        if eachLine.find('_chem_comp_angle') > -1:
            eachLine.strip()
            chem_comp_angle.append(eachLine)
            chem_comp_angle_head = 'yes'

        #####################
        # 3. Chiral centers #
        #####################

        if ciftag == 'yes' and extract_chiral_header == 'done':
            read_chiral = 'done'

        if ciftag == 'no':

            # Turn on reading once header mapping is extracted

            if chem_comp_chiral_head == 'yes' and extract_chiral_header != 'done':
                extract_chiral_header = 'yes'
                read_chiral = 'yes'

            # Extract fields for chiral center list mapping

            if extract_chiral_header == 'yes':

                num_chem_comp_chiral = len(chem_comp_chiral)

                count = 0
                while count < num_chem_comp_chiral:

                    tag = chem_comp_chiral[count]

                    if tag.find('_chem_comp_chir.id') > -1:
                        chiral_id_point = count

                    if tag.find('_chem_comp_chir.atom_id_centre') > -1:
                        chiral_atom_id_point = count

                    if tag.find('_chem_comp_chir.atom_id_1') > -1:
                        chiral_atom_id_1_point = count

                    if tag.find('_chem_comp_chir.atom_id_2') > -1:
                        chiral_atom_id_2_point = count

                    if tag.find('_chem_comp_chir.atom_id_3') > -1:
                        chiral_atom_id_3_point = count

                    if tag.find('_chem_comp_chir.volume_sign') > -1:
                        chiral_sign_point = count

                    count = count + 1

                    extract_chiral_header = 'done'

            # Extract data into lists using pointers

            if read_chiral == 'yes':

                parseLine = eachLine.split()
                chiral_id = parseLine[chiral_id_point]
                chiral_atom_id = parseLine[chiral_atom_id_point]
                chiral_atom_id_1 = parseLine[chiral_atom_id_1_point]
                chiral_atom_id_2 = parseLine[chiral_atom_id_2_point]
                chiral_atom_id_3 = parseLine[chiral_atom_id_3_point]
                chiral_sign = parseLine[chiral_sign_point]

                aList_chiral_id.append(chiral_id)
                aList_chiral_atom.append(chiral_atom_id) 
                aList_chiral_atom1.append(chiral_atom_id_1)
                aList_chiral_atom2.append(chiral_atom_id_2)
                aList_chiral_atom3.append(chiral_atom_id_3)
                aList_chiral_sign.append(chiral_sign)

        # Identify and collect _chem_comp_angle tags

        if eachLine.find('_chem_comp_chir') > -1:
            eachLine.strip()
            chem_comp_chiral.append(eachLine)
            chem_comp_chiral_head = 'yes'

        #############
        # 4. Planes #
        #############

        if ciftag == 'yes' and extract_plane_header == 'done':
            read_plane = 'done'

        if ciftag == 'no':

            # Turn on reading once header mapping is extracted

            if chem_comp_plane_head == 'yes' and extract_plane_header != 'done':
                extract_plane_header = 'yes'
                read_plane = 'yes'

            # Extract fields for bond list mapping

            if extract_plane_header == 'yes':

                num_chem_comp_plane = len(chem_comp_plane)

                count = 0
                while count < num_chem_comp_plane:

                    tag = chem_comp_plane[count]

                    if tag.find('_chem_comp_plane_atom.comp_id') > -1:
                        plane_id_point = count

                    if tag.find('_chem_comp_plane_atom.atom_id') > -1:
                        plane_atom_id_point = count                    

                    count = count + 1

                    extract_plane_header = 'done'

            # Extract data into lists using pointers

            if read_plane == 'yes':

                parseLine = eachLine.split()
                plane_id = parseLine[plane_id_point]
                plane_atom_id = parseLine[plane_atom_id_point]

                aList_plane_id.append(plane_id)
                aList_plane_atom.append(plane_atom_id)

        # Identify and collect _chem_comp_angle tags

        if eachLine.find('_chem_comp_plane_atom') > -1:
            eachLine.strip()
            chem_comp_plane.append(eachLine)
            chem_comp_plane_head = 'yes'

    # Entry list counts

    length_atomlist = len(aList_atom_id)
    length_atomelementlist = len(aList_atom_symbol)
    length_bondlist = len(aList_atom1_bond)
    length_anglelist = len(aList_atom1_angle)
    length_chirallist = len(aList_chiral_id)
    length_planelist = len(aList_plane_id)
    length_xyz = len(aList_atom_x)

    ################################
    # Write restraints SHELX style #
    ################################

    if refprogram == 'shelx' or refprogram == 'SHELX':

        outputfile = os.path.join(workingdir,'mi_shelx.lib')

        print '\nGenerating dictionary data for entity:',entityname
        print 'Job-ID:',job_id
        print 'Writing SHELX restraints file: ',outputfile

        # Check if dictionary file already exists for possible append mode

        fileexists = os.path.exists(outputfile)
        if fileexists != 0:

            print '\nA SHELX dictionary file already exists'

            file = open(outputfile,'r')
            allLines = file.readlines()
            file.close()

            search_entity = 'DFIX_' + entityname

            for eachLine in allLines:           
                if eachLine.find(search_entity) > -1:
                    write_mode = 'exists'

            if write_mode == 'exists':
                print 'Exiting - entity is already in this dictionary'
                time.sleep(4)
                return 1

            if write_mode == 'new':
                print 'Entity will be appended to this current SHELX dictionary'
                write_mode = 'append'

        # Convert the angle restraint to a 1-3 distance restraint

        count = 0
        while count < length_anglelist:

            atom1 = aList_atom1_angle[count]
            atom2 = aList_atom2_angle[count]
            atom3 = aList_atom3_angle[count]

            length12 = 0.0
            length23 = 0.0

            # Find 12 bond length

            count_bond = 0
            while count_bond < length_bondlist:

                if atom1 == aList_atom1_bond[count_bond] and atom2 == aList_atom2_bond[count_bond]:
                    length12 = aList_atom12_dist_bond[count_bond]

                if atom2 == aList_atom1_bond[count_bond] and atom1 == aList_atom2_bond[count_bond]:
                    length12 = aList_atom12_dist_bond[count_bond]

                count_bond = count_bond + 1

            # Find 23 bond length

            count_bond = 0
            while count_bond < length_bondlist:

                if atom2 == aList_atom1_bond[count_bond] and atom3 == aList_atom2_bond[count_bond]:
                    length23 = aList_atom12_dist_bond[count_bond]

                if atom3 == aList_atom1_bond[count_bond] and atom2 == aList_atom2_bond[count_bond]:
                    length23 = aList_atom12_dist_bond[count_bond]

                count_bond = count_bond + 1

            # Compute distance

            if length12 != 0.0 and length23 != 0.0:

                length12 = float(length12)
                length23 = float(length23)

                angle13 = aList_atom13_dist_angle[count]
                angle13 = float(angle13)
                angle13 = angle13/twopi
                cos_angle13 = math.cos(angle13)

                length13 = length12*length12 + length23*length23 + 2.0*length12*length23*cos_angle13
                length13 = math.sqrt(length13)

                length13_round = round(length13,3)
                length13_round = str(length13_round)

                aList_angle13.append(length13_round)

            else:

                print 'Could not evaluate an angle 1-3 distance'
                time.sleep(4)
                return 1

            count = count + 1

        # Write dictionary

        bondtag = 'DFIX_' + entityname
        angletag = 'DANG_' + entityname
        chiraltag = 'CHIV_' + entityname
        planetag = 'FLAT_' + entityname

        if write_mode == 'append':
            file = open('mi_shelx.lib','a')
        else:
            file = open('mi_shelx.lib','w')

        file.write('\n')

        # Write bonds

        count= 0
        while count < length_bondlist:

            bondout = bondtag + ' ' + \
                      aList_atom1_bond[count] + ' ' + aList_atom2_bond[count] + \
                      ' ' + aList_atom12_dist_bond[count]
            file.write(bondout)
            file.write('\n')

            count = count + 1

        file.write('\n')

        # Write angles

        count = 0
        while count < length_anglelist:

            angleout = angletag + ' ' + \
                       aList_atom1_angle[count] + ' ' + aList_atom3_angle[count] + \
                       ' ' + aList_angle13[count]

            file.write(angleout)
            file.write('\n')

            count = count + 1

        file.write('\n')

        # Write chiral centers

        count = 0
        while count < length_chirallist:

            chiralout = chiraltag + ' ' + \
                        aList_chiral_atom[count]

            file.write(chiralout)
            file.write('\n')

            count = count + 1

        file.write('\n')

        # Write planes

        prev_plane_id = 'no_id'

        count = 0
        while count < length_planelist:

            plane_id = aList_plane_id[count]

            if plane_id != prev_plane_id:

                file.write('\n')
                planeout = planetag + ' '
                file.write(planeout)

            atomout = aList_plane_atom[count] + ' '
            file.write(atomout)

            prev_plane_id = plane_id

            count = count + 1

        file.write('\n')
        file.close()

    ##################################
    # Write restraints CNS/CNX style #
    ##################################

    if refprogram == 'cns' or refprogram == 'cnx' or refprogram == 'CNS' or refprogram == 'CNX':

        top_file = entityname + '.top'
        par_file = entityname + '.par'

        outputfile_top = os.path.join(workingdir,top_file)
        outputfile_par = os.path.join(workingdir,par_file)

        print '\nGenerating dictionary data for entity:',entityname
        print 'Job-ID:',job_id
        print 'Writing CNS/CNX topology file:  ',outputfile_top
        print 'Writing CNX/CNX parameter file: ',outputfile_par


        #################################################
        # Build chemical atom type list and mass list   #
        #################################################

        count = 0
        while count < length_atomlist:

            atom_element = aList_atom_symbol[count]
            str_count = str(count)

            # Assign mass and VDV 6-12 for elements found in ligands

            mass = '12.00'
            vdv = '0.1200  3.7418    0.1000  3.7418'

            if atom_element == 'H':
                mass = '1.00'
                vdv = '0.0498 1.4254 0.0498 1.4254'

            if atom_element == 'C':
                mass = '12.00'
                vdv = '0.120 3.7418 0.120 3.7418'

            if atom_element == 'N':
                mass = '14.00'
                vdv = '0.2384 2.8509 0.2384 2.8509'

            if atom_element == 'O':
                mass = '16.00'
                vdv = '0.1591 2.8508 0.1591 2.8508'

            if atom_element == 'F':
                mass = '19.00'
                vdv = '0.10 3.029 0.10 3.029'

            if atom_element == 'P':
                mass = '31.00'
                vdv = '0.5849 3.3854 0.5849 3.3854'

            if atom_element == 'S':
                mass = '32.00'
                vdv = '0.043 3.3676 0.043 3.3676'

            if atom_element == 'CL':
                mass = '35.5'
                vdv = '0.26 3.671 0.26 3.671'

            if atom_element == 'SE':
                mass = '79.0'
                vdv = '0.043 3.510 0.043 3.510'

            if atom_element == 'BR':
                mass = '80.0'
                vdv = '0.320 3.884 0.320 3.884'

            if atom_element == 'I':
                mass = '127.0'
                vdv = '0.80 3.635 0.80 3.635'

            aList_atom_mass.append(mass)
            aList_atom_vdv.append(vdv)

            # Assign a chemical atom type (limited to 4 characters)

            length_element = len(atom_element)

            if length_element == 1:
                out_atom_element = atom_element
                atom_type = atom_element + '_' + str_count
            else:
                out_atom_element = atom_element + str_count

            aList_atom_type.append(atom_type)

            # end of loop

            count = count + 1

        #############################################
        # Establish improper dihedrals for planes  #
        #############################################

        # Build list of quartets of connected atoms

        count1 = 0
        while count1 < length_anglelist:

            atom1 = aList_atom1_angle[count1]
            atom2 = aList_atom2_angle[count1]
            atom3 = aList_atom3_angle[count1]

            count2 = 0
            while count2 < length_anglelist:

                atom4 = aList_atom1_angle[count2]
                atom5 = aList_atom2_angle[count2]
                atom6 = aList_atom3_angle[count2]

                if atom2 == atom4 and atom3 == atom5:
                    aList_quartet.append(atom1)
                    aList_quartet.append(atom2)
                    aList_quartet.append(atom3)
                    aList_quartet.append(atom6)

                if atom2 == atom6 and atom3 == atom5:
                    aList_quartet.append(atom1)
                    aList_quartet.append(atom2)
                    aList_quartet.append(atom3)
                    aList_quartet.append(atom4)

                if atom2 == atom4 and atom1 == atom5:
                    aList_quartet.append(atom3)
                    aList_quartet.append(atom2)
                    aList_quartet.append(atom1)
                    aList_quartet.append(atom6)

                if atom3 == atom5 and atom2 == atom6:
                    aList_quartet.append(atom4)
                    aList_quartet.append(atom5)
                    aList_quartet.append(atom6)
                    aList_quartet.append(atom1)                

                count2 = count2 + 1

            count1 = count1 + 1

        length_quartet = len(aList_quartet) - 4

        # Look to see if each quartet exists in a plane

        match = 0
        count2 = 0
        while count2 < length_quartet:

            # Get a quartet of 4 connected atoms 

            count2_0 = count2
            count2_1 = count2 + 1
            count2_2 = count2 + 2
            count2_3 = count2 + 3

            atom1 = aList_quartet[count2_0]
            atom2 = aList_quartet[count2_1]
            atom3 = aList_quartet[count2_2]
            atom4 = aList_quartet[count2_3]

            # Loop through the planes to see if this quartet is within a plane

            if length_planelist > 1:
                plane_id_prev = aList_plane_id[0]

            match = 0
            count1 = 0
            while count1 < length_planelist:

                plane_id = aList_plane_id[count1]
                atom0 = aList_plane_atom[count1]

                if atom0 == atom1 or atom0 == atom2 or atom0 == atom3 or atom0 == atom4:
                    match = match + 1

                if plane_id != plane_id_prev:
                    match = 0

                if match == 4:
                    improper_out =  ' ' + atom1 + ' ' + atom2 + ' ' + atom3 + ' ' + atom4 + ' '
                    aList_improper_planes.append(improper_out)
                    match = 0

                plane_id_prev = plane_id

                count1 = count1 + 1  

            count2 = count2 + 4

        length_improper_planes = len(aList_improper_planes)

        ########################
        # Write topology file  #
        ########################

        file = open(top_file,'w')

        # Header

        file.write('REMARKS ')
        file.write(top_file)
        file.write('\nREMARKS Created by Molecular Images script mi_convertlib.py\n')
        file.write('\n')
        file.write(' set echo=false end \n')
        file.write('\n')

        # Mass list

        count = 0
        while count < length_atomlist:

            atom_type = aList_atom_type[count]
            atom_mass = aList_atom_mass[count]

            file.write(' MASS ')
            file.write(atom_type)
            file.write('  ')
            file.write(atom_mass)
            file.write('\n')

            count = count + 1

        file.write('\n')
        file.write(' autogenerate angles=true end\n')
        file.write(' \n')
        file.write('RESIdue ')
        file.write(entityname)
        file.write('\n\n')

        # Group assignment list

        file.write('GROUP\n')

        count = 0
        while count < length_atomlist:

            atom_type = aList_atom_type[count]
            atom_name = aList_atom_id[count]

            length_atom_name = len(atom_name)
            if length_atom_name == 3:
                atom_name_out = atom_name + ' '
            if length_atom_name == 2:
                atom_name_out = atom_name + '  '
            if length_atom_name == 1:
                atom_name_out = atom_name + '   '

            file.write(' ATOM ')
            file.write(atom_name_out)
            file.write(' TYPE ')
            file.write(atom_type)
            file.write(' CHARge  0.0 END\n')

            count = count + 1

        # Bond list

        file.write('\n')

        line_count = 0
        count = 0

        while count < length_bondlist:

            atom1 = aList_atom1_bond[count]
            atom2 = aList_atom2_bond[count]

            length_atom1_name = len(atom1)
            if length_atom1_name == 3:
                atom1_name_out = atom1 + ' '
            if length_atom1_name == 2:
                atom1_name_out = atom1 + '  '
            if length_atom1_name == 1:
                atom1_name_out = atom1 + '   '

            length_atom2_name = len(atom2)
            if length_atom2_name == 3:
                atom2_name_out = atom2 + ' '
            if length_atom2_name == 2:
                atom2_name_out = atom2 + '  '
            if length_atom2_name == 1:
                atom2_name_out = atom2 + '   '

            file.write(' BOND ')
            file.write(atom1_name_out)
            file.write(' ')
            file.write(atom2_name_out)
            file.write('   ')

            line_count = line_count + 1

            if line_count == 4:
                file.write('\n')
                line_count = 0

            count = count + 1

        # IMPRoper list for chirality

        file.write('\n\n')

        count = 0
        while count < length_chirallist:

            atom0 = aList_chiral_atom[count]
            atom1 = aList_chiral_atom1[count]
            atom2 = aList_chiral_atom2[count]
            atom3 = aList_chiral_atom3[count]

            file.write(' IMPRoper  ')
            file.write(atom0)
            file.write(' ')
            file.write(atom1)
            file.write(' ')
            file.write(atom2)
            file.write(' ')
            file.write(atom3)
            file.write('\n')    

            count = count + 1

        # IMPROPer list for planes

        count = 0
        while count < length_improper_planes:

            improper_out = aList_improper_planes[count]

            file.write(' IMPRoper ')
            file.write(improper_out)
            file.write('\n')

            count = count + 1

        # Close out

        file.write('\nEND { RESIdue ')
        file.write(entityname)
        file.write(' }\n')

        file.close()

        #########################
        # Write parameter file  #
        #########################

        file = open(par_file,'w')

        # Header

        file.write('REMARKS ')
        file.write(par_file)
        file.write('\nREMARKS Created by Molecular Images script mi_convertlib.py\n')
        file.write('\n')
        file.write(' set echo=false end \n')
        file.write('\n')   

        # Write parameters for bond types

        count = 0
        while count < length_bondlist:

            atom1 = aList_atom1_bond[count]
            atom2 = aList_atom2_bond[count]
            bond_distance = aList_atom12_dist_bond[count]

            # Look up the chemical atom type based on atom name

            count1 = 0
            while count1 < length_atomlist:

                atom0 = aList_atom_id[count1]

                if atom0 == atom1:
                    atom1_type = aList_atom_type[count1]

                if atom0 == atom2:
                    atom2_type = aList_atom_type[count1]

                count1 = count1 + 1

            file.write(' BOND ')
            file.write(atom1_type)
            file.write(' ')
            file.write(atom2_type)
            file.write(' ')
            file.write(bond_weight)
            file.write(' ')
            file.write(bond_distance)
            file.write('\n')

            count = count + 1

        file.write('\n')

        # Write parameters for angle types

        count = 0
        while count < length_anglelist:

            atom1 = aList_atom1_angle[count]
            atom2 = aList_atom2_angle[count]
            atom3 = aList_atom3_angle[count]        
            angle_degrees = aList_atom13_dist_angle[count]

            # Look up chemical atom type based on atom name

            count1 = 0
            while count1 < length_atomlist:

                atom0 = aList_atom_id[count1]

                if atom0 == atom1:
                    atom1_type = aList_atom_type[count1]

                if atom0 == atom2:
                    atom2_type = aList_atom_type[count1]

                if atom0 == atom3:
                    atom3_type = aList_atom_type[count1]            

                count1 = count1 + 1

            file.write(' ANGLe ')
            file.write(atom1_type)
            file.write(' ')
            file.write(atom2_type)
            file.write(' ')
            file.write(atom3_type)
            file.write(' ')        
            file.write(angle_weight)
            file.write(' ')
            file.write(angle_degrees)
            file.write('\n')

            count = count + 1

        file.write('\n')

        # Write parameters for chiral center impropers

        count = 0
        while count < length_chirallist:

            atom1 = aList_chiral_atom[count]
            atom2 = aList_chiral_atom1[count]
            atom3 = aList_chiral_atom2[count]
            atom4 = aList_chiral_atom3[count]
            chiral_sign = aList_chiral_sign[count]

            # Look up chemical atom type based on atom name

            count1 = 0
            while count1 < length_atomlist:

                atom0 = aList_atom_id[count1]

                if atom0 == atom1:
                    atom1_type = aList_atom_type[count1]

                    if length_xyz > 0:
                        atom1_x = aList_atom_x[count1]
                        atom1_y = aList_atom_y[count1]
                        atom1_z = aList_atom_z[count1]

                if atom0 == atom2:
                    atom2_type = aList_atom_type[count1]

                    if length_xyz > 0:
                        atom2_x = aList_atom_x[count1]
                        atom2_y = aList_atom_y[count1]
                        atom2_z = aList_atom_z[count1] 

                if atom0 == atom3:
                    atom3_type = aList_atom_type[count1]

                    if length_xyz > 0:
                        atom3_x = aList_atom_x[count1]
                        atom3_y = aList_atom_y[count1]
                        atom3_z = aList_atom_z[count1] 

                if atom0 == atom4:
                    atom4_type = aList_atom_type[count1]

                    if length_xyz > 0:
                        atom4_x = aList_atom_x[count1]
                        atom4_y = aList_atom_y[count1]
                        atom4_z = aList_atom_z[count1] 

                count1 = count1 + 1

            # Establish chirality sign from given 

            if chiral_sign == 'negativ':
                improper_chiral_angle = '-35.264'

            if chiral_sign == 'positiv':
                improper_chiral_angle = '35.264'

            # If chirality not given compute angle between 123 and 234 atom planes from coordinates 

            if chiral_sign == 'none' and length_xyz > 0:

                atom1_x = float(atom1_x)
                atom1_y = float(atom1_y)
                atom1_z = float(atom1_z)
                atom2_x = float(atom2_x)
                atom2_y = float(atom2_y)
                atom2_z = float(atom2_z)
                atom3_x = float(atom3_x)
                atom3_y = float(atom3_y)
                atom3_z = float(atom3_z)
                atom4_x = float(atom4_x)
                atom4_y = float(atom4_y)
                atom4_z = float(atom4_z)

                rx_12 = atom1_x - atom2_x
                ry_12 = atom1_y - atom2_y
                rz_12 = atom1_z - atom2_z
                rx_23 = atom2_x - atom3_x
                ry_23 = atom2_y - atom3_y
                rz_23 = atom2_z - atom3_z
                rx_34 = atom3_x - atom4_x
                ry_34 = atom3_y - atom4_y
                rz_34 = atom3_z - atom4_z

                ax = ry_12*rz_34 - rz_12*ry_23
                ay = rz_12*rx_34 - rx_12*rz_23
                az = rx_12*ry_34 - ry_12*rx_23
                bx = ry_23*rz_34 - ry_34*rz_23
                by = rz_23*rx_34 - rz_34*rx_23
                bz = rx_23*ry_34 - rx_34*ry_23
                cx = ry_23*az - rz_23*ay
                cy = rz_23*ax - rx_23*az
                cz = rx_23*ay - ry_23*ax

                rar = 1.0/math.sqrt(ax*ax + ay*ay + az*az)
                rbr = 1.0/math.sqrt(bx*bx + by*by + bz*bz)
                rcr = 1.0/math.sqrt(cx*cx + cy*cy + cz*cz)

                ax = ax*rar
                ay = ay*rar
                az = az*rar
                bx = bx*rbr
                by = by*rbr
                bz = bz*rbr
                cx = cx*rcr
                cy = cy*rcr
                cz = cz*rcr

                cp = ax*bx + ay*by + az*bz
                sp = cx*bx + cy*by + cz*bz

                # angle is acos(cp), sign convention is opposite of sign of sp

                if sp < 0.0:
                    chiral_sign = 'positiv'
                    improper_chiral_angle = '35.264'
                else:
                    chiral_sign = 'negativ'
                    improper_chiral_angle = '-35.264'

            if chiral_sign == 'none':
                print 'WARNING - sign of chiral center could not be determined'
                print '+/- 35.264 ambiguity in CNS/CNX dictionary file'

            # Write

            file.write(' DIHEdral ')
            file.write(atom1_type)
            file.write(' ')
            file.write(atom2_type)
            file.write(' ')
            file.write(atom3_type)
            file.write(' ')
            file.write(atom4_type)
            file.write(' ')       
            file.write(improper_weight)
            file.write(' ')
            file.write(improper_chiral_angle)
            file.write('\n')

            count = count + 1

        # Write plane impropers here

        improper_plane_angle = '0.0'

        count = 0
        while count < length_improper_planes:

            improper_quartet = aList_improper_planes[count]
            parseLine = improper_quartet.split()

            atom1 = parseLine[0]
            atom2 = parseLine[1]
            atom3 = parseLine[2]
            atom4 = parseLine[3]

            # Look up chemical atom type based on atom name

            count1 = 0
            while count1 < length_atomlist:

                atom0 = aList_atom_id[count1]

                if atom0 == atom1:
                    atom1_type = aList_atom_type[count1]

                if atom0 == atom2:
                    atom2_type = aList_atom_type[count1]

                if atom0 == atom3:
                    atom3_type = aList_atom_type[count1]            

                if atom0 == atom4:
                    atom4_type = aList_atom_type[count1] 

                count1 = count1 + 1

            file.write(' DIHEdral ')
            file.write(atom1_type)
            file.write(' ')
            file.write(atom2_type)
            file.write(' ')
            file.write(atom3_type)
            file.write(' ')
            file.write(atom4_type)
            file.write(' ')       
            file.write(improper_weight)
            file.write(' ')
            file.write(improper_plane_angle)
            file.write('\n')

            count = count + 1

        file.write('\n')

        # Write non-bonded (VDV) parameters

        count = 0
        while count < length_atomlist:

            atom_type = aList_atom_type[count]
            atom_vdv = aList_atom_vdv[count]

            file.write(' NONBonded ')
            file.write(atom_type)
            file.write(' ')
            file.write(atom_vdv)
            file.write('\n')

            count = count + 1

        # Finish

        file.write('\nset echo=true end\n')
        file.close()

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
    file.write('\nInput mmCIF library: ')
    file.write(ciffile)

    if refprogram == 'shelx' or refprogram == 'SHELX':
        file.write('\nOutput SHELX library: ')
        file.write(outputfile)

        if write_mode == 'new':
            file.write('\nWrite mode: new\n')
        else:
            file.write('\nWrite mode: append\n')

    if refprogram == 'cns' or refprogram == 'cnx' or refprogram == 'CNS' or refprogram == 'CNX':
        file.write('\nOutput CNS/CNX topology file: ')
        file.write(outputfile_top)
        file.write('\nOutput CNS/CNX parameter file: ')
        file.write(outputfile_par)
        file.write('\nWrite mode: new\n')

    file.write('---------------\n')
    file.close()

    #
    return 0

if __name__ == "__main__":
    sys.exit(Run())
