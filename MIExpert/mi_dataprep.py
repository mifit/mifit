#################################################################
#                                                               #
# Refinement data preparation script                            #
#                                                               #
# Copyright: Molecular Images   2005                            #
#                                                               #
# This script is distributed under the same conditions as MIFit #
#                                                               #
#################################################################

import sys
import os
import time
import string
import getopt 
import ccp4check

def Usage():
    print "Usage: %s [Parameters]" % sys.argv[0]
    print "Parameters:"
    print "  -h,--hklin=FILE             input HKL file"
    print "  -d,--workdir=DIR            working directory to use"
    print "  -g,--spacegroup=NUM         spacegroup number, default from hkl file"
    print "  -m,--reference_mtz=MTZFILE  the reference mtz file"
    print "  -?,--help                   for this help"


def Run(argv=None):
    if argv is None:
        argv=sys.argv
    # Initialize

    hklin = 'none'
    mtzout = 'none'
    workingdir = 'none'
    spacegroup_no = 'none'
    referencemtz = 'none'

    ilabel = 'none'
    sigilabel = 'none'
    flabel = 'none'
    sigflabel = 'none'
    rfreelabel = 'none'
    spacegroup_ip = 'none'
    wavelength_ip = 'none'

    temp_refdata = 'mi_refdata.mtz'
    temp_refsort = 'mi_reference_sort.mtz'
    temp_data_name = 'mi_data.hkl'
    temp_mtz_name = 'mi_data.mtz'
    temp_mtz_fname = 'mi_fdata.mtz'
    temp_mtz_sortname = 'mi_sortdata.mtz'
    temp_mtz_final = 'mi_finaldata.mtz'
    temp_mtz_add = 'mi_add.mtz'
    temp_combined_sort = 'mi_data_combined.mtz'
    temp_scaled = 'mi_data_scaled.mtz'

    aList = []
    labelList = []
    colList = []
    reindex_hkl = []
    rfactorList = []
    list_length = 0
    runid = '1'
    runid_int = 0
    projectlog = 'project_history.txt'
    job_prefix = 'dataprep_'
    reindex_flag = 'no'
    output_reindexed = 'no'
    anisotropic_data_flag = 'no'
    spacegroup_name = 'none'

    ###########################################################
    # Read local input file for info on files and key params  #
    ###########################################################
    # Read args
    number_of_args = len(argv)
    args = argv[1:]
    optlist, args = getopt.getopt(
        args,'h:d:g:r:?',
        ['hklin=','workdir=','spacegroup=','reference_mtz=','help'])
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
            if arg_value == '-h' or arg_value=='--hklin':
                hklin = param_value
            elif arg_value == '-d' or arg_value=='--workdir':
                workingdir = param_value
            elif arg_value == '-g' or arg_value=='--spacegroup':
                spacegroup_no = param_value
            elif arg_value == '-m' or arg_value=='--reference_mtz':
                referencemtz = param_value
        count=count+1

    ccp4,error = ccp4check.ccp4check()
    if not ccp4:
      print '\n' + error + '\n'
      time.sleep(4)
      return 1

    fileexists = os.path.exists(hklin)
    if fileexists == 0:
        print 'The integrated intensity data file was not found ',hklin
        time.sleep(4)
        return 1

    fileexists = os.path.exists(workingdir)
    if fileexists == 0:
        print 'The working directory was not found ',workingdir
        time.sleep(4)
        return 1

    fileexists = os.path.exists(referencemtz)
    if fileexists == 0 and referencemtz != 'none':
        print 'The reference mtz file was not found ',referencemtz
        time.sleep(4)
        return 1

    os.chdir(workingdir)
    pwd = os.getcwd()

    # Identify input file type

    if hklin.find('.mtz') > -1:

        # Extension identifies SCALA type

        datatype = 'scala'

        file = open(hklin,'rb')
        allLines = file.readlines()
        file.close()

        file = open(temp_mtz_name,'wb')
        file.writelines(allLines)
        file.close()

    else:

        # Check file for d*trek or SCALEPACK type by finding d*trek keywords

        file = open(hklin,'r')
        allLines = file.readlines()
        file.close()

        for eachLine in allLines:

            if eachLine.find('CRYSTAL_SPACEGROUP') > -1:
                aList = eachLine.split('=')
                spacegroup_ip = aList[1]
                spacegroup_ip = spacegroup_ip.replace(';','')

            if eachLine.find('SOURCE_WAVELENGTH') > -1:
                aList = eachLine.split()
                wavelength_ip = aList[2]
                wavelength_ip = wavelength_ip.replace(';','')

        if spacegroup_ip != 'none':    
            datatype = 'dstartrek'
        else:
            datatype = 'scalepack'

            # Extract spacegroup name

            count = 0
            for eachLine in allLines:

                count = count + 1

                if count == 3:
                    aLine = eachLine.split()
                    line_length = len(aLine)
                    if line_length == 7:
                        spacegroup_name = aLine[6]
                        spacegroup_name = spacegroup_name.upper()

        file = open(temp_data_name,'w')
        file.writelines(allLines)
        file.close() 

    # Check space group dependencies of input formats

    if datatype == 'dstartrek' and spacegroup_no == 'none':
        spacegroup_no = spacegroup_ip

    if datatype == 'scalepack' and spacegroup_no == 'none' and spacegroup_name == 'none':
        print '\nA space group name was expected as item six on line three in the SCALEPACK data\n'
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
                aList = eachLine.split('_')
                runid = aList[1]
                runid_int = int(runid)

        runid_int = runid_int + 1
        runid = str(runid_int)

    job_id = job_prefix + runid

    # Set default output mtz based on input hkl. 

    hkl_basename = os.path.basename(hklin)
    aList = hkl_basename.split('.')
    hkl_rootname = aList[0]

    if datatype != 'scala':
        mtzout = hkl_rootname + '.mtz'
    else:
        mtzout = hkl_rootname + '_truncate.mtz'

    mtzout_full = os.path.join(workingdir,mtzout)

    # Report

    print '\nPreparing amplitude refinement data in MTZ format'
    print 'Job-ID:',job_id
    print 'Output data file:',mtzout_full

    if datatype == 'scala':

        print 'Input file appears to be from SCALA'

        ######################################################
        # For SCALA/MTZ no data conversion needed            #
        # assuming default labels from SCALA or just I,sigI  #
        ######################################################

        # Check labels

        file = open('mi_mtzdump.inp','w')
        file.write('HEADER\n')
        file.write('END\n')
        file.close()

        runmtz = 'mtzdump HKLIN ' + temp_mtz_name + ' < mi_mtzdump.inp > mi_mtzdump.log'
        os.system(runmtz)

        file = open('mi_mtzdump.log','r')
        allLines = file.readlines()
        file.close()

        os.remove('mi_mtzdump.log')
        os.remove('mi_mtzdump.inp')

        read_columns = 'no'
        read_labels = 'no'

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

        list_length = len(labelList)
        count = 0

        while count < list_length:
            if labelList[count] == 'J' and ilabel == 'none':
                ilabel = colList[count]

            if labelList[count] == 'Q' and sigilabel == 'none':
                sigilabel = colList[count]

            count = count + 1

        if ilabel == 'none' or sigilabel == 'none':
            print 'MTZ labels for I,sd(I) could not be established'
            time.sleep(4)
            return 1

    if datatype == 'scalepack':

        ##################################
        # Run SCALEPACK2MTZ conversion   #
        ##################################

        print 'Input file appears to be from SCALEPACK'

        filename_inp = 'mi_runscalepack2mtz.inp'
        filename_log = 'mi_scalepack2mtz.log'

        file = open(filename_inp,'w')
        file.write('name project noname crystal 1 dataset 1\n')

        file.write('symm ')
        if spacegroup_no == 'none':
            file.write(spacegroup_name)
        else:
            file.write(spacegroup_no)

        file.write('\n')
        file.write('wave 0.0\n')
        file.write('end\n')
        file.close()

        runscalepack2mtz = 'scalepack2mtz hklin ' + temp_data_name + ' hklout ' + temp_mtz_name + ' < ' + filename_inp + ' > ' + filename_log

        os.system(runscalepack2mtz)

        fileexists = os.path.exists(temp_mtz_name)
        if fileexists != 0:
            os.remove(temp_data_name)
            os.remove(filename_inp)
            os.remove(filename_log)

        else:
            print 'Output mtz file from SCALEPACK2MTZ was not found'
            time.sleep(4)
            return 1    

    if datatype == 'dstartrek':

        ###########################
        # Run D*TREK conversion   #
        ###########################

        print 'Input file appears to be from D*TREK'

        filename_inp = 'mi_rundstartrek.inp'
        filename_log = 'mi_rundstartrek.log'

        file = open(filename_inp,'w')
        file.write('title noname\n')
        file.write('symm ')
        file.write(spacegroup_no)
        file.write('\n')

        if wavelength_ip != 'none':
            file.write('WAVE ')
            file.write(wavelength_ip)
            file.write('\n')

        file.write('end\n')
        file.close()

        rundtrek2mtz = 'dtrek2mtz hklin ' + temp_data_name + ' hklout ' + temp_mtz_name + ' < ' + filename_inp + ' > ' + filename_log

        os.system(rundtrek2mtz)

        fileexists = os.path.exists(temp_mtz_name)
        if fileexists != 0:
            os.remove(temp_data_name)
            os.remove(filename_inp)
            os.remove(filename_log)
        else:
            print 'Output mtz file from DTREK2MTZ was not found'
            time.sleep(4)
            return 1

    ######################################
    # Run TRUNCATE to get F's from I's   #
    ######################################

    print 'Running TRUNCATE to reduce data from I to F'

    filename_inp = 'mi_truncate.inp'
    truncate_log = 'truncate_' + runid + '.log'
    truncate_log_full = os.path.join(pwd,truncate_log)

    fileexists = os.path.exists(temp_mtz_fname)
    if fileexists != 0:
        os.remove(temp_mtz_fname)

    fileexists = os.path.exists(truncate_log)
    if fileexists != 0:
        os.remove(truncate_log)

    file = open(filename_inp,'w')
    file.write('title noname\n')
    file.write('truncate yes\n')
    file.write('NRESIDUE 500\n')

    # If only I,sig(I) were given input use the given labels. Otherwise assume SCALA defaults

    if datatype == 'scala' and list_length == 5:
        file.write('LABIN IMEAN=')
        file.write(ilabel)
        file.write(' SIGIMEAN=')
        file.write(sigilabel)
        file.write('\n')
        file.write('LABOUT F=FP SIGF=SIGFP\n')
    else:
        file.write('LABOUT F=FP SIGF=SIGFP DANO=DANO SIGDANO=SIGDANO\n')

    file.write('NOHARVEST\n')
    file.write('END\n')
    file.close()

    runtruncate = 'truncate hklin ' + temp_mtz_name + ' hklout ' + temp_mtz_fname + ' < ' + filename_inp + ' > ' + truncate_log

    os.system(runtruncate)

    # Determine input data resolution from log 

    read_resolution = 'no'
    resolution_mtz = 'none'

    fileexists = os.path.exists(truncate_log)
    if fileexists != 0:
        file = open(truncate_log,'r')
        allLines = file.readlines()
        file.close()

        for eachLine in allLines:

            line_length = len(eachLine)

            if read_resolution == 'yes' and line_length > 1:
                parseLine = eachLine.split()
                resolution_mtz = parseLine[5]
                read_resolution = 'no'

            if eachLine.find('*  Resolution Range :') > -1:
                read_resolution = 'yes'

    print 'Resolution limit for initial input to TRUNCATE:',resolution_mtz

    # Check that TRUNCATE ran

    fileexists = os.path.exists(temp_mtz_fname)
    if fileexists != 0:
        os.remove(temp_mtz_name)
        os.remove(filename_inp)
    else:
        print 'Output mtz file from TRUNCATE was not found'
        time.sleep(4)
        return 1

    # Check log for space group number and warnings about problem data

    fileexists = os.path.exists(truncate_log)
    if fileexists != 0:
        file = open(truncate_log,'r')
        allLines = file.readlines()
        file.close()

        for eachLine in allLines:

            if eachLine.find('Beware-serious ANISOTROPY') > -1:
                print '\nWarning - data is seriously anisotropic\n'
                anisotropic_data_flag = 'yes'

            if eachLine.find('* Space group =') > -1:
                aLine = eachLine.split('number')
                spacegroup_no = aLine[1]
                spacegroup_no = spacegroup_no.replace(')','')
                spacegroup_no = spacegroup_no.strip()

    #########################################
    # Fix sort/asymmetric unit to standard  #
    #########################################

    print 'Setting standard sort order'

    filename_inp = 'mi_cad.inp'
    filename_log = 'mi_cad.log'

    file = open(filename_inp,'w')

    if datatype == 'scala' and list_length == 5:
        file.write('LABIN FILE_NUMBER 1 E1=FP E2=SIGFP\n')
        file.write('LABOUT FILE_NUMBER 1 E1=FP E2=SIGFP\n')
    else:
        file.write('LABIN FILE_NUMBER 1 E1=FP E2=SIGFP E3=DANO E4=SIGDANO\n')
        file.write('LABOUT FILE_NUMBER 1 E1=FP E2=SIGFP E3=DANO E4=SIGDANO\n')

    file.write('SORT H K L \n')
    file.write('END\n')
    file.close()

    runcad = 'cad HKLIN1 ' + temp_mtz_fname + ' HKLOUT ' + temp_mtz_sortname + ' < ' + filename_inp + ' > ' + filename_log

    os.system(runcad)

    fileexists = os.path.exists(temp_mtz_sortname)
    if fileexists != 0:
        os.remove(temp_mtz_fname)
        os.remove(filename_inp)
        os.remove(filename_log)
    else:
        print 'Output mtz file from CAD sort was not found'
        time.sleep(4)
        return 1

    ##########################
    # Analyse reference data #
    ##########################

    if referencemtz != 'none':

        print 'Checking reference data'

        # Obtain reference data file locally

        file = open(referencemtz,'rb')
        allLines = file.readlines()
        file.close()

        file = open(temp_refdata,'wb')
        file.writelines(allLines)
        file.close()  

        # Run CAD to force reference data to standard sort/au

        filename_inp = 'mi_cad.inp'
        filename_log = 'mi_cad.log'

        file = open(filename_inp,'w')
        file.write('LABIN FILE 1 ALL\n')
        file.write('SORT H K L \n')
        file.write('END\n')
        file.close()

        runcad = 'cad HKLIN1 ' + temp_refdata + ' HKLOUT ' + temp_refsort + ' < ' + filename_inp + ' > ' + filename_log

        os.system(runcad)

        # Obtain data information

        fileexists = os.path.exists(temp_refsort)
        if fileexists != 0:
            os.remove(temp_refdata)

            file = open(filename_log,'r')
            allLines = file.readlines()
            file.close()        

            read_columns = 'no'
            read_labels = 'no'

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

            list_length = len(labelList)

            count = 0
            while count < list_length:
                if labelList[count] == 'F' and flabel == 'none':
                    flabel = colList[count]

                if labelList[count] == 'Q' and sigflabel == 'none':
                    sigflabel = colList[count]

                if labelList[count] == 'I' and rfreelabel == 'none':
                    rfreelabel = colList[count]

                count = count + 1

            if flabel == 'none' or sigflabel == 'none' or rfreelabel == 'none':
                print 'MTZ labels for F,sd(F),Rfree were not established from reference data'
                time.sleep(4)
                return 1
            else:
                os.remove(filename_inp)
                os.remove(filename_log)

        else:
            print 'Output mtz file from CAD sort of reference was not found'
            time.sleep(4)
            return 1

    #################################
    # Reindexing per reference data #
    #################################

    fileexists = os.path.exists(temp_refsort)
    if fileexists != 0:

        reindex = 'h, k, l'
        reindex_hkl.append(reindex)

        # P4i/related 4i groups  kh-l

        if spacegroup_no == '75' or spacegroup_no == '76' or spacegroup_no == '77' or spacegroup_no == '78' or\
           spacegroup_no == '79' or spacegroup_no == '80':

            reindex = 'k, h, -l'
            reindex_hkl.append(reindex)

        # P3i and R3 groups   -h-kl, kh-l, -k-h-l

        if spacegroup_no == '143' or spacegroup_no == '144' or spacegroup_no == '145' or spacegroup_no == '146':

            reindex = '-h, -k, l'
            reindex_hkl.append(reindex)

            reindex = 'k, h, -l'
            reindex_hkl.append(reindex)

            reindex = '-k, -h, -l'
            reindex_hkl.append(reindex)

        # P3il2  k,h,-l groups

        if spacegroup_no == '149' or spacegroup_no == '151' or spacegroup_no == '153':

            reindex = 'k, h, -l'
            reindex_hkl.append(reindex)

        # P3i2l and R32

        if spacegroup_no == '150' or spacegroup_no == '152' or spacegroup_no == '154' or spacegroup_no == '155':

            reindex = '-h, -k, l'
            reindex_hkl.append(reindex)

        # P6i groups

        if spacegroup_no == '168' or spacegroup_no == '169' or spacegroup_no == '170' or spacegroup_no == '171'\
           or spacegroup_no == '172' or spacegroup_no == '173':

            reindex = 'k, h, -l'
            reindex_hkl.append(reindex)

        # P2i3 and related 2i3 groups

        if spacegroup_no == '195' or spacegroup_no == '196' or spacegroup_no == '197' or spacegroup_no == '198'\
           or spacegroup_no == '199':

            reindex = 'k, h, -l'
            reindex_hkl.append(reindex)

        # Flag for reindexing

        length_reindex_hkl = len(reindex_hkl)
        if length_reindex_hkl > 1:
            reindex_flag = 'yes'

    ########################################
    # Loop through indexing possibilities  #
    ########################################

    if reindex_flag == 'yes':

        print 'Checking for indexing versus reference data'

        count = 0
        while count < length_reindex_hkl:

            # Run REINDEX 

            filename_inp = 'mi_reindex.inp'
            filename_log = 'mi_reindex.log'

            file = open(filename_inp,'w')
            file.write('REINDEX ')
            reindex = reindex_hkl[count]
            file.write(reindex)

            file.write('\nEND\n')
            file.close()

            file_index = count + 1
            file_index = str(file_index)            
            temp_reindexed = 'mi_reindexdata' + file_index + '.mtz'
            temp_reindexed_sort = 'mi_sortdata' + file_index + '.mtz'

            runreindex = 'reindex HKLIN ' + temp_mtz_sortname + ' HKLOUT ' + temp_reindexed + ' < ' + filename_inp + ' > ' + filename_log
            os.system(runreindex)

            os.remove(filename_inp)
            os.remove(filename_log)           

            # Run CAD to get standard sort/au

            fileexists = os.path.exists(temp_reindexed)
            if fileexists != 0:

                filename_inp = 'mi_cad.inp'
                filename_log = 'mi_cad.log'

                file = open(filename_inp,'w')
                file.write('LABIN FILE 1 ALL\n')
                file.write('SORT H K L \n')
                file.write('END\n')
                file.close()

                runcad = 'cad HKLIN1 ' + temp_reindexed + ' HKLOUT ' + temp_reindexed_sort + ' < ' + filename_inp + ' > ' + filename_log
                os.system(runcad)

                os.remove(temp_reindexed)
                os.remove(filename_inp)
                os.remove(filename_log)

            else:

                print 'REINDEX run failed'
                time.sleep(4)
                return 1

            # Run CAD to combine with reference

            fileexists = os.path.exists(temp_reindexed_sort)
            if fileexists != 0:

                filename_inp = 'mi_cad.inp'
                filename_log = 'mi_cad.log'

                file = open(filename_inp,'w')
                file.write('RESOLUTION OVERALL 10.0 5.0\n')
                file.write('LABIN FILE_NUMBER 1 E1=FP E2=SIGFP\n')
                file.write('LABIN FILE_NUMBER 2 E1=')
                file.write(flabel)
                file.write(' E2=')
                file.write(sigflabel)
                file.write('\nLABOUT FILE_NUMBER 2 E1=FREF E2=SIGFREF\n')
                file.write('END\n')
                file.close()

                runcad = 'cad HKLIN1 ' + temp_reindexed_sort + ' HKLIN2 ' + temp_refsort + ' HKLOUT ' + temp_combined_sort + \
                         ' < ' + filename_inp + ' > ' + filename_log

                os.system(runcad)

                os.remove(filename_inp)
                os.remove(filename_log)

            else:

                print 'CAD to sort reindexed data failed'
                time.sleep(4)
                return 1           

            # Run SCALEIT to get agreement

            fileexists = os.path.exists(temp_combined_sort)
            if fileexists != 0:

                filename_inp = 'mi_scaleit.inp'
                filename_log = 'mi_scaleit.log'

                file = open(filename_inp,'w')
                file.write('LABIN FP=FP SIGFP=SIGFP FPH1=FREF SIGFPH1=SIGFREF\n')
                file.write('END\n')
                file.close()

                runscaleit = 'scaleit HKLIN ' + temp_combined_sort + ' HKLOUT ' + temp_scaled + ' < ' + filename_inp \
                             + ' > ' + filename_log

                os.system(runscaleit)

            else:

                print 'CAD to combine reindexed data with reference failed'
                time.sleep(4)
                return 1

            # Parse SCALEIT log

            fileexists = os.path.exists(temp_scaled)
            if fileexists != 0:

                os.remove(temp_scaled)
                os.remove(temp_combined_sort)

                file = open(filename_log,'r')
                allLines = file.readlines()
                file.close()

                os.remove(filename_inp)
                os.remove(filename_log)

                for eachLine in allLines:

                    if eachLine.find('THE TOTALS') > -1:
                        aList = eachLine.split()

                        rfactor = aList[7]
                        rfactorList.append(rfactor)

            else:

                print 'SCALEIT to test data set agreement failed'
                time.sleep(4)
                return 1

            # End of loop over indexing possibilities

            count = count + 1

        # Find the best indexing solution and switch to correct sortdata file to o/p

        length_rfactor = len(rfactorList)
        rfactor_prev = 999.00

        count = 0
        while count < length_rfactor:

            rfactor = rfactorList[count]
            rfactor = float(rfactor)

            if rfactor < rfactor_prev:
                index_solution = count + 1
                rfactor_prev = rfactor 

            count = count + 1

        # Report

        file_index = str(index_solution)

        print 'Using reindexed solution ',file_index,' with R(agree) 10-5A:',rfactor_prev
        temp_mtz_sortname = 'mi_sortdata' + file_index + '.mtz'
        output_reindexed = 'yes'

    ########################
    # Attach R-free flags  #
    ########################

    # If reference data is available combine with those Rfree flags

    if referencemtz != 'none':

        print 'Using reference data to establish Rfree flags'

        filename_inp = 'mi_cad.inp'
        filename_log = 'mi_cad.log'

        file = open(filename_inp,'w')
        file.write('labi file 1 E1=FP E2=SIGFP\n')
        file.write('labi file 2 E1=')
        file.write(rfreelabel)
        file.write('\nEND\n')
        file.close()

        runcad = 'cad hklin1 ' + temp_mtz_sortname + ' hklin2 ' + temp_refsort + ' hklout ' + temp_mtz_add + ' < ' + filename_inp + ' > ' + filename_log

        os.system(runcad)  

        fileexists = os.path.exists(temp_mtz_add)
        if fileexists != 0:
            os.remove(temp_mtz_sortname)
            os.remove(temp_refsort)
            os.remove(filename_inp)
            os.remove(filename_log)
            os.rename(temp_mtz_add,temp_mtz_sortname)
        else:
            print 'Output mtz file from CAD combine was not found'
            time.sleep(4)
            return 1

    else:

        print 'Establishing Rfree flags'

    # Attach flags

    filename_inp = 'mi_rfree.inp'
    filename_log = 'mi_rfree.log'

    file = open(filename_inp,'w')

    if referencemtz != 'none':
        file.write('COMPLETE FREE=FreeR_flag\n')

    file.write('FREERFRAC 0.05\n')
    file.write('END\n')
    file.close()

    runfreer = 'freerflag hklin ' + temp_mtz_sortname + ' hklout ' + temp_mtz_final + ' < ' + filename_inp + ' > ' + filename_log

    os.system(runfreer)

    fileexists = os.path.exists(temp_mtz_final)
    if fileexists != 0:
        os.remove(temp_mtz_sortname)
        os.remove(filename_inp)
        os.remove(filename_log)
    else:
        print 'Output mtz file from FREERFLAG was not found'
        time.sleep(4)
        return 1

    fileexists = os.path.exists(mtzout)
    if fileexists != 0:
        os.remove(mtzout)

    os.rename(temp_mtz_final,mtzout)

    # Clean up

    count = 0
    while count < 4:

        file_index = count + 1
        file_index = str(file_index)            
        temp_reindexed_sort = 'mi_sortdata' + file_index + '.mtz'

        fileexists = os.path.exists(temp_reindexed_sort)
        if fileexists != 0:
            os.remove(temp_reindexed_sort)

        count = count + 1

    ########################
    # Append project log   #
    ########################

    print 'Writing project log'

    runtime = time.ctime(time.time())

    file = open(projectlog,'a')
    file.seek(0,2)
    file.write('Job ID: ')
    file.write(job_id)
    file.write('\nDate: ')
    file.write(runtime)
    file.write('\nInput data: ')
    file.write(hklin)
    file.write('\nInput reference data: ')
    file.write(referencemtz)
    file.write('\nOutput TRUNCATE log: ')
    file.write(truncate_log_full)
    file.write('\nOutput mtz data: ')
    file.write(mtzout_full)
    file.write('\nSpace group number: ')
    file.write(spacegroup_no)
    file.write('\nSummary:')
    if output_reindexed == 'no':
        if anisotropic_data_flag == 'no':
            file.write('no reindexing done\n')
        else:
            file.write('highly anisotropic data, no reindexing done\n')
    else:
        if anisotropic_data_flag == 'no':    
            file.write('reindexed to reference\n')
        else:
            file.write('highly anisotropic data, reindexed to reference\n')        

    file.write('---------------\n')
    file.close()

    #
    return 0

if __name__ == "__main__":
    sys.exit(Run())

