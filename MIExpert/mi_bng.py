######################################################################
#                                                                    #
# Runs 'bind-n-grind' automation including refinement data           #
# preparation, MR, refinement, rigid-ligand fitting and MIFit launch #
#                                                                    #
# Copyright: Molecular Images   2005                                 #
#                                                                    #
# This script is distributed under the same conditions as MIFit      #
#                                                                    #
######################################################################

import sys
import os
import time
import string
import dircache
import getopt
import ccp4check

def Usage():
    print "Usage: %s [options]" % sys.argv[0]
    print "Options are:"
    print "  --hklin=FILE              Data file (option may be repeated)"
    print "  --workdir=DIR             Working dir (option may be repeated)"
    print "  --spacegroup_no=NUM       The spacegroup number to use"
    print "  --reference_mtz=FILE      "
    print "  --pdbin=FILE              The input coordinate file."
    print "  --multi_search=yes or no  Try all pdb files in pdbin's directory? Default: no"   
    print "  --libfile                 The input library file.  Default: no file"
    print "  --fragfit=FILE            The fragfit file. Default: no file"
    print "  --chemfit=FILE            The chemistry structure file. Default: no file"    
    print "  --frag_center=\"x y z\"     Center of the fragment. Default: none"
    print "  --mlwfile=FILE            Input session file for viewpoint. Default: none"
    print "  --pdbviewfile=FILE        Input PDB marker file for viewpoint. Default: none"
    print "  --mifit=yes or no         Launch mifit when done.  Default: no"
    print "  --bngsummary=DIR          Path for bng summary html file. Default: no dir, do not produce"
    print "  --sg_search=yes or no     Do a spacegroup search.  Default: no"
    print "  --molimagehome=DIR        Path to MIFit."
    print "  --writemap=yes or no      Write map around target point.  Default: no"
    print "  --detector_constants=FILE Detector constants file.  Default: no file"
    print "  --process_engine=type     One of dstartrek or mosflm.  Default: dstartrek"
    print "  -?,--help                 This help."
    print ""
    print "Note:"
    print "  The number of hklin options should equal the number of workingdir options"

def Run(argv=None):
    if argv is None:
        argv=sys.argv

    # Initialization and defaults

    integrate = 'mi_runintegrate.txt'
    dataprep = 'mi_rundataprep.txt'
    molrep = 'mi_runmr.txt'
    refine = 'mi_runrefine.txt'
    inputfile = 'mi_runbng.txt'

    mifit_root = 'none'
    symlibinstall = 'none'

    mlw_file = 'bng_milaunch.mlw'
    phs_file = 'bng_mlmap.phs'
    bngsummaryname = 'bng_jobsummary_1.htm'
    bngsummaryroot = 'bng_jobsummary_'

    hklin = 'none'
    workingdir = 'none'
    spacegroup_no = 'none'
    referencemtz = 'none'
    pdbin = 'none'
    multi_search = 'no'    
    libfile = 'none'
    fragfile = 'none'
    chemfile = 'none'
    fragcenter = 'none'
    mlwfilein = 'none'
    pdbviewfile = 'none'
    launch_mifit = 'no'
    bngsummary = 'none'
    pdbfileout = 'none'
    sg_search = 'no'
    write_map = 'no'
    detector_constants = 'none'
    mr_sg = 'none'
    mr_sg_best = 'none'
    process_engine = 'dstartrek'

    fragview = '1.0000 0.0000 0.0000 0.0000 1.0000 0.0000 0.0000 0.0000 1.0000'
    zoom = '30.00'
    frontclip = '3.00'
    backclip = '-3.00'
    contourradius = '12.000000'
    aList_contourlevels1 = []
    aList_contourleveldefault1 = []
    aList_contourcolor1 = []
    aList_color1 = []
    aList_contourlevels2 = []
    aList_contourleveldefault2 = []
    aList_contourcolor2 = []
    aList_color2 = []
    first_map = 'no'
    second_map = 'no'
    border = 10.0

    water_radius = 15.0
    water_radius = water_radius * water_radius

    aList = []
    aList_hklin = []
    aList_workingdir = []
    aList_dir = []
    aList_sg = []

    quote = """'"""

    # parse args
    number_of_args = len(argv)
    args = argv[1:]
    optlist, args = getopt.getopt(
        args,'?',
        ['hklin=','workdir=','spacegroup_no=','reference_mtz=',
         'pdbin=','multi_search=','libfile=','fragfit=','chemfit=',
         'frag_center=','mlwfile=','pdbviewfile=','mifit=','bngsummary=',
         'sg_search=','molimagehome=','writemap=',
         'detector_constants=','process_engine=','help'])
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
            if arg_value=='--hklin':
                hklin = param_value
                aList_hklin.append(hklin)
            if arg_value=='--workdir':
                workingdir = param_value
                aList_workingdir.append(workingdir)
            if arg_value=='--spacegroup_no':
                spacegroup_no = param_value
            if arg_value=='--reference_mtz':
                referencemtz = param_value
            if arg_value=='--pdbin':
                pdbin = param_value
            if arg_value=='--multi_search':
                multi_search = param_value
            if arg_value=='--libfile':
                libfile = param_value
            if arg_value=='--fragfit':
                fragfile = param_value
            if arg_value=='--chemfit':
                chemfile = param_value
            if arg_value=='--frag_center':
                fragcenter = param_value
            if arg_value=='--mlwfile':
                mlwfilein = param_value
            if arg_value=='--pdbviewfile':
                pdbviewfile = param_value                
            if arg_value=='--mifit':
                launch_mifit = param_value
            if arg_value=='--bngsummary':
                bngsummary = param_value
            if arg_value=='--sg_search':
                sg_search = param_value
            if arg_value=='--molimagehome':
                mifit_root = param_value
            if arg_value=='--writemap':
                write_map = param_value
            if arg_value=='--detector_constants':
                detector_constants = param_value
            if arg_value=='--process_engine':
                process_engine = param_value
        count=count+1

    ccp4,error = ccp4check.ccp4check()
    if not ccp4:
      print '\n' + error + '\n'
      time.sleep(4)
      return 1

    # Capture initial CCP4 scratch space because subscripts may reset it

    ccp4_scratch = os.environ['CCP4_SCR']

    # Check MIFit installation to access various files (default and direct path)

    if mifit_root == 'none':
        print '\nMIFit root directory was not given\n'
        time.sleep(4)
        return 1
        
    if not os.path.exists(mifit_root):
        print '\nMIFit root directory was not found\n'
        time.sleep(4)
        return 1

    # Set paths to MIFit executable and symmetry library

    mifit_root_data = os.path.join(mifit_root,'data')    
    symlibinstall = os.path.join(mifit_root_data,'symlib')

    test_platform = sys.platform
    if test_platform.find('win') > -1:
        mifitinstall = os.path.join(mifit_root,'MIFit.exe')
    else:
        mifitinstall = os.path.join(mifit_root,'MIFit')

    if not os.path.exists(symlibinstall):
        print '\nThe MI symmetry library was not located\n'
        time.sleep(4)
        return 1

    # Check paired datasets/working directories

    number_datasets = len(aList_hklin)
    number_workingdirs = len(aList_workingdir)

    if number_datasets != number_workingdirs and number_workingdirs != 0:
        print '\nThe number of working directories must equal the number of datasets - stopping !\n'
        time.sleep(4)
        return 1

    runcount = 0
    while runcount < number_datasets:
        hklin = aList_hklin[runcount]

        if not os.path.exists(hklin):
            time.sleep(4)
            print '\nAn input data file was not found',hklin,'\n'
            time.sleep(4)
            return 1

        runcount = runcount + 1

    runcount = 0
    while runcount < number_workingdirs:
        workingdir = aList_workingdir[runcount]

        if not os.path.exists(workingdir):
            time.sleep(4)
            print '\nAn input working directory file was not found',workingdir,'\n'
            time.sleep(4)
            return 1        

        runcount = runcount + 1

    # Check all common input file paths

    if not os.path.exists(pdbin):
        print '\nThe input coordinate file was not found\n'
        time.sleep(4)
        return 1

    if not os.path.exists(libfile) and os.path.basename(libfile) != 'none':
        print '\nThe input library file was not found\n'
        time.sleep(4)
        return 1

    if not os.path.exists(referencemtz) and os.path.basename(referencemtz) != 'none':
        print '\nThe input reference data file was not found\n'
        time.sleep(4)
        return 1

    if not os.path.exists(fragfile) and os.path.basename(fragfile) != 'none':
        print '\nThe input fragment coordinate file was not found\n'
        time.sleep(4)
        return 1

    if not os.path.exists(chemfile) and os.path.basename(chemfile) != 'none':
        print '\nThe input chemistry structure file was not found\n'
        time.sleep(4)
        return 1

    if not os.path.exists(mlwfilein) and os.path.basename(mlwfilein) != 'none':
        print '\nThe input session file was not found\n'
        time.sleep(4)
        return 1

    if not os.path.exists(bngsummary) and os.path.basename(bngsummary) != 'none':
        print '\nThe directory for the BNG job summary was not found\n'
        time.sleep(4)
        return 1

    # Check the image processing engine parameter is recognized

    if process_engine != 'none' and process_engine != 'dstartrek' and process_engine != 'mosflm':
        print '\nThe image processing engine must be set to one of none/dstartrek/mosflm'
        time.sleep(4)
        return 1

    # Check and obtain fragment/view center

    if fragcenter != 'none':

        aList = fragcenter.split()
        number_args = len(aList)

        if number_args != 3:
            print '\nThe fragment center must be three numbers for x,y,z\n'
            time.sleep(4)
            return 1
        else:
            x_center = aList[0]
            y_center = aList[1]
            z_center = aList[2]
            x_center = float(x_center)
            y_center = float(y_center)
            z_center = float(z_center)

    # Or parse view view point from first coordinate in PDB marker file (ATOM/HETATM record)

    if pdbviewfile != 'none':

        file = open(pdbviewfile,'r')
        allLines = file.readlines()
        file.close()

        found_marker = 'no'

        for eachLine in allLines:

            tag = eachLine[0:6]
            tag = tag.strip()

            if tag == 'ATOM' or tag == 'HETATM':
                if found_marker == 'no':

                    x_center = eachLine[30:38]
                    y_center = eachLine[38:46]
                    z_center = eachLine[46:54]
                    fragcenter = x_center + ' ' + y_center + ' ' + z_center
                    x_center = float(x_center)
                    y_center = float(y_center)
                    z_center = float(z_center)

                    found_marker = 'yes'

    ##########################################
    # Parse view definitions from mlw file   #
    ##########################################

    if mlwfilein != 'none':

        file = open(mlwfilein,'r')
        allLines = file.readlines()
        file.close()

        for eachLine in allLines:

            # Get view point

            if eachLine.find('translation') > -1:
                aLine = eachLine.split()
                number_args = len(aLine)

                if number_args == 4:
                    x_center = aLine[1]
                    y_center = aLine[2]
                    z_center = aLine[3]
                    fragcenter = x_center + ' ' + y_center + ' ' + z_center
                    x_center = float(x_center)
                    y_center = float(y_center)
                    z_center = float(z_center)

            if eachLine.find('rotation') > -1:
                aLine = eachLine.split()
                number_args = len(aLine)

                if number_args == 10:
                    fragview = aLine[1] + ' ' + aLine[2] + ' ' + aLine[3] + ' ' \
                                + aLine[4] + ' ' + aLine[5] + ' ' + aLine[6] + ' ' \
                                + aLine[7] + ' ' + aLine[8] + ' ' + aLine[9]

            if eachLine.find('zoom') > -1:
                aLine = eachLine.split()
                number_args = len(aLine)

                if number_args == 2:
                    zoom = aLine[1]

            if eachLine.find('frontclip') > -1:
                aLine = eachLine.split()
                number_args = len(aLine)

                if number_args == 2:
                    frontclip = aLine[1]

            if eachLine.find('backclip') > -1:
                aLine = eachLine.split()
                number_args = len(aLine)

                if number_args == 2:
                    backclip = aLine[1]


            if eachLine.find('contourradius') > -1:
                aLine = eachLine.split()
                number_args = len(aLine)

                if number_args == 2:
                    contourradius = aLine[1]

            # Get colors and contors

            if eachLine.find('maptocont 1') > -1:
                first_map = 'yes'

            if eachLine.find('maptocont 2') > -1:
                second_map = 'yes'

            if second_map == 'no':

                if eachLine.find('contourlevels') > -1:
                    aList = eachLine.split()
                    aList_contourlevels1.append(aList[1])

                if eachLine.find('contourleveldefault') > -1:
                    aLine = eachLine[19:100]
                    aLine = aLine.strip()
                    aList_contourleveldefault1.append(aLine)

                if eachLine.find('contourcolor') > -1:
                    aList = eachLine.split()
                    aList_contourcolor1.append(aList[1])

                if eachLine.find('color') > -1:
                    aList = eachLine.split()
                    if aList[0] == 'color':
                        aList_color1.append(aList[1])

            if second_map == 'yes':

                if eachLine.find('contourlevels') > -1:
                    aList = eachLine.split()
                    aList_contourlevels2.append(aList[1])

                if eachLine.find('contourleveldefault') > -1:
                    aLine = eachLine[19:100]
                    aLine = aLine.strip()
                    aList_contourleveldefault2.append(aLine)

                if eachLine.find('contourcolor') > -1:
                    aList = eachLine.split()
                    aList_contourcolor2.append(aList[1])

                if eachLine.find('color') > -1:
                    aList = eachLine.split()
                    if aList[0] == 'color':
                        aList_color2.append(aList[1])

    # Check that for fragment input a center is also given

    if fragfile != 'none' and fragcenter == 'none':
        print '\nA fragment center must be available if a ligand file is given\n'
        time.sleep(4)
        return 1

    ##############################    
    # Prepare summary HTML file  #
    ##############################

    file_tag_max = 0

    if os.path.basename(bngsummary) != 'none':

        bngsummaryfile = os.path.join(bngsummary,bngsummaryname)

        # If file already exists increment name

        if os.path.exists(bngsummaryfile):

            aList_dir = os.listdir(bngsummary)
            number_files = len(aList_dir)

            count = 0
            while count < number_files:

                test_file = aList_dir[count]

                if test_file.find(bngsummaryroot) > -1:

                    aList = test_file.split('_')

                    number_args = len(aList)

                    if number_args > 2:

                        file_tag = aList[2]
                        file_tag = file_tag.replace('.htm','')

                        if file_tag.isdigit() == 1:
                            file_tag = int(file_tag)

                            if file_tag > file_tag_max:
                                file_tag_max = file_tag

                count = count + 1

            file_tag_max = int(file_tag_max)
            file_tag_max = file_tag_max + 1
            file_tag_max = str(file_tag_max)

            bngsummaryname = bngsummaryroot + file_tag_max + '.htm'
            bngsummaryfile = os.path.join(bngsummary,bngsummaryname)

        print '\nRUN SUMMARY DATA LOG:',bngsummaryfile

        runtime = time.ctime(time.time())

        # HTML header

        filename = bngsummaryfile
        file = open(filename,'w')

        file.write('<html>\n')
        file.write('<head><title>BNG Run Summary</title></head>\n')
        file.write('<body bgcolor = "white">\n')
        file.write('<h2><center>BNG Run Summary</center></h2>\n')
        file.write('<p>\n')
        file.write('<b>Job start: ')
        file.write(runtime)
        file.write('</b>\n')
        file.write('<p>\n')

        file.write('<table border=1>\n')

        file.write('<tr>\n')
        file.write('<tr bgcolor = "yellow">\n')

        file.write('<td>History</td>')
        file.write('<td>Res. (&#197)</td>')
        file.write('<td>R<sub>work</sub></td>')
        file.write('<td>R<sub>free</sub></td>')
        file.write('<td>Error List</td>')
        file.write('<td>Working Directory</td>')
        file.write('</tr>\n')

        file.write('<tr>\n')

        file.close()

    ######################################################
    # Loop over multiple dataset/working directory pairs #
    ######################################################

    runcount = 0
    while runcount < number_datasets:

        list_runcount = runcount + 1
        list_runcount = str(list_runcount)

        print '\nSTARTING BNG PROCESS ON DATASET',list_runcount
        print '\nUsing MIFit installation:',mifit_root

        image_data_processed = 'yes'

        hklin = aList_hklin[runcount]

        # Default working directory path is dataset location

        if number_workingdirs == 0:
            workingdir = os.path.dirname(hklin)
        else:
            workingdir = aList_workingdir[runcount]            

        dir_list = os.listdir(workingdir)
        number_files = len(dir_list)

        os.chdir(workingdir)

        ##############################################################
        # Data processing option (invoked via image file extensions) #
        ##############################################################

        if hklin.find('.img') > -1 or hklin.find('.osc') > -1:

            print '\nDATA PROCESSING'

            image_data_processed = 'no'

            # Use beam.mask file in image directory or assume Rigaku CCD (only used for d*TREK)

            beam_mask_path = os.path.join(workingdir,'beam.mask')
            if not os.path.exists(beam_mask_path):
                beam_mask_path = 'rigakuccd'

            # Check space group is set

            if spacegroup_no == 'none':
                print '\nThe space group number must be given for image data processing\n'
                time.sleep(4)
                return 1

            # Establish a working directory (BNG) for merged data file and structure solution process

            bng_workingdir = os.path.join(workingdir,'BNG')

            if not os.path.exists(bng_workingdir):
                os.mkdir(bng_workingdir)

            # set args for image data processing
            tmpargs=[]
            tmpargs.append("integrate")
            tmpargs.append("--template_image")
            tmpargs.append(hklin)            
            tmpargs.append("--spacegroup")
            tmpargs.append(spacegroup_no)
            tmpargs.append('--workdir')
            tmpargs.append(bng_workingdir)
            if detector_constants != 'none':
                tmpargs.append('--detector_constants')
                tmpargs.append(detector_constants)

            # Execute image data processing 
            if process_engine == 'dstartrek':
                tmpargs.append('--beammask')
                tmpargs.append(beam_mask_path)
                import mi_integrate
                if mi_integrate.Run(tmpargs)!=0:
                    runcount = runcount + 1
                    continue

                hklin = os.path.join(bng_workingdir,'ScalAverage_1.ref')

            elif process_engine == 'mosflm':
                try:
                    import mi_integrate_mosflm
                except:
                    print '\nMOSFLM data processing is not enabled in this release\n'
                    time.sleep(4)
                    return 1   
                if mi_integrate_mosflm.Run(tmpargs)!=0:
                    runcount = runcount + 1
                    continue

                hklin = os.path.join(bng_workingdir,'ScalAverage_1.mtz')            

            # Check image data processing succeeded

            if not os.path.exists(hklin):
                print '\nImage integration for failed\n'
                time.sleep(4)
                runcount = runcount + 1
                continue
            else:
                workingdir = bng_workingdir
                image_data_processed = 'yes'

        ######################
        # Structure Solution #
        ######################

        if image_data_processed == 'yes':

            ################
            # Data import  #
            ################

            print '\nDATA SETUP'

            # set args
            tmpargs=[]
            tmpargs.append("dataprep")
            tmpargs.append('--hklin')
            tmpargs.append(hklin)
            tmpargs.append('--workdir')
            tmpargs.append(workingdir)
            tmpargs.append('--spacegroup')
            tmpargs.append(spacegroup_no)
            tmpargs.append('--reference_mtz')
            tmpargs.append(referencemtz)

            # Execute
            import mi_dataprep
            if mi_dataprep.Run(tmpargs)!=0:
                runcount = runcount + 1
                continue # try next dataset

            ############
            # Run MR   #
            ############

            print '\nMR CALCULATIONS'

            # Select mtz file name from the last data setup job and check space group

            file = open('project_history.txt','r')
            allLines = file.readlines()
            file.close()

            read_mtz_name = 'no'
            mtzout = 'none'
            spacegroup_check = 'none'

            for eachLine in allLines:

                if eachLine.find('Job ID:') > -1 and eachLine.find('dataprep_') > -1:
                    read_mtz_name = 'yes'

                if read_mtz_name == 'yes' and eachLine.find('Output mtz data:') > -1:        
                    mtzout = eachLine[16:200]
                    mtzout = mtzout.strip()
                    mtzout_full=mtzout
                    read_mtz_name = 'no'

                if eachLine.find('Space group number:') > -1:
                    spacegroup_check = eachLine[19:100]
                    spacegroup_check = spacegroup_check.strip()

            if mtzout == 'none':
                print '\nThere is no mtz file to use for MR !\n'
                time.sleep(4)
                return 1

            # Check user options are possible given space group indexing possibilities and available reference data  

            permutable = 'no'
            if spacegroup_check == '75' or spacegroup_check == '76' or spacegroup_check == '77' or spacegroup_check == '78' \
               or spacegroup_check == '79' or spacegroup_check == '80' or spacegroup_check == '143' or spacegroup_check == '144' \
               or spacegroup_check == '145' or spacegroup_check == '146' or spacegroup_check == '149' \
               or spacegroup_check == '151' or spacegroup_check == '153' or spacegroup_check == '150' or spacegroup_check == '152' \
               or spacegroup_check == '154' or spacegroup_check == '155' or spacegroup_check == '168' or spacegroup_check == '169' \
               or spacegroup_check == '170' or spacegroup_check == '171' or spacegroup_check == '172' or spacegroup_check == '173' \
               or spacegroup_check == '195' or spacegroup_check == '196' or spacegroup_check == '197' or spacegroup_check == '198' \
               or spacegroup_check == '199':
                permutable = 'yes'

            if permutable == 'yes' and referencemtz == 'none':
                print '\nWarning: this space group is reindexable!'
                print 'A reference data set is needed for target site water exclusion, ligand-fitting, molecule repositioning\n'

                mlwfilein = 'none'
                fragcenter = 'none'
                fragfile = 'none'

            # Write run file
            tmpargs=[]
            tmpargs.append("molrep")
            tmpargs.append('--pdbfile')
            tmpargs.append(pdbin)
            tmpargs.append('--mtzfile')
            tmpargs.append(mtzout)
            tmpargs.append('--fixed_pdb')
            tmpargs.append('none')
            tmpargs.append('--workdir')
            tmpargs.append(workingdir)
            tmpargs.append('--multi_search')
            tmpargs.append(multi_search)

            tmpargs.append('--match_pdbin')
            if permutable == 'yes' and referencemtz == 'none':
                tmpargs.append('no')
            else:
                tmpargs.append('yes')

            tmpargs.append('--sg_search')
            if sg_search == 'yes':
                tmpargs.append('yes')
            else:
                tmpargs.append('no')       
            # Execute
            import mi_molrep
            if mi_molrep.Run(tmpargs)!=0:
                runcount = runcount + 1
                continue # try next dataset

            os.environ['CCP4_SCR'] = ccp4_scratch

            ##########################
            # Run initial refinement #
            ##########################

            print '\nPRELIMINARY REFINEMENT'

            read_r = 'no'
            aList_pdb = []
            aList_rvalue = []
            pdbin_best_full = 'none'
            pdbfile = 'none'

            # Select the best coordinate file from the project log

            file = open('project_history.txt','r')
            allLines = file.readlines()
            file.close()

            for eachLine in allLines:

                if eachLine.find('Job ID:') > -1 and eachLine.find('molrep_') > -1:
                    read_r = 'yes'
                    
                if eachLine.find('Output atoms:') > -1 and eachLine.find('molrep_') > -1 and read_r=='yes':
                    pdbfile = eachLine[14:]
                    pdbfile = pdbfile.strip()
                    aList_pdb.append(pdbfile)

                if eachLine.find('MR space group:') > -1 and read_r == 'yes':
                    aList = eachLine.split(':')
                    mr_sg = aList[1]
                    aList_sg.append(mr_sg.strip())

                if eachLine.find('Summary:') > -1 and read_r == 'yes':
                    aList = eachLine.split(':')
                    rvalue = aList[1]
                    rvalue = rvalue.replace('R=','')
                    aList_rvalue.append(rvalue)
                    read_r = 'no'

            number_models = len(aList_rvalue)
            rvalue_best = 999.00

            if number_models == 0:
                print '\nThere is no MR model to refine !\n'
                time.sleep(4)
                return 1

            count = 0
            while count < number_models:

                rvalue = aList_rvalue[count]

                if rvalue.isalpha() == 1:
                    print '\nMR R-value in history file is not a number - value is:',rvalue,'\n'
                    time.sleep(4)
                    return 1
                else:
                    rvalue = float(rvalue)

                if rvalue < rvalue_best:
                    pdbin_best = aList_pdb[count]
                    rvalue_best = rvalue
                    mr_sg_best = aList_sg[count]

                count = count + 1

            pdbin_best_full = pdbin_best.strip()

            if rvalue_best > 0.65:
                print '\nMR R-value is too high - stopping !\n'
                time.sleep(4)
                return 1

            # Change the space group number in the data file if space group search option was invoked

            if sg_search == 'yes':

                if os.path.exists('mi_temp_sg_new.mtz'):
                    os.remove('mi_temp_sg_new.mtz')

                file = open(mtzout_full,'rb')
                allLines = file.readlines()
                file.close()

                file = open('mi_temp_sg.mtz','wb')
                file.writelines(allLines)
                file.close()

                filename_inp = 'mi_cad_sg.inp'
                filename_log = 'mi_cad_sg.log'

                file = open(filename_inp,'w')

                file.write('LABIN FILE_NUMBER 1 ALL\n')
                file.write('SYMMETRY ')
                file.write(mr_sg_best)
                file.write('\n')
                file.write('SORT H K L \n')
                file.write('END\n')
                file.close()

                runcad = 'cad HKLIN1 mi_temp_sg.mtz HKLOUT mi_temp_sg_new.mtz < ' + filename_inp + ' > ' + filename_log

                os.system(runcad)

                if os.path.exists('mi_temp_sg_new.mtz'):
                    os.remove(filename_inp)
                    os.remove(filename_log)
                    os.remove('mi_temp_sg.mtz')
                    os.remove(mtzout_full)
                    os.rename('mi_temp_sg_new.mtz',mtzout_full)

            # Write run file

            tmpargs=[]
            tmpargs.append("refine")
            tmpargs.append('--pdbfile')
            tmpargs.append(pdbin_best_full)
            tmpargs.append('--mtzfile')
            tmpargs.append(mtzout_full)
            tmpargs.append('--workdir')
            tmpargs.append(workingdir)
            tmpargs.append('--libfile')
            tmpargs.append(libfile)
            tmpargs.append('--engine')
            tmpargs.append('refmac5')
            tmpargs.append('--weight')
            tmpargs.append('0.1')
            tmpargs.append('--max_res')
            tmpargs.append('none')
            tmpargs.append('--bref_type')
            tmpargs.append('none')
            tmpargs.append('--cycles')
            tmpargs.append('5')
            tmpargs.append('--water_cycles')
            tmpargs.append('0')
            tmpargs.append('--mifithome')
            tmpargs.append(mifit_root)

            # Execute
            import mi_refine
            if mi_refine.Run(tmpargs)!=0:
                runcount = runcount + 1
                continue # try next dataset

            os.environ['CCP4_SCR'] = ccp4_scratch

            #####################################
            # Run refinement with water-picking #
            #####################################

            print '\nEXTENDED REFINEMENT'

            # Select last coordinate file to continue refinement

            pdbfile = 'none'

            file = open('project_history.txt','r')
            allLines = file.readlines()
            file.close()

            for eachLine in allLines:

                if eachLine.find('Output atoms:') > -1:
                    pdbfile = eachLine[13:200]
                    pdbfile = pdbfile.strip()

            if pdbfile == 'none':
                print '\nPDB file for extended refinement was not found !\n'
                time.sleep(4)
                return 1

            # Write run args
            tmpargs=[]
            tmpargs.append("refine")
            tmpargs.append('--pdbfile')
            tmpargs.append(pdbfile)
            tmpargs.append('--mtzfile')
            tmpargs.append(mtzout_full)
            tmpargs.append('--workdir')
            tmpargs.append(workingdir)
            tmpargs.append('--libfile')
            tmpargs.append(libfile)
            tmpargs.append('--engine')
            tmpargs.append('refmac5')
            tmpargs.append('--weight')
            tmpargs.append('0.1')
            tmpargs.append('--max_res')
            tmpargs.append('none')
            tmpargs.append('--bref_type')
            tmpargs.append('none')
            tmpargs.append('--cycles')
            tmpargs.append('5')
            tmpargs.append('--water_cycles')
            tmpargs.append('3')
            tmpargs.append('--mifithome')
            tmpargs.append(mifit_root)

            # Execute
            import mi_refine
            if mi_refine.Run(tmpargs)!=0:
                runcount = runcount + 1
                continue # try next dataset

            os.environ['CCP4_SCR'] = ccp4_scratch

            #############################################################
            # Run water deletion from target site if viewpoint was set  #
            #############################################################

            if fragcenter != 'none':

                print '\nRECOMPUTING MODEL AND MAP DATA WITH TARGET SITE WATERS REMOVED'

                # Find coordinates and delete waters

                pdbfile = 'none'

                file = open('project_history.txt','r')
                allLines = file.readlines()
                file.close()

                for eachLine in allLines:

                    if eachLine.find('Output atoms:') > -1:
                        pdbfile = eachLine[13:200]
                        pdbfile = pdbfile.strip()

                if pdbfile == 'none':
                    print '\nPDB file for omit map calculations was not found !\n'
                    time.sleep(4)
                    return 1

                file = open(pdbfile,'r')
                allLines = file.readlines()
                file.close()

                pdbfile_omit = pdbfile.replace('.pdb','_omit.pdb')
                file = open(pdbfile_omit,'w')

                for eachLine in allLines:

                    tag = eachLine[0:6]
                    tag = tag.strip()

                    write_record = 'yes'

                    if tag == 'ATOM' or tag == 'HETATM':

                        if eachLine.find('HOH') > -1:

                            x = eachLine[30:38]
                            y = eachLine[38:46]
                            z = eachLine[46:54]
                            x = float(x)
                            y = float(y)
                            z = float(z)

                            dist = (x_center - x) ** 2 + (y_center - y)**2 + (z_center - z)** 2

                            if dist < water_radius:
                                write_record = 'no'

                    if write_record == 'yes':
                        file.write(eachLine)

                file.close()

                # Recalculate map data without target waters 
                # Write run args
                tmpargs=[]
                tmpargs.append("refine")
                tmpargs.append('--pdbfile')
                tmpargs.append(pdbfile_omit)
                tmpargs.append('--mtzfile')
                tmpargs.append(mtzout_full)
                tmpargs.append('--workdir')
                tmpargs.append(workingdir)
                tmpargs.append('--libfile')
                tmpargs.append(libfile)
                tmpargs.append('--engine')
                tmpargs.append('refmac5')
                tmpargs.append('--weight')
                tmpargs.append('0.1')
                tmpargs.append('--max_res')
                tmpargs.append('none')
                tmpargs.append('--bref_type')
                tmpargs.append('none')
                tmpargs.append('--cycles')
                tmpargs.append('0')
                tmpargs.append('--water_cycles')
                tmpargs.append('0')
                tmpargs.append('--mifithome')
                tmpargs.append(mifit_root)
            
                import mi_refine
                if mi_refine.Run(tmpargs)!=0:
                    runcount = runcount + 1
                    continue # try next dataset
                
            os.environ['CCP4_SCR'] = ccp4_scratch

            ###############################################################
            # Automated ligand-fitting subscripts may be invoked here     #
            ###############################################################

            # Example script for this option uses FFFEAR to perform 6D search with a single i/p ligand conformation
            
            if fragcenter != 'none' and fragfile != 'none':

                print '\nRIGID-BODY LIGAND FITTING'
                print 'This option uses the example script: mi_ligandfit.py'

                # Obtain paths to the last coordinate file and last data file

                file = open('project_history.txt','r')
                allLines = file.readlines()
                file.close()

                for eachLine in allLines:

                    if eachLine.find('Output atoms:') > -1:
                        pdbfile = eachLine[13:200]
                        pdbfile = pdbfile.strip()

                    if eachLine.find('Output phased data:') > -1:
                        mtzout = eachLine[19:200]
                        mtzout = mtzout.strip()

                if pdbfile == 'none':
                    print '\nPDB file for rigid-body ligand fitting was not found\n'
                    time.sleep(4)
                    return 1

                if mtzout == 'none':
                    print '\nMTZ file for rigid-body ligand fitting was not found\n'
                    time.sleep(4)
                    return 1

                # Set output name for protein with ligand

                pdbfileout = pdbfile.replace('.pdb','_ligand.pdb')

                ############################################################################
                # Run the ligand-fitting routine                                           #
                # Ligand-fitting routines usually need:                                    #
                #  pdbfile = the path to current protein model                             #
                #  mtzout = the path to phased diffraction data from REFMAC                #
                #  workingdir = the path to the working directory                          #
                #  fragcenter = approx x,y,z coordinates for ligand                        #
                #  fragfile = the path to an input 3D model of ligand                      #
                #  pdbfileout = the path to o/p coordinates of protein with fitted ligand  #
                ############################################################################

                tmpargs=[]
                tmpargs.append("ligandfit")
                tmpargs.append('-f')
                tmpargs.append(fragfile)
                tmpargs.append('-m')
                tmpargs.append(mtzout)  
                tmpargs.append('-p')
                tmpargs.append(pdbfile)      
                tmpargs.append('-c')
                tmpargs.append(fragcenter)
                tmpargs.append('-d')
                tmpargs.append(workingdir)        
                tmpargs.append('-o')
                tmpargs.append(pdbfileout)
                import mi_ligandfit
                mi_ligandfit.Run(tmpargs)

            # User may supply script/technology here, with inputs patterned on the rigid-body exemple        

            if fragcenter != 'none' and chemfile != 'none':

                print '\nFLEXIBLE LIGAND FITTING'
                print 'This option requires a user-supplied script to be called from mi_bng.py'

            #########################################################
            # Setup crystal and data files for interactive graphics #
            #########################################################

            pdbfile = 'none'
            mtzout = 'none'

            print '\nFILE CREATION FOR MIFIT'
            print 'View center:',fragcenter

            # Obtain paths to the last coordinate file and last data file

            file = open('project_history.txt','r')
            allLines = file.readlines()
            file.close()

            for eachLine in allLines:

                if eachLine.find('Output atoms:') > -1:
                    pdbfile = eachLine[13:200]
                    pdbfile = pdbfile.strip()

                if eachLine.find('Output phased data:') > -1:
                    mtzout = eachLine[19:200]
                    mtzout = mtzout.strip()

            if pdbfile == 'none':
                print '\nPDB file for session file was not found\n'
                time.sleep(4)
                return 1

            if mtzout == 'none':
                print '\nMTZ file for session file was not found\n'
                time.sleep(4)
                return 1

            # Replace pdb file to version with ligand if automated-ligand fitting was done

            if fragcenter != 'none' and fragfile != 'none' and pdbfileout != 'none':
                fileexists = os.path.exists(pdbfileout)
                if fileexists != 0:
                    pdbfile = pdbfileout
                
            # Need local paths

            pdbfile_local = os.path.basename(pdbfile)
            mtzout_local = os.path.basename(mtzout)

            #   
            # Write a minimal session (mlw) file to launch the display
            #

            # Obtain a model center for translation keyword if no fragment center was given

            if fragcenter == 'none':

                xmean = 0.0
                ymean = 0.0
                zmean = 0.0
                number_atoms = 0.0

                file = open(pdbfile,'r')
                allLines = file.readlines()
                file.close()

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
                        xmean = x + xmean
                        ymean = y + ymean
                        zmean = z + zmean
                        number_atoms = number_atoms + 1.0

                xmean = xmean/number_atoms
                ymean = ymean/number_atoms
                zmean = zmean/number_atoms
                xmean = round(xmean,3)
                ymean = round(ymean,3)
                zmean = round(zmean,3)
                xmean = str(xmean)
                ymean = str(ymean)
                zmean = str(zmean)

                fragcenter_use = ' ' + xmean + ' ' + ymean + ' ' + zmean

            else:

                fragcenter_use = fragcenter

            # Write session file

            print 'Session file:',mlw_file

            file = open(mlw_file,'w')
            file.write('LoadPDB 1 ')
            file.write(pdbfile_local)
            file.write('\n')

            # Standard likelihood weighted map

            file.write('MapColumns FO=FWT PHI=PHWT\n')
            file.write('LoadMapPhase 1 ')
            file.write(mtzout)
            file.write('\n')
            file.write('silentmode\n')
            file.write('coefficients Direct FFT\n')
            file.write('fftapply\n')
            file.write('maptocont 1\n')
            file.write('maplinewidth 1.000000\n')
            file.write('contourradius ')
            file.write(contourradius)
            file.write('\n')

            # Write standard or user-defined colors

            if first_map == 'no':

                file.write('contourlevels 4\n')
                file.write('contourleveldefault 50.000000 100.000000 50.000000 200.000000 250.000000\n')
                file.write('color 21\n')
                file.write('contourcolor 1\n')
                file.write('color 22\n')
                file.write('contourcolor 2\n')
                file.write('color 23\n')
                file.write('contourcolor 3\n')
                file.write('color 24\n')
                file.write('contourcolor 4\n')
                file.write('color 25\n')
                file.write('contourcolor 5\n')

            else:

                contourlevels1 = aList_contourlevels1[0]
                contourleveldefault1 = aList_contourleveldefault1[0]

                aLine = contourleveldefault1.split()
                num_contours = len(aLine)

                file.write('contourlevels ')
                file.write(contourlevels1)
                file.write('\n')
                file.write('contourleveldefault ')
                file.write(contourleveldefault1)
                file.write('\n')

                count = 0
                while count < num_contours:

                    file.write('color ')
                    file.write(aList_color1[count])
                    file.write('\n')
                    file.write('contourcolor ')
                    file.write(aList_contourcolor1[count])
                    file.write('\n')

                    count = count + 1

            file.write('contourmap 1\n')

            # Likelihood weighted difference map

            file.write('MapColumns FO=DELFWT PHI=PHDELFWT\n')
            file.write('LoadMapPhase 2 ')
            file.write(mtzout)
            file.write('\n')
            file.write('silentmode\n')
            file.write('coefficients Direct FFT\n')
            file.write('fftapply\n')
            file.write('maptocont 2\n')
            file.write('maplinewidth 1.000000\n')
            file.write('contourradius ')
            file.write(contourradius)
            file.write('\n')

            # Write standard or user-defined colors

            if second_map == 'no':

                file.write('contourlevels 6\n')
                file.write('contourleveldefault -200.000000 -150.000000 150.000000 200.000000 250.000000\n')
                file.write('color 24\n')
                file.write('contourcolor 1\n')
                file.write('color 25\n')
                file.write('contourcolor 2\n')
                file.write('color 21\n')
                file.write('contourcolor 3\n')
                file.write('color 22\n')
                file.write('contourcolor 4\n')
                file.write('color 23\n')
                file.write('contourcolor 5\n')

            else:

                contourlevels2 = aList_contourlevels2[0]
                contourleveldefault2 = aList_contourleveldefault2[0]

                aLine = contourleveldefault2.split()
                num_contours = len(aLine)

                file.write('contourlevels ')
                file.write(contourlevels2)
                file.write('\n')
                file.write('contourleveldefault ')
                file.write(contourleveldefault2)
                file.write('\n')

                count = 0
                while count < num_contours:

                    file.write('color ')
                    file.write(aList_color2[count])
                    file.write('\n')
                    file.write('contourcolor ')
                    file.write(aList_contourcolor2[count])
                    file.write('\n')

                    count = count + 1

            file.write('contourmap 2\n')

            # View parameters

            file.write('translation ')    
            file.write(fragcenter_use)
            file.write('\nrotation ')
            file.write(fragview)
            file.write('\n')
            file.write('zoom ')
            file.write(zoom)
            file.write('\n')
            file.write('perspective 0.000\n')
            file.write('frontclip ')
            file.write(frontclip)
            file.write('\n')
            file.write('backclip ')
            file.write(backclip)
            file.write('\n')
            file.write('transform\n')
            file.write('stereo off\n')
            file.close()

            #
            # Option to write map around target point
            #

            if write_map == 'yes':

                mapout_local = mtzout_local.replace('.mtz','_site.map')

                # Compute complete map with CCP4/FFT

                file = open('mi_fft.inp','w')
                file.write('LABIN F1=FWT PHI=PHWT\n')
                file.write('END\n')
                file.close()

                runfft = 'fft HKLIN ' + mtzout_local + ' MAPOUT mi_2ff.map < mi_fft.inp > mi_fft.log'
                os.system(runfft)

                if not os.path.exists('mi_2ff.map'):
                    print 'FFT for map display failed'
                    time.sleep(4)
                    return 1
                else:
                    os.remove('mi_fft.inp')
                    os.remove('mi_fft.log')

                if fragcenter == 'none':
                    xyz_limits = 'BORDER ' + str(int(border))
                else:

                    # Obtain box coordinates

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

                    xyz_limits = 'XYZLIM ' + xmin + ' ' + xmax + ' ' + ymin + ' ' + ymax + ' ' + zmin + ' ' + zmax

                # Use CCP4/MAPMASK to box out required region

                file = open('mi_mapmask.inp','w')
                file.write(xyz_limits)
                file.write('\n')
                file.write('EXTEND XTAL\n')
                file.write('END\n')
                file.close()

                runmapmask = 'mapmask MAPIN mi_2ff.map XYZIN ' + pdbfile_local + ' MAPOUT ' + mapout_local + ' < mi_mapmask.inp > mi_mapmask.log'
                os.system(runmapmask)

                if not os.path.exists(mapout_local):
                    print 'CCP4/MAPMASK for target site failed'
                    time.sleep(4)
                    return 1
                else:
                    os.remove('mi_mapmask.inp')
                    os.remove('mi_mapmask.log')
                    os.remove('mi_2ff.map')

                    print 'Map of target region',mapout_local

            ###################################
            # Append content to HTML summary  #
            ###################################

            if bngsummary != 'none':

                # Learn information from project history

                read_ref_summary = 'no'
                refine_resolution = '?'
                refine_rwork = '?'
                refine_rfree = '?'
                error_list_path= '?'

                project_history_path = os.path.join(workingdir,'project_history.txt')

                file = open(project_history_path,'r')
                allLines = file.readlines()
                file.close()

                for eachLine in allLines:

                    # Parse for refinement information

                    if eachLine.find('Job ID:') > -1 and eachLine.find('refine_') > -1:
                        read_ref_summary = 'yes'

                    if eachLine.find('Summary:') > -1 and read_ref_summary == 'yes':
                        refine_summary = eachLine[16:100]

                        aLine = refine_summary.split()
                        number_args = len(aLine)

                        if number_args > 3:
                            refine_rwork = aLine[0]
                            refine_rwork = refine_rwork.replace('Rwork=','')
                            refine_rfree = aLine[1]
                            refine_rfree = refine_rfree.replace('Rfree=','')
                            refine_resolution = aLine[3]
                            refine_resolution = refine_resolution.replace('Resolution=','')

                        read_ref_summary = 'no'

                    if eachLine.find('Output error list:') > -1 and read_ref_summary == 'yes':
                        error_list_path = eachLine[18:200]
                        error_list_path = error_list_path.strip()

                # Setup summary data and links

                run_id_link = runcount + 1
                run_id_link = str(run_id_link)

                filename = bngsummaryfile
                file = open(filename,'a')

                # Link to history file

                file.write('<td><a href = "file:///')
                file.write(project_history_path)
                file.write('">')
                file.write(run_id_link)
                file.write('</a></td>\n')

                # Resolution

                file.write('<td>')
                file.write(refine_resolution)
                file.write('</td>\n')

                # Rw

                file.write('<td>')
                file.write(refine_rwork)
                file.write('</td>\n')

                # Rf

                file.write('<td>')
                file.write(refine_rfree)
                file.write('</td>\n')

                # Error list link

                file.write('<td><a href = "file:///')
                file.write(error_list_path)
                file.write('">')
                file.write(run_id_link)
                file.write('</a></td>\n')

                # Working directory name and link

                file.write('<td><a href = "file:///')
                file.write(workingdir)
                file.write('">')
                file.write(workingdir)
                file.write('</a></td></tr>\n')

                file.close()

            # End of structure solution 

        # End of loop

        runcount = runcount + 1

    # Close out HTML summary

    if bngsummary != 'none':

        runtime = time.ctime(time.time())

        filename = bngsummaryfile
        file = open(filename,'a')

        file.write('</table>')
        file.write('<p>\n')
        file.write('<b>Job End: ')
        file.write(runtime)
        file.write('</b><p>\n')
        file.write('</body>\n')
        file.write('</html>\n')
        file.close()

    #######################################
    # Launch MIFIT option for single runs #
    #######################################

    if number_datasets == 1 and launch_mifit == 'yes':

        os.chdir(workingdir)

        print '\nMIFIT LAUNCH'
        os.execl(mifitinstall,mlw_file,mlw_file)

    #
    return 0

if __name__ == "__main__":
    sys.exit(Run())


