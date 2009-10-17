#################################################################
#                                                               #
# Integrate images with MOSFLM and merge with SCALA             #
#                                                               #
# Copyright: Molecular Images   2007                            #
#                                                               #
# This script is distributed under the same conditions as MIFit #
#                                                               #
#################################################################

import sys
import os
import getopt
import time
import string
import dircache
import ccp4check

def Usage():
    print "Usage: %s [options]" % sys.argv[0]
    print "Options are:"
    print "  -t,--template_image=FILE:         name of template image file"
    print "  -s,--spacegroup=NUM               the space group number"
    print "  -f,--first_image=NUM              first image number to process, has default"
    print "  -l,--last_image=NUM               last image number to process, has default"
    print "  -i,--integrate_resolution=STRING  integrate resolution, if any"
    print "  -m,--merge_resolution=STRING      merging resolution, if any"
    print "  -g,--batch_prefix=NUM             group number, prefix for batch. default: 1"
    print "  -o,--scale_only=no or yes         only scale the images. default: no"
    print "  -w,--workdir=DIR                  the working directory. default: image directory"
    print "  -d,--detector_constants=FILE      file for detector constants: default: none"
    print "  -?,--help                         this help message"

def Run(argv=None):
    if argv is None:
        argv=sys.argv
    # Initialize

    first_image = 'none'
    last_image = 'none'
    dt_spacegroup = 'none'
    image_name = 'none'
    final_workdir = 'none'
    integrate_res = 'none'
    merging_res = 'none'
    scale_only = 'no'
    integrate_res = 'none'
    merging_res = 'none'
    batch_prefix = 'none'
    detector_constants = 'none'
    beam_x_image = 'none'
    beam_y_image = 'none'
    ioversigi_limit = 1.0
    image_extension = 4
    gain = '1.0'

    second_index = 45
    second_index_prev = second_index - 1

    quote = '''"'''

    ##################
    # Parse args     #
    ##################
    args = argv[1:]
    optlist, args = getopt.getopt(
        args,'t:s:w:f:l:i:m:g:o:d:?',
        ['template_image=','spacegroup=','first_image=',
         'last_image=','integrate_resolution=','merge_resolution=',
         'batch_prefix=','scale_only=','workdir=','detector_constants=',
         'help'])
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
        if arg_value=='-?' or arg_value=='--help':
            Usage()
            return
        if number_of_list_inputs >=2:
            param_value = aList[1]
            if arg_value == '-t' or arg_value=='--template_image':
                image_name = param_value
            elif arg_value == '-s' or arg_value=='--spacegroup':
                dt_spacegroup = param_value
            elif arg_value == '-f' or arg_value=='--first_image':
                first_image = param_value
            elif arg_value == '-l' or arg_value=='--last_image':
                last_image = param_value
            elif arg_value == '-i' or arg_value=='--integrate_resolution':
                integrate_res = param_value
            elif arg_value == '-m' or arg_value=='--merge_resolution':
                merging_res = param_value
            elif arg_value == '-g' or arg_value=='--batch_prefix':
                batch_prefix = param_value
            elif arg_value == '-o' or arg_value=='--scale_only':
                scale_only = param_value
            elif arg_value == '-w' or arg_value=='--workdir':
                final_workdir = param_value
            elif arg_value == '-d' or arg_value=='--detector_constants':
                detector_constants = param_value                
        count = count + 1        

    #######################
    # Checks and defaults #
    #######################        

    ccp4,error = ccp4check.ccp4check()
    if not ccp4:
        print '\n' + error + '\n'
        time.sleep(4)
        return 1

    # Check MOSFLM is in place in standard install directory in CCP4.
    # Path may be changed for independent installations.

    test_platform = sys.platform
    if test_platform.find('win') > -1:
        ipmosflm_path = os.path.join(ccp4.bin,'ipmosflm.exe')        
    else:   
        ipmosflm_path = os.path.join(ccp4.bin,'ipmosflm')

    fileexists = os.path.exists(ipmosflm_path)
    if fileexists == 0:
        print '\nMOSFLM executable was not found as',ipmosflm_path
        time.sleep(4)
        return 1 

    # Check for image directory

    fileexists = os.path.exists(image_name)
    if fileexists == 0:
        print 'The image directory/template data file was not found:',image_name
        time.sleep(4)
        return 1

    # Check for processed file directory

    filexists = os.path.exists(final_workdir)
    if fileexists == 0 and final_workdir != 'none':
        print 'The final directory for processed data was not found:',final_workdir
        time.sleep(4)
        return 1

    # Check for space group

    if dt_spacegroup == 'none':
        print 'The space group number was not given'
        time.sleep(4)
        return 1

    # Merging resolution limits

    if merging_res != 'none':
        aLine = merging_res.split()
        num = len(aLine)
        if num != 2:
           print 'Merging resolution should be two numbers'
           time.sleep(4)
           return 1

    # Set Data run identification

    if batch_prefix == 'none':
        data_code_number = '1'
    else:
        data_code_number = batch_prefix

    # Set path to images depending on whether template is directory or a file

    if os.path.isfile(image_name) == 0:
        image_dir = image_name
        image_name = 'none'
    else:
        image_dir = os.path.dirname(image_name)

    # Find image root from files in folder

    aList_image_number = []
    aList_image_number_low = []
    aList_image_number_high = []

    aList_image_root = []
    aList_image_root_count = []
    aList_image_file = []
    aList_image_extension = []   

    image_number_low = 99999
    image_number_high = -99999
    image_span = 0
    image_span_prev = 0

    aList_dir = dircache.listdir(image_dir)
    number_files = len(aList_dir)

    # Ensure image ranges are set

    if first_image != 'none':
        image_number_low = int(first_image)
    if last_image != 'none':
        image_number_high = int(last_image)

    #

    count = 0
    while count < number_files:
        imagefile = aList_dir[count]

        # Detect image files using extension .img or .osc preceeded by a 3 or 4 digit number

        if imagefile.find('.img') > -1 or imagefile.find('.osc') > -1:
            imagefile_split = imagefile.split('.')
            imagefile_root = imagefile_split[0]
            i_end = len(imagefile_root)
            i_start = i_end - 4
            image_number = imagefile_root[i_start:i_end]

            if image_number.isdigit() == 1:               
                image_number = int(image_number)
                image_extension = 4
            else:
                i_start = i_end - 3
                image_number = imagefile_root[i_start:i_end]
                image_extension = 3

                if image_number.isdigit() == 1:
                    image_number = int(image_number)
                else:
                    print '\nImage file names must contain 3 or 4 digits preceeding extension .osc/.img\n'
                    time.sleep(4)
                    return 1

            # Collect the image number, file name and file root for each image

            aList_image_number.append(image_number)
            aList_image_file.append(imagefile)
            aList_image_extension.append(image_extension)
                
            true_imagefile_root = imagefile_root[0:i_start]               
            aList_image_root.append(true_imagefile_root)

        count = count + 1

    number_of_images = len(aList_image_number)

    if number_of_images < 2:
        print 'Number of images found was only ',number_of_images
        time.sleep(4)
        return 1

    # Get the correct file root - the one with most instances

    count = 0
    while count < number_of_images:

        test_root = aList_image_root[count]
        root_count = aList_image_root.count(test_root)
        aList_image_root_count.append(root_count)

        count = count + 1

    test_count_prev = 0
    test_count_max = 0
      
    count = 0
    while count < number_of_images:

        test_count = aList_image_root_count[count]

        if test_count > test_count_prev:
            root_index = count
            test_count_prev = test_count

        count = count + 1

    true_root_image = aList_image_root[root_index]
    image_extension = aList_image_extension[root_index]

    # Automated determination of image range

    if first_image == 'none' or last_image == 'none':

        # Find the largest contiguous set of frames with the correct file root

        image_number_target = 0
        count = 0
        while count < number_of_images:

            test_image_root = aList_image_root[count]

            if test_image_root == true_root_image and image_number_target == 0:
                image_number_target = aList_image_number[count]

            count = count + 1

        count = 0
        while count < number_of_images:

            test_image_root = aList_image_root[count]

            if test_image_root == true_root_image:

                image_number = aList_image_number[count]

                # Capture and reset if the image numbers break

                if image_number != image_number_target:
                    aList_image_number_low.append(image_number_low)
                    aList_image_number_high.append(image_number_high)                
                    image_number_low = image_number
                    image_number_high = image_number

                if image_number < image_number_low:
                    image_number_low = image_number

                if image_number > image_number_high:
                    image_number_high = image_number

                image_number_target = image_number + 1

            count = count + 1

        aList_image_number_low.append(image_number_low)
        aList_image_number_high.append(image_number_high)

        # Now get largest segment from lists of low/high values

        number_of_spans = len(aList_image_number_low)

        count = 0
        while count < number_of_spans:
            image_number_low_temp = aList_image_number_low[count]
            image_number_high_temp = aList_image_number_high[count]

            image_span = image_number_high_temp - image_number_low_temp

            if image_span > image_span_prev:
                image_number_low = image_number_low_temp
                image_number_high = image_number_high_temp
                image_span_prev = image_span

            count = count + 1

        if first_image == 'none':
            first_image = str(image_number_low)

        if last_image == 'none':
            last_image = str(image_number_high)

    #############################
    # Setup program parameters  #
    #############################

    # Total number of images

    number_images = image_number_high - image_number_low + 1

    # Establish the indexing and refinement images

    index_first = image_number_low
    refine_segment1_first = image_number_low
    refine_segment1_last = image_number_low + 3
    cell_refine_images_segment1 = str(refine_segment1_first) + ' to ' + str(refine_segment1_last)

    if number_images > second_index:
        index_second = image_number_low + second_index_prev
        refine_segment2_last = image_number_low + second_index_prev
 
    else:
        index_second = image_number_high
        refine_segment2_last = image_number_high

    refine_segment2_first = refine_segment2_last - 3
    cell_refine_images_segment2 = str(refine_segment2_first) + ' to ' + str(refine_segment2_last)   

    image_seq_find = str(index_first) + ' ' + str(index_second)

    # Set template image name to first image if not defined

    if image_name == 'none':
    
        count = 0
        while count < number_of_images:

            test_image_number = aList_image_number[count]
            test_image_number = str(test_image_number)
            
            test_image_root = aList_image_root[count]

            if test_image_root == true_root_image and test_image_number == first_image:
                imagefile = aList_image_file[count]
                image_name = os.path.join(image_dir,imagefile)

            count = count + 1    

        if image_name == 'none':
            print 'Automated identification of the template image failed'
            time.sleep(4)
            return 1

    # Establish MOSFLM style templates

    num = len(image_name)

    if image_extension == 3:
        num = num - 7
        hashes = '###'
    else:
        num = num - 8
        hashes = '####'

    if image_name.find('.img') > -1:
        mosflm_template_full = image_name[0:num] + hashes + '.img'

    if image_name.find('.osc') > -1:
        mosflm_template_full = image_name[0:num] + hashes + '.osc'     

    mosflm_template = os.path.basename(mosflm_template_full)

    # Go to frame directory

    os.chdir(image_dir)

    # Collect user-defined beamcenter or distance data from a special constants file

    beam_x = 'none'
    beam_y = 'none'
    xtal_detector_distance = 'none'

    fileexists = os.path.exists(detector_constants)
    if fileexists != 0 and detector_constants != 'none':

        file = open(detector_constants,'r')
        allLines = file.readlines()
        file.close()

        for eachLine in allLines:

            if eachLine.find('beam_center') > -1:
                aLine = eachLine.split()
                number_args = len(aLine)
                if number_args > 2:
                    beam_x = aLine[1]
                    beam_y = aLine[2] 
                else:
                    print 'There should be two numbers on the beam_center line'
                    time.sleep(4)
                    return 1

            if eachLine.find('xtal_detector_distance') > -1:
                aLine = eachLine.split()
                number_args = len(aLine)

                if number_args > 1:
                    xtal_detector_distance = aLine[1]
                else:
                    print 'There should be one number on the xtal_detector_distance line'
                    time.sleep(4)
                    return 1

            if eachLine.find('detector_gain') > -1:
                aLine = eachLine.split()
                number_args = len(aLine)

                if number_args > 1:
                    detector_gain = aLine[1]
                else:
                    print 'There should be one number on the detector_gain line'
                    time.sleep(4)
                    return 1

    ##########
    # Start  #
    ##########

    print '\nAutomated integration and merging\n'
    print 'Image directory:',image_dir
    print 'Template image:',image_name
    print 'Batch number:',batch_prefix
    print 'First image to process:',first_image
    print 'Last image to process:',last_image
    print 'Images for initial indexing:',image_seq_find
    print 'Integration resolution limits:',integrate_res
    print 'Merging resolution limits:',merging_res
    print 'Expected space group number:',dt_spacegroup

    time.sleep(5)

    if beam_x != 'none':
        print 'Using input beam center:',beam_x,beam_y

    if xtal_detector_distance != 'none':
        print 'Using input detector distance:',xtal_detector_distance

    # Process log

    runtime = time.ctime(time.time())

    file = open('autoprocess.log','w')
    file.write('Processing start time       : ')
    file.write(runtime)
    file.write('\n')
    file.close()

    if scale_only != 'yes':

        image_name_base = os.path.basename(image_name)

        fileexists = os.path.exists('SUMMARY')
        if fileexists != 0:
            os.remove('SUMMARY')

        fileexists = os.path.exists('NEWMAT')
        if fileexists != 0:
            os.remove('NEWMAT')

        fileexists = os.path.exists('NEWMAT_REFINED')
        if fileexists != 0:
            os.remove('NEWMAT_REFINED')

        ###############
        # Auto index  #
        ###############

        # may need to search over threshold and frame nos

        runtime = time.ctime(time.time())
        print 'Date:',runtime
        print 'Autoindexing'

        file = open('mi_mosflm_index.inp','w')

        file.write('TITLE Autoindexing\n')
        file.write('TEMPLATE ')
        file.write(mosflm_template)
        file.write('\n') 
        file.write('DIRECTORY "')
        file.write(image_dir)
        file.write('"\n')

        # Apply user specified beam center if given

        if beam_x != 'none':
            file.write('BEAM ')
            file.write(beam_x)
            file.write(' ')
            file.write(beam_y)
            file.write('\n')

        # Apply user specified distance if given

        if xtal_detector_distance != 'none':
            file.write('DISTANCE ')
            file.write(xtal_detector_distance)
            file.write('\n')

        # Symmetry and indexing

        file.write('SYMM ')
        file.write(dt_spacegroup)
        file.write('\n')

        file.write('SEPARATION CLOSE\n')
        file.write('FINDSPOTS THRESHOLD 10\n')

        file.write('AUTOINDEX DPS IMAGES ')
        file.write(image_seq_find)
        file.write('\n')

        file.write('GO\n')
        file.close()

        # Execute

        runmosflm = ipmosflm_path + ' < mi_mosflm_index.inp > mi_mosflm_index.log'
        os.system(runmosflm)

        fileexists = os.path.exists('NEWMAT')
        if fileexists != 0:

            os.remove('mi_mosflm_index.inp')

            runtime = time.ctime(time.time())

            file = open('autoprocess.log','a')
            file.write('Autoindexing done           : ')
            file.write(runtime)
            file.write('\n')
            file.close()

            # Obtain beam center

            file = open('mi_mosflm_index.log')   
            allLines = file.readlines()
            file.close()

            for eachLine in allLines:
                if eachLine.find('Beam coordinates of') > -1 and eachLine.find('have been refined') > -1:
                    aLine = eachLine.split()
                    beam_x_image = aLine[9]
                    beam_y_image = aLine[10]

            if beam_x_image == 'none' or beam_y_image == 'none':
                print 'Parsing of beam center seems to have failed'
                time.sleep(4)
                return 1

        else:

            print 'MOSFLM Autoindexing seems to have failed'
            time.sleep(4)
            return 1

        fileexists = os.path.exists('COORDS')
        if fileexists != 0:
            os.remove('COORDS')

        fileexists = os.path.exists('SUMMARY')
        if fileexists != 0:
            os.remove('SUMMARY')

        #####################
        # Cell Refinement   #
        #####################

        runtime = time.ctime(time.time())
        print 'Date:',runtime
        print 'Cell refinement' 

        file = open('mi_mosflm_refine.inp','w')

        file.write('TITLE Cell refinement\n')
        file.write('TEMPLATE ')
        file.write(mosflm_template)
        file.write('\n') 
        file.write('DIRECTORY "')
        file.write(image_dir)
        file.write('"\n')

        file.write('MATRIX NEWMAT\n')

        # Apply user specified beam center if given

        if beam_x != 'none':
            beam_x_apply = beam_x
            beam_y_apply = beam_y
        else:
            beam_x_apply = beam_x_image
            beam_y_apply = beam_y_image

        file.write('BEAM ')
        file.write(beam_x_apply)
        file.write(' ')
        file.write(beam_y_apply)
        file.write('\n')
        file.write('BACKSTOP CENTRE ')
        file.write(beam_x_apply)
        file.write(' ')
        file.write(beam_y_apply)
        file.write(' RADIUS 5.00\n')

        # Apply user specified distance if given

        if xtal_detector_distance != 'none':
            file.write('DISTANCE ')
            file.write(xtal_detector_distance)
            file.write('\n')

        file.write('SYMM ')
        file.write(dt_spacegroup)
        file.write('\n')

        # Machine and crystal default

        file.write('MOSAIC 0.70\n')
        file.write('SEPARATION CLOSE\n')
        file.write('GAIN ')
        file.write(gain)
        file.write('\n')
        file.write('OVERLOAD CUTOFF 65500\n')
        file.write('DISTORTION YSCALE 1.0000 TILT 0 TWIST 0\n')

        # Refinement

        file.write('NEWMATRIX NEWMAT_REFINED\n')
        file.write('POSTREF SEGMENT 2 MAXRESIDUAL 1.3 SHIFTFAC 3.0 MAXSHIFT 0.1 RESOLUTION 4.0\n')
        file.write('PROCESS ')
        file.write(cell_refine_images_segment1)
        file.write('\nGO\n')
        file.write('PROCESS ')
        file.write(cell_refine_images_segment2)
        file.write('\nGO\n')
        file.close()

        # Execute

        runmosflm = ipmosflm_path +  ' < mi_mosflm_refine.inp > mi_mosflm_refine.log'
        os.system(runmosflm)

        fileexists = os.path.exists('NEWMAT_REFINED')
        if fileexists != 0:

            os.remove('mi_mosflm_refine.inp')

            runtime = time.ctime(time.time())

            file = open('autoprocess.log','a')
            file.write('Cell refinement done        : ')
            file.write(runtime)
            file.write('\n')
            file.close()

        else:

            print 'MOSFLM cell refinement seems to have failed'
            time.sleep(4)
            return 1

        fileexists = os.path.exists('SUMMARY')
        if fileexists != 0:
            os.remove('SUMMARY')

        fileexists = os.path.exists('GENFILE')
        if fileexists != 0:
            os.remove('GENFILE')

        fileexists = os.path.exists('NEWMAT')
        if fileexists != 0:
            os.remove('NEWMAT')

        ###########################
        # Integrate               #
        ###########################

        fileexists = os.path.exists('mi_integrate.mtz')
        if fileexists != 0:
            os.remove('mi_integrate.mtz')

        runtime = time.ctime(time.time())
        print 'Date:',runtime
        print 'Image integration' 

        file = open('mi_mosflm_integrate.inp','w')

        file.write('TITLE Cell refinement\n')
        file.write('TEMPLATE ')
        file.write(mosflm_template)
        file.write('\n') 
        file.write('DIRECTORY "')
        file.write(image_dir)
        file.write('"\n')

        file.write('MATRIX NEWMAT_REFINED\n')
        file.write('GENFILE GENFILE\n')
        file.write('HKLOUT mi_integration.mtz\n')

        # Apply user specified beam center if given

        if beam_x != 'none':
            beam_x_apply = beam_x
            beam_y_apply = beam_y
        else:
            beam_x_apply = beam_x_image
            beam_y_apply = beam_y_image

        file.write('BEAM ')
        file.write(beam_x_apply)
        file.write(' ')
        file.write(beam_y_apply)
        file.write('\n')
        file.write('BACKSTOP CENTRE ')
        file.write(beam_x_apply)
        file.write(' ')
        file.write(beam_y_apply)
        file.write(' RADIUS 5.00\n')

        # Apply user specified distance if given

        if xtal_detector_distance != 'none':
            file.write('DISTANCE ')
            file.write(xtal_detector_distance)
            file.write('\n')

        file.write('SYMM ')
        file.write(dt_spacegroup)
        file.write('\n')

        # Machine and crystal default

        file.write('MOSAIC 0.70\n')
        file.write('GAIN ')
        file.write(gain)
        file.write('\n')
        file.write('OVERLOAD CUTOFF 65500\n')
        file.write('DISTORTION YSCALE 1.0000 TILT 0 TWIST 0\n')

        # Refinement

        if integrate_res != 'none':
            file.write('RESOLUTION ')
            file.write(integrate_res)
            file.write('\n') 

        file.write('POSTREF FIX ALL\n')
        file.write('PROCESS ')
        file.write(first_image)
        file.write(' TO ')
        file.write(last_image)
        file.write('\nGO\n')
        file.close()

        # Execute

        runmosflm = ipmosflm_path + ' < mi_mosflm_integrate.inp > mi_mosflm_integrate.log' 
        os.system(runmosflm)

        fileexists = os.path.exists('mi_integration.mtz')
        if fileexists != 0:
            os.remove('mi_mosflm_integrate.inp')
        else:
            print 'MOSFLM integration seems to have failed'
            time.sleep(4)
            return 1

        fileexists = os.path.exists('SUMMARY')
        if fileexists != 0:
            os.remove('SUMMARY')

        fileexists = os.path.exists('GENFILE.gen')
        if fileexists != 0:
            os.remove('GENFILE.gen')

        fileexists = os.path.exists('NEWMAT_REFINED')
        if fileexists != 0:
            os.remove('NEWMAT_REFINED')

        #################################
        # Need to sort prior to merging #
        #################################

        fileexists = os.path.exists('mi_integration_sorted.mtz')
        if fileexists != 0:
            os.remove('mi_integration_sorted.mtz')

        file = open('mi_sortmtz.inp','w')

        file.write('H K L M/ISYM BATCH\n')
        file.close()

        run_sortmtz = 'sortmtz HKLIN mi_integration.mtz HKLOUT mi_integration_sorted.mtz < mi_sortmtz.inp > mi_sortmtz.log'

        os.system(run_sortmtz)

        fileexists = os.path.exists('mi_integration_sorted.mtz')
        if fileexists == 0:
            print 'Sorting process failed'
            time.sleep(4)
            return 1
        else:
            os.remove('mi_integration.mtz')
            os.rename('mi_integration_sorted.mtz','mi_integration.mtz')
            os.remove('mi_sortmtz.inp')
            os.remove('mi_sortmtz.log')

    ###############################################################
    # Scale and average the integrated and profile-fitted reflns  #
    ###############################################################

    runtime = time.ctime(time.time())

    print 'Date:',runtime
    print '\nMerging'

    fileexists = os.path.exists('mi_integration.mtz')
    if fileexists == 0:
        print 'File mi_integration.mtz is not available for merging'
        time.sleep(4)
        return 1
    else:
        runtime = time.ctime(time.time())

        file = open('autoprocess.log','a')
        file.write('Integration done            : ')
        file.write(runtime)
        file.write('\n')
        file.close()

    output_ref = 'ScalAverage_' + data_code_number + '.mtz'
    output_log = 'scala_scaleaverage_' + data_code_number + '.log'
    output_ref_initial = 'ScalAverage_' + data_code_number + '_initial.mtz'
    output_log_initial = 'scala_scaleaverage_' + data_code_number + '_initial.log'

    file = open('mi_scala.inp','w')

    file.write('TITLE First pass merging\n')
    file.write('RUN 1 ALL\n')
    file.write('INTENSITIES PARTIAL\n')
    file.write('CYCLES 20\n')
    file.write('ANOMALOUS OFF\n')
    file.write('SDCORRECTION 1.3 0.02\n')  

    # Protocol allows for some decay

    file.write('SCALES ROTATION SPACING 5 SECONDARY 6 TAILS BFACTOR ON BROTATION SPACING 20\n')
    file.write('TIE BFACTOR 0.5\n')
    file.write('REJECT 6 6 ALL -8 -8\n')   
    file.write('EXCLUDE EMAX 10\n')

    if merging_res != 'none':
        file.write('RESOLUTION ')
        file.write(merging_res)
        file.write('\n')

    file.close()

    run_scala = 'scala HKLIN mi_integration.mtz HKLOUT ' + output_ref +\
                ' SCALES mi_scales.txt ROGUES mi_rogues.txt NORMPLOT mi_normplot.txt PLOT mi_plot.txt < mi_scala.inp > ' \
                 + output_log

    os.system(run_scala)

    fileexists = os.path.exists(output_ref)
    if fileexists == 0:    
        print 'Process SCALA seems to have failed'
        time.sleep(4)
        return 1
    else:

        os.remove('mi_scala.inp')

        runtime = time.ctime(time.time())

        file = open('autoprocess.log','a')
        file.write('Merging done                : ')
        file.write(runtime)
        file.write('\n')
        file.close()

    fileexists = os.path.exists('mi_scales.txt')
    if fileexists == 1:
        os.remove('mi_scales.txt')

    fileexists = os.path.exists('mi_rogues.txt')
    if fileexists == 1:
        os.remove('mi_rogues.txt')

    fileexists = os.path.exists('mi_normplot.txt')
    if fileexists == 1:
        os.remove('mi_normplot.txt')

    fileexists = os.path.exists('mi_plot.txt')
    if fileexists == 1:
        os.remove('mi_plot.txt')

    fileexists = os.path.exists('fort.10')
    if fileexists == 1:
        os.remove('fort.10')

    fileexists = os.path.exists('COORDS')
    if fileexists == 1:
        os.remove('COORDS')

    ####################################################################
    # (Re) Scale and average the integrated and profile-fitted reflns  #
    # This is approximate - may be better criteria available           # 
    ####################################################################

    # Capture table data to see how initial merging went (note: depend on precise SCALA file format)

    aList_res_high = []
    aList_ioversigi = []
    table_length = 0
    read_table = 'no'

    file = open(output_log,'r')
    allLines = file.readlines()
    file.close()

    for eachLine in allLines:

        if eachLine.find('N 1/d^2 Dmin(A) Rmrg  Rfull   Rcum  Ranom  Nanom') > -1:
            read_table = 'yes'

        if read_table == 'yes':

            aLine = eachLine.split()
            number_items = len(aLine)

            if number_items == 18:

                res_high = aLine[2]
                ioversigi = aLine[12]

                if res_high.find('.') > -1 and ioversigi.find('.') > -1:
                    aList_res_high.append(res_high)
                    aList_ioversigi.append(ioversigi)

        if eachLine.find('Overall:') > -1:
            read_table = 'no'

    # 

    table_length = len(aList_res_high)
    max_table_length_index = table_length - 1

    if table_length > 1:

        count = 0
        while count < table_length:

            res_high = aList_res_high[count]
            ioversigi = aList_ioversigi[count]
            res_high = float(res_high)
            ioversigi = float(ioversigi)

            if ioversigi > ioversigi_limit:
                res_to_process = res_high

            count = count + 1

        # Remerge at lower resolution if required by I/sigI criteria

        res_upper = aList_res_high[max_table_length_index]
        res_upper = float(res_upper)

        if res_to_process > res_upper:

            res_to_process = str(res_to_process)

            runtime = time.ctime(time.time())

            print 'Date:',runtime
            print '\nRemerging to',res_to_process,'angstroms resolution'

            fileexists = os.path.exists(output_ref_initial)
            if fileexists != 0:
                os.remove(output_ref_initial)

            fileexists = os.path.exists(output_log_initial)
            if fileexists != 0:
                os.remove(output_log_initial)       

            os.rename(output_ref,output_ref_initial)
            os.rename(output_log,output_log_initial)

            res_to_process = str(res_to_process)

            file = open('mi_scala.inp','w')

            file.write('TITLE First pass remerging\n')
            file.write('RUN 1 ALL\n')
            file.write('INTENSITIES PARTIAL\n')
            file.write('CYCLES 20\n')
            file.write('ANOMALOUS OFF\n')
            file.write('SDCORRECTION 1.3 0.02\n')  
            file.write('SCALES ROTATION SPACING 5 SECONDARY 6 TAILS BFACTOR ON BROTATION SPACING 20\n')
            file.write('TIE BFACTOR 0.5\n')
            file.write('REJECT 6 6 ALL -8 -8\n')
            file.write('EXCLUDE EMAX 10\n')
            file.write('RESOLUTION HIGH ')
            file.write(res_to_process)
            file.write('\n')
            file.close()

            run_scala = 'scala HKLIN mi_integration.mtz HKLOUT ' + output_ref +\
                        ' SCALES mi_scales.txt ROGUES mi_rogues.txt NORMPLOT mi_normplot.txt PLOT mi_plot.txt < mi_scala.inp > ' \
                        + output_log

            os.system(run_scala)

            fileexists = os.path.exists(output_ref)
            if fileexists == 0:    
                print 'Process SCALA seems to have failed'
                time.sleep(4)
                return 1
            else:

                os.remove('mi_scala.inp')

                runtime = time.ctime(time.time())

                file = open('autoprocess.log','a')
                file.write('Remerging done              : ')
                file.write(runtime)
                file.write('\n')
                file.close()

            fileexists = os.path.exists('mi_scales.txt')
            if fileexists == 1:
                os.remove('mi_scales.txt')

            fileexists = os.path.exists('mi_rogues.txt')
            if fileexists == 1:
                os.remove('mi_rogues.txt')

            fileexists = os.path.exists('mi_normplot.txt')
            if fileexists == 1:
                os.remove('mi_normplot.txt')

            fileexists = os.path.exists('mi_plot.txt')
            if fileexists == 1:
                os.remove('mi_plot.txt')

            fileexists = os.path.exists('fort.10')
            if fileexists == 1:
                os.remove('fort.10')

            fileexists = os.path.exists('COORDS')
            if fileexists == 1:
                os.remove('COORDS')

    else:

        print 'Warning - parsing to test for remerging failed'

    ######################
    # Tail end processes #
    ######################

    # Transfer the final reflection and output files to defined space or leave in image directory

    if final_workdir != 'none' and final_workdir != image_dir:

        output_ref_destination = os.path.join(final_workdir,output_ref)
        output_log_destination = os.path.join(final_workdir,output_log)

        fileexists = os.path.exists(output_ref_destination)
        if fileexists != 0:
            os.remove(output_ref_destination)

        fileexists = os.path.exists(output_log_destination)
        if fileexists != 0:
            os.remove(output_log_destination)

        os.rename(output_ref,output_ref_destination)
        os.rename(output_log,output_log_destination)

    else:

        output_ref_destination = os.path.join(image_dir,output_ref)
        output_log_destination = os.path.join(image_dir,output_log)

    # Log and clean-up

    runtime = time.ctime(time.time())

    print 'Date:',runtime
    print '\nIntegrated and merged intensity file:',output_ref_destination
    print 'Merging statistics file:',output_log_destination

    time.sleep(4)
    return 0

if __name__ == "__main__":
    sys.exit(Run())


