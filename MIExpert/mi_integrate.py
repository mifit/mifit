#################################################################
#                                                               #
# Integrate images and merge with D*TREK                        #
#                                                               #
# Loosely based on Jim Pflugrath's shell scripts for automated  #
# data processing                                               #
#                                                               #
# Copyright: Molecular Images   2005                            #
#                                                               #
#################################################################

import sys
import os
import getopt
import time
import string
import dircache

def Usage():
    print "Usage: %s [options]" % sys.argv[0]
    print "Options are:"
    print "  -t,--template_image=DIR           name of image folder or template file"
    print "  -s,--spacegroup=NUM               the space group number"
    print "  -b,--beammask=FILE                none (default), rigakuccd, or filename"
    print "  -f,--first_image=NUM              first image number to process, has default"
    print "  -l,--last_image=NUM               last image number to process, has default"
    print "  -i,--integrate_resolution=STRING  integrate resolution, if any"
    print "  -m,--merge_resolution=STRING      merging resolution, if any"
    print "  -h,--start_header=FILE            header image file. Default: template image"
    print "  -g,--batch_prefix=NUM             group number, prefix for batch. default: 1"
    print "  -o,--scale_only=no or yes         only scale the images. default: no"
    print "  -w,--workdir=DIR                  the working directory. Default, image directory"
    print "  -d,--detector_constants=FILE      file for detector constants: default: none"
    print "  -?,--help                         this help message"

def Run(argv=None):
    if argv is None:
        argv=sys.argv
    # Initialize

    first_image = 'none'
    second_image = 'none'
    last_image = 'none'
    dt_spacegroup = 'none'
    image_name = 'none'
    start_head = 'none'
    final_workdir = 'none'
    integrate_res = 'none'
    merging_res = 'none'
    scale_only = 'no'
    beammask = 'none'
    integrate_res = 'none'
    merging_res = 'none'
    batch_prefix = 'none'
    detector_constants = 'none'

    ioversigi_limit = 1.5
    image_extension = 4

    ##################
    # parse args     #
    ##################
    args = argv[1:]
    optlist, args = getopt.getopt(
        args,'t:s:b:f:l:i:m:h:g:o:w:d:?',
        ['template_image=','spacegroup=','beammask=','first_image=',
         'last_image=','integrate_resolution=','merge_resolution=',
         'start_header=','batch_prefix=','scale_only=','workdir=','detector_constants=',
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
        if arg_value == '-h' or arg_value=='--help':
            Usage()
            return
        if number_of_list_inputs >=2:
            param_value = aList[1]
            if arg_value == '-t' or arg_value=='--template_image':
                image_name = param_value
            elif arg_value == '-s' or arg_value=='--spacegroup':
                dt_spacegroup = param_value
            elif arg_value == '-b' or arg_value=='--beammask':
                beammask = param_value
            elif arg_value == '-f' or arg_value=='--first_image':
                first_image = param_value
            elif arg_value == '-l' or arg_value=='--last_image':
                last_image = param_value
            elif arg_value == '-i' or arg_value=='--integrate_resolution':
                integrate_res = param_value
            elif arg_value == '-m' or arg_value=='--merge_resolution':
                merging_res = param_value
            elif arg_value == '-h' or arg_value=='--start_header':
                start_head = param_value
            elif arg_value == '-g' or arg_value=='--batch_prefix':
                batch_prefix = param_value
            elif arg_value == '-o' or arg_value=='--scale_only':
                scale_only = param_value
            elif arg_value == '-w' or arg_value=='--workdir':
                final_workdir = param_value
            elif arg_value == '-d' or arg_value=='--detector_constants':
                detector_constants = param_value                
        count = count + 1        

    # Tip for faster execution 

    os.environ['DTREK_REFLN_BINARY'] = 'YES'

    # Check the D*TREK environment
    env_vars = os.environ.keys()
    find_dtrek = env_vars.count('DTREK_ROOT')
    if find_dtrek == 0:
        print '\nDTREK_ROOT was not found - is the d*trek environment established ?'
        time.sleep(4)
        return 1

    # Set DSTARTREK program paths including OS dependence

    dtrek_path = os.environ['DTREK_ROOT']

    test_platform = sys.platform
    if test_platform.find('win') > -1:
        dtrek_path_exec = os.path.join(dtrek_path,'Bin_Win')
    else:
        find_dtrek_bin_suffix = env_vars.count('DTREK_BINSUFFIX')
        if find_dtrek_bin_suffix == 0:
            print '\nDTREK_BINSUFFIX was not found'
            time.sleep(4)
            return 1
        else:
            dtrek_bin_suffix = os.environ['DTREK_BINSUFFIX']
            dtrek_bin = 'bin' + dtrek_bin_suffix
            dtrek_path_exec = os.path.join(dtrek_path,dtrek_bin)

    dtheadermerge_path = os.path.join(dtrek_path_exec,'dtheadermerge')
    dtfind_path = os.path.join(dtrek_path_exec,'dtfind')
    dtindex_path = os.path.join(dtrek_path_exec,'dtindex')
    dtrefine_path = os.path.join(dtrek_path_exec,'dtrefine')
    dtpredict_path = os.path.join(dtrek_path_exec,'dtpredict')
    dtintegrate_path = os.path.join(dtrek_path_exec,'dtintegrate')
    dtscaleaverage_path = os.path.join(dtrek_path_exec,'dtscaleaverage')
    dtheaderedit_path = os.path.join(dtrek_path_exec,'dtheaderedit')

    quote = '''"'''
    dtheadermerge = quote + dtheadermerge_path + quote
    dtfind = quote + dtfind_path + quote
    dtindex = quote + dtindex_path + quote
    dtrefine = quote + dtrefine_path + quote
    dtpredict = quote + dtpredict_path + quote
    dtintegrate = quote + dtintegrate_path + quote
    dtscaleaverage = quote + dtscaleaverage_path + quote
    dtheaderedit = quote + dtheaderedit_path + quote

    #######################
    # Checks and defaults #
    #######################        

    # Check for image directory

    fileexists = os.path.exists(image_name)
    if fileexists == 0:
        print 'The image directory/template data file was not found:',image_name
        time.sleep(4)
        return 1

    # Check beammask

    if beammask == 'NONE' or beammask == 'None':
        beammask = 'none'

    if beammask == 'RIGAKUCCD' or beammask == 'RigakuCCD':
        beammask = 'rigakuccd'

    if beammask != 'rigakuccd' and beammask != 'none':
        fileexists = os.path.exists(beammask)
        if fileexists == 0:
            print 'The beam.mask file was not found:',beammask
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

    # Check header image

    filexists = os.path.exists(start_head)
    if fileexists == 0 and start_head != 'none':
        print 'The header image was not found:',start_head
        time.sleep(4)
        return 1

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

    # Set values for indexing and processing images

    image_seq_int = first_image + ' ' + last_image

    first_image = int(first_image)
    second_image = first_image + 1
    first_image = str(first_image)
    second_image = str(second_image)
    image_seq_find = first_image + ' ' + second_image

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

    # Set header image to template image if not defined

    if start_head == 'none':
        start_head = image_name

    # Set any integration resolution limits

    if integrate_res != 'none':
        dt_integreso = '-reso ' + integrate_res
    else:
        dt_integreso = ''

    # Set any merging resolution limits

    if merging_res != 'none':
        dt_scalereso = '-reso ' + merging_res
    else:
        dt_scalereso = ''

    # Set the default beam mask environment variables

    if beammask == 'rigakuccd':
        os.environ['CCD_NONUNF_INFO'] = 'FirstScanImage'
        os.environ['CCD_NONUNF_TYPE'] = 'Simple_mask'
        os.environ['RX_NONUNF_INFO'] = 'None'
        os.environ['RX_NONUNF_TYPE'] = 'None'
        os.environ['A4_NONUNF_INFO'] = 'None'
        os.environ['A4_NONUNF_TYPE'] = 'None'

    if beammask == 'none':
        os.environ['CCD_NONUNF_INFO'] = 'None'
        os.environ['CCD_NONUNF_TYPE'] = 'None'    
        os.environ['RX_NONUNF_INFO'] = 'None'
        os.environ['RX_NONUNF_TYPE'] = 'None'
        os.environ['A4_NONUNF_INFO'] = 'None'
        os.environ['A4_NONUNF_TYPE'] = 'None'

    if beammask != 'rigakuccd' and beammask != 'none':
        os.environ['CCD_NONUNF_INFO'] = beammask
        os.environ['RX_NONUNF_INFO'] = beammask
        os.environ['A4_NONUNF_INFO'] = beammask
        os.environ['CCD_NONUNF_TYPE'] = 'Simple_mask'
        os.environ['RX_NONUNF_TYPE'] = 'Simple_mask'
        os.environ['A4_NONUNF_TYPE'] = 'Simple_mask'

    # Set image batch identification

    if batch_prefix == 'none':
        batch_prefix = '1'

    dtrek_prefix = batch_prefix + '_'

    # Specifically set environment variable SCAN_TEMPLATE for 3 digit image numbers to avoid d*TREK problems 

    if image_extension == 3:
        i_end = len(image_name)
        i_start = i_end - 3

        image_file_extension = image_name[i_start:i_end]

        i_start = i_end - 7
        image_file_root = image_name[0:i_start]

        scan_template = image_file_root + '???' + '.' + image_file_extension

        os.environ['SCAN_TEMPLATE'] = scan_template

    # Set final file names based on batch prefix input  

    data_code_number = batch_prefix

    # Go to frame directory

    os.chdir(image_dir)

    # Delete any residual junk from failed integration

    aList_dir = os.listdir(image_dir)
    number_files = len(aList_dir)

    count = 0
    while count < number_files:
        target_file = aList_dir[count]

        if target_file.find('dtintegrate.head.') > -1:
            os.remove(target_file)

        if target_file.find('dtscaleaverage.head.') > -1:
            os.remove(target_file) 

        count = count + 1

    # Collect user-defined beamcenter or distance data from a special constants file

    beam_x = 'none'
    beam_y = 'none'
    xtal_detector_distance = 'none'
    set_beam_center = 'none'
    set_gonio_values = 'none'

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

    ##########
    # Start  #
    ##########

    print '\nAutomated integration and merging\n'
    print 'Image directory:',image_dir
    print 'Template image:',image_name
    print 'Header image:',start_head
    print 'Beam mask:',beammask
    print 'Batch number:',batch_prefix
    print 'First image to process:',first_image
    print 'Last image to process:',last_image
    print 'Images for initial indexing:',image_seq_find
    print 'Integration resolution limits:',integrate_res
    print 'Merging resolution limits:',merging_res
    print 'Expected space group number:',dt_spacegroup

    if beam_x != 'none' and beam_y != 'none':
        input_beamcenter = ' -beamcenter ' + beam_x + ' ' + beam_y + ' '
        print 'Using input beam center:',beam_x,beam_y
    else:
        input_beamcenter = ''
        
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
        start_head_base = os.path.basename(start_head) 

        ##########################################################################################################
        # Template code for changing image header hardware machine constant parameters - currently only distance #
        ##########################################################################################################

        if xtal_detector_distance != 'none':           

            runtime = time.ctime(time.time())
            print 'Date:',runtime
            print 'Correcting file header'

            new_start_head_base = 'new_' + start_head_base
            backup_start_head_base = start_head_base.replace('.img','.backup_img')
            backup_start_head_base = backup_start_head_base.replace('.osc','.backup_osc')       

            file = open('dtheader_edit.inp','w')
            file.write(start_head_base)
            file.write('\n')
            file.write(new_start_head_base)
            file.write('\n')

            if xtal_detector_distance != 'none':
                set_gonio_values = 'CCD_GONIO_VALUES=0.0000 0.0000 0.0000 0.0000 0.0000 ' + xtal_detector_distance
                file.write(set_gonio_values)
                file.write(';\n')

                set_gonio_values = 'A4_GONIO_VALUES=0.0000 0.0000 0.0000 0.0000 0.0000 ' + xtal_detector_distance
                file.write(set_gonio_values)
                file.write(';\n')

            file.close()

            run_dtheaderedit = dtheaderedit + ' < dtheader_edit.inp > dtheader_edit.log'
            os.system(run_dtheaderedit)

            fileexists = os.path.exists(new_start_head_base)
            if fileexists == 0:
                print 'dtheaderedit run to change header data failed'
                time.sleep(4)
                return 1
            else:

                # Keep original backup of this image

                fileexists = os.path.exists(backup_start_head_base)
                if fileexists == 0:
                    os.rename(start_head_base,backup_start_head_base)

                fileexists = os.path.exists(start_head_base)
                if fileexists != 0:
                    os.remove(start_head_base)
                    print 'Retained original backup of image with initial header\n',backup_start_head_base
                    time.sleep(4)

                os.rename(new_start_head_base,start_head_base)

                os.remove('dtheader_edit.inp')
                os.remove('dtheader_edit.log')

        ###################################################################
        # Merge info from the image header into the starting .head file   #
        ###################################################################

        runtime = time.ctime(time.time())
        print 'Date:',runtime
        print 'Setting up starting head file'

        fileexists = os.path.exists('newmerge.head')
        if fileexists != 0:
            os.remove('newmerge.head')

        run_dtheadermerge = dtheadermerge + ' ' + start_head_base + ' ' + image_name_base + \
                            ' +SCAN\* +ROTATION\* +\*_GONIO_VALUES +SOURCE_WAVELENGTH newmerge.head > ' + dtrek_prefix + 'dtheadermerge.log'

        os.system(run_dtheadermerge)

        fileexists = os.path.exists('newmerge.head')
        if fileexists == 0:
            print 'Process dtheadermerge seems to have failed'
            time.sleep(4)
            return 1
        else:
            runtime = time.ctime(time.time())

            file = open('autoprocess.log','a')
            file.write('Header merge done           : ')
            file.write(runtime)
            file.write('\n')
            file.close()

        ##################################################
        # Find spots in some images - default parameters #
        ##################################################

        print 'Finding spots'

        fileexists = os.path.exists('dtfind.head')
        if fileexists != 0:
            os.remove('dtfind.head')

        fileexists = os.path.exists('dtfind.ref')
        if fileexists != 0:
            os.remove('dtfind.ref')

        run_dtfind = dtfind + ' newmerge.head -seq ' + image_seq_find + input_beamcenter + \
                     ' -sigma 3 -min 50 -filter 6 -window 0 0 -display -out dtfind.head > ' + dtrek_prefix + 'dtfind.log'

        os.system(run_dtfind)

        fileexists = os.path.exists('dtfind.ref')
        if fileexists == 0:
            print 'Process dtfind seems to have failed'
            time.sleep(4)
            return 1
        else:
            runtime = time.ctime(time.time())

            file = open('autoprocess.log','a')
            file.write('Spot find done              : ')
            file.write(runtime)
            file.write('\n')
            file.close()

        ############################################################
        # Autoindex the spots that were found - default parameters #
        ############################################################

        print 'Autoindexing'

        fileexists = os.path.exists('dtindex.head')
        if fileexists != 0:
            os.remove('dtindex.head')

        run_dtindex = dtindex + ' dtfind.head dtfind.ref -maxresid 3.0 -sigma 5 -spacegroup ' + dt_spacegroup + ' > ' + dtrek_prefix + 'dtindex.log'

        os.system(run_dtindex)

        ####################################################
        # Autoindex failure as indicated by no header file #
        ####################################################

        fileexists = os.path.exists('dtindex.head')
        if fileexists == 0:

            print 'Initial (default) indexing did not succeed - entering search' 

            count = 1
            count_max = 6
            while count < count_max:

                new_maxcell = 100.0 + count * 50.0
                new_maxcell = str(new_maxcell)

                print 'Trying maximum cell:',new_maxcell

                run_dtindex = dtindex + ' dtfind.head dtfind.ref -maxresid 3.0 -sigma 5 -maxcell ' + new_maxcell + ' -spacegroup ' + dt_spacegroup + ' > ' + dtrek_prefix + 'dtindex.log'

                os.system(run_dtindex)

                fileexists = os.path.exists('dtindex.head')
                if fileexists != 0:
                    count = count_max

                count = count + 1

        # Final check for success of indexing trials

        fileexists = os.path.exists('dtindex.head')
        if fileexists == 0:
            print 'Multiple attempts with process dtindex have failed to find a solution'
            time.sleep(4)
            return 1

        else:        
            runtime = time.ctime(time.time())

            file = open('autoprocess.log','a')
            file.write('Autoindexing done           : ')
            file.write(runtime)
            file.write('\n')
            file.close()

        ##############################################################
        # Refine the autoindexing solution with the dtfind spot list #
        ##############################################################

        print 'Index refinement using spot list'

        fileexists = os.path.exists('dtrefine.head')
        if fileexists != 0:
            os.remove('dtrefine.head')    

        run_dtrefine = dtrefine + ' dtindex.head dtfind.ref ' + \
                       ' +CrysAll +DetAll -sigma 5.0 -rej 0.5 0.5 0.5 -cycles 500 -display -go -go -go ' + \
                       '> ' + dtrek_prefix + 'dtrefine1.log'

        os.system(run_dtrefine)

        fileexists = os.path.exists('dtrefine.head')
        if fileexists == 0:
            print 'Process dtrefine seems to have failed'
            time.sleep(4)
            return 1
        else:
            runtime = time.ctime(time.time())

            file = open('autoprocess.log','a')
            file.write('Refinement of spot list done: ')
            file.write(runtime)
            file.write('\n')
            file.close()

        ##################################################
        # Refine the autoindexing solution from images   #
        ##################################################

        print 'Index refinement using images'

        fileexists = os.path.exists('dtrefine.head.1')
        if fileexists != 0:
            os.remove('dtrefine.head.1')

        run_dtrefine = dtrefine + ' dtrefine.head -seq ' + \
                       image_seq_find + \
                       ' +CrysAll +DetAll -sigma 5.0 -rej 0.5 0.5 0.5 -cycles 500 -display -go -go -go ' + \
                       ' > ' + dtrek_prefix + 'dtrefine2.log'

        os.system(run_dtrefine)

        fileexists = os.path.exists('dtrefine.head.1')
        if fileexists == 0:
            print 'Process dtrefine seems to have failed'
            time.sleep(4)
            return 1
        else:
            runtime = time.ctime(time.time())

            file = open('autoprocess.log','a')
            file.write('Refinement from images done : ')
            file.write(runtime)
            file.write('\n')
            file.close()

        ##########################################################
        # Predict to set mosaicity where we want it to start at  #
        ##########################################################

        print 'Predict mosaicity'

        fileexists = os.path.exists('dtpredict.head')
        if fileexists != 0:
            os.remove('dtpredict.head')

        fileexists = os.path.exists('dtpredict.ref')
        if fileexists != 0:
            os.remove('dtpredict.ref')    

        run_dtpredict = dtpredict + ' dtrefine.head -seq ' + image_seq_find + \
                        ' -mosaicity 1.5000 -display > ' + dtrek_prefix + 'dtpredict.log'

        os.system(run_dtpredict)

        fileexists = os.path.exists('dtpredict.head')
        if fileexists == 0:
            print 'Process dtpredict seems to have failed'
            time.sleep(4)
            return 1
        else:
            runtime = time.ctime(time.time())

            file = open('autoprocess.log','a')
            file.write('Mosaicity prediction done   : ')
            file.write(runtime)
            file.write('\n')
            file.close()

        ##############################################################
        # Integrate the images (option values for protein crystals)  #
        ##############################################################

        fileexists = os.path.exists('dtprofit.ref')
        if fileexists != 0:
            os.remove('dtprofit.ref')

        fileexists = os.path.exists('dtintegrate.head')
        if fileexists != 0:
            os.remove('dtintegrate.head')

        print 'Integrating'

        run_dtintegrate = dtintegrate + ' dtpredict.head -window 0 0 -pad 1 ' + \
                          dt_integreso + ' -mosaicitymodel 1.000 0.000 -profit 50 7 -batch 1 4 -batchprefix ' + \
                          batch_prefix + ' -display -prerefine 2 -seq ' + image_seq_int + \
                      ' > ' + dtrek_prefix + 'dtintegrate.log'

        os.system(run_dtintegrate)

        #############
        # Clean-up  #
        #############

        fileexists = os.path.exists('dtintpartials.ref')
        if fileexists != 0:
            os.remove('dtintpartials.ref')

        fileexists = os.path.exists('dtintrejects.ref')
        if fileexists != 0:
            os.remove('dtintrejects.ref')

        fileexists = os.path.exists('dtpredict.ref')
        if fileexists != 0:
            os.remove('dtpredict.ref')

        fileexists = os.path.exists('dtfind.ref')
        if fileexists != 0:
            os.remove('dtfind.ref')   

        fileexists = os.path.exists('dtintegrate.ref')
        if fileexists != 0:
            os.remove('dtintegrate.ref')

        fileexists = os.path.exists('dtfind.head')
        if fileexists != 0:
            os.remove('dtfind.head')   

        aList_dir = os.listdir(image_dir)
        number_files = len(aList_dir)

        count = 0
        while count < number_files:
            target_file = aList_dir[count]

            if target_file.find('dtintegrate.head.') > -1:
                os.remove(target_file)

            if target_file.find('dtscaleaverage.head.') > -1:
                os.remove(target_file) 

            if target_file.find('tmp.ref') > -1:
                os.remove(target_file)      

            count = count + 1

    ###############################################################
    # Scale and average the integrated and profile-fitted reflns  #
    ###############################################################

    runtime = time.ctime(time.time())

    print 'Date:',runtime
    print '\nMerging'

    fileexists = os.path.exists('dtscale.ref')
    if fileexists != 0:
        os.remove('dtscale.ref')

    fileexists = os.path.exists('dtintegrate.head')
    if fileexists == 0:
        print 'File dtintegrate.head is not available for merging'
        time.sleep(4)
        return 1

    fileexists = os.path.exists('dtprofit.ref')
    if fileexists == 0:
        print 'File dtprofit.ref is not available for merging'
        time.sleep(4)
        return 1
    else:
        runtime = time.ctime(time.time())

        file = open('autoprocess.log','a')
        file.write('Integration done            : ')
        file.write(runtime)
        file.write('\n')
        file.close()

    output_ref = 'ScalAverage_' + data_code_number + '.ref'
    output_log = 'dtscaleaverage_' + data_code_number + '.log'
    output_ref_initial = 'ScalAverage_' + data_code_number + '_initial.ref'
    output_log_initial = 'dtscaleaverage_' + data_code_number + '_initial.log'

    run_dtscaleaverage = dtscaleaverage + ' dtintegrate.head dtprofit.ref ' + dt_scalereso + \
                         ' -sigma 5.0 -anom -errormodel -reject .0075 -batchscale -reqab spherical 4 3 ' + \
                         output_ref + ' > ' + output_log

    os.system(run_dtscaleaverage)

    fileexists = os.path.exists(output_ref)
    if fileexists == 0:    
        print 'Process dtscaleaverage seems to have failed'
        time.sleep(4)
        return 1
    else:
        runtime = time.ctime(time.time())

        file = open('autoprocess.log','a')
        file.write('Merging done                : ')
        file.write(runtime)
        file.write('\n')
        file.close()

    ####################################################################
    # (Re) Scale and average the integrated and profile-fitted reflns  #
    ####################################################################

    # Capture table data to see how initial merging went (note: depend on precise d*trek file format)

    aList_res_low = []
    aList_res_high = []
    aList_ioversigi = []
    table_length = 0
    read_table = 'no'

    file = open(output_log,'r')
    allLines = file.readlines()
    file.close()

    for eachLine in allLines:

        if eachLine.find('Rmerge vs Resolution') > -1:
            read_table = 'yes'

        if read_table == 'yes':

            aLine = eachLine.split()
            number_items = len(aLine)

            if eachLine.find('---') == -1 and number_items == 12:

                res_low = aLine[0]
                res_high = aLine[2]
                ioversigi = aLine[7]

                if res_low.find('.') > -1 and res_high.find('.') > -1 and ioversigi.find('.') > -1:
                    aList_res_low.append(res_low)
                    aList_res_high.append(res_high)
                    aList_ioversigi.append(ioversigi)

        if eachLine.find('Summary of data collection statistics') > -1:
            read_table = 'no'

    # Note that table includes overall statistics in final entry

    table_length = len(aList_res_low)
    table_length = table_length - 1
    max_table_length_index = table_length - 1
    valid_data = 'no'

    if table_length > 1:

        count = 0
        while count < table_length:

            res_high = aList_res_high[count]
            ioversigi = aList_ioversigi[count]
            res_high = float(res_high)
            ioversigi = float(ioversigi)

            if ioversigi > ioversigi_limit:
                res_to_process = res_high
                valid_data = 'yes'

            count = count + 1

        # Remerge at lower resolution if required by I/sigI criteria

        if valid_data == 'yes':

            res_lower = aList_res_low[0]
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

                dt_scalereso = '-reso ' + res_lower + ' ' + res_to_process

                run_dtscaleaverage = dtscaleaverage + ' dtintegrate.head dtprofit.ref ' + dt_scalereso + \
                                     ' -sigma 5.0 -anom -errormodel -reject .0075 -batchscale -reqab spherical 4 3 ' + \
                                     output_ref + ' > ' + output_log

                os.system(run_dtscaleaverage)

                fileexists = os.path.exists(output_ref)
                if fileexists == 0:    
                    print 'Process dtscaleaverage to remerge seems to have failed'
                    time.sleep(4)
                    return 1
                else:
                    runtime = time.ctime(time.time())

                    file = open('autoprocess.log','a')
                    file.write('Remerging done              : ')        
                    file.write(runtime)
                    file.write('\n')
                    file.close()   
        else:

            print 'Warning - parsing to test for remerging failed'

    ######################
    # Tail end processes #
    ######################

    # Clean-up processing

    fileexists = os.path.exists('zone.ref')
    if fileexists != 0:
        os.remove('zone.ref')

    fileexists = os.path.exists('dtscaleaverage_rejects.ref')
    if fileexists != 0:
        os.remove('dtscaleaverage_rejects.ref')

    fileexists = os.path.exists('newmerge.head')
    if fileexists != 0:
        os.remove('newmerge.head')

    fileexists = os.path.exists('dtpredict.head')
    if fileexists != 0:
        os.remove('dtpredict.head')

    fileexists = os.path.exists('dtrefine.head')
    if fileexists != 0:
        os.remove('dtrefine.head')

    fileexists = os.path.exists('dtrefine.head.1')
    if fileexists != 0:
        os.remove('dtrefine.head.1')

    fileexists = os.path.exists('dtscaleaverage.head')
    if fileexists != 0:
        os.remove('dtscaleaverage.head')

    fileexists = os.path.exists('dtindex.head')
    if fileexists != 0:
        os.remove('dtindex.head')

    # Transfer the final reflection and output files to predefined space or push to default location folder

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

        aList = os.path.split(image_dir)
        image_folder_name = aList[1]
        destination_folder = image_folder_name + '_processed'
        destination_folder_path = os.path.join(image_dir,destination_folder)

        fileexists = os.path.exists(destination_folder_path)
        if fileexists == 0:
            os.mkdir(destination_folder_path)

        output_ref_destination = os.path.join(destination_folder_path,output_ref)
        output_log_destination = os.path.join(destination_folder_path,output_log)

        fileexists = os.path.exists(output_ref_destination)
        if fileexists != 0:
            os.remove(output_ref_destination)

        fileexists = os.path.exists(output_log_destination)
        if fileexists != 0:
            os.remove(output_log_destination)

        os.rename(output_ref,output_ref_destination)
        os.rename(output_log,output_log_destination)

    # Log and clean-up

    runtime = time.ctime(time.time())

    print 'Date:',runtime
    print '\nIntegrated and merged intensity file:',output_ref_destination
    print 'Merging statistics file:',output_log_destination

    time.sleep(4)
    return 0


if __name__ == "__main__":
    sys.exit(Run())
