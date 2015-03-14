# 8 Running external crystallographic software #

MIFit has the capability of running external crystallographic programs (principally from the CCP4 suite) through a command-line library of python 'smart scripts' contained in the MIExpert distribution. Since the Python scripts are fully exposed, these applications may be customized or adapted at individual sites and may be easily updated if the underlying CCP4 software changes between MIExpert releases.

MIExpert protects the user from common problems with using the CCP4 suite (for example, handling paths containing spaces on Windows systems, some atom name justification issues in PDB files). MIExpert also manages various data management tasks by daisy-chaining different programs to perform common tasks. For example, the automated structure solution application includes both molecular replacement and refinement steps and may be initiated directly from image data or integrated intensity data.

All of these applications may be launched from the GUIs within the MIFit **Job**menu once MIExpert is installed. The interfaces themselves are written in Python, using the pyQt module and might be altered locally if desired. For example, it is perfectly possible to add additional parameters to an interface or even provide entirely new major functions. For example, different molecular replacement programs could be added as options to the molecular replacement interface.

**8.1 Accessing the CCP4 suite**

MIExpert uses the CCP4 software suite as a basis for SAD phased map calculation, molecular replacement and structure refinement. To execute the MIExpert scripts the users operating environment should contain the CCP4 environment variables and the CCP4 programs should be in the users path.

On the Windows operating system the relevant environment variables should always be automatically available once the CCP4 suite is installed. User environment variables may be checked by (i) a right mouse click on the **My Computer**icon on the desktop, (ii) selecting **Properties**, (iii) selecting **Advanced**, (iv) selecting **Environment Variables**.

On Linux systems the window in which MIFit was launched will need to know the CCP4 environment. This will not be the case unless execution of the CCP4 setup process (_i.e._ by sourcing the CCP4 setup script) is incorporated into the users login process. You can check to see which environment variables are established in a particular window by using the Linux/UNIX _printenv_command.

**8.2 The Job List menu**

** **

MIFit interface contains a **Job List**that may be accessed via the **Jobs** tab at the bottom of the navigation tree. This list contains all the jobs that have been run in the current MIFit session. Each job is identified by a serial number. A right mouse click on a job provide a set of options, **Delete Job**, **Job Properties**, **Show Log File**, **Clean Successful**, **Clean All**and **Detach Job**. Jobs are color coded by yellow (running), green (successfully completed) or red (failed).

The **Show Log File** option shows all of the diagnostic information from the job that would normally be printed to a terminal. Note that this command is only accessible after the job is completed in order to avoid interfering with a running job. If a job should unexpectedly fail it is often possible to find the cause (or at least, the point of failure) by checking this log.

Occasionally a job is killed or otherwise fails but MIFit continues to list the job as active. This is potentially a problem as the number of simultaneous jobs is limited through the **File/Preferences/General**menu. The **Detach Job**command may be used to close this job.

The **Open Results**option provides a convenient way to open certain files from a completed job. The dialog box lists available session, PDB, and MTZ files in reverse chronological order, since in most cases the most recent files are the results of the job. The user may toggle on or off the types of files to be loaded in to MIFit with the check boxes on the left of the dialog box.

**8.3 The project history file**

Many of the scripted applications described in chapters 9 and 10 employ the concept of a history file (called _project\_history.txt_) that is created in and accessed from the working directory in which the calculations were performed. This file logs the key input and output files from each application run together with a summary of the results. This file may make a useful record, for example, of a series of molecular replacement calculations starting from different models.

To maintain accurate and consistent book-keeping, standard file root names and serial numbers are used to name the output files resulting from these scripted applications. Besides providing an indexed record of the script runs, the project history file ensures that the outputs from program runs are named sequentially and that files are not overwritten. Information is appended to the project history file upon job completion _i.e._ failed jobs do not result in an entry.

**8.4 Running custom scripts**

![http://mifit.googlecode.com/svn/wiki/images/image080.jpg](http://mifit.googlecode.com/svn/wiki/images/image080.jpg)

Figure 8.1 The custom job interface

The MIFit **Job/Run Custom Job**interface was developed to provide a very simple but general mechanism for running external crystallographic applications and reloading any resulting coordinate (pdb) or diffraction data (mtz) files in to a new MIFit session. This interface provides the user a route for applying a custom script-driven procedure to a model without leaving the MIFit interface.

The **Job Name**parameter is used to provide a temporary identification of the job.

The **Command**parameter specifies the command-line command that the user wishes to execute. The exact form of the command will depend on the operating system and type of script. For example, for a C-shell script on a Linux system running the C-shell the parameter would just be the path to the script. For a python script on Linux it might be necessary to include the path to the python interpreter in the command.

The **Arguments**parameter provides the input to the script. By supplying the arguments _$model_and _$data_paths to the model (supplied by the **Model**parameter) and the diffraction data file (supplied by the **Data**parameter) this information will be passed as the first two command-line arguments of the script. (The quotes are included here to handle paths containing spaces and might be unnecessary on a Linux file system.). For the **Model**parameter either the current active model loaded into MIFit or a model from a file may be passed to the script. Additional arguments could also be included as part of the **Arguments**parameter, as needed by the script.

The **Working Directory**parameter specifies the directory in which calculations will be performed and the output files are expected to appear.

The parameters used by the **Job/Run Custom Job**interface are preserved across MIFit sessions.