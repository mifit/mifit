# Tutorial Lesson: Automated co-crystal structure determination #

(This tutorial lesson requires that MIExpert is installed and the system has access to the CCP4 software suite)

_Scenario: You are collecting many data sets on related protein:ligand complexes as part of ligand discovery and optimization project. You wish to obtain pre-refined structures and maps to screen these data as rapidly and as conveniently as possible._

1. Copy the files automation\_start.pdb and automation\_start.ref from the examples directory of your MIFit installation into a new directory. File automation\_start.pdb is the initial search model that will be used for the structure determination  it is a well refined example of this protein structure. File automation\_start.ref is a file of intensity data (in this case, processed using the d\*TREK program).

2. Select Job/Cocrystal Solution from the MIFit interface. You will see the Cocrystal Solution menu appear.

3. Click on the Add button in the lower left of the menu and use the browser to locate and load the file of intensity data, automation\_start.ref. Note that data merged using any of d\*TREK, SCALEPACK or CCP4/SCALA may be used by. If you had multiple data sets corresponding to the same type of crystal then multiple co-crystal solution jobs could be run in series by adding them all into this selection.

4. Use the Browse... button to load the model _automation\_start.pdb_as the **Model (pdb)**parameter.

5. An option that is useful when multiple data sets are loaded into the co-crystal solution menu is HTML Summary. This provides a summary table of all the structure determination jobs. To see what this looks like, select the HTML Summary checkbox and use the associated Browse... button to determine the directory in which this file will be deposited - select the directory established in (1) for the data.

6. Deselect the **Reference Data (mtz)** and **Dictionary (cif)** menu items as there is no reference data and the model does not contain any ligands requiring user-defined parameters.

6. A useful option is to set the position on which the display will be centered in the resulting MIFit session file. To do this, check the Viewpoint selection, make sure that Coordinates is selected and enter values 61, 33 and 26 as X, Y and Z parameters. These numbers are the position of the target site on which to center the display in angstroms. If a previous session file were available that specified this view it could have been entered through the Session File (mlw) parameter. Besides setting the viewpoint, this option will also remove waters from the target site area (potentially occupied by a ligand) in order to provide more useable difference map for ligand fitting.

7. Select OK. The job runs molecular replacement and several cycles of refinement and water-picking. On a 2GHz laptop the job takes about 3.5 minutes.

8. Select the Jobs tab near the bottom of the navigation tree and use a right mouse button click on the bottom entry in the Jobs List to select Show Log File. This action will provide a view of summary information from the structure solution process. Note that this command is grayed-out (inaccessible) until the job completes. Completion is indicated by the co-crystal solution job icon changing from yellow to green.

9. Use the File/Open models, data, maps, etc... command to load the session file _bng\_milaunch.mlw_. The refined model from this run will be loaded as well as the standard likelihood-weighted electron density map and a likelihood weighted difference map. In this case default MIFit contour levels and colors are used but if the **Session File (mlw)** option had been used the display would have assumed the values from that file. You may wish to right-click on the map icons in the tree on the left and select **Show/Hide**to hide one or other of these maps. Looking at these maps you will conclude that there are several ordered water molecules but no bound ligand. This process allowed you to evaluate this data set with minimum of work. (Alternatively you may load the final model and mtz data file as described in the previous lesson on refinement).

10. Look at the directory specified in (1). Files specified by root refine\_3 correspond to the last refinement run. File refine\_3\_errors.txt lists amino acids that may be in error (as indicated by abnormal geometries or a significant degree of mismatch to the elec-tron density). Note that the job history and this error list are linked to the HTML summary file.

11. You may select File/Close or File/Exit to close this session or shut down MIFit.