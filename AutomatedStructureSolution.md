# Automated Structure Solution #

**10.1 Overview**

An automated structure determination process has been developed as a special application within _MIExpert_ to provide pre-refined structures and maps for co-crystal structure analysis in as efficient and convenient a manner as possible.

**Process Summary**

** **

This process prepares refinement data from files of integrated intensities (D\*TREK, SCALEPACK and SCALA (MTZ) formats are supported), runs molecular replacement with CCP4/MOLREP and performs preliminary refinement calculations (including water picking) with CCP4/REFMAC5 to obtain a pre-refined model of the co-crystal structure. The interface and associated script allow for the entry of multiple data sets in a single job.

A MIFit session file is created at the completion of the refinement process with pre-computed map data and with a view centered on the expected ligand site. An option is available to provide an HTML report summarizing results from a set of co-crystal solution jobs on multiple structures. Summary information on the structure solution process is provided in the file _project\_history.txt_ in each working directory. A convenient list of probable structure errors is also provided for each structure. The entire automated process will normally complete within a few minutes per structure. For example, on a 2GHz laptop, it takes 3.5 minutes to complete this process on a 2 resolution data for hsp90.

The application may also start with unprocessed image data and will perform the preliminary image data processing and merging. Obviously, image data processing is a relatively long operation and will typically add ~10 minutes to the run time.

**Reference data and structure transformations**

** **

If a reference data set is provided the process will re-index the new data sets (for applicable space groups) to find the indexing that is most consistent with the reference. Following molecular replacement, the script will also apply symmetry operations to place the resulting models as close as possible to the input model. The aim of these operations is to ensure that all structures within a series of co-crystals are conveniently placed at the same position in the crystal cell. The cross-validation flags in the reference data set are preserved in the setup of refinement data for the new structure.

**Molecular replacement and refinement**

If the structure is known to exist in multiple conformational states, a set of possible starting models may be automatically input to the molecular replacement process. In the case of searches involving multiple trial models, the model that gives the lowest R-factor following molecular replacement is carried forward for further refinement.

The model refinement is run in two stages, initially with the model directly from molecular replacement and then by four short runs that are interspersed with water-picking calculations. All data are used for refinement calculations together with the 'mask' bulk solvent correction and isotropic B-factor refinement. Should a special restraint library be needed for the refinement (_i.e._ if the MR model contains non-standard components) it may be supplied as part of the input.

**Ligand placement**

An interface option is available to support ligand placement. A rigid ligand 6D search may be used if the ligand conformer can be anticipated (for example, in the case of fragments or ligands closely related to previous examples) although this is somewhat slow. The flexible search option may be used if other ligand-fitting technology is available and has been implemented in the MIExpert_mi\_bng.py_ driver script.

**10.2 The automated co-crystal structure determination interface**

![http://mifit.googlecode.com/svn/wiki/images/image094.jpg](http://mifit.googlecode.com/svn/wiki/images/image094.jpg)

Figure 10.1 The **Job/Cocrystal Solution** interface

The **Job/Cocrystal Solution**interface automates co-crystal structure solution. The minimum requirements for running this application are a molecular replacement search model and a file of intensity data.

The molecular replacement search model may be entered as the **Model (pdb)**parameter.

The filenames for the data set(s) to be analyzed may be entered under the **Intensity Data**parameter using the **Add**button. An important and useful feature of this interface is that a series of related data sets may be added to the **Intensity Data**list and the structure solution operations will work on each of them in turn. If a data set is incorrectly included in the processing list it may be removed with the **Remove**button.

For this application the working directory is assumed to correspond to the location of the data set_. i.e._ it is best to place each set of intensity data within a separate directory. If the extension for the data set is _.img_ or _.osc_ then it is assumed that this corresponds to an image data file and the application will attempt to process all image data in that directory. A subdirectory called 'BNG' is created as the working directory for all subsequent steps in this case.

A set of optional menu items are available for these structure solution jobs.

The **Image Data**option may be set to **Preprocessed**if the input file is merged integrated data or **Process**if the user wishes to process the image data as part of the automation process. Automated image data processing is somewhat less robust than the rest of the structure solution process although it is practical for good quality data.

The **Detector Constants**option is only applicable if the structure solution job includes automated image data processing and information (usually the beam center or crystal to detector distance) is incorrect. The browse button associated with this option may be used to select a text file (_.txt_) that may contain values for the beam center (beam\_center) in mm and the crystal to detector distance (xtal\_detector\_distance) in pixels.

For example,

beam\_center 1013.0000 1009.0000

xtal\_detector\_distance 50.0000

The **Spacegroup**option may be set at **Unchanged**to use the space group encoded in the header of the diffraction data file. If this value is incorrect but the data was processed in the same point group the correct **Space Group Number**may be entered. The space group number must be provided if this application is to run image data processing.



Although (arguably) not strictly necessary, it will often be preferred to provide a **Reference Data (mtz)**parameter to apply the indexing conventions and R-free flags used in this reference data to the new data set. In some space groups the provision of this data set is necessary to ensure that the final model is placed in the crystal at approximately the same position as the input search model.

A potentially useful option is the provision of a **Search Multiple Models**checkbox which, when toggled on, will perform MR calculations using all PDB files in the same directory as the model specified by the **Model (pdb)**parameter. The model that gave the best MR solution (lowest R-factor) is carried forward for subsequent refinement. This capability is useful when a co-crystal is known to exist in multiple crystal forms (say, an open and a closed conformation) as the MR process will select the most appropriate model for subsequent refinement and model-fitting.

The **HTML Summary**option is a path to a directory in which a summary of the results from a set of structure solution jobs will appear.

The **Dictionary (cif)**parameter may be used to provide any CCP4/REFMAC5 dictionary that might be required to run refinement on the input model. For example, this would be the case if the molecular replacement model incorporated coordinates for a novel ligand or cofactor.

The **View Point**option with the **Points**parameter toggled on requires **X**, **Y**, **Z**as parameters corresponding to the view point in angstroms. This option provides an input for centering the view in the resulting MIFit session file. Alternatively, with the **Session**parameter toggled on, the required input is a previous session file. In that case both the view center and other display attributes (view direction, slab thickness, map colors etc) from that file are applied to the new session file. These options also eliminate water molecules from the vicinity of the target point in order to facilitate the generation of clearer difference density maps.

The **Place Ligand (pdb)**parameter may be used to provide ligand placement near the **View Point**target as a final step in the structure solution process. Unless ligand docking software is installed and implemented the available option is purely a density docking process by 6D search _i.e._ it adds the ligand to the protein but does not attempt to refine it. For this reason, no dictionary is required for the ligand. The rigid-body docking is mainly useful is the ligand conformation can be anticipated in advance or the ligand has few degrees of freedom (_i.e._ for fragment screening). The ligand coordinates are taken from files specified by the browsers for the two ligand search options.

Once the **OK**button is clicked the automated structure solution begins.

The model resulting from this process is _refine\_3.pdb_with associated data file _refine\_3.mtz_. If the 6D ligand placement option was applied a coordinate file _refine\_3\_ligfit.pdb_, which contains the fitted ligand, is also created. A MIFit session file (_bng\_milaunch.mlw_) that displays the final model in the context of the final electron density map and the final difference map is available for facile review of results. Local abnormalities in the final structure are listed in file _refine\_3\_errors.txt._ Crystal symmetry operations are applied to superimpose the final refined model at the same point in the crystal as the input model.

As co-crystal solution jobs complete the emerging session files may be viewed using MIFit. On a moderately powerful laptop computer, session files can usually be reviewed concurrent with ongoing structure solution processes.