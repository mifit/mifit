# Interfaces to CCP4 software #

**9.1 The Job/Integrate interface**

![http://mifit.googlecode.com/svn/wiki/images/image082.jpg](http://mifit.googlecode.com/svn/wiki/images/image082.jpg)

Figure 9.1 The **Job/Integrate** interface

Automated processing of image data may be performed using the **Job/Integrate**command, which runs MOSFLM and CCP4/SCALA. This interface is not intended to replace traditional image data processing interfaces for analyzing data on completely new or difficult-to-process crystals; it is intended to automate processing of routine reasonable data sets obtained in high throughput co-crystallography and fragment screening modes.

The **Intensity Data**option is used to select an image file (usually the first image and with file extension _.img_ or _.osc_) from the set of images to process. This image serves as the template for determining which image files belong in the processing set.



Since this GUI is intended to relieve the crystallographer from the task of processing data from previously characterized crystals the expected space group is known and should be entered as the **Space Group**parameter. No other information is needed to auto-process image data.

In some cases the hardware parameters for the beam center or crystal to detector distance encoded in the headers of the image files are not correct. These values may be rectified by providing a parameter file using the **Reset Hardware Parameters**option. The format of this file is simply a keyword for the parameter that needs to be changed followed by the user defined values.

For example,

beam\_center 1013.0000 1009.0000

xtal\_detector\_distance 50.0000

In some cases it may be necessary to process only a subset of the available images (for example, if the crystal became severely radiation damaged during data collection) or process to a predetermined resolution limit (for example, if the detector area is much larger than the extent of the diffraction pattern). For these purposes the **First Image**, **Last Image**and **Resolution range**parameters are available. Image numbers may be entered in the **First Image**and **Last Image**parameters. The **Resolution range**parameters are low and high resolution limits in angstroms.

The **OK** button is used to launch data integration and merging.

Data integration and merging usually takes several minutes. After completion the standard output files _ScalAverage\_1.mtz_and _scala\_scaleaverage\_1.log_contain the integrated merged intensity data and merging statistics.

**9.2 The Job/Molecular Replacement interface**

![http://mifit.googlecode.com/svn/wiki/images/image084.jpg](http://mifit.googlecode.com/svn/wiki/images/image084.jpg)

Figure 9.2 The **Job/Molecular Replacement** interface

Molecular replacement calculations may be run using the **Job/Molecular Replacement**command. This application employs the CCP4/MOLREP and CCP4/PHASER programs to carry out molecular replacement calculations involving complete rotation-translation searches.

Required parameters for these calculations are the **Working directory**parameter, which defines where the molecular replacement calculation outputs will be created, **Model (PDB)**, which defines the input search model and **Data**, which defines the input structure factor data. Note that data may be entered as structure factor amplitudes in mtz format or as merged intensities in d\*TREK (._ref),_ SCALEPACK (._sca_) or mtz (_.mtz_) formats.

The **Fixed model**checkbox allows the user to enter a partial model to hold fixed while carrying out a search with the model defined by the **Model (PDB)**parameter. This might be useful for solving a structure with two different types of molecule in the crystal asymmetric unit where it is necessary to locate the two molecules independently.

If there are multiple molecules of the same type, the **Copies** option may be used to set the number of copies that should be searched for. This required number of copies might be estimated using a Matthews coefficient analysis.



The **Space group**option may be set at **Default**to use the space group encoded in the header of the diffraction data file or a **Space Group Number**may be entered.

The **Method**option may be set to **MolRep**or **Phaser**to select between CCP4/MOLREP and CCP4/PHASER as the molecular replacement engine. At the present time, CCP4/PHASER is a more sensitive and reliable program for solving new structures but CCP4/MOLREP is faster and might be used on simple cases, involving known structures.



If the **Search multiple models**checkbox is selected molecular replacement processes will be run using each of the PDB files that are found in the same directory as the model defined by the **Model (pdb)**parameter. This option may be useful when a new structure is being solved and it is unclear which of several homologous structures would be the best candidate for solving the structure by MR. The results from multiple searches can be assessed by looking at the project history file.

The **Match input position**checkbox may be used to apply crystal symmetry operations to place the model resulting from the MR process as closely as possible to the input model. This option is useful when working on co-crystal or other projects where there are precedent structures in the same space group since it is usually convenient to have all structures placed at the same position in the crystal cell.

The **OK**button is used to launch the molecular replacement application. Output files resulting from this application are generated in the working directory and identified by file root _molrep_#_where the # is the serial number corresponding to this molecular replacement run. The project history file may be consulted in order to identify the next available sequence number and find the input, output and summary of each molecular replacement run._

**9.3 The Job/SAD Phasing interface**

![http://mifit.googlecode.com/svn/wiki/images/image086.jpg](http://mifit.googlecode.com/svn/wiki/images/image086.jpg)

Figure 9.3 The **Job/SAD** phasing interface

The calculation of refined electron density maps with CCP4/PHASER, given intensity data and a set of sites (for example, previously located with SHELXD) is supported through the **Job/SAD Phasing** interface.

The **Working directory**parameter is the directory in which you wish the SAD phasing outputs to appear. This directory may be located or created via the associated **Browse**button.

The **Intensity data**parameter should contain the path to a file containing merged intensity data, in which the Bijovoet reflection mates have been kept separate. Input data from CCP4/SCALA, SCALEPACK or d\*TREK should be automatically recognized. After performing an internal conversion to the CCP4 MTZ file format, these data are processed through the CCP4/TRUNCATE program prior to carrying out phasing calculations.

The **Scatterer sites**parameter is the file containing anomalous scattering sites. The **Files of type**filter in the browser contains a pull down for PDB (._pdb)_ or RES (_.res_) input files, which correspond to PDB and SHELX formats.

The **Scatterer type**parameter may be used to set the element code for the anomalous scattering centers. The default value is S (sulfur) since this phasing system has been successfully tested on examples of sulfur-SAD phasing from a chromium source. SE (selenium) would also be a common setting.

The **Number of sites**parameter is the number of scattering sites that will be used.

The **Solvent fraction**parameter is the expected solvent content of the crystal and is used for the density modification (phase refinement) step for map improvement.

The **Switch site hand**option is used to generate phases the inverted constellation of scatterer sites, since this is initially unknown and both possible hands may be tried.

The **Change spacegroup to**parameter may be used if it is necessary to change the space group from the value given in the input file of intensity data. It may be necessary to use this parameter to test phasing in enantiomorphic space groups or for cases where the processing only determined the point group.



The **OK**button launches the phasing procedure. It will typically take a few minutes to phase a protein.

**9.4 The Job/Refinement interface**

![http://mifit.googlecode.com/svn/wiki/images/image088.jpg](http://mifit.googlecode.com/svn/wiki/images/image088.jpg)

Figure 9.4 The **Job/Refinement** interface

Automated refinement jobs with error detection may be run using the **Job/Refinement**command.

The **Working directory**parameter specifies the location where the output files resulting from the refinement will be generated.

The **Model (PDB)**parameter specifies the input model and the **Data**parameter specifies the refinement data file in mtz format.

The **Refinement Type**options include **Rigid Body**, which performs rigid-body refinement of the independent chains within the coordinate set, and options to perform refinement on all atomic parameters with either **Refmac5**or **SHELX**.

Individual atomic temperature factors may be refined with either **Isotropic**or **Anisotropic**restraints.

The refinement weight applied to the x-ray term in restrained refinement may be specified using the **Weight**parameter. In practice, the optimal weight tends to larger values with higher resolution data and also tends to diminish as the refinement approaches convergence.

The **Number of cycles**parameter specifies the number of refinement cycles within a refinement run. A value of 5-15 cycles is often sufficient for preliminary refinement of a co-crystal structure.

The **Water picking cycles**parameter intersperses water-picking with refinement cycles. For example, if this parameter is set to a value of 3 then electron density maps will be scanned for water molecules following each of three refinement runs and the process will then be completed by one further refinement run.

The optional items that may be supplied through this interface include **Max resolution**, which allows the user to truncate the upper resolution of data used in the refinement.

The **TLS Specification**parameter may be used to input specifications of bodies to use for TLS refinement (see the TLSIN format as specified by the CCP4 documentation).

The **Dictionary**parameter may be used to provide a restraint dictionary for an input structure that contains a novel small molecule ligand.

The **OK**button is used to launch the refinement application.

Output files resulting from this application are generated in the working directory and identified by file root _refine_#_where the # is the serial number corresponding to this refinement run. The project history file may be consulted in order to identify the next available sequence number and find the input, output and summary of each refinement run. Output files retained from the refinement process using CCP4/REFMAC5 include coordinates, phased reflection data, the standard CCP4/REFMAC5 summary files in free text and mmCIF format and the error list file._

**Note: Automated structure validation and the 'error list'**

Structure refinements run via the **Job/Refinement**and the **Job/Cocrystal Solution**interfaces automatically create an 'error list' text file that reports certain global structure quality values and a list of amino acids in which a severe abnormality was detected.

This file may be loaded as an annotation for a given structure to provide a display on the main canvas of the error type and site. Select the **Display**tab in the navigation tree, right click on the **Annotations**icon and select **Import from error list**to load this list. In addition the errors are encoded in the header of the PDB file and the command **Auto-show imported errors** may be used to display/conceal them from the main canvas. (Often you may prefer to conceal this display!)

Global quality metrics include the standard crystallographic indices R<sub>work </sub>and R<sub>free </sub>as calculated by CCP4/REFMAC5, the percentage of residues outside the core area of the Richardsons Ramachandran analysis (treating general, Pro and Gly amino acids as separate cases), the percentage of residues with severely abnormal χ1 angles (using the Richardsons side chain torsion angle data as implemented in CCP4/ROTAMER program) and the total number of amino acids that appear to contain local errors. These metrics give a good indication of the quality of a structure model at a particular point in the refinement. The text file is concise enough to be printed and included in a laboratory notebook.

More useful for achieving model improvement is the information from the list of local errors, which may be used as targets for model examination and rebuilding efforts. Although proteins do contain some genuine structural anomalies (with interesting energetic rationales) and there is a grey zone in assessing what degree of misfit between data and map is acceptable, the calculation routes and thresholds used to identify most of these errors types have been tested on very large numbers of structures and will usually indicate portions of the model that need to be checked and corrected. The statistical data used for these tests and calculation routes are more sophisticated and modern than standards embodied in the PROCHECK (1993) program and the density checks supply a means of detecting geometrically correct features that misfit the data.

More specifically, the criteria for flagging structure errors are:

Geometry errors  bonds length deviations greater than 6x the refinement sigma or bond angles greater than 8x the refinement sigma used by REFMAC5

Van der Waals errors  close contact deviations greater than 4.25x the refinement sigma used by REFMAC5

Omega errors  deviations of 4x the true deviation (5.6) from the peak (178.9)

Phi-Psi errors  lies outside a prescribed area containing 99.95% of correct amino acids in the general area, 99.8% of correct amino acids for GLY and 99.8% of data for PRO in the Richardson tabulations

Cis peptide errors  any non-PRO amino acid flagged as cis by REFMAC5

Rotamer errors in χ-1  deviations of more than 45 in χ-1 from the nearest conformer as calculated by ROTAMER, using rotamer data from the Richardson lab

Density errors  identification of a 5 sigma peak/hole within 1.5 of any model atom. This threshold may need adjustment.

With the possible exception of density errors, a fully refined structure will have very few residual abnormalities.

**9.5 The Job/NCS Modeling interface**

![http://mifit.googlecode.com/svn/wiki/images/image090.jpg](http://mifit.googlecode.com/svn/wiki/images/image090.jpg)

Figure 9.6 The **Job/NCS Modeling** interface

When building models of protein crystals in which there is more than one subunit per asymmetric unit it is usually the case that model fitting efforts are concentrated on one copy that is defined by a single chain id in the PDB coordinate file. In the early building stages it would be desirable to replicate this rebuilt model across all copies of the protein. The **Job/NCS Modeling**interface is designed to provide a convenient mechanism for this task.

This interface is available when a protein model is loaded into the MIFit canvas. Required parameters are the **Primary Chain ID**, which defines the template chain (the model that the user wishes to replicate over the less completely modeled chains in the structure) and the **Working Directory**where the rebuilt model will be placed.

The optional **Model Exclusions**parameter may be used provide the path to a file that contains data lines for chain-id, start-amino-acid-number, end-amino-acid-number over which the NCS symmetry will not be applied. For example, B 120 130, would exclude chain B, amino acids 120-130 from reconstruction.

When launched by selecting **OK**the rebuilt model is returned to the MIFit canvas.

Optionally, a diffraction data file may be entered with the **Data (mtz)**parameter and the density map may then be averaged within the model envelope to provide phases according one of the **Average Calculated Phases**, **Average Experimental Phases**and **Average Combined Phases**selections. These are all somewhat experimental options. The process will automatically calculate phases but expects the input file to contain Hendrickson-Lattman ABCD coefficients to make use of experimental phases.

The **Mask Additions**option provides a specification (a file containing lines with chain-id, amino acid-number) for including atomic groups outside the primary chain in the subunit mask that is used for averaging. The output mtz file is named after the input pdb file but prefixed with _ncs__._

**9.6 Automated multiple protein:ligand structure superposition**

** **

**http://mifit.googlecode.com/svn/wiki/images/image092.jpg*****

Figure 9.6 The **Job/Cocrystal Superposition** interface

A very common scenario in protein crystallography is that several closely related structures (with the same sequence) have been solved and the structure analyst wants to compare the binding modes of many different ligands within the active site. To perform an accurate comparison of ligand coordinates, only atoms around the ligand binding sites should be used for the structure superposition. The method used here is to base the structure superposition on Cα atoms within 15 of an active site coordinate within the reference structure.

The **Job/Cocrystal Superposition**command is used to run the structure comparison.

This command is only accessible when a model is loaded into the MIFit main canvas and it will often be convenient to first load the reference structure that will be used as a basis for the structure superposition.

The **Working Directory**parameter defines the path to the directory in which the results from the structure comparison will appear.

The **Structure Directory**should contain all the coordinate files (PDB format) for the protein:ligand structures that will be compared. These structures must have consistent sequence naming conventions. _i.e._ residue A 22 in one structure corresponds to residue A 22 in all the others.

The **Target File**parameter defines the particular coordinate file that will be used as the reference structure for comparisons. This file is usually one of the structures in the **Structure Directory**but need not be.

The **Target Coordinates**are the **X**, **Y**and **Z**positions in angstroms in the **Target File**for the center of the ligand binding site. This position is used to select the surrounding Cα atoms for structure superposition.

After clicking on **OK**the superposition process will return a coordinate file called _allligands.pdb_to the **Working Directory**and load this file into the MIFit main canvas. This file contains all of the ligands with the structure set superimposed together and distinguished by unique chain ids.

**9.8 The Job/Refmac5 restraints interface**

Protein structures that contain cofactors or other small molecule ligands may require the input of restraint information for their refinement. The CCP4/REFMAC5 installation contains dictionaries in mmCIF format for most of the small molecule entities found in public domain structures. Many of these dictionaries are in a form such that the ligand stereochemistry is automatically known to CCP4/REFMAC5 (_i.e._ in the same way that the standard amino acids are known) but others are in a minimal form that and require preprocessing to a complete description before being used. In addition, some structures contain several novel ligands, or a combination of known and novel ligands. It will not always be obvious before attempting to run a refinement whether all of the entities in a coordinate file are accounted for and that the refinement will be able to proceed.

In order to check a coordinate file and, if necessary, combine various sources of restraint information, the **Job/Set Refmac5 restraints**command is available. The input for this command is the coordinate file that you wish to use for refinement and, in the same directory as that file, any mmCIF dictionary files that are considered necessary for the refinement. The dictionary files should be named according to the three character residue code in the coordinate file with file extension _.cif_. For example, if the structure contains a novel ligand MMM there should be dictionary file available called _MMM.cif_.

The results from running this pre-refinement check are reported in the project history file. If a restraint dictionary was created that contains the information from several inputs it is named _restraints_#.cif_, where # is a serial number reflecting the number of times this command was executed in the current working directory._

Outputs from this process are most conveniently viewed via the **Job List**control in the document tree pane