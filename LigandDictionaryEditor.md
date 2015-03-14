# The Ligand Dictionary Editor #

The MIFit software contains specialized tools to simplify the task of generating refinement restraints for novel ligands. Ligand restraint files may be exported for use with the CCP4/REFMAC5 refinement software or converted to other formats for use with SHELX or CNS/CNX.

**6.1 Entering ligand data**

The **Dictionary/Import Ligand**command provides several paths for entering small molecule information into the MIFit dictionary. This command branches into **CIF**, **Mol**, **Pdb**and **Smiles**options, which may be used to enter small molecule information in any of these formats.

In order to develop information for a novel ligand an effective approach will often be to draw the ligand with a chemistry program containing a molecular sketcher (for example, ChemDraw or ChemSketch, which will efficiently manage bond order information) and load the resulting model into the dictionary editor as a MOL or SMILES file. Alternatively, a set of template coordinates may already be available from the Protein Data Bank web sites. Useful resources for finding those small molecules that have already been found associated with protein structures are now available from the Protein Data Bank web sites.

The **Dictionary/Import Ligand/CIF**, **Mol**and **Pdb**options provide browsers for loading ligand coordinate information in mmCIF (__.cif_), MOL (.mol_, __.sdf_, .sd_) or PDB (**_.pdb_) formats.**

From the **Mol**option, the user must also provide a 3 character entity code after entering the ligand data file since, unlike the PDB and mmCIF formats, this information is not a part of that format. To avoid confusion and operational conflicts, it is best not to use three character codes that have already used by the Protein Data Bank for structures in the public domain. These entry codes may be searched from the PDB web site; currently about 5% of the available name space appears to have been taken. In particular, the CCP4/REFMAC5 software includes mmCIF dictionary entries for the majority of these ligands and it is convenient to take advantage of these files.

Both 2D and 3D mol file types are supported, depending on the 2D/3D parameter setting embedded in the mol file format.

The **Smiles**option leads to an interface that provides three different methods for importing the SMILES string.

The **File**option may be used with the associated **Browse**button to enter a file with extension **_.smi_ that contains a SMILES string. Alternatively, the**Smiles **option provides a text field in which to cut-and-paste a SMILES string.**

The **Database**option may be used if it is possible to provide a script that is capable of returning a SMILES string as standard output from a ligand identification number (reg number). To utilize this capability the script and the command to execute it should have been provided as the **File/Preferences/Environment/Smiles Database Command**parameter. This capability is mainly intended for users with access to a corporate small molecule data base since these data bases typically store SMILES data for chemical entries.

The ID Code (3 letters) field is required to provide the Dictionary Editor with a 3 character entity code for the ligand since these codes are not part of the SMILES description.

**6.2 The Dictionary Editor interface**

![http://mifit.googlecode.com/svn/wiki/images/image078.jpg](http://mifit.googlecode.com/svn/wiki/images/image078.jpg)

Figure 6.1 Image of a ligand loaded into the Dictionary Editor

The Dictionary Editor interface appears after completing the entry of ligand data. This application contains a window for viewing the molecule and displaying the current restraints. The left and right mouse buttons control the rotation and translation of the molecule; by holding Shift while pressing the left or right mouse button the size of the molecule may be adjusted. The keyboard short cut | controls the display of split screen stereo, as for the main canvas.

A right mouse button click may be used to provide access to a menu that allows changing the drawing style. Currently supported styles are **Line**, **Ball and Line**, **Stick**and **Ball and Stick**.

The **Show**check boxes on the left hand side of the interface provides options for displaying the various refinement restraints established for this molecule and for voiding the display of hydrogen atoms and atom labels. For example, checking the **Planes**box will display green meshes across those sets of atoms that a refinement program will restrain to lie in the same planes.

The **Refine**button optimizes the ligand atom positions towards consistency with the current set of restraints. This option is most often used when it is necessary to change a plane definition or a bond order.

The **Atoms**pull down contains options for removing individual selected atoms from the ligand. The **Atoms/Remove Hydrogens**option may be used to eliminate all hydrogen atoms from the entity. Hydrogen atoms should be stripped if you intend to use the main menu semi-automated ligand fitting option (**Refine/Find Ligand Fit and Conformer**).

The **Angles**, **Bonds**, **Chirals**and **Planes**pull down menus at the top of the interface window are access commands for changing the restraint definitions. These tools work in conjunction with interactive atom picking from the molecular image. For example, to change a bond order the **Bonds** checkbox should be selected and an atom bond should be picked. The resulting dialog box may be used to change the bond order.

The most common problem with automatically generated restraints is that the plane group definitions are only partially correct. The software attempts to deduce the correct subsets of atoms to combine into planes from the input coordinate data but this is a relatively difficult problem. In order to add or remove atoms from a plane first click on the plane identifier (for example, a 1 marking plane-1). The plane mesh will then change color to red, indicating that the selected plane is active. The relevant options under the **Planes**pull down menu are now active. Initially, only the **Remove Plane**option is available. After clicking on an atom in the model, options to add or remove atoms also become available and may be used to change the atomic composition of the plane.

It is important to note that the reliability of the automatically derived restraints depends on the quality of the input coordinates. Energy minimized ligand models will usually give correct restraints whereas inaccurate coordinate sets are likely to miss some restraints and may require adjustment.

For input from mmCIF files, which directly encode the required restraint information, the definitions of chiral centers, planar groups etc are read from the mmCIF file rather than generated by the MIFit code. This behavior is useful for checking the mmCIF dictionary files prepared for use with the CCP4/REFMAC5 program.

Although normally used to simple to provide coordinate data for single conformers, the capability generating and exposing multiple conformers for the ligand is available via the **Dictionary/Generate Conformers**command. The ligand conformations are generated in torsion angle space and are subsequently stored in the MIFit dictionary file as explicit coordinates, in a way analogous to the storage of amino acid side chains conformers. The conformers may be viewed using the **Conformer** counter.

**6.3 Exporting ligand information**

Once the ligand is parameterized by a correct set of refinement restraints the **Export to File**button may be used to write ligand dictionary information in the mmCIF format for refinement with CCP4/REFMAC5 or to write a coordinate file in PDB or MOL format. The **Save as type**selection in the file dialog is used to specify the output file format.

After selecting the **Save and Exit**button in the main Dictionary Editor window the ligand coordinates and associated restraints are added to the MIFit dictionary.

**6.4 Restraint dictionary files for SHELX and CNS/CNX**

The mmCIF (REFMAC) dictionary files output by the Dictionary Editor may also be converted to formats useful for refinement with either SHELX or the CNS/CNX family of programs. Although MIFit maintains its own dictionary format, the mmCIF format is used as the most standardized vehicle for information exchange. In addition, the mmCIF dictionary files available within the CCP4 installation provide a valuable resource and the conversion tools allow the utilization of these data with the other refinement programs. The twin commands **Job/Convert Cif to SHELX**and **Job/Convert Cif to CNS/CNX**may be used perform this format conversion.

Once the application for SHELX dictionary conversion has executed a dictionary file for the SHELX program will appear with default name _mi\_shelx.lib_in the same directory as the input mmCIF data file. If the command is subsequently repeated with a different mmCIF dictionary files then the additional restraints will be appended to the file _mi\_shelx.lib_.

The output files from this application that are intended for use with the CNS/CNX program have names based on the input ligand entity name. For example, if a ligand is described with entity code MMM in the input mmCIF file then the output is a parameter file, _MMM.par_, and a topology file, _MMM.top_. As with the files produced for SHELX, these files will appear in the same directory as the input mmCIF file.