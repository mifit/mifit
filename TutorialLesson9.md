# Tutorial Lesson: Molecular replacement with CCP4/MOLREP or CCP4/PHASER #

(This tutorial lesson requires that MIExpert is installed and the system has access to the CCP4 software suite)

_Scenario: You have just obtained data on a new crystal form of a previously solved struc-ture. In contrast to the previously solved structure, which contained a monomer in the crystal asymmetric unit, the new crystal form may contain a dimer in the crystal asymmetric unit._

1. Copy the files 1A28\_mr.pdb and r1a28sf.mtz from the examples directory of your MIFit installation into a new directory. File 1A28\_mr.pdb is the initial search model (a monomer) that will be used for the structure determination. File r1a28sf.mtz con-tains structure factor data for the new crystal form.

2. Select Job/Molecular Replacement from the MIFit interface and the Molecular Replacement menu will appear.

3. Use the Browse... button associated with the Working directory parameter to find and load the directory defined in (1).

4. Use the Browse... button associated with the Model (PDB) parameter to load file _1A28\_mr.pdb_. Use the Browse... button associated with the **Data**parameter to load file r1a28sf.mtz.

5. Set the **Copies** parameter to 2 since crystal volume analysis suggest that there are two molecules in the crystal asymmetric unit.

6. The options to input a Fixed model, Search multiple models and Match input position should remain unchecked. The Fixed model parameter would be used if a molecular replacement search had been able to establish one component of the crystal structure (perhaps a single molecule from a complex or a domain of a large flexible structure). The Search multiple models parameter (a directory containing many MR search models) can be used in difficult molecular replacement cases to conveniently test a number of different models. The Match input position is useful when MR is applied to structures with same space group and cell as previously solved structures (_i.e._ co-crystals). It will adjust the MR solution by symmetry operations to match the initial search model.

7. The Spacegroup option may remain set to Default to use the space group specified by the diffraction data file header. The Method option may be set to MolRep to use the CCP4/MOLREP program to perform the molecular replacement as this is an easy case.

Select OK to launch the calculation. This process takes about 1 minute on a 2GHz laptop.

8. Select the Jobs tab at the bottom of the navigation tree. With a right mouse button click select the bottom entry in the Jobs List and select Show Log File. This action will provide a view of summary information from the molecular replacement process. Note that this command is grayed-out (inaccessible) until the job completes. Completion is indicated by the molecular replacement job icon changing from yellow to green.

9. Use the File/Open models, data, maps, etc... command to load the resulting model molrep\_1.pdb from the directory defined in (1). You will notice that two molecules were found by the molecular replacement process (indicated by A and B chains in the **Models List** section of the navigation tree) to fit the diffraction data and that these were automatically packed together.

10. Look at the files in the directory defined in (1). The project\_history.txt file logs the automation process. File molrep\_1.log is the log file from the molecular replacement run.

11. Although the default behavior of most molecular replacement programs is to exclude solutions in which molecules grossly overlap in the crystal it is always useful to check the packing of the solution.

12. Select Show/Symmetry Atoms/Show symmetry atoms as atoms to generate all symmetry related atoms in the volume displayed in the main MIFit canvas. You may wish to zoom out or move around the model and regenerate symmetry related positions to get a full sense of the crystal arrangement.

13. Since this appears to be a correct solution (R &lt; 0.35 with good packing) you may also wish to view this model in the context of an electron density map. A convenient way to do this is to perform a quick refinement of the current structure and then load the resulting refined likelihood weighted electron density map. The Job/Refinement command is described in the next tutorial lesson.

14. You may select File/Close or File/Exit to close this session or shut down MIFit.