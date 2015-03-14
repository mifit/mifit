# Tutorial Lesson: Superimposing ligands from multiple related structures #

(This tutorial lesson requires that MIExpert is installed and the system has access to the CCP4 software suite)

_Scenario: You have solved the structures of six HIV protease:ligand complexes and wish to display and make an accurate comparison of the binding positions of the ligands in each of them._

1. Copy files 1d4h.pdb, 1d4i.pdb, 1d4y.pdb, 1eby.pdb, 1ebw.pdb, 1d4j.pdb from the examples directory of your MIFit installation into a new directory. Each of these files contains the HIV protein with a different ligand. This step is required because the superposition application will superimpose all PDB files from within a specified directory.

2. Create a working directory (different from the directory established in (1)) into which files containing the superimposed ligands will be placed.

3. Use the MIFit command File/Open models, data, maps, etc... to load coordinate file 1d4h.pdb from the directory in which you placed it. The structure will appear in the main canvas.

4. In the **Models List**section of the navigation tree select segment find and select entity BEH 501. This action is performed so as to center the display on the bound ligand BEH in the target structure but does not affect the superposition calculation.

5. Select the **Job/Cocrystal Superposition** menu.

6. Set Working Directory to the name of the directory specified in (2) above. Set Structure Directory to the name of the directory specified in (1), in which you stored the set of six coordinate files. Set Target File to structure _1d4h.pdb_.

7. Set the X, Y and Z parameters in the **Target Coordinates**field to a point 14, 23, 5. These are coordinates in angstroms that are close to atom C01, a central atom in the ligand that is used to define the binding site volume for the structure superposition. Select **OK**.

8. After a few seconds a coordinate file, _allligands.pdb_, containing all superimposed ligands will appear in the working directory. These ligands were identified and superimposed using the ligand binding site Cα atoms from the target structure. You may load this file using the **File/Open models, data, maps, etc...**command.

9. The superimposed ligands in this file are identified by separate chains. If you wish, you can select individual ligand _Segmt_ entries from the **Models List** in the navigation tree with a right mouse button click and use the **Color** command to highlight individual ligands.

10. You may select File/Close or File/Exit to close this session or shut down MIFit.