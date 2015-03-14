# Tutorial Lesson: Establishing restraints for structure refinement with CCP4/REFMAC5 #

_Scenario: You have completed the preliminary refinement of a new protein structure and now wish to create a dictionary for a small molecule ligand molecule within the structure for use with CCP4/REFMAC5. You have used a chemical drawing program to obtain a SMILES string describing the molecular structure of the ligand._

1. Copy the file ligand\_example.smi from the examples directory of your MIFit installation into a new directory. This file contains a SMILES string describing the chemical structure of the ligand. Note that the Dictionary Editor also supports entry of ligand coordinate information in the form of MOL, PDB and mmCIF formats.

2. Select the command Dictionary/Import Ligand/Smiles. Make sure that the File option is toggled on and use the Browse... button to enter the name of the file containing the SMILES string. Enter a three character entity code (for example, LIG) for the ligand in the ID Code field. Select Ok.

3. You may use mouse button clicks on the Chirals and Planes checkboxes to confirm that the interpretation of the ligand in terms of restraints was correct. You will find that there will be no chiral center restraints and one plane covering the aromatic ring and atoms immediately connected to it. Incorrect restraint assignments may be changed using the pull down menus at the top of the Dictionary Editor window.

4. To write out the refinement dictionary, select the Export to File... button. Enter a file name, select Refmac mmCIF Dictionary File (.cif) in the Save as type field, and click on OK. Many crystallographers prefer to ignore the relatively weak and often poorly defined torsion restraints from ligand refinement and select No in the resulting Write Torsions? dialog box.

5. Select the Save and Exit button to exit from the Ligand Dictionary Editor.

6. You may select File/Exit if you wish to shut down MIFit. When you exit from MIFit you will be asked if you wish to write the new ligand entry into the MIFit dictionary file so as to retain it for future use.