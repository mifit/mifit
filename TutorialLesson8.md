# Tutorial Lesson: Automated ligand density docking #

_Scenario: You have an almost completely refined structure and now wish to fit a ligand to residual difference density within the target site._

1. Copy the files ligand\_example.pdb, ligand\_example.mtz and ligand\_example.smi from the examples directory of your MIFit installation into a new directory. File ligand.pdb is the parent protein model, file ligand\_example.mtz contains data from the last cycles of refinement in MTZ format and file ligand\_example.smi contains a SMILES string corresponding to the ligand that is to be fit.

2. To establish the ligand dictionary parameters select the Dictionary/Import Ligand/Smiles command. With the File option toggled on, use the Browse... button to select file _li-gand\_example.smi_from the dictionary defined in (1). Enter a three character code (for example, LIG) as the ID Code (3 letters) parameter and select **OK**. You should see the Dictionary Editor interface appear.

3. See the previous lesson if you wish to confirm that the ligand restraints were correctly defined. Click on Save and Exit to load the ligand into the dictionary.

4. Next, load the protein model and map using the File/Open models, data, maps, etc... dialog box to locate file ligand\_example.pdb and ligand\_example.mtz from (1). This lesson only re-quires the difference map using the DELFWT and PHDELFWT columns so for clarity you may wish to deselect the calculation of the normal density map. If multiple maps are loaded then you should click on the required map in the **Models List** section of the navigation tree to make it the fitting target.

5. Use the left (rotate) and right (translate) mouse buttons to make sure that the large density feature near residue A 163 is centered on the cross hairs.

6. To add the ligand to the map select the Model/Add residue command. Use the Residue Type pull down to select the LIG entity. Select End of Model for the Insert Position parameter. Make sure that Screen Center is toggled on for the Put At parameter. You may enter a new chain-id W in the **Chain Id** field. Select OK to add the ligand to the model. It should appear in the center of the canvas.

7. Click on any atom in the ligand and select Fit/Fit Residue. Alternatively you may use the keyboard shortcut f. You should see the ligand turn green, indicating that it is active.

8. Select the Refine/Find Ligand Fit and Conformer command to automatically fit the ligand to the difference density. The ligand docking search will continue until the maximum number of trials is reached or a sufficiently good fit is established. At the completion of this process a good fit of the ligand will usually be obtained. Since this search starts from a randomized starting point, fully reproducible results are not obtained on every run and occasionally the search will result in a plausible misfit.

9. You may sometimes be able to optimize the fit by selecting the Refine/Rigid-Body Refine Current Atoms option. If just one or two torsion groups are incorrect then they may be fixed by clicking on the both defining the torsion group near the move-able atoms and rotating around that bond with the right mouse button.

10. Once you have achieved a satisfactory fit you may accept the fit by clicking on the icon from the toolbar.

11. Ligand fits to density may also be optimized using real-space refinement. Click an atom on the ligand and then select Refine/Refine Residue. The ligand color will change to blue to indicate that it is active for refinement. Hold down the space bar to refine the ligand to convergence. Select Refine/Accept Refine to accept the optimized solution.

12. You may select File/Close or File/Exit to close this session or shut down MIFit.