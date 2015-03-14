# Tutorial Lesson: Refitting individual residues #


_Scenario: You have been refining a protein structure which is correct except for the positions of a few residues which still misfit the electron density map after the last cycles of refinement. You wish to refit these amino acids. The lesson illustrates methods for refitting individual amino acids._

1. Copy files automation\_misfit.pdb and automation\_out.mtz from the examples directory of your MIFit installation into a new directory. These files are a coordinate file that contains a few model fitting errors and file of phased x-ray data in mtz format from a refinement with CCP4/REFMAC5.

2. Select the MIFit File/Open models, data, maps, etc command to load file automa-tion\_misfit.pdb and automation\_out.mtz from the directory established in (1). (Hold down the Shift keyboard button as you use the mouse to select both files.)

3. A dialog box will appear that allow selection of map types and MTZ columns labels. In this example we will create the electron density map from pre-computed likelihood-weighted structure factor coefficients and phases that were calculated by the REFMAC5 program. _i.e._ we will simplify the calculation route by utilizing pre-computed map coefficients with the **Direct FFT** option. Since MIFit automatically recognizes MTZ files with standard REFMAC5 column labels, the dialog will automatically be set to load these values. Simply click on OK.

4. Scroll down to residue X 142 in the Residues List section of the navigation tree and make a double mouse click on it. This action will re-center the model and map on the Cα atom in residue X 142 Tyr. You will notice that the side chain of this amino acid does not fit the density.

5. Select the residue to be fit (X 142 Tyr) by clicking on any atom in it in the main canvas. Now select either the Fit/Fit Residue command or the icon in the toolbar or the keyboard shortcut f. The residue should turn green to indicate that it is active for fitting.

6. Try translating the residue by either clicking on the Fit/Translate command or the translate icon in the toolbar. Using the right mouse button you can now move the selected residue around the canvas.

7. Try rotating the residue by either clicking on the Fit/Rotate command or the rotate icon in the toolbar. Using the right mouse button you can now rotate the selected residue around the canvas.

8. Reset the position of the amino acid back to the starting position but remain in fitting mode by selecting the Fit/Reset Fit command.

9. To activate rotation of the side chain about the χ1 angle (_1_ along the bond between the Cα and C atoms) enter 1 on the keyboard. You should see a grey arrow that points from the Cα to the C atom appear on the model in the main canvas window. Use the right mouse button to rotate the side chain about χ1 into the electron density.

10. The side chain also needs to be adjusted slightly about the χ2 angle (_i.e._ the bond between the C and Cγ atoms). Enter 2 on the keyboard and you should see the grey arrow shift

to point from C to Cγ. Use the right mouse button to rotate the side chain about χ2 into the density.

11. In a real fitting session you would now accept this fit either by selecting the Fit/Accept Fit command, by clicking on the green check icon on the toolbar or by using the keyboard shortcut ;. Here, select Fit/Reset Fit so that the lesson can continue to show other fitting methods. Select the cancel fit icon (a red cross) to deactivate this residue.

12. A convenient way to refit side chains is to use the Model/Replace and fit command (or the keyboard shortcut r). Try this by making a mouse click on any atom in X 142 Tyr in the main canvas and entering 'r' from the keyboard. A dialog box will pop up with a pull down menu containing all entity types in the MIFit dictionary. Since we want to retain the default amino acid type (Tyr) select OK. You will see that the Tyr side chain move to fit the electron density. To save this fit you could select the green check icon from the toolbar. Alternatively, to cancel this fit click, select the cancel fit icon from the toolbar.

13. We will now correct a common type of main chain fitting error, a peptide flip, in which the carbonyl oxygen is rotated by approximately 180 degrees from the correct orientation. Make sure that Chain:X in the **Models List** part of the navigation tree is still selected (icon is colored blue), scroll up the Residues **List** to residue X 37 and make a double mouse click on it. This action will re-center the model and map on the Cα atom in X 37 Phe. You will notice that the main chain of this amino acid does not fit the density.

14. In the main canvas, click on the carbonyl (O atom) in X 37 Phe with the mouse. Select the Fit/Fix Backbone/Flip Peptide command. You will see the peptide plane flip over so that the carbonyl now properly fits the electron density.

15. You will now use a pentamer library to check a part of the protein structure. Make sure that Chain:X in the **Models List** part of the navigation tree is still selected (icon is colored blue), scroll down to residue X 192 and double click on it. This action will re-center the model and map on the Cα atom in X 192 Glu.

16. Click on any atom in X 192 Glu in the main canvas. Select Fit/Fix Backbone/Suggest Backbone Match. You will see a pentamer match to the backbone, starting at X 192 Glu, appear on the screen in purple. In this case the protein conformation is similar to the pentamer and there is no reason to change it. If the atomic model dims you may turn of the dimming by unselecting the Render/Dim Non-active models command. The image of the pentamer match can be removed with the Fit/Fix Backbone/Clear Backbone Match command. In cases where the backbone should be refit to match the pentamer the four Fit/Fix Backbone/Replace commands could be used to replace parts of the structure.

17. You may select File/Close or File/Exit to close this session or shut down MIFit.