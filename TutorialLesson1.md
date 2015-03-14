# Tutorial Lesson:  Basic manipulation of the model display #

_This lesson introduces some of the mechanisms for controlling the view and representation of an atomic model._

1. Copy the file automation\_out.pdb from the examples directory of your MIFit installation into a new directory.

2. Select the File/Open models, data, maps, etc... command to load a coordinate file in PDB format. Use the Look in pull down menu in the resulting browser window to locate the file in the directory used for (1). Load the file (by a mouse clicking on it) and select Open. You should see the protein model appear in the main canvas window.

3. Try rotating the model by dragging the left mouse button in the canvas window.

4. Try zooming in and out by dragging the left mouse button while holding down the Shift key. Dragging the left mouse button up zooms in and dragging the mouse button down zooms out.

5. Try panning the model by dragging with right mouse button in the canvas window.

6. Notice that if the cursor is positioned in the extreme top of the main canvas the function of the left mouse button changes so that dragging it left or right rotates the view about the Z-axis. The function of the right mouse button also changes so that it controls the slab.

7. Try mouse clicks on the two slab icons from the toolbar to slab the canvas view in and out.

8. Try to re-center the display by a double mouse click on an atom in the main canvas.

9. Try centering on a residue by a double mouse click on an amino acid in the Residues List section of the navigation tree on the left.

10. Try centering on a residue by entering a chain-id, residue name pair (say, X 35) in the Go to residue parameter box at the top of the navigation tree on the left and hitting the Enter key.

11. Try centering on a residue by a double mouse click on an amino acid in the sequence window on the right.

12. Try a single mouse click on an atom in the main canvas. Notice that this atom is placed in the stack in the lower left of the canvas. The stack is used to identify targets for commands that operate on specific atoms or residues.

13. Try coloring a residue by clicking on a residue, clicking on the paint can icon in the toolbar and selecting the Color Last Picked Residue option. In the resulting dialog box select a color by clicking on it and then use the Method pull down menu to select to All atoms. When you select OK you should see the color of all atoms in the selected residue change. Try some of the other color options. You can also perform most of the same commands from the Render/Color options.

14. Try changing the rendering style by selecting the different options from the Render pull down menu. Many people prefer to use the Sticks option for speed and simplicity when model fitting.

15. Try the command Viewpoint/Center model on screen to zoom and slab the model automatically to fill the screen.

16. Try changing the model to a backbone only representation by selecting the Show/Backbone/Show backbone as CA trace and the Show/Sidechains/Hide sidechain atoms commands. You should now see only connected CA atoms of each residue (and any ligands, waters or prosthetic groups, if they were present in the input model). It is simpler to see the overall fold of the protein in this representation

17. Try the command Show/Secondary Structure/Make Ribbon to make a ribbon representation. A popup dialog box asks whether you wish to retain the original representation of the model (the CA trace) in addition to building the ribbon diagram. You may select Yes to remove the CA trace. The resulting model will be represented in different colors for helix, sheet and coil portions of the molecule.

18. Try changing the ribbon colors by selecting the Show/Secondary Structure/Ribbon Colors command. Clicking on each of the Helix Color..., Sheet Color... and Random Coil Color... buttons will provide a color palette for selecting new colors for each of these secondary structure types. Select red for the helix, yellow for the sheet and green for random coil. Click OK to see the model re-colored with these choices.

19. Select Show/Secondary Structure/Clear Ribbon to remove the ribbon display. The canvas may now appear empty as the model is hidden.

20. Select Show/Secondary Structure/Show Tube Secondary Structure to draw the molecule as a glossy tube. You may answer either Yes or No to the dialog box that asks Do you want to hide residue atoms. Select Show/Hide Secondary Structure to remove the tube representation.

21. Select Show/Secondary Structure/Show Schematic Secondary Structure to draw the molecule with schematic secondary structure rendering. You may answer either Yes or No to the dialog box that asks Do you want to hide residue atoms.

22. You may select File/Close or File/Exit to close this session or shut down MIFit.