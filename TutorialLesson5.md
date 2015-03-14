# Tutorial Lesson: Displaying molecular surfaces #

_Scenario: You have completed the refinement of a protein:ligand structure and you wish to examine the ligand binding volume and the protein surfaces._

1. Copy file 1A28.pdb from the examples directory of your MIFit installation into a new directory. This file contains a protein structure with a bound ligand (progesterone).

2. Use the MIFit command File/Open models, data, maps, etc... to load coordinate file 1A28.pdb from the directory in which you placed it. The structure will appear in the main canvas.

3. Use the mouse to select chain B in the Models List pane of the navigation tree. Scroll down to the end of the **Residues List** and perform a double mouse click on entity STR 2 in the Residues List section of the tree. These actions will center the canvas display on the bound ligand STR 2.

4. Use the mouse to click on any atom in the main canvas.

5. Make sure that the Show/Solid Surface/Molecular surface **mode** option is checked and select the Show/Solid Surface/Build Surface command. After a few seconds you should see the solvent exposed surface including a completely enclosed cavity that contained the STR molecule, _i.e._ in this protein the bound ligand is not accessible from the outside.

6. You may wish to click on the zoom and slab out icons or use the mouse controls to obtain a more global view of the protein surface.

7. You may select File/Close or File/Exit to close this session or shut down MIFit.