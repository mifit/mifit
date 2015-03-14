# Model building and refitting #



MIFit contains a complete set of tools for building protein models into electron density maps. Capabilities include baton style Cα chain-tracing, conversion of Cα-traces to poly-alanine chain segments, the application of pentamer libraries to refit main chain atoms, convenient side chain mutation and density refitting, automated and manual water picking, ligand placement and adjustment.

**5.1 Chain tracing**

Most structure models are either built into experimentally phased electron density maps using automated procedures (for example, arp/wARP) or were derived from molecular replacement searches. In these scenarios the starting point for interactive model building is a model that is largely complete but where it is necessary to manually trace some portions of the map that were either not fit by these methods or were misfit.

The option to provide smart baton-style placement of Cα markers (residue type MRK) is available via the **Model/Add MRK Before**and the **Model/Add MRK After**commands. Convenient key-board short cuts equivalent to these commands are &lt; and &gt;. The function of these commands is to predict (based on density and geometric constraints) and add a Cα atom either before or after the current Cα atom. If the predicted position appears to be incorrect then the right mouse button may be used to adjust the new Cα site at fixed distance of 3.8 from the previous Cα.

The **Fit/Accept Fit**command may be used to accept the current set of positions and cancel the chain tracing.

Selection of the end points of the sequence of MRK points followed by selection of the **Model/Poly-ALA Range command**will convert the Cα trace to a polyalanine trace. Alternatively, if the entire chain consists of MRK points, the **Model/Poly-ALA Chain**command may be used to carry out the conversion to a full atom representation of the chain trace. The accuracy of the polyalanine trace depends on the accuracy of the Cα placement. The endpoints of the trace tend to be the least reliable because the atomic positions are less constrained than those in the central region. In the absence of a following residue, the positions of atoms in the C-terminal amino acid are not fully defined so this amino acid remains a MRK residue.

**5.2 Adding and removing individual residues**

![http://mifit.googlecode.com/svn/wiki/images/image072.jpg](http://mifit.googlecode.com/svn/wiki/images/image072.jpg)

Figure 5.1 Menu for adding a new residue to extend a chain or to add a ligand entity

The **Model/Add residue**command may be used to add new entities to the model. The pull down menu corresponding to the **Residue Type**parameter in the resulting dialog box is used to specify the type of entity that is to be added to the model. The **Insert Position**options describe where in the coordinate set the new entity will be placed. The **Chain id**parameter box may be filled if you wish to set a new chain id for the inserted residue. (This **Chain id** option is not normally used for protein model building but may be useful for adding ligands.)

The **Put At**options may be used to simply add the new residue at the screen center, extend an existing chain in **Alpha helix**or **Beta sheet**conformation or extend the existing chain with **Best Fit**of the new residue to the electron density.

The **Model/Delete residue**(keyboard shortcut **D**) may be used to remove a selected entity from the structure, _i.e._ after using this command the entity disappears from both the graphics display and the tree control and is no longer a part of the model.

The **Model/Rename residue**may be used to change the residue number for an individual entity.

The **Model/Delete atom**command (keyboard shortcut **Delete**) may be used to delete an individual atom from the structure _i.e._ the last atom selected from the canvas disappears from the canvas display and the tree control and is no longer part of the model. This command is mostly useful for pruning disordered side chains.

**5.3 Refitting backbone atoms**

The MIFit software contains commands to check and refit protein backbones based on comparison with a pentamer fragment library. A tool is also available for flipping individual amino acids so that the carbonyl oxygen is rotated about the peptide plane.

To employ the pentamer library the user must select an atom from the model displayed in the main canvas. After selecting the **Fit/Fix Backbone/Suggest Backbone Match**command the pentamer that best fits the current structure will be displayed in purple in the main canvas. The first Cα in the pentamer will belong to the amino acid associated with the selected atom.

Portions of the protein backbone that overlap the pentamer may be adjusted to improve the match using the **Fit/Fix Backbone**commands **Replace Middle 3**, **Replace First 4**, **Replace Last 4**or **Replace All 5**. If none of these commands appears to be useful then the **Fit/Fix Backbone/Clear Backbone Match**command may be used to eliminate the pentamer from the display and the tree control.

To flip an individual peptide the user should click on an atom in the peptide plane and then select the **Fit/Fix Backbone/Flip Peptide**command. The peptide will then be observed to rotate so that the carbonyl atom appears on the opposite side of the peptide plane.

**5.4 Refitting and changing side chains**

![http://mifit.googlecode.com/svn/wiki/images/image074.jpg](http://mifit.googlecode.com/svn/wiki/images/image074.jpg)

Figure 5.2 Menu for mutating residues

Side chain groups may be changed from one type to another by selecting an atom in the main canvas and then using the **Model/Replace residue**command. The action of this command is simply to mutate one amino acid type to another without reference to any electron density map. After selecting this command a dialog window appears from which the new amino acid type may be selected.

An extremely useful command for both model building and model correction (_i.e._ adjusting the position of a side chain to improve the fit to the electron density) is the **Model/Replace and Fit**command. This command also generates a dialog selection for mutating an amino acid, with default value set to the original amino acid type. The action of this command is to replace the amino acid side chain but also to make the best fit to the electron density. An alternative method for selecting this command is the keyboard shortcut r. The proposed fit (with the amino acid remaining colored green to indicate that it is still live) may be accepted or rejected using the green check or red cross toolbar icons.

It is also possible to cycle through the possible preferred conformers for a particular amino acid type by first selecting an atom in the main canvas and then using either the **Model/Next Conformer**command or the keyboard shortcut c.

A fitting method that allows arbitrary interactive adjustment of side chain torsion angles is to select the amino acid and then choose either the **Fit/Fit residue**command or the keyboard shortcut f. The numeric keyboard keys (1, 2 ) may then be used to select a torsion angle (χ1, χ2 ) for rotation. A grey arrow appears on the screen to indicate the selected torsion and the right-mouse button controls the torsion about the indicated bond. The green check and red cross toolbar icons and may be used to accept or reject the fit.

**5.5 General refitting**

Any entity (amino acid, water etc) in the current model may be activated for fitting by using either the **Fit/Fit residue**command or the keyboard shortcut f. Active entities are colored green in the main canvas. Entities may be translated by selecting the translate icon and adjusted using the right mouse button (the left mouse button retains control of the view orientation). The rotation icon may be selected to change the right mouse button to control the rotation of the active entity. Clicking on a bond in the active entity sets an arbitrary torsion with the rotating atoms at the arrow head marking the bond. The atom closest to selected position is marked with the arrow head which indicates the moving portion. The torsion icon may be selected to set the right mouse button control to the torsion.

Similarly, atomic groups other than individual amino acids may be activated for fitting. The **Fit**commands **Fit Atom**, **Fit Residues**, **Fit Residue Range**, **Fit Atoms**and **Fit Molecule**may be used in conjunction with the selection of atoms from the main canvas to activate portions of the structure for fitting. Note that for multiple selections the list in the bottom left hand corner of the main canvas records the entities currently in the stack. The **Show/Stack**commands may be used to manipulate the contents of the stack.

To accept a fit the **Fit/Accept Fit**command or the icon may be used. To cancel a fit the **Fit/Cancel Fit**command or the icon may be used. If you do not wish to deactivate the model but do wish to return to the original position of the active entity then the **Fit/Reset Fit**command may be selected.

**5.5 Modeling discrete disorder**

Discrete disorder (the presence of multiple well defined conformations, usually for side chains) is very simple to model with MIFit and is handled transparently by the refinement program CCP4/REFMAC5. Well refined structures at resolutions better than 2 will usually show some abnormal side chain densities that may be reasonably interpreted in terms of two major side chain conformations.

To set up an amino acid for modeling discrete disorder, select an atom within the amino acid in which the disorder needs to be modeled and activate it with **Fit/Fit residue**or the keyboard shortcut 'f'. The **Fit/Disorder/Split Fit**command may then be used to duplicate the entire amino acid, with the duplicate copy slightly offset from the original copy for clarity. The duplicated copy may be manipulated with the interactive fitting tools _i.e._ the translation, rotation and torsion adjustment options available from the toolbar.

More commonly, the required modeling requires splitting just a side chain or part of a side chain into multiple conformations. To do this, specify a side chain torsion angle using one of the methods previously described. The **Fit/Disorder/Split Fit**command may then be used to split the side chain at the specified torsion and the new copy may be rotated into position. As usual, the toolbar icons and may be used to accept or reject the fit.

**5.6 Local structure refinement**

The entire molecule or portions of the molecule may be optimized using commands beneath the **Refine**menu. The **Refine/Refine Options**command shows the weights used for the optimization. It will be noted that the refinement target includes a term for matching the electron density as well as tethers for the Cα positions at the ends of the refinement zone.

An active region for structure optimization is colored pale blue. Following optimization, the changed structure may be accepted with the **Refine/Accept Refine**command, canceled with the **Refine/Cancel Refine**command or allowed to remain active but reset to the unrefined position with the **Refine/Reset Refine**commands.

Once a region is active, holding down the keyboard space bar results in a continuous refinement of the selected region to convergence. This is often a powerful optimization technique, resulting in large shifts in the atomic coordinates to fit the electron density.

The entire molecule may be optimized using the **Refine/Refine Molecule**option. More usually, an amino acid is selected for optimization. The individual selected amino acid may be optimized using the **Refine/Refine Residue**command or the amino acid the amino acids on either side of it may be included in the optimization by using the **Refine/Refine Local Region**command.

Selecting two atoms by clicking on them in the main canvas and then applying the **Re-fine/Refine Range**command will result in all amino acids between the first and the last selected atom being included in the optimization.

**5.7 Adding water molecules**

The MIFit software contains several mechanisms for adding water molecules to an atomic model. One method is to fold water picking into automated refinement cycles through the **Job/Refinement**interface.

![http://mifit.googlecode.com/svn/wiki/images/image076.jpg](http://mifit.googlecode.com/svn/wiki/images/image076.jpg)

Figure 5.3 Menu for global addition of ordered water molecules to a model

MIFit also provides an independent method for water picking using the **Model/Add Waters...**command. This command is intended to pick waters from a map loaded into the main canvas. This will usually be a difference map in which densities for ordered water molecules are visible as peaks around the protein. The resulting dialog box from this command contains criteria for selecting waters within a border around the protein (via the **Minimum distance from protein**and **Maximum distance from protein**parameters). The **Minimum map level**parameter sets the threshold level for density to be considered as a possible water site. For a difference map this would ordinarily be set to a value of 3-4 sigma; for a conventional map a value of ~1 sigma might be appropriate. The **Maximum number of waters to add**parameter is useful to avoid overloading the map with water molecules. Since the optimal density map level parameter is somewhat problem dependent this is useful to prevent adding too many waters if that value is initially set too low. The **Asymmetric unit**parameter need not ordinarily be changed  it will normally correspond to a unique region of the crystal cell.

Upon selecting **OK**the water fitting process run and the new waters will be displayed in the main canvas. The number of waters added will be reported in the lower left of the status bar.

These methods serve to fit the majority of ordered water molecules in a protein structure. To delete an individual water molecule the most efficient method is to select it via a single mouse click in the main canvas and then using the **Delete**button on the keyboard.

Individual waters may either be added using the **Model/Add residue**command, setting the **Residue Type**pull down menu to HOH (the standard code for a water molecule). A potentially faster method is to use the keyboard shortcut **Shift+W**. This option caused the new water to appear at the mouse position in the main canvas and optimized to fit the density (_i.e._ not at the canvas center).

**5.8 Interactive ligand fitting**

After loading a protein model into the main canvas, the **Model/Add residue**command may be used to add any ligand in the MIFit dictionary to the model. (See the chapter on the Dictionary Editor for information on adding ligands to the dictionary).



The **Residue Type**parameter is a pull down menu that allows selection of any entity in the MIFit dictionary. The **Insert Position**parameter controls where in the coordinate file the ligand atom record will appear. The **Put At**parameter controls position of the ligand in the structure  to simply add a ligand to the model the ligand density should be moved to the screen center and the **Screen Center**option should be selected. The ligand will appear in arbitrary orientation on the ligand density.

MIFit does not require dictionary data to simply manipulate (rotate and translate) a ligand molecule or to rotate a part of the ligand about an arbitrary torsion angle. In a scenario where the atomic model (including the ligand) and electron density map are loaded in the main canvas and you wish to fit the ligand to a density feature: (i) select the ligand for fitting with a single mouse click on any atom in the ligand (ii) activate the ligand with a mouse click on the **Fit/Fit Residue**menu item or by using the keyboard shortcut f. These actions will highlight the ligand molecule in green. The left mouse button now controls rotation of the entire canvas view and the right mouse button may be set to rotate or translate the ligand molecule. As already described for fitting amino acids, the two icons to the right of the Cancel Fit icon are used to toggle between translational and rotational fitting modes.

A mouse click on a bond in an active residue (ligand) establishes an arbitrary torsion angle within the ligand. This torsion angle is highlighted as a grey arrow along the rotatable bond. The atomic group closest to the point that was clicked will be the part of the ligand that will move. If the icon to the right of the icons that control the translation and rotation operations is selected, a portion of the ligand may be rotated about this bond using the right mouse button. It will usually be fairly obvious to a human with minimal chemistry training which torsion angles within a ligand are allowed relatively free rotation.



Taken together, these options allow the crystallographer to not only quickly orient a ligand molecule in density but also to fit groups of atoms within the ligand by twisting them about rotatable bonds.

**5.9 Semi-automated ligand fitting**

MIFit contains a command (**Refine/Find Ligand Fit and Conformer**) for docking ligands to difference density during model building sessions. This command performs rotation and translation searches using a genetic algorithm (D. E. McRee, _Acta Cryst_. **D\*60, 2276-2279, 2004) on ligands that have been arbitrarily overlapped on to the ligand difference density, including selection of the best fitting conformer. The ligand must be present in the MIFit dictionary to use this option.**

The ligand fitting method works best on difference maps that show strong, well defined ligand densities. The interface to this application currently limits the map data to better than 2.5 resolution in order for the method to be applicable. Trials suggest that this is the resolution at which the shape of the density will usually be sufficiently distinctive to correctly fit a ligand. The fitting process will usually work best if stray electron density features (usually corresponding to ordered water molecules near the ligand) are already modeled so as to remove them from the difference map. Large and flexible ligands tend to be the more difficult to fit than small rigid examples; molecules that contain more than 60 non-hydrogen atoms will probably be too large to fit and are currently disallowed by the interface. On the other hand, compact molecules (say consisting of 3-4 planar groups separated by torsion angles) are often fit quickly and correctly.

When the ligand density has been identified in a difference density map, the ligand molecule is first added to the protein coordinate data by selecting the required entity using the **Model/Add residue**command. The **Residue Type**pull down may be used to select the ligand from the list of entities in the MIFit dictionary. The **Insert Position**parameter should usually be **End of Model** and the **Put At**parameter will be set to **Screen Center**. The **Chain id**field may be used to override the default chain selection by entering another single character chain id  usually when a user wishes to enter a new, unique chain-id for the ligand. Clicking on **OK**will place the ligand in an arbitrary orientation at the screen center, which should have been adjusted to overlap the ligand difference density.

The ligand should be selected as 'active' using the **Fit/Fit Residue**command (or keyboard shortcut 'f') and the fitting process launched using the **Refine/Find Ligand Fit and Conformer**command. The process will terminate when either the maximum allowed number of generations is reached in the genetic algorithm. Scores and progress for a running process are displayed in the log window pane and the right segment of the status bar. It should be noted that the position of the ligand that is displayed in the graphics during a ligand docking process is a current position and is not the best fit to that point. The best fit is loaded when the process completes. This procedure typically takes ~1 minute and is directly proportional to the number of conformations sampled (often 200-300).

In some cases the result is a ligand fit that is basically correct but is only an approximate fit to the electron density. When this occurs a useful possibility is to use the **Refine/Refine Residue**option to apply real space refinement to optimize the fit of the ligand to the map.
