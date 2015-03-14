# Tutorial Lesson: Chain tracing #


_Scenario: You have just obtained an electron density map and wish to build a model into the electron density. (Note: automated model-fitting programs are now usually used to build initial models into experimentally phased map. In other cases, molecular replacement methods are used to obtain an approximate model. In both these situations the manual building process is often limited to tracing loops and model correction.)_

1. Copy files start\_trace.pdb and automation\_out.mtz from the examples directory of your MIFit installation into a new directory. The file start\_trace.pdb just contains two Cα marker atoms (MRK residues) and is used to initiate the model-building process and guide the tutorial. The file automation\_out.mtz contains phased structure factor data that will be used to compute an electron density map into which the model will be built.

2. Use the command File/Open models, data, maps, etc... to load the coordinate file _start\_trace.pdb_and data file automation\_out.mtz from the directory in which you placed it. In the phase file import dialog, deselect loading the difference map (map coefficients DELFWT and PHDELFWT) and click on **OK**. The two marker structure and map should appear in the main canvas.

3. Make sure that the region of the map has appropriate volume and contour levels for model building. With the right mouse button click on the map icon in the Models List section of the navigation tree and select the Contour options... command. In the resulting dialog box, move the Radius slider to about 10 in order to set a 10 box for the map volume. You may also wish to set the Preset Map Styles option to Blue Map, 1,2,3,4,5 sigma and then uncheck the display of contour levels 2 and 4 to remove some of the contour lines from the map display. Select OK and you will see the map view change to these settings.

4. In the Go to residue parameter box at the top of the navigation tree enter X 150. This action will center the display on Cα marker X 150.

5. Perform a mouse click on the (centered) MRK 150 entity in the main canvas and then select the Model/Add MRK After command. You should see a green line appearing to connect that marker to the expected position for the next Cα atom.

6. By continuing to select Model/Add MRK After or by using the keyboard short cut SHIFT &gt; trace out four more Cα positions.

7. Select Refine/Accept Refine and Fit/Accept Fit to accept this trace.

8. Select Model/Poly-Ala Chain to convert this Cα chain trace to a polyalanine chain trace.

9. The C-terminal residue is not fully converted to ALA (it remains as a MRK residue) and may be deleted by making a mouse button click on the MRK atom in the main canvas. Then, make a right mouse button click in the main canvas and select Delete Residue from the pop-up menu.

10. Since this is a high resolution map showing atomic detail it is useful to optimize the poly-ALA trace before adding side chains. Perform a mouse click on an atom at both the N- and C-terminal ends the chain. Select Refine/Refine Range to activate this chain for refinement (it will change color to blue) and then hold down the space bar to drive the refinement to convergence. Accept the result by selecting Refine/Accept Refine.

11. As an exercise you could now add side chains by selecting a residue with a mouse button click in the main canvas and then using keyboard short cut 'r' to open the 'replace and fit' dialog box. After selecting the correct residue type from the residue pull down menu it is replaced with optimal matching to the electron density map. (The true sequence for this polypeptide trace is TVITKHN)

12. You may select File/Close or File/Exit to close this session or shut down MIFit.