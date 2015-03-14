# Tutorial Lesson: Loading and changing electron density map displays #

_This lesson introduces some of the mechanisms for loading and adjusting the display of electron density maps._

1. Copy files automation\_out.pdb, bng\_mlmap.fcf and bng\_milaunch.mlw from the examples directory of your MIFit installation into a new directory. These files are a co-ordinate file, a file of phased x-ray data and a MIFit session file for loading the coordinates and data with preset view information.

2. Use the File/Open models, data, maps, etc... command to load the file bng\_milaunch.mlw from the directory (1) in which it was placed. You should see the electron density map appear around a protein model in the main canvas. Note that space group and cell information (the text beside the crystal icons) are independently read from both the input coordinate file and the x-ray data file. In the case of an _.fcf_ format file the space group information is absent but the crystal symmetry operators are encoded.

3. Try using the right mouse button to pan across the map. You should see that the map is automatically re-contoured as you scroll across the main canvas and you move to the edge of the current box of map density.

4. The map size and contour levels can be changed using the Contour Options command in which may be accessed via a right mouse click on the map icon in the **Models List** tree. Move the Radius scrollbar to change the size of the display density to a 9 box. From the pull down menu in the Preset Map Styles list box make sure that the first option - Blue Map 1,2,3,4 5 sigma is selected. The lowest contour level for displaying this map is set to one sigma (Crystallographers refer to the root-mean-square density fluctuation of the map in terms of sigma, although this may be confusing as it does not relate to the error in the electron density for this type of map). MIFit internally scales maps so that one sigma is set to 50 units. The other density levels are initially set to 2 sigma, 3 sigma etc for up to 5 contour levels.

5. Try turning off contour levels 2 and 4 by removing the checks from the associated Show parameters. Click on OK and you should see these two intermediate contour levels disappear from the display, leaving a less cluttered image.

6. Try changing the resolution of the map by selecting the FFT Phases command via a right mouse click on the map icon in the **Models List** tree. Set the Min Resolution parameter to 3.0 and select OK. The operating default for the FFT had been set to create the map using the full resolution limits of the data but you should now see the map at 3 resolution showing much less detail.

7. Try changing the grid spacing for the map contours by selecting the FFT Phases option again and changing the Grid option from the default (Medium Grid) to Fine Grid. Click on OK and you should see that the map is returned with a much smoother appearance. The Medium Grid setting is used for most model building activities since extra contour lines impede graphics performance and may make the model hard to see but presentation images sometimes look better when over-contoured on a fine grid.

8. Zoom out the model and select the Show/Symmetry Atoms/Show symmetry atoms as atoms command. You should see some symmetry related atoms (purple) appear around the molecule.

9. You may select File/Close or File/Exit to close this session or shut down MIFit.