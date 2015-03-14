# Displaying structure data #

**7.1 Setting model display modes**

MIFit provides many options for displaying protein structures. These options may be selected from the **Render**menu. The basic choices for displaying the atomic model in the main canvas are listed at the top of the menu as:

> Sticks

> Knob and Stick

> Ball and Cylinder

> CPK

The **Render/Sticks**option represents the atomic model by drawing lines between the atomic centers. The **Render/Knob and Stick**option adds to the stick representation by marking the atomic positions with small spheres. The **Render/Ball and Cylinder**option is a version of the **Render/Knob and Stick**option that presents a more three-dimensional appearance. The **Render/CPK**option provides a fully space-filling representation of the molecule.

For the **Render/Ball and Cylinder**display mode the user is able to specify the diameter of the balls and cylinders representing the atoms in the atomic model. The **Render/Set Ball/Cylinder Size**command provides access to a dialog window in which the **Ball Diameter**field allows specification of the ball size as a percentage of the CPK diameter and the **Cylinder Diameter**field allows control of the cylinder size as a percentage diameter of the ball diameter. Setting the **Cylinder Diameter**parameter to ~90% of the **Ball Diameter**parameter provides a way to make thick sticks with a full three-dimensional appearance.

The **Render/Smooth Lines**option is usually selected as this improves the appearance of lines by anti-aliasing. However, this option also slows down the interactive aspects of the model display and, depending on the speed of the host computer, it may be necessary to set **Render/No Line Smoothing**for general model fitting work.

The **Render/Depthcue Lines**and **Render/Depthcue Colors**options both provide a sense of three-dimensionality in the z-direction in the absence of stereo display and these options are normally toggled on. Sometimes it may be more useful to maintain a brighter image by turning off the depth cueing of colors. Depth cueing of only lines still provides some sense of a three dimensional image.

The **Render/Line Thickness**option is used to change the line thickness (in pixels) for the model display. This option is usually set so that lines have a thickness of one pixel but it may occasionally be useful to brighten the image by making the lines slightly thicker.

All model display styles allow the concurrent display of electron density contours, although the electron density surfaces may not be very clearly visible when the **Render/CPK**option is applied. However, the CPK representation combined with the display of difference density contours is a good way to visualize both the crystallographic aspects of ligand binding sites and surface cavities.

All of the model display options are available stereo displays as well as for mono viewing. If several model documents are open, each open model document may be set to different model display parameters and/or stereo mode.

When displaying multiple models during model-fitting operation it may be helpful for the non-active (reference) models to be less displayed with less intensity than the working model. Conversely, if displaying two models in order to create a figure to illustrate their similarities, the models would be displayed with equal visibility. For this reason a **Dim Non-active models**option is available and the extent of the diming may be controlled with the **Set Amount to Dim Non-active Models**parameter.

**7.2 Displaying model surfaces**

A different type of model display involves the representation of molecular surfaces. These representations are somewhat akin to the display of electron maps in that they may appear concurrently with a normal model display option. However, unlike electron density maps, which require structure factor data for their calculation, surface calculations do not require any additional information beyond the atomic model.



**Dot surfaces**

Options related to dot surface displays are selected from either the Show/Dot Surface menu or from inside the canvas window by performing a right mouse button click to pop up the Quick Menu. Selecting the Show/Dot Surface/van der Waal Surface option creates a dialog window that contains options for controlling the calculation of the dot surface. Only one dot surface can exist in MIFit at any given time. The maximum number of dots/Angstrom is 10 and for a high dot density the display somewhat resembles a CPK display. When the surface extends over large parts of a protein structure the response of the interactive display may be slowed down.

The majority of the Show/Dot Surface options result in the same type of dialog boxes except that in some cases the bounding box is an unnecessary parameter. A particularly useful surface is the solvent exposed surface that may be generated using the Show/Dot Surface/Solvent Exposed Surface command. This is a relatively smooth surface that provides a means of visualizing niches and cavities in the protein that are capable of binding small molecules.

The Show/Dot Surface/Clear Surface command may be used to eliminate an existing dot surface.

**Solid surfaces**

** **

For some purposes it is an advantage to create opaque, solid surfaces rather than semi-transparent dot surfaces. The solid surfaces hide the clutter of the atomic model representation.

Commands related to displaying solid surfaces may be accessed from the Show/Solid Surface menu. The top set of commands Build Surface, Color Surface, Color Surface by atom type and Clear Surface are 'greyed out' (inaccessible) until an atom is selected.

Build Surface renders the surface in either Molecular Surface Mode or Accessible Surface Mode. The surface produced by Molecular Surface Mode is similar in extent to a Van der Waals surface while the surface produced in the Accessible Surface Mode is expanded by the radius of a nominal water probe molecule.

Surfaces may be built based on an Atom Selection, a Single Residue Selection, a Multiple Residue Selection, Peptide Selection or a Molecule Selection. The Peptide Selection option is the most commonly used command since this will create a surface based on all the protein (amino acid) component of a model, leaving cavities containing ligands unfilled. Note that it is possible to build multiple surfaces so, for example, it is possible to render surfaces for two different models showing component of a molecular assembly in the same document.

Solid surfaces are often used to provide a quick visual assessment of whether a designed ligand molecule is able to bind to a particular site without clashing with the protein. Surface depictions are also frequently seen in journal articles, providing an overview of ligand binding cavities and protein surface topography.

Software requirements of fast interactive rendering versus creating an accurate and attractive appearance are somewhat contradictory so the user has two types of option to fine tune the appearance of the surface display. The Set Smooth Level command allows the user to provide increasing numbers of surface smoothing iterations (on a scale of 0-5). Note that the smoothing changes the exact dimensions of the surface. The Standard Quality, High Quality and Ultra High Quality settings change the size of the grid on which the surface is initially computed. These settings improve the accuracy of the surface at cost of some computer time in the initial surface computation and slowing down graphics speed for interactive displays.

**7.3 Backbone and ribbon displays**

Displaying all of the atoms in a model shows the greatest amount of detail but may also make it difficult to visualize the overall features of a structure. An understanding of the fold and major topological features of the entire structure may be facilitated by displaying just the Cα carbon trace or the chain trace as a ribbon diagram, a display mode popularized by Jane Richardson in the 1980s.

The **Show/Backbone**command may be set to **Show backbone as atoms**(_i.e._ the default all atom representation), **Show backbone as CA trace**or **Hide Backbone**in order to control the representation of the backbone. Similarly, the **Show/Sidechains**command contains options **Show sidechain atoms**and **Hide sidechain atoms**to manage the display of all side chains. In order to create an alpha-carbon trace you can select **Show/Backbone/Show backbone as CA trace**and then select **Show/Backbone/Hide sidechain atoms**to hide the side chains.

Similarly, the **Show/Secondary Structure/Make Ribbon**command creates a ribbon image of the current model. The **Show/Secondary Structure/Clear Ribbon**command removes the ribbon image and the **Show/Secondary Structure/Ribbon Colors**command spawns a dialog window that allows the user to change the default colors for the helix, sheet and random coil secondary structure elements.

MIFit also contain options for creating worm and schematic cartoon displays of the protein chain trace in user-managed colors with the **Show/Secondary Structure/Show Tube Secondary Structure**and the **Show/Secondary Structure/Show Schematic Secondary Structure**commands. The **Show/Secondary Structure/Hide Secondary Structure**command may be used to remove these displays.

**7.4 Choosing a background color**

A black background is the standard (default) choice for interactive work with atomic models. However, if you wish to print the display canvas on overhead film or capture the image for a PowerPoint presentation, black is usually a poor choice. An image that appears attractive and displays with good contrast on a computer screen will often appear dark and hard to see when viewed as an overhead or on a PowerPoint slide. Although these issues are to some extent a matter of personal preference, light (often white) backgrounds often look best for these types of presentation.

Printing on paper those images that contain black backgrounds also has several negative consequences. One practical issue, especially when using plain paper in ink-jet printers, is that what emerges from the printer is a soggy welled piece of paper as this soaks up a lot of black ink. The paper takes quite a while to dry and become useable. Another side effect is on the lifetime of the ink cartridge. Manufacturer-specified ink cartridge life times are based on the assumption that roughly 5% of the paper is actually covered with ink, which is a reasonable estimate when printing text.

The background color in the main MIFit canvas may be changed using the **Render/Set Background Color**command.

This command provides a color palette from which a color may be selected by clicking on it. After selecting **OK**the background color for the main MIFit canvas will change to the selected color.

**7.5 Capturing images from the MIFit main canvas**

MIFit provides the **File/Export Image As...**option to write the canvas image as a file in one of a variety of standard image formats _(_.jpeg_, .png_ or **_.tiff_).**

Alternatively, it is possible to capture the screen view directly for pasting into an open document (perhaps a Word or PowerPoint file) using the **File/Copy Canvas**command.

While displaying a structure with MIFit on the computer screen one might simply want to produce a quick hard copy of the MIFit canvas contents and the **File/Print**command is available for that purpose.