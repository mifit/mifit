# Tutorial Lesson: Refinement with CCP4/REFMAC5 #

(This tutorial lesson requires that MIExpert is installed and the system has access to the CCP4 software suite)

_Scenario: You have a reasonably good protein model and wish to continue to finalize the structure with further refinement including water-picking and automated error detection._

1. Copy the files _automation\_out.pdb_and _automation\_out.mtz_from the _examples_directory of your MIFit installation into a new directory. File _automation\_out.pdb_is the model that will be refined. File _automation\_out.mtz_is the refinement data in MTZ format.

2. Select the **Job/Refinement**command from the MIFit interface. In the resulting interface use the **Browse...**button to select the directory defined in (1) as the **Working directory**. Use the associated **Browse...**buttons to load _automation\_out.pdb_as **Model (PDB)**and _automation\_out.mtz_as the **Data**parameter.

3. Leave the **Refinement type**parameter as **Refmac5** and the **B-Factor treatment**as **Isotropic**.

4. The **Weight**parameter controls the relative contributions of the diffraction data and the stereochemical restraints to the refinement process, with a lower value indicating a lower contribution from the data. This is a relatively high resolution structure (which allows larger weights) but is also near the end of the refinement (for which smaller weights are use) so set the weight to 0.1.

5. We are going to include water picking in this refinement so leave the **Number of cycles**parameter as 5. Select the **Water picking cycles**check boxand set the value to 2. With these settings, three refinement runs will be performed, with water picking between runs (_i.e._ the sequence refinement/water-pick/refinement/water-pick/refinement).

6. Leave the **Optional items**alone and click on the **OK**button. The refinement run will take a few minutes. You can check the progress of the refinement job by selecting the **Jobs** tab from the base of the navigation tree. Completion is indicated by the refinement job icon changing from yellow to green.

7. In a real-life application the next step would be load the output coordinates (_refine\_1.pdb_) and phased structure factor data (_refine\_1.mtz_) with the **File/Open models, data, maps, etc...**menu command or by using the **Open Results**... menu item available via a right mouse button click on the green refinement job icon under the **Jobs** tab. Although metrics for identifying structure misfits are somewhat fuzzy it is interesting to look through the amino acids in the associated error list file, _refine\_1\_errors.txt_. These items are also imported in to MIFit via special annotations in the PDB file but this behavior may be suppressed by selecting the **Display** tab at the base of the navigation tree and deselecting the **Annotations/Auto-show imported errors**

9. You may select **File/Close**or **File/Exit**to close this session or shut down MIFit.