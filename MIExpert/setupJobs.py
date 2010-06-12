import mifit

MIExpert_dir = mifit.directory() + "/MIExpert"

mifit.addJob("&Integrate with d*TREK",
             "Integrate",
             "python",
             [ MIExpert_dir + "/mi_integrate_ui.py" ],
             "")

mifit.addJob("&SAD Phasing",
             "SAD Phasing",
             "python",
             [ MIExpert_dir + "/mi_sadphase_ui.py" ],
             "")

mifit.addJob("&Molecular Replacement",
             "Molecular Replacement",
             "python",
             [ MIExpert_dir + "/mi_molrep_ui.py" ],
             "")

mifit.addJob("&Refinement",
             "Refinement",
             "python",
             [ MIExpert_dir + "/mi_refine_ui.py" ],
             "")

mifit.addJob("C&ocrystal Solution",
             "Cocrystal Solution",
             "python",
             [ MIExpert_dir + "/mi_bng_ui.py" ],
             "")

mifit.addJob("Re&port",
             "Report",
             "python",
             [ MIExpert_dir + "/mi_deposit3d_ui.py" ],
             "")

mifit.addJob("&NCS Modeling",
             "NCS Modeling",
             "python",
             [ MIExpert_dir + "/mi_ncsmodeler_ui.py" ],
             "")

mifit.addJob("Cocr&ystal Superposition",
             "Cocrystal Superposition",
             "python",
             [ MIExpert_dir + "/mi_ligandoverlap_ui.py" ],
             "")

mifit.addJob("Convert CIF to &SHELX",
             "Convert CIF to SHELX",
             "python",
             [ MIExpert_dir + "/mi_convertlib_ui.py", "--refprogram", "shelx" ],
             "")

mifit.addJob("Convert CIF to &CNX",
             "Convert CIF to CNX",
             "python",
             [ MIExpert_dir + "/mi_convertlib_ui.py", "--refprogram", "cns" ],
             "")

mifit.addJob("R&efmac5 restraints",
             "Refmac5 restraints",
             "python",
             [ MIExpert_dir + "/mi_restraints_ui.py" ],
             "")

