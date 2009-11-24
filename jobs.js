var jobs = [
    {
        menuName: "&Integrate with d*TREK",
        jobName: "Integrate",
        executable: "python",
        arguments: [ mifit.directory() + "/MIExpert/mi_integrate_ui.py" ]
    },
    {
        menuName: "&SAD Phasing",
        jobName: "SAD Phasing",
        executable: "python",
        arguments: [ mifit.directory() + "/MIExpert/mi_sadphase_ui.py" ]
    },
    {
        menuName: "&Molecular Replacement",
        jobName: "Molecular Replacement",
        executable: "python",
        arguments: [ mifit.directory() + "/MIExpert/mi_molrep_ui.py" ]
    },
    {
        menuName: "&Refinement",
        jobName: "Refinement",
        executable: "python",
        arguments: [ mifit.directory() + "/MIExpert/mi_integrate_ui.py" ]
    },
    {
        menuName: "C&ocrystal Solution",
        jobName: "Cocrystal Solution",
        executable: "python",
        arguments: [ mifit.directory() + "/MIExpert/mi_bng_ui.py" ]
    },
    {
        menuName: "Re&port",
        jobName: "Report",
        executable: "python",
        arguments: [ mifit.directory() + "/MIExpert/mi_deposit3d_ui.py" ]
    },
    {
        menuName: "&NCS Modeling",
        jobName: "NCS Modeling",
        executable: "python",
        arguments: [ mifit.directory() + "/MIExpert/mi_ncsmodeler_ui.py" ]
    },
    {
        menuName: "Cocr&ystal Superposition",
        jobName: "Cocrystal Superposition",
        executable: "python",
        arguments: [ mifit.directory() + "/MIExpert/mi_ligandoverlap_ui.py" ]
    },
    {
        menuName: "Convert CIF to &SHELX",
        jobName: "Convert CIF to SHELX",
        executable: "python",
        arguments: [ mifit.directory() + "/MIExpert/mi_convertlib_ui.py", "--refprogram", "shelx" ]
    },
    {
        menuName: "Convert CIF to &CNX",
        jobName: "Convert CIF to CNX",
        executable: "python",
        arguments: [ mifit.directory() + "/MIExpert/mi_convertlib_ui.py", "--refprogram", "cns" ]
    }
];

for (var i = 0; i < jobs.length; ++i)
{
    var job = jobs[i]
    mifit.addJob(job.menuName, job.jobName, job.executable, job.arguments)
}
