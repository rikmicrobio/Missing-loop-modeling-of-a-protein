**MODUS-PRIME-RGLAB**

**Installation and First-Run Guide**

*Native PDB protein modelling pipeline*

# 1. Overview

Modus-Prime-RGLAB starts from a local PDB-format structure. Run the
master script, select a PDB file in the file picker (when a graphical
desktop is available), choose a chain, and let the pipeline extract
sequences, identify gaps, align sequences, and build/rank models. The
original input PDB is not modified.

# 2. Requirements

  ---------------------------------------------------------------------
  Component                          Purpose / notes
  ---------------------------------- ----------------------------------
  Python                             Use a Python version supported by
                                     the installed MODELLER release.
                                     Keep the pipeline in the same
                                     environment as MODELLER.

  Conda                              Recommended for creating and
                                     managing an isolated environment.

  MODELLER                           Required for building and
                                     assessing protein models. A valid
                                     MODELLER license/key may be
                                     required; follow the official Sali
                                     Lab installation and licensing
                                     instructions.

  Clustal Omega (clustalo)           Required for automatic alignment.
                                     It can be skipped if a compatible
                                     pre-aligned FASTA file is supplied
                                     with \--aligned.

  Tkinter / desktop display          Needed for the graphical
                                     file-selection dialog. On headless
                                     Linux/WSL, use the command-line
                                     \--pdb option or configure GUI
                                     support.

  Input structure                    A PDB-format .pdb file. This
                                     pipeline does not accept mmCIF as
                                     its native template input.
  ---------------------------------------------------------------------

# 3. Create an environment

Open a terminal with Conda available. Create an environment using a
Python version supported by your MODELLER package:

> conda create -n modus-prime python=3.10 -y
>
> conda activate modus-prime

Python 3.10 is an example starting point, not a guarantee for every
MODELLER release. If Conda reports a dependency conflict, choose a
Python version supported by the MODELLER build available for your
operating system.

# 4. Install MODELLER

Use the official MODELLER installation instructions and the package
channel/version appropriate for your platform. For Conda installations,
a common starting command is:

> conda install -c salilab modeller

If that package is unavailable or incompatible, consult the official
MODELLER documentation rather than mixing packages from unrelated Python
environments. Complete any required license/key setup.

# 5. Install Clustal Omega

With the modus-prime environment activated, install Clustal Omega
through Conda:

> conda install -c conda-forge -c bioconda clustalo -y

Clustal Omega is not needed when you provide an aligned FASTA file using
the \--aligned option.

# 6. Verify the installation

> python \--version
>
> python -c \"from modeller import Environ; Environ(); print(\'MODELLER
> OK\')\"
>
> clustalo \--version

If the MODELLER test fails, check that the active Conda environment is
the one where MODELLER was installed and licensed. If clustalo is not
found, confirm its installation and that the environment\'s executable
directory is on PATH.

# 7. Run the pipeline

Extract the Modus-Prime-RGLAB ZIP archive. In a terminal, move into the
folder containing run_pipeline.py and activate the environment:

> cd path/to/modeler_automated
>
> conda activate modus-prime
>
> python run_pipeline.py

The console displays the Modus-Prime-RGLAB branding. When a graphical
desktop is available, select the local .pdb file in the file picker.

# 8. Command-line alternatives

Select a specific file and chain without using the file picker:

> python run_pipeline.py \--pdb path/to/native_structure.pdb \--chain A
>
> python run_pipeline.py \--pdb path/to/native_structure.pdb \--chain A
> \--loop \--models 5 \--loop-models 20
>
> python run_pipeline.py \--pdb path/to/native_structure.pdb \--chain A
> \--target full_length.fasta

Use \--target when the supplied PDB is missing unresolved residues but a
full-length FASTA sequence is available. If the PDB has SEQRES records,
the pipeline can use them for the full sequence; otherwise it extracts
the observed sequence from the structure. Only one chain is processed
per run.

# 9. Output files

  ---------------------------------------------------------------------
  Output                             Description
  ---------------------------------- ----------------------------------
  \*\_run/                           Output folder for the selected
                                     structure and run-specific files.

  \*.seq / \*\_seqres.fasta          Extracted template sequence and,
                                     when available, full sequence from
                                     SEQRES.

  \*\_alignment.ali                  MODELLER alignment used for model
                                     building.

  \*\_gaps.json                      Detected missing regions/gaps.

  \*\_modelled_ranking.tsv           Model ranking and available DOPE /
                                     GA341 scores.

  \*\_modelled_best.pdb              Best-ranked model produced by the
                                     pipeline.

  logs/                              Logs for individual pipeline
                                     steps.
  ---------------------------------------------------------------------

# 10. Troubleshooting

-   **No file-selection window:** Run on a desktop-enabled environment.
    On headless Linux/WSL, provide the path with \--pdb.

-   **No module named \'modeller\':** Activate the environment where
    MODELLER is installed; confirm package and Python compatibility.

-   **clustalo not found:** Install Clustal Omega in the active
    environment, or provide an aligned FASTA with \--aligned.

-   **PDB has no SEQRES records:** The pipeline can extract observed
    residues, but unresolved regions require a full-length sequence
    supplied with \--target.

-   **Long or missing loops:** Modelling uncertainty increases for long
    gaps, especially gaps longer than roughly 10--15 residues; inspect
    the output and validate models before downstream use.

# 11. Scientific-use note

Generated structures are computational predictions, not experimentally
verified structures. Review logs, alignments, gap definitions, model
scores, and structural quality before using models for docking,
publication, or other downstream analyses.

*Modus-Prime-RGLAB \| Installation Note*

## Optional: create the environment from the repository file

From the repository root, you can also try:

```bash
conda env create -f environment.yml
conda activate modus-prime
```

Conda package availability and MODELLER compatibility vary by operating system. If environment creation fails, follow the manual installation steps above and choose a Python version supported by the MODELLER build available for your platform. MODELLER licensing still needs to be configured separately.
