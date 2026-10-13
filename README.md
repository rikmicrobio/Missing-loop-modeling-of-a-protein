<div align="center">

# Modus-Prime-RGLAB
### Native PDB Protein Modelling Pipeline

**Select a native PDB structure → identify missing regions → build and rank protein models**

</div>

---

## Table of Contents

- [Overview](#overview)
- [Pipeline workflow](#pipeline-workflow)
- [Requirements](#requirements)
- [Installation](#installation)
- [Configure MODELLER](#configure-modeller)
- [Run the pipeline](#run-the-pipeline)
- [Using your own full-length sequence](#using-your-own-full-length-sequence)
- [Optional settings](#optional-settings)
- [Understanding the output](#understanding-the-output)
- [Troubleshooting](#troubleshooting)
- [Scientific considerations](#scientific-considerations)

## Overview

**Modus-Prime-RGLAB** is a Python workflow for modelling protein structures from a native/template **PDB-format** structure. It extracts the observed template sequence, uses the PDB's `SEQRES` records when available (or a user-provided FASTA sequence), aligns sequences, identifies missing regions, builds models with MODELLER, and ranks the resulting structures.

The pipeline displays the Modus-Prime-RGLAB ASCII banner in green when launched. The original input structure is not modified; run files are written to a separate working directory.

> **Input format:** PDB (`.pdb` or `.ent`) with `ATOM` records. mmCIF files are not accepted as input by this workflow.

## Pipeline workflow

```text
Native PDB structure
        │
        ▼
Select protein chain
        │
        ▼
Extract observed template sequence
        │
        ├── Use PDB SEQRES sequence, when available
        └── Or provide a full-length FASTA sequence
        │
        ▼
Sequence alignment (Clustal Omega or supplied alignment)
        │
        ▼
Identify missing regions / gaps
        │
        ▼
Build comparative models (MODELLER)
        │
        ├── Optional: refine missing regions with LoopModel
        ▼
Rank models and save output structures
```

## Requirements

| Requirement | Purpose |
|---|---|
| Python | Runs the pipeline scripts; use a version compatible with your MODELLER installation |
| Conda | Recommended for managing the environment |
| MODELLER | Required for comparative modelling; follow the official installation and licence instructions |
| Clustal Omega (`clustalo`) | Used for sequence alignment unless you provide an aligned FASTA file |
| Tkinter / desktop GUI | Opens the PDB file picker; if unavailable, the script prompts for the file path |
| PDB structure | A local `.pdb`/`.ent` file or a valid four-character RCSB PDB ID |

A GPU is **not required** for the standard workflow.

## Installation

### 1. Download the repository

Clone the GitHub repository:

```bash
git clone <YOUR-GITHUB-REPOSITORY-URL>
cd <YOUR-REPOSITORY-FOLDER>
```

Alternatively, download the repository as a ZIP from GitHub and extract it. Change into the directory that contains `run_pipeline.py`.

### 2. Create and activate a Conda environment

Start with a dedicated environment. Python 3.10 is an example; the correct version depends on the MODELLER release available for your operating system.

```bash
conda create -n modus-prime python=3.10 -y
conda activate modus-prime
```

If MODELLER cannot be installed in this environment, use a Python version supported by the MODELLER package you are installing. Avoid mixing packages across unrelated Python environments.

### 3. Install MODELLER

Install MODELLER using the official instructions for your operating system and release. A Conda installation may be available with:

```bash
conda install -c salilab modeller
```

If that command is unavailable or produces dependency conflicts, consult the official [MODELLER documentation](https://salilab.org/modeller/) for the appropriate package and installation route.

**Licence setup:** MODELLER requires a valid licence key, including for academic use. Obtain and configure the key as instructed by MODELLER. The pipeline checks that MODELLER can start before beginning the modelling steps.

### 4. Install Clustal Omega

With the `modus-prime` environment active:

```bash
conda install -c conda-forge -c bioconda clustalo -y
```

Check that it is available:

```bash
clustalo --version
```

Clustal Omega can be omitted if you supply a compatible aligned FASTA file with `--aligned`.

### 5. Check Tkinter / file selection

Test Tkinter:

```bash
python -c "import tkinter; print('Tkinter OK')"
```

If the test fails on Ubuntu/WSL, the operating-system package may be required:

```bash
sudo apt update
sudo apt install python3-tk
```

This installs Tkinter for the system Python; a Conda Python may need separate configuration. On systems without a graphical desktop, the pipeline can prompt you to enter the PDB path in the terminal instead.

### 6. Verify the main dependencies

Run these commands in the activated environment:

```bash
python --version
python -c "from modeller import Environ; Environ(); print('MODELLER OK')"
clustalo --version
python -c "import tkinter; print('Tkinter OK')"
```

The MODELLER test should only be considered successful if the licence is configured and `Environ()` starts without an error. Tkinter is optional when you use the command-line `--pdb` option.

## Run the pipeline

### Option A — Select a local PDB interactively

From the directory containing `run_pipeline.py`:

```bash
conda activate modus-prime
python run_pipeline.py
```

1. The green Modus-Prime-RGLAB banner appears in the terminal.
2. Select your local `.pdb` or `.ent` file in the file picker.
3. The pipeline uses the first available chain by default.
4. The pipeline runs the sequence, alignment, gap-identification, modelling, and ranking steps.
5. Results are saved in a run directory, normally named `<PDB_BASENAME>_run/`.

If a graphical file picker cannot open, enter the PDB path when prompted.

### Option B — Provide the PDB path and chain

```bash
python run_pipeline.py --pdb path/to/native_structure.pdb --chain A
```

Replace `path/to/native_structure.pdb` with your actual file path and `A` with a chain present in the structure.

On Windows, quote paths containing spaces:

```powershell
python run_pipeline.py --pdb "C:\Users\YourName\Desktop\protein structure.pdb" --chain A
```

### Option C — Use a PDB ID

If the supplied value is not an existing local file and matches a four-character PDB ID, the pipeline attempts to download the structure from RCSB PDB. Internet access is required.

```bash
python run_pipeline.py --pdb 6I07 --chain A
```

## Using your own full-length sequence

A PDB structure may not contain coordinates for every residue. If the PDB has `SEQRES` records, the pipeline uses them as the full sequence by default. If `SEQRES` records are absent, it extracts the observed sequence from the structure; unresolved residues cannot be recovered from coordinates alone.

To model against a full-length sequence, prepare a FASTA file, for example `full_length.fasta`:

```fasta
>target_protein
MKK...REPLACE_WITH_THE_FULL_LENGTH_AMINO_ACID_SEQUENCE...K
```

Run:

```bash
python run_pipeline.py \
  --pdb path/to/native_structure.pdb \
  --chain A \
  --target full_length.fasta
```

On Windows PowerShell, the same command can be entered on one line:

```powershell
python run_pipeline.py --pdb "C:\path\to\native_structure.pdb" --chain A --target "C:\path\to\full_length.fasta"
```

Use the complete amino-acid sequence for the intended target. Confirm that the selected chain and target sequence correspond to the same protein.

## Optional settings

### Build multiple models and refine loops

The pipeline builds five models by default. To enable loop refinement:

```bash
python run_pipeline.py \
  --pdb path/to/native_structure.pdb \
  --chain A \
  --target full_length.fasta \
  --loop \
  --models 5 \
  --loop-models 20
```

This requests five comparative models and 20 loop-refinement models per model. Loop refinement can increase runtime.

### Common command-line options

| Option | Description |
|---|---|
| `--pdb FILE_OR_ID` | Local PDB file or four-character RCSB PDB ID. If omitted, a file picker opens |
| `--chain A` | Select the chain to process; defaults to the first available protein chain |
| `--target FILE.fasta` | Full-length target sequence in FASTA format |
| `--aligned FILE` | Use a supplied aligned FASTA (template sequence first); skips Clustal Omega |
| `--workdir DIR` | Choose the output/work directory |
| `--models N` | Number of comparative models to build (default: `5`) |
| `--loop` | Enable missing-region refinement with MODELLER `LoopModel` |
| `--loop-models N` | Loop-refinement models per model (default: `10`) |
| `--md-level LEVEL` | Modelling refinement level: `very_fast`, `fast`, `slow`, or `very_slow` |
| `--flank N` | Additional residues to refine on each side of a gap (default: `2`) |
| `--include-termini` | Also refine missing N- and C-terminal tails |
| `--start-at STEP` | Resume at `getseq`, `align`, or `model`; required earlier output files must already exist |
| `--install-missing` | Ask the script to install Clustal Omega with Conda if it is missing |
| `--verbose` | Display full output from each pipeline step |
| `--soap` | Also calculate SOAP-Protein-OD; requires its supporting library |

To view the available options:

```bash
python run_pipeline.py --help
```

Individual steps also support help, for example:

```bash
python getseq.py --help
python make_alignment.py --help
python model-single.py --help
```

## Understanding the output

For a local structure such as `native_structure.pdb`, the default work directory is `native_structure_run/`. The exact file prefix depends on the selected PDB file or PDB ID.

| Output | Description |
|---|---|
| `<code>.pdb` | Staged working copy of the input structure; the original is left unchanged |
| `<code>.seq` | Observed template sequence |
| `<code>_seqres.fasta` | Full sequence extracted from PDB `SEQRES` records, when available |
| `<code>_alignment.ali` | MODELLER alignment used to build models |
| `<code>_gaps.json` | Detected missing regions/gaps |
| `<code>_modelled.B*.pdb` | Individual generated models |
| `<code>_modelled_ranking.tsv` | Model ranking, including available DOPE and GA341 scores |
| `<code>_modelled_best.pdb` | Best-ranked model according to the pipeline's ranking rule |
| `logs/` | Step-specific logs, useful for diagnosing errors |

**Model ranking:** lower DOPE scores are generally preferred when comparing models of the same target under comparable conditions. GA341 values closer to 1 can support confidence in the fold, but neither score replaces structural validation.

## Troubleshooting

### `Modeller could not start`

- Confirm that `modus-prime` is activated.
- Confirm that MODELLER is installed in that exact environment.
- Check that your MODELLER licence key is correctly configured.
- Run the verification command from the [installation section](#6-verify-the-main-dependencies).

### `clustalo not found`

Install it in the active environment:

```bash
conda install -c conda-forge -c bioconda clustalo -y
```

Alternatively, supply a compatible aligned FASTA file with `--aligned`.

### The file picker does not open

- Check Tkinter with `python -c "import tkinter"`.
- Use a desktop-enabled terminal session, or enter the PDB path when prompted.
- For headless sessions, bypass the picker:

```bash
python run_pipeline.py --pdb /absolute/path/to/native_structure.pdb --chain A
```

### `chain not in PDB`

Check the chain identifiers present in the structure and use one of those values with `--chain`. The pipeline processes one chain per run.

### Missing residues or gaps are not reconstructed as expected

The coordinates in a PDB cannot reveal the identity of unresolved residues. Supply the full-length sequence with `--target full_length.fasta`. Inspect the alignment and `*_gaps.json` before interpreting the model.

### A pipeline step fails

Read the relevant file in the run directory's `logs/` folder. The terminal error should identify the failed step and its log path. Fix the underlying issue before restarting or resuming.

## Scientific considerations

- Models are computational predictions, not experimentally determined structures.
- Long loops and unresolved terminal regions can be especially uncertain; inspect them carefully.
- Review sequence alignment, gap definitions, model scores, stereochemistry, clashes, and other structural-quality metrics before using models for docking, publication, or downstream analyses.
- Run one chain at a time. Repeat with a different `--chain` value to process another chain.

---

<div align="center">

**Modus-Prime-RGLAB**  
*Native PDB Modelling Pipeline*

</div>
