# Free Energy Pipeline for Polymer Monomers

A closed loop for learning the **size-independent thermodynamic stability** of polymer
monomers, so that a generative model can be steered toward molecules that are actually
stable.

```
pool (2.06M monomers)
   └─ seed selection ──> 3,915 molecules
          └─ DFT (ORCA) ──> ΔG
                 └─ FS5 + Ridge ──> ΔG residual, the size-independent signal
                        └─ GIN ensemble ──> prediction + uncertainty over the whole pool
                               └─ active learning ──> next batch ──┐
                                                                    │
                        ┌───────────────────────────────────────────┘
                        └─> back to DFT
```

Raw formation free energy scales with molecule size, so it cannot be compared across
molecules and is useless as a guidance signal — it would only push a generator toward
smaller or larger structures. Subtracting a composition baseline leaves the part that
reflects structure, which is what the loop above learns.

The DFT half (`scripts/`) converts SMILES to 3D structures, runs two-step geometry
optimisation and frequencies in [ORCA](https://orcaforum.kofo.mpg.de/), and computes ΔG
and its residual. The learning half (`ml/`) builds the candidate pool, selects what to
compute, trains the GNN, and runs acquisition.

> The name says "polyimide" for historical reasons: the first dataset was 1,077
> polyimides. That set is 19.3% imide against 1.8% in the pool it was meant to explore,
> covered only 31% of it, and has been retired. The current pool spans all four polymer
> databases.

---

## Pipeline Overview

```
┌─────────────────────────────────────────────────────────────────────────┐
│                     Free Energy Calculation Pipeline                    │
├─────────────────────────────────────────────────────────────────────────┤
│                                                                         │
│  ┌──────────────┐    ┌──────────────┐    ┌──────────────────────────┐   │
│  │  SMILES CSV   │───▶│  RDKit 3D    │───▶│  Step 1: PBE/def2-SVP   │   │
│  │ (miss_point)  │    │  XYZ coords  │    │  Geometry Optimization  │   │
│  └──────────────┘    └──────────────┘    └───────────┬──────────────┘   │
│                                                       │                  │
│  ┌──────────────────────────────────────────────────────────────────┐   │
│  │                              ▼                                    │   │
│  │  ┌───────────────────────────┐    ┌──────────────────────────┐   │   │
│  │  │  Step 2: B3LYP-D3BJ/     │───▶│  Frequency Calculation   │   │   │
│  │  │  def2-TZVP Optimization  │    │  (298.15 K)              │   │   │
│  │  └───────────────────────────┘    └───────────┬──────────────┘   │   │
│  └──────────────────────────────────────────────────────────────────┘   │
│                                                       │                  │
│  ┌───────────────────────────────────────────────────┐                  │
│  │                        ▼                           │                  │
│  │  ┌──────────────────┐    ┌──────────────────────┐ │                  │
│  │  │  Extract Gibbs   │───▶│  Compute ΔG with     │ │                  │
│  │  │  Free Energy     │    │  Atomization Ref     │ │                  │
│  │  └──────────────────┘    └───────────┬──────────┘ │                  │
│  └────────────────────────────────────────┼──────────┘                  │
│                                            ▼                             │
│                        ┌────────────────────────────────────┐            │
│                        │  Step 6: FS5 + Ridge baseline      │            │
│                        │  ΔG_residual = ΔG − f(composition) │            │
│                        │  → SIZE-INDEPENDENT stability      │            │
│                        └────────────────────────────────────┘            │
│                                                                         │
└─────────────────────────────────────────────────────────────────────────┘
```

**Why Step 6 matters.** Raw ΔG scales ~linearly with molecule size, so it cannot be
compared across molecules of different size, and it is useless as a reward/guidance
signal (it would just favour bigger or smaller molecules). Step 6 subtracts a linear
composition baseline (FS5: element counts + bond-element-pair × bond-type counts +
ring size × aromaticity, fitted with Ridge), leaving a **size-independent measure of
thermodynamic stability**. On the full dataset the residual has r(size) ≈ −0.001
while retaining a σ ≈ 8.9 kcal/mol structural signal. See the module docstring of
`scripts/compute_residual_deltaG.py` for the feature-set comparison that selected FS5
and for why nonlinear baselines are rejected.

---

## Repository Structure

```
Free_Energy/
├── README.md                          # This file
├── .gitignore
│
├── submit_free_energy.csh             # SGE job script (runs full pipeline for ONE molecule)
├── full_submit_free_energy.sh         # Batch submission wrapper (submits all molecules)
├── submit_atom_ref.csh                # SGE job: single-atom reference energies
│
├── scripts/                           # Python scripts for each pipeline step
│   ├── generate_xyz.py                # Step 0: SMILES → 3D XYZ coordinates
│   ├── generate_step1_opt_inp.py      # Step 1: Generate PBE optimization inputs
│   ├── generate_step2_opt_inp_from_xyz.py  # Step 2: Generate B3LYP optimization inputs
│   ├── generate_freq_inp.py           # Step 3: Generate frequency calculation inputs
│   ├── extract_thermo.py              # Step 4: Extract Gibbs energies from output
│   ├── compute_deltaG.py              # Step 5: Compute ΔG with atomization reference
│   ├── compute_residual_deltaG.py     # Step 6: FS5+Ridge → SIZE-INDEPENDENT ΔG residual
│   │
│   ├── atom_ref.csv                   # Atomic reference energies USED BY STEP 5 (required)
│   ├── generate_atom_ref_inp.py       # Generate single-atom ORCA inputs (13 elements × 2 methods)
│   ├── parse_atom_ref.py              # Parse those outputs → atom_ref.csv
│   ├── compare_atom_ref_methods.py    # Compare B3LYP-D3BJ vs ωB97X-D3 reference sets
│   │
│   ├── bond_count.py                  # Utility: count backbone/total bonds
│   ├── plot.py                        # Analysis: scatter plots of ΔG vs atom count
│   ├── tsne.py                        # Analysis: t-SNE visualization of molecules
│   └── tsne_all.py                    # Analysis: t-SNE with background polymer set
│
├── ml/                                # Learning half: pool, selection, GNN, acquisition
│   ├── build_clean_pool.py            # Chemistry + element filter over the raw pool
│   ├── seed_select.py                 # Iteration-0 seed: 4 components (see seeds/README)
│   ├── active_learning_v3.py          # Acquisition: uncertainty-weighted max-coverage
│   ├── gin_residual_v2_cv.py          # Train the GIN ensemble on the ΔG residual
│   ├── predict_dataset_gin.py         # Ensemble inference over the whole pool
│   ├── tau_calibration.py             # Where tau=0.4 comes from (measured, not chosen)
│   ├── run_residual.sh                # Run step 6 against a dataset outside this repo
│   └── sync_check.sh                  # Detect drift against the analysis workspace
│
├── pool/
│   ├── pool_molecules.csv.gz          # 2,057,755 candidate monomers (14 MB)
│   └── README.md                      # What was filtered out and why
│
├── seeds/
│   ├── seed_v9.csv                    # 3,915 molecules: the current DFT campaign
│   └── README.md                      # Why four components, with the measurements
│
└── data/                              # Data directory (see "Data Directory Guide" below)
    ├── input_molecules.csv            # THE INPUT LIST: columns PID, smiles
    ├── atom_ref/                      # Atomic reference energies
    ├── xyz/                           # RDKit-generated 3D structures
    ├── opt_inp/                       # ORCA optimization input/output files
    ├── opt_out/                       # ORCA optimization log outputs
    ├── freq_inp/                      # ORCA frequency input files
    └── freq_out/                      # ORCA frequency output files
```

## Running a campaign

```bash
# 1. Turn a seed set into the pipeline's input list
python -c "import pandas as pd; d=pd.read_csv('seeds/seed_v9.csv'); \
           d['PID']='SD'+d['rank'].astype(str); \
           d[['PID','smiles']].to_csv('data/input_molecules.csv', index=False)"

# 2. 3D structures (parallel; conformer search is on by default)
python scripts/generate_xyz.py                      # N_WORKERS=16 to go faster

# 3. Submit one ORCA job per molecule
NUMBER=4000 ./full_submit_free_energy.sh

# 4. When they finish: thermo -> ΔG -> residual
python scripts/extract_thermo.py
python scripts/compute_deltaG.py
python scripts/compute_residual_deltaG.py
```

`data/input_molecules.csv` is the one name the DFT half agrees on: `generate_xyz.py`
reads it, `full_submit_free_energy.sh` submits from it, and every output file is named
after its `PID`, so the two stay in step. Override with `INPUT_CSV=...` for a
second batch rather than editing either script.

**Hold out the validation block.** 1,000 of the 3,915 seed rows have
`is_validation=True`. They are a uniform random draw, which makes them the only
unbiased measure of progress across active-learning rounds — every other component is
chosen by a criterion correlated with the model or the coverage. Keep them out of GNN
*training*; they should still be included when *fitting the baseline*, where each point
carries only ~0.025 leverage.

---

## Data Directory Guide

The `data/` directory is shipped empty. Below is a description of what files each subdirectory should contain after running the pipeline.

| Subdirectory | Contents | Generated By | Example Files |
|---|---|---|---|
| `data/atom_ref/` | Atomic reference energies CSV | **User-provided** | `atom_ref.csv` (columns: `atom`, `energy` in Hartree) |
| `data/xyz/` | 3D molecular coordinates | `scripts/generate_xyz.py` | `PI1.xyz`, `PI2.xyz`, ... |
| `data/opt_inp/` | ORCA optimization input files + ORCA-generated intermediate files | `scripts/generate_step1_opt_inp.py`, `scripts/generate_step2_opt_inp_from_xyz.py` | `PI1_step1_pbe_opt.inp`, `PI1_step2_b3lyp_opt.inp`, `PI1_step1_pbe_opt.xyz` (optimized geometry) |
| `data/opt_out/` | ORCA optimization log outputs | ORCA (via `submit_free_energy.csh`) | `PI1_step1_pbe_opt.out`, `PI1_step2_b3lyp_opt.out` |
| `data/freq_inp/` | ORCA frequency input files + intermediate files | `scripts/generate_freq_inp.py` | `PI1_freq.inp` |
| `data/freq_out/` | ORCA frequency output files | ORCA (via `submit_free_energy.csh`) | `PI1_freq.out` |

**Additionally, the following CSV files will be generated in `data/`:**

| File | Description | Generated By |
|---|---|---|
| `miss_point.csv` | **Input file** — monomer IDs and SMILES strings | **User-provided** (columns: `monomer_ID`, `smiles`) |
| `deltaG_raw.csv` | Extracted raw Gibbs free energies | `scripts/extract_thermo.py` (columns: `mol`, `Gibbs_Eh`, `G_minus_Eel`) |
| `failed_monomers.csv` | SMILES that failed 3D conversion | `scripts/generate_xyz.py` (only if failures occur) |

---

## Prerequisites

### Software

| Software | Version | Purpose |
|---|---|---|
| [ORCA](https://orcaforum.kofo.mpg.de/) | 6.1.0 | DFT calculations (geometry optimization + frequency) |
| Python | 3.8+ | Pipeline scripts |
| SGE (Sun Grid Engine) | — | Job scheduler (for HPC clusters) |

### Python Packages

```bash
pip install pandas rdkit numpy scikit-learn matplotlib seaborn
```

Or create a conda environment:

```bash
conda create -n free_energy python=3.10 pandas rdkit numpy scikit-learn matplotlib seaborn
conda activate free_energy
```

---

## Step-by-Step Usage

### 0. Prepare Input Data

Create `data/miss_point.csv` with your monomer SMILES:

```csv
monomer_ID,smiles
PI1,*CC*
PI2,*C(C*)C
PI3,*C(C*)CC
...
```

> **Note:** The `*` symbols represent polymer repeat-unit attachment points and will be automatically replaced with `C` (methyl caps) during processing.

Atomic reference energies live in **`scripts/atom_ref.csv`** (this is the path
`compute_deltaG.py` reads) with columns `atom,energy` in Hartree:

```csv
atom,energy
C,-37.838153520821
H,-0.498764293374
O,-75.066887071071
N,-54.579005649305
...
```

These are shipped with the repo for all 13 supported elements
(H, C, N, O, F, Si, P, S, Cl, Ge, Br, Sn, I). **They must be computed with the same
functional/basis as the molecules** (`B3LYP D3BJ def2-TZVP`), otherwise ΔG is
meaningless — the atoms are subtracted from the molecular Gibbs energy.

To regenerate them (e.g. to add an element, or to switch functional):

```bash
# 1. Write ORCA inputs for every atom (ground-state multiplicities are built in:
#    H=2, C=3, N=4, O=3, F=2, Si=3, P=4, S=3, Cl=2, Ge=3, Br=2, Sn=3, I=2)
python scripts/generate_atom_ref_inp.py

# 2. Run them all (serial -- a single atom has ~10-30 basis functions, so MPI
#    parallelism is counter-productive; the inputs deliberately contain no %pal)
qsub submit_atom_ref.csh

# 3. Parse -> data/atom_ref/atom_ref_{b3lyp,wb97x}.csv, then write the production file
python scripts/parse_atom_ref.py --write-production --use E
```

`parse_atom_ref.py` prints an old-vs-new comparison and refuses to write if any
energy failed to parse. It defaults to the **electronic** energy (`--use E`) to stay
consistent with the existing table; `--use G` would switch to Gibbs and shift every
molecule's ΔG by a per-element constant.

`scripts/compare_atom_ref_methods.py` quantifies the difference between two
functional sets (a `wB97X-D3/def2-TZVPPD` set is generated alongside B3LYP for
comparison — **do not mix the two**).

### 1. Generate 3D XYZ Coordinates

```bash
python scripts/generate_xyz.py
```

- **Input:** `data/miss_point.csv`
- **Output:** `data/xyz/{monomer_ID}.xyz`
- Uses RDKit ETKDGv3 conformer generation + UFF force field optimization
- Failed molecules are saved to `data/failed_monomers.csv`

### 2. Generate Step 1 Optimization Inputs (PBE/def2-SVP)

```bash
python scripts/generate_step1_opt_inp.py
```

- **Input:** `data/xyz/*.xyz`
- **Output:** `data/opt_inp/{monomer_ID}_step1_pbe_opt.inp`
- ORCA keywords: `! PBE def2-SVP TightSCF Opt RIJCOSX def2/J`
- Uses 8 parallel cores (`%pal nprocs 8 end`)

### 3. Run the Full Pipeline (Per Molecule via SGE)

Submit a single molecule:

```bash
qsub -v MOL="PI1" -N "FE_PI1" submit_free_energy.csh
```

Or submit all molecules in batch:

```bash
bash full_submit_free_energy.sh
```

**What `submit_free_energy.csh` does for each molecule:**

1. **Step 1 — PBE/def2-SVP optimization** (fast pre-optimization)
2. **Step 2 — B3LYP-D3BJ/def2-TZVP optimization** (accurate geometry)
3. **Frequency calculation** at 298.15 K (thermodynamic corrections)
4. **Extract thermodynamic data** (Gibbs free energy)
5. **Compute ΔG** relative to atomic reference energies

### 4. Post-Processing & Analysis

After all jobs complete:

```bash
# Step 4 — extract thermodynamic data from frequency outputs
#          → scripts/merged_G_raw.csv
python scripts/extract_thermo.py

# Step 5 — atomization ΔG for all molecules (needs scripts/atom_ref.csv)
#          → scripts/final_data_with_deltaG.csv
python scripts/compute_deltaG.py

# Step 6 — size-independent residual (FS5 + Ridge). THIS is the quantity to
#          compare across molecules or feed to an ML / generative model.
#          → scripts/final_data_with_residual_deltaG.csv + scripts/residual_plots/
python scripts/compute_residual_deltaG.py

# Generate scatter plots
python scripts/plot.py

# Generate t-SNE chemical space visualization
python scripts/tsne.py
python scripts/tsne_all.py  # with background polymer set
```

---

## Computational Details

### Two-Step Geometry Optimization

| Step | Method | Basis Set | Dispersion | Aux Basis | Purpose |
|---|---|---|---|---|---|
| 1 | PBE | def2-SVP | — | def2/J (RIJCOSX) | Fast rough optimization |
| 2 | B3LYP | def2-TZVP | D3BJ | def2/J (RIJCOSX) | Accurate final geometry |

### Frequency Calculation

| Property | Value |
|---|---|
| Method | B3LYP-D3BJ/def2-TZVP |
| Temperature | 298.15 K |
| Aux Basis | def2/J (RIJCOSX) |
| SCF Convergence | TightSCF |

### ΔG Calculation

The atomization free energy is computed as:

```
ΔG(molecule) = G(molecule) − Σ nᵢ × E(atomᵢ)
```

Where:
- `G(molecule)` = Gibbs free energy from the frequency calculation (in Hartree)
- `nᵢ` = number of atom type `i` in the molecule
- `E(atomᵢ)` = DFT reference energy for isolated atom `i` (from `atom_ref.csv`)

### Normalized ΔG Metrics

The pipeline also computes several normalized ΔG values:

| Metric | Description |
|---|---|
| `DeltaG_per_heavy_atom` | ΔG divided by number of heavy (non-H) atoms |
| `DeltaG_per_atom` | ΔG divided by total number of atoms |
| `DeltaG_per_backbone_bond` | ΔG divided by number of heavy-atom bonds |
| `DeltaG_per_bond` | ΔG divided by total number of bonds |
| `DeltaG_per_CH` | ΔG divided by number of C–H bonds |

---

## Analysis Scripts

### `scripts/plot.py`

Generates scatter plots of normalized ΔG values vs. total atom count, useful for identifying size-dependent trends in thermodynamic stability.

### `scripts/tsne.py`

Performs t-SNE dimensionality reduction on Morgan fingerprints (radius=2, 1024 bits) and generates colored scatter plots for each ΔG metric. This helps visualize chemical space coverage.

### `scripts/tsne_all.py`

Similar to `tsne.py`, but overlays the computed molecules on a background set of polymers (`data/all_polymer.csv`) to show where the target polyimides sit within the broader polymer chemical space.

### `scripts/bond_count.py`

Utility script for counting backbone bonds, total bonds, and C–H bonds from SMILES strings. Useful for understanding normalization denominators.

---

## Output Files

After a complete pipeline run, the final results are saved to:

| File | Description |
|---|---|
| `scripts/merged_G_raw.csv` | Step 4 output: `mol, smiles, Gibbs_Eh, G_minus_Eel` (intermediate) |
| `scripts/final_data_with_deltaG.csv` | Step 5 output: SMILES, atom/bond counts, Gibbs energy, **ΔG** and all normalized metrics |
| `scripts/final_data_with_residual_deltaG.csv` | **Step 6 output — the main result.** Adds `Delta_G_predicted` (FS5 Ridge baseline) and **`Delta_G_residual`**, the size-independent stability signal, in kcal/mol. Net-charged molecules get `NaN` residual (the DFT pipeline assumes neutral singlets) so downstream `dropna()` excludes them |
| `scripts/residual_plots/` | Step 6 diagnostics: residual vs atom count / MW (both should be ≈ flat), parity plots, residual histogram |
| `data/deltaG_raw.csv` | Raw extracted Gibbs energies (intermediate) |
| `tsne_plots/` | t-SNE scatter plots colored by each ΔG metric |
| `tsne_plots_with_background/` | t-SNE plots with background polymer set |
| `deltaG_scatter_plot.png` | Scatter plot of energy vs. atom count |

---

## Customization

### Changing the DFT Method

Edit the ORCA input templates in:
- `scripts/generate_step1_opt_inp.py` — Step 1 method (default: PBE/def2-SVP)
- `scripts/generate_step2_opt_inp_from_xyz.py` — Step 2 method (default: B3LYP-D3BJ/def2-TZVP)
- `scripts/generate_freq_inp.py` — Frequency method (default: B3LYP-D3BJ/def2-TZVP)

### Changing Parallelization

All ORCA inputs use `%pal nprocs 8 end`. To change core count:
1. Update the template strings in the Python scripts
2. Update `#$ -pe smp 8` in `submit_free_energy.csh`

### Changing Temperature

Edit the `temp` variable in `scripts/generate_freq_inp.py` (default: 298.15 K).

---

## HPC Job Submission

The pipeline is designed for **SGE (Sun Grid Engine)** clusters:

```bash
# Edit full_submit_free_energy.sh to set:
#   - file_path: path to your input CSV
#   - start: row to start from (2 = skip header)
#   - number: max number of molecules to submit

bash full_submit_free_energy.sh
```

Each molecule runs as an independent job requesting:
- 8 SMP cores
- Long queue

---

## License

This project is part of ongoing research. Please contact the authors before use.
