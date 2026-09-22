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
│  │  SMILES CSV  │───▶│  RDKit 3D    │───▶│  Step 1: PBE/def2-SVP   │   │
│  │ input_molec. │    │  XYZ coords  │    │  Geometry Optimization  │   │
│  └──────────────┘    └──────────────┘    └───────────┬──────────────┘   │
│                                                       │                  │
│  ┌──────────────────────────────────────────────────────────────────┐   │
│  │                              ▼                                    │   │
│  │  ┌───────────────────────────┐    ┌──────────────────────────┐   │   │
│  │  │  Step 2: wB97X-D3/       │───▶│  Frequency Calculation   │   │   │
│  │  │  def2-TZVP Optimization  │    │  (298.15 K, DefGrid3)    │   │   │
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
│   ├── generate_step2_opt_inp_from_xyz.py  # Step 2: Generate wB97X-D3 optimization inputs
│   ├── generate_freq_inp.py           # Step 3: Generate frequency calculation inputs
│   ├── extract_thermo.py              # Step 4: Extract Gibbs energies from output
│   ├── compute_deltaG.py              # Step 5: Compute ΔG with atomization reference
│   ├── compute_residual_deltaG.py     # Step 6: FS5+Ridge → SIZE-INDEPENDENT ΔG residual
│   │
│   ├── atom_ref.csv                   # Atomic reference energies USED BY STEP 5 (required)
│   ├── orca_settings.py               # SINGLE SOURCE OF TRUTH for the level of theory
│   ├── generate_atom_ref_inp.py       # Generate single-atom ORCA inputs (13 elements)
│   ├── parse_atom_ref.py              # Parse those outputs → atom_ref.csv
│   ├── compare_atom_ref_methods.py    # Compare two atomic-reference sets
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

Every file is named after its molecule's `PID`, which is how the stages stay in step:
`generate_xyz.py` writes `<PID>.xyz`, ORCA writes `<PID>_freq.out`, and
`extract_thermo.py` matches them back up by that id. Examples below use `SD1` for a
molecule from the seed set.

| Subdirectory | Contents | Generated By | Example Files |
|---|---|---|---|
| `data/atom_ref/` | Single-atom ORCA inputs/outputs, when regenerating the reference table | `scripts/generate_atom_ref_inp.py` | `C_atom.inp`, `C_atom.out` |
| `data/xyz/` | 3D molecular coordinates | `scripts/generate_xyz.py` | `SD1.xyz`, `SD2.xyz`, ... |
| `data/opt_inp/` | ORCA optimization input files + ORCA-generated intermediate files | `scripts/generate_step1_opt_inp.py`, `scripts/generate_step2_opt_inp_from_xyz.py` | `SD1_step1_pbe_opt.inp`, `SD1_step2_opt.inp`, `SD1_step1_pbe_opt.xyz` (optimized geometry) |
| `data/opt_out/` | ORCA optimization log outputs | ORCA (via `submit_free_energy.csh`) | `SD1_step1_pbe_opt.out`, `SD1_step2_opt.out` |
| `data/freq_inp/` | ORCA frequency input files + intermediate files | `scripts/generate_freq_inp.py` | `SD1_freq.inp` |
| `data/freq_out/` | ORCA frequency output files | ORCA (via `submit_free_energy.csh`) | `SD1_freq.out` |

> The **production** atomic reference table is `scripts/atom_ref.csv`, not
> `data/atom_ref/` — that is the path `compute_deltaG.py` reads. `data/atom_ref/` only
> holds the raw single-atom ORCA runs when you regenerate it.

**Additionally, the following CSV files will be generated in `data/`:**

| File | Description | Generated By |
|---|---|---|
| `input_molecules.csv` | **Input file** — columns `PID`, `smiles`. Override the path with `INPUT_CSV=...` | **User-provided**, usually derived from `seeds/` |
| `deltaG_raw.csv` | Extracted Gibbs free energies **plus per-molecule diagnostics** | `scripts/extract_thermo.py` (`mol`, `Gibbs_Eh`, `G_minus_Eel`, `n_imag_saddle`, `terminated_normally`, `level`, `usable`, ...) |
| `failed_monomers.csv` | SMILES that failed 3D conversion, with the reason | `scripts/generate_xyz.py` (only if failures occur) |

---

## Prerequisites

### Software

| Software | Version | Purpose |
|---|---|---|
| [ORCA](https://orcaforum.kofo.mpg.de/) | 6.1.0 | DFT calculations (geometry optimization + frequency) |
| Python | 3.8+ | Pipeline scripts |
| SGE (Sun Grid Engine) | — | Job scheduler (for HPC clusters) |

### Python Packages

The DFT half (`scripts/`) needs:

```bash
conda create -n free_energy python=3.10 pandas rdkit numpy scipy scikit-learn matplotlib
conda activate free_energy
```

The learning half (`ml/`) additionally needs PyTorch and PyTorch Geometric. A GPU is not
required for selection — `seed_select.py` and `active_learning_v3.py` fall back to CPU —
but the Tanimoto similarity that dominates both is a dense matmul, so a GPU is worth
having: scoring 80,000 candidates against the whole 2.06M pool takes **138 s on an
NVIDIA A10 versus ~3.8 h on CPU**. Training the GIN ensemble effectively requires one.

```bash
conda install pytorch pytorch-cuda -c pytorch -c nvidia
pip install torch-geometric
```

> Pin your RDKit version if you intend to rebuild the pool. Canonical SMILES and
> aromaticity perception can change between releases, which would silently shift the
> canonical forms and the fingerprints. `pool/README.md` records the versions this pool
> was built with (RDKit 2025.09.2).

---

## Step-by-Step Usage

### 0. Prepare Input Data

Create `data/input_molecules.csv` with two columns, `PID` and `smiles`:

```csv
PID,smiles
SD1,*CC*
SD2,*C(C*)C
SD3,*C(C*)CC
...
```

Any other columns are ignored, so a seed file can be trimmed down directly:

```bash
python -c "import pandas as pd; d=pd.read_csv('seeds/seed_v9.csv'); \
           d['PID']='SD'+d['rank'].astype(str); \
           d[['PID','smiles']].to_csv('data/input_molecules.csv', index=False)"
```

> **Note:** The `*` symbols represent polymer repeat-unit attachment points and are
> replaced with `C` (methyl caps) during processing, so `*CC*` is computed as `CCC`.
> Capping is consistent across the whole dataset, so the caps' contribution is a
> constant that the composition baseline absorbs — but for very small repeat units the
> caps dominate the molecule, which is why the pool applies a 14-atom lower bound.

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
functional, basis, grid and integral approximation as the molecules**, otherwise ΔG
is meaningless — the atoms are subtracted from the molecular Gibbs energy and nothing
cancels. That is why the keyword line is not written in any generator: every script
imports it from `scripts/orca_settings.py`, and `parse_atom_ref.py` re-reads the level
of theory echoed in each ORCA output and **refuses to write `atom_ref.csv` on a
mismatch**. (This is not hypothetical: the atoms previously ran at `def2-TZVPPD`
while the molecules ran at `def2-TZVP`.)

They are **Gibbs free energies**, not electronic energies, because the molecular side
of the subtraction is a Gibbs free energy. A free atom has no vibrations and no
rotations, so its *G* is just *E*(elec) plus the translational term and the electronic
degeneracy — small, but not zero and not constant across elements.

To regenerate them (e.g. to add an element, or to switch functional):

```bash
# 1. Write ORCA inputs for every atom (ground-state multiplicities are built in:
#    H=2, C=3, N=4, O=3, F=2, Si=3, P=4, S=3, Cl=2, Ge=3, Br=2, Sn=3, I=2)
python scripts/generate_atom_ref_inp.py

# 2. Run them all (serial -- a single atom has ~10-30 basis functions, so MPI
#    parallelism is counter-productive; the inputs deliberately contain no %pal)
qsub submit_atom_ref.csh

# 3. Parse -> data/atom_ref/atom_ref_production.csv, then write the production file
python scripts/parse_atom_ref.py --write-production --use G
```

`parse_atom_ref.py` prints an old-vs-new comparison, checks SCF convergence and
multiplicity, verifies the level of theory against `orca_settings.py`, and refuses to
write if anything is wrong. `--use G` (Gibbs) is the default and the correct choice;
`--use E` exists only to reproduce the historical pre-2026-09 ΔG values.

> The old pipeline subtracted **electronic** atom energies from **Gibbs** molecular
> energies, so its ΔG was not the free energy of any process. It survived in practice
> only because the mismatch is a per-element constant × atom count, which the
> composition-level Ridge baseline absorbs exactly (R² = 1.0000000000 against the
> correction term) — the residual target, i.e. the actual training signal, was
> unaffected. It is fixed now regardless.

### 1. Generate 3D XYZ Coordinates

```bash
N_WORKERS=16 python scripts/generate_xyz.py
```

- **Input:** `data/input_molecules.csv` (or `INPUT_CSV=...`)
- **Output:** `data/xyz/{PID}.xyz`; failures land in `data/failed_monomers.csv`
- Parallel, with `chunksize=1`: per-molecule cost spans ~0.02 s to ~50 s, so work is
  handed out one at a time rather than pre-sliced. Output is byte-identical to a serial
  run. 3,915 molecules take ~25 min on 16 workers.

**This step decides the answer, not just the starting point.** ORCA's geometry
optimisation is a *local* minimiser: it relaxes into the nearest torsional minimum and
never crosses a barrier, so whichever conformer it is handed determines the final
energy. Measured on four molecules that were submitted twice under different SMILES
spellings of the same structure:

| starting conformer | final Gibbs differed by |
|---|---|
| different | 0.78 and 3.14 kcal/mol |
| identical | 0.000 and 0.015 kcal/mol |

Four out of four, no exceptions. So this script:

1. **canonicalises the SMILES before embedding.** ETKDG is deterministic for a fixed
   seed only for a fixed *atom order*, so the same molecule written two ways embedded
   differently and got two different ΔG values. That was the actual cause of the
   discrepancies above.
2. **generates many conformers and keeps the lowest UFF energy**, scaling the count with
   rotatable bonds (20 / 50 / 100). A fixed budget is statistically thin for floppy
   molecules: a 34-rotatable-bond PEG has ~3³⁴ torsional minima.

Against ~15 h of DFT per molecule this costs ~12 s, and across 100 pool molecules it
improved the starting UFF energy for 96, tied 4, and regressed none — a median gain of
4.7 kcal/mol, rising to 8.2 for molecules with 13+ rotatable bonds. Set `N_CONFS=1` to
restore the old single-conformer behaviour.

### 2. Generate Step 1 Optimization Inputs (PBE/def2-SVP)

```bash
python scripts/generate_step1_opt_inp.py
```

- **Input:** `data/xyz/*.xyz`
- **Output:** `data/opt_inp/{PID}_step1_pbe_opt.inp`
- ORCA keywords: `! PBE def2-SVP TightSCF Opt RIJCOSX def2/J`
- Cores: `%pal nprocs` is filled in from `$NSLOTS`, so it always matches the
  `-pe smp` the scheduler actually granted (ORCA ignores the request otherwise)
- Use `--mol PID` to generate one molecule's input instead of the whole directory

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
2. **Step 2 — wB97X-D3/def2-TZVP optimization** (accurate geometry, `TightOpt`)
3. **Frequency calculation** at 298.15 K, same functional/basis/grid (`VeryTightSCF`)
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
| 2 | ωB97X-D3 | def2-TZVP | D3 (built in) | def2/J (RIJCOSX, DefGrid3) | Accurate final geometry (`TightOpt`) |

### Frequency Calculation

| Property | Value |
|---|---|
| Method | ωB97X-D3/def2-TZVP |
| Temperature | 298.15 K |
| Aux Basis | def2/J (RIJCOSX) |
| Integration grid | DefGrid3 |
| SCF Convergence | VeryTightSCF |

Chosen by measurement, not by habit — see
[`docs/orca_benchmark_2026-09-22.md`](docs/orca_benchmark_2026-09-22.md). In short:
RIJCOSX is the *only* option (ORCA refuses an analytic Hessian with RIJK, and exact
4-centre integrals cost 17.7× at 29 atoms), it costs 0.13 kcal/mol in Gibbs, and
DefGrid3 is what removes spurious imaginary frequencies.

### ΔG Calculation

The atomization free energy is computed as:

```
ΔG(molecule) = G(molecule) − Σ nᵢ × G(atomᵢ)
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
# All three are environment variables -- no need to edit the script:
#   INPUT_CSV  input list, molecule id in the first column (default data/input_molecules.csv)
#   START      first line to submit; 2 skips the header. Raise it to resume a batch.
#   NUMBER     cap on jobs submitted (default: no cap)

INPUT_CSV=data/input_molecules.csv START=2 NUMBER=4000 ./full_submit_free_energy.sh
```

Each molecule runs as an independent job requesting:
- 8 SMP cores
- Long queue

---

## License

This project is part of ongoing research. Please contact the authors before use.
