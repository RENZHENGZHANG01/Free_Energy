"""
Size-Independent ΔG via FS5 + Ridge Regression
==============================================
Produces `Delta_G_residual = ΔG − Ridge_FS5(composition)`: a size-independent
thermodynamic-stability signal used to guide the diffusion / generative model.
Raw formation free energy scales ~linearly with atom count, so an un-referenced
ΔG would merely push the generator toward smaller molecules — removing the
composition baseline is mandatory, not an accuracy trick.

FS5 — the ONLY feature set (adopted 2026-06-25; earlier FS1/FS2 removed)
------------------------------------------------------------------------
  (1) N_<element>            element counts, 13 elements incl. Ge/Sn   →  13
  (2) N_{single,double,triple,aromatic}   bond-type totals             →   4
      BP_<el>-<el>_<bondtype>  bond-element-pair × bond-type counts    →  ~34 (dynamic)
  (3) N_rings, N_aromatic_rings           ring totals                  →   2
      RING_<size>_<arom|aliph> ring counts by size & aromaticity       →   ~7 (dynamic)
                                                                        ≈ 60 columns
BP_/RING_ keys are DYNAMIC (absent pair/ring ⇒ 0); ring sizes > 8 bucket to
"big" so macrocycles are never out-of-vocabulary.

Why exactly this — and why nothing richer
-----------------------------------------
Formation free energy is ~additive over BONDS (not atoms), so bond-resolved
counts are the physically-correct dominant descriptor; ring size captures strain
(strained 3-ring vs stable 6-ring). Everything stays COMPOSITION-level — never
topological fingerprints (Morgan) — because a baseline that encodes structure
would absorb the stability signal and drive the residual to 0, killing the
guidance signal. Two constraints must hold simultaneously: strong enough to
remove size dependence, weak enough to leave the structural signal.

Validated alternatives, all rejected (out-of-fold GroupKFold new-subset MAE,
kcal/mol; sources `baseline_comparison.py` + a 15× repeated group-split sweep,
both deleted 2026-07-28 once FS5 was final):
  Ridge + FS1 (element+bond-type+ring totals, the original)  54.6
  Ridge + FS2 (element × hybridization)                      45.9
  Ridge + FS4 (bond-additivity, element-pair × bondtype)     ~24
  Ridge + FS5 (FS4 + ring size × aromatic)          20.5 ± 2.4  ← ADOPTED
  Ridge + FS6 (FS5 + atom environments)              22.6 ± 4.1  overfits
  HistGradientBoosting / RandomForest + FS1         244 / 372  ← see below
  LASSO / ElasticNet                                 worse than Ridge
NONLINEAR baselines are decisively worse AND re-introduce size dependence
(r(size) ≈ −0.2 vs FS5's −0.001), because ΔG is a linear function of composition
counts (linear Ridge reaches ΔG-R² = 0.9999) and trees cannot extrapolate a
linear function to out-of-range compositions. Do not switch to nonlinear, and
do not predict raw ΔG directly.

Known limit: any count baseline has a fixed vocabulary, so a generated molecule
containing an unseen bond/ring/element silently drops that contribution. Not
fixable by featurization — only by data coverage. Flag OOV molecules downstream.

Input / output
--------------
Input  `Delta_G` is in HARTREE and is converted up front; the Ridge fit,
predictions, and residual are therefore all in KCAL/MOL. (Anything consuming
this residual must use kcal/mol — a stale Hartree-era scale constant in
`predict_dataset_gin.py` once caused a silent 1.77× error.)
Net-charged species are excluded (the DFT pipeline runs everything as a neutral
singlet, so their ΔG is invalid): rows are kept with NaN residual so downstream
`dropna()` drops them naturally.

Current state: n=1077, residual mean −0.0000, sd 8.8706 kcal/mol,
r(residual, total_atoms) = −0.0009  ⇒ size-independent, signal intact.
"""

import os
import warnings
import numpy as np
import pandas as pd
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from rdkit import Chem
from rdkit.Chem import Descriptors, rdMolDescriptors, rdchem
from sklearn.linear_model import Ridge
from sklearn.model_selection import cross_val_score, train_test_split
from sklearn.preprocessing import StandardScaler
from sklearn.metrics import r2_score, mean_absolute_error
from scipy import stats

warnings.filterwarnings("ignore", category=DeprecationWarning)

# ─── paths ───────────────────────────────────────────────────
# This file is Step 6 of the ORCA pipeline and lives with it, but the ML workspace
# one level up runs it on a different (larger) dataset. Defaults therefore point at
# files NEXT TO THIS SCRIPT (standalone pipeline use); override via env vars to run
# it against another dataset without copying the script:
#
#   RESIDUAL_INPUT=/groups/tluo/FFE_Renzheng/final_data_with_deltaG.csv \
#   RESIDUAL_OUTPUT=/groups/tluo/FFE_Renzheng/final_data_with_residual_deltaG.csv \
#   python Renzheng/scripts/compute_residual_deltaG.py
#
# Plots default to a residual_plots/ next to the OUTPUT, so they follow the dataset.
DIR = os.path.dirname(os.path.abspath(__file__))
INPUT   = os.environ.get("RESIDUAL_INPUT",  os.path.join(DIR, "final_data_with_deltaG.csv"))
OUTPUT  = os.environ.get("RESIDUAL_OUTPUT", os.path.join(DIR, "final_data_with_residual_deltaG.csv"))
PLOTDIR = os.environ.get("RESIDUAL_PLOTDIR",
                         os.path.join(os.path.dirname(os.path.abspath(OUTPUT)), "residual_plots"))
os.makedirs(PLOTDIR, exist_ok=True)

# The DFT Delta_G in the input CSV is in Hartree. We convert it to kcal/mol up
# front so the Ridge baseline, predictions and residual are all in kcal/mol.
HARTREE_TO_KCAL = 627.509474

# ─── plotting style ──────────────────────────────────────────
plt.rcParams.update({
    "font.family": "DejaVu Sans", "font.size": 12,
    "axes.labelsize": 14, "axes.titlesize": 16,
    "xtick.labelsize": 11, "ytick.labelsize": 11,
    "figure.dpi": 150,
})

# ─── feature extraction ──────────────────────────────────────
# Element vocabulary for the count-based baseline. Must cover every element in
# the dataset, otherwise out-of-vocab atoms are silently dropped from the feature
# vector and the baseline (hence the residual) is wrong for those molecules.
ATOM_TYPES = ["C", "H", "O", "N", "F", "Cl", "Br", "S", "I", "P", "Si", "Ge", "Sn"]


def extract_features(smiles_clean):
    """Extract the FS5 feature vector from a cleaned SMILES (see module docstring)."""
    if smiles_clean is None or (isinstance(smiles_clean, float) and np.isnan(smiles_clean)):
        return None

    mol = Chem.MolFromSmiles(str(smiles_clean))
    if mol is None:
        return None

    mol_H = Chem.AddHs(mol)

    # FS5 (1): element counts
    atom_counts = {a: 0 for a in ATOM_TYPES}
    for atom in mol_H.GetAtoms():
        sym = atom.GetSymbol()
        if sym in atom_counts:
            atom_counts[sym] += 1

    # FS5 (2): bond-type totals + bond-element-pair × bond-type counts (the FS5 core:
    # formation free energy is ~additive over BONDS, so bond-resolved counts dominate).
    bond_counts = {"single": 0, "double": 0, "triple": 0, "aromatic": 0}
    bond_pairs = {}
    for bond in mol_H.GetBonds():
        bt = bond.GetBondType()
        if bt == Chem.rdchem.BondType.SINGLE:
            bond_counts["single"] += 1
        elif bt == Chem.rdchem.BondType.DOUBLE:
            bond_counts["double"] += 1
        elif bt == Chem.rdchem.BondType.TRIPLE:
            bond_counts["triple"] += 1
        elif bt == Chem.rdchem.BondType.AROMATIC:
            bond_counts["aromatic"] += 1
        e1, e2 = sorted([bond.GetBeginAtom().GetSymbol(), bond.GetEndAtom().GetSymbol()])
        bp_key = f"BP_{e1}-{e2}_{str(bt)}"
        bond_pairs[bp_key] = bond_pairs.get(bp_key, 0) + 1

    # FS5 (3): ring topology — totals + counts resolved by (size, aromatic), which
    # captures ring strain (a strained 3-ring vs a stable 6-ring). Sizes > 8 are
    # bucketed as "big" so arbitrarily large macrocycles never become OOV.
    ri = mol.GetRingInfo()
    n_rings = ri.NumRings()
    n_aromatic_rings = rdMolDescriptors.CalcNumAromaticRings(mol)
    ring_size_counts = {}
    for ring in ri.AtomRings():
        sz = len(ring) if len(ring) <= 8 else "big"
        is_arom = all(mol.GetAtomWithIdx(i).GetIsAromatic() for i in ring)
        rk = f"RING_{sz}_{'arom' if is_arom else 'aliph'}"
        ring_size_counts[rk] = ring_size_counts.get(rk, 0) + 1

    # Assemble, DROPPING the columns that are exact linear combinations of others.
    # Keeping them made the design matrix rank-deficient by construction (measured on
    # 3957 molecules: 107 columns, rank 95). That is harmless for prediction -- Ridge
    # simply splits the coefficient between the redundant columns -- but it costs two
    # things that matter here:
    #   * coefficients stop being identifiable. Adding d to N_single while subtracting
    #     d from all 49 BP_*_SINGLE columns gives identical predictions, so "the energy
    #     of a C-C single bond" has no unique answer.
    #   * any D-optimal / leverage-driven selection degenerates: det(X'X) is identically
    #     zero, so the criterion is driven purely by the ridge term and goes chasing
    #     null-space directions that no amount of data can ever pin down. That was the
    #     real reason the seed selector kept reaching for exotic molecules.
    #
    # Dropped, each with the identity that makes it redundant (all verified to exactly
    # 0.000000 across 3957 molecules):
    #   N_single / N_double / N_triple / N_aromatic = sum of the BP_*_<bondtype> columns
    #   N_rings                                     = sum of the RING_* columns
    #   N_H / N_F / N_Cl / N_Br / N_I               = number of bonds containing that
    #       atom, because hydrogen and the halogens are monovalent: each atom
    #       contributes exactly one bond. (Hypervalent iodine breaks this identity,
    #       which is one more reason to keep such structures out of the pool.)
    #
    # The TOTALS are dropped rather than the detail columns: keeping BP_C-H, BP_N-H and
    # BP_O-H lets each bond type carry its own energy, whereas keeping only N_H would
    # collapse them into a single meaningless "hydrogen count".
    #
    # N_aromatic_rings is KEPT. It is NOT redundant with the RING_*_arom columns --
    # RDKit's aromatic-ring perception and the all-atoms-aromatic test used above
    # disagree (measured maximum difference: 3 rings).
    _REDUNDANT = {"N_single", "N_double", "N_triple", "N_aromatic", "N_rings",
                  "N_H", "N_F", "N_Cl", "N_Br", "N_I"}

    feat = {}
    feat.update({f"N_{k}": v for k, v in atom_counts.items()})
    feat.update({f"N_{k}": v for k, v in bond_counts.items()})
    feat["N_rings"] = n_rings
    feat["N_aromatic_rings"] = n_aromatic_rings
    feat.update(bond_pairs)
    feat.update(ring_size_counts)
    for _k in _REDUNDANT:
        feat.pop(_k, None)

    return feat


# ─── helper: parity plot ─────────────────────────────────────
def save_scatter(x, y_vals, xlabel, ylabel, title, fname, color):
    fig, ax = plt.subplots(figsize=(9, 6))
    ax.scatter(x, y_vals, s=18, alpha=0.55, c=color, edgecolors="white", linewidths=0.3)
    mask = np.isfinite(x) & np.isfinite(y_vals)
    if mask.sum() > 2:
        z = np.polyfit(x[mask], y_vals[mask], 1)
        p = np.poly1d(z)
        xf = np.linspace(x[mask].min(), x[mask].max(), 200)
        ax.plot(xf, p(xf), "--", color="#1e293b", lw=1.5, alpha=0.7,
                label=f"y = {z[0]:.6f}x + {z[1]:.4f}")
        ss_r = np.sum((y_vals[mask] - p(x[mask])) ** 2)
        ss_t = np.sum((y_vals[mask] - np.mean(y_vals[mask])) ** 2)
        r2 = 1 - ss_r / ss_t if ss_t > 0 else 0
        ax.legend(title=f"R² = {r2:.4f}", loc="best", fontsize=10,
                  title_fontsize=11, framealpha=0.9, edgecolor="#cbd5e1")
    ax.set_xlabel(xlabel, fontweight="bold")
    ax.set_ylabel(ylabel, fontweight="bold")
    ax.set_title(title, fontweight="bold", pad=12)
    ax.grid(True, alpha=0.2)
    ax.set_axisbelow(True)
    for sp in ax.spines.values():
        sp.set_color("#94a3b8"); sp.set_linewidth(0.8)
    fig.tight_layout()
    fig.savefig(fname, bbox_inches="tight")
    plt.close(fig)
    print(f"  📊 {fname}")


# ─── main ─────────────────────────────────────────────────────
def main():
    print("=" * 60)
    print("  Size-Independent ΔG  —  FS5 + Ridge Regression")
    print("  (with proper train/test split)")
    print("=" * 60)

    df = pd.read_csv(INPUT)
    print(f"\nLoaded {len(df)} rows from {INPUT}")

    # keep only rows with Delta_G
    mask_valid = df["Delta_G"].notna()
    print(f"Rows with Delta_G: {mask_valid.sum()}")

    # Exclude net-charged species. The DFT pipeline runs every molecule as a
    # neutral singlet (charge=0, mult=1), so any net-charged SMILES has an
    # invalid ΔG and a wildly wrong residual. Drop them from the baseline /
    # training target (rows are kept in the CSV but with NaN residual, so the
    # GNN scripts' dropna() naturally excludes them).
    def _net_charge(smi):
        m = Chem.MolFromSmiles(str(smi))
        return None if m is None else Chem.GetFormalCharge(m)
    charges = df["smiles_clean"].apply(_net_charge)
    mask_neutral = (charges == 0)
    n_charged = int((mask_valid & ~mask_neutral).sum())
    print(f"Excluding {n_charged} net-charged / unparseable molecules (DFT assumes neutral)")
    mask_valid = mask_valid & mask_neutral
    print(f"Rows kept (neutral, with Delta_G): {mask_valid.sum()}")

    # ── extract features ──
    print("\nExtracting FS5 features (elements + bond-pairs + ring size/aromaticity)...")
    feat_list = []
    for i, smi in enumerate(df["smiles_clean"]):
        if i % 100 == 0:
            print(f"  → row {i}/{len(df)}")
        feat_list.append(extract_features(smi))

    feat_df = pd.DataFrame(feat_list)
    # Bond-pair / ring-size keys are dynamic (a molecule lacking a given pair/ring
    # simply has 0 of it) → fill missing keys with 0. Rows from a failed parse are
    # all-NaN; keep them NaN so the finite-feature check below drops them.
    failed_rows = feat_df.isna().all(axis=1)
    feat_df = feat_df.fillna(0.0)
    feat_df.loc[failed_rows, :] = np.nan
    feature_names = list(feat_df.columns)
    print(f"\nFeature columns ({len(feature_names)}): "
          f"{sum(c.startswith('BP_') for c in feature_names)} bond-pair, "
          f"{sum(c.startswith('RING_') for c in feature_names)} ring-size, "
          f"{sum(c.startswith('N_') for c in feature_names)} atom/bond/ring count")

    # attach features to main df
    for col in feature_names:
        df[col] = feat_df[col].values

    # ── prepare regression data ──
    df_reg = df[mask_valid].copy()
    X_all = df_reg[feature_names].values.astype(float)
    # Convert DFT Delta_G from Hartree to kcal/mol; everything downstream
    # (Ridge fit, predictions, residual, MAE, plots) is therefore in kcal/mol.
    y_all = df_reg["Delta_G"].values.astype(float) * HARTREE_TO_KCAL

    # check for NaN in features
    valid_feat = np.all(np.isfinite(X_all), axis=1) & np.isfinite(y_all)
    X_all = X_all[valid_feat]
    y_all = y_all[valid_feat]
    df_reg = df_reg[valid_feat].copy()
    print(f"\nTotal regression samples (after NaN removal): {len(y_all)}")

    # ══════════════════════════════════════════════════════════
    #  STEP 1: Train/Test Split (85/15)
    # ══════════════════════════════════════════════════════════
    X_train, X_test, y_train, y_test, idx_train, idx_test = train_test_split(
        X_all, y_all, np.arange(len(y_all)),
        test_size=0.15, random_state=42
    )
    print(f"\n  Train set: {len(y_train)}")
    print(f"  Test set:  {len(y_test)}")

    # ── standardize (fit on train only) ──
    scaler = StandardScaler()
    X_train_scaled = scaler.fit_transform(X_train)
    X_test_scaled = scaler.transform(X_test)

    # ══════════════════════════════════════════════════════════
    #  STEP 2: Ridge CV on TRAIN SET only → select best alpha
    # ══════════════════════════════════════════════════════════
    alphas = [0.01, 0.1, 1.0, 10.0, 100.0]
    best_alpha, best_r2 = None, -np.inf
    print("\nRidge CV (5-fold on TRAIN set only):")
    for alpha in alphas:
        model = Ridge(alpha=alpha)
        scores = cross_val_score(model, X_train_scaled, y_train, cv=5, scoring="r2")
        mean_r2 = scores.mean()
        print(f"  α = {alpha:>6.2f}  →  R² = {mean_r2:.6f} ± {scores.std():.6f}")
        if mean_r2 > best_r2:
            best_r2, best_alpha = mean_r2, alpha

    print(f"\n  ✅ Best α = {best_alpha}, CV R² (train) = {best_r2:.6f}")

    # ══════════════════════════════════════════════════════════
    #  STEP 3: Fit on TRAIN, evaluate on TEST
    # ══════════════════════════════════════════════════════════
    model_eval = Ridge(alpha=best_alpha)
    model_eval.fit(X_train_scaled, y_train)

    y_train_pred = model_eval.predict(X_train_scaled)
    y_test_pred = model_eval.predict(X_test_scaled)

    train_r2 = r2_score(y_train, y_train_pred)
    test_r2 = r2_score(y_test, y_test_pred)
    train_mae = mean_absolute_error(y_train, y_train_pred)
    test_mae = mean_absolute_error(y_test, y_test_pred)

    print(f"\n{'=' * 60}")
    print(f"  EVALUATION RESULTS")
    print(f"{'=' * 60}")
    print(f"  Train R²:  {train_r2:.6f}    Train MAE: {train_mae:.6f} kcal/mol")
    print(f"  Test  R²:  {test_r2:.6f}    Test  MAE: {test_mae:.6f} kcal/mol")
    print(f"{'=' * 60}")

    # ── feature importance (from train-only model) ──
    print("\n  Feature coefficients (scaled, from train model):")
    coef_order = np.argsort(np.abs(model_eval.coef_))[::-1]
    for idx in coef_order:
        print(f"    {feature_names[idx]:>20s}  {model_eval.coef_[idx]:>+10.6f}")

    # ── Test set parity plot ──
    print("\nGenerating plots...")
    fig, ax = plt.subplots(figsize=(8, 8))
    # plot train points (faded)
    ax.scatter(y_train, y_train_pred, s=14, alpha=0.25, c="#94a3b8",
               edgecolors="white", linewidths=0.2, label=f"Train (n={len(y_train)})")
    # plot test points (bold)
    ax.scatter(y_test, y_test_pred, s=30, alpha=0.8, c="#DC2626",
               edgecolors="white", linewidths=0.4, label=f"Test (n={len(y_test)})")
    lims = [min(y_all.min(), y_train_pred.min()) - 0.5,
            max(y_all.max(), y_train_pred.max()) + 0.5]
    ax.plot(lims, lims, "--", color="#1e293b", lw=1.5, alpha=0.7, label="Ideal (y=x)")
    ax.set_xlim(lims); ax.set_ylim(lims)
    ax.set_xlabel("DFT ΔG (kcal/mol)", fontweight="bold")
    ax.set_ylabel("Ridge Predicted ΔG (kcal/mol)", fontweight="bold")
    ax.set_title(f"Ridge: Predicted vs Actual ΔG\n"
                 f"Train R² = {train_r2:.4f} | Test R² = {test_r2:.4f} | Test MAE = {test_mae:.4f}",
                 fontweight="bold", pad=12)
    ax.legend(fontsize=11, framealpha=0.9, loc="upper left")
    ax.grid(True, alpha=0.2); ax.set_aspect("equal")
    for sp in ax.spines.values():
        sp.set_color("#94a3b8"); sp.set_linewidth(0.8)
    fig.tight_layout()
    fig.savefig(os.path.join(PLOTDIR, "predicted_vs_actual_test.png"), bbox_inches="tight")
    plt.close(fig)
    print(f"  � {os.path.join(PLOTDIR, 'predicted_vs_actual_test.png')}")

    # ══════════════════════════════════════════════════════════
    #  STEP 4: Refit on ALL data → produce residuals for CSV
    # ══════════════════════════════════════════════════════════
    print("\n--- Refitting on ALL data for final residual CSV ---")
    scaler_all = StandardScaler()
    X_all_scaled = scaler_all.fit_transform(X_all)

    model_final = Ridge(alpha=best_alpha)
    model_final.fit(X_all_scaled, y_all)
    y_all_pred = model_final.predict(X_all_scaled)
    residuals_all = y_all - y_all_pred

    all_r2 = r2_score(y_all, y_all_pred)
    all_mae = mean_absolute_error(y_all, y_all_pred)
    print(f"  All-data R²: {all_r2:.6f} | MAE: {all_mae:.6f}")

    # ── write residuals back ──
    df_reg["Delta_G_predicted"] = y_all_pred
    df_reg["Delta_G_residual"] = residuals_all

    df["Delta_G_predicted"] = np.nan
    df["Delta_G_residual"] = np.nan
    df.loc[df_reg.index, "Delta_G_predicted"] = df_reg["Delta_G_predicted"].values
    df.loc[df_reg.index, "Delta_G_residual"] = df_reg["Delta_G_residual"].values

    df.to_csv(OUTPUT, index=False)
    print(f"\n💾 Saved to: {OUTPUT}")

    # ══════════════════════════════════════════════════════════
    #  Additional validation plots (using all-data residuals)
    # ══════════════════════════════════════════════════════════
    total_atoms = df_reg["total_atoms"].values
    mw = df_reg["MW"].values

    # 1) Residual vs Total Atoms
    save_scatter(total_atoms, residuals_all,
                 "Total Atoms", "ΔG Residual (kcal/mol)",
                 "ΔG Residual vs Total Atoms (should be ≈ 0 correlation)",
                 os.path.join(PLOTDIR, "residual_vs_total_atoms.png"), "#2563EB")

    # 2) Residual vs MW
    save_scatter(mw, residuals_all,
                 "Molecular Weight (g/mol)", "ΔG Residual (kcal/mol)",
                 "ΔG Residual vs Molecular Weight (should be ≈ 0 correlation)",
                 os.path.join(PLOTDIR, "residual_vs_MW.png"), "#7C3AED")

    # 3) All-data predicted vs actual (labeled as "all data")
    save_scatter(y_all_pred, y_all,
                 "Predicted ΔG (kcal/mol)", "Actual ΔG (kcal/mol)",
                 "Predicted vs Actual ΔG (all data, for residual computation)",
                 os.path.join(PLOTDIR, "predicted_vs_actual_alldata.png"), "#059669")

    # 4) Original Delta_G vs total_atoms
    save_scatter(total_atoms, y_all,
                 "Total Atoms", "ΔG (kcal/mol)",
                 "Original ΔG vs Total Atoms (before correction)",
                 os.path.join(PLOTDIR, "original_deltaG_vs_total_atoms.png"), "#DC2626")

    # 5) Residual histogram
    fig, ax = plt.subplots(figsize=(9, 6))
    ax.hist(residuals_all, bins=60, color="#6366f1", edgecolor="white", linewidth=0.5, alpha=0.8)
    mean_r, std_r = residuals_all.mean(), residuals_all.std()
    ax.axvline(mean_r, color="#dc2626", linestyle="--", linewidth=2, label=f"Mean = {mean_r:.4f}")
    ax.set_xlabel("ΔG Residual (kcal/mol)", fontweight="bold")
    ax.set_ylabel("Count", fontweight="bold")
    ax.set_title(f"Distribution of ΔG Residual (σ = {std_r:.4f})", fontweight="bold", pad=12)
    ax.legend(fontsize=11, framealpha=0.9)
    ax.grid(True, alpha=0.2); ax.set_axisbelow(True)
    for sp in ax.spines.values():
        sp.set_color("#94a3b8"); sp.set_linewidth(0.8)
    fig.tight_layout()
    hist_path = os.path.join(PLOTDIR, "residual_histogram.png")
    fig.savefig(hist_path, bbox_inches="tight")
    plt.close(fig)
    print(f"  📊 {hist_path}")

    # ── Pearson correlation check ──
    r_atoms, p_atoms = stats.pearsonr(total_atoms, residuals_all)
    r_mw, p_mw = stats.pearsonr(mw, residuals_all)
    print(f"\n  Pearson correlation (residual vs total_atoms): r = {r_atoms:.6f}, p = {p_atoms:.4e}")
    print(f"  Pearson correlation (residual vs MW):          r = {r_mw:.6f}, p = {p_mw:.4e}")

    print("\n🎉 DONE!")


if __name__ == "__main__":
    main()
