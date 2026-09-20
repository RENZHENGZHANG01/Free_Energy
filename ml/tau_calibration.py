#!/usr/bin/env python3
"""
Calibrate the learnability radius tau from REAL GIN CV results.

Question this answers (for AL v2 design): "Does having ONE labeled neighbor
within Tanimoto tau make a molecule predictable by the GNN?" If yes, what tau?

For each labeled point we already have the GNN's out-of-fold |error| (from
gin_cv_plots/cv_predictions.csv). We compute that point's nearest-neighbor
Tanimoto to the REST of the labeled set (excluding its own canonical-SMILES
group, since GroupKFold groups by canonical SMILES -> no leakage). Then we bin
|error| by that NN-Tanimoto. If error falls sharply once NN >= tau, that tau is
the radius at which a single labeled neighbor confers learnability.
"""
import os, warnings, numpy as np, pandas as pd
from rdkit import Chem, RDLogger
from rdkit.Chem import AllChem, DataStructs
warnings.filterwarnings("ignore"); RDLogger.DisableLog("rdApp.*")

DIR = os.path.dirname(os.path.abspath(__file__))
NBITS = 2048

def canon(smi):
    m = Chem.MolFromSmiles(str(smi))
    return Chem.MolToSmiles(m) if m else None

def to_fp(smi):
    m = Chem.MolFromSmiles(str(smi))
    if m is None: return None
    a = np.zeros(NBITS, np.float32)
    DataStructs.ConvertToNumpyArray(AllChem.GetMorganFingerprintAsBitVect(m, 2, nBits=NBITS), a)
    return a

def tani(A, B):
    inter = A @ B.T; a = A.sum(1)[:, None]; b = B.sum(1)[None, :]
    return inter / (a + b - inter + 1e-9)

def main():
    cv = pd.read_csv(os.path.join(DIR, "gin_cv_plots/cv_predictions.csv"))
    meta = pd.read_csv(os.path.join(DIR, "final_data_with_residual_deltaG.csv"),
                       usecols=["mol", "smiles_clean"])
    df = cv.merge(meta, on="mol", how="left").dropna(subset=["smiles_clean"]).reset_index(drop=True)

    fps, grp, keep = [], [], []
    for s in df["smiles_clean"]:
        f = to_fp(s); c = canon(s)
        keep.append(f is not None and c is not None)
        fps.append(f if f is not None else np.zeros(NBITS, np.float32))
        grp.append(c)
    df["grp"] = grp; df = df[keep].reset_index(drop=True)
    fp = np.array([fps[i] for i in range(len(keep)) if keep[i]])
    grp = df["grp"].values

    # nearest-neighbor Tanimoto to the rest of the labeled set, excluding same group
    S = tani(fp, fp)
    for i in range(len(S)):
        S[i, grp == grp[i]] = -1.0          # mask self + same canonical SMILES (GroupKFold leakage)
    nn = S.max(1)
    df["nn_tani"] = nn

    print("=" * 72)
    print("  tau calibration: GNN out-of-fold |error| vs nearest train-neighbor")
    print("=" * 72)
    print(f"  labeled points: {len(df)}  | new={int(df['is_new'].sum())}  old={int((~df['is_new']).sum())}")
    print(f"  residual target std: {df['y_true'].std():.2f} kcal/mol\n")

    bins = [(-0.01, 0.2), (0.2, 0.3), (0.3, 0.4), (0.4, 0.5), (0.5, 0.6), (0.6, 0.8), (0.8, 1.01)]
    print(f"  {'NN-Tanimoto bin':>16s} {'n':>5s} {'%new':>5s} {'median|err|':>12s} {'mean|err|':>10s} {'P90|err|':>9s}")
    for lo, hi in bins:
        m = (nn > lo) & (nn <= hi)
        if m.sum() == 0:
            print(f"  ({lo:.2f},{hi:.2f}]{'':>4s} {0:5d}"); continue
        e = df.loc[m, "abs_error"].values
        pn = 100 * df.loc[m, "is_new"].mean()
        print(f"  ({lo:.2f},{hi:.2f}]{'':>4s} {m.sum():5d} {pn:5.0f} "
              f"{np.median(e):12.2f} {e.mean():10.2f} {np.percentile(e,90):9.2f}")

    # correlation + a simple "learnable" threshold scan
    from scipy.stats import spearmanr
    rho = spearmanr(nn, df["abs_error"]).correlation
    print(f"\n  Spearman(NN-Tanimoto, |error|) = {rho:.3f}  (negative => closer neighbor, lower error)")
    print(f"\n  {'tau':>5s} {'frac labeled w/ NN>=tau':>24s} {'median|err| above':>18s} {'median|err| below':>18s}")
    for tau in [0.3, 0.35, 0.4, 0.45, 0.5, 0.55, 0.6]:
        ab = df.loc[nn >= tau, "abs_error"]; be = df.loc[nn < tau, "abs_error"]
        print(f"  {tau:5.2f} {100*np.mean(nn>=tau):23.1f}% "
              f"{(np.median(ab) if len(ab) else np.nan):18.2f} {(np.median(be) if len(be) else np.nan):18.2f}")
    print("\n  DONE!")

if __name__ == "__main__":
    main()
