import os
import pandas as pd
import numpy as np
from rdkit import Chem
from rdkit.Chem import AllChem
from sklearn.manifold import TSNE
import matplotlib.pyplot as plt
import seaborn as sns

# The DFT dataset is regenerated per campaign and is not in the repo. Point at it
# explicitly rather than failing with a bare FileNotFoundError, because the usual
# reason it is absent is that the level of theory changed and the old one was
# retired -- in which case silently picking up a leftover file would be worse.
def _require(path, what):
    import os as _os, sys as _sys
    if not _os.path.exists(path):
        _sys.exit(
            f"\n  ABORT: {what} not found:\n    {path}\n\n"
            "  It is produced by the DFT pipeline:\n"
            "      python scripts/extract_thermo.py --expect-level\n"
            "      python scripts/compute_deltaG.py\n"
            "      python scripts/compute_residual_deltaG.py\n\n"
            "  Set the matching environment variable to use a different path.")
    return path


# =========================================================
# Compute Morgan Fingerprint
# =========================================================
def fp_from_smiles(smiles, nBits=1024):
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return np.zeros(nBits)
    fp = AllChem.GetMorganFingerprintAsBitVect(mol, radius=2, nBits=nBits)
    return np.array(fp)


# =========================================================
# Main TSNE Function
# =========================================================
def generate_tsne_plots():

    # compute_deltaG.py writes the dataset next to itself, in scripts/ -- not in the
    # project root, which is where this used to look.
    here = os.path.dirname(os.path.abspath(__file__))
    root = os.path.dirname(here)

    df_path = _require(os.environ.get(
        "DELTAG_CSV", os.path.join(here, "final_data_with_deltaG.csv")),
        "the Delta_G dataset")
    print("Loading:", df_path)
    df = pd.read_csv(df_path)

    # =====================================================
    # Compute fingerprints
    # =====================================================
    print("Computing Morgan fingerprints...")
    fps = np.array([fp_from_smiles(smi) for smi in df["smiles_clean"]])

    # =====================================================
    # Run t-SNE
    # =====================================================
    print("Running t-SNE...")
    tsne = TSNE(
        n_components=2,
        perplexity=30,
        learning_rate=200,
        # n_iter=2000,
        metric="cosine",
        random_state=42,
        verbose=1,
    )
    tsne_result = tsne.fit_transform(fps)

    df["tsne_x"] = tsne_result[:, 0]
    df["tsne_y"] = tsne_result[:, 1]

    # =====================================================
    # Values to color with
    # =====================================================
    color_keys = [
        "Delta_G",
        "DeltaG_per_heavy_atom",
        "DeltaG_per_atom",
        "DeltaG_per_backbone_bond",
        "DeltaG_per_bond",
    ]

    out_dir = os.path.join(root, "tsne_plots")
    os.makedirs(out_dir, exist_ok=True)

    # =====================================================
    # Plot each TSNE
    # =====================================================
    for key in color_keys:
        plt.figure(figsize=(7, 6))
        sc = plt.scatter(df["tsne_x"], df["tsne_y"],
                         c=df[key], cmap="Spectral_r",
                         s=22, edgecolor="none")

        plt.colorbar(sc, label=key)
        plt.title(f"t-SNE Colored by {key}", fontsize=14)
        plt.xlabel("t-SNE 1")
        plt.ylabel("t-SNE 2")

        save_path = os.path.join(out_dir, f"tsne_{key}.png")
        plt.savefig(save_path, dpi=300, bbox_inches='tight')
        plt.close()

        print("Saved:", save_path)

    print("\n🎉 All t-SNE plots done! Saved in:", out_dir)


if __name__ == "__main__":
    generate_tsne_plots()
