import os
import pandas as pd
import numpy as np
from rdkit import Chem
from rdkit.Chem import AllChem
from sklearn.manifold import TSNE
import matplotlib.pyplot as plt

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



def preprocess_smiles(smi):
    return smi.replace("*", "C") if isinstance(smi, str) else smi


def fp_from_smiles(smiles, nBits=1024):
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return np.zeros(nBits)
    fp = AllChem.GetMorganFingerprintAsBitVect(mol, radius=2, nBits=nBits)
    return np.array(fp)


def generate_tsne_plots():

    # compute_deltaG.py writes the dataset next to itself, in scripts/ -- not in the
    # project root, which is where this used to look.
    here = os.path.dirname(os.path.abspath(__file__))
    root = os.path.dirname(here)

    # -----------------------------------------
    # Load main dataset (with ΔG)
    # -----------------------------------------
    df_main_path = _require(os.environ.get(
        "DELTAG_CSV", os.path.join(here, "final_data_with_deltaG.csv")),
        "the Delta_G dataset")
    print("Loading foreground dataset:", df_main_path)
    df_main = pd.read_csv(df_main_path)
    df_main["smiles_clean"] = df_main["smiles_clean"].apply(preprocess_smiles)

    # -----------------------------------------
    # Load background polymer set
    # -----------------------------------------
    df_bg_path = os.path.join(root, "data", "all_polymer.csv")
    print("Loading background dataset:", df_bg_path)
    df_bg = pd.read_csv(df_bg_path)
    df_bg["smiles_clean"] = df_bg["smiles"].astype(str).apply(preprocess_smiles)

    # -----------------------------------------
    # Combine both sets for joint TSNE
    # -----------------------------------------
    df_main["is_main"] = 1
    df_bg["is_main"]   = 0

    df_all = pd.concat([df_main, df_bg], ignore_index=True)

    print("Computing Morgan fingerprints for ALL polymers...")
    fps = np.array([fp_from_smiles(smi) for smi in df_all["smiles_clean"]])

    # -----------------------------------------
    # t-SNE
    # -----------------------------------------
    print("Running t-SNE...")
    tsne = TSNE(
        n_components=2,
        perplexity=30,
        learning_rate=200,
        metric="cosine",
        random_state=42,
        verbose=1,
    )
    tsne_result = tsne.fit_transform(fps)

    df_all["tsne_x"] = tsne_result[:, 0]
    df_all["tsne_y"] = tsne_result[:, 1]

    df_main_tsne = df_all[df_all["is_main"] == 1].reset_index(drop=True)
    df_bg_tsne   = df_all[df_all["is_main"] == 0].reset_index(drop=True)

    # -----------------------------------------
    # What values to plot for colored points
    # -----------------------------------------
    color_keys = [
        "Delta_G",
        "DeltaG_per_heavy_atom",
        "DeltaG_per_atom",
        "DeltaG_per_backbone_bond",
        "DeltaG_per_bond",
    ]

    out_dir = os.path.join(root, "tsne_plots_with_background")
    os.makedirs(out_dir, exist_ok=True)

    print("Plotting...")

    # -----------------------------------------
    # Make each t-SNE plot
    # -----------------------------------------
    for key in color_keys:

        plt.figure(figsize=(7, 6))

        # background points first (grey)
        plt.scatter(
            df_bg_tsne["tsne_x"], df_bg_tsne["tsne_y"],
            c="lightgrey", s=10, alpha=0.35, edgecolors="none"
        )

        # foreground colored points
        sc = plt.scatter(
            df_main_tsne["tsne_x"], df_main_tsne["tsne_y"],
            c=df_main_tsne[key], cmap="Spectral_r",
            s=22, edgecolors="none"
        )

        plt.colorbar(sc, label=key)
        plt.title(f"t-SNE with background — colored by {key}", fontsize=14)
        plt.xlabel("t-SNE 1")
        plt.ylabel("t-SNE 2")

        save_path = os.path.join(out_dir, f"tsne_with_bg_{key}.png")
        plt.savefig(save_path, dpi=300, bbox_inches='tight')
        plt.close()

        print("Saved:", save_path)

    print("\n🎉 All t-SNE plots saved at:", out_dir)


if __name__ == "__main__":
    generate_tsne_plots()
