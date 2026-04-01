import os
import pandas as pd
import numpy as np
from rdkit import Chem
from rdkit.Chem import AllChem
from sklearn.manifold import TSNE
import matplotlib.pyplot as plt
import seaborn as sns

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

    root = os.path.dirname(os.path.dirname(__file__))

    df_path = os.path.join(root, "final_data_with_deltaG.csv")
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
