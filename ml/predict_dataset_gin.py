#!/usr/bin/env python3
"""
Generalized GIN Ensemble Prediction Pipeline (v2 Architecture)
For OMG / omics and any other candidate datasets.
"""

import os
import argparse
import pandas as pd
import numpy as np
import torch
import warnings

import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D
from sklearn.preprocessing import StandardScaler
from sklearn.model_selection import GroupKFold
from sklearn.manifold import TSNE
import json

from torch_geometric.data import Data
from torch_geometric.loader import DataLoader
from rdkit import Chem, RDLogger
from rdkit.Chem import AllChem, rdchem
from torch.nn import Linear, Sequential, BatchNorm1d, ReLU, Dropout
import torch.nn.functional as F
from torch_geometric.nn import GINEConv, global_mean_pool, global_max_pool, global_add_pool

warnings.filterwarnings('ignore')
RDLogger.DisableLog("rdApp.*")

# ----------------- GLOBALS & CONSTANTS -----------------
DIR = os.path.dirname(os.path.abspath(__file__))
MODEL_DIR = os.path.join(DIR, "gin_cv_models")
TRAIN_CSV = os.path.join(DIR, "final_data_with_residual_deltaG.csv")
XYZ_DIR = os.path.join(DIR, "xyz_unified")
HARTREE_TO_KCAL = 627.509474

# NOTE: must match gin_residual_v2_cv.py EXACTLY — the gin_cv_models checkpoints were
# trained with the 13-element vocab (+Ge,Sn) => NODE_DIM=34. An 11-element/32-dim
# featurizer fails to load them (node_encoder shape mismatch 34 vs 32).
COMMON_SYMBOLS = ['C', 'H', 'N', 'O', 'S', 'F', 'Cl', 'P', 'Si', 'Br', 'I', 'Ge', 'Sn']
ELECTRONEGATIVITY = {
    'H': 2.20, 'C': 2.55, 'N': 3.04, 'O': 3.44, 'F': 3.98,
    'P': 2.19, 'S': 2.58, 'Cl': 3.16, 'Si': 1.90, 'Br': 2.96, 'I': 2.66,
    'Ge': 2.01, 'Sn': 1.96, '*': 2.0
}
VDW_RADIUS = {
    'H': 1.20, 'C': 1.70, 'N': 1.55, 'O': 1.52, 'F': 1.47,
    'P': 1.80, 'S': 1.80, 'Cl': 1.75, 'Si': 2.10, 'Br': 1.85, 'I': 1.98,
    'Ge': 2.11, 'Sn': 2.17, '*': 1.5
}
HYBRIDIZATION_TYPES = [
    rdchem.HybridizationType.SP, rdchem.HybridizationType.SP2,
    rdchem.HybridizationType.SP3, rdchem.HybridizationType.SP3D,
    rdchem.HybridizationType.SP3D2,
]
RING_SIZES = [3, 4, 5, 6, 7, 8]
BOND_STEREO_TYPES = [
    rdchem.BondStereo.STEREONONE, rdchem.BondStereo.STEREOZ,
    rdchem.BondStereo.STEREOE, rdchem.BondStereo.STEREOANY,
]
NODE_DIM = 34
EDGE_DIM = 10

# ----------------- FEATURIZATION -----------------
def get_atom_features(atom):
    symbol = atom.GetSymbol()
    sym_feat = [1.0 if symbol == s else 0.0 for s in COMMON_SYMBOLS]
    sym_feat.append(1.0 if symbol not in COMMON_SYMBOLS else 0.0)
    en = ELECTRONEGATIVITY.get(symbol, 2.0)
    vdw = VDW_RADIUS.get(symbol, 1.5)
    mass = atom.GetMass() * 0.01
    hyb = atom.GetHybridization()
    hyb_feat = [1.0 if hyb == h else 0.0 for h in HYBRIDIZATION_TYPES]
    features = (
        sym_feat + [en, vdw, mass] + hyb_feat
        + [float(atom.GetDegree()), float(atom.GetFormalCharge()),
           float(atom.GetIsAromatic()), float(atom.GetTotalNumHs()),
           float(atom.GetNumRadicalElectrons()), float(atom.IsInRing())]
        + [float(atom.IsInRingSize(s)) for s in RING_SIZES]
    )
    return np.array(features, dtype=np.float32)

def get_bond_features(bond):
    bt = bond.GetBondType()
    bond_type = [float(bt == rdchem.BondType.SINGLE), float(bt == rdchem.BondType.DOUBLE),
                 float(bt == rdchem.BondType.TRIPLE), float(bt == rdchem.BondType.AROMATIC)]
    stereo = bond.GetStereo()
    stereo_feat = [float(stereo == s) for s in BOND_STEREO_TYPES]
    return np.array(bond_type + [float(bond.GetIsConjugated()), float(bond.IsInRing())]
                    + stereo_feat, dtype=np.float32)

def smiles_to_graph(smiles):
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return None
    atom_feats = [get_atom_features(atom) for atom in mol.GetAtoms()]
    x = torch.tensor(np.array(atom_feats), dtype=torch.float)
    edge_list, edge_feat_list = [], []
    for bond in mol.GetBonds():
        i, j = bond.GetBeginAtomIdx(), bond.GetEndAtomIdx()
        bf = get_bond_features(bond)
        edge_list.extend([[i, j], [j, i]])
        edge_feat_list.extend([bf, bf])
    if edge_list:
        edge_index = torch.tensor(edge_list, dtype=torch.long).t().contiguous()
        edge_attr = torch.tensor(np.array(edge_feat_list), dtype=torch.float)
    else:
        edge_index = torch.empty((2, 0), dtype=torch.long)
        edge_attr = torch.empty((0, EDGE_DIM), dtype=torch.float)
    return Data(x=x, edge_index=edge_index, edge_attr=edge_attr,
                y=torch.tensor([0.0], dtype=torch.float).view(1, -1))

def get_morgan_fingerprint(smiles, radius=2, nbits=2048):
    mol = Chem.MolFromSmiles(smiles)
    if mol is None: return np.zeros(nbits)
    return np.array(AllChem.GetMorganFingerprintAsBitVect(mol, radius, nBits=nbits), dtype=np.float32)

# ----------------- ARCHITECTURE -----------------
class GINEResidualModel(torch.nn.Module):
    def __init__(self, node_dim, edge_dim, hidden_dim=128, num_layers=4, dropout=0.15):
        super().__init__()
        self.num_layers = num_layers
        self.dropout = dropout
        self.node_encoder = Sequential(Linear(node_dim, hidden_dim), ReLU())
        self.convs = torch.nn.ModuleList()
        self.batch_norms = torch.nn.ModuleList()
        for _ in range(num_layers):
            mlp = Sequential(Linear(hidden_dim, hidden_dim * 2),
                             BatchNorm1d(hidden_dim * 2), ReLU(),
                             Linear(hidden_dim * 2, hidden_dim))
            self.convs.append(GINEConv(mlp, edge_dim=edge_dim))
            self.batch_norms.append(BatchNorm1d(hidden_dim))
        self.predictor = Sequential(
            Linear(hidden_dim * 3, hidden_dim * 2), BatchNorm1d(hidden_dim * 2),
            ReLU(), Dropout(dropout),
            Linear(hidden_dim * 2, hidden_dim), ReLU(), Dropout(dropout * 0.5),
            Linear(hidden_dim, hidden_dim // 2), ReLU(),
            Linear(hidden_dim // 2, 1),
        )

    def forward(self, data):
        x, edge_index, edge_attr, batch = data.x, data.edge_index, data.edge_attr, data.batch
        x = self.node_encoder(x)
        for i in range(self.num_layers):
            x_id = x
            x = self.convs[i](x, edge_index, edge_attr)
            x = self.batch_norms[i](x)
            x = F.relu(x)
            x = F.dropout(x, p=self.dropout, training=self.training)
            x = x + x_id
        h_graph = torch.cat([global_mean_pool(x, batch), global_max_pool(x, batch),
                              global_add_pool(x, batch)], dim=1)
        return self.predictor(h_graph)


# ----------------- MAIN PIPELINE -----------------
def main():
    parser = argparse.ArgumentParser()
    parser.add_argument('--input_csv', required=True)
    parser.add_argument('--out_dir', required=True)
    parser.add_argument('--dataset_name', required=True)
    args = parser.parse_args()

    os.makedirs(args.out_dir, exist_ok=True)
    device = torch.device('cuda' if torch.cuda.is_available() else 'cpu')
    print(f"Device: {device}")

    # ---------------------------------------------------------
    # Target scaling, DERIVED FROM THE TRAINING DATA (never hardcoded).
    #
    # History of the bug this replaces: the old code hardcoded
    #     scale_=0.02501143583517879, mean_=0.00015726585024300533
    # and then multiplied by HARTREE_TO_KCAL. That was CORRECT in the era when
    # Delta_G_residual was in HARTREE with the FS2 baseline (training log: "Mean:
    # 0.000157 +/- 0.025011 Hartree"). Since then the target changed twice --
    # units Hartree -> kcal/mol, and baseline FS2 -> FS5 -- so the retrained
    # gin_cv_models use sd = 8.8665 kcal/mol. The hardcoded constants were not
    # updated, and the now-redundant x627.509 stayed in, so uncertainties came out
    # 0.025011*627.509/8.8665 = 1.77x too large. Deriving the scaler from the same
    # CSV the models were trained on makes that class of error impossible.
    #
    # Also: training used a PER-FOLD StandardScaler (gin_residual_v2_cv.py fits it on
    # each fold's train+val split), so each fold's 10 bootstrap models live on a
    # slightly different scale (fold sd spread 8.61-9.07, i.e. 5.1%). We therefore
    # invert each model with ITS OWN fold scaler BEFORE averaging, instead of applying
    # one global scaler to the pooled predictions.
    # ---------------------------------------------------------
    # Preferred source: a sidecar written next to the checkpoints at training time.
    # This is what makes drift impossible -- the scaler travels WITH the models instead
    # of being re-derived from whatever CSV happens to be on disk later.
    sidecar = os.path.join(MODEL_DIR, "target_scaler.json")
    fold_scalers, meta_src = None, None
    if os.path.exists(sidecar):
        with open(sidecar) as fh:
            meta = json.load(fh)
        fold_scalers = [(float(m), float(s)) for m, s in meta["fold_scalers"]]
        meta_src = f"{os.path.basename(sidecar)} (train={meta.get('train_csv')}, " \
                   f"n={meta.get('n')}, units={meta.get('units')})"
    else:
        train_csv = os.environ.get("TRAIN_CSV", TRAIN_CSV)
        if not os.path.isabs(train_csv):
            train_csv = os.path.join(DIR, train_csv)
        tdf = pd.read_csv(train_csv).dropna(subset=['Delta_G_residual', 'smiles', 'mol'])
        y_train = tdf['Delta_G_residual'].values.astype(float)
        groups_tr = np.array([
            (lambda mm: Chem.MolToSmiles(mm) if mm else str(s))(Chem.MolFromSmiles(str(s)))
            for s in tdf['smiles_clean'].astype(str)])
        fold_scalers = []
        for tv, _ in GroupKFold(n_splits=5).split(np.arange(len(y_train)), y_train, groups_tr):
            fold_scalers.append((float(y_train[tv].mean()), float(y_train[tv].std())))
        meta_src = f"{os.path.basename(train_csv)} (n={len(y_train)}, " \
                   f"sd={y_train.std():.4f}) -- NO sidecar, verify this is the CSV the " \
                   f"checkpoints were trained on!"

    sds = np.array([s for _, s in fold_scalers])
    print(f"Target scaler <- {meta_src}")
    print("  per-fold (mean, sd): " + ", ".join(f"({m:+.3f},{s:.3f})" for m, s in fold_scalers))
    if not (1.0 < sds.mean() < 100.0):
        print(f"  WARNING: mean sd {sds.mean():.4g} is outside the expected kcal/mol range "
              f"-- wrong units? Predictions would be mis-scaled.")

    # Read target Dataset
    print(f"Reading dataset: {args.input_csv}")
    df_target = pd.read_csv(args.input_csv)
    
    # Auto-detect SMILES col
    tmp_target = next((c for c in ['SMILES', 'smiles_list', 'smiles'] if c in df_target.columns), None)
    smi_col = tmp_target if tmp_target else df_target.columns[0]
        
    df_target['smiles_clean'] = df_target[smi_col].astype(str).str.replace('*', 'C', regex=False)
    
    pyg_list = []
    valid_idx = []
    for idx, row in df_target.iterrows():
        g = smiles_to_graph(row['smiles_clean'])
        if g is not None:
            pyg_list.append(g)
            valid_idx.append(idx)
            
    df_valid = df_target.iloc[valid_idx].reset_index(drop=True)
    print(f"Successfully featurized {len(df_valid)} graph molecules.")

    # Load the 50 pre-trained GIN models, REMEMBERING each one's fold (its scaler differs)
    all_models = []
    for fold in range(5):
        for mi in range(10):
            ckpt = os.path.join(MODEL_DIR, f"gin_fold{fold}_model{mi}.pt")
            if not os.path.exists(ckpt): continue

            model = GINEResidualModel(node_dim=NODE_DIM, edge_dim=EDGE_DIM, hidden_dim=128, num_layers=4, dropout=0.15).to(device)
            model.load_state_dict(torch.load(ckpt, map_location=device, weights_only=True))
            model.eval()
            all_models.append((fold, model))

    if not all_models:
        print("No models found. Exiting.")
        return
    print(f"Loaded {len(all_models)} models successfully")

    # Inference. Each model is inverted with ITS OWN fold scaler into kcal/mol BEFORE
    # the ensemble statistics are taken -- the models live on slightly different scales
    # (per-fold StandardScaler at training time), so pooling raw scaled outputs and
    # applying one global scaler afterwards is not correct.
    loader = DataLoader(pyg_list, batch_size=64, shuffle=False)
    all_preds = []
    for k, (fold, model) in enumerate(all_models):
        model_preds = []
        with torch.no_grad():
            for data in loader:
                data = data.to(device)
                pred = model(data)
                model_preds.append(pred.cpu().numpy())
        p_scaled = np.vstack(model_preds).flatten()
        f_mean, f_sd = fold_scalers[fold]
        all_preds.append(p_scaled * f_sd + f_mean)          # -> kcal/mol
        if (k + 1) % 10 == 0:
            print(f"  Processed {k + 1}/{len(all_models)} models.")

    preds_stack = np.stack(all_preds)                        # already kcal/mol
    y_mean_kcal = preds_stack.mean(axis=0)
    y_std_kcal = preds_stack.std(axis=0)                     # ensemble sd = uncertainty

    df_valid['Delta_G_residual_kcal_mol'] = y_mean_kcal
    df_valid['Uncertainty_kcal_mol'] = y_std_kcal
    # Hartree columns kept for backward compatibility (derived, not primary)
    df_valid['Delta_G_residual_Hartree'] = y_mean_kcal / HARTREE_TO_KCAL
    df_valid['Uncertainty_Hartree'] = y_std_kcal / HARTREE_TO_KCAL
    print(f"  predicted residual: mean={y_mean_kcal.mean():.3f} sd={y_mean_kcal.std():.3f} kcal/mol")
    print(f"  uncertainty:        median={np.median(y_std_kcal):.3f} "
          f"mean={y_std_kcal.mean():.3f} max={y_std_kcal.max():.3f} kcal/mol")

    # Unique Identifier assignment if PID doesn't exist
    if 'PID' not in df_valid.columns:
        df_valid['PID'] = [f"{args.dataset_name}_{i}" for i in range(len(df_valid))]

    out_csv = os.path.join(args.out_dir, f"{args.dataset_name}_predictions.csv")
    df_valid.to_csv(out_csv, index=False)
    print(f"Base predictions stored to {out_csv}")

    # Plot 1: Prediction vs Uncertainty Distribution
    fig, axes = plt.subplots(1, 2, figsize=(14, 6))
    ax = axes[0]
    ax.hist(df_valid['Delta_G_residual_kcal_mol'], bins=50, color='#3B82F6', edgecolor='white')
    ax.set_title(f"Predicted $\\Delta G$ Residuals ({args.dataset_name})")
    ax.set_xlabel("kcal/mol")
    ax = axes[1]
    ax.hist(df_valid['Uncertainty_kcal_mol'], bins=50, color='#F59E0B', edgecolor='white')
    ax.set_title(f"Ensemble Uncertainty Distribution")
    ax.set_xlabel("Uncertainty (kcal/mol)")
    plt.tight_layout()
    plt.savefig(os.path.join(args.out_dir, f"{args.dataset_name}_gin_prediction_distribution.png"), dpi=300)
    plt.close()

    # Plot 2: Prediction vs Uncertainty Scatter
    fig, ax = plt.subplots(figsize=(10, 7))
    scatter = ax.scatter(df_valid['Delta_G_residual_kcal_mol'], df_valid['Uncertainty_kcal_mol'], 
                         c=df_valid['Uncertainty_kcal_mol'], cmap='RdYlGn_r',
                         s=12, alpha=0.6, edgecolors='none')
    cbar = plt.colorbar(scatter, ax=ax, label='Uncertainty (kcal/mol)')
    ax.set_xlabel('Predicted ΔG Residual (kcal/mol)', fontsize=13, fontweight='bold')
    ax.set_ylabel('Ensemble Uncertainty (kcal/mol)', fontsize=13, fontweight='bold')
    ax.set_title(f'Prediction vs Uncertainty — {args.dataset_name} (N={len(df_valid)})',
                 fontsize=14, fontweight='bold')
    ax.grid(True, alpha=0.2)
    plt.tight_layout()
    plt.savefig(os.path.join(args.out_dir, f"{args.dataset_name}_prediction_vs_uncertainty.png"), dpi=300)
    plt.close()

    # Plot 3: 3D t-SNE Plot
    print("\n  Computing 3D t-SNE")
    n_for_tsne = min(len(df_valid), 5000)
    if n_for_tsne < len(df_valid):
        tsne_idx = np.random.RandomState(42).choice(len(df_valid), n_for_tsne, replace=False)
    else:
        tsne_idx = np.arange(len(df_valid))
    
    fps = np.array([get_morgan_fingerprint(df_valid['smiles_clean'].iloc[i]) for i in tsne_idx])
    tsne = TSNE(n_components=3, random_state=42, perplexity=30, max_iter=1000)
    coords_3d = tsne.fit_transform(fps)
    
    pred_subset = df_valid['Delta_G_residual_kcal_mol'].values[tsne_idx]
    std_subset = df_valid['Uncertainty_kcal_mol'].values[tsne_idx]

    fig = plt.figure(figsize=(22, 9))
    ax1 = fig.add_subplot(121, projection='3d')
    sc1 = ax1.scatter(coords_3d[:, 0], coords_3d[:, 1], coords_3d[:, 2],
                      c=pred_subset, cmap='coolwarm', s=8, alpha=0.7, edgecolors='none')
    fig.colorbar(sc1, ax=ax1, shrink=0.5, pad=0.08, label='ΔG Residual (kcal/mol)')
    ax1.set_title('Predicted ΔG Residual', fontsize=13, fontweight='bold')
    ax1.view_init(elev=25, azim=45)
    
    ax2 = fig.add_subplot(122, projection='3d')
    sc2 = ax2.scatter(coords_3d[:, 0], coords_3d[:, 1], coords_3d[:, 2],
                      c=std_subset, cmap='YlOrRd', s=8, alpha=0.7, edgecolors='none')
    fig.colorbar(sc2, ax=ax2, shrink=0.5, pad=0.08, label='Uncertainty (kcal/mol)')
    ax2.set_title('Prediction Uncertainty', fontsize=13, fontweight='bold')
    ax2.view_init(elev=25, azim=45)
    
    fig.suptitle(f'3D t-SNE Chemical Space — {args.dataset_name} Predictions', fontsize=15, fontweight='bold', y=1.01)
    fig.tight_layout()
    plt.savefig(os.path.join(args.out_dir, f"{args.dataset_name}_tsne_3d_combined.png"), dpi=300, bbox_inches='tight')
    plt.close()

if __name__ == '__main__':
    main()
