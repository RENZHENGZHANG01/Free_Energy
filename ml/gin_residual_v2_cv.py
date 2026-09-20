"""
GIN v2 — 5-Fold Cross-Validation on All Available 2D Data
==========================================================

This script is based on gin_residual_v2_cv.py with the following changes:
  1. Removed the 3D xyz-coordinate filtering step.
  2. Uses all valid molecules in final_data_with_residual_deltaG.csv.
  3. Treats Delta_G_residual as already being in kcal/mol.
  4. Uses publication-style plotting settings inspired by Nature-style figures.

Core model/training logic is kept from gin_residual_v2_cv.py:
  - 2D molecular graph from SMILES without explicit hydrogen expansion.
  - 5-fold random KFold cross-validation.
  - Per-fold target standardization using training/validation fold only.
  - 10-model bootstrap ensemble per fold.
"""

import os
import pandas as pd
import numpy as np
import torch
import torch.nn.functional as F
import random
from torch.nn import Linear, Sequential, ReLU, BatchNorm1d, Dropout
from torch_geometric.data import Data
from torch_geometric.loader import DataLoader
from torch_geometric.nn import GINEConv, global_mean_pool, global_max_pool, global_add_pool
from rdkit import Chem
from rdkit.Chem import rdchem
from sklearn.model_selection import KFold, train_test_split, GroupKFold, GroupShuffleSplit
from sklearn.preprocessing import StandardScaler
from sklearn.metrics import r2_score, mean_absolute_error
from torch.optim.lr_scheduler import CosineAnnealingWarmRestarts
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt

# =====================================================================
#  Publication-style plotting settings
# =====================================================================
plt.rcParams.update({
    "font.family": "sans-serif",
    "font.sans-serif": ["Arial", "Helvetica", "DejaVu Sans"],
    "font.size": 18,
    "axes.labelsize": 20,
    "axes.titlesize": 22,
    "xtick.labelsize": 18,
    "ytick.labelsize": 18,
    "legend.fontsize": 15,
    "figure.dpi": 300,
    "savefig.dpi": 600,
    "axes.linewidth": 1.5,
    "axes.edgecolor": "black",
    "axes.spines.top": True,
    "axes.spines.right": True,
    "xtick.direction": "in",
    "ytick.direction": "in",
    "xtick.top": True,
    "ytick.right": True,
    "xtick.major.width": 1.5,
    "ytick.major.width": 1.5,
    "xtick.minor.width": 1.2,
    "ytick.minor.width": 1.2,
    "pdf.fonttype": 42,
    "ps.fonttype": 42,
})

DIR = os.path.dirname(os.path.abspath(__file__))
CSV_PATH = os.path.join(DIR, "final_data_with_residual_deltaG.csv")
PLOT_DIR = os.path.join(DIR, "gin_cv_plots")
MODEL_DIR = os.path.join(DIR, "gin_cv_models")
os.makedirs(PLOT_DIR, exist_ok=True)
os.makedirs(MODEL_DIR, exist_ok=True)

SEED = 42
N_FOLDS = 5
random.seed(SEED)
np.random.seed(SEED)
torch.manual_seed(SEED)


# =====================================================================
#  Feature Engineering (same as gin_residual_v2.py)
# =====================================================================
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
# 14 element one-hot (13 symbols + "other") + 3 (EN, vdW, mass) + 5 hybridization
# + 6 atom flags + 6 ring-size flags = 34
NODE_DIM = 34
EDGE_DIM = 10


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


def smiles_to_graph(smiles, target):
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
                y=torch.tensor([target], dtype=torch.float).view(1, -1))


# =====================================================================
#  Model (same as gin_residual_v2.py)
# =====================================================================
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


# =====================================================================
#  Ensemble (modified for CV)
# =====================================================================
class GNNEnsemble:
    def __init__(self, n_models, fold_id, node_dim, edge_dim,
                 hidden_dim=128, num_layers=4, dropout=0.15):
        self.n_models = n_models
        self.fold_id = fold_id
        self.device = torch.device('cuda' if torch.cuda.is_available() else 'cpu')
        self.models = [GINEResidualModel(node_dim, edge_dim, hidden_dim,
                                         num_layers, dropout).to(self.device)
                       for _ in range(n_models)]
        n_params = sum(p.numel() for p in self.models[0].parameters())
        print(f"    Parameters per model: {n_params:,}")

    def train_ensemble(self, train_dataset, val_dataset, max_epochs=300):
        for i, model in enumerate(self.models):
            print(f"\n    Ensemble {i+1}/{self.n_models}", end=" — ")

            rng = random.Random(SEED + self.fold_id * 10000 + i * 1000)
            bootstrap = rng.choices(train_dataset, k=len(train_dataset))
            train_loader = DataLoader(bootstrap, batch_size=32, shuffle=True)
            val_loader = DataLoader(val_dataset, batch_size=32)

            optimizer = torch.optim.AdamW(model.parameters(), lr=1e-3, weight_decay=1e-4)
            scheduler = CosineAnnealingWarmRestarts(optimizer, T_0=50, T_mult=2, eta_min=1e-6)
            criterion = torch.nn.SmoothL1Loss()

            best_val_loss = float('inf')
            patience_counter = 0
            patience_limit = 30
            ckpt_path = os.path.join(MODEL_DIR, f"gin_fold{self.fold_id}_model{i}.pt")

            for epoch in range(1, max_epochs + 1):
                model.train()
                epoch_loss, n_samples = 0.0, 0
                for data in train_loader:
                    data = data.to(self.device)
                    optimizer.zero_grad()
                    pred = model(data)
                    loss = criterion(pred, data.y)
                    loss.backward()
                    torch.nn.utils.clip_grad_norm_(model.parameters(), max_norm=1.0)
                    optimizer.step()
                    epoch_loss += loss.item() * data.num_graphs
                    n_samples += data.num_graphs
                epoch_loss /= n_samples

                model.eval()
                val_loss, n_val = 0.0, 0
                with torch.no_grad():
                    for data in val_loader:
                        data = data.to(self.device)
                        val_loss += criterion(model(data), data.y).item() * data.num_graphs
                        n_val += data.num_graphs
                val_loss /= n_val
                scheduler.step()

                if val_loss < best_val_loss:
                    best_val_loss = val_loss
                    torch.save(model.state_dict(), ckpt_path)
                    patience_counter = 0
                else:
                    patience_counter += 1

                if patience_counter >= patience_limit:
                    break

            model.load_state_dict(torch.load(ckpt_path, weights_only=True))
            print(f"stopped ep {epoch}, best_val={best_val_loss:.6f}")

    def predict(self, dataset):
        loader = DataLoader(dataset, batch_size=32)
        per_model = []
        for model in self.models:
            model.eval()
            preds = []
            with torch.no_grad():
                for data in loader:
                    data = data.to(self.device)
                    preds.append(model(data).cpu().numpy())
            per_model.append(np.vstack(preds))
        stacked = np.stack(per_model)
        return stacked.mean(axis=0).flatten(), stacked.std(axis=0).flatten()


# =====================================================================
#  Main — 5-Fold Cross-Validation
# =====================================================================
def main():
    print("=" * 60)
    print("  GIN v2 — 5-Fold Cross-Validation")
    print("=" * 60)

    # ── Load CSV ──
    df = pd.read_csv(CSV_PATH)
    target_col = 'Delta_G_residual'
    df = df.dropna(subset=[target_col, 'smiles', 'mol']).reset_index(drop=True)

    # ── Use all available molecules; no 3D xyz-coordinate filtering ──
    print(f"\n  Total usable CSV rows after dropping missing values: {len(df)}")

    y_raw = df[target_col].values.astype(float)
    print(f"  Target: {target_col}")
    print(f"  Unit: kcal/mol")
    print(f"  Range: [{y_raw.min():.6f}, {y_raw.max():.6f}] kcal/mol")
    print(f"  Mean:  {y_raw.mean():.6f} ± {y_raw.std():.6f} kcal/mol")

    # ── Reference set of ORIGINAL molecule IDs (for old-vs-new subset eval) ──
    # Molecules present in the pre-active-learning dataset are "old"; everything
    # else is a "new" AL point. Used to report subset-resolved CV metrics.
    REF_CSV = os.path.join(DIR, "final_data_with_residual_deltaG_june24.csv")
    if os.path.exists(REF_CSV):
        original_mol_ids = set(pd.read_csv(REF_CSV)['mol'].astype(str))
        print(f"  Reference: {len(original_mol_ids)} original mol IDs "
              f"from {os.path.basename(REF_CSV)} (old-vs-new subset eval enabled)")
    else:
        original_mol_ids = None
        print(f"  ⚠️ {os.path.basename(REF_CSV)} not found — old/new subset eval disabled")

    # ── Build graphs (with raw targets) ──
    print("\n  Converting SMILES to molecular graphs...")
    all_smiles = df['smiles'].tolist()
    all_smiles_clean = df['smiles_clean'].tolist()
    graphs_raw = []
    valid_idx = []
    groups = []
    for i, smi in enumerate(all_smiles):
        g = smiles_to_graph(smi, y_raw[i])  # store raw target
        if g is not None:
            graphs_raw.append(g)
            valid_idx.append(i)
            # Group key for leakage-free CV (aligned with MACE): canonical SMILES
            # of the cleaned structure, so identical structures never split across
            # train / test. Falls back to the raw SMILES if parsing fails.
            mol_clean = Chem.MolFromSmiles(str(all_smiles_clean[i]))
            groups.append(Chem.MolToSmiles(mol_clean) if mol_clean else smi)
    y_raw = y_raw[valid_idx]
    groups = np.array(groups)

    # Track molecule IDs (parallel to graphs) to resolve old/new subsets later.
    valid_mols = df['mol'].astype(str).values[valid_idx]
    if original_mol_ids is not None:
        is_new_all = np.array([m not in original_mol_ids for m in valid_mols])
        print(f"  Subset sizes — original (old): {(~is_new_all).sum()} "
              f"| new (AL): {int(is_new_all.sum())}")

    # Optional size descriptor for diagnostic plotting.
    # Prefer an existing total_atoms column; otherwise compute total atoms from RDKit.
    if 'total_atoms' in df.columns:
        total_atoms_all = df['total_atoms'].values[valid_idx].astype(float)
    else:
        total_atoms_all = []
        for smi in np.array(all_smiles, dtype=object)[valid_idx]:
            mol = Chem.MolFromSmiles(str(smi))
            total_atoms_all.append(float(sum(atom.GetTotalNumHs() + 1 for atom in mol.GetAtoms())) if mol is not None else np.nan)
        total_atoms_all = np.array(total_atoms_all, dtype=float)

    print(f"  Valid graphs: {len(graphs_raw)}")
    print(f"  Node features: {NODE_DIM}, Edge features: {EDGE_DIM}")

    # ══════════════════════════════════════════════════════════════════
    #  5-Fold Cross-Validation
    # ══════════════════════════════════════════════════════════════════
    gkf = GroupKFold(n_splits=N_FOLDS)
    indices = np.arange(len(graphs_raw))

    fold_results = []
    all_y_true_kcal = []
    all_y_pred_kcal = []
    all_y_std_kcal = []
    all_total_atoms = []
    all_mols = []

    saved_fold_scalers = []          # persisted next to the checkpoints (see below)

    for fold_idx, (train_val_idx, test_idx) in enumerate(gkf.split(indices, y_raw, groups)):
        print(f"\n{'═' * 60}")
        print(f"  FOLD {fold_idx + 1}/{N_FOLDS}")
        print(f"{'═' * 60}")

        # ── Per-fold target scaling ──
        y_fold_train_val = y_raw[train_val_idx]
        scaler = StandardScaler()
        scaler.fit(y_fold_train_val.reshape(-1, 1))
        y_all_scaled = scaler.transform(y_raw.reshape(-1, 1)).flatten()
        saved_fold_scalers.append([float(scaler.mean_[0]), float(scaler.scale_[0])])

        # Rebuild graphs with scaled targets
        graphs_scaled = []
        for idx in range(len(graphs_raw)):
            g_orig = graphs_raw[idx]
            g_new = Data(x=g_orig.x, edge_index=g_orig.edge_index,
                         edge_attr=g_orig.edge_attr,
                         y=torch.tensor([y_all_scaled[idx]], dtype=torch.float).view(1, -1))
            graphs_scaled.append(g_new)

        # Split train_val into train + val (grouped, to prevent leakage into early stopping)
        gss = GroupShuffleSplit(n_splits=1, test_size=0.15, random_state=SEED + fold_idx)
        train_sub_idx, val_sub_idx = next(gss.split(train_val_idx, groups=groups[train_val_idx]))
        train_idx = train_val_idx[train_sub_idx]
        val_idx = train_val_idx[val_sub_idx]
        train_set = [graphs_scaled[i] for i in train_idx]
        val_set = [graphs_scaled[i] for i in val_idx]
        test_set = [graphs_scaled[i] for i in test_idx]

        print(f"  Train: {len(train_set)} | Val: {len(val_set)} | Test: {len(test_set)}")

        # ── Train ensemble ──
        ensemble = GNNEnsemble(
            n_models=10, fold_id=fold_idx,
            node_dim=NODE_DIM, edge_dim=EDGE_DIM,
            hidden_dim=128, num_layers=4, dropout=0.15,
        )
        ensemble.train_ensemble(train_set, val_set, max_epochs=300)

        # ── Evaluate ──
        y_pred_s, y_std_s = ensemble.predict(test_set)
        y_true_s = np.array([g.y.item() for g in test_set])

        y_pred = scaler.inverse_transform(y_pred_s.reshape(-1, 1)).flatten()
        y_true = scaler.inverse_transform(y_true_s.reshape(-1, 1)).flatten()
        y_std = y_std_s * scaler.scale_[0]

        # Delta_G_residual is already in kcal/mol.
        y_pred_kcal = y_pred
        y_true_kcal = y_true
        y_std_kcal_fold = y_std

        r2 = r2_score(y_true, y_pred)
        mae_kcal = mean_absolute_error(y_true_kcal, y_pred_kcal)

        fold_results.append({'fold': fold_idx + 1, 'r2': r2, 'mae_kcal': mae_kcal})
        all_y_true_kcal.extend(y_true_kcal.tolist())
        all_y_pred_kcal.extend(y_pred_kcal.tolist())
        all_y_std_kcal.extend(y_std_kcal_fold.tolist())
        all_total_atoms.extend(total_atoms_all[test_idx].tolist())
        all_mols.extend(valid_mols[test_idx].tolist())

        print(f"\n  Fold {fold_idx+1} Results: R²={r2:.4f}, MAE={mae_kcal:.2f} kcal/mol")

        del ensemble
        torch.cuda.empty_cache() if torch.cuda.is_available() else None

    # ══════════════════════════════════════════════════════════════════
    #  Persist the target scaler NEXT TO the checkpoints
    # ══════════════════════════════════════════════════════════════════
    # Inference must undo exactly this scaling. Writing it here (rather than letting
    # predict_dataset_gin.py hardcode or re-derive it) is what prevents scaler drift:
    # the target has already changed units (Hartree->kcal/mol) and baseline (FS2->FS5)
    # once, which silently mis-scaled predictions by 1.77x until 2026-07-24.
    import json as _json
    _meta = {
        "train_csv": os.path.basename(CSV_PATH),
        "n": int(len(y_raw)),
        "units": "kcal/mol",
        "target": target_col,
        "global_mean": float(np.mean(y_raw)),
        "global_sd": float(np.std(y_raw)),
        "fold_scalers": saved_fold_scalers,        # [[mean, sd], ...] per fold
        "note": "Written by gin_residual_v2_cv.py. predict_dataset_gin.py reads this "
                "and inverts each fold's models with its own scaler.",
    }
    with open(os.path.join(MODEL_DIR, "target_scaler.json"), "w") as _fh:
        _json.dump(_meta, _fh, indent=2)
    print(f"\n  💾 Saved target scaler -> {os.path.join(MODEL_DIR, 'target_scaler.json')}"
          f"  (sd={_meta['global_sd']:.6f} {_meta['units']})")

    # ══════════════════════════════════════════════════════════════════
    #  Aggregate Results
    # ══════════════════════════════════════════════════════════════════
    r2_values = [r['r2'] for r in fold_results]
    mae_values = [r['mae_kcal'] for r in fold_results]

    print(f"\n{'═' * 60}")
    print(f"  5-FOLD CV RESULTS — GIN v2 (2D)")
    print(f"{'═' * 60}")
    for r in fold_results:
        print(f"  Fold {r['fold']}: R²={r['r2']:.4f}, MAE={r['mae_kcal']:.2f} kcal/mol")
    print(f"  {'─' * 50}")
    print(f"  Mean R²:  {np.mean(r2_values):.4f} ± {np.std(r2_values):.4f}")
    print(f"  Mean MAE: {np.mean(mae_values):.2f} ± {np.std(mae_values):.2f} kcal/mol")
    print(f"{'═' * 60}")

    # ══════════════════════════════════════════════════════════════════
    #  Old (original) vs New (active-learning) subset evaluation
    # ══════════════════════════════════════════════════════════════════
    all_y_true_arr = np.array(all_y_true_kcal)
    all_y_pred_arr = np.array(all_y_pred_kcal)
    all_y_std_arr  = np.array(all_y_std_kcal)
    all_atoms_arr  = np.array(all_total_atoms, dtype=float)
    all_mols_arr   = np.array(all_mols, dtype=object)

    if original_mol_ids is not None:
        is_new_pred = np.array([m not in original_mol_ids for m in all_mols_arr])
    else:
        is_new_pred = np.zeros(len(all_mols_arr), dtype=bool)

    # Save per-molecule out-of-fold CV predictions for downstream analysis
    pred_df = pd.DataFrame({
        'mol': all_mols_arr,
        'y_true': all_y_true_arr,
        'y_pred': all_y_pred_arr,
        'y_std': all_y_std_arr,
        'abs_error': np.abs(all_y_pred_arr - all_y_true_arr),
        'total_atoms': all_atoms_arr,
        'is_new': is_new_pred,
    })
    pred_csv = os.path.join(PLOT_DIR, 'cv_predictions.csv')
    pred_df.to_csv(pred_csv, index=False)
    print(f"\n  📄 Per-molecule CV predictions saved: {pred_csv}")

    if original_mol_ids is not None and is_new_pred.any():
        print(f"\n{'═' * 60}")
        print(f"  SUBSET EVALUATION — original (old) vs new (AL) points")
        print(f"{'═' * 60}")
        for label, mask in [("Original (old)", ~is_new_pred),
                            ("New (AL)", is_new_pred),
                            ("All", np.ones_like(is_new_pred))]:
            n = int(mask.sum())
            if n < 2:
                print(f"  {label:16s}: n={n} (too few to score)")
                continue
            yt, yp = all_y_true_arr[mask], all_y_pred_arr[mask]
            print(f"  {label:16s}: n={n:5d}  R²={r2_score(yt, yp):.4f}  "
                  f"MAE={mean_absolute_error(yt, yp):.2f} kcal/mol")
        print(f"{'═' * 60}")
        print(f"  → If 'Original (old)' R² is still ~0.94, adding AL data did NOT")
        print(f"    degrade the original chemical space; the global drop is driven")
        print(f"    by the harder, sparsely-sampled new region.")

        # Parity plot colored by subset
        fig, ax = plt.subplots(figsize=(8, 8))
        for mask, color, lab in [(~is_new_pred, '#2563EB', 'Original (old)'),
                                 (is_new_pred, '#dc2626', 'New (AL)')]:
            ax.scatter(all_y_true_arr[mask], all_y_pred_arr[mask], s=34, alpha=0.65,
                       color=color, edgecolors='white', linewidths=0.4, label=lab)
        lims = [min(all_y_true_arr.min(), all_y_pred_arr.min()) - 2,
                max(all_y_true_arr.max(), all_y_pred_arr.max()) + 2]
        ax.plot(lims, lims, color='black', linestyle='--', lw=2.0, zorder=5)
        ax.set_xlim(lims); ax.set_ylim(lims)
        ax.set_xlabel(r'Ground-truth $\Delta G_{\mathrm{res}}$ (kcal mol$^{-1}$)')
        ax.set_ylabel(r'Predicted $\Delta G_{\mathrm{res}}$ (kcal mol$^{-1}$)')
        ax.set_title('CV Parity — Original vs New (AL) subsets', pad=12)
        ax.legend(frameon=True, framealpha=0.9, loc='upper left')
        ax.grid(False); ax.set_box_aspect(1)
        fig.tight_layout()
        fig.savefig(os.path.join(PLOT_DIR, 'parity_cv_by_subset.png'), bbox_inches='tight')
        fig.savefig(os.path.join(PLOT_DIR, 'parity_cv_by_subset.pdf'), bbox_inches='tight')
        plt.close(fig)
        print(f"  📊 {os.path.join(PLOT_DIR, 'parity_cv_by_subset.png')}")

    # ══════════════════════════════════════════════════════════════════
    #  Plots
    # ══════════════════════════════════════════════════════════════════
    all_y_true_kcal = np.array(all_y_true_kcal)
    all_y_pred_kcal = np.array(all_y_pred_kcal)
    all_y_std_kcal = np.array(all_y_std_kcal)
    all_total_atoms = np.array(all_total_atoms, dtype=float)

    # --- Parity ---
    fig, ax = plt.subplots(figsize=(8, 8))
    ax.errorbar(
        all_y_true_kcal, all_y_pred_kcal,
        yerr=all_y_std_kcal,
        fmt='o', alpha=0.72,
        ecolor='lightgray', elinewidth=0.8, capsize=2,
        markersize=6, markerfacecolor='#2563EB',
        markeredgecolor='white', markeredgewidth=0.4,
        linestyle='none'
    )
    lims = [
        min(all_y_true_kcal.min(), all_y_pred_kcal.min()) - 2,
        max(all_y_true_kcal.max(), all_y_pred_kcal.max()) + 2,
    ]
    ax.plot(lims, lims, color='#dc2626', linestyle='--', lw=2.2, zorder=5)
    ax.set_xlim(lims)
    ax.set_ylim(lims)
    ax.set_xlabel(r'Ground-truth $\Delta G_{\mathrm{res}}$ (kcal mol$^{-1}$)')
    ax.set_ylabel(r'Predicted $\Delta G_{\mathrm{res}}$ (kcal mol$^{-1}$)')

    stats_text = (
        f'$R^2$ = {np.mean(r2_values):.3f} $\\pm$ {np.std(r2_values):.3f}\n'
        f'MAE = {np.mean(mae_values):.3f} $\\pm$ {np.std(mae_values):.3f} kcal mol$^{{-1}}$'
    )
    ax.text(
        0.05, 0.95, stats_text, transform=ax.transAxes,
        fontsize=15, va='top', ha='left',
        bbox=dict(boxstyle='round,pad=0.35', facecolor='white', edgecolor='0.75', alpha=0.85)
    )
    ax.set_title('GIN v2 — 5-Fold Cross-Validation Parity', pad=12)
    ax.grid(False)
    ax.set_box_aspect(1)
    fig.tight_layout()
    fig.savefig(os.path.join(PLOT_DIR, 'parity_cv.png'), bbox_inches='tight')
    fig.savefig(os.path.join(PLOT_DIR, 'parity_cv.pdf'), bbox_inches='tight')
    plt.close(fig)
    print(f"\n  📊 {os.path.join(PLOT_DIR, 'parity_cv.png')}")
    print(f"  📄 {os.path.join(PLOT_DIR, 'parity_cv.pdf')}")

    # --- Per-fold R² ---
    fig, ax = plt.subplots(figsize=(8, 5.5))
    folds_x = np.arange(1, len(fold_results) + 1)
    ax.bar(
        folds_x, r2_values,
        color='#2563EB', alpha=0.85,
        edgecolor='white', linewidth=1.4
    )
    ax.axhline(
        y=np.mean(r2_values), color='#dc2626', linestyle='--', lw=2.2,
        label=f'Mean = {np.mean(r2_values):.3f}'
    )
    ax.set_xticks(folds_x)
    ax.set_xticklabels([f'Fold {i}' for i in folds_x])
    ax.set_ylabel(r'$R^2$')
    ax.set_title(r'Per-Fold $R^2$ Scores', pad=12)
    ax.set_ylim(0, 1)
    ax.legend(frameon=True, framealpha=0.9, edgecolor='0.8')
    ax.grid(False)
    fig.tight_layout()
    fig.savefig(os.path.join(PLOT_DIR, 'fold_r2.png'), bbox_inches='tight')
    fig.savefig(os.path.join(PLOT_DIR, 'fold_r2.pdf'), bbox_inches='tight')
    plt.close(fig)
    print(f"  📊 {os.path.join(PLOT_DIR, 'fold_r2.png')}")
    print(f"  📄 {os.path.join(PLOT_DIR, 'fold_r2.pdf')}")

    # --- Uncertainty calibration ---
    all_errors = np.abs(all_y_true_kcal - all_y_pred_kcal)
    fig, ax = plt.subplots(figsize=(8, 8))
    ax.scatter(
        all_y_std_kcal, all_errors,
        s=36, alpha=0.65, c='#D97706',
        edgecolors='white', linewidths=0.4
    )
    max_val = max(all_y_std_kcal.max(), all_errors.max())
    ax.plot([0, max_val], [0, max_val], color='#dc2626', linestyle='--', lw=2.2, alpha=0.75)
    ax.set_xlim(left=0)
    ax.set_ylim(bottom=0)
    ax.set_xlabel(r'Ensemble standard deviation (kcal mol$^{-1}$)')
    ax.set_ylabel(r'Absolute error (kcal mol$^{-1}$)')
    ax.set_title('Uncertainty Calibration', pad=12)
    ax.grid(False)
    ax.set_box_aspect(1)
    fig.tight_layout()
    fig.savefig(os.path.join(PLOT_DIR, 'uncertainty_cv.png'), bbox_inches='tight')
    fig.savefig(os.path.join(PLOT_DIR, 'uncertainty_cv.pdf'), bbox_inches='tight')
    plt.close(fig)
    print(f"  📊 {os.path.join(PLOT_DIR, 'uncertainty_cv.png')}")
    print(f"  📄 {os.path.join(PLOT_DIR, 'uncertainty_cv.pdf')}")

    # --- Error distribution ---
    errors = all_y_pred_kcal - all_y_true_kcal
    fig, ax = plt.subplots(figsize=(8, 6))
    ax.hist(errors, bins=30, color='#2563EB', edgecolor='white', linewidth=0.6, alpha=0.85)
    ax.axvline(x=0, color='#dc2626', linestyle='--', lw=2.2)
    ax.set_xlabel(r'Prediction error (kcal mol$^{-1}$)')
    ax.set_ylabel('Count')
    ax.set_title(
        f'Error Distribution\nMean = {errors.mean():.2f}, Std = {errors.std():.2f} kcal mol$^{{-1}}$',
        pad=12
    )
    ax.grid(False)
    fig.tight_layout()
    fig.savefig(os.path.join(PLOT_DIR, 'error_dist_cv.png'), bbox_inches='tight')
    fig.savefig(os.path.join(PLOT_DIR, 'error_dist_cv.pdf'), bbox_inches='tight')
    plt.close(fig)
    print(f"  📊 {os.path.join(PLOT_DIR, 'error_dist_cv.png')}")
    print(f"  📄 {os.path.join(PLOT_DIR, 'error_dist_cv.pdf')}")

    # --- Predicted residual vs molecular size diagnostic ---
    valid_atom_mask = np.isfinite(all_total_atoms)
    if valid_atom_mask.sum() >= 3:
        try:
            from scipy.stats import pearsonr
            r_val, p_val = pearsonr(all_total_atoms[valid_atom_mask], all_y_pred_kcal[valid_atom_mask])
            stats_text = f'Pearson $r$ = {r_val:.3f}\n$p$ = {p_val:.2e}'
        except Exception:
            r_val = np.corrcoef(all_total_atoms[valid_atom_mask], all_y_pred_kcal[valid_atom_mask])[0, 1]
            stats_text = f'Pearson $r$ = {r_val:.3f}'

        fig, ax = plt.subplots(figsize=(8, 8))
        ax.scatter(
            all_total_atoms[valid_atom_mask], all_y_pred_kcal[valid_atom_mask],
            alpha=0.72, color='#2563EB',
            edgecolors='white', linewidths=0.4, s=42
        )
        ax.axhline(y=0, color='#dc2626', linestyle='--', lw=2.2, alpha=0.75)
        ax.set_xlabel('Total atoms')
        ax.set_ylabel(r'Predicted $\Delta G_{\mathrm{res}}$ (kcal mol$^{-1}$)')
        ax.text(
            0.05, 0.95, stats_text, transform=ax.transAxes,
            fontsize=15, va='top', ha='left',
            bbox=dict(boxstyle='round,pad=0.35', facecolor='white', edgecolor='0.75', alpha=0.85)
        )
        ax.set_title('Predicted Residual vs Molecular Size', pad=12)
        ax.grid(False)
        ax.set_box_aspect(1)
        fig.tight_layout()
        fig.savefig(os.path.join(PLOT_DIR, 'residual_vs_atoms_cv.png'), bbox_inches='tight')
        fig.savefig(os.path.join(PLOT_DIR, 'residual_vs_atoms_cv.pdf'), bbox_inches='tight')
        plt.close(fig)
        print(f"  📊 {os.path.join(PLOT_DIR, 'residual_vs_atoms_cv.png')}")
        print(f"  📄 {os.path.join(PLOT_DIR, 'residual_vs_atoms_cv.pdf')}")

    print(f"\n🎉 DONE!")


if __name__ == "__main__":
    main()
