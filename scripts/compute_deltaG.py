import os
import pandas as pd
from rdkit import Chem

# =========================================================
# Utility Functions
# =========================================================

def preprocess_smiles(smi: str):
    """Replace '*' with 'C' and return valid RDKit SMILES."""
    if smi is None:
        return None
    return smi.replace("*", "C")


def count_atoms(smiles):
    """Return atom_counts dict including H."""
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return None

    mol_H = Chem.AddHs(mol)
    atom_counts = {}

    for atom in mol_H.GetAtoms():
        symbol = atom.GetSymbol()
        atom_counts[symbol] = atom_counts.get(symbol, 0) + 1

    return atom_counts


def count_bonds(smiles):
    """Return backbone bonds, total bonds, heavy atoms, total atoms."""
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return None, None, None, None, None

    backbone_bonds = mol.GetNumBonds()          # heavy atom bonds

    mol_H = Chem.AddHs(mol)
    total_bonds = mol_H.GetNumBonds()           # total bonds with H

    heavy_atoms = mol.GetNumHeavyAtoms()
    total_atoms = mol_H.GetNumAtoms()

    CH_bonds = total_bonds - backbone_bonds

    return backbone_bonds, total_bonds, heavy_atoms, total_atoms, CH_bonds


# =========================================================
# Main ΔG Computation Pipeline
# =========================================================

def compute_deltaG():

    # script/ → project_root/
    root = os.path.dirname(os.path.dirname(__file__))

    # --- YOUR INPUT FILE IS HERE ---
    input_csv = os.path.join(root, "data", "PI1070.csv")

    atom_ref_csv = os.path.join(root, "data/atom_ref", "atom_ref.csv")
    deltaG_raw_csv = os.path.join(root, "data", "deltaG_raw.csv")

    print("Reading input SMILES file:", input_csv)
    df = pd.read_csv(input_csv)

    print("Reading atom reference energy:", atom_ref_csv)
    atom_ref = pd.read_csv(atom_ref_csv)
    ref_dict = dict(zip(atom_ref["atom"], atom_ref["energy"]))

    print("Reading Gibbs energies:", deltaG_raw_csv)
    df_G = pd.read_csv(deltaG_raw_csv)

    if "Gibbs_Eh" not in df_G.columns:
        raise ValueError("deltaG_raw.csv missing column: Gibbs_Eh")

    # --- FIX: sort PI order numerically ---
    df_G["mol_index"] = df_G["mol"].apply(lambda x: int(x.replace("PI", "")))
    df_G = df_G.sort_values("mol_index").reset_index(drop=True)
    df_G = df_G.drop(columns=["mol_index"])

    # Attach sorted Gibbs energy back to main df
    df["Gibbs_Eh"] = df_G["Gibbs_Eh"]


    # =====================================================
    # Clean SMILES
    # =====================================================
    df["smiles_clean"] = df["smiles"].astype(str).apply(preprocess_smiles)

    # =====================================================
    # Count atoms and bonds
    # =====================================================
    atom_results = []
    backbone_list, total_bond_list, heavy_atom_list, total_atom_list, CH_list = [], [], [], [], []

    print("Counting atoms + bonds...")

    for i, smi in enumerate(df["smiles_clean"]):
        if i % 20 == 0:
            print(f"  → row {i}/{len(df)}")

        atom_counts = count_atoms(smi)
        b1, b2, ha, ta, ch = count_bonds(smi)

        atom_results.append(atom_counts)
        backbone_list.append(b1)
        total_bond_list.append(b2)
        heavy_atom_list.append(ha)
        total_atom_list.append(ta)
        CH_list.append(ch)

    df["atom_counts"] = atom_results
    df["backbone_bonds"] = backbone_list
    df["total_bonds"] = total_bond_list
    df["heavy_atoms"] = heavy_atom_list
    df["total_atoms"] = total_atom_list
    df["CH_bonds"] = CH_list

    # =====================================================
    # Compute reference energy sum
    # =====================================================
    def compute_reference_energy(atom_dict):
        total_ref_E = 0.0
        for atom, n in atom_dict.items():
            if atom not in ref_dict:
                raise ValueError(f"Missing reference energy for atom '{atom}' in atom_ref.csv")
            total_ref_E += n * ref_dict[atom]
        return total_ref_E

    df["E_ref_sum"] = df["atom_counts"].apply(compute_reference_energy)

    # =====================================================
    # Compute ΔG
    # =====================================================
    df["Delta_G"] = df["Gibbs_Eh"] - df["E_ref_sum"]

    # =====================================================
    # Normalized ΔG
    # =====================================================
    df["DeltaG_per_heavy_atom"] = df["Delta_G"] / df["heavy_atoms"]
    df["DeltaG_per_atom"] = df["Delta_G"] / df["total_atoms"]
    df["DeltaG_per_backbone_bond"] = df["Delta_G"] / df["backbone_bonds"]
    df["DeltaG_per_bond"] = df["Delta_G"] / df["total_bonds"]
    df["DeltaG_per_CH"] = df["Delta_G"] / df["CH_bonds"]

    # =====================================================
    # Save output
    # =====================================================
    output_csv = os.path.join(root, "final_data_with_deltaG.csv")
    df.to_csv(output_csv, index=False)

    print("🎉 DONE! Saved to:", output_csv)
    print(df.head())

    return df


# Run
if __name__ == "__main__":
    compute_deltaG()
