import os
import pandas as pd
from rdkit import Chem
from rdkit.Chem import Descriptors

# =========================================================
# Utility Functions
# =========================================================

def preprocess_smiles(smi: str):
    """Replace '*' with 'C' and return valid RDKit SMILES."""
    if smi is None or (isinstance(smi, float) and pd.isna(smi)):
        return None
    return str(smi).replace("*", "C")


def count_atoms(smiles):
    """Return atom_counts dict including H."""
    if smiles is None:
        return None
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
    """Return backbone bonds, total bonds, heavy atoms, total atoms, CH_bonds."""
    if smiles is None:
        return None, None, None, None, None
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


def compute_molecular_weight(smiles):
    """Compute molecular weight from SMILES string."""
    if smiles is None:
        return None
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return None
    return Descriptors.MolWt(mol)


# =========================================================
# Main ΔG Computation Pipeline
# =========================================================

def compute_deltaG():

    # All files are in the same directory as this script
    script_dir = os.path.dirname(os.path.abspath(__file__))

    # --- Input: the reordered file with mol, smiles, Gibbs_Eh, G_minus_Eel ---
    input_csv = os.path.join(script_dir, "merged_G_raw.csv")
    atom_ref_csv = os.path.join(script_dir, "atom_ref.csv")

    print("Reading input file:", input_csv)
    df = pd.read_csv(input_csv)

    print("Reading atom reference energy:", atom_ref_csv)
    atom_ref = pd.read_csv(atom_ref_csv)
    ref_dict = dict(zip(atom_ref["atom"], atom_ref["energy"]))

    # =====================================================
    # Clean SMILES (* → C)
    # =====================================================
    df["smiles_clean"] = df["smiles"].apply(preprocess_smiles)

    # =====================================================
    # Count atoms and bonds
    # =====================================================
    atom_results = []
    backbone_list, total_bond_list, heavy_atom_list, total_atom_list, CH_list = [], [], [], [], []

    print("Counting atoms + bonds...")

    for i, smi in enumerate(df["smiles_clean"]):
        if i % 100 == 0:
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
    # Compute molecular weight
    # =====================================================
    print("Computing molecular weights...")
    df["MW"] = df["smiles_clean"].apply(compute_molecular_weight)

    # =====================================================
    # Compute reference energy sum
    # =====================================================
    def compute_reference_energy(atom_dict):
        if atom_dict is None:
            return None
        total_ref_E = 0.0
        for atom, n in atom_dict.items():
            if atom not in ref_dict:
                raise ValueError(f"Missing reference energy for atom '{atom}' in atom_ref.csv")
            total_ref_E += n * ref_dict[atom]
        return total_ref_E

    df["E_ref_sum"] = df["atom_counts"].apply(compute_reference_energy)

    # =====================================================
    # Compute ΔG (only for rows that have Gibbs_Eh data)
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
    df["Delta_G_per_MW"] = df["Delta_G"] / df["MW"]

    # =====================================================
    # Select output columns (drop intermediate atom_counts dict)
    # =====================================================
    output_cols = [
        "mol", "smiles", "smiles_clean", "MW",
        "Gibbs_Eh", "G_minus_Eel", "E_ref_sum", "Delta_G",
        "heavy_atoms", "total_atoms", "backbone_bonds", "total_bonds", "CH_bonds",
        "DeltaG_per_heavy_atom", "DeltaG_per_atom",
        "DeltaG_per_backbone_bond", "DeltaG_per_bond", "DeltaG_per_CH",
        "Delta_G_per_MW",
    ]
    df_out = df[output_cols]

    # =====================================================
    # Save output
    # =====================================================
    output_csv = os.path.join(script_dir, "final_data_with_deltaG.csv")
    df_out.to_csv(output_csv, index=False)

    print(f"\n🎉 DONE! Saved to: {output_csv}")
    print(f"   Total rows: {len(df_out)}")
    print(f"   Rows with Delta_G: {df_out['Delta_G'].notna().sum()}")
    print(f"   Rows without Gibbs data: {df_out['Delta_G'].isna().sum()}")
    print(df_out.head())

    return df_out


# Run
if __name__ == "__main__":
    compute_deltaG()
