import os
import sys
import pandas as pd
from rdkit import Chem
from rdkit.Chem import Descriptors

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import orca_settings as S


def check_level_consistency(df, atom_ref_csv):
    """Refuse to combine molecular and atomic energies from different methods.

    Delta_G = G(molecule) - SUM n_i * G(atom_i) is a difference of ABSOLUTE
    energies. If the two sides came from different functionals, basis sets or
    grids, nothing cancels and every Delta_G is wrong by tens of kcal/mol per
    atom -- with no symptom in the output. The 2026-09 switch from
    B3LYP-D3BJ/def2-TZVP to wB97X-D3/def2-TZVP made that failure reachable by
    simply re-running this script against an old thermo table, which is exactly
    the kind of silent corruption worth one loud check.

    extract_thermo.py writes a 'level' column holding the keyword line ORCA
    echoed for each molecule. A table without that column predates the check and
    cannot be verified, so it is refused unless ALLOW_UNVERIFIED_LEVEL=1.
    """
    expect = " ".join(S.FREQ.lower().split())

    if "level" not in df.columns:
        if os.environ.get("ALLOW_UNVERIFIED_LEVEL") == "1":
            print("  WARNING: input has no 'level' column; level of theory NOT verified.")
            return
        raise SystemExit(
            "\n  ABORT: the input thermo table has no 'level' column, so the level of\n"
            "  theory it was computed at cannot be checked against the atomic\n"
            f"  references in {atom_ref_csv}.\n\n"
            f"  Expected: {S.FREQ}\n\n"
            "  Regenerate it with the current extract_thermo.py, which records the\n"
            "  level of theory per molecule:\n"
            "      python scripts/extract_thermo.py --expect-level\n\n"
            "  If you are certain the input predates that column AND was computed at\n"
            "  the level above, set ALLOW_UNVERIFIED_LEVEL=1.")

    found = sorted({" ".join(str(x).lower().split()) for x in df["level"].dropna()} - {""})
    wrong = [lv for lv in found if lv != expect]
    if wrong:
        raise SystemExit(
            "\n  ABORT: molecular energies and atomic references are at different\n"
            "  levels of theory. Delta_G would be meaningless.\n\n"
            f"  Expected: {S.FREQ}\n"
            + "".join(f"  Found:    {lv}\n" for lv in wrong)
            + "\n  Either recompute the molecules at the expected level, or point\n"
              "  orca_settings.py at the level the molecules actually used and\n"
              "  regenerate atom_ref.csv to match.")
    print(f"  Level of theory verified: {S.FREQ}")

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

    # --- Input: written directly by extract_thermo.py, no manual merge step ---
    root = os.path.dirname(script_dir)
    input_csv = os.environ.get("THERMO_CSV",
                               os.path.join(root, "data", "deltaG_raw.csv"))
    atom_ref_csv = os.path.join(script_dir, "atom_ref.csv")

    print("Reading input file:", input_csv)
    if not os.path.exists(input_csv):
        raise SystemExit(
            f"\n  ABORT: {input_csv} not found.\n"
            "  Produce it from the ORCA frequency outputs first:\n"
            "      python scripts/extract_thermo.py --expect-level")
    df = pd.read_csv(input_csv)

    check_level_consistency(df, atom_ref_csv)

    # Drop the runs extract_thermo.py flagged as unusable (no normal termination,
    # SCF failure, no frequencies). Their Gibbs_Eh is already NaN; dropping them
    # here keeps the count honest instead of emitting a row of NaNs.
    if "usable" in df.columns:
        n_bad = int((~df["usable"].astype(bool)).sum())
        if n_bad:
            print(f"  Dropping {n_bad} unusable calculation(s) flagged by extract_thermo.py")
            df = df[df["usable"].astype(bool)].reset_index(drop=True)
    if "smiles" not in df.columns:
        raise SystemExit(
            f"\n  ABORT: {input_csv} has no 'smiles' column, so atoms cannot be\n"
            "  counted. Re-run extract_thermo.py with --smiles-csv pointing at the\n"
            "  campaign input CSV.")
    n_nosmi = int(df["smiles"].isna().sum())
    if n_nosmi:
        print(f"  Dropping {n_nosmi} molecule(s) with no SMILES")
        df = df[df["smiles"].notna()].reset_index(drop=True)

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
