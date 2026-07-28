import os
import random
import pandas as pd
from rdkit import Chem
from rdkit.Chem import AllChem

# ===== Path setup =====
input_csv = "data/round1.csv"  
output_dir = "data/xyz"
os.makedirs(output_dir, exist_ok=True)

df = pd.read_csv(input_csv)

success, failed = 0, 0
failed_records = []

def _write_xyz(mol, smiles, mol_name):
    conf = mol.GetConformer()
    xyz_path = os.path.join(output_dir, f"{mol_name}.xyz")
    with open(xyz_path, "w") as f:
        f.write(f"{mol.GetNumAtoms()}\n")
        f.write(f"Generated from SMILES: {smiles}\n")
        for atom in mol.GetAtoms():
            pos = conf.GetAtomPosition(atom.GetIdx())
            f.write(f"{atom.GetSymbol():2s} {pos.x:.6f} {pos.y:.6f} {pos.z:.6f}\n")

def _try_embed(mol, **kwargs):
    """Try ETKDGv3 embedding with given keyword args. Returns (success, mol)."""
    m = Chem.RWMol(mol)
    p = AllChem.ETKDGv3()
    kwargs.setdefault("randomSeed", 42)
    kwargs.setdefault("maxIterations", 200)
    for k, v in kwargs.items():
        try:
            setattr(p, k, v)
        except AttributeError:
            if k == "maxIterations":
                try: setattr(p, "maxAttempts", v)
                except AttributeError: pass
            pass
    return AllChem.EmbedMolecule(m, p) != -1, m

def smiles_to_3d_xyz(smiles, mol_name, seed=42):
    try:
        mol = Chem.MolFromSmiles(smiles)
        if mol is None:
            return False, "invalid_smiles"

        mol = Chem.AddHs(mol)

        # Strategy 1: standard ETKDGv3
        p = AllChem.ETKDGv3()
        p.randomSeed = seed
        try: p.numThreads = 1
        except AttributeError: pass
        try: p.maxIterations = 200
        except AttributeError:
            try: p.maxAttempts = 200
            except AttributeError: pass
        
        ok = AllChem.EmbedMolecule(mol, p) != -1

        # Strategy 2: ETKDGv3 + random coords + no chirality enforcement
        if not ok:
            ok, mol = _try_embed(mol, randomSeed=seed, enforceChirality=False,
                                  useRandomCoords=True, maxIterations=200)

        # Strategy 3: same but also ignore smoothing failures (helps for
        # unusual atom types like hypervalent S, P, etc.)
        if not ok:
            ok, mol = _try_embed(mol, randomSeed=seed, enforceChirality=False,
                                  useRandomCoords=True, ignoreSmoothingFailures=True,
                                  maxIterations=200)

        # Strategy 4: fallback to 2D + small Z perturbation so ORCA can
        # optimize from a chemically reasonable (flat) starting geometry
        if not ok:
            AllChem.Compute2DCoords(mol)
            conf = mol.GetConformer()
            rng = random.Random(seed)
            for i in range(mol.GetNumAtoms()):
                pos = conf.GetAtomPosition(i)
                conf.SetAtomPosition(i, (pos.x, pos.y, rng.uniform(-0.1, 0.1)))
            ok = True  # 2D-based geometry written; ORCA will optimize it

        if not ok:
            return False, "embed_failed"

        # UFF optimization (best-effort; skip if atom types unsupported)
        try:
            AllChem.UFFOptimizeMolecule(mol, maxIters=200)
        except Exception:
            pass

        _write_xyz(mol, smiles, mol_name)
        return True, "ok"

    except Exception as e:
        return False, str(e)

# ============================
# 主循环：对所有 monomer 生成 xyz
# ============================
for idx, row in df.iterrows():
    monomer_id = str(row["PID"])
    smiles_raw = str(row["smiles"])

    # Cap open polymer ends: * → C
    smiles_fixed = smiles_raw.replace("*", "C")

    ok, reason = smiles_to_3d_xyz(smiles_fixed, monomer_id)
    if ok:
        success += 1
    else:
        failed += 1
        failed_records.append((monomer_id, smiles_raw, smiles_fixed, reason))

# ============================
# 保存失败记录
# ============================
if failed_records:
    fail_df = pd.DataFrame(
        failed_records,
        columns=["PID", "raw_SMILES", "converted_SMILES", "reason"]
    )
    fail_df.to_csv("data/failed_monomers.csv", index=False)

print("\n========== SUMMARY ==========")
print(f"Total monomers: {len(df)}")
print(f"Successfully generated: {success}")
print(f"Failed:                 {failed}")
if failed > 0:
    print("Saved failed SMILES to: data/failed_monomers.csv")
print("XYZ files saved in:", output_dir)
