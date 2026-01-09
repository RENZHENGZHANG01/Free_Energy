import os
import pandas as pd
from rdkit import Chem
from rdkit.Chem import AllChem

# ===== Path setup =====
input_csv = "data/PI1070.csv"  
output_dir = "data/xyz"           # 统一放到 data/xyz 目录
os.makedirs(output_dir, exist_ok=True)

df = pd.read_csv(input_csv)

success, failed = 0, 0
failed_records = []  # 存失败的 monomer_ID & smiles

def smiles_to_3d_xyz(smiles, mol_name, seed=42):
    try:
        mol = Chem.MolFromSmiles(smiles)
        if mol is None:
            return False

        mol = Chem.AddHs(mol)

        # 用 ETKDGv3 生成 3D 构型
        params = AllChem.ETKDGv3()
        params.randomSeed = seed
        AllChem.EmbedMolecule(mol, params)

        # 用 UFF 优化几何
        AllChem.UFFOptimizeMolecule(mol, maxIters=200)

        conf = mol.GetConformer()
        xyz_path = os.path.join(output_dir, f"{mol_name}.xyz")

        with open(xyz_path, "w") as f:
            # 第一行：原子数
            f.write(f"{mol.GetNumAtoms()}\n")
            # 第二行：注释行
            f.write(f"Generated from SMILES: {smiles}\n")
            # 后面每一行：元素 + 坐标
            for atom in mol.GetAtoms():
                pos = conf.GetAtomPosition(atom.GetIdx())
                f.write(f"{atom.GetSymbol():2s} {pos.x:.6f} {pos.y:.6f} {pos.z:.6f}\n")

        return True

    except Exception:
        return False

# ============================
# 主循环：对所有 monomer 生成 xyz
# ============================
for idx, row in df.iterrows():
    monomer_id = str(row["monomer_ID"])
    smiles_raw = str(row["smiles"])

    # 封端：* → C 
    smiles_fixed = smiles_raw.replace("*", "C")

    ok = smiles_to_3d_xyz(smiles_fixed, monomer_id)
    if ok:
        success += 1
    else:
        failed += 1
        failed_records.append((monomer_id, smiles_raw, smiles_fixed))

# ============================
# 保存失败记录
# ============================
if failed_records:
    fail_df = pd.DataFrame(
        failed_records,
        columns=["monomer_ID", "raw_SMILES", "converted_SMILES"]
    )
    fail_df.to_csv("data/failed_monomers.csv", index=False)

print("\n========== SUMMARY ==========")
print(f"Total monomers: {len(df)}")
print(f"Successfully generated: {success}")
print(f"Failed:                 {failed}")
if failed > 0:
    print("Saved failed SMILES to: data/failed_monomers.csv")
print("XYZ files saved in:", output_dir)
