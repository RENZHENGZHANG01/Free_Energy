from rdkit import Chem

def count_bonds(smiles):
    mol = Chem.MolFromSmiles(smiles)

    # Backbone bonds = only heavy atoms (no H)
    backbone_bonds = mol.GetNumBonds()

    # Total bonds including H
    mol_H = Chem.AddHs(mol)
    total_bonds = mol_H.GetNumBonds()

    return backbone_bonds, total_bonds


# Example monomers
monomers = {
    "M1": "CC(CC)CC",
    "M2": "CC(C)CCC(CC)CC",
    "M3": "CC(C)CCC(C)CCC(CC)CC",
    "M4": "CC(C)CCC(C)CCC(C)CCC(CC)CC",
    "M8": "CC(C)CCC(C)CCC(C)CCC(C)CCC(C)CCC(C)CCC(C)CCC(CC)CC"
}

for monomer_id, smi in monomers.items():
    backbone, total = count_bonds(smi)
    print(f"{monomer_id}: Backbone = {backbone}, Total = {total}, C–H = {total - backbone}")
