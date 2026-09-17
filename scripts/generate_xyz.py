import os
import random
import multiprocessing as mp
import pandas as pd
from rdkit import Chem
from rdkit.Chem import AllChem, rdMolDescriptors

# Worker processes. Conformer generation is embarrassingly parallel (each molecule
# writes its own file) and per-molecule cost spans ~0.02 s to ~50 s, so the pool is
# fed with chunksize=1 to keep the slow floppy molecules from serialising a worker.
# Default is deliberately modest: this often gets run on a shared login node. For a
# full seed batch, submit it to a compute node and raise N_WORKERS.
def _default_workers():
    # sched_getaffinity respects cpuset/cgroup limits (what a scheduler actually
    # grants); cpu_count() reports the whole node and would oversubscribe.
    try:
        avail = len(os.sched_getaffinity(0))
    except AttributeError:
        avail = mp.cpu_count() or 1
    return max(1, min(8, avail))


N_WORKERS = int(os.environ.get("N_WORKERS", _default_workers()))

# ===== Path setup =====
input_csv = "data/round1.csv"
output_dir = "data/xyz"
os.makedirs(output_dir, exist_ok=True)

# Number of ETKDG conformers to generate per molecule before force-field ranking.
# ORCA's geometry optimisation is a LOCAL minimiser: it relaxes into the nearest
# minimum and never crosses a torsional barrier, so whichever conformer we hand it
# decides the final energy. Measured on the duplicate molecules in this dataset
# (same molecule submitted twice under two different SMILES spellings):
#   different starting conformer -> final Gibbs differed by 0.78 and 3.14 kcal/mol
#   identical starting conformer -> final Gibbs differed by 0.000 and 0.015 kcal/mol
# Screening the conformers here is FREE (RDKit + UFF, ~0.1-10 s) compared with the
# ~15 h of DFT each molecule costs, so there is no reason to hand ORCA a random one.
# Set N_CONFS=1 to restore the old single-conformer behaviour.
#
# The count is scaled by the number of rotatable bonds, because a fixed budget is
# statistically thin for floppy molecules: a 34-rotatable-bond PEG has ~3^34 torsional
# minima, so "best of 20" is barely better than one random draw (measured: it was the
# only molecule out of 100 where a fixed 20 lost to a single conformer). Thresholds
# follow the usual conformer-generation heuristic (Ebejer et al.).
N_CONFS = int(os.environ["N_CONFS"]) if "N_CONFS" in os.environ else None


def _n_confs_for(mol):
    """Conformers to generate: fixed if N_CONFS is set, else scaled by flexibility."""
    if N_CONFS is not None:
        return N_CONFS
    nrot = rdMolDescriptors.CalcNumRotatableBonds(mol)
    if nrot <= 7:
        return 20
    if nrot <= 12:
        return 50
    return 100

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

        # Canonicalise the atom ordering. ETKDG is deterministic for a fixed seed
        # only for a fixed atom ORDER, so the same molecule written two different
        # ways (e.g. *C(C(C*)C)(C)C vs *C(C*)C(C)(C)C, both 2,2,3-trimethylpentane)
        # would otherwise embed to different conformers and end up with different
        # DFT energies. Without this line the same molecule appearing in two source
        # databases gets two different ΔG values.
        _canon = Chem.MolFromSmiles(Chem.MolToSmiles(mol))
        if _canon is not None:
            mol = _canon

        mol = Chem.AddHs(mol)

        # Strategy 1: standard ETKDGv3
        n_confs = _n_confs_for(mol)

        p = AllChem.ETKDGv3()
        p.randomSeed = seed
        try: p.numThreads = 1
        except AttributeError: pass
        # maxIterations is a TOTAL embedding-attempt budget, so it must scale with the
        # number of conformers requested. Left at the single-embedding value of 200, a
        # floppy molecule asking for 100 conformers gets only ~3 -- i.e. the conformer
        # search silently does nothing for exactly the molecules that need it most.
        try: p.maxIterations = max(200, 30 * n_confs)
        except AttributeError:
            try: p.maxAttempts = max(200, 30 * n_confs)
            except AttributeError: pass

        # Generate n_confs conformers at once; the lowest-UFF-energy one is picked
        # below. Strategies 2-4 stay single-conformer: they only run for molecules
        # that fail normal embedding, where getting ANY geometry is the hard part.
        cids = list(AllChem.EmbedMultipleConfs(mol, numConfs=n_confs, params=p))
        ok = len(cids) > 0

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

        # UFF optimisation + conformer selection (best-effort; skip if atom types
        # unsupported). With several conformers, relax them all and keep the lowest
        # in energy -- this is what stops ORCA from being handed an arbitrary one.
        try:
            if mol.GetNumConformers() > 1:
                # maxIters must be generous: long floppy chains (e.g. multi-arm PEG)
                # leave a third of the conformers unconverged at 500, which makes the
                # reported energies -- and therefore the pick -- noise. UFF steps are
                # cheap relative to the embedding, so this costs almost nothing.
                res = AllChem.UFFOptimizeMoleculeConfs(mol, maxIters=5000)
                # res is [(not_converged, energy), ...] in conformer order
                order = [c.GetId() for c in mol.GetConformers()]
                best_id = min(zip(order, res), key=lambda t: t[1][1])[0]
                best_conf = Chem.Conformer(mol.GetConformer(best_id))
                mol.RemoveAllConformers()
                mol.AddConformer(best_conf, assignId=True)
            else:
                AllChem.UFFOptimizeMolecule(mol, maxIters=200)
        except Exception:
            # UFF unavailable for these atom types: fall back to the first
            # conformer, which is exactly the pre-existing behaviour.
            if mol.GetNumConformers() > 1:
                keep = Chem.Conformer(mol.GetConformer(mol.GetConformers()[0].GetId()))
                mol.RemoveAllConformers()
                mol.AddConformer(keep, assignId=True)

        _write_xyz(mol, smiles, mol_name)
        return True, "ok"

    except Exception as e:
        return False, str(e)

# ============================
# 主循环：对所有 monomer 生成 xyz (并行)
# ============================
def _one(task):
    """Worker: build one molecule's xyz. Returns a record for the summary."""
    monomer_id, smiles_raw = task
    # Cap open polymer ends: * → C
    smiles_fixed = smiles_raw.replace("*", "C")
    ok, reason = smiles_to_3d_xyz(smiles_fixed, monomer_id)
    return monomer_id, smiles_raw, smiles_fixed, ok, reason


def main():
    global success, failed
    tasks = [(str(r["PID"]), str(r["smiles"])) for _, r in df.iterrows()]
    print(f"Generating {len(tasks)} structures with {N_WORKERS} worker(s) ...",
          flush=True)

    if N_WORKERS > 1:
        pool = mp.Pool(N_WORKERS)
        # chunksize=1: per-molecule cost varies by ~1000x, so hand out work one at
        # a time rather than pre-slicing (a static split strands workers on the
        # floppy molecules that need 100 conformers).
        it = pool.imap_unordered(_one, tasks, chunksize=1)
    else:
        pool = None
        it = map(_one, tasks)

    try:
        for i, (monomer_id, smiles_raw, smiles_fixed, ok, reason) in enumerate(it, 1):
            if ok:
                success += 1
            else:
                failed += 1
                failed_records.append((monomer_id, smiles_raw, smiles_fixed, reason))
            if i % 100 == 0 or i == len(tasks):
                print(f"  {i}/{len(tasks)}  ok={success} failed={failed}", flush=True)
    finally:
        if pool is not None:
            pool.close()
            pool.join()


main()

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
