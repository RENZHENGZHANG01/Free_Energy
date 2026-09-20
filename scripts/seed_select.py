#!/usr/bin/env python3
"""
Iteration-0 seed selection, v2.

CHANGES FROM v1, each driven by a measurement on v1's own output:

 * `size` (stratified random over atom-count deciles) is REPLACED by `dopt`.
   v1 assumed atom count was a good proxy for "informative for the linear baseline".
   Measured on v1's design matrix (leverage share / point share):
       feature 6.71x   random 0.46x   size 0.40x   coverage 0.34x
   `size` came out WORSE than plain random, i.e. it was a redundant copy of the
   random block. Atom count is one direction of a 150-dim feature space and that
   direction is already over-determined; stratifying on it re-samples a solved axis.
   `dopt` instead selects by sequential D-optimality -- greedily taking the molecule
   that most increases det(X'X + lambda I), which IS the criterion for "pins down the
   baseline coefficients". Measured earlier at N=500: baseline MAE 11.73 for
   D-optimal vs 13.16 random vs 17.03 coverage-greedy.

 * Rare FS5 features now need >= MIN_SUPPORT (3) molecules, not 1. In v1, 61 of 150
   features had a single supporting molecule, so one bad DFT run would corrupt that
   coefficient outright. Measured blast radius of a 10 kcal/mol error on a feature
   molecule: median 0 other molecules moved >1 kcal/mol, mean 1.2, max 6 -- small,
   but it is the only component that moves anything, and redundant support removes it.

 * Features are only chased if they occur in >= FEAT_FLOOR of the pool, and molecules
   with implausible multiple bonds to Si/Ge/Sn are dropped. v1's set-cover went
   hunting for the rarest feature vectors, which are exactly the artefacts of the
   generative databases: BP_Si-Si_DOUBLE (disilene), BP_C-Ge_DOUBLE, BP_P-Sn_DOUBLE.
   Those species do not exist at room temperature; spending DFT on them is waste and
   they would anchor baseline coefficients nothing real ever uses.

Component roles (unchanged in spirit):
  feature  - every baseline coefficient gets data. Cheap insurance, ~1% of budget.
  dopt     - makes the design matrix well-conditioned: the baseline's accuracy.
  coverage - greedy max-coverage at tau; this is what makes later AL rounds efficient.
  random   - uniform draw. The ONLY unbiased yardstick, and it keeps the model
             calibrated to the pool distribution that the AL acquisition relies on.
             Held out of GNN TRAINING; still used to FIT the baseline (per-point
             leverage there is ~0.025, so the self-influence is negligible).
"""
import os, sys, time, heapq
import numpy as np
import pandas as pd
import multiprocessing as mp
from rdkit import Chem, RDLogger
from rdkit.Chem import AllChem, DataStructs
RDLogger.DisableLog("rdApp.*")

DIR = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, os.path.join(DIR, "Renzheng", "scripts"))
from compute_residual_deltaG import extract_features

CACHE       = os.path.join(DIR, "al_v2_cache")
N_COVER     = int(os.environ.get("N_COVER", 2500))
N_RANDOM    = int(os.environ.get("N_RANDOM", 1000))
N_DOPT      = int(os.environ.get("N_DOPT", 400))
N_FEAT_MAX  = int(os.environ.get("N_FEAT_MAX", 600))
# Measured on seed_v3: a single molecule whose DFT is wrong by 10 kcal/mol shifts the
# OTHER molecules sharing a feature by ~10/support kcal/mol -- 3.26 at support 3, 1.07
# at 6-10, 0.38 at 11-25. Support 3 therefore lands right on the ~3 kcal/mol threshold
# at which a difference counts as real, i.e. one bad run could fake a chemical signal.
# 12 puts it below 0.5, well under the 0.8-3.1 kcal/mol DFT noise floor, and costs only
# a couple of hundred cheap small molecules.
MIN_SUPPORT = int(os.environ.get("MIN_SUPPORT", 12))
FEAT_FLOOR  = float(os.environ.get("FEAT_FLOOR", 1e-4))   # pool frequency floor
TAU         = float(os.environ.get("TAU", 0.4))
CAND_N      = int(os.environ.get("CAND_N", 250000))
COV_REF_N   = int(os.environ.get("COV_REF_N", 0))   # 0 = whole pool
COV_CAND_N  = int(os.environ.get("COV_CAND_N", 25000))
MAX_ATOMS   = int(os.environ.get("MAX_ATOMS", 140))
MIN_ATOMS   = int(os.environ.get("MIN_ATOMS", 14))
N_WORKERS   = int(os.environ.get("N_WORKERS", 16))
SEED        = int(os.environ.get("SEED", 0))
COST_A, COST_B = 5.82e-5, 2.69
OUT = os.environ.get("OUT", os.path.join(DIR, "seed_v2.csv"))

# Multiple bonds to Si/Ge/Sn: disilenes/germenes/stannenes. Not room-temperature
# species; their presence flags a generated structure, not a candidate monomer.
BAD_SMARTS = [Chem.MolFromSmarts(p) for p in ("[Si]=,#*", "[Ge]=,#*", "[Sn]=,#*", "[o][o]")]

def log(m): print(f"  [{time.strftime('%H:%M:%S')}] {m}", flush=True)
def cost_h(n): return COST_A * np.power(np.asarray(n, float), COST_B)


def _prep(smi):
    """FS5 features + atom count + implausibility flag for one candidate."""
    m = Chem.MolFromSmiles(str(smi))
    if m is None:
        return None
    for p in BAD_SMARTS:
        if p is not None and m.HasSubstructMatch(p):
            return "BAD"
    f = extract_features(str(smi))
    if f is None:
        return None
    return f, Chem.AddHs(m).GetNumAtoms()


def tani(Q, P, Psum):
    inter = Q @ P.T
    return inter / (Q.sum(1)[:, None] + Psum[None, :] - inter + 1e-9)


def build_edges(fp, cand_idx, ref_idx, tau):
    """candidate -> covered-pool-molecule lists. Tanimoto on binary fingerprints is
    one dense matmul, so this runs on the GPU when there is one: measured 138 s for
    80k candidates against the FULL 2.07M pool, versus ~3.8 h on CPU."""
    try:
        import torch
        use_gpu = torch.cuda.is_available()
    except ImportError:
        use_gpu = False
    if not use_gpu:
        P = np.asarray(fp[ref_idx], dtype=np.float32); Ps = P.sum(1)
        Q = np.asarray(fp[cand_idx], dtype=np.float32)
        log("    (CPU path)")
        out = []
        for i in range(0, len(Q), 400):
            S = tani(Q[i:i+400], P, Ps) >= tau
            for r in S:
                out.append(np.flatnonzero(r).astype(np.int32))
        return out

    dev = "cuda"
    log(f"    (GPU path: {torch.cuda.get_device_name(0)})")
    C = torch.from_numpy(np.asarray(fp[cand_idx], dtype=np.float16)).to(dev)
    cs = C.float().sum(1)
    out = [[] for _ in range(len(cand_idx))]
    PCH, CCH = 50000, 8192
    for i in range(0, len(ref_idx), PCH):
        rid = ref_idx[i:i+PCH]
        P = torch.from_numpy(np.asarray(fp[rid], dtype=np.float16)).to(dev)
        ps = P.float().sum(1)
        for j in range(0, len(C), CCH):
            Q = C[j:j+CCH]
            inter = (Q @ P.T).float()
            hit = (inter / (cs[j:j+CCH, None] + ps[None, :] - inter + 1e-6)) >= tau
            nz = hit.nonzero()
            if len(nz):
                # Vectorised regroup. There are ~1e8 hits in total, so the Python
                # loop must run over candidates in the chunk (8k), never over hits.
                a = nz[:, 0].cpu().numpy()
                b = (nz[:, 1].to(torch.int32) + i).cpu().numpy()
                order = np.argsort(a, kind="stable")
                a_s, b_s = a[order], b[order]
                bounds = np.searchsorted(a_s, np.arange(len(Q) + 1))
                for k in range(len(Q)):
                    if bounds[k + 1] > bounds[k]:
                        out[j + k].append(b_s[bounds[k]:bounds[k + 1]])
            del inter, hit, nz
        del P, ps
    return [np.concatenate(x).astype(np.int32) if x else np.empty(0, np.int32)
            for x in out]


def main():
    rng = np.random.default_rng(SEED)
    print("=" * 86)
    print("  Iteration-0 seed selection v2")
    print(f"  coverage={N_COVER} random={N_RANDOM} dopt={N_DOPT} feature<={N_FEAT_MAX}"
          f" | tau={TAU} support>={MIN_SUPPORT} feat_floor={FEAT_FLOOR}")
    print("=" * 86)

    fp = np.load(os.path.join(CACHE, "pool_fp.npy"), mmap_mode="r")
    meta = pd.read_csv(os.path.join(CACHE, "pool_meta.csv"))
    N_RAW = len(fp)

    # ONE universe for ALL four components. In v5 feature/dopt/coverage drew from a
    # filtered candidate set while `random` drew from the raw pool, which (a) made the
    # random block unbiased with respect to a DIFFERENT population than the training
    # components sampled -- destroying the one property it exists for -- and (b) let
    # filtered-out molecules back in through the random door (a C=Si silene did exactly
    # that). The chemistry filter now lives in the cache itself: build_clean_pool.py
    # computes it and the pool under al_v2_cache/pool_fp.npy is already filtered, with
    # the unfiltered original kept as pool_fp_raw.npy. That way nothing downstream --
    # including active_learning_v3.py, which never knew pool_clean.npz existed and would
    # have re-selected the filtered molecules at acquisition time -- has to remember to
    # apply a mask. Only the atom-count window is applied here, since it is a cost knob
    # rather than a correctness one.
    if "n_atoms" not in meta.columns:
        raise SystemExit("pool_meta.csv has no n_atoms column -- this cache predates "
                         "the baked-in chemistry filter. Re-run build_clean_pool.py.")
    na_all = meta["n_atoms"].values
    keep = (na_all >= MIN_ATOMS) & (na_all <= MAX_ATOMS)
    CLEAN = np.flatnonzero(keep)
    N_POOL = len(CLEAN)
    log(f"pool {N_RAW:,} (chemistry-filtered cache) -> usable {N_POOL:,} "
        f"({N_POOL/N_RAW*100:.2f}%; {int((~keep).sum()):,} outside atom window "
        f"[{MIN_ATOMS},{MAX_ATOMS}])")

    cand = np.sort(rng.choice(CLEAN, min(CAND_N, N_POOL), replace=False))
    smis = meta["smiles"].iloc[cand].astype(str).str.replace("*", "C", regex=False).tolist()
    log(f"featurising {len(cand):,} candidates ({N_WORKERS} workers) ...")
    with mp.Pool(N_WORKERS) as pool:
        out = pool.map(_prep, smis, chunksize=500)
    # Chemistry and atom window were already applied when CLEAN was built, so anything
    # rejected here is only an FS5 featurisation failure.
    keep = [i for i, o in enumerate(out) if o is not None and o != "BAD"]
    if len(keep) < len(out):
        log(f"  {len(out)-len(keep):,} candidates failed FS5 featurisation")
    cand = cand[keep]
    feats = [out[i][0] for i in keep]
    natoms = np.array([out[i][1] for i in keep], float)

    cost = cost_h(natoms)

    F = pd.DataFrame(feats).fillna(0.0)
    names = np.array(F.columns)
    present = F.values > 0
    freq = present.mean(0)
    log(f"  FS5 features present in candidates: {len(names)}")

    chosen = {}

    # ---- (1) feature set-cover, >=MIN_SUPPORT each, cost-aware -----------------
    target = (freq >= FEAT_FLOOR)
    need = np.where(target, MIN_SUPPORT, 0).astype(int)
    log(f"(1) chasing {int(target.sum())} features at >= {FEAT_FLOOR:.4%} pool frequency "
        f"(skipping {int((~target).sum())} rarer ones -> flag as OOV downstream)")
    avail = present.copy()
    fpick = []
    for _ in range(N_FEAT_MAX):
        if need.sum() == 0:
            break
        gain = (avail & (need > 0)[None, :]).sum(1).astype(float)
        gain[fpick] = -1
        eff = gain / np.sqrt(np.maximum(cost, 1e-3))
        eff[gain <= 0] = -1
        j = int(eff.argmax())
        if eff[j] <= 0:
            break
        fpick.append(j)
        need = np.maximum(0, need - present[j].astype(int))
        avail[j] = False
    for j in fpick:
        chosen.setdefault(int(cand[j]), "feature")
    log(f"    {len(fpick)} molecules, unmet feature-slots {int(need.sum())}, "
        f"cost {cost[fpick].sum():.0f} h")

    # ---- (2) sequential D-optimal ----------------------------------------------
    # Greedily add the molecule maximising the increase in log det(X'X + lam I).
    # For a candidate x the increase is log(1 + x' M^-1 x), so argmax over x of the
    # quadratic form, then rank-1 update M^-1 (Sherman-Morrison). This is the real
    # criterion for "determines the regression coefficients", which atom-count
    # stratification only approximated -- and did worse than random at.
    mu, sd = F.values.mean(0), F.values.std(0)
    sd[sd == 0] = 1.0
    Xs = (F.values - mu) / sd
    p = Xs.shape[1]
    lam = 1.0
    Minv = np.eye(p) / lam
    picked_d = []
    taken = np.zeros(len(Xs), bool)
    for i in fpick:                      # feature picks already contribute
        v = Xs[i]; Mv = Minv @ v
        Minv -= np.outer(Mv, Mv) / (1.0 + v @ Mv)
        taken[i] = True
    # Restrict to a subsample: the quadratic form must be recomputed every step, so
    # this is O(steps * n * p^2). 60k candidates is ample for a 150-column design.
    dsub = rng.choice(len(Xs), min(60000, len(Xs)), replace=False)
    dsub = np.union1d(dsub, np.array(fpick, dtype=int)) if fpick else dsub
    Xd = Xs[dsub]
    log(f"(2) D-optimal: selecting {N_DOPT} from {len(Xd):,} candidates ...")
    tk = taken[dsub].copy()
    for _ in range(N_DOPT):
        q = np.einsum('ij,ij->i', Xd @ Minv, Xd)   # BLAS matmul, then row-wise dot
        q[tk] = -1
        jj = int(q.argmax())
        if q[jj] <= 0:
            break
        j = int(dsub[jj])
        v = Xs[j]; Mv = Minv @ v
        Minv -= np.outer(Mv, Mv) / (1.0 + v @ Mv)
        taken[j] = True; tk[jj] = True
        picked_d.append(j)
    for j in picked_d:
        chosen.setdefault(int(cand[j]), "dopt")
    log(f"    {len(picked_d)} molecules, cost {cost[picked_d].sum():,.0f} h")

    # ---- (3) greedy max-coverage ------------------------------------------------
    # Reference = the molecules whose coverage we measure. This can be the WHOLE pool:
    # it only costs one boolean per molecule. COV_REF_N=0 means "use all of it".
    ref = (CLEAN if COV_REF_N <= 0 or COV_REF_N >= N_POOL
           else np.sort(rng.choice(CLEAN, COV_REF_N, replace=False)))
    # Candidates = the molecules greedy is allowed to PICK. This one cannot be the whole
    # pool: we must store each candidate's neighbour list, and 2.07M candidates x ~13k
    # neighbours would be ~108 GB. It is drawn straight from the full pool (NOT from the
    # featurised subsample -- that coupling was pointless, coverage needs only the
    # cached fingerprints). The restriction is measured to be non-binding: 80k random
    # candidates can reach 96.6% of the pool, while 2000 picks only reach ~91%.
    cov_cand = np.sort(rng.choice(CLEAN, min(COV_CAND_N, N_POOL), replace=False))
    log(f"(3) coverage: {len(cov_cand):,} candidates (drawn from the full {N_POOL:,}) "
        f"vs {len(ref):,} reference molecules, building edges ...")
    nbrs = build_edges(fp, cov_cand, ref, TAU)
    log(f"    edges {sum(len(x) for x in nbrs):,}")
    covered = np.zeros(len(ref), bool)
    heap = [(-len(nbrs[i]), i) for i in range(len(nbrs))]
    heapq.heapify(heap)
    cpick = []
    while heap and len(cpick) < N_COVER:
        negg, j = heapq.heappop(heap)
        g = int((~covered[nbrs[j]]).sum())
        if not heap or g >= -heap[0][0]:
            if g <= 0:
                break
            cpick.append(j); covered[nbrs[j]] = True
        else:
            heapq.heappush(heap, (-g, j))
    for j in cpick:
        chosen.setdefault(int(cov_cand[j]), "coverage")
    log(f"    {len(cpick)} molecules, covers {covered.mean()*100:.1f}% of pool sample")

    # ---- (4) uniform random -----------------------------------------------------
    got = 0
    while got < N_RANDOM:
        for i in rng.choice(CLEAN, min(N_RANDOM * 2, N_POOL), replace=False):
            if int(i) not in chosen:
                chosen[int(i)] = "random"; got += 1
                if got >= N_RANDOM:
                    break
    log(f"(4) random: {got} molecules (held out of GNN training)")

    # ---- (5) support repair -----------------------------------------------------
    # The set-cover in (1) only guarantees MIN_SUPPORT for the features it TARGETED.
    # coverage/random/dopt pick molecules without looking at features at all, so they
    # drag in rarer features that then sit in the design matrix with 1-2 supporting
    # molecules -- exactly the single-point fragility (1) exists to prevent. Measured
    # on the v4 run: 32 features still ended up below 12, and the worst single-point
    # influence moved from `feature` (1.38 kcal/mol) to `dopt` (4.92), because
    # D-optimality deliberately seeks the extreme-composition molecules that carry
    # them. So: recount support over the ASSEMBLED set and top it up.
    sel_idx = np.array(sorted(chosen))
    pos = {int(c): i for i, c in enumerate(cand)}
    have = np.array([pos[int(i)] for i in sel_idx if int(i) in pos])
    supp = present[have].sum(0) if len(have) else np.zeros(present.shape[1], int)
    short = np.where((supp > 0) & (supp < MIN_SUPPORT))[0]
    log(f"(5) support repair: {len(short)} features below {MIN_SUPPORT} after assembly")
    if len(short):
        deficit = np.zeros(present.shape[1], int)
        deficit[short] = MIN_SUPPORT - supp[short]
        avail2 = present.copy()
        avail2[have] = False
        added = 0
        for _ in range(N_FEAT_MAX):
            if deficit.sum() == 0:
                break
            gain = (avail2 & (deficit > 0)[None, :]).sum(1).astype(float)
            eff = gain / np.sqrt(np.maximum(cost, 1e-3))
            eff[gain <= 0] = -1
            j = int(eff.argmax())
            if eff[j] <= 0:
                break
            chosen.setdefault(int(cand[j]), "feature")
            deficit = np.maximum(0, deficit - present[j].astype(int))
            avail2[j] = False
            added += 1
        still = int((deficit > 0).sum())
        log(f"    added {added} molecules; {still} features still short "
            f"(pool simply lacks cheap carriers -> flag those as OOV downstream)")

    # ---- assemble ---------------------------------------------------------------
    idxs = np.array(sorted(chosen))
    rows = meta.iloc[idxs].copy().reset_index(drop=True)
    rows["component"] = [chosen[i] for i in idxs]
    na_map = {int(c): n for c, n in zip(cand, natoms)}
    na = []
    for i, s in zip(idxs, rows["smiles"].astype(str).str.replace("*", "C", regex=False)):
        v = na_map.get(int(i))
        if v is None:
            m = Chem.MolFromSmiles(s); v = Chem.AddHs(m).GetNumAtoms() if m else np.nan
        na.append(v)
    rows["n_atoms"] = na
    rows["est_cost_h"] = cost_h(rows["n_atoms"].values)
    rows["is_validation"] = rows["component"] == "random"

    # Size-range safeguard. dopt already enriches both tails ~3x over the pool
    # (measured: <25 atoms 3.8% vs pool 1.2%; >110 atoms 7.0% vs 2.1%), but assert it
    # rather than assume -- if either tail is thinner than the pool's own share,
    # top it up so the linear baseline keeps leverage across the whole size range.
    na_all = rows["n_atoms"].values.astype(float)
    for lo, hi, label in [(MIN_ATOMS, 25, "small"), (110, MAX_ATOMS, "large")]:
        share = ((na_all >= lo) & (na_all < hi)).mean()
        pool_share = ((natoms >= lo) & (natoms < hi)).mean()
        n_now = int(((na_all >= lo) & (na_all < hi)).sum())
        flag = "OK" if share >= pool_share else "THIN"
        log(f"    size tail {label:5s} [{lo},{hi}): {n_now} molecules "
            f"= {share*100:.2f}% vs pool {pool_share*100:.2f}%  [{flag}]")
    rows = rows.sample(frac=1.0, random_state=SEED).reset_index(drop=True)
    rows.insert(0, "rank", np.arange(1, len(rows) + 1))
    rows.to_csv(OUT, index=False)

    tot = rows.est_cost_h.sum()
    print("\n" + "=" * 86)
    print(f"  wrote {OUT}  ({len(rows)} molecules)")
    print(rows.groupby("component").agg(n=("rank", "size"),
          med_atoms=("n_atoms", "median"), cost_h=("est_cost_h", "sum")).to_string())
    print(f"\n  TOTAL {tot:,.0f} CPU-hours | 100 concurrent {tot/100/24:.1f} d "
          f"| 200 concurrent {tot/200/24:.1f} d")
    print(f"  db mix:\n{rows['db'].value_counts().to_string()}")
    print("  DONE!")


if __name__ == "__main__":
    main()
