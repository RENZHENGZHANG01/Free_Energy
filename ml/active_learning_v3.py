#!/usr/bin/env python3
"""
Active Learning v3 — WHOLE-POOL uncertainty-weighted greedy max-coverage
========================================================================
v2 chose points by pure geometric coverage (facility-location on Morgan FPs) and
ignored the GIN ensemble. v3 keeps the coverage backbone (which prevents the
outlier-chasing that made v1 cover only 4%) but STEERS the budget with the
ensemble's per-molecule uncertainty.

Objective — choose N_SELECT anchors maximizing the blended-need-weighted pool mass
that gets a learnable neighbor within tau, over the WHOLE in-domain pool:

    gain(c) = Σ_{i ∈ (not-yet-captured ∩ Ball_tau(c))} w_i
    w_i     = (EPS_COV + cov_gap_i) · (EPS_UNC + û_i)

  û_i      = rank-normalized ensemble uncertainty in [0,1]  (robust to heavy tail)
  cov_gap_i= 1 if i is uncovered by the current train (NN-to-train < tau) else 0

Why whole-pool works without wasting budget on resolved old space: û_i is the
ensemble's epistemic state AFTER training on the current labeled set, so a
near-train point the model already nailed has low û (ignored), while a near-train
point it is still unsure about has high û (targeted — the user's ask). The cov_gap
factor keeps NEW chemical space in play even where the model is confidently wrong.
Pure coverage (v2) is the corner EPS_UNC→∞; this strictly generalizes it.

Uncertainty source: the 4 fresh <db>_predictions.csv (predict_dataset_gin.py on the
CURRENT gin_cv_models). Only uncertainty RANKS are used, so the known stale absolute
scaler in predict_dataset_gin.py does not matter here.

Env knobs: TAU, N_SELECT, K_CLUSTERS, CAND_SAMPLE, EPS_COV, EPS_UNC, N_WORKERS,
           POOL_CAP, FP_CHUNK, TRAIN_CSV, OUT_CSV, KMEANS_FIT_CAP.
Reuses al_v2_cache/pool_fp.npy (fingerprints identical to v2).
"""
import os, sys, time, warnings, multiprocessing as mp
import numpy as np, pandas as pd
from rdkit import Chem, RDLogger
from rdkit.Chem import AllChem, DataStructs
from sklearn.cluster import MiniBatchKMeans
from scipy.stats import rankdata
warnings.filterwarnings("ignore"); RDLogger.DisableLog("rdApp.*")

DIR = os.path.dirname(os.path.abspath(__file__))
NBITS = 2048
TAU         = float(os.environ.get("TAU", 0.4))
N_SELECT    = int(os.environ.get("N_SELECT", 200))
K_CLUSTERS  = int(os.environ.get("K_CLUSTERS", 5000))
CAND_SAMPLE = int(os.environ.get("CAND_SAMPLE", 300000))   # weighted-sample size for candidates
EPS_COV     = float(os.environ.get("EPS_COV", 0.15))       # floor so near-train stays eligible
EPS_UNC     = float(os.environ.get("EPS_UNC", 0.15))       # floor so coverage still breaks ties
N_WORKERS   = int(os.environ.get("N_WORKERS", max(1, mp.cpu_count())))
POOL_CAP    = int(os.environ.get("POOL_CAP", 0))
FP_CHUNK    = int(os.environ.get("FP_CHUNK", 20000))
TRAIN_CSV   = os.environ.get("TRAIN_CSV", "final_data_with_residual_deltaG.csv")
# Retired with the 2026-09 functional change along with the checkpoints trained on
# it; a residual target from one level of theory is not valid for another.
OUT_CSV     = os.environ.get("OUT_CSV", "al_v3_selected.csv")
CACHE = os.path.join(DIR, "al_v2_cache"); os.makedirs(CACHE, exist_ok=True)

FS5_ELEMENTS = {"C", "H", "N", "O", "S", "F", "Cl", "P", "Si", "Br", "I", "Ge", "Sn"}

DBS = [("PolyInfo", "Polyinfo/PolyInfo.csv",     "SMILES"),
       ("omics",    "omics/omics.csv",           "smiles_list"),
       ("OMG",      "OMG/OMG.csv",               "smiles"),
       ("PI1M",     "PI1M/PI1M_predictions.csv", "smiles_clean")]

# fresh ensemble-uncertainty predictions (predict_dataset_gin.py outputs)
PRED = {"PolyInfo": "Polyinfo/PolyInfo_predictions.csv",
        "omics":    "omics/omics_predictions.csv",
        "OMG":      "OMG/OMG_predictions.csv",
        "PI1M":     "PI1M/PI1M_predictions.csv"}


def _rss_gb():
    import resource
    return resource.getrusage(resource.RUSAGE_SELF).ru_maxrss / 1e6

def log(msg): print(f"  [{time.strftime('%H:%M:%S')}] (peakRSS {_rss_gb():.1f}G) {msg}", flush=True)


def _process(smi_orig):
    s = str(smi_orig).replace("*", "C")
    m = Chem.MolFromSmiles(s)
    if m is None: return (None, None)
    if Chem.GetFormalCharge(m) != 0: return ("__CHG__", None)
    for at in m.GetAtoms():
        if at.GetSymbol() not in FS5_ELEMENTS: return ("__OOV__", None)
    canon = Chem.MolToSmiles(m)
    bv = AllChem.GetMorganFingerprintAsBitVect(m, 2, nBits=NBITS)
    a = np.zeros(NBITS, np.uint8); DataStructs.ConvertToNumpyArray(bv, a)
    return (canon, a)

def _canon(s):
    m = Chem.MolFromSmiles(str(s))
    return Chem.MolToSmiles(m) if m else None


def build_pool():
    cache_fp = os.path.join(CACHE, f"pool_fp{'_cap'+str(POOL_CAP) if POOL_CAP else ''}.npy")
    cache_meta = os.path.join(CACHE, f"pool_meta{'_cap'+str(POOL_CAP) if POOL_CAP else ''}.csv")
    if os.path.exists(cache_fp) and os.path.exists(cache_meta):
        log(f"loading cached pool from {os.path.basename(cache_fp)}")
        meta = pd.read_csv(cache_meta)
        # The cache is the CHEMISTRY-FILTERED pool built by build_clean_pool.py: element
        # whitelist without Ge/Sn, plus eleven substructure rules for generative-database
        # artefacts (deuterium, hypervalent iodine, Si/Ge/Sn multiple bonds, I-P, I-N, ...).
        # The n_atoms column is the marker that the filter has been applied.
        if "n_atoms" not in meta.columns:
            raise SystemExit(
                f"{os.path.basename(cache_meta)} has no n_atoms column, so it predates the "
                "chemistry filter. Rebuild with build_clean_pool.py rather than letting the "
                "fallback path below regenerate an unfiltered pool.")
        return np.load(cache_fp), meta

    # NOTE: this fallback rebuilds with the ELEMENT filter only -- it does not apply the
    # substructure rules. It exists for bootstrapping from scratch, and the result must be
    # passed through build_clean_pool.py before use. Refuse to run it silently, because a
    # quietly-unfiltered pool would put all 17,398 filtered molecules back into acquisition
    # with nothing in the logs to say so.
    if os.environ.get("ALLOW_UNFILTERED_REBUILD") != "1":
        raise SystemExit(
            "No pool cache found. Rebuilding here produces an UNFILTERED pool (element "
            "whitelist only). Either restore al_v2_cache/pool_fp.npy + pool_meta.csv, or "
            "set ALLOW_UNFILTERED_REBUILD=1 and then run build_clean_pool.py on the result.")

    rows = []
    for name, path, col in DBS:
        s = pd.read_csv(os.path.join(DIR, path), usecols=[col])[col].dropna().astype(str).drop_duplicates()
        if POOL_CAP and len(s) > POOL_CAP:
            s = s.sample(POOL_CAP, random_state=42)
        for smi in s: rows.append((name, smi))
        log(f"{name:9s}: {len(s):>8d} raw unique SMILES")
    df_raw = pd.DataFrame(rows, columns=["db", "smiles"]).drop_duplicates("smiles").reset_index(drop=True)
    log(f"raw pool (dedup across DBs): {len(df_raw)}")

    log(f"parse + in-domain filter + fingerprint  ({N_WORKERS} workers) ...")
    db_arr = df_raw["db"].tolist(); smi_arr = df_raw["smiles"].tolist()
    if N_WORKERS > 1:
        pool = mp.Pool(N_WORKERS, maxtasksperchild=100000)
        it = pool.imap(_process, smi_arr, chunksize=2000)
    else:
        pool, it = None, (_process(s) for s in smi_arr)
    fp_all = np.empty((len(df_raw), NBITS), np.uint8)
    dbs, origs, canons, seen = [], [], [], set()
    n_unp = n_chg = n_oov = w = 0
    for i, (c, f) in enumerate(it):
        if f is None:
            n_unp += (c is None); n_chg += (c == "__CHG__"); n_oov += (c == "__OOV__"); continue
        if c in seen: continue
        seen.add(c); fp_all[w] = f; w += 1
        dbs.append(db_arr[i]); origs.append(smi_arr[i]); canons.append(c)
        if w % 200000 == 0: log(f"  ... {w} in-domain unique so far")
    if pool is not None: pool.close(); pool.join()
    fp = fp_all[:w]
    log(f"dropped: unparseable={n_unp}  charged={n_chg}  OOV-element={n_oov}")
    meta = pd.DataFrame({"db": dbs, "smiles": origs, "canon": canons})
    log(f"in-domain unique pool: {len(meta)}  (FP array {fp.nbytes/1e9:.2f} GB)")
    np.save(cache_fp, fp); meta.to_csv(cache_meta, index=False)
    return fp, meta


def load_uncertainty(meta):
    """Per-pool-row ensemble uncertainty (kcal/mol). Join by smiles_clean: both the
    pool and predict_dataset_gin.py derive smiles_clean = smiles.replace('*','C') from
    the SAME source strings, so string-matching is EXACT and needs no RDKit canon pass.
    Missing rows -> pool-median imputed."""
    frames = []
    for name, path in PRED.items():
        fpath = os.path.join(DIR, path)
        if not os.path.exists(fpath):
            log(f"  WARN: missing predictions {path} (uncertainty for {name} will impute)"); continue
        d = pd.read_csv(fpath, usecols=lambda c: c in ("smiles_clean", "Uncertainty_kcal_mol"))
        frames.append(d.dropna(subset=["smiles_clean", "Uncertainty_kcal_mol"]))
    if not frames:
        log("  WARN: no predictions found -> uniform uncertainty (v3 degrades to v2-like)")
        return np.ones(len(meta), np.float64), 1.0
    allp = pd.concat(frames, ignore_index=True).drop_duplicates("smiles_clean")
    umap = dict(zip(allp["smiles_clean"].astype(str), allp["Uncertainty_kcal_mol"].astype(float)))
    key = meta["smiles"].astype(str).str.replace("*", "C", regex=False)
    u = key.map(umap).values.astype(float)
    miss = np.isnan(u); u[miss] = np.nanmedian(u)
    log(f"  uncertainty joined by smiles_clean: {100*(1-miss.mean()):.1f}% matched, "
        f"{100*miss.mean():.1f}% imputed  ({len(umap)} pred keys)")
    return u, float(miss.mean())


def tani_block(A, B):
    A = A.astype(np.float32); B = B.astype(np.float32)
    inter = A @ B.T
    return inter / (A.sum(1)[:, None] + B.sum(1)[None, :] - inter + 1e-9)

def max_sim_to(ref_fp, pool_fp, chunk=FP_CHUNK):
    out = np.empty(len(pool_fp), np.float32)
    reff = ref_fp.astype(np.float32); rs = reff.sum(1)[None, :]
    for i in range(0, len(pool_fp), chunk):
        P = pool_fp[i:i+chunk].astype(np.float32)
        inter = P @ reff.T
        out[i:i+chunk] = (inter / (P.sum(1)[:, None] + rs - inter + 1e-9)).max(1)
    return out


def main():
    t0 = time.time()
    print("=" * 80)
    print(f"  Active Learning v3 — WHOLE-POOL uncertainty-weighted coverage @ tau={TAU}, select {N_SELECT}")
    print(f"  EPS_COV={EPS_COV} EPS_UNC={EPS_UNC} K={K_CLUSTERS} cand_sample={CAND_SAMPLE} "
          f"workers={N_WORKERS} pool_cap={POOL_CAP or 'FULL'}")
    print("=" * 80)

    pool_fp, meta = build_pool()
    N = len(meta)

    tr = pd.read_csv(os.path.join(DIR, TRAIN_CSV))["smiles_clean"].dropna().astype(str)
    train_fp = np.stack([f for _, f in (_process(s) for s in tr) if f is not None]).astype(np.uint8)
    log(f"train labeled set: {len(train_fp)} in-domain FPs")

    # geometric coverage gap vs current train
    log("computing covered_by_train ...")
    sim_tr = max_sim_to(train_fp, pool_fp)
    covered = sim_tr >= TAU
    cov_gap = (~covered).astype(np.float64)
    base_cov = covered.mean()
    log(f"train covers {base_cov*100:.1f}% of pool at tau={TAU}  |  uncovered = {int(cov_gap.sum())} "
        f"({100*cov_gap.mean():.1f}%)")

    # ensemble uncertainty -> rank-normalized û -> blended need weight w
    u_raw, miss = load_uncertainty(meta)
    u_hat = (rankdata(u_raw) - 1) / (len(u_raw) - 1)              # rank-normalize to [0,1] (robust)
    w = (EPS_COV + cov_gap) * (EPS_UNC + u_hat)
    log(f"weight w: mean={w.mean():.3f} p50={np.median(w):.3f} p99={np.percentile(w,99):.3f} "
        f"| û p50={np.median(u_hat):.2f} | uncovered share of Σw = {w[cov_gap>0].sum()/w.sum()*100:.1f}%")

    # ── candidates: weighted subsample (∝ w) -> cluster -> real-point medoids ──
    #    concentrates candidates where the weighted (uncertain × uncovered) mass is,
    #    incl. near-train-uncertain pockets (whole-pool eligibility).
    K = min(K_CLUSTERS, N)
    S = min(CAND_SAMPLE, N)
    log(f"sampling {S} candidate seeds ∝ w, clustering into K={K} medoids ...")
    samp = np.random.default_rng(0).choice(N, size=S, replace=False, p=w / w.sum())
    samp_X = pool_fp[samp].astype(np.float32)
    km = MiniBatchKMeans(n_clusters=K, n_init=3, random_state=0, batch_size=4096, max_iter=200).fit(samp_X)
    centers = km.cluster_centers_; slabels = km.labels_
    medoid_pool_idx = []
    for c in range(K):
        mem = np.where(slabels == c)[0]
        if len(mem) == 0: continue
        mloc = mem[int(np.argmin(np.linalg.norm(samp_X[mem] - centers[c], axis=1)))]
        medoid_pool_idx.append(int(samp[mloc]))
    medoid_pool_idx = np.array(sorted(set(medoid_pool_idx)))
    del samp_X
    log(f"candidate medoids: {len(medoid_pool_idx)}")

    # ── coverage sets over the WHOLE pool (which pool pts each medoid covers) ──
    log("computing candidate -> WHOLE-pool coverage sets ...")
    med_fp = pool_fp[medoid_pool_idx].astype(np.float32); med_sum = med_fp.sum(1)[None, :]
    nC = len(medoid_pool_idx)
    er, ec = [], []
    for i in range(0, N, FP_CHUNK):
        Pc = pool_fp[i:i+FP_CHUNK].astype(np.float32)
        sim = Pc @ med_fp.T
        sim = sim / (Pc.sum(1)[:, None] + med_sum - sim + 1e-9)
        rows, cols = np.where(sim >= TAU)
        er.append(rows + i); ec.append(cols)
    er = np.concatenate(er); ec = np.concatenate(ec)
    log(f"candidate->pool edges: {len(er)} (~{len(er)/max(nC,1):.0f}/candidate); grouping ...")
    order = np.argsort(ec, kind="stable"); ec_s, er_s = ec[order], er[order]
    bounds = np.searchsorted(ec_s, np.arange(nC + 1))
    cover_sets = [er_s[bounds[c]:bounds[c+1]] for c in range(nC)]

    # ── weighted greedy max-coverage (lazy) ──
    log(f"weighted greedy: selecting {N_SELECT} ...")
    import heapq
    captured = np.zeros(N, bool)
    gains = np.array([w[cs].sum() for cs in cover_sets])
    heap = [(-g, i) for i, g in enumerate(gains)]; heapq.heapify(heap)
    selected, sel_wgain, sel_ncov = [], [], []
    while len(selected) < N_SELECT and heap:
        neg_g, i = heapq.heappop(heap)
        cs = cover_sets[i]
        true_gain = float(w[cs][~captured[cs]].sum()) if len(cs) else 0.0
        if true_gain < -neg_g - 1e-9:
            heapq.heappush(heap, (-true_gain, i)); continue
        if true_gain <= 0: break
        n_new = int((~captured[cs]).sum())
        selected.append(int(medoid_pool_idx[i])); sel_wgain.append(true_gain); sel_ncov.append(n_new)
        captured[cs] = True
    log(f"selected {len(selected)} | captured weighted mass = {w[captured].sum()/w.sum()*100:.1f}% of Σw "
        f"| newly-covered pool pts = {int(captured.sum())} ({100*captured.sum()/N:.1f}%)")

    # ── report + write ──
    sel_fp = pool_fp[selected]
    new_cov = (max_sim_to(np.vstack([train_fp, sel_fp]), pool_fp) >= TAU).mean()
    pick_u = u_raw[selected]; pick_nn = sim_tr[selected]
    near_unc = int(((pick_nn >= TAU) & (u_hat[selected] >= 0.8)).sum())
    log(f"PROJECTED geometric combined coverage: {base_cov*100:.1f}% -> {new_cov*100:.1f}% "
        f"(+{(new_cov-base_cov)*100:.1f}%)")
    log(f"picks' uncertainty: median={np.median(pick_u):.2f} vs pool median={np.median(u_raw):.2f} kcal/mol "
        f"| {near_unc}/{len(selected)} picks are near-train(>=τ) AND high-û (whole-pool effect)")

    out = meta.iloc[selected].copy()
    out["uncertainty_kcal_mol"] = pick_u
    out["nn_to_train"] = pick_nn
    out["weighted_gain"] = sel_wgain
    out["n_covered"] = sel_ncov
    out = out.sort_values("weighted_gain", ascending=False).reset_index(drop=True)
    out.insert(0, "rank", np.arange(1, len(out) + 1))
    out.to_csv(os.path.join(DIR, OUT_CSV), index=False)
    log(f"wrote {OUT_CSV}  ({len(out)} molecules for DFT)")
    print(f"\n  DB breakdown of picks:\n{out['db'].value_counts().to_string()}")
    print(f"\n  Done in {time.time()-t0:.0f}s.  DONE!")


if __name__ == "__main__":
    import traceback
    try:
        main()
    except Exception:
        print("\n!!! EXCEPTION — full traceback:", flush=True)
        traceback.print_exc(); sys.stdout.flush(); raise
