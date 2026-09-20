#!/usr/bin/env python3
"""
One-off pass over the whole pool: mark chemically implausible molecules and record
atom counts, then cache. Every selection component then draws from the SAME cleaned
universe.

Why this has to happen at POOL level rather than per component: in the v5 run the
feature/dopt/coverage components drew from a filtered candidate set while `random`
drew straight from the raw pool. That makes the random block's claim to be an
unbiased validation sample false -- it was unbiased with respect to a DIFFERENT
population than the one the other components sampled. It also let filtered-out
molecules back in through the random door (confirmed: a C=Si silene appeared in the
v5 random block that the filter would have rejected).

Rule set, with pool hit rates measured on a 120k sample before adopting:
    hypervalent iodine   0.022%    Si multiple bond   0.017%
    Ge multiple bond     0.008%    Sn multiple bond   0.004%
    C=I / I multiple     0.001%    aromatic O-O       0.001%
    2H/3H isotopes       0.000%
    -> ~0.05% of the pool, about 1,000 molecules.
DELIBERATELY NOT FILTERED: P=[!O] hits 0.854% (~17,700 molecules) because P=N is the
backbone of polyphosphazenes, a real polymer family. Dropping it would delete a
legitimate chemical class, not database noise.
"""
import os, numpy as np, pandas as pd, multiprocessing as mp
from rdkit import Chem, RDLogger
RDLogger.DisableLog("rdApp.*")

DIR = "/groups/tluo/FFE_Renzheng"
CACHE = os.path.join(DIR, "al_v2_cache")
N_WORKERS = int(os.environ.get("N_WORKERS", 16))

# Element whitelist. The downstream goal is a physics engine for an ORGANIC polymer
# generative model, which will never emit Ge or Sn, so modelling them is pure overhead
# -- and it is EXPENSIVE overhead, because the two highest-leverage components chase
# exactly the rarest features. Measured on seed_v6 (pool share -> component share):
#     Ge  0.08% -> dopt 4.00%, feature 9.69%   (50x / 121x enrichment)
#     Sn  0.11% -> dopt 4.40%, feature 7.75%   (40x /  70x)
# i.e. ~132 of 4015 seed molecules were spent on 0.19% of the chemistry, taken out of
# the budget of the components that carry 79% of the baseline's leverage.
# KEPT DELIBERATELY: Si (4.30%, polysiloxanes/silicones) and P (2.11%, polyphosphazenes)
# are backbone elements of real polymer families, not database noise. Dropping them
# because they look "inorganic" would delete two whole materials classes.
ALLOWED_ELEMENTS = set("C H O N F Cl Br S P Si I".split())

# Each rule is (SMARTS, a SMILES that MUST match, description). The positive example
# is not decoration: a SMARTS can parse cleanly and then never match anything, failing
# silently forever. That is exactly what "[2H,3H]" did here -- it parsed, reported zero
# isotope hits across 120k molecules, and let deuterated structures straight through.
# Split into [2H] and [3H] it works. Anything without a verifiable example is marked
# UNVERIFIED rather than quietly trusted.
BAD_RULES = [
    ("[Si]=,#*",              "C=[Si](C)C",           "Si multiple bond (disilene etc.)"),
    ("[Ge]=,#*",              "C=[Ge](C)C",           "Ge multiple bond (germene)"),
    ("[Sn]=,#*",              "C=[Sn](C)C",           "Sn multiple bond (stannene)"),
    ("[I]=,#*",               "CCCN(C)C(C)=I",        "I multiple bond (C=I)"),
    ("[I;v2,v3,v4,v5,v6,v7]", "CCC(C)(C)[IH](F)(F)F", "hypervalent iodine"),
    ("[2H]",                  "CC([2H])(Cl)C",        "deuterium"),
    ("[3H]",                  "CC([3H])(Cl)C",        "tritium"),
    ("[I][P]",                "CCC(C)P(=O)(Cl)I",     "I-P bond"),
    ("[I][N]",                "CCC(C)NI",             "I-N bond"),
    ("[c]:[p]",               "c1ccpcc1",             "aromatic C-P"),
    ("[o][o]",                None,                   "aromatic peroxide (UNVERIFIED)"),
]
BAD_SMARTS = [r[0] for r in BAD_RULES]


def _selftest():
    """Fail loudly if a rule cannot match its own positive example."""
    bad = []
    for sma, pos, desc in BAD_RULES:
        p = Chem.MolFromSmarts(sma)
        if p is None:
            bad.append(f"{desc}: SMARTS does not parse ({sma})")
            continue
        if pos is None:
            continue                      # no constructible example; kept on trust
        m = Chem.MolFromSmiles(pos)
        if m is None or not m.HasSubstructMatch(p):
            bad.append(f"{desc}: SMARTS {sma} does not match its example {pos}")
    if bad:
        raise SystemExit("SMARTS self-test failed:\n  " + "\n  ".join(bad))
    print(f"SMARTS self-test: {sum(1 for r in BAD_RULES if r[1])}/{len(BAD_RULES)} "
          f"rules verified against a positive example", flush=True)


_PATS = None


def _init():
    global _PATS
    _PATS = [Chem.MolFromSmarts(p) for p in BAD_SMARTS]


def _check(smi):
    """-> (ok, n_atoms). ok=False for unparseable or implausible."""
    m = Chem.MolFromSmiles(str(smi).replace("*", "C"))
    if m is None:
        return False, 0
    for p in _PATS:
        if p is not None and m.HasSubstructMatch(p):
            return False, 0
    for a in m.GetAtoms():
        if a.GetSymbol() not in ALLOWED_ELEMENTS:
            return False, 0
    return True, Chem.AddHs(m).GetNumAtoms()


if __name__ == "__main__":
    _selftest()
    meta = pd.read_csv(os.path.join(CACHE, "pool_meta.csv"))
    smis = meta["smiles"].astype(str).tolist()
    print(f"pool {len(smis):,} molecules, {N_WORKERS} workers", flush=True)
    with mp.Pool(N_WORKERS, initializer=_init) as pool:
        res = pool.map(_check, smis, chunksize=2000)
    ok = np.array([r[0] for r in res], dtype=bool)
    na = np.array([r[1] for r in res], dtype=np.int32)
    out = os.path.join(CACHE, "pool_clean.npz")
    print(f"element whitelist: {sorted(ALLOWED_ELEMENTS)}")
    np.savez_compressed(out, ok=ok, n_atoms=na)
    print(f"\nclean  : {ok.sum():,} ({ok.mean()*100:.3f}%)")
    print(f"dropped: {(~ok).sum():,} ({(~ok).mean()*100:.3f}%)")
    print(f"atoms  : median {np.median(na[ok]):.0f}  p99 {np.percentile(na[ok],99):.0f}  max {na[ok].max()}")
    print(f"wrote {out}")
