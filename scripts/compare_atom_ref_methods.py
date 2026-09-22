#!/usr/bin/env python3
"""
Quantify what switching the atomic reference set from one level of theory to another
does -- e.g. legacy_b3lyp (B3LYP-D3BJ/def2-TZVP) -> production (wB97X-D3/def2-TZVP).

Two very different questions, answered separately:

(1) Per-ATOM reference energies differ by a LOT (different functional => different
    absolute energy). That number alone is NOT the error you would fix -- it mostly
    cancels in Delta_G, because Delta_G subtracts the atoms from the molecule.

(2) What actually matters is the effect on Delta_G, and specifically on the RESIDUAL
    (the training signal). Since
        Delta_G = G(mol) - SUM n_i E_atom(i)
    changing only the atom references shifts Delta_G by a per-element constant times
    atom count -- which the composition-level Ridge baseline absorbs EXACTLY. So an
    atom-reference change alone cannot change the residual at all. Only recomputing
    the MOLECULES with the new functional changes the science.

This script therefore reports:
  a) the raw per-atom differences (both E_elec and G)
  b) the induced Delta_G shift for the real dataset composition (per molecule)
  c) an explicit statement of what is/ isn't affected, so the cost/benefit is clear.

The two sets are named by their subdirectory under data/atom_ref/ (as written by
generate_atom_ref_inp.py and parsed by parse_atom_ref.py).

Usage:
  python compare_atom_ref_methods.py [--old legacy_b3lyp] [--new production] \
                                     [--dataset ../../final_data_with_residual_deltaG.csv]
"""
import os, argparse, csv
import pandas as pd, numpy as np

H2K = 627.509474
HERE = os.path.dirname(os.path.abspath(__file__))


def load(tag):
    p = os.path.join(HERE, "..", "data", "atom_ref", f"atom_ref_{tag}.csv")
    if not os.path.exists(p):
        return None
    d = pd.read_csv(p)
    return {r["atom"]: r for _, r in d.iterrows()}


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--dataset", default=os.path.join(HERE, "..", "..",
                                                      "final_data_with_residual_deltaG.csv"))
    ap.add_argument("--old", default="legacy_b3lyp", help="baseline set name")
    ap.add_argument("--new", default="production", help="comparison set name")
    args = ap.parse_args()

    b, w = load(args.old), load(args.new)
    if not b or not w:
        have = sorted(os.path.basename(f)[len("atom_ref_"):-4]
                      for f in __import__("glob").glob(
                          os.path.join(HERE, "..", "data", "atom_ref", "atom_ref_*.csv")))
        print(f"  Need both atom_ref_{args.old}.csv and atom_ref_{args.new}.csv "
              "(run parse_atom_ref.py after the ORCA job).")
        print(f"  Sets available: {have or 'none'}")
        return

    print("=" * 88)
    print(f"  (a) PER-ATOM reference energies: {args.old}  vs  {args.new}")
    print("=" * 88)
    print(f"  {'atom':5s} {'mult':>4s} | {args.old[:15]:>15s} {args.new[:15]:>15s} "
          f"{'diff(Eh)':>11s} {'diff(kcal/mol)':>14s}")
    print("  " + "-" * 84)
    rows = []
    for a in b:
        if a not in w:
            print(f"  {a:5s}  -- missing in wB97X set"); continue
        eb, ew = b[a]["E_elec_Eh"], w[a]["E_elec_Eh"]
        gb, gw = b[a]["G_Eh"], w[a]["G_Eh"]
        d = ew - eb
        rows.append((a, eb, ew, d, gb, gw))
        print(f"  {a:5s} {int(b[a]['mult']):>4d} | {eb:15.6f} {ew:15.6f} "
              f"{d:11.6f} {d*H2K:14.1f}")
    print("\n  NOTE: these absolute shifts are huge but LARGELY IRRELEVANT on their own --")
    print("        they cancel against the molecular energy in Delta_G.")

    # (b) induced Delta_G shift for the real dataset, using actual atom counts
    ds = os.path.abspath(args.dataset)
    if os.path.exists(ds):
        df = pd.read_csv(ds)
        cnt_cols = {c[2:]: c for c in df.columns if c.startswith("N_") and c[2:] in b}
        if cnt_cols:
            print("\n" + "=" * 88)
            print("  (b) Induced Delta_G shift IF ONLY the atom references were swapped")
            print("      (molecules still B3LYP) -- per-element constant x atom count")
            print("=" * 88)
            shift = np.zeros(len(df))
            for el, col in cnt_cols.items():
                d = w[el]["E_elec_Eh"] - b[el]["E_elec_Eh"]
                shift += df[col].fillna(0).values * d
            shift_k = shift * H2K
            print(f"  molecules: {len(df)}   elements matched: {sorted(cnt_cols)}")
            print(f"  Delta_G shift  mean={shift_k.mean():+.1f}  sd={shift_k.std():.1f}  "
                  f"min={shift_k.min():+.1f}  max={shift_k.max():+.1f} kcal/mol")
            resid = df["Delta_G_residual"].dropna()
            print(f"  (for scale: current residual sd = {resid.std():.2f} kcal/mol)")
            print("\n  >>> This shift is EXACTLY a linear function of atom counts, and the Ridge")
            print("      baseline uses those very counts as features -> it is absorbed with ZERO")
            print("      change to the residual target. Swapping atom refs alone buys NOTHING.")
            print("      Real gains require recomputing the MOLECULES with wB97X-D3/def2-TZVPPD.")
    else:
        print(f"\n  (dataset not found at {ds}; skipped part b)")

    # (c) thermal terms, for reference
    print("\n" + "=" * 88)
    print("  (c) G - E_elec per atom (thermal/entropic term; single atoms have no ZPE)")
    print("=" * 88)
    print(f"  {'atom':5s} | {'B3LYP G-E':>12s} {'wB97X G-E':>12s}   (kcal/mol)")
    for a, eb, ew, d, gb, gw in rows:
        if gb is None or gw is None or (isinstance(gb, float) and np.isnan(gb)):
            print(f"  {a:5s} |  (no Freq thermochemistry parsed)"); continue
        print(f"  {a:5s} | {(gb-eb)*H2K:12.2f} {(gw-ew)*H2K:12.2f}")


if __name__ == "__main__":
    main()
