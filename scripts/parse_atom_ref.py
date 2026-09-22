#!/usr/bin/env python3
"""
Parse isolated-atom ORCA outputs -> atom_ref CSV used by compute_deltaG.py.

WHICH ENERGY GOES IN atom_ref.csv?  ->  G_Eh  (Gibbs free energy). Default --use G.

    Delta_G = G(molecule) - SUM_i n_i * G(atom_i)

Both sides must be the same kind of quantity. extract_thermo.py pulls "Final Gibbs
free energy" for the molecules, so the atoms must contribute Gibbs free energies
too. A free atom has no vibrations and no rotations, so its G is E_elec plus the
translational term and the electronic degeneracy; small, but not zero, and not
constant across elements.

    Historical note. Until 2026-09 this file defaulted to --use E and subtracted
    ELECTRONIC atom energies from MOLECULAR Gibbs energies. That is not a free
    energy of any process. It survived in practice only because the mismatch is a
    per-element constant times atom count, which the composition-level Ridge
    baseline in compute_residual_deltaG.py absorbs exactly (R^2 = 1.0000000000
    against the correction term), leaving the residual target -- the actual
    training signal -- unchanged. Delta_G itself was still uninterpretable.
    --use E is kept only for reproducing those historical numbers.

CONSISTENCY IS CHECKED, NOT ASSUMED
    Each .out echoes its own keyword line. For the "production" set this script
    compares that line against orca_settings.ATOM_REF and refuses to write
    atom_ref.csv if they differ -- which is what catches a stale directory of
    outputs left over from a previous functional. That exact failure (atoms at
    def2-TZVPPD, molecules at def2-TZVP) was live in this pipeline before.

Usage:
    python parse_atom_ref.py                                  # report only
    python parse_atom_ref.py --write-production               # writes atom_ref.csv (Gibbs)
    python parse_atom_ref.py --write-production --use E       # historical/electronic
    python parse_atom_ref.py --set legacy_b3lyp --write-production --use E
"""
import os, re, sys, argparse, csv

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import orca_settings as S

HARTREE_KCAL = 627.509474


def parse_one(path):
    txt = open(path, errors="ignore").read()

    def last(pat):
        m = re.findall(pat, txt)
        return float(m[-1]) if m else None

    # ORCA echoes the input file as "|  1> ! <keywords>". Used to verify that these
    # outputs were actually produced with the level of theory we think they were.
    kw = re.search(r"\|\s*\d+>\s*(!.*)", txt)

    return {
        "E_elec_Eh": last(r"FINAL SINGLE POINT ENERGY\s+(-?\d+\.\d+)"),
        "G_Eh":      last(r"Final Gibbs free energy\s*\.*\s*(-?\d+\.\d+)"),
        "mult":      last(r"Multiplicity\s+Mult\s*\.*\s*(\d+)"),
        "keywords":  kw.group(1).strip() if kw else None,
        "converged": "ORCA TERMINATED NORMALLY" in txt,
        "scf_fail":  ("SCF NOT CONVERGED" in txt)
                     or ("SCF ITERATIONS DID NOT CONVERGE" in txt),
    }


def norm_kw(s):
    """Compare keyword lines ignoring case and whitespace, not order."""
    return " ".join(s.lower().split()) if s else None


def read_set(d):
    rows = []
    for fn in sorted(os.listdir(d)):
        # STRICT "<El>_atom.out". ORCA emits auxiliary per-ECP-element files such as
        # "I_atom_atom53.out" / "Sn_atom_atom50.out"; those are not reference
        # energies and must never be treated as elements.
        m = re.match(r"^([A-Z][a-z]?)_atom\.out$", fn)
        if m:
            rows.append((m.group(1), parse_one(os.path.join(d, fn))))
    return rows


def report(tag, rows, expect_kw):
    print(f"\n  === {tag} ===   {len(rows)} atoms")
    print(f"  {'atom':5s} {'E_elec_Eh':>18s} {'G_Eh':>18s} {'G-E(kcal)':>10s} {'mult':>4s}  status")
    problems = []
    for sym, r in rows:
        bad = []
        if r["E_elec_Eh"] is None: bad.append("no E")
        if r["G_Eh"] is None:      bad.append("no G")
        if not r["converged"]:     bad.append("no normal termination")
        if r["scf_fail"]:          bad.append("SCF not converged")
        if expect_kw and norm_kw(r["keywords"]) != norm_kw(expect_kw):
            bad.append("keyword mismatch")
        gap = ((r["G_Eh"] - r["E_elec_Eh"]) * HARTREE_KCAL
               if (r["G_Eh"] is not None and r["E_elec_Eh"] is not None) else None)
        print(f"  {sym:5s} {str(r['E_elec_Eh']):>18s} {str(r['G_Eh']):>18s} "
              f"{(f'{gap:10.2f}' if gap is not None else ' ' * 10)} "
              f"{int(r['mult'] or 0):>4d}  {'** ' + ', '.join(bad) + ' **' if bad else 'ok'}")
        if bad:
            problems.append((sym, bad, r["keywords"]))
    return problems


def main():
    ap = argparse.ArgumentParser()
    here = os.path.dirname(os.path.abspath(__file__))
    ap.add_argument("--dir", default=os.path.join(here, "..", "data", "atom_ref"))
    ap.add_argument("--set", default="production",
                    help="subdirectory of --dir to use (default: production)")
    ap.add_argument("--write-production", action="store_true",
                    help="overwrite scripts/atom_ref.csv from the chosen set")
    ap.add_argument("--use", choices=["G", "E"], default="G",
                    help="energy column for atom_ref.csv. Default G (Gibbs), which "
                         "is what matches the molecular side; E is historical only.")
    ap.add_argument("--force", action="store_true",
                    help="write even if the keyword-consistency check fails")
    args = ap.parse_args()

    base = os.path.abspath(args.dir)
    d = os.path.join(base, args.set)
    if not os.path.isdir(d):
        avail = [x for x in sorted(os.listdir(base))
                 if os.path.isdir(os.path.join(base, x))] if os.path.isdir(base) else []
        raise SystemExit(f"no such set: {d}\n  available: {avail}")

    rows = read_set(d)
    if not rows:
        raise SystemExit(f"no <El>_atom.out files in {d} -- has the job finished?")

    # Only the production set is pinned to orca_settings; a legacy set is by
    # definition at a different level of theory.
    expect_kw = S.ATOM_REF if args.set == "production" else None
    if expect_kw:
        print(f"  expecting: {expect_kw}")

    problems = report(args.set, rows, expect_kw)

    out_csv = os.path.join(base, f"atom_ref_{args.set}.csv")
    with open(out_csv, "w", newline="") as fh:
        w = csv.writer(fh)
        w.writerow(["atom", "E_elec_Eh", "G_Eh", "mult", "ok", "keywords"])
        for sym, r in rows:
            ok = (r["converged"] and not r["scf_fail"]
                  and r["E_elec_Eh"] is not None and r["G_Eh"] is not None)
            w.writerow([sym, r["E_elec_Eh"], r["G_Eh"], int(r["mult"] or 0),
                        "yes" if ok else "NO", r["keywords"] or ""])
    print(f"\n  full table -> {out_csv}")

    if not args.write_production:
        return

    if problems and not args.force:
        print("\n  ABORT: not writing atom_ref.csv. Problems:")
        for sym, bad, kw in problems:
            print(f"    {sym:3s} {', '.join(bad)}")
            if "keyword mismatch" in bad:
                print(f"        found:    {kw}")
                print(f"        expected: {expect_kw}")
        print("  Fix the runs (or pass --force if you are certain).")
        raise SystemExit(1)

    key = "G_Eh" if args.use == "G" else "E_elec_Eh"
    vals = {s: r[key] for s, r in rows}
    if any(v is None for v in vals.values()):
        raise SystemExit("\n  ABORT: some energies missing; not writing atom_ref.csv")

    prod = os.path.join(here, "atom_ref.csv")
    old = {}
    if os.path.exists(prod):
        with open(prod) as fh:
            for row in csv.DictReader(fh):
                old[row["atom"]] = float(row["energy"])
        # Back up first: this file drives every molecule's Delta_G.
        import shutil, time as _t
        bak = prod + ".bak_" + _t.strftime("%Y%m%d_%H%M%S")
        shutil.copy2(prod, bak)
        print(f"\n  backed up old file -> {os.path.basename(bak)}")

    with open(prod, "w", newline="") as fh:
        w = csv.writer(fh)
        w.writerow(["atom", "energy"])
        for s, v in sorted(vals.items()):
            w.writerow([s, f"{v:.12f}"])

    print(f"\n  WROTE {prod}")
    print(f"  column 'energy' = {args.set} {key}")
    print(f"  level of theory = {rows[0][1]['keywords']}")
    print(f"\n  {'atom':5s} {'old':>18s} {'new':>18s} {'diff(Eh)':>12s} {'diff(kcal/mol)':>15s}")
    for s, v in sorted(vals.items()):
        o = old.get(s)
        if o is None:
            print(f"  {s:5s} {'(new)':>18s} {v:18.9f}")
        else:
            dd = v - o
            print(f"  {s:5s} {o:18.9f} {v:18.9f} {dd:12.6f} {dd * HARTREE_KCAL:15.2f}")
    moved = [s for s, v in vals.items()
             if s in old and abs(v - old[s]) * HARTREE_KCAL > 0.01]
    print(f"\n  Elements that moved >0.01 kcal/mol: {len(moved)}")
    if moved:
        print("  => Every molecule's Delta_G changes. The molecular dataset must be")
        print("     recomputed at the same level of theory before this file is used:")
        print("       scripts/extract_thermo.py -> compute_deltaG.py -> compute_residual_deltaG.py")


if __name__ == "__main__":
    main()
