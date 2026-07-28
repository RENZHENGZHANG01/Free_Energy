#!/usr/bin/env python3
"""
Parse isolated-atom ORCA outputs -> atom_ref CSVs.

Emits BOTH energies per atom so the choice is explicit:
  E_elec_Eh : "FINAL SINGLE POINT ENERGY"  (electronic energy only)
  G_Eh      : "Final Gibbs free energy"    (what extract_thermo.py pulls for MOLECULES)

Which one belongs in atom_ref.csv?  --> E_elec_Eh  (verified 2026-07-24)
  The EXISTING atom_ref.csv holds ELECTRONIC energies: a fresh
  B3LYP-D3BJ/def2-TZVP/mult=3 run of the C atom reproduces
  FINAL SINGLE POINT ENERGY = -37.838153520618 vs the stored -37.83815352 (exact match),
  while that atom's Gibbs value is -37.85270000. So the 11 original elements were
  computed correctly with this very setup, and Ge/Sn must be added the SAME way
  (electronic energy) or atom_ref.csv becomes internally inconsistent and every
  molecule's Delta_G shifts.

  Note Delta_G = Gibbs_Eh(molecule) - SUM n_i * E_elec(atom_i) is therefore not a
  strict formation FREE energy (molecule side has thermal terms, atom side does not).
  That mismatch is a per-element CONSTANT x atom count, which the composition-level
  Ridge baseline (N_C, N_H, ... features) absorbs exactly -- so the RESIDUAL target
  (the actual training signal) is unaffected. Switching to G_Eh would change every
  historical Delta_G, so DON'T, unless you recompute the whole dataset.

Also flags SCF convergence and multiplicity so a silently-wrong reference is caught.
Usage:  python parse_atom_ref.py [--dir ../data/atom_ref] [--write-production]
"""
import os, re, argparse, csv

def parse_one(path):
    txt = open(path, errors="ignore").read()
    def last(pat):
        m = re.findall(pat, txt)
        return float(m[-1]) if m else None
    rec = {
        "E_elec_Eh": last(r"FINAL SINGLE POINT ENERGY\s+(-?\d+\.\d+)"),
        "G_Eh":      last(r"Final Gibbs free energy\s*\.*\s*(-?\d+\.\d+)"),
        "mult":      last(r"Multiplicity\s+Mult\s*\.*\s*(\d+)"),
        "converged": ("SUCCESS" in txt) or ("ORCA TERMINATED NORMALLY" in txt),
        "scf_fail":  ("SCF NOT CONVERGED" in txt) or ("SCF ITERATIONS DID NOT CONVERGE" in txt),
    }
    return rec


def main():
    ap = argparse.ArgumentParser()
    here = os.path.dirname(os.path.abspath(__file__))
    ap.add_argument("--dir", default=os.path.join(here, "..", "data", "atom_ref"))
    ap.add_argument("--write-production", action="store_true",
                    help="overwrite scripts/atom_ref.csv from the B3LYP set")
    ap.add_argument("--use", choices=["G", "E"], default="E",
                    help="energy column for the production atom_ref.csv. Default E "
                         "(electronic) -- matches the existing 11 entries; see module docstring.")
    args = ap.parse_args()

    base = os.path.abspath(args.dir)
    summary = {}
    for tag in ("b3lyp", "wb97x"):
        d = os.path.join(base, tag)
        if not os.path.isdir(d):
            continue
        rows = []
        for fn in sorted(os.listdir(d)):
            # STRICT match "<El>_atom.out" only. ORCA emits extra per-ECP-element files
            # like "I_atom_atom53.out" / "Sn_atom_atom50.out" (auxiliary ECP atom runs);
            # those are NOT reference energies and must not be treated as elements.
            m = re.match(r"^([A-Z][a-z]?)_atom\.out$", fn)
            if not m:
                continue
            sym = m.group(1)
            r = parse_one(os.path.join(d, fn))
            rows.append((sym, r))
        if not rows:
            print(f"  {tag}: no .out files yet"); continue
        out_csv = os.path.join(base, f"atom_ref_{tag}.csv")
        with open(out_csv, "w", newline="") as fh:
            w = csv.writer(fh); w.writerow(["atom", "E_elec_Eh", "G_Eh", "mult", "ok"])
            for sym, r in rows:
                w.writerow([sym, r["E_elec_Eh"], r["G_Eh"], int(r["mult"] or 0),
                            "yes" if (r["converged"] and not r["scf_fail"]) else "NO"])
        summary[tag] = rows
        print(f"\n  === {tag} ===  -> {out_csv}")
        print(f"  {'atom':5s} {'E_elec_Eh':>16s} {'G_Eh':>16s} {'mult':>4s}  status")
        for sym, r in rows:
            bad = (not r["converged"]) or r["scf_fail"] or r["E_elec_Eh"] is None
            print(f"  {sym:5s} {str(r['E_elec_Eh']):>16s} {str(r['G_Eh']):>16s} "
                  f"{int(r['mult'] or 0):>4d}  {'** CHECK **' if bad else 'ok'}")
        missing = [s for s, r in rows if r["E_elec_Eh"] is None]
        if missing: print(f"  !! no energy parsed for: {missing}")

    # optional: write the production file used by compute_deltaG.py
    if args.write_production and "b3lyp" in summary:
        rows = summary["b3lyp"]
        key = "G_Eh" if args.use == "G" else "E_elec_Eh"
        vals = {s: r[key] for s, r in rows}
        if any(v is None for v in vals.values()):
            print("\n  ABORT: some B3LYP energies missing; not writing production atom_ref.csv")
            return
        prod = os.path.join(here, "atom_ref.csv")
        old = {}
        if os.path.exists(prod):
            with open(prod) as fh:
                for row in csv.DictReader(fh):
                    old[row["atom"]] = float(row["energy"])
        # back up before overwriting (the file drives every molecule's Delta_G)
        if os.path.exists(prod):
            import shutil, time as _t
            bak = prod + ".bak_" + _t.strftime("%Y%m%d_%H%M%S")
            shutil.copy2(prod, bak)
            print(f"\n  backed up old file -> {os.path.basename(bak)}")
        with open(prod, "w", newline="") as fh:
            w = csv.writer(fh); w.writerow(["atom", "energy"])
            for s, v in vals.items(): w.writerow([s, f"{v:.12f}"])
        print(f"\n  WROTE production {prod}  (column 'energy' = B3LYP {key})")
        print(f"  {'atom':5s} {'old':>16s} {'new':>16s} {'diff(Eh)':>12s} {'diff(kcal/mol)':>15s}")
        for s, v in vals.items():
            o = old.get(s)
            if o is None:
                print(f"  {s:5s} {'(new)':>16s} {v:16.9f} {'':>12s} {'':>15s}")
            else:
                d = v - o
                print(f"  {s:5s} {o:16.9f} {v:16.9f} {d:12.6f} {d*627.509474:15.2f}")
        n_chg = sum(1 for s, v in vals.items()
                    if s in old and abs(v - old[s]) * 627.509474 > 0.01)
        print(f"\n  Elements whose value moved >0.01 kcal/mol: {n_chg}")
        if n_chg == 0:
            print("  => the 11 original entries reproduced EXACTLY; only Ge/Sn are new.")
            print("     Existing Delta_G values stay valid; just re-run compute_deltaG.py")
            print("     (+ compute_residual_deltaG.py) so Ge/Sn molecules stop erroring.")
        else:
            print("  => WARNING: existing entries CHANGED. Every molecule's Delta_G shifts;")
            print("     you must re-run compute_deltaG.py + compute_residual_deltaG.py")
            print("     and retrain. Investigate before accepting.")


if __name__ == "__main__":
    main()
