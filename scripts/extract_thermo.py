#!/usr/bin/env python3
"""
Parse ORCA frequency outputs -> thermochemistry table (data/deltaG_raw.csv).

This used to pull two numbers and trust them. It now records whether the
calculation is trustworthy at all, because the two ways a frequency job goes wrong
are both silent in the old output:

  1. THE JOB DID NOT FINISH. A run killed by the queue wall-clock, or one whose SCF
     never converged, still leaves a .out file. If it got far enough to print a
     "Final Gibbs free energy" from an earlier step, that number was read and used.
     Those rows are now written with Gibbs_Eh = NaN so compute_deltaG.py drops them,
     and listed by name at the end.

  2. IMAGINARY FREQUENCIES. The harmonic thermochemistry ORCA prints assumes every
     mode is a real vibration. An imaginary mode means the geometry is not a
     minimum, and the entropy term is computed from a frequency that does not
     physically exist. Two very different causes get the same symbol:

       |v| < IMAG_NOISE_CM (default 20)   numerical noise -- a nearly-flat torsion
           on a not-quite-converged geometry, or COSX grid noise. Measured on the
           earlier B3LYP dataset: 89.7% of all imaginary modes were in this range.
           Its effect on G is a fraction of a kcal/mol. Usually acceptable.

       |v| >= IMAG_NOISE_CM                a genuine saddle point -- the optimiser
           stopped on a transition state, not a minimum. 0.5% of that dataset.
           These need re-optimisation from a geometry displaced along the imaginary
           mode; the Gibbs value as printed is not the molecule's.

     Both are reported. Neither is discarded by default (the classification depends
     on a threshold, and that is the user's call) -- use --drop-saddle to blank the
     genuine ones.

  It also records the level of theory echoed in each .out. Mixing functionals or
  basis sets inside one dataset makes every Delta_G meaningless, and this is the
  cheapest place to notice it: --expect-level compares against orca_settings.FREQ.

Usage:
    python extract_thermo.py
    python extract_thermo.py --expect-level          # fail on a mixed-method dataset
    python extract_thermo.py --drop-saddle           # also blank real saddle points
    python extract_thermo.py --imag-noise-cm 50      # looser noise threshold
"""
import os, re, sys, argparse
import pandas as pd

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import orca_settings as S

HARTREE_KCAL = 627.509474

# "     6:      -8.38 cm**-1  ***imaginary mode***"
FREQ_LINE = re.compile(r"^\s*(\d+):\s*(-?\d+\.\d+)\s*cm\*\*-1(.*)$")


def _last_float(txt, pat):
    m = re.findall(pat, txt)
    return float(m[-1]) if m else None


def parse_freq_out(path):
    txt = open(path, errors="ignore").read()

    rec = {
        "terminated":    "ORCA TERMINATED NORMALLY" in txt,
        "scf_fail":      ("SCF NOT CONVERGED" in txt)
                         or ("SCF ITERATIONS DID NOT CONVERGE" in txt),
        "E_elec_Eh":     _last_float(txt, r"FINAL SINGLE POINT ENERGY\s+(-?\d+\.\d+)"),
        "ZPE_Eh":        _last_float(txt, r"Zero point energy\s*\.*\s*(-?\d+\.\d+)\s*Eh"),
        "entropy_Eh":    _last_float(txt, r"Final entropy term\s*\.*\s*(-?\d+\.\d+)\s*Eh"),
        "thermal_Eh":    _last_float(txt, r"Total thermal correction\s*\.*\s*(-?\d+\.\d+)\s*Eh"),
        "Gibbs_Eh":      _last_float(txt, r"Final Gibbs free energy\s*\.*\s*(-?\d+\.\d+)"),
        "G_minus_Eel":   _last_float(txt, r"G-E\(el\)\s*\.*\s*(-?\d+\.\d+)\s*Eh"),
        "n_atoms":       _last_float(txt, r"Number of atoms\s*\.*\s*(\d+)"),
    }

    m = re.search(r"\|\s*\d+>\s*(!.*)", txt)
    rec["level"] = m.group(1).strip() if m else ""

    # Imaginary frequencies. Take the LAST "VIBRATIONAL FREQUENCIES" block: a
    # restarted or multi-step job prints more than one, and only the final geometry
    # counts. The first six entries are the projected-out translations and rotations
    # and are printed as exactly 0.00; they are not modes and are skipped.
    #
    # Two counts are kept because they genuinely differ. ORCA only appends
    # "***imaginary mode***" above an internal threshold of about 1 cm^-1, so a run
    # can print "-0.45 cm**-1" with no marker at all. Counting sign is the complete
    # answer; counting markers is what a grep of the output would give. Reporting
    # only one of them makes the table disagree with the .out file it came from.
    starts = [m.start() for m in re.finditer(r"^VIBRATIONAL FREQUENCIES", txt, re.M)]
    neg, flagged = [], []
    if starts:
        block = txt[starts[-1]:].split("NORMAL MODES")[0]
        for line in block.splitlines():
            m = FREQ_LINE.match(line)
            if not m:
                continue
            v = float(m.group(2))
            marked = "imaginary" in m.group(3).lower()
            if v < 0.0:
                neg.append(v)
            if marked:
                flagged.append(v)
    rec["neg_freqs"] = neg
    rec["n_neg"] = len(neg)
    rec["n_imag_flagged"] = len(flagged)
    rec["min_freq_cm"] = min(neg) if neg else None
    rec["has_freq_block"] = bool(starts)
    return rec


def main():
    ap = argparse.ArgumentParser()
    here = os.path.dirname(os.path.abspath(__file__))
    root = os.path.dirname(here)
    ap.add_argument("--freq-dir", default=os.environ.get(
        "FREQ_OUT_DIR", os.path.join(root, "data", "freq_out")))
    ap.add_argument("--out", default=os.environ.get(
        "THERMO_CSV", os.path.join(root, "data", "deltaG_raw.csv")))
    ap.add_argument("--imag-noise-cm", type=float, default=20.0,
                    help="|v| below this is treated as numerical noise (default 20)")
    ap.add_argument("--drop-saddle", action="store_true",
                    help="also blank Gibbs_Eh for genuine saddle points")
    ap.add_argument("--expect-level", action="store_true",
                    help="fail if any output's keyword line differs from "
                         "orca_settings.FREQ")
    args = ap.parse_args()

    files = sorted(f for f in os.listdir(args.freq_dir) if f.endswith("_freq.out"))
    if not files:
        raise SystemExit(f"no *_freq.out in {args.freq_dir}")

    rows = []
    for fn in files:
        name = fn[:-len("_freq.out")]
        r = parse_freq_out(os.path.join(args.freq_dir, fn))

        # Classification is by magnitude, not by ORCA's marker: a -0.45 cm^-1 mode
        # is noise whether or not ORCA chose to label it.
        n_noise = sum(1 for v in r["neg_freqs"] if abs(v) < args.imag_noise_cm)
        n_saddle = r["n_neg"] - n_noise

        # Unusable: the job did not finish, the SCF failed, or no frequencies were
        # ever printed. The Gibbs value, if any, belongs to an intermediate step.
        unusable = (not r["terminated"]) or r["scf_fail"] or (not r["has_freq_block"])
        if args.drop_saddle and n_saddle > 0:
            unusable = True

        rows.append({
            "mol": name,
            "Gibbs_Eh": None if unusable else r["Gibbs_Eh"],
            "G_minus_Eel": None if unusable else r["G_minus_Eel"],
            "E_elec_Eh": r["E_elec_Eh"],
            "ZPE_Eh": r["ZPE_Eh"],
            "entropy_term_Eh": r["entropy_Eh"],
            "thermal_corr_Eh": r["thermal_Eh"],
            "n_atoms": int(r["n_atoms"]) if r["n_atoms"] else None,
            "terminated_normally": r["terminated"],
            "scf_converged": not r["scf_fail"],
            "n_neg_freq": r["n_neg"],
            "n_imag_orca_flagged": r["n_imag_flagged"],
            "n_imag_noise": n_noise,
            "n_imag_saddle": n_saddle,
            "min_freq_cm": r["min_freq_cm"],
            "level": r["level"],
            "usable": not unusable,
        })

    df = pd.DataFrame(rows)
    os.makedirs(os.path.dirname(os.path.abspath(args.out)), exist_ok=True)
    df.to_csv(args.out, index=False)

    n = len(df)
    print(f"Parsed {n} frequency outputs -> {args.out}\n")
    print(f"  usable                     {int(df['usable'].sum()):5d} / {n}")
    print(f"  did not terminate normally {int((~df['terminated_normally']).sum()):5d}")
    print(f"  SCF not converged          {int((~df['scf_converged']).sum()):5d}")
    print(f"  clean (no negative modes)  {int((df['n_neg_freq'] == 0).sum()):5d}")
    print(f"  imaginary, |v| < {args.imag_noise_cm:g} cm-1  "
          f"{int((df['n_imag_noise'] > 0).sum()):5d}   (numerical noise)")
    print(f"  imaginary, |v| >= {args.imag_noise_cm:g} cm-1 "
          f"{int((df['n_imag_saddle'] > 0).sum()):5d}   (genuine saddle points)")
    print(f"  (ORCA's own ***imaginary mode*** marker fires on "
          f"{int((df['n_imag_orca_flagged'] > 0).sum())} of these; it has an internal "
          f"~1 cm-1 threshold)")

    bad = df[~df["terminated_normally"] | ~df["scf_converged"]]
    if len(bad):
        print(f"\n  UNUSABLE ({len(bad)}), Gibbs_Eh blanked:")
        for m in bad["mol"].head(30):
            print(f"    {m}")
        if len(bad) > 30:
            print(f"    ... and {len(bad) - 30} more")

    sad = df[df["n_imag_saddle"] > 0].sort_values("min_freq_cm")
    if len(sad):
        print(f"\n  SADDLE POINTS ({len(sad)}) -- re-optimise from a geometry "
              f"displaced along the imaginary mode:")
        for _, r in sad.head(30).iterrows():
            print(f"    {r['mol']:20s} {r['n_imag_saddle']} mode(s), "
                  f"lowest {r['min_freq_cm']:.2f} cm-1")
        if len(sad) > 30:
            print(f"    ... and {len(sad) - 30} more")

    levels = sorted(set(df["level"]) - {""})
    print(f"\n  level(s) of theory found: {len(levels)}")
    for lv in levels:
        print(f"    {int((df['level'] == lv).sum()):5d}  {lv}")
    if len(levels) > 1:
        print("  WARNING: more than one level of theory in this dataset. Delta_G is")
        print("  a difference of absolute energies -- mixing methods makes it meaningless.")
    if args.expect_level:
        wrong = [lv for lv in levels
                 if " ".join(lv.lower().split()) != " ".join(S.FREQ.lower().split())]
        if wrong:
            print(f"\n  ABORT (--expect-level): expected\n    {S.FREQ}")
            raise SystemExit(1)


if __name__ == "__main__":
    main()
