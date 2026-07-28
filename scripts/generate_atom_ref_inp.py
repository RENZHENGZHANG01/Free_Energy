#!/usr/bin/env python3
"""
Generate ORCA inputs for isolated-atom reference energies (atom_ref).

WHY: Delta_G = Gibbs_Eh(molecule) - SUM_i n_i * E(atom_i)  [compute_deltaG.py:144].
The atom references MUST use the SAME functional/basis as the molecules, otherwise
Delta_G is meaningless. The molecules were computed with
    ! B3LYP D3BJ def2-TZVP TightSCF Opt/Freq RIJCOSX def2/J
so the B3LYP set below is the ONE that is consistent with the existing 1277 molecules.
The wB97X-D3/def2-TZVPPD set is generated in parallel for FUTURE use only -- do NOT
mix it with B3LYP molecular energies.

Also fixes the missing Ge/Sn references (the old atom_ref.csv had only 11 elements,
while the AL in-domain vocabulary allows 13 -> compute_reference_energy() would
raise ValueError on any Ge/Sn molecule).

CRITICAL -- ground-state multiplicities (2S+1) of the free atoms. Getting these
wrong silently produces a wrong reference energy (and hence wrong Delta_G for every
molecule containing that element):
    H  1s1        2S   -> 2      Si [Ne]3s2 3p2  3P  -> 3
    C  2s2 2p2    3P   -> 3      P  [Ne]3s2 3p3  4S  -> 4
    N  2s2 2p3    4S   -> 4      S  [Ne]3s2 3p4  3P  -> 3
    O  2s2 2p4    3P   -> 3      Cl [Ne]3s2 3p5  2P  -> 2
    F  2s2 2p5    2P   -> 2      Ge/Sn  ns2 np2  3P  -> 3
                                 Br/I   ns2 np5  2P  -> 2

Usage:  python generate_atom_ref_inp.py [--outdir ../data/atom_ref]
"""
import os, argparse

# element -> ground-state multiplicity (2S+1). See docstring; do not "simplify" these.
ATOM_MULT = {
    "H": 2, "C": 3, "N": 4, "O": 3, "F": 2,
    "Si": 3, "P": 4, "S": 3, "Cl": 2,
    "Ge": 3, "Br": 2, "Sn": 3, "I": 2,
}

# keyword lines. "Freq" is included so we also get the Gibbs free energy, matching
# what extract_thermo.py pulls for molecules ("Final Gibbs free energy").
METHODS = {
    # CONSISTENT with the existing 1277 molecules -> use this one for production.
    "b3lyp": "! B3LYP D3BJ def2-TZVP TightSCF Freq RIJCOSX def2/J",
    # Newer / more accurate (range-separated hybrid + diffuse basis). FUTURE USE ONLY:
    # switching to it requires recomputing ALL molecules with the same method.
    "wb97x": "! wB97X-D3 def2-TZVPPD TightSCF Freq RIJCOSX def2/J",
}

# NOTE: run these SERIAL (no %pal). A single atom has only ~10-30 basis functions, so
# MPI communication overhead dwarfs the work: with "%pal nprocs 8" the job spun for
# 20 min on MPI "Read -1" errors and finished only 2/13, whereas serial takes seconds
# per atom. %pal affects speed only, never the energy.
TEMPLATE = """{keywords}

* xyz 0 {mult}
{sym}  0.0  0.0  0.0
*
"""


def main():
    ap = argparse.ArgumentParser()
    here = os.path.dirname(os.path.abspath(__file__))
    ap.add_argument("--outdir", default=os.path.join(here, "..", "data", "atom_ref"))
    args = ap.parse_args()

    n = 0
    for tag, kw in METHODS.items():
        d = os.path.abspath(os.path.join(args.outdir, tag))
        os.makedirs(d, exist_ok=True)
        for sym, mult in ATOM_MULT.items():
            path = os.path.join(d, f"{sym}_atom.inp")
            with open(path, "w") as fh:
                fh.write(TEMPLATE.format(keywords=kw, mult=mult, sym=sym))
            n += 1
        print(f"  {tag:6s}: {len(ATOM_MULT)} inputs -> {d}")
        print(f"          {kw}")
    print(f"\n  Generated {n} atom inputs "
          f"({len(ATOM_MULT)} elements x {len(METHODS)} methods).")
    print("  Multiplicities used:", ", ".join(f"{k}={v}" for k, v in ATOM_MULT.items()))


if __name__ == "__main__":
    main()
