#!/usr/bin/env python3
"""
Generate ORCA inputs for the isolated-atom reference energies.

WHY THIS EXISTS
    Delta_G = G(molecule) - SUM_i n_i * G(atom_i)      [compute_deltaG.py]

The atom references MUST use the same functional, basis, grid and integral
approximation as the molecules, or the subtraction cancels nothing. The keyword
line therefore is not written here at all -- it is imported from orca_settings.py,
the same object the molecular frequency inputs use. They cannot drift apart.

    production   = orca_settings.ATOM_REF  (== orca_settings.FREQ)
    legacy_b3lyp = the old B3LYP/def2-TZVP set, emitted only with --legacy-b3lyp,
                   and useful for one thing only: reproducing the historical
                   Delta_G values of the pre-2026-09 dataset. Never mix it with
                   wB97X-D3 molecular energies.

    (The previous version of this script generated a "wb97x" set at def2-TZVPPD
    while the molecules ran at def2-TZVP. That mismatch would have silently
    corrupted every Delta_G, which is exactly why the string now comes from
    one shared module.)

WHY GIBBS, NOT ELECTRONIC ENERGY
    The molecule side of the subtraction is a Gibbs free energy, so the atom side
    has to be one too. "Freq" is in the keyword line so ORCA prints
    "Final Gibbs free energy" for the atoms as well. A free atom has no vibrations
    and no rotations, so its G is just E_elec plus the translational term and the
    electronic degeneracy -- cheap, but it has to be there.
    (The old pipeline subtracted ELECTRONIC atom energies from MOLECULAR Gibbs
    energies. The mismatch is a per-element constant times atom count, which the
    composition-level Ridge baseline absorbed exactly, so the residual target
    survived -- but Delta_G itself was not a free energy of anything.)

CRITICAL -- ground-state multiplicities (2S+1) of the free atoms. Getting one wrong
silently produces a wrong reference, and hence a wrong Delta_G for every molecule
containing that element:
    H  1s1        2S   -> 2      Si [Ne]3s2 3p2  3P  -> 3
    C  2s2 2p2    3P   -> 3      P  [Ne]3s2 3p3  4S  -> 4
    N  2s2 2p3    4S   -> 4      S  [Ne]3s2 3p4  3P  -> 3
    O  2s2 2p4    3P   -> 3      Cl [Ne]3s2 3p5  2P  -> 2
    F  2s2 2p5    2P   -> 2      Ge/Sn  ns2 np2  3P  -> 3
                                 Br/I   ns2 np5  2P  -> 2

Usage:  python generate_atom_ref_inp.py [--outdir ../data/atom_ref] [--legacy-b3lyp]
"""
import os, sys, argparse

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import orca_settings as S

# element -> ground-state multiplicity (2S+1). See docstring; do not "simplify" these.
ATOM_MULT = {
    "H": 2, "C": 3, "N": 4, "O": 3, "F": 2,
    "Si": 3, "P": 4, "S": 3, "Cl": 2,
    "Ge": 3, "Br": 2, "Sn": 3, "I": 2,
}

LEGACY_B3LYP = "! B3LYP D3BJ def2-TZVP TightSCF Freq RIJCOSX def2/J"

# NOTE: no %pal. A single atom has ~10-30 basis functions, so MPI communication
# dwarfs the work: an earlier "%pal nprocs 8" attempt spun 20 minutes on MPI
# "Read -1" errors and finished 2 of 13, where serial takes seconds per atom.
# %pal affects speed only, never the energy.
TEMPLATE = """{keywords}

* xyz 0 {mult}
{sym}  0.0  0.0  0.0
*
"""


def write_set(outdir, tag, keywords):
    d = os.path.abspath(os.path.join(outdir, tag))
    os.makedirs(d, exist_ok=True)
    for sym, mult in ATOM_MULT.items():
        with open(os.path.join(d, f"{sym}_atom.inp"), "w") as fh:
            fh.write(TEMPLATE.format(keywords=keywords, mult=mult, sym=sym))
    print(f"  {tag:13s}: {len(ATOM_MULT)} inputs -> {d}")
    print(f"  {'':13s}  {keywords}")
    return len(ATOM_MULT)


def main():
    ap = argparse.ArgumentParser()
    here = os.path.dirname(os.path.abspath(__file__))
    ap.add_argument("--outdir", default=os.path.join(here, "..", "data", "atom_ref"))
    ap.add_argument("--legacy-b3lyp", action="store_true",
                    help="also emit the old B3LYP set (only for reproducing "
                         "historical Delta_G; do not mix with production data)")
    args = ap.parse_args()

    n = write_set(args.outdir, "production", S.ATOM_REF)
    if args.legacy_b3lyp:
        n += write_set(args.outdir, "legacy_b3lyp", LEGACY_B3LYP)

    print(f"\n  Generated {n} atom inputs.")
    print("  Multiplicities:", ", ".join(f"{k}={v}" for k, v in ATOM_MULT.items()))
    print("\n  Next:  qsub submit_atom_ref.csh")
    print("         python scripts/parse_atom_ref.py --write-production")


if __name__ == "__main__":
    main()
