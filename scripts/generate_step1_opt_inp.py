#!/usr/bin/env python3
"""
Step 1: RDKit geometry (data/xyz/<MOL>.xyz) -> cheap PBE pre-optimisation input.

The level of theory comes from scripts/orca_settings.py -- see that file for why.
This step is deliberately cheap and its functional is irrelevant to the final
numbers: step 2 re-optimises from scratch, so all this has to do is clean up the
force-field geometry enough that the expensive optimisation starts somewhere sane.

Usage:
    python generate_step1_opt_inp.py              # every molecule in data/xyz
    python generate_step1_opt_inp.py --mol SD123  # just one (what the job script does)
"""
import os, sys, argparse

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import orca_settings as S

TEMPLATE = """{keywords}

{pal}* xyz {charge} {mult}
{coords}*
"""


def read_coords(path):
    with open(path) as f:
        lines = f.readlines()
    if len(lines) < 3:
        raise ValueError(f"XYZ too short: {path}")
    n = int(lines[0].split()[0])
    coords = lines[2:2 + n]
    if len(coords) < n:
        raise ValueError(f"XYZ header says {n} atoms, found {len(coords)}: {path}")
    return "".join(coords)


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--mol", default=None,
                    help="single molecule id; default is every .xyz in --xyz-dir")
    ap.add_argument("--xyz-dir", default=os.environ.get("XYZ_DIR", "data/xyz"))
    ap.add_argument("--out-dir", default=os.environ.get("OPT_INP_DIR", "data/opt_inp"))
    ap.add_argument("--charge", type=int, default=0)
    ap.add_argument("--mult", type=int, default=1)
    args = ap.parse_args()

    os.makedirs(args.out_dir, exist_ok=True)

    if args.mol:
        names = [args.mol]
    else:
        names = sorted(f[:-4] for f in os.listdir(args.xyz_dir) if f.endswith(".xyz"))

    n = 0
    for name in names:
        xyz = os.path.join(args.xyz_dir, f"{name}.xyz")
        if not os.path.exists(xyz):
            raise SystemExit(f"missing geometry: {xyz}")
        out = os.path.join(args.out_dir, f"{name}_{S.STEP1_TAG}.inp")
        with open(out, "w") as fw:
            fw.write(TEMPLATE.format(keywords=S.STEP1_OPT, pal=S.pal_block(),
                                     charge=args.charge, mult=args.mult,
                                     coords=read_coords(xyz)))
        n += 1

    print(f"Step 1 inputs written: {n} -> {args.out_dir}")
    print(f"  {S.STEP1_OPT}")


if __name__ == "__main__":
    main()
