#!/usr/bin/env python3
"""
Step 2: PBE-optimised geometry -> production wB97X-D3 optimisation input.

The level of theory comes from scripts/orca_settings.py -- see that file for the
benchmark that chose it. In short: wB97X-D3/def2-TZVP with RIJCOSX and DefGrid3,
TightOpt because a loose minimum turns into spurious imaginary frequencies in the
Hessian computed at the next step.

Usage:
    python generate_step2_opt_inp_from_xyz.py              # every optimised geometry
    python generate_step2_opt_inp_from_xyz.py --mol SD123  # just one

NOTE the --mol flag. submit_free_energy.csh has always passed it, but earlier
versions of this script did not parse any arguments and silently rebuilt the input
for EVERY molecule in the directory -- so with a few thousand jobs in flight, each
one rewrote every other one's input file while they were being read.
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
    bad = [ln for ln in coords if len(ln.split()) < 4]
    if bad:
        raise ValueError(f"malformed coordinate line in {path}: {bad[0]!r}")
    return "".join(coords)


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--mol", default=None)
    ap.add_argument("--opt-dir", default=os.environ.get("OPT_INP_DIR", "data/opt_inp"))
    ap.add_argument("--charge", type=int, default=0)
    ap.add_argument("--mult", type=int, default=1)
    args = ap.parse_args()

    suffix = f"_{S.STEP1_TAG}.xyz"
    if args.mol:
        names = [args.mol]
    else:
        names = sorted(f[:-len(suffix)] for f in os.listdir(args.opt_dir)
                       if f.endswith(suffix))

    n = 0
    for name in names:
        xyz = os.path.join(args.opt_dir, f"{name}{suffix}")
        if not os.path.exists(xyz):
            raise SystemExit(f"step 1 did not produce a geometry: {xyz}")
        out = os.path.join(args.opt_dir, f"{name}_{S.STEP2_TAG}.inp")
        with open(out, "w") as fw:
            fw.write(TEMPLATE.format(keywords=S.STEP2_OPT, pal=S.pal_block(),
                                     charge=args.charge, mult=args.mult,
                                     coords=read_coords(xyz)))
        n += 1

    print(f"Step 2 inputs written: {n} -> {args.opt_dir}")
    print(f"  {S.STEP2_OPT}")


if __name__ == "__main__":
    main()
