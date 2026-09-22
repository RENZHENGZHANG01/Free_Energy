#!/usr/bin/env python3
"""
Step 3: production-optimised geometry -> frequency / thermochemistry input.

The level of theory comes from scripts/orca_settings.py and is the SAME functional,
basis and grid as the step 2 optimisation. This is not a style preference: a Hessian
evaluated at a different level than the geometry it sits on is not a Hessian at a
stationary point, and every such mismatch shows up as imaginary frequencies that
are pure artefact.

VeryTightSCF here (rather than TightSCF) because the Hessian differentiates the SCF
solution, so SCF noise is amplified straight into the low-frequency modes -- which
are exactly the modes the entropy term is most sensitive to.

Usage:
    python generate_freq_inp.py              # every step-2 geometry
    python generate_freq_inp.py --mol SD123  # just one (what the job script does)
"""
import os, sys, argparse

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import orca_settings as S

TEMPLATE = """{keywords}

{pal}%freq
  Temp {temp}
end

* xyz {charge} {mult}
{coords}*
"""


def read_coords(path):
    """Coordinate block of an XYZ file, with the header checked.

    ORCA fails in confusing ways on a truncated XYZ (an optimisation killed by the
    queue leaves one behind), so the atom count is verified against the body rather
    than trusted.
    """
    with open(path) as f:
        lines = f.readlines()
    if len(lines) < 3:
        raise ValueError(f"XYZ too short: {path}")
    try:
        n = int(lines[0].split()[0])
    except Exception:
        raise ValueError(f"XYZ first line is not an atom count: {path}\n{lines[0]!r}")
    coords = lines[2:2 + n]
    if len(coords) < n:
        raise ValueError(f"XYZ incomplete: {path} header says {n}, found {len(coords)}")
    bad = [ln for ln in coords if len(ln.split()) < 4]
    if bad:
        raise ValueError(f"malformed coordinate line in {path}: {bad[0]!r}")
    return "".join(coords)


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--mol", default=None)
    ap.add_argument("--opt-dir", default=os.environ.get("OPT_INP_DIR", "data/opt_inp"))
    ap.add_argument("--freq-dir", default=os.environ.get("FREQ_INP_DIR", "data/freq_inp"))
    ap.add_argument("--charge", type=int, default=0)
    ap.add_argument("--mult", type=int, default=1)
    ap.add_argument("--temp", type=float, default=S.TEMP_K)
    args = ap.parse_args()

    os.makedirs(args.freq_dir, exist_ok=True)

    suffix = f"_{S.STEP2_TAG}.xyz"
    if args.mol:
        names = [args.mol]
    else:
        names = sorted(f[:-len(suffix)] for f in os.listdir(args.opt_dir)
                       if f.endswith(suffix))

    n = 0
    for name in names:
        xyz = os.path.join(args.opt_dir, f"{name}{suffix}")
        if not os.path.exists(xyz):
            raise SystemExit(f"step 2 did not produce a geometry: {xyz}")
        out = os.path.join(args.freq_dir, f"{name}_{S.FREQ_TAG}.inp")
        with open(out, "w") as fw:
            fw.write(TEMPLATE.format(keywords=S.FREQ, pal=S.pal_block(),
                                     temp=args.temp, charge=args.charge,
                                     mult=args.mult, coords=read_coords(xyz)))
        n += 1

    print(f"Frequency inputs written: {n} -> {args.freq_dir}")
    print(f"  {S.FREQ}")


if __name__ == "__main__":
    main()
