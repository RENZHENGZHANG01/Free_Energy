# import os, re

# xyz_dir = "data/opt_inp"     # READ OPTIMIZED XYZ HERE!!
# freq_inp_dir = "data/freq_inp"
# os.makedirs(freq_inp_dir, exist_ok=True)

# # ORCA Frequency template
# template = """! B3LYP D3BJ def2-TZVP TightSCF Freq

# %freq
#   Temp 298.15
# end

# * xyz 0 1
# {coords}
# *
# """

# # Loop through all optimized XYZ files
# for file in os.listdir(xyz_dir):
#     if file.endswith("_opt.xyz"):
#         name = file[:-8]  # e.g., M2_opt.xyz → M2
#         with open(os.path.join(xyz_dir, file), "r") as f:
#             lines = f.readlines()
#             coords = "".join(lines[2:])  # skip atom_count + title lines

#         # Write frequency input
#         out_path = os.path.join(freq_inp_dir, f"{name}_freq.inp")
#         with open(out_path, "w") as fw:
#             fw.write(template.format(coords=coords))

# print(f"✔ Frequency inputs saved in {freq_inp_dir}")


#!/usr/bin/env python3
import os

xyz_dir = "data/opt_inp"     
freq_inp_dir = "data/freq_inp"
os.makedirs(freq_inp_dir, exist_ok=True)

charge = 0
mult = 1
temp = 298.15

template = """! B3LYP D3BJ def2-TZVP TightSCF Freq RIJCOSX def2/J

%freq
  Temp {temp}
end

* xyz {charge} {mult}
{coords}
*
"""

def read_xyz_coords(xyz_path: str) -> str:
    """Read XYZ file and return coordinate lines (skip first 2 lines).
    Also validates XYZ header to avoid ORCA xyzfile pitfalls.
    """
    with open(xyz_path, "r") as f:
        lines = f.readlines()

    if len(lines) < 3:
        raise ValueError(f"XYZ too short: {xyz_path}")

    # Validate first line is atom count integer
    try:
        n_atoms = int(lines[0].strip())
    except Exception:
        raise ValueError(f"XYZ first line is not an integer atom count: {xyz_path}\nFirst line: {lines[0]!r}")

    coord_lines = lines[2:]
    if len(coord_lines) < n_atoms:
        raise ValueError(
            f"XYZ seems incomplete: {xyz_path}\n"
            f"Header says {n_atoms} atoms but only {len(coord_lines)} coordinate lines found."
        )

    # Basic sanity: each coord line should have at least 4 tokens
    bad = [ln for ln in coord_lines[:n_atoms] if len(ln.split()) < 4]
    if bad:
        raise ValueError(f"XYZ coordinate line format looks wrong in {xyz_path}\nExample bad line: {bad[0]!r}")

    return "".join(coord_lines[:n_atoms])

count = 0
for file in sorted(os.listdir(xyz_dir)):
    if not file.endswith("_step2_b3lyp_opt.xyz"):
        continue

    name = file.replace("_step2_b3lyp_opt.xyz", "")
    xyz_path = os.path.join(xyz_dir, file)

    coords = read_xyz_coords(xyz_path)

    out_path = os.path.join(freq_inp_dir, f"{name}_freq.inp")
    with open(out_path, "w") as fw:
        fw.write(template.format(temp=temp, charge=charge, mult=mult, coords=coords))

    count += 1

print(f"✔ Frequency inputs generated: {count} (saved in {freq_inp_dir})")
