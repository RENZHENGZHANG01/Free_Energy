# import os

# xyz_dir = "data/xyz"
# out_dir = "data/opt_inp"
# os.makedirs(out_dir, exist_ok=True)

# template = """! B3LYP D3BJ def2-TZVP TightSCF TightOpt OPT

# * xyz 0 1
# {coords}
# *
# """

# for file in os.listdir(xyz_dir):
#     if file.endswith(".xyz"):
#         name = file[:-4]
#         with open(f"{xyz_dir}/{file}") as f:
#             lines = f.readlines()[2:]
#         with open(f"{out_dir}/{name}_opt.inp", "w") as fw:
#             fw.write(template.format(coords="".join(lines)))

# print("Done: ORCA OPT input created.")


# import os

# xyz_dir = "data/xyz"
# out_dir = "data/opt_inp"
# os.makedirs(out_dir, exist_ok=True)

# ###############################################
# # Step 1: PBE/SVP + RIJCOSX (fast rough opt)
# ###############################################
# template_pbe = """! PBE def2-SVP TightSCF Opt RIJCOSX def2/J

# * xyz 0 1
# {coords}
# *
# """

# ###############################################
# # Step 2: B3LYP/TZVP + RIJCOSX (accurate opt)
# ###############################################
# # NOTE: This references the PBE-optimized geometry file
# template_b3lyp = """! B3LYP D3BJ def2-TZVP TightSCF Opt RIJCOSX def2/J

# * xyzfile 0 1 {pbe_xyz}
# *
# """

# ###############################################
# # Main Loop
# ###############################################
# for file in os.listdir(xyz_dir):
#     if file.endswith(".xyz"):
#         name = file[:-4]

#         # read XYZ coordinates (skip 1st two lines)
#         with open(f"{xyz_dir}/{file}") as f:
#             lines = f.readlines()[2:]
#         coords = "".join(lines)

#         # ---- Step 1: write PBE input ----
#         inp1 = f"{out_dir}/{name}_step1_pbe_opt.inp"
#         with open(inp1, "w") as fw:
#             fw.write(template_pbe.format(coords=coords))

#         # ---- Step 2: write B3LYP input ----
#         # This step uses xyzfile from Step 1 results
#         pbe_xyz_name = f"{name}_step1_pbe_opt.xyz"
#         inp2 = f"{out_dir}/{name}_step2_b3lyp_opt.inp"
#         with open(inp2, "w") as fw:
#             fw.write(template_b3lyp.format(pbe_xyz=pbe_xyz_name))

# print("Done: Two-step ORCA OPT input files created.")


import os

xyz_dir = "data/xyz"
out_dir = "data/opt_inp"
os.makedirs(out_dir, exist_ok=True)

charge = 0
mult = 1

template_pbe = """! PBE def2-SVP TightSCF Opt RIJCOSX def2/J

* xyz {charge} {mult}
{coords}
*
"""

for file in os.listdir(xyz_dir):
    if not file.endswith(".xyz"):
        continue

    name = file[:-4]

    with open(os.path.join(xyz_dir, file)) as f:
        lines = f.readlines()

    coords = "".join(lines[2:])

    with open(os.path.join(out_dir, f"{name}_step1_pbe_opt.inp"), "w") as fw:
        fw.write(template_pbe.format(
            charge=charge,
            mult=mult,
            coords=coords
        ))

print("✔ Step 1 PBE OPT inputs generated")
