import os

xyz_dir = "data/xyz"
out_dir = "data/opt_inp"
os.makedirs(out_dir, exist_ok=True)

charge = 0
mult = 1

template_pbe = """! PBE def2-SVP TightSCF Opt RIJCOSX def2/J

%pal nprocs 8 end

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
