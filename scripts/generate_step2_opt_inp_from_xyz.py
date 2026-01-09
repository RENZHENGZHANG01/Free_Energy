import os

opt_dir = "data/opt_inp"

charge = 0
mult = 1

template_b3lyp = """! B3LYP D3BJ def2-TZVP TightSCF Opt RIJCOSX def2/J

* xyz {charge} {mult}
{coords}
*
"""

for file in os.listdir(opt_dir):
    if not file.endswith("_step1_pbe_opt.xyz"):
        continue

    name = file.replace("_step1_pbe_opt.xyz", "")

    with open(os.path.join(opt_dir, file)) as f:
        lines = f.readlines()

    # skip atom count + comment
    coords = "".join(lines[2:])

    with open(os.path.join(opt_dir, f"{name}_step2_b3lyp_opt.inp"), "w") as fw:
        fw.write(template_b3lyp.format(
            charge=charge,
            mult=mult,
            coords=coords
        ))

    print(f"✔ Generated step2 B3LYP OPT input for {name}")
