import os, re, pandas as pd

freq_dir = "data/freq_out"
records = []

for file in os.listdir(freq_dir):
    if file.endswith("_freq.out"):
        name = file[:-9]
        G, thermal_corr = None, None

        with open(os.path.join(freq_dir,file)) as f:
            for line in f:
                if "Final Gibbs free energy" in line:
                    # Example: Final Gibbs free energy            ...   -236.99213046 Eh
                    G = float(re.search(r"(-?\d+\.\d+)", line).group(1))

                if "G-E(el)" in line:
                    # Example: G-E(el) ... 0.15632109 Eh
                    thermal_corr = float(re.search(r"(-?\d+\.\d+)", line).group(1))

        records.append([name, G, thermal_corr])

df = pd.DataFrame(records, columns=["mol", "Gibbs_Eh", "G_minus_Eel"])
df.to_csv("data/deltaG_raw.csv", index=False)

print("Saved: data/deltaG_raw.csv")
print(df.head())
