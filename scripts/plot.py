import pandas as pd
import matplotlib.pyplot as plt

# Load final ΔG dataset
df = pd.read_csv("final_data_with_deltaG.csv")

# X-axis
x_col = "total_atoms"

# Y columns: energy-related features
y_cols = [
    # "Delta_G",
    "DeltaG_per_heavy_atom",
    "DeltaG_per_atom",
    "DeltaG_per_backbone_bond",
    "DeltaG_per_bond"
]

# Colors and markers for scatter points
colors = ['red', 'green', 'purple', 'orange']
markers = ['o', 's', '^', 'D']

plt.figure(figsize=(10, 6))

# Plot scatter for each energy column
for i, col in enumerate(y_cols):
    plt.scatter(
        df[x_col],
        df[col],
        label=col.replace("_", " ").title(),
        color=colors[i],
        marker=markers[i],
        s=60,              # Adjust marker size
        alpha=0.8,         # Slight transparency
        edgecolor='black', # Make points more visible
        linewidth=0.5
    )

# Optional: label points with monomer name (if exists)
if "Monomer" in df.columns:
    for i in df.index:
        plt.text(
            df[x_col][i] + 0.5,
            df["DeltaG_per_atom"][i],
            str(df["Monomer"][i]),
            fontsize=8
        )

# Axes labels
plt.xlabel("Total Number of Atoms", fontsize=12)
plt.ylabel("Energy (a.u.)", fontsize=12)
plt.title("Scatter Plot: Total Atoms vs Energy Components", fontsize=14)

plt.legend(title="Energy Type", fontsize=10)
plt.grid(True, linestyle='--', alpha=0.5)

plt.savefig("deltaG_scatter_plot.png", dpi=300, bbox_inches='tight')
print("Saved scatter plot as deltaG_scatter_plot.png")
