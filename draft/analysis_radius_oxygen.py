import re
import pandas as pd
import matplotlib.pyplot as plt

# ---- settings ----
CSV_PATH = "./output/data_dependent_on_radius_tumor.csv"
MINUTE = 300  # change this to the minute you want to plot
ATTRIBUTE = "average_oxygen_all_cells_radius"  # change this to the attribute you want to plot
# other examples: "average_oncoprotein_radius", "average_oxygen_cancer_cells_radius",
# "num_alive_cells_radius", "num_alive_tumor_cells_radius", etc.

# ---- load data ----
df = pd.read_csv(CSV_PATH)
row = df[df["total_minutes"] == MINUTE].iloc[0]

# ---- find columns for the chosen attribute ----
pattern = re.compile(rf"^{re.escape(ATTRIBUTE)}_(\d+)_to_(\d+)$")
points = []
for col in df.columns:
    m = pattern.match(col)
    if m and pd.notna(row[col]):
        x_start, x_end = int(m.group(1)), int(m.group(2))
        radius = (x_start + x_end) / 2
        points.append((radius, row[col]))

if not points:
    raise SystemExit(f"No columns found matching attribute '{ATTRIBUTE}'")

points.sort()
radii = [p[0] for p in points]
values = [p[1] for p in points]

print(f"Attribute: {ATTRIBUTE}")
print(f"Minute: {MINUTE}")
print(f"Radius range with data: {radii[0]} to {radii[-1]}")
print(f"Number of points: {len(radii)}")

# ---- plot ----
plt.figure(figsize=(8, 5))
plt.plot(radii, values, marker="o")
plt.xlabel("Radius (um)")
plt.ylabel(ATTRIBUTE)
plt.title(f"{ATTRIBUTE} vs radius - minute {MINUTE}")
plt.grid(True, alpha=0.3)
plt.tight_layout()
plt.savefig("./draft/oxygen_vs_radius.png", dpi=150)
print("Saved plot to oxygen_vs_radius.png")