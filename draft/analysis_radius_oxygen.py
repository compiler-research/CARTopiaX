import re
import pandas as pd
import matplotlib.pyplot as plt

# ---- settings ----
CSV_PATH = "./output/data_dependent_on_radius_tumor.csv"
MINUTE = 3780  # change this to the minute you want to plot

# ---- load data ----
df = pd.read_csv(CSV_PATH)
print(df.head())
row = df[df["total_minutes"] == MINUTE].iloc[0]

# ---- find radius columns ----
pattern = re.compile(r"^average_oxygen_all_cells_radius_(\d+)_to_(\d+)$")
points = []
for col in df.columns:
    m = pattern.match(col)
    if m and pd.notna(row[col]):
        x_start, x_end = int(m.group(1)), int(m.group(2))
        radius = (x_start + x_end) / 2
        points.append((radius, row[col]))

points.sort()
radii = [p[0] for p in points]
oxygen = [p[1] for p in points]

print(f"Minute: {MINUTE}")
print(f"Radius range with data: {radii[0]} to {radii[-1]}")
print(f"Number of points: {len(radii)}")

# ---- plot ----
plt.figure(figsize=(8, 5))
plt.plot(radii, oxygen, marker="o")
plt.xlabel("Radius (um)")
plt.ylabel("Average oxygen (all cells)")
plt.title(f"Average oxygen vs radius - minute {MINUTE}")
plt.grid(True, alpha=0.3)
plt.tight_layout()
plt.savefig("oxygen_vs_radius.png", dpi=150)
plt.show()