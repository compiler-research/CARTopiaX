import re
import pandas as pd
import matplotlib.pyplot as plt
import re
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt

# ---- settings ----
CSV_PATH = "./output/data_dependent_on_radius_tumor.csv"
# MINUTE = 4320 # change this to the minute you want to plot
MINUTE = 4320 # change this to the minute you want to plot
# ATTRIBUTE = "average_oxygen_all_cells_radius"  # change this to the attribute you want to plot
ATTRIBUTE = "tumor_cells_type5_dead_radius"  # change this to the attribute you want to plot
# DIVIDING_FACTOR = 585  # change this to the factor you want to divide by. Use 525 to pass from mmHg to mol/m3
DIVIDING_FACTOR = 1  # change this to the factor you want to divide by. Use 525 to pass from mmHg to mol/m3
# other examples: "average_oncoprotein_radius", "average_oxygen_cancer_cells_radius",
# "num_alive_cells_radius", "num_alive_tumor_cells_radius", etc.
APPLY_RADIAL_NORMALIZATION = True  # change this to True if you want to apply radial normalization to the attribute

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
        points.append((radius, x_start, x_end, row[col]))

if not points:
    raise SystemExit(f"No columns found matching attribute '{ATTRIBUTE}'")

points.sort(key=lambda p: p[0])
radii = [p[0] for p in points]
raw_values = np.array([p[3] for p in points], dtype=float)

if APPLY_RADIAL_NORMALIZATION:
    r_in = np.array([p[1] for p in points], dtype=float)
    r_out = np.array([p[2] for p in points], dtype=float)
    ring_area = np.pi * (r_out**2 - r_in**2)
    raw_values = raw_values / ring_area

values = raw_values / DIVIDING_FACTOR

print(f"Attribute: {ATTRIBUTE}")
print(f"Minute: {MINUTE}")
print(f"Radius range with data: {radii[0]} to {radii[-1]}")
print(f"Number of points: {len(radii)}")
print(f"Radial normalization applied: {APPLY_RADIAL_NORMALIZATION}")

# ---- plot ----
plt.figure(figsize=(8, 5))
plt.plot(radii, values, marker="o")
plt.xlabel("Radius (um)")
ylabel = ATTRIBUTE + (" (per um^2)" if APPLY_RADIAL_NORMALIZATION else "")
plt.ylabel(ylabel)
plt.title(f"{ATTRIBUTE} vs radius - minute {MINUTE}")
plt.grid(True, alpha=0.3)
plt.tight_layout()
plt.savefig("./draft/" + ATTRIBUTE + "_vs_radius.png", dpi=150)
print("Saved plot to " + ATTRIBUTE + "_vs_radius.png")