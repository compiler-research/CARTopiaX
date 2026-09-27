"""Run one CARTopiaX simulation in-process from Python.

Parameters are set directly on a SimParam (no params.json), and the summary
statistics come back in memory as SummaryRow objects, with the same columns
as output/final_data.csv. Build the project first, then run from the
repository root.
"""

import os
#single threaded by default, so the final results are bit reproducible.
os.environ.setdefault("OMP_NUM_THREADS", "1")

import hashlib
import cppyy

cppyy.add_include_path(os.path.abspath("src"))
cppyy.load_library("build/libCARTopiaX.so")
cppyy.include("params/hyperparams.h")
cppyy.include("cart_tumor.h")


cppyy.cppdef('const char* g_cartopiax_argv[] = {"CARTopiaX"};')

params = cppyy.gbl.bdm.SimParam()
params.__python_owns__ = False
params.seed = 42
params.output_performance_statistics = False
params.total_minutes_to_simulate = 1440   
params.initial_tumor_radius = 150.0
params.bounded_space_length = 1000
params.resolution_grid_substances = 50
params.dt_substances = 0.01
params.dt_step = 0.1
params.treatment.clear()
params.treatment[0] = 3957
params.ComputeDerived()
params.output_csv_interval = 600           

rows = cppyy.gbl.std.vector[cppyy.gbl.bdm.SummaryRow]()
cppyy.gbl.bdm.Simulate(1, cppyy.gbl.g_cartopiax_argv, params, rows)

:if len(rows) == 0:
    raise RuntimeError("Simulate returned no rows")
last = rows.back()
print("rows:", len(rows), "| tumour cells:", last.num_tumor_cells,
      "| radius:", last.tumor_radius, flush=True)

h = hashlib.sha256(open("output/final_data.csv", "rb").read()).hexdigest()
print("sha256:", h, flush=True)

