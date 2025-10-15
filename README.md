Rendezvous (SCP + QCQP) ??Quick Guide

Overview
- MATLAB-based orbit-rendezvous experiments using the Gim?“Alfriend STM (with J2) in LVLH.
- Includes: Successive Convexification (SCP) with keep-out zones (KOZ), and QCQP/QP/LP baselines.

Requirements
- MATLAB (R2021a+ recommended)
  - Optimization Toolbox (quadprog, fmincon)
  - Symbolic Math Toolbox (optional: enables exact grads for implicit surfaces; falls back to numeric diffs)
- Optional (for QCQP): Gurobi with MATLAB interface installed and on the MATLAB path
- Optional (helpers): Python 3 + NumPy

Repo Layout
- `GA-STM/` MATLAB classes and scripts (GimAlfriendSTM, SCP, baselines)
- `startup.m` Adds `GA-STM/` recursively to MATLAB path
- `run_qcqp_custom.m` QCQP-like run with fallback to `fmincon`
- `custom_states.csv` Two rows `[x vx y vy z vz]` in meters (start/goal)
- `prepare_custom_states.py` Helper to build `custom_states.csv` from `.npy`

Getting Started (MATLAB)
1) Open MATLAB in this repo folder (or `cd` here) and run `startup` to add paths.
2) Choose one of the pipelines below.

Run: SCP with Implicit KOZ
- Script: `GA-STM/SCP_KOZ_QP.m`
- What it does:
  - Builds discrete STM from chief elements, sets a time grid, and runs SCP.
  - Imports multiple implicit KOZ level sets `F(x,y,z)-1=0` from external files and linearizes them with gradients.
  - Saves `GA-STM/scp_last.mat` and `GA-STM/scp_koz_results.png`.
- Before running:
  - Ensure the expression directory exists (edit in the script if needed):
    - Default: `C:\Users\98kim\Desktop\Acta-Astronautica\Funcs_ISS_expr`
  - If Symbolic Toolbox is missing, it will use numeric gradients automatically.
- Run: open the script and press Run, or call `SCP_KOZ_QP`.

Run: SCP with Ellipsoidal KOZ (simple)
- Script: `GA-STM/SCP_Ellipsoid.m`
- Reads `custom_states.csv` for start/goal, enforces an ellipsoidal KOZ only.
- Produces `GA-STM/scp_ellipsoid_last.mat` and figures.

Run: MIQP with AABB KOZ (same domain as SCP)
- Script: `GA-STM/SCP_KOZ_MIQP.m` (defaults to quadprog QP-relaxation; optional Gurobi)
- What it does:
  - Same time grid and start/goal as `SCP_KOZ_QP` by default (`dt=1, N=60`, Â±50 m on T-axis).
  - L2 objective (QP-like) with mixed-integer big?‘M disjunctions to stay outside AABB approximations of implicit KOZ.
  - Saves `GA-STM/scp_miqp_last.mat` and plots.
- Configure in-file:
  - `expr_dir` (default: `C:\Users\98kim\Desktop\Acta-Astronautica\Funcs_ISS_expr`)
  - domain `[-100,100] m`, `grid_res`, and big?‘M/clearance.

Run: QCQP (custom start/goal)
- Script: `run_qcqp_custom.m`
- Input: `custom_states.csv` (two rows `[x vx y vy z vz]` in meters)
- Behavior:
  - If Gurobi is not available, falls back to `fmincon` with identical objective/constraints.
  - Saves `qcqp_custom_results.png` with trajectory/control plots and summary.
- Tip (start/goal from .npy):
  - `python prepare_custom_states.py <track.npy> custom_states.csv`
  - Prints `N_suggest=<samples>`; you can adjust `N_override` or `dt` in the script if needed.
- Note about Gurobi path and solver function:
  - `run_qcqp_custom.m` calls `solveQCQP_L2` when Gurobi is on the path.
  - `solveQCQP_L2` is currently defined as a local function inside `GA-STM/OptimizationBaselineQCQP.m`.
  - Options if you want Gurobi here:
    1) Run `GA-STM/OptimizationBaselineQCQP.m` directly; or
    2) Refactor `solveQCQP_L2` into its own `solveQCQP_L2.m` file on the MATLAB path.
  - Otherwise, leave Gurobi off the MATLAB path to use the `fmincon` fallback.

Python Helpers
- `prepare_custom_states.py <npy_path> <out_csv>`
  - Loads `.npy` (RTN positions in km), builds `[x vx y vy z vz]` rows in meters, saves CSV.
- `read_npy.py <npy_path>` Quick inspect of shape and first/last samples.
- `check_F_sign.py` Evaluates a single implicit expression file; edit the absolute path inside if needed.

Local Git
- A git repo already exists here (`.git`). To create your first commit:
  - `git add .`
  - `git commit -m "Initial import and README"`
  - Optional default branch rename: `git branch -M main`

Notes
- Large figure asset: `ISS_data3_mesh120_opac05_dull.fig` is used for optional overlay.
- Images such as `koz_surfaces_with_path.png`, `scp_koz_results.png`, `compare_qp_scp.png` land in `GA-STM/`.

PowerShell (Terminal) Commands
- Set working directory:
  - `Set-Location C:\Users\98kim\Desktop\SCP\Rendezvous`
- MIQP (QP relaxation, default):
  - `matlab -batch "startup; SCP_KOZ_MIQP"`
- MIQP with Gurobi (exact MIQP; set solver='gurobi' in script):
  - `matlab -batch "startup; addpath('C:\gurobi1202\win64\matlab'); SCP_KOZ_MIQP"`
- SCP (implicit KOZ, successive convexification):
  - `matlab -batch "startup; SCP_KOZ_QP"`
- SCP (ellipsoidal KOZ):
  - `matlab -batch "startup; SCP_Ellipsoid"`
- QCQP (custom start/goal, fmincon fallback):
  - `matlab -batch "startup; run_qcqp_custom"`
- KOZ surfaces + path visualization:
  - `matlab -batch "startup; plot_implicit_koz_and_path"`
- QP vs SCP comparison plot:
  - `matlab -batch "startup; compare_QP_SCP_plot"`
- Absolute MATLAB path example:
  - `& "C:\Program Files\MATLAB\R2023b\bin\matlab.exe" -batch "startup; SCP_KOZ_MIQP"`
