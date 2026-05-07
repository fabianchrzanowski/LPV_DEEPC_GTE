# What Each File Does — Detailed Guide

## File 1: `interp_matrix_eig.m` — The Core Algorithm

**Where:** `Library/utilities/`
**Needs Simulink?** ❌ No. Pure maths.

### What it does in plain English

You have ~32 identified A matrices, one per Wf operating point. You want a smooth A matrix at *any* Wf — not just the 32 you identified at. This function does that interpolation.

**The problem it solves:** If you just interpolate each of the 16 entries of the 4×4 A matrix separately (which is what `interp_matrix.m` does), the eigenvalues of the result can accidentally leave the unit circle between operating points → the model goes unstable at some Wf values even though it was stable at all 32 identified points.

### How it works step by step

1. **Eigen-decompose** each of the 32 identified A matrices: `A_k = V_k · diag(λ_k) · V_k⁻¹`
2. **Sort eigenvalues consistently** across all operating points (by magnitude then angle) so eigenvalue #1 at Wf=0.5 corresponds to eigenvalue #1 at Wf=0.55 etc.
3. **Polar interpolation** — convert each eigenvalue to magnitude + angle, unwrap the angles to avoid 2π jumps, interpolate both with `pchip` across the fine ρ grid
4. **Reconstruct** — pick eigenvectors from the midpoint operating point (V_nom), and rebuild: `A(ρ) = V_nom · diag(λ(ρ)) · V_nom⁻¹`
5. **Flag** any ρ where |λ| > 1 as a warning (doesn't hide it, doesn't clamp)

### Inputs/outputs

```matlab
% IN:  rho points (32 values), A matrices (4×4×32), query grid (1600 values)
% OUT: A matrices (4×4×1600), diagnostic info struct
[A_interp, info] = interp_matrix_eig(rho, A_raw, rho_fine, opts);
```

### Limitations (to acknowledge in thesis)

- **Fixed eigenvector basis** introduces reconstruction error at operating points far from the midpoint. The full SMILE approach (De Caigny et al., 2011) would interpolate V(ρ) too, but is significantly more complex.
- **Sorting-based tracking** may swap eigenvalue identity at mode veering points. More sophisticated approaches like MAC-based matching (Al-Jiboory & Zhu, 2019) exist but add complexity.

---

## File 2: `AGTF30_sys_id_4.m` (modified) — Builds the LUT

**Where:** `Library/system_identification/`
**Needs Simulink?** ❌ No. Just loads `.mat` files and does interpolation.

### What it does

Takes the raw identified matrices from `sys_id_2` (stored in `Params/tables/matrices_DE_wf.mat`) and interpolates them onto a fine ρ grid:

| Matrix | Interpolation method | Why |
|--------|---------------------|-----|
| **A** | `interp_matrix_eig` (eigenvalue) | A's eigenvalues determine stability |
| **B, C, D** | `interp_matrix` (element-wise) | These don't affect stability |
| **Offsets** (x₀, y₀, u₀, z₀) | `interp1` with pchip | Smooth scalar/vector values |

**This is the file that turns your 32 identified models into a single continuous LPV model.**

Saves the result as a lookup table in:
`Params/tables/models/AGTF30_model_tables__1_input.mat`

Also generates comparison plots: eigenvalues before vs after interpolation.

---

## File 3: `AGTF30_compare_interpolation_methods.m` — Dissertation Figures

**Where:** `Library/system_identification/`
**Needs Simulink?** ⚠️ **Partially.**

### What it does

Builds **two** A-matrix interpolations from a deliberately **sparse** subset of operating points (every 5th → only ~7 instead of 32):

1. **Element-wise** — standard `interp_matrix` (what you had before)
2. **Eigenvalue** — `interp_matrix_eig`

Then generates 3 figures:

| Figure | Needs Simulink | What it shows |
|--------|:-:|---|
| **Fig 1**: eigenvalue magnitudes + complex plane | ❌ | Where element-wise goes unstable with sparse data |
| **Fig 2**: time traces vs Simulink plant | ✅ | Prediction accuracy comparison |
| **Fig 3**: error histograms | ✅ | Statistical accuracy comparison |

### To run ONLY Fig 1 (no Simulink)

Add `return` after the Fig 1 section (~line 155):
```matlab
fprintf('Stopping before Simulink validation.\n'); return;
```

### Tuning the stress test

Change `sparse_step` at line 46:
- `sparse_step = 1` → all points, both methods work well (baseline)
- `sparse_step = 5` → stress test, element-wise more likely to go unstable
- `sparse_step = 7` → aggressive, strongest difference between methods

---

## File 4: `AGTF30_eig_lpv_mpc_trajectory.m` — The MPC Controller

**Where:** `Library/controllers/`
**Needs Simulink?** ✅ **Yes**

### What it does

Closed-loop MPC controller. At each time step:

1. **Reads current state** from the Simulink plant (N_H or N_L measurement)
2. **Looks up A, B, C, D** from the eigenvalue-interpolated LUT at current ρ (= current fuel flow)
3. **Solves a QP** (quadratic program via YALMIP) to find optimal fuel flow for the next N_pred steps, subject to:
   - Fuel flow bounds: `u_min ≤ Wf ≤ u_max`
   - Rate constraint: `|ΔWf| ≤ du_max`
   - Tracking objective: minimise `(output - reference)²`
4. **Applies first fuel command** to the Simulink plant
5. **Repeats**

### Trajectory options (line ~74)

```matlab
trajectory_type = 'takeoff_profile';
```

| Option | Description |
|--------|-------------|
| `'sinusoidal'` | Smooth 300 rpm oscillation around nominal |
| `'step_sequence'` | Big jumps: idle → cruise → max → cruise → idle |
| `'takeoff_profile'` | Realistic: idle → spool up → climb → cruise |
| `'custom'` | Your own — edit the code block |

---

## How to Test WITHOUT Simulink

| Test | Command | What you see |
|------|---------|-------------|
| Eigenvalue plots | `AGTF30_sys_id_4` | Smooth eigenvalue curves across ρ |
| Stability comparison | `AGTF30_compare_interpolation_methods` (Fig 1 only) | Element-wise vs eigenvalue |
| MPC on linear model | `AGTF30_simple_mpc_yalmip` (lines 77-113) | MPC tracking without Simulink plant |

---

## How to Run in Simulink

The dependency chain:

```
Your script → AGTF30_initial_conditions(MWS)
    → sim('AGTF30SysSS')
        → AGTF30_eng.mdl
            → S-function "Ambient_TMATS"
                → T-MATS toolbox must be on MATLAB path
```

### Fix the T-MATS path

```matlab
% 1. Open the project (this should set up paths):
openProject('path/to/GTE-master - Copy')

% 2. Check if T-MATS is loaded:
which Ambient_TMATS
% If "not found" → T-MATS didn't copy over with the project

% 3. Check the ORIGINAL GTE-master folder:
%    T-MATS may be a referenced project that didn't copy

% 4. Find and add manually:
addpath(genpath('C:\path\to\T-MATS'))

% 5. Verify:
which Ambient_TMATS   % should now show a path
```

---

## Summary: What Needs What

| File | Simulink | YALMIP | T-MATS | Output |
|------|:---:|:---:|:---:|---|
| `interp_matrix_eig.m` | ❌ | ❌ | ❌ | (utility function) |
| `sys_id_4` | ❌ | ❌ | ❌ | LUT `.mat` file + eigenvalue plots |
| `compare...` (Fig 1) | ❌ | ❌ | ❌ | Stability comparison figure |
| `compare...` (Figs 2-3) | ✅ | ❌ | ✅ | Accuracy comparison figures |
| `eig_lpv_mpc_trajectory` | ✅ | ✅ | ✅ | Closed-loop tracking results |
