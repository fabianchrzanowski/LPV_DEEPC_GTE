# AGTF30 Engine Controller — LPV-MPC & LPV-DeePC

MATLAB/Simulink project implementing gain-scheduled predictive controllers (LPV-MPC and LPV-DeePC) on the AGTF30 twin-spool turbofan engine model (NASA). Requires T-MATS solvers on the MATLAB path.

---

## Prerequisites

- [YALMIP](https://yalmip.github.io/) + MOSEK (for LPV-MPC)
- [CasADi](https://web.casadi.org/) + IPOPT (for LPV-DeePC)
- Open `GTE.prj` before running any script — all scripts resolve paths through `proj.RootFolder`.

System identification must be completed first. Run scripts in `Library/system_identification/` in order (`AGTF30_sys_id_1.m` → `AGTF30_sys_id_4.m`) to generate the LPV model tables in `Params/tables/`.

---

## Project Structure

```
Library/
├── controllers/
│   ├── setup_simulation_params.m                          ← master config
│   ├── AGTF30_lpv_mpc_yalmip_closed_loop_nonlinear.m
│   └── AGTF30_lpv_deepc_casadi_fast_closed_loop_nonlinear.m
└── Plotting/
    ├── plot_options.m                                     
    ├── plot_controller_comparison.m                       ← main comparison dashboard
    ...
```

---

## Simulation Parameters — `setup_simulation_params.m`

All controllers share this single config file.
**Scenario options:**

- `'QUICK_TEST'` — short two-segment step, useful for debugging
- `'MID_TEST'` — ~600-step trajectory: idle → ramp → cruise → high power → toggles
- `'LONG_TEST'` — aggressive ramps, rapid steps, sinusoidal inputs

---

## Critical Controller Scripts

### LPV-MPC — `AGTF30_lpv_mpc_yalmip_closed_loop_nonlinear.m`

Gain-scheduled MPC using YALMIP + MOSEK. At each time step, selects the nearest LPV operating point and solves a QP.

---

### LPV-DeePC — `AGTF30_lpv_deepc_casadi_fast_closed_loop_nonlinear.m`

Data-driven predictive controller. Builds a separate Hankel matrix per Wf operating point and schedules between them at runtime. Uses CasADi + IPOPT.

Key tuning parameters:

```matlab
u_min = 0.35;  u_max = 2.0;   du_max = 0.15;
q_track = 15;  r_u = 1;       r_du = 1;
lambda_g     = 50;    % regularisation on Hankel coefficients
lambda_sigma = 10;    % slack for DC offset mismatch
max_cols_per_wf = 3000;  % Hankel columns per operating point
```
---

## Running a Simulation

1. Open `GTE.prj` in MATLAB.
2. Optionally set variables before running:
   ```matlab
   SCENARIO = 'MID_TEST'; % or other trajectory
   % and possibly other variables
   ```
3. Run the desired controller script — it calls `setup_simulation_params.m`.
4. Results are saved to `Results/results_<controller_name>.mat`.

---

## Plotting

All plotting scripts are in `Library/Plotting/`.

Run `plot_controller_comparison.m` after any controller script — it lists available results in `Results/` and generates overlay plots interactively.
All other scripts explain what they plot in the header comment.
