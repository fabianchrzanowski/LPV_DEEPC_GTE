# Dissertation: Retrofitting Data-Driven Feedback into Feedforward Controllers

## The Narrative Arc

Your AGTF30 engine already has a **feedforward controller** — it's in `setup_Controller.m`. It maps:

```
PLA (throttle) → fuel schedule lookup tables → Wf
```

plus limiters for N1_max, N3_max, T45_max, Ps3 bounds, and Wf/Ps3 ratio limits. This is a classic **open-loop schedule with safety limits**. It has no model of the plant — it just follows pre-computed tables.

**Your dissertation adds a feedback layer on top using data-driven models.** The story is:

```
Old: PLA → Feedforward schedule → Wf → Engine
New: PLA → Feedforward schedule → Wf_nominal
                                      ↓
     Reference → [MPC using data-driven LPV model] → ΔWf correction
                                      ↓
                              Wf_nominal + ΔWf → Engine → Measurement → (back to MPC)
```

This is the "retrofit" — you're not replacing the feedforward, you're augmenting it with feedback.

---

## What You Already Have

| Component | Script | Status |
|-----------|--------|--------|
| **Feedforward controller** | `setup_Controller.m` (Pablo's) | ✅ Exists, built into Simulink model |
| **System identification** | `sys_id_1` → `sys_id_2` | ✅ Working |
| **LPV model construction** | `sys_id_4` + `interp_matrix_eig` | ✅ Working (eigenvalue interp) |
| **Single-point MPC** | `AGTF30_mpc_yalmip_closed_loop_nonlinear` | ✅ Pablo's, works at one Wf |
| **LPV MPC (element-wise)** | `AGTF30_lpv_mpc_yalmip_closed_loop_nonlinear` | ✅ Pablo's, gain-scheduled |
| **LPV MPC (eigenvalue)** | `AGTF30_eig_lpv_mpc_trajectory` | ✅ Ours, needs Simulink |
| **DeePC** | `AGTF30_deepc_yalmip_closed_loop_nonlinear` | ✅ Pablo's, model-free alternative |
| **Interpolation comparison** | `AGTF30_compare_interpolation_methods` | ✅ Ours, dissertation figures |

---

## Suggested Dissertation Structure

### Chapter 1: Introduction
- Gas turbine control problem: why feedforward schedules are insufficient
- Motivation: disturbance rejection, constraint handling, adaptation
- Thesis statement: "Can data-driven models enable feedback control retrofit without physical modelling?"

### Chapter 2: Background / Literature Review
- Feedforward fuel scheduling (the current industrial practice)
- System identification for nonlinear systems (LPV modelling)
- MPC fundamentals
- Eigenvalue interpolation for LPV stability (cite De Caigny, Al-Jiboory & Zhu)

### Chapter 3: System Identification
- PRBS experiment design (`sys_id_1`)
- Local LTI model identification via regression (`sys_id_2`)
- **Key contribution**: eigenvalue-based interpolation with consistent tracking vs element-wise
- Validation against nonlinear plant (`sys_id_3`, `sys_id_5`)
- Figures: eigenvalue comparison plots, prediction accuracy

### Chapter 4: MPC Controller Design
- **Single-point MPC**: works at one operating point (baseline)
- **LPV gain-scheduled MPC**: updates model at each step based on current Wf
- Constraint formulation: fuel bounds, rate limits
- (Optional) DeePC comparison: model-free vs model-based

### Chapter 5: Results & Comparison
This is where the "retrofit" story comes together:

| Experiment | What it shows |
|-----------|--------------|
| **Feedforward only** | Engine follows schedule but can't reject disturbances or track arbitrary references |
| **Feedforward + single-point MPC** | Works near one operating point, degrades elsewhere |
| **Feedforward + LPV MPC (element-wise)** | Works across range, but potential stability issues with sparse data |
| **Feedforward + LPV MPC (eigenvalue)** | Works across range with stability transparency |
| **(Optional) DeePC** | Model-free alternative — how does it compare? |

### Chapter 6: Discussion
- Stability-accuracy tradeoff (fixed V limitation)
- Data efficiency: how many operating points are needed?
- Computational feasibility: can MPC run in real-time at Ts = 15ms?
- Limitations and future work

---

## What's Missing — Ranked by Impact

### 1. Feedforward-only baseline (HIGH IMPACT, EASY)

You need to show what the engine does **without** your MPC — just the built-in feedforward controller responding to PLA changes. This is your "before" picture.

Run the Simulink model with a PLA profile and record the outputs. No MPC code needed — just `setup_everything` and `sim('AGTF30SysDyn')`.

### 2. Get closed-loop MPC running on Simulink (HIGH IMPACT, MEDIUM)

The T-MATS path issue. Try running Pablo's existing scripts first:
1. `AGTF30_mpc_yalmip_closed_loop_nonlinear` — if this works, Simulink is fine
2. `AGTF30_lpv_mpc_yalmip_closed_loop_nonlinear` — Pablo's LPV version
3. `AGTF30_eig_lpv_mpc_trajectory` — our eigenvalue version

### 3. Side-by-side comparison table (HIGH IMPACT, MEDIUM)

Run all controllers on the **same reference trajectory** and tabulate:

| Metric | FF only | Single MPC | LPV MPC (elem) | LPV MPC (eig) |
|--------|---------|-----------|-----------------|----------------|
| RMSE (N_H) | — | ? | ? | ? |
| Max error | — | ? | ? | ? |
| Wf usage | — | ? | ? | ? |
| Solve time/step | 0 | ? | ? | ? |
| Stability guaranteed | N/A | Local | ❌ | Transparent |

### 4. Sparse data experiment (MEDIUM IMPACT, EASY)

Already built into the comparison script. Run with `sparse_step = 1, 3, 5, 7` and show how element-wise degrades while eigenvalue interpolation maintains quality. This demonstrates data efficiency.

### 5. (OPTIONAL) Disturbance rejection test (MEDIUM IMPACT, MEDIUM)

Change altitude or temperature mid-trajectory:
```matlab
MWS.In.Alt = timeseries([0; 0; 5000; 5000], [0, 5, 5.1, 20], 'Name', 'Alt');
```
Feedforward can't handle this; MPC should adapt.

### 6. (OPTIONAL) Surge margin constraints (LOW IMPACT, EASY)

Add to MPC constraints:
```matlab
% z_pred indices 14,15 are SMFan and SMHPC
constraints = [constraints, z_pred(14) >= 5, z_pred(15) >= 5];
```
Shows MPC can enforce safety limits that feedforward limiters handle less elegantly.

### 7. (OPTIONAL) DeePC comparison (LOW IMPACT, HARD)

Pablo already wrote the DeePC script. If it works, run it on the same trajectory for a model-free vs model-based comparison. Note: the DeePC script says "takes forever, doesn't work that great" in its comments (line 4), so this may not be worth the effort.

---

## How to Run Everything (Order)

```
1. Double-click GTE.prj                              → loads paths
2. AGTF30_sys_id_4                                    → builds LUT (no Simulink)
3. AGTF30_compare_interpolation_methods (Fig 1 only)  → dissertation figures (no Simulink)
4. AGTF30_mpc_yalmip_closed_loop_nonlinear            → test Simulink works (Pablo's)
5. AGTF30_eig_lpv_mpc_trajectory                      → our eigenvalue MPC (Simulink)
6. Tabulate results from 4 and 5
```
