# Mathematical Formulations of the Control Architectures

This document details the exact mathematics governing the four control architectures compared in your dissertation. All methods aim to minimize a quadratic tracking cost function over a prediction horizon $N$, but they differ fundamentally in how they predict the future behaviour of the nonlinear AGTF30 engine.

---

## 1. The Core Problem: Nonlinearity
The real AGTF30 engine is a highly nonlinear system. Its dynamics can be represented generally as:
$$ x_{k+1} = f(x_k, u_k) $$
$$ y_k = h(x_k, u_k) $$

Where $u$ is the fuel flow ($W_f$), $y$ are the measured outputs (e.g., $N_H, N_L$), and $x$ are the internal states. Because optimal control of a generic nonlinear system is computationally intractable for real-time engine control, we rely on linear approximations.

### Linearisation and Eigenvalues

To apply linear control theory, we must linearise the nonlinear dynamics around an equilibrium point $(x_{eq}, u_{eq})$. Mathematically, this is done using a Taylor series expansion, retaining only the first-order Jacobian matrices:

$$ A = \left. \frac{\partial f}{\partial x} \right|_{(x_{eq}, u_{eq})}, \quad B = \left. \frac{\partial f}{\partial u} \right|_{(x_{eq}, u_{eq})} $$

This gives us the local Linear Time-Invariant (LTI) model in deviation variables:
$$ \delta x_{k+1} = A \delta x_k + B \delta u_k $$

**Eigenvalues and Local Stability:**
The mathematical key to understanding whether this local model is stable lies in the **eigenvalues** of the state transition matrix $A$. By solving the characteristic equation $\det(A - \lambda I) = 0$, we find the eigenvalues $\lambda_i$. 

Because your controllers operate in discrete-time (sampled at $T_s = 0.015s$), the mathematical condition for local asymptotic stability is that all eigenvalues must lie strictly inside the unit circle on the complex plane:
$$ |\lambda_i| < 1 \quad \text{for all } i $$

If the magnitude of any eigenvalue exceeds 1 (i.e., $\max(|eig(A)|) > 1$), the local linear model will diverge to infinity. In your engine system identification process, it's normal for numerical inaccuracies or nonlinear quirks at certain fuel flows to produce marginally unstable local models (e.g., $|\lambda| = 1.01$). This is why we mathematically filter or ignore these specific operating points when designing the controller or simulating the LPV plant.

---

## 2. Frozen Model Predictive Control (Linear MPC)

**Concept:** 
Assume the engine operates near a single, fixed equilibrium point (e.g., Cruise, $W_{f,eq} = 1.2$ pps). We linearise the engine equations via a Taylor expansion at this point to yield a single Linear Time-Invariant (LTI) model.

**Prediction Model:**
The future states are predicted using constant matrices $A, B, C, D$:
$$ \delta x_{k+1} = A \delta x_k + B \delta u_k $$
$$ \delta y_k = C \delta x_k + D \delta u_k $$

Where the variables are deviations from the equilibrium:
- $\delta u_k = u_k - u_{eq}$
- $\delta y_k = y_k - y_{eq}$

**Optimisation Problem (QP):**
At every time step $t$, the controller solves:
$$ \min_{\delta u} \sum_{i=1}^{N} \Big( q \| y_{t+i} - r_{t+i} \|^2_2 + r_{du} \| \Delta u_{t+i} \|^2_2 \Big) $$
Subject to:
- Model dynamics: $x_{k+1} = Ax_k + Bu_k$
- Input limits: $u_{min} \leq u_k \leq u_{max}$
- Rate limits: $-du_{max} \leq \Delta u_k \leq du_{max}$

> [!WARNING]
> **The Flaw:** $A$ and $B$ are only valid near $u_{eq}$. When you command Max thrust ($W_f = 1.8$), the physical engine dynamics change dramatically (e.g., time constants shrink, DC gain shifts). The Frozen MPC is entirely blind to this, predicting the future using incorrect cruise dynamics, which leads to massive tracking errors.

---

## 3. Linear Parameter-Varying MPC (LPV MPC)

**Concept:**
Instead of a single LTI model, the engine is modeled as a family of linear models scheduled by a measurable parameter $\rho$ (in this case, $\rho = W_f$). 

**Prediction Model:**
$$ x_{k+1} = A(\rho_k) x_k + B(\rho_k) u_k $$
$$ y_k = C(\rho_k) x_k + D(\rho_k) u_k $$

In your implementation, the scheduling variable is the *current* fuel flow: $\rho_k = u_{k-1}$. At the start of every control step, the algorithm interpolates the look-up tables to find the exact $A, B, C, D$ and equilibrium offsets for the *current* operating point.

**Optimisation Problem:**
The QP is identical to Frozen MPC, but the matrices $A, B$ are updated at every time step $t$. 

> [!NOTE]
> **Implementation detail:** To keep the QP convex and solvable in milliseconds, your `AGTF30_lpv_mpc` assumes the matrices $A(\rho_{t}), B(\rho_{t})$ remain constant *across the prediction horizon $N$*. The matrices change at every real-time step $t$, but during the $N$-step prediction, the solver assumes "frozen-at-current-step" dynamics. This is a standard and highly effective LPV approximation.

---

## 4. Frozen Data-Enabled Predictive Control (Linear DeePC)

**Concept:**
DeePC skips the system identification step entirely. Instead of identifying $A$ and $B$ matrices, it uses raw input/output data directly to predict the future. This relies on **Willems' Fundamental Lemma**, which states that for a linear, controllable system, any valid trajectory of length $L$ can be expressed as a linear combination of previously recorded data.

**The Hankel Matrix:**
We collect a long dataset of length $T_d$ by applying a persistently exciting input $u_d$ to the engine and measuring $y_d$. We stack this data into a block Hankel matrix, which is then partitioned into "past" (length $T_{ini}$) and "future" (length $N$) blocks:

$$ \begin{bmatrix} U_p \\ U_f \\ Y_p \\ Y_f \end{bmatrix} g = \begin{bmatrix} u_{ini} \\ u_{pred} \\ y_{ini} \\ y_{pred} \end{bmatrix} $$

Where $g$ is a vector of multipliers. By forcing the past data $U_p, Y_p$ to match our *actual recent history* $u_{ini}, y_{ini}$, the future data $U_f, Y_f$ naturally forms a valid prediction of the future.

**Optimisation Problem:**
$$ \min_{g, \sigma_y} \sum_{i=1}^{N} \Big( q \| y_{pred, i} - r_{i} \|^2_2 + r_{du} \| \Delta u_{pred, i} \|^2_2 \Big) + \lambda_g \|g\|^2_2 + \lambda_\sigma \|\sigma_y\|^2_2 $$

Subject to:
- Data constraint: $\begin{bmatrix} U_p \\ Y_p \end{bmatrix} g = \begin{bmatrix} u_{ini} \\ y_{ini} + \sigma_y \end{bmatrix}$
- Future extraction: $u_{pred} = U_f g$, $y_{pred} = Y_f g$
- Standard input/rate limits on $u_{pred}$

*Note: $\sigma_y$ is a slack variable that absorbs measurement noise so the equality constraint doesn't become infeasible.*

> [!WARNING]
> **The Flaw:** Willems' Lemma *only holds for LTI systems*. If you build the Hankel matrix using data from a single operating point (Frozen DeePC), it will fail during transients exactly like Frozen MPC. If you build the Hankel matrix by concatenating data from *all* operating points, you mathematically violate the lemma by mixing incompatible nonlinear dynamics, leading to the "flat" unresponsive behavior you saw earlier.

---

## 5. LPV DeePC (Gain-Scheduled DeePC)

**Concept:**
This is the state-of-the-art solution that bridges data-driven control and nonlinear adaptability. Since Willems' Lemma requires LTI data, we cannot mix data from different power settings. Instead, we generate *multiple* independent Hankel matrices, one for each steady-state operating point $W_f$.

**The Mechanism:**
1. **Offline:** Identify a discrete set of Hankel matrices $\mathcal{H}_1, \mathcal{H}_2, \dots, \mathcal{H}_m$ at various fuel flows. Formulate and pre-compile a separate CasADi QP solver for each Hankel matrix.
2. **Online:** At time step $t$, measure the current fuel flow $W_f$.
3. **Schedule:** Select the pre-compiled solver $j$ whose base $W_f$ is closest to the current state.
4. **Solve:** Pass the recent history $(u_{ini}, y_{ini})$ to solver $j$ to compute the optimal $u_{pred}$.

**Optimisation Problem:**
Identical to Linear DeePC, but the underlying data matrices $(U_p, U_f, Y_p, Y_f)_j$ switch dynamically based on the scheduling parameter.

> [!TIP]
> **Why it's brilliant:** LPV-DeePC completely eliminates the need to identify complex parameter-varying matrices $A(\rho), B(\rho)$. It naturally captures the local system dynamics (including unmodeled linear dynamics that parametric system ID might miss) while adapting to the global nonlinearity of the engine envelope.

---

## 6. Fast Surrogate Mode (Plant Simulation)

When evaluating controller performance, the true AGTF30 Simulink model acts as the real-world nonlinear engine ("the plant"). However, simulating a 1-hour fluid-dynamics engine model for every single step of a receding-horizon controller takes an immense amount of computational time.

To accelerate validation, we employ a **Fast Surrogate Mode**. 

### The Mathematics of the Surrogate

Instead of calling the high-fidelity Simulink physics engine to find out how the engine reacts to a fuel command, we simulate the engine mathematically using the set of identified Linear Parameter-Varying (LPV) state-space matrices: $\mathcal{A}(\rho), \mathcal{B}(\rho), \mathcal{C}(\rho), \mathcal{D}(\rho)$.

For any given step $t$, we observe the *actual* fuel flow from the previous step $u_{t-1}$. We use this to interpolate the global nonlinear engine dynamics into a localized linear matrix.

The surrogate plant updates its internal "fast state" $x_{fast}$ using the exact state-space difference equations:

$$ x_{fast}[t+1] = A_p(u_{t-1}) \cdot x_{fast}[t] + B_p(u_{t-1}) \cdot (u_{now} - u_{0,p}) $$
$$ y_{meas}[t] = z_{0,p} + C_p(u_{t-1}) \cdot x_{fast}[t] + D_p(u_{t-1}) \cdot (u_{prev} - u_{0,p}) $$

Where:
- $A_p, B_p, C_p, D_p$ are the plant matrices evaluated at the closest equilibrium point to the current fuel flow.
- $u_{0,p}, z_{0,p}$ are the steady-state trim offsets at that operating point.
- $u_{now}$ is the control action just decided by the controller.
- $y_{meas}$ is the simulated output (e.g., $N_H$ RPM) fed back into the controller.

### Why it works

Because we are doing pure matrix multiplications in MATLAB instead of spinning up a Simulink solver, the computation time drops from roughly 30 minutes to just **3 seconds**. 

This allows us to rapidly prototype, debug, and tune the controllers (like adjusting $N_{pred}$ or tracking weights $Q$) in near real-time. Once the controller is perfectly tuned and proven stable on the mathematical surrogate, `FAST_SURROGATE_MODE` can be turned to `false` for a final, high-fidelity verification run against the true Simulink physics engine.
