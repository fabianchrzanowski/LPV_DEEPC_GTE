# Eigenvalue Interpolation for Stable LPV Scheduling

## The Problem: Why Element-wise Interpolation Breaks

Your current code in [AGTF30_sys_id_4.m](file:///c:/Users/fafa2/OneDrive/Dokumenty/MATLAB/MATLAB%20Drive/YEAR_4/dissertation/GTE-master%20-%20Copy/Library/system_identification/AGTF30_sys_id_4.m) and [LPV2mat.m](file:///c:/Users/fafa2/OneDrive/Dokumenty/MATLAB/MATLAB%20Drive/YEAR_4/dissertation/GTE-master%20-%20Copy/Library/utilities/LPV2mat.m) does **element-wise** `interp1` on the A matrix entries across ρ (fuel flow). This gives smooth curves for each `a₁₁, a₁₂, a₂₁, a₂₂`, but the **eigenvalues** of the interpolated matrix can leave the unit circle between operating points even if every identified model is stable at its own operating point.

**Why?** Eigenvalues are a nonlinear function of the matrix entries. Linearly interpolating the entries does not linearly interpolate the eigenvalues. Two stable A matrices can produce an unstable convex combination.

---

## Your Supervisor's Approach: Interpolate Eigenvalues Instead

### The Characteristic Polynomial (2×2 case)

For your 2×2 system matrix `A` (which combines your 2 differential + 2 algebraic states into a 4×4, but the eigenvalue argument works for any size — let's use 2×2 to illustrate):

```
det(λI - A) = det | λ - a₁₁   -a₁₂  |
                  | -a₂₁    λ - a₂₂ |

             = (λ - a₁₁)(λ - a₂₂) - a₁₂·a₂₁

             = λ² - (a₁₁ + a₂₂)λ + (a₁₁·a₂₂ - a₁₂·a₂₁)

             = λ² - tr(A)·λ + det(A)
```

The eigenvalues are:
```
λ₁,₂ = tr(A)/2  ±  √( (tr(A)/2)² - det(A) )
```

So the eigenvalues are determined by **two quantities**: `tr(A)` and `det(A)`.

### The Key Idea

Instead of interpolating the four matrix entries `a₁₁, a₁₂, a₂₁, a₂₂` independently (which gives you no control over where the eigenvalues land), you:

1. **Interpolate the eigenvalues directly** (or equivalently, interpolate `tr(A)` and `det(A)`)
2. **Reconstruct** the matrix entries from the interpolated eigenvalues

This guarantees that at every intermediate ρ value, the eigenvalues are inside the unit circle (stable).

---

## How to Apply It: Concrete Recipe

### Step 1: Decompose each identified A(ρₖ)

At each operating point ρₖ, you already have an identified A matrix. Decompose it:

```matlab
[V, Lambda] = eig(A_k);   % V = eigenvectors,  Lambda = diag(eigenvalues)
% so A_k = V * Lambda * V⁻¹
```

### Step 2: Interpolate the eigenvalues

Across your operating points, interpolate the eigenvalues (not the matrix entries):

```matlab
% λ₁(ρ) and λ₂(ρ) interpolated separately
lambda1_interp = interp1(rho_points, lambda1_values, rho_query, 'pchip');
lambda2_interp = interp1(rho_points, lambda2_values, rho_query, 'pchip');
```

> [!IMPORTANT]
> If eigenvalues are complex conjugate pairs, interpolate their **magnitude** and **angle** separately (polar form) to avoid crossing the real axis incorrectly.

### Step 3: Reconstruct the matrix

For a **2×2 system**, you can recover a matrix from its eigenvalues + eigenvectors. But the eigenvectors may also change across ρ. Two practical approaches:

---

### Approach A: Fix Eigenvectors, Interpolate Eigenvalues Only (Simplest)

If the system doesn't change structure much across ρ (often true for gain scheduling around a trim line):

```matlab
% Use eigenvectors from a nominal point (e.g., middle of range)
[V_nom, ~] = eig(A_nominal);

% At any query rho:
Lambda_interp = diag([lambda1_interp, lambda2_interp]);
A_interp = V_nom * Lambda_interp / V_nom;   % V * Λ * V⁻¹
```

This guarantees the eigenvalues are exactly where you put them. If you clamp the interpolated eigenvalues inside the unit circle, the system stays stable.

---

### Approach B: Interpolate via Trace and Determinant (2×2 only, no eigenvector needed)

For a 2×2 A, the characteristic polynomial is fully determined by `tr(A)` and `det(A)`. So:

```matlab
% At each operating point k, compute:
tr_k  = trace(A_k);   % = a₁₁ + a₂₂ = λ₁ + λ₂
det_k = det(A_k);     % = a₁₁*a₂₂ - a₁₂*a₂₁ = λ₁ * λ₂

% Interpolate these two scalars:
tr_interp  = interp1(rho_points, tr_values, rho_query, 'pchip');
det_interp = interp1(rho_points, det_values, rho_query, 'pchip');
```

Now you need to recover a 2×2 matrix with the desired trace and determinant. One natural way: keep the off-diagonal elements `a₁₂, a₂₁` interpolated element-wise (they affect eigenvectors/mode shapes more than stability), and adjust the diagonal:

```matlab
% Interpolate off-diagonals normally
a12 = interp1(rho_points, a12_values, rho_query);
a21 = interp1(rho_points, a21_values, rho_query);

% Solve for diagonals from tr and det constraints:
%   a₁₁ + a₂₂ = tr_interp
%   a₁₁ * a₂₂ = det_interp + a₁₂ * a₂₁
% This is a quadratic in a₁₁:
%   a₁₁² - tr_interp * a₁₁ + (det_interp + a₁₂*a₂₁) = 0
p = det_interp + a12*a21;
a11 = (tr_interp + sqrt(tr_interp^2 - 4*p)) / 2;
a22 = tr_interp - a11;

A_interp = [a11, a12; a21, a22];
```

> [!TIP]
> **Stability check**: For a discrete-time 2×2 system, stability ⟺ all eigenvalues inside the unit circle. The conditions on `tr(A)` and `det(A)` are:
> 1. `|det(A)| < 1`
> 2. `|tr(A)| < 1 + det(A)`
>
> You can **clamp** the interpolated trace and determinant to satisfy these before reconstructing A.

---

## For Your 4×4 System (2 diff + 2 alg states)

Your actual A matrix is 4×4. The same idea generalizes — the characteristic polynomial is:

```
det(λI - A) = λ⁴ - c₃λ³ + c₂λ² - c₁λ + c₀
```

where the coefficients `c₀...c₃` are functions of the matrix entries (related to traces of powers of A). You can:

1. Compute eigenvalues at each operating point
2. Interpolate the 4 eigenvalues across ρ (keeping them inside the unit circle)
3. Reconstruct using **Approach A** (fixed nominal eigenvectors, varying eigenvalues)

Approach A scales cleanly to any dimension. Approach B (trace/det) is elegant but only practical for 2×2.

---

## Suggested Modification to Your Code

Replace the element-wise interpolation in [AGTF30_sys_id_4.m](file:///c:/Users/fafa2/OneDrive/Dokumenty/MATLAB/MATLAB%20Drive/YEAR_4/dissertation/GTE-master%20-%20Copy/Library/system_identification/AGTF30_sys_id_4.m) lines 81-84 with eigenvalue-based interpolation for the A matrix only (B, C, D don't affect stability directly):

```matlab
% --- Current (element-wise, can go unstable) ---
A = interp_matrix(rho, models.A_dt, rho_es);

% --- Proposed (eigenvalue-interpolated, stability guaranteed) ---
A = interp_matrix_eig(rho, models.A_dt, rho_es);
```

Where `interp_matrix_eig` would:
1. Eigen-decompose each `A(:,:,k)`
2. Interpolate eigenvalues (in polar form if complex)
3. Use a nominal eigenvector basis to reconstruct

> [!NOTE]
> B, C, D matrices don't have a stability constraint, so element-wise interpolation is fine for them.

---

## Summary

| Method | What's interpolated | Stability guarantee? | Complexity |
|--------|-------------------|---------------------|------------|
| Current (`interp_matrix`) | Each `aᵢⱼ` independently | ❌ No | Simple |
| Eigenvalue interp (Approach A) | λ values + fixed V | ✅ Yes (if λ clamped) | Medium |
| Trace/det interp (Approach B) | tr(A), det(A) | ✅ Yes (if clamped) | 2×2 only |
