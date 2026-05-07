% ========================================================================
% EIGENVALUE-BASED MATRIX INTERPOLATION
%
% Interpolates a 3-D page matrix (n x n x p) across the third dimension
% using eigenvalue decomposition.
%
% Pipeline:
%   1. Eigen-decompose each page A(:,:,k)
%   2. Sort eigenvalues consistently across pages (by magnitude, then angle)
%   3. Interpolate eigenvalues in polar form (magnitude, angle) with pchip
%   4. Reconstruct using fixed nominal eigenvector basis (midpoint)
%   5. Flag (don't clamp) any |lambda| > 1 as a stability warning
%
% INPUTS
%   x     : scheduling parameter values for the input pages (length p)
%   M     : matrix page object (n x n x p)
%   xi    : query points for interpolation
%   opts  : (optional) struct with fields:
%             .nom_idx - index into x for nominal eigenvectors (default: middle)
%             .method  - interpolation method (default: 'pchip')
%
% OUTPUTS
%   Mi    : interpolated matrix page object (n x n x length(xi))
%   info  : struct with diagnostic information
%
% ========================================================================

function [Mi, info] = interp_matrix_eig(x, M, xi, opts)

    % --- Parse optional arguments ---
    if nargin < 4; opts = struct(); end
    if ~isfield(opts,'method'); opts.method = 'pchip'; end

    [n, m, p] = size(M);
    assert(n == m, 'interp_matrix_eig: matrix must be square (n x n x p)');
    assert(numel(x) == p, 'interp_matrix_eig: length(x) must equal size(M,3)');

    if ~isfield(opts,'nom_idx')
        opts.nom_idx = round(p/2);
    end

    Nq = numel(xi);

    % =====================================================================
    % STEP 1: Eigen-decompose each page
    % =====================================================================
    all_lambda = zeros(n, p);
    all_V      = zeros(n, n, p);

    for k = 1:p
        [V, Lambda] = eig(M(:,:,k));
        lam = diag(Lambda);

        % Sort by magnitude (descending), then by angle for consistency
        [~, sort_idx] = sortrows([abs(lam), angle(lam)], [-1, -2]);
        lam = lam(sort_idx);
        V   = V(:, sort_idx);

        all_lambda(:, k) = lam;
        all_V(:,:,k)     = V;
    end

    % =====================================================================
    % STEP 2: Refine sorting — greedy nearest-neighbour matching
    %
    % After initial sorting by magnitude/angle, refine by minimising
    % jumps between adjacent operating points. This prevents accidental
    % swaps when two eigenvalues have similar magnitudes.
    % =====================================================================
    for k = 2:p
        lam_prev = all_lambda(:, k-1);
        lam_curr = all_lambda(:, k);
        V_curr   = all_V(:,:, k);

        used = false(n, 1);
        perm = zeros(n, 1);
        for ii = 1:n
            dists = abs(lam_curr - lam_prev(ii));
            dists(used) = inf;
            [~, jj] = min(dists);
            perm(ii) = jj;
            used(jj) = true;
        end

        all_lambda(:, k) = lam_curr(perm);
        all_V(:,:, k)    = V_curr(:, perm);
    end

    % =====================================================================
    % STEP 3: Interpolate eigenvalues in polar form
    % =====================================================================
    eig_mag   = abs(all_lambda);
    eig_angle = angle(all_lambda);

    % Unwrap angles to avoid 2*pi jumps
    for ii = 1:n
        eig_angle(ii, :) = unwrap(eig_angle(ii, :));
    end

    mag_interp   = zeros(n, Nq);
    angle_interp = zeros(n, Nq);

    for ii = 1:n
        mag_interp(ii, :)   = interp1(x, eig_mag(ii,:),   xi, opts.method, 'extrap');
        angle_interp(ii, :) = interp1(x, eig_angle(ii,:), xi, opts.method, 'extrap');
    end

    % Magnitude must be non-negative
    mag_interp = max(mag_interp, 0);

    % Reconstruct complex eigenvalues
    lambda_interp = mag_interp .* exp(1i * angle_interp);

    % =====================================================================
    % STEP 4: Flag instability (DO NOT clamp)
    % =====================================================================
    unstable_mask = mag_interp > 1.0;
    n_unstable    = sum(unstable_mask(:));

    if n_unstable > 0
        [~, rho_idx] = find(unstable_mask);
        rho_unstable = xi(unique(rho_idx));

        warning('interp_matrix_eig:unstable', ...
            ['%d eigenvalue point(s) have |lambda| > 1 (unstable).\n' ...
             '  rho range: [%.4f, %.4f]\n' ...
             '  Max |lambda|: %.6f\n' ...
             '  Investigate whether this is an identification artefact\n' ...
             '  or genuine physical instability.'], ...
            n_unstable, min(rho_unstable), max(rho_unstable), max(mag_interp(:)));
    end

    % =====================================================================
    % STEP 5: Reconstruct using fixed nominal eigenvectors
    % =====================================================================
    V_nom = all_V(:,:, opts.nom_idx);

    cond_V = cond(V_nom);
    if cond_V > 1e8
        warning('interp_matrix_eig:illcond', ...
            'Nominal eigenvector matrix is poorly conditioned (cond = %.1e).', cond_V);
    end

    V_nom_inv = V_nom \ eye(n);

    Mi = zeros(n, n, Nq);
    for k = 1:Nq
        lam_k = enforce_conjugate_pairs(lambda_interp(:, k));
        Mi(:,:,k) = real(V_nom * diag(lam_k) * V_nom_inv);
    end

    % =====================================================================
    % Diagnostic output
    % =====================================================================
    if nargout > 1
        % Element-wise interpolation for comparison
        M_elemwise = interp_matrix_elemwise(x, M, xi);
        eig_mag_elemwise = zeros(n, Nq);
        for k = 1:Nq
            eig_mag_elemwise(:, k) = sort(abs(eig(M_elemwise(:,:,k))), 'descend');
        end

        info.eig_mag_elemwise = eig_mag_elemwise;
        info.eig_mag_eig      = mag_interp;
        info.nom_idx          = opts.nom_idx;
        info.cond_V_nom       = cond_V;
        info.lambda_interp    = lambda_interp;
        info.lambda_tracked   = all_lambda;
        info.V_nom            = V_nom;
        info.n_unstable       = n_unstable;
        info.unstable_mask    = unstable_mask;
    end

end


%% =========================================================================
%  Enforce complex conjugate pairing
%  =========================================================================
function lam = enforce_conjugate_pairs(lam)
    n = numel(lam);
    used = false(n, 1);
    for ii = 1:n
        if used(ii); continue; end
        if abs(imag(lam(ii))) < 1e-12
            lam(ii) = real(lam(ii));
            used(ii) = true;
        else
            diffs = abs(lam - conj(lam(ii)));
            diffs(ii) = inf;
            diffs(used) = inf;
            [~, jj] = min(diffs);
            if jj > 0 && ~used(jj)
                lam_avg = (lam(ii) + conj(lam(jj))) / 2;
                lam(ii) = lam_avg;
                lam(jj) = conj(lam_avg);
                used(ii) = true;
                used(jj) = true;
            end
        end
    end
end


%% =========================================================================
%  Standard element-wise interpolation (for comparison)
%  =========================================================================
function Mi = interp_matrix_elemwise(x, M, xi)
    M = permute(M, [3, 1, 2]);
    Mi = interp1(x, M, xi);
    Mi = permute(Mi, [2, 3, 1]);
end
