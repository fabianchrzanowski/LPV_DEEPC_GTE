% ========================================================================
% Nonlinearity Demonstration via Operating-Point-Dependent Linear Models
%
% This script loads the identified discrete-time state-space models at each
% fuel-flow operating point and demonstrates that:
%   (a) The global engine map Wf → (N_L, N_H, ...) is NONLINEAR
%   (b) At each Wf, the system is well-approximated by a LINEAR model
%       (a "localised linear system")
%   (c) These local linear models (A, B, C, D) change with Wf in a
%       nonlinear fashion — eigenvalues, gains, and time constants all
%       vary as functions of the operating point.
%
% Figures produced (all dissertation-quality with LaTeX formatting):
%   1. Step response family — same input perturbation, different dynamics
%   2. Frequency response (Bode magnitude) variation across Wf
%   3. Eigenvalue migration in the complex plane
%   4. Steady-state gain map — the "global nonlinearity"
%   5. Local linearity validation — one-step prediction vs measured data
%
% ========================================================================
% Fabian Chrzanowski — University of Sheffield — 2026
% ========================================================================

clearvars; close all; clc

plot_options;

proj = currentProject;
cd(proj.RootFolder)

%% ===== Load identified models =====
id_file = fullfile(proj.RootFolder,'Params','tables','matrices_DE_wf.mat');
wf_file = fullfile(proj.RootFolder,'Params','test_files_fuel_only','test_identification','wf_values.mat');
data_folder = fullfile(proj.RootFolder,'Params','test_files_fuel_only','test_identification','AGTF30');

if ~isfile(id_file);  error('Missing ID model file: %s', id_file); end
if ~isfile(wf_file);  error('Missing Wf value file: %s', wf_file); end

models = load(id_file);
wf_data = load(wf_file);

A_dt = models.A_dt;      % (nx, nx, n_op)
B_dt = models.B_dt;      % (nx, nu, n_op)
C_dt = models.C_dt;      % (nz, nx, n_op)
D_dt = models.D_dt;      % (nz, nu, n_op)
wf_values = wf_data.Wf_values(:)';

n_op = size(A_dt, 3);
nx   = size(A_dt, 1);
nz   = size(C_dt, 1);

Ts = 0.015;  % sampling time [s]

% Output names from the identified system
output_names = models.sys_eng.OutputName;

% Offsets (equilibrium values at each operating point)
x0_all = models.offset.x;  % (nx, n_op)
z0_all = models.offset.z;  % (nz, n_op)
u0_all = models.offset.u;  % (nu, n_op)

fprintf('Loaded %d operating points across Wf = [%.2f, %.2f] pps\n', ...
    n_op, wf_values(1), wf_values(end));
fprintf('State dimension: %d, Output dimension: %d\n', nx, nz);

%% ===== Select representative operating points for clarity =====
% Use a subset for plotting (too many curves = visual clutter)
n_show = min(8, n_op);
idx_show = unique(round(linspace(1, n_op, n_show)));
n_show = numel(idx_show);

% Colormap for operating points
cmap = turbo(n_show);

% Identify N_H and N_L output indices
idx_NH = find(contains(lower(string(output_names)), 'n_h'), 1);
idx_NL = find(contains(lower(string(output_names)), 'n_l'), 1);
if isempty(idx_NH); idx_NH = 2; end
if isempty(idx_NL); idx_NL = 1; end

%% =========================================================================
%  FIGURE 1: Step Response Family
%  Same δWf step applied at each operating point → different transients
% =========================================================================
N_step = 200;                    % number of simulation steps
delta_wf = 0.05;                 % step amplitude [pps] (small perturbation)

figure('Name','Step Response Family');
tiledlayout(2, 1, 'TileSpacing', 'compact', 'Padding', 'compact');

% --- N_L subplot ---
ax1 = nexttile;
hold(ax1, 'on');

% --- N_H subplot ---
ax2 = nexttile;
hold(ax2, 'on');

valid_legs = {};
for ii = 1:n_show
    k = idx_show(ii);
    A = A_dt(:,:,k);
    B = B_dt(:,:,k);
    C = C_dt(:,:,k);
    D = D_dt(:,:,k);

    % Skip unstable models (spectral radius > 1)
    rho_k = max(abs(eig(A)));
    if rho_k > 1.01
        fprintf('  Skipping Wf=%.2f (unstable, |lambda|=%.3f)\n', wf_values(k), rho_k);
        continue
    end

    % Simulate step response (deviation form)
    x = zeros(nx, 1);
    y_dev = zeros(N_step, nz);
    for t = 1:N_step
        y_dev(t,:) = (C * x + D * delta_wf)';
        x = A * x + B * delta_wf;
        % Truncate if diverging
        if norm(x) > 1e6; y_dev(t+1:end,:) = NaN; break; end
    end

    t_vec = (0:N_step-1)' * Ts;

    plot(ax1, t_vec, y_dev(:,idx_NL), '-', 'Color', cmap(ii,:), 'LineWidth', 1.5);
    plot(ax2, t_vec, y_dev(:,idx_NH), '-', 'Color', cmap(ii,:), 'LineWidth', 1.5);

    valid_legs{end+1} = sprintf('$W_f = %.2f$ pps', wf_values(k)); %#ok<AGROW>
end

ylabel(ax1, '$\Delta N_L$ [rpm]');
title(ax1, sprintf('Step response to $\\delta W_f = %.3f$ pps at different operating points', delta_wf));
legend(ax1, valid_legs, 'Location', 'best', 'FontSize', 9);

xlabel(ax2, 'Time [s]');
ylabel(ax2, '$\Delta N_H$ [rpm]');
title(ax2, 'Same perturbation $\rightarrow$ different dynamics = \textbf{nonlinear system}');

%% =========================================================================
%  FIGURE 2: Bode-like Frequency Response Variation
%  |G(e^{jω})| at each operating point — gain changes with Wf
% =========================================================================
N_freq = 500;
omega = logspace(-1, log10(pi/Ts), N_freq);  % [rad/s]

figure('Name','Frequency Response Variation');
tiledlayout(2, 1, 'TileSpacing', 'compact', 'Padding', 'compact');

ax_bode_NL = nexttile;
hold(ax_bode_NL, 'on');

ax_bode_NH = nexttile;
hold(ax_bode_NH, 'on');

bode_legs = {};
for ii = 1:n_show
    k = idx_show(ii);
    A = A_dt(:,:,k);
    B = B_dt(:,:,k);
    C = C_dt(:,:,k);
    D = D_dt(:,:,k);

    % Skip unstable models
    if max(abs(eig(A))) > 1.01; continue; end

    mag_NL = zeros(N_freq, 1);
    mag_NH = zeros(N_freq, 1);

    for jj = 1:N_freq
        z = exp(1j * omega(jj) * Ts);
        G = C * ((z * eye(nx) - A) \ B) + D;  % transfer matrix at z
        mag_NL(jj) = min(abs(G(idx_NL, 1)), 1e6);  % cap magnitude
        mag_NH(jj) = min(abs(G(idx_NH, 1)), 1e6);
    end

    semilogx(ax_bode_NL, omega, 20*log10(mag_NL), '-', 'Color', cmap(ii,:), 'LineWidth', 1.3);
    semilogx(ax_bode_NH, omega, 20*log10(mag_NH), '-', 'Color', cmap(ii,:), 'LineWidth', 1.3);
    bode_legs{end+1} = sprintf('$W_f = %.2f$ pps', wf_values(k)); %#ok<AGROW>
end

ylabel(ax_bode_NL, '$|G_{N_L}(e^{j\omega T_s})|$ [dB]');
title(ax_bode_NL, 'Frequency response: $W_f \rightarrow N_L$');
legend(ax_bode_NL, bode_legs, 'Location', 'best', 'FontSize', 9);

xlabel(ax_bode_NH, 'Frequency $\omega$ [rad/s]');
ylabel(ax_bode_NH, '$|G_{N_H}(e^{j\omega T_s})|$ [dB]');
title(ax_bode_NH, 'Frequency response: $W_f \rightarrow N_H$');

%% =========================================================================
%  FIGURE 3: Eigenvalue Migration in the Complex Plane
%  Poles move as Wf changes — the "system DNA" is operating-point dependent
% =========================================================================
figure('Name','Eigenvalue Migration');

% Unit circle
th = linspace(0, 2*pi, 400);
plot(cos(th), sin(th), 'k--', 'LineWidth', 1.2);
hold on; axis equal; grid on;

% Plot eigenvalues coloured by Wf
cmap_full = turbo(n_op);
for k = 1:n_op
    ev = eig(A_dt(:,:,k));
    scatter(real(ev), imag(ev), 40, cmap_full(k,:), 'filled', ...
        'MarkerEdgeColor', 'k', 'LineWidth', 0.3);
end

% Add arrows showing direction of eigenvalue travel (connecting adjacent points)
% Track eigenvalues with consistent ordering
ev_tracked = nan(nx, n_op);
ev_tracked(:,1) = sort(eig(A_dt(:,:,1)));
for k = 2:n_op
    ev_k = eig(A_dt(:,:,k));
    ev_tracked(:,k) = align_poles_local(ev_tracked(:,k-1), ev_k);
end

for jj = 1:nx
    plot(real(ev_tracked(jj,:)), imag(ev_tracked(jj,:)), '-', ...
        'Color', [0.5 0.5 0.5 0.4], 'LineWidth', 0.8);
    % Arrow at midpoint
    mid = round(n_op/2);
    quiver(real(ev_tracked(jj,mid)), imag(ev_tracked(jj,mid)), ...
        real(ev_tracked(jj,mid+1) - ev_tracked(jj,mid)), ...
        imag(ev_tracked(jj,mid+1) - ev_tracked(jj,mid)), ...
        0.5, 'k', 'MaxHeadSize', 2, 'LineWidth', 1.5);
end

cb = colorbar;
colormap(turbo);
clim([wf_values(1) wf_values(end)]);
ylabel(cb, '$W_f$ [pps]');
xlabel('$\mathrm{Re}(\lambda)$');
ylabel('$\mathrm{Im}(\lambda)$');
title('Eigenvalue migration: poles move nonlinearly as $W_f$ changes');

%% =========================================================================
%  FIGURE 4: Steady-State Gain Map (The Global Nonlinearity)
%  Plot Wf vs steady-state N_H, N_L — this is the "operating map"
%  The nonlinear curvature motivates the need for LPV / DeePC
% =========================================================================
figure('Name','Steady-State Gain Map (Global Nonlinearity)');
tiledlayout(1, 2, 'TileSpacing', 'compact', 'Padding', 'compact');

% Steady-state output at each operating point (from offsets)
NL_ss = z0_all(idx_NL, :);
NH_ss = z0_all(idx_NH, :);

% Also compute DC gain at each point: y_ss = (C(I-A)^{-1}B + D) * δu
dc_gain_NL = zeros(1, n_op);
dc_gain_NH = zeros(1, n_op);
for k = 1:n_op
    % Skip if (I-A) is near-singular (pole at z=1 → infinite DC gain)
    IminA = eye(nx) - A_dt(:,:,k);
    if rcond(IminA) < 1e-10
        dc_gain_NL(k) = NaN;
        dc_gain_NH(k) = NaN;
    else
        G_dc = C_dt(:,:,k) * (IminA \ B_dt(:,:,k)) + D_dt(:,:,k);
        dc_gain_NL(k) = G_dc(idx_NL, 1);
        dc_gain_NH(k) = G_dc(idx_NH, 1);
    end
end

% --- Panel A: Operating map ---
ax_map = nexttile;
yyaxis(ax_map, 'left');
plot(ax_map, wf_values, NL_ss, '-o', 'MarkerSize', 4, 'LineWidth', 1.5);
ylabel(ax_map, '$N_L$ [rpm]');
yyaxis(ax_map, 'right');
plot(ax_map, wf_values, NH_ss, '-s', 'MarkerSize', 4, 'LineWidth', 1.5);
ylabel(ax_map, '$N_H$ [rpm]');
xlabel(ax_map, '$W_f$ [pps]');
title(ax_map, 'Steady-state operating map');

% Overlay tangent lines at a few points to show local linearity
yyaxis(ax_map, 'left');
hold(ax_map, 'on');
for ii = [1, round(n_show/2), n_show]
    k = idx_show(ii);
    if k > 1 && k < n_op
        slope = (NL_ss(k+1) - NL_ss(k-1)) / (wf_values(k+1) - wf_values(k-1));
        wf_local = linspace(wf_values(max(1,k-3)), wf_values(min(n_op,k+3)), 20);
        NL_tangent = NL_ss(k) + slope * (wf_local - wf_values(k));
        plot(ax_map, wf_local, NL_tangent, '--', 'Color', cmap(ii,:), ...
            'LineWidth', 1.0);
    end
end

% --- Panel B: DC gain variation ---
ax_gain = nexttile;
yyaxis(ax_gain, 'left');
plot(ax_gain, wf_values, dc_gain_NL, '-o', 'MarkerSize', 4, 'LineWidth', 1.5);
ylabel(ax_gain, '$G_{DC}^{N_L}$ [rpm/pps]');
yyaxis(ax_gain, 'right');
plot(ax_gain, dc_gain_NH, '-s', 'MarkerSize', 4, 'LineWidth', 1.5);
ylabel(ax_gain, '$G_{DC}^{N_H}$ [rpm/pps]');
xlabel(ax_gain, '$W_f$ [pps]');
title(ax_gain, 'DC gain $G(1) = C(I-A)^{-1}B + D$ varies with $W_f$');

%% =========================================================================
%  FIGURE 5: Local Linearity Validation
%  Show that WITHIN a neighbourhood of each Wf, the linear model is a good
%  approximation (small one-step prediction error), but ACROSS operating
%  points it fails (large error = the system is globally nonlinear).
% =========================================================================
if isfolder(data_folder)
    filelist = dir(fullfile(data_folder, '*.mat'));
    [~, idx_sort] = sort({filelist.name});
    filelist = filelist(idx_sort);

    if numel(filelist) >= n_op
        tstart = 19.5;
        ridge_lambda = 1e-6;

        % Pick 3 representative operating points: low, mid, high
        idx_demo = [1, round(n_op/2), n_op];

        figure('Name','Local Linearity Validation');
        tiledlayout(numel(idx_demo), 3, 'TileSpacing', 'compact', 'Padding', 'compact');

        for ii = 1:numel(idx_demo)
            k = idx_demo(ii);

            % Load measured data at this operating point
            sample = load(fullfile(filelist(k).folder, filelist(k).name));
            out_dyn = get_out_dyn_local(sample);
            [psi_meas, u_meas] = extract_cloud_local(out_dyn, Ts, tstart);

            N_data = size(psi_meas, 1);

            % --- Column 1: model at own operating point (good fit) ---
            A_own = A_dt(:,:,k);
            B_own = B_dt(:,:,k);
            psi0_own = [models.offset.x(:,k); models.offset.y(:,k)].';
            u0_own = models.offset.u(:,k);
            psi_pred_own = one_step_predict_local(A_own, B_own, psi_meas, u_meas, psi0_own, u0_own);

            % --- Column 2: model from a DIFFERENT operating point (poor fit) ---
            k_wrong = idx_demo(mod(ii, numel(idx_demo)) + 1);  % use next demo point
            A_wrong = A_dt(:,:,k_wrong);
            B_wrong = B_dt(:,:,k_wrong);
            psi0_wrong = [models.offset.x(:,k_wrong); models.offset.y(:,k_wrong)].';
            u0_wrong = models.offset.u(:,k_wrong);
            psi_pred_wrong = one_step_predict_local(A_wrong, B_wrong, psi_meas, u_meas, psi0_wrong, u0_wrong);

            t_data = (0:N_data-1)' * Ts;

            % Column 1: N_H measured vs own-model prediction
            ax = nexttile;
            plot(ax, t_data, psi_meas(:,2), 'k-', 'LineWidth', 1.0);
            hold(ax, 'on');
            plot(ax, t_data, psi_pred_own(:,2), '-', 'Color', [0 0.45 0.74], 'LineWidth', 1.0);
            ylabel(ax, '$N_H$ [rpm]');
            if ii == 1; title(ax, 'Local model (correct $W_f$)'); end
            if ii == numel(idx_demo); xlabel(ax, 'Time [s]'); end
            text(ax, 0.02, 0.95, sprintf('$W_f = %.2f$', wf_values(k)), ...
                'Units', 'normalized', 'VerticalAlignment', 'top', 'FontSize', 10);

            % Column 2: N_H measured vs wrong-model prediction
            ax = nexttile;
            plot(ax, t_data, psi_meas(:,2), 'k-', 'LineWidth', 1.0);
            hold(ax, 'on');
            plot(ax, t_data, psi_pred_wrong(:,2), '-', 'Color', [0.85 0.33 0.1], 'LineWidth', 1.0);
            ylabel(ax, '$N_H$ [rpm]');
            if ii == 1; title(ax, sprintf('Wrong model ($W_f = %.2f$)', wf_values(k_wrong))); end
            if ii == numel(idx_demo); xlabel(ax, 'Time [s]'); end

            % Column 3: prediction error comparison
            err_own   = abs(psi_meas(:,2) - psi_pred_own(:,2));
            err_wrong = abs(psi_meas(:,2) - psi_pred_wrong(:,2));

            ax = nexttile;
            plot(ax, t_data, err_own, '-', 'Color', [0 0.45 0.74], 'LineWidth', 1.0);
            hold(ax, 'on');
            plot(ax, t_data, err_wrong, '-', 'Color', [0.85 0.33 0.1], 'LineWidth', 1.0);
            ylabel(ax, '$|e|$ [rpm]');
            if ii == 1
                title(ax, 'Prediction error');
                legend(ax, 'Correct $W_f$', 'Wrong $W_f$', 'Location', 'best', 'FontSize', 9);
            end
            if ii == numel(idx_demo); xlabel(ax, 'Time [s]'); end
        end

        sgtitle('Local linearity: correct model fits well $\rightarrow$ wrong model fails $\rightarrow$ system is globally nonlinear');
    else
        fprintf('Not enough data files for Figure 5 (have %d, need %d).\n', ...
            numel(filelist), n_op);
    end
else
    fprintf('Data folder not found — skipping Figure 5.\n');
end

%% =========================================================================
%  FIGURE 6: Summary — A matrix entry variation (nonlinear in ρ)
%  Show how individual entries of A change across operating points,
%  emphasising the smooth but nonlinear variation.
% =========================================================================
figure('Name','A Matrix Entry Variation across Wf');
nr = nx; nc = nx;
tiledlayout(nr, nc, 'TileSpacing', 'compact', 'Padding', 'compact');

for r = 1:nr
    for c = 1:nc
        nexttile;
        entries = squeeze(A_dt(r,c,:));
        plot(wf_values, entries, '-o', 'MarkerSize', 3, 'LineWidth', 1.3, ...
            'Color', [0 0.45 0.74]);
        title(sprintf('$a_{%d%d}$', r, c));
        if r == nr; xlabel('$W_f$ [pps]'); end
    end
end
sgtitle('$A$ matrix entries vary \textbf{nonlinearly} with operating point $W_f$');

fprintf('\n===== Summary =====\n');
fprintf('All 6 figures demonstrate that:\n');
fprintf('  - The engine is a NONLINEAR system (Figs 1,3,4,6)\n');
fprintf('  - At each Wf, a LOCAL LINEAR model is accurate (Fig 5)\n');
fprintf('  - Dynamics (eigenvalues, gains, time constants) change with Wf (Figs 2,3,4)\n');
fprintf('  - This motivates LPV / gain-scheduled control (different A,B,C,D per region)\n');

%% ===== Helper functions =====
function ev2_aligned = align_poles_local(ev_ref, ev2)
    n = numel(ev_ref);
    ev2_aligned = zeros(size(ev2));
    used = false(n,1);
    for i = 1:n
        d = abs(ev2 - ev_ref(i));
        d(used) = inf;
        [~,idx] = min(d);
        ev2_aligned(i) = ev2(idx);
        used(idx) = true;
    end
end

function out_dyn = get_out_dyn_local(data)
    if isstruct(data) && isfield(data,'out_Dyn')
        out_dyn = data.out_Dyn;
        return
    end
    error('Could not find out_Dyn in loaded data file.');
end

function [psi, u] = extract_cloud_local(out_dyn, Ts, tstart)
    t_raw = signal_time_local(out_dyn.eng.S0.W);
    t = (tstart:Ts:t_raw(end)).';

    x1 = signal_data_local(out_dyn.eng.Shaft.N_Fan);
    x2 = signal_data_local(out_dyn.eng.Shaft.N_HPC);
    y1 = signal_data_local(out_dyn.eng.S21.Pt);
    y2 = signal_data_local(out_dyn.eng.S25.Pt);
    u_raw = signal_data_local(out_dyn.cntrl.Wfact);

    x1_i = interp1(t_raw, x1, t, 'linear', 'extrap');
    x2_i = interp1(t_raw, x2, t, 'linear', 'extrap');
    y1_i = interp1(t_raw, y1, t, 'linear', 'extrap');
    y2_i = interp1(t_raw, y2, t, 'linear', 'extrap');
    u_i  = interp1(t_raw, u_raw, t, 'linear', 'extrap');

    psi = [x1_i x2_i y1_i y2_i];
    u = u_i;
end

function psi_hat = one_step_predict_local(A, B, psi_meas, u_meas, psi0, u0)
    N = size(psi_meas, 1);
    nx_local = size(psi_meas, 2);
    psi_hat = nan(N, nx_local);
    psi_hat(1,:) = psi_meas(1,:);

    for k = 1:N-1
        xk = psi_meas(k,:).' - psi0(:);
        uk = u_meas(k) - u0;
        xkp1 = psi0(:) + A * xk + B * uk;
        psi_hat(k+1,:) = xkp1.';
    end
end

function d = signal_data_local(sig)
    if isa(sig,'timeseries')
        d = sig.Data(:);
    elseif isstruct(sig) && isfield(sig,'Data')
        d = sig.Data(:);
    else
        d = sig(:);
    end
    d = double(d);
end

function t = signal_time_local(sig)
    if isa(sig,'timeseries')
        t = sig.Time(:);
    elseif isstruct(sig) && isfield(sig,'Time')
        t = sig.Time(:);
    else
        error('Signal time vector not found.');
    end
    t = double(t);
end
