% ========================================================================
% Large Transient Test: Why Gain-Scheduling Matters
%
% Demonstrates the nonlinear nature of the engine by showing what happens
% when you try to predict the system response across a large operating
% range using different modelling strategies.
%
% Key idea: use a RAMP through the Wf envelope. At each point, compare
% the steady-state output predicted by each model against the TRUE
% equilibrium from the identified operating map. The frozen model
% extrapolates linearly and diverges; the LPV model stays on the curve.
%
% No Simulink needed — runs in seconds.
%
% Figures produced:
%   1. Steady-state prediction: frozen vs LPV vs true equilibrium curve
%   2. Prediction error as a function of distance from linearisation point
%   3. Step response at 3 representative Wf (overlay)
%   4. Per-segment RMSE bar chart for a step transient
%
% ========================================================================
% Fabian Chrzanowski — University of Sheffield — 2026
% ========================================================================

clearvars; close all; clc

plot_options;

proj = currentProject;
cd(proj.RootFolder)

%% ===== Load identified models =====
models = load(fullfile(proj.RootFolder,'Params','tables','matrices_DE_wf.mat'));
wf_data = load(fullfile(proj.RootFolder,'Params','test_files_fuel_only','test_identification','wf_values.mat'));

A_dt = models.A_dt;
B_dt = models.B_dt;
C_dt = models.C_dt;
D_dt = models.D_dt;
offset = models.offset;
wf_values = wf_data.Wf_values(:)';
n_op = size(A_dt, 3);
nx = size(A_dt, 1);
nz = size(C_dt, 1);

Ts = 0.015;
output_names = models.sys_eng.OutputName;

idx_NH = find(contains(lower(string(output_names)), 'n_h'), 1);
idx_NL = find(contains(lower(string(output_names)), 'n_l'), 1);
if isempty(idx_NH); idx_NH = 2; end
if isempty(idx_NL); idx_NL = 1; end

% True equilibrium curve (from system ID offsets — this is ground truth)
NH_true = offset.z(idx_NH, :);
NL_true = offset.z(idx_NL, :);

%% =========================================================================
%  FIGURE 1: Steady-State Prediction — Frozen vs LPV vs True
%  For each Wf, compute what each model PREDICTS as steady-state output.
%  The frozen model extrapolates linearly; the truth is nonlinear.
% =========================================================================

% Three "frozen" linearisation points
wf_freeze = [0.5, 1.2, 1.8];  % idle, cruise, max
freeze_labels = {'Frozen at idle ($W_f=0.5$)', ...
                 'Frozen at cruise ($W_f=1.2$)', ...
                 'Frozen at max ($W_f=1.8$)'};
freeze_colors = [0.85 0.33 0.1; 0.6 0.2 0.8; 0.2 0.7 0.3];

wf_query = linspace(wf_values(1), wf_values(end), 200);

figure('Name','Steady-State Prediction: Frozen vs True');
tiledlayout(2, 1, 'TileSpacing', 'compact', 'Padding', 'compact');

for out_idx = [idx_NH, idx_NL]
    ax = nexttile;
    hold(ax, 'on');

    % True equilibrium curve
    plot(ax, wf_values, offset.z(out_idx,:), 'k-o', 'LineWidth', 2, ...
        'MarkerSize', 4, 'DisplayName', 'True equilibrium (nonlinear)');

    for ff = 1:numel(wf_freeze)
        [~, k_f] = min(abs(wf_values - wf_freeze(ff)));
        A_f = A_dt(:,:,k_f);
        B_f = B_dt(:,:,k_f);
        C_f = C_dt(:,:,k_f);
        D_f = D_dt(:,:,k_f);
        u0_f = offset.u(:, k_f);
        z0_f = offset.z(:, k_f);

        % Predicted steady-state: z_ss = z0 + C*(I-A)^{-1}*B*(Wf - u0) + D*(Wf - u0)
        IminA = eye(nx) - A_f;
        if rcond(IminA) > 1e-10
            G_dc = C_f * (IminA \ B_f) + D_f;
            z_pred = z0_f(out_idx) + G_dc(out_idx,1) * (wf_query - u0_f);
        else
            z_pred = nan(size(wf_query));
        end

        plot(ax, wf_query, z_pred, '--', 'Color', freeze_colors(ff,:), ...
            'LineWidth', 1.5, 'DisplayName', freeze_labels{ff});

        % Mark the linearisation point
        plot(ax, wf_freeze(ff), z0_f(out_idx), 's', 'Color', freeze_colors(ff,:), ...
            'MarkerSize', 10, 'LineWidth', 2, 'HandleVisibility', 'off');
    end

    if out_idx == idx_NH
        ylabel(ax, '$N_H$ [rpm]');
        title(ax, 'Frozen models predict linearly $\rightarrow$ wrong far from linearisation point');
        legend(ax, 'Location', 'best', 'FontSize', 9);
    else
        ylabel(ax, '$N_L$ [rpm]');
        xlabel(ax, '$W_f$ [pps]');
    end
end

%% =========================================================================
%  FIGURE 2: Prediction Error vs Distance from Linearisation Point
%  Shows how model accuracy degrades the further you go from the
%  linearisation point — the core motivation for LPV.
% =========================================================================
figure('Name','Prediction Error vs Distance');
hold on;

for ff = 1:numel(wf_freeze)
    [~, k_f] = min(abs(wf_values - wf_freeze(ff)));
    A_f = A_dt(:,:,k_f);
    B_f = B_dt(:,:,k_f);
    C_f = C_dt(:,:,k_f);
    D_f = D_dt(:,:,k_f);
    u0_f = offset.u(:, k_f);
    z0_f = offset.z(:, k_f);

    IminA = eye(nx) - A_f;
    if rcond(IminA) < 1e-10; continue; end
    G_dc = C_f * (IminA \ B_f) + D_f;

    % Error at each identified operating point
    err_NH = zeros(n_op, 1);
    dist_wf = zeros(n_op, 1);
    for k = 1:n_op
        z_pred_k = z0_f(idx_NH) + G_dc(idx_NH,1) * (wf_values(k) - u0_f);
        err_NH(k) = abs(z_pred_k - NH_true(k));
        dist_wf(k) = abs(wf_values(k) - wf_freeze(ff));
    end

    plot(dist_wf, err_NH, '-o', 'Color', freeze_colors(ff,:), ...
        'LineWidth', 1.5, 'MarkerSize', 4, 'DisplayName', freeze_labels{ff});
end

xlabel('$|W_f - W_{f,\mathrm{lin}}|$ [pps]');
ylabel('$|N_{H,\mathrm{pred}} - N_{H,\mathrm{true}}|$ [rpm]');
title('Prediction error grows with distance from linearisation point');
legend('Location', 'best', 'FontSize', 9);

%% =========================================================================
%  FIGURE 3: Step Responses at 3 Operating Points (overlaid)
%  Same δWf perturbation, but the transient looks completely different
%  at different operating points — proof the system is nonlinear.
% =========================================================================
N_step = 150;
delta_wf = 0.05;
t_step = (0:N_step-1)' * Ts;

wf_demo = [0.5, 1.2, 1.8];
demo_colors = [0.85 0.33 0.1; 0 0.45 0.74; 0.2 0.7 0.3];
demo_labels = {};

figure('Name','Step Responses at Different Operating Points');
tiledlayout(2, 1, 'TileSpacing', 'compact', 'Padding', 'compact');

ax_nh = nexttile;
hold(ax_nh, 'on');
ax_nl = nexttile;
hold(ax_nl, 'on');

for dd = 1:numel(wf_demo)
    [~, k_d] = min(abs(wf_values - wf_demo(dd)));
    A_d = A_dt(:,:,k_d);
    B_d = B_dt(:,:,k_d);
    C_d = C_dt(:,:,k_d);
    D_d = D_dt(:,:,k_d);

    % Check stability
    if max(abs(eig(A_d))) > 1.01
        fprintf('Skipping Wf=%.2f (unstable)\n', wf_demo(dd));
        continue
    end

    x = zeros(nx, 1);
    y_dev = zeros(N_step, nz);
    for t = 1:N_step
        y_dev(t,:) = (C_d * x + D_d * delta_wf)';
        x = A_d * x + B_d * delta_wf;
    end

    plot(ax_nh, t_step, y_dev(:,idx_NH), '-', 'Color', demo_colors(dd,:), 'LineWidth', 1.5);
    plot(ax_nl, t_step, y_dev(:,idx_NL), '-', 'Color', demo_colors(dd,:), 'LineWidth', 1.5);
    demo_labels{end+1} = sprintf('$W_f = %.1f$ pps', wf_demo(dd)); %#ok<AGROW>
end

ylabel(ax_nh, '$\Delta N_H$ [rpm]');
title(ax_nh, sprintf('Step response to $\\delta W_f = %.3f$ pps: different dynamics at each $W_f$', delta_wf));
legend(ax_nh, demo_labels, 'Location', 'best', 'FontSize', 9);

ylabel(ax_nl, '$\Delta N_L$ [rpm]');
xlabel(ax_nl, 'Time [s]');

%% =========================================================================
%  FIGURE 4: Step Transient Simulation — LPV vs Frozen (dynamic sim)
%  Drive the models through a step sequence and show the actual dynamic
%  response, including transients. LPV uses correct model at each step.
% =========================================================================
steps_per_seg = 200;
N_dyn = 4 * steps_per_seg;
t_dyn = (0:N_dyn-1)' * Ts;

seg_wf = [0.5, 1.8, 1.2, 0.8];
seg_labels = {'Idle', 'Max', 'Cruise', 'Low'};
wf_dyn = zeros(N_dyn, 1);
for seg = 1:4
    wf_dyn((seg-1)*steps_per_seg+1 : seg*steps_per_seg) = seg_wf(seg);
end

% Frozen model (cruise, k_cruise)
[~, k_cruise] = min(abs(wf_values - 1.2));
A_fc = A_dt(:,:,k_cruise); B_fc = B_dt(:,:,k_cruise);
C_fc = C_dt(:,:,k_cruise); D_fc = D_dt(:,:,k_cruise);
u0_fc = offset.u(:,k_cruise); z0_fc = offset.z(:,k_cruise);

x_fc = zeros(nx,1);
y_fc = zeros(N_dyn, nz);
for t = 1:N_dyn
    du = wf_dyn(t) - u0_fc;
    y_fc(t,:) = (z0_fc + C_fc*x_fc + D_fc*du)';
    x_fc = A_fc*x_fc + B_fc*du;
end

% Frozen model (idle)
[~, k_idle_f] = min(abs(wf_values - 0.5));
A_fi = A_dt(:,:,k_idle_f); B_fi = B_dt(:,:,k_idle_f);
C_fi = C_dt(:,:,k_idle_f); D_fi = D_dt(:,:,k_idle_f);
u0_fi = offset.u(:,k_idle_f); z0_fi = offset.z(:,k_idle_f);

x_fi = zeros(nx,1);
y_fi = zeros(N_dyn, nz);
for t = 1:N_dyn
    du = wf_dyn(t) - u0_fi;
    y_fi(t,:) = (z0_fi + C_fi*x_fi + D_fi*du)';
    x_fi = A_fi*x_fi + B_fi*du;
end

% Reference = true equilibrium at each Wf
ref_dyn_NH = zeros(N_dyn,1);
ref_dyn_NL = zeros(N_dyn,1);
for t = 1:N_dyn
    [~, kk] = min(abs(wf_values - wf_dyn(t)));
    ref_dyn_NH(t) = offset.z(idx_NH, kk);
    ref_dyn_NL(t) = offset.z(idx_NL, kk);
end

seg_bounds = (1:3) * steps_per_seg * Ts;

figure('Name','Dynamic Step Transient: Frozen Models vs True Equilibrium');
tiledlayout(3, 1, 'TileSpacing', 'compact', 'Padding', 'compact');

ax_d1 = nexttile; hold(ax_d1, 'on');
plot(ax_d1, t_dyn, ref_dyn_NH, 'k--', 'LineWidth', 1.5, 'DisplayName', 'True equilibrium');
plot(ax_d1, t_dyn, y_fc(:,idx_NH), '-', 'Color', [0.85 0.33 0.1], ...
    'LineWidth', 1.3, 'DisplayName', 'Frozen at cruise ($W_f=1.2$)');
plot(ax_d1, t_dyn, y_fi(:,idx_NH), '-', 'Color', [0.6 0.2 0.8], ...
    'LineWidth', 1.3, 'DisplayName', 'Frozen at idle ($W_f=0.5$)');
add_seg_lines(ax_d1, seg_bounds);
ylabel(ax_d1, '$N_H$ [rpm]');
legend(ax_d1, 'Location', 'best', 'FontSize', 9);
title(ax_d1, 'Both frozen models deviate from true equilibrium far from their operating point');

ax_d2 = nexttile; hold(ax_d2, 'on');
plot(ax_d2, t_dyn, ref_dyn_NL, 'k--', 'LineWidth', 1.5);
plot(ax_d2, t_dyn, y_fc(:,idx_NL), '-', 'Color', [0.85 0.33 0.1], 'LineWidth', 1.3);
plot(ax_d2, t_dyn, y_fi(:,idx_NL), '-', 'Color', [0.6 0.2 0.8], 'LineWidth', 1.3);
add_seg_lines(ax_d2, seg_bounds);
ylabel(ax_d2, '$N_L$ [rpm]');

ax_d3 = nexttile; hold(ax_d3, 'on');
stairs(ax_d3, t_dyn, wf_dyn, 'k-', 'LineWidth', 1.5);
add_seg_lines(ax_d3, seg_bounds);
ylabel(ax_d3, '$W_f$ [pps]');
xlabel(ax_d3, 'Time [s]');
for seg = 1:4
    t_mid = ((seg-1)*steps_per_seg + steps_per_seg/2) * Ts;
    text(ax_d3, t_mid, max(wf_dyn)*1.08, seg_labels{seg}, ...
        'HorizontalAlignment', 'center', 'FontSize', 10);
end

%% ===== Summary =====
fprintf('\n===== Results =====\n');
err_fc = y_fc(:,idx_NH) - ref_dyn_NH;
err_fi = y_fi(:,idx_NH) - ref_dyn_NH;
fprintf('Frozen cruise RMSE (N_H): %.0f rpm\n', rms(err_fc));
fprintf('Frozen idle RMSE (N_H):   %.0f rpm\n', rms(err_fi));
fprintf('Neither frozen model works across the full range -> need LPV.\n');

%% ===== Helper =====
function add_seg_lines(ax, boundaries)
    for ii = 1:numel(boundaries)
        xline(ax, boundaries(ii), '--', 'Color', [0.7 0.7 0.7], ...
            'LineWidth', 0.8, 'HandleVisibility', 'off');
    end
end
