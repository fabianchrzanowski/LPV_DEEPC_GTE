% ========================================================================
% LPV-DeePC Dissertation Plotter: Control Diagnostics & Mathematical Plots
%
% This script loads data from saved .mat closed-loop simulation files 
% and generates 4 figures:
%   1. Cost Function Comparison
%   2. Singular Value Decomposition of the Hankel Matrix
%   3. Slack Variable (\sigma) Norm over Time
%   4. g-Vector Weight Distribution (Histogram & Stem)
% ========================================================================
clearvars; close all; clc

try proj = currentProject; catch; proj = openProject('GTE.prj'); end
cd(proj.RootFolder)

% Load formatting defaults if available
if isfile(fullfile(proj.RootFolder, 'Library', 'Plotting', 'plot_options.m'))
    run('plot_options.m');
end

%% ===== File Selection Logic =====
results_dir = fullfile(pwd, 'Results');
mat_files = dir(fullfile(results_dir, '*.mat'));
if isempty(mat_files)
    error('No .mat files found in %s', results_dir);
end

disp('Available Result Files:');
for i = 1:length(mat_files)
    fprintf('[%d] %s\n', i, mat_files(i).name);
end

fprintf('\nSelect exactly ONE LPV-DeePC file for diagnostics (e.g., 1): ');
selection_deepc = input('');
if isempty(selection_deepc) || selection_deepc < 1 || selection_deepc > length(mat_files)
    error('Invalid selection.');
end
file_deepc = fullfile(results_dir, mat_files(selection_deepc).name);
data_deepc = load(file_deepc);

fprintf('Select exactly ONE LPV-MPC file for cost comparison (or press ENTER to skip): ');
selection_mpc = input('');
if ~isempty(selection_mpc) && selection_mpc >= 1 && selection_mpc <= length(mat_files)
    file_mpc = fullfile(results_dir, mat_files(selection_mpc).name);
    data_mpc = load(file_mpc);
    has_mpc = true;
else
    has_mpc = false;
end

%% ===== Check for Required Diagnostic Data =====
if ~isfield(data_deepc, 'cost_track_hist') || ~isfield(data_deepc, 'g_hist_matrix') || ~isfield(data_deepc, 'sigma_hist')
    error('The selected DeePC file does not contain diagnostic variables (g_hist, sigma_hist, etc.). Please re-run the DeePC controller simulation with the updated scripts to generate them!');
end

Ts = 0.015;
T_ini = 20;

cost_total_deepc = data_deepc.cost_hist(T_ini+1:end);
cost_track_deepc = data_deepc.cost_track_hist(T_ini+1:end);
cost_du_deepc    = data_deepc.cost_du_hist(T_ini+1:end);
cost_g_deepc     = data_deepc.cost_g_hist(T_ini+1:end);
cost_sigma_deepc = data_deepc.cost_sigma_hist(T_ini+1:end);

g_hist_matrix   = data_deepc.g_hist_matrix(T_ini+1:end, :);
sigma_hist      = data_deepc.sigma_hist(T_ini+1:end, :);

N_sim = length(cost_total_deepc);
time_sim = (0:N_sim-1)'*Ts;

sigma_norm_hist = sqrt(sum(sigma_hist.^2, 2)); % 2-norm of rows

%% ===== Load Hankel Matrices for SVD Plot =====
fprintf('Loading training data for Hankel SVD...\n');
training_folder = fullfile(proj.RootFolder,'Params','test_files_fuel_only','test_identification','AGTF30');
wf_data = load(fullfile(proj.RootFolder,'Params','test_files_fuel_only','test_identification','wf_values.mat'));
wf_values = wf_data.Wf_values(:)';

% L and N_sys from theory
L = 30; % N_ini + N_pred
n_sys = 4;
rank_limit = L + n_sys;

filelist = dir(fullfile(training_folder,'*.mat'));
[~,idx_sort] = sort({filelist.name});
filelist = filelist(idx_sort);
n_op = min(numel(filelist), numel(wf_values));

hankel_U = cell(n_op, 1);
hankel_Y = cell(n_op, 1);
hankel_wf = wf_values(1:n_op);

for k = 1:n_op
    data_train = load(fullfile(filelist(k).folder, filelist(k).name));
    out_dyn = fetch_out_dyn(data_train);
    u = extract_input_trace(out_dyn);
    y = measure_tracking_series(out_dyn, data_deepc.track_name);
    n_common = min(numel(u), numel(y));
    
    u_mean = mean(u(1:n_common));
    y_mean = mean(y(1:n_common));
    du = u(1:n_common) - u_mean;
    dy = y(1:n_common) - y_mean;
    
    hankel_U{k} = block_hankel(du(:), L);
    hankel_Y{k} = block_hankel(dy(:), L);
end

% Remove empty entries
valid = ~cellfun(@isempty, hankel_U);
hankel_U = hankel_U(valid);
hankel_Y = hankel_Y(valid);
hankel_wf = hankel_wf(valid);
n_valid = sum(valid);


%% ===== Figure 1: Cost Function comparison =====
fig1 = figure('Name', 'Cost Function Comparison');
t1 = tiledlayout(fig1, 2, 1, 'TileSpacing', 'compact', 'Padding', 'compact');

ax1a = nexttile(t1);
if has_mpc
    cost_total_mpc = data_mpc.cost_hist(T_ini+1:end);
    % Pad or truncate MPC to match length
    if length(cost_total_mpc) > N_sim
        cost_total_mpc = cost_total_mpc(1:N_sim);
    elseif length(cost_total_mpc) < N_sim
        cost_total_mpc = [cost_total_mpc; nan(N_sim - length(cost_total_mpc), 1)];
    end
    
    % Shift MPC forward by 1 step to match DeePC alignment
    cost_total_mpc = [cost_total_mpc(2:end); cost_total_mpc(end)];
    
    if isfield(data_mpc, 'cost_track_hist')
        cost_track_mpc = data_mpc.cost_track_hist(T_ini+1:end);
        cost_du_mpc = data_mpc.cost_du_hist(T_ini+1:end);
        cost_track_mpc = [cost_track_mpc(2:end); cost_track_mpc(end)];
        cost_du_mpc = [cost_du_mpc(2:end); cost_du_mpc(end)];
    else
        cost_track_mpc = [];
    end
    plot(ax1a, time_sim, cost_total_mpc, '-', 'Color', [0.85 0.35 0.10], 'LineWidth', 1.6, 'DisplayName', 'LPV-MPC Total Cost');
end
hold(ax1a, 'on');
plot(ax1a, time_sim, cost_total_deepc, '-', 'Color', [0.15 0.65 0.25], 'LineWidth', 1.8, 'DisplayName', 'LPV-DeePC Total Cost');
set(ax1a, 'TickLabelInterpreter', 'latex', 'FontSize', 11);
ylabel(ax1a, 'Total Cost $J(t)$', 'Interpreter', 'latex', 'FontSize', 12);
title(ax1a, '\textbf{LPV-DeePC vs LPV-MPC Stage Cost Comparison}', 'Interpreter', 'latex', 'FontSize', 13);
legend(ax1a, 'Location', 'best', 'Interpreter', 'latex');
grid(ax1a, 'on');

% Breakdown of DeePC cost
ax1b = nexttile(t1);
hold(ax1b, 'on');
plot(ax1b, time_sim, cost_track_deepc, '-', 'Color', [0.15 0.65 0.25], 'LineWidth', 1.2, 'DisplayName', 'DeePC Tracking Cost ($J_y$)');
plot(ax1b, time_sim, cost_du_deepc, '-', 'Color', [0.15 0.65 0.25], 'LineWidth', 1.2, 'LineStyle', ':', 'DisplayName', 'DeePC Control Slew Cost ($J_{\Delta u}$)');
plot(ax1b, time_sim, cost_g_deepc, '-', 'LineWidth', 1.2, 'DisplayName', 'DeePC Regularization ($J_g$)');
plot(ax1b, time_sim, cost_sigma_deepc, '-', 'LineWidth', 1.2, 'DisplayName', 'DeePC Slack Cost ($J_\sigma$)');

if has_mpc && ~isempty(cost_track_mpc)
    plot(ax1b, time_sim, cost_track_mpc, '-', 'Color', [0.85 0.35 0.10], 'LineWidth', 1.2, 'DisplayName', 'MPC Tracking Cost ($J_y$)');
    plot(ax1b, time_sim, cost_du_mpc, '-', 'Color', [0.85 0.35 0.10], 'LineWidth', 1.2, 'LineStyle', ':', 'DisplayName', 'MPC Control Slew Cost ($J_{\Delta u}$)');
end
set(ax1b, 'TickLabelInterpreter', 'latex', 'FontSize', 11);
xlabel(ax1b, 'Time [s]', 'Interpreter', 'latex', 'FontSize', 12);
ylabel(ax1b, 'Cost Components', 'Interpreter', 'latex', 'FontSize', 12);
title(ax1b, '\textbf{LPV-DeePC Stage Cost Components Breakdown}', 'Interpreter', 'latex', 'FontSize', 12);
legend(ax1b, 'Location', 'best', 'Interpreter', 'latex', 'FontSize', 9);
grid(ax1b, 'on');

%% ===== Figure 2: Hankel Matrix Singular Values =====
fig2 = figure('Name', 'Hankel Singular Values');
ax2 = axes(fig2);
hold(ax2, 'on');
set(ax2, 'YScale', 'log');

% Set up colormap based on Wf range
colormap_name = 'turbo';
colormap(ax2, colormap_name);
wf_min = min(hankel_wf);
wf_max = max(hankel_wf);
colors = turbo(n_valid);

% Sort operating points by Wf to ensure smooth color transition
[sorted_wf, sort_idx] = sort(hankel_wf);

for idx = 1:n_valid
    k = sort_idx(idx);
    Hu_k = hankel_U{k};
    Hy_k = hankel_Y{k};
    H_k = [Hu_k; Hy_k];
    
    s_vals = svd(H_k);
    c = colors(idx, :);
    plot(ax2, s_vals, '-', 'Color', c, 'LineWidth', 1.2, 'HandleVisibility', 'off');
end

% Overlay the midpoint operating point with black markers for clear reference
% mid_idx = sort_idx(round(n_valid / 2));
% H_mid = [hankel_U{mid_idx}; hankel_Y{mid_idx}];
% s_vals_mid = svd(H_mid);
% semilogy(ax2, s_vals_mid, 'o--', 'Color', 'k', 'LineWidth', 1.8, ...
%     'MarkerFaceColor', 'k', 'MarkerSize', 4, ...
%     'DisplayName', sprintf('Midpoint OP ($W_f = %.2f$ pps)', hankel_wf(mid_idx)));

% Plot numerical rank cutoff
xline(ax2, rank_limit, 'r--', 'LineWidth', 2.0, ...
    'DisplayName', sprintf('Willems Rank Limit ($L+n_x = %d$)', rank_limit));

% Add colorbar to explain color range
cb = colorbar(ax2);
if exist('clim', 'file')
    clim(ax2, [wf_min, wf_max]);
else
    caxis(ax2, [wf_min, wf_max]);
end
ylabel(cb, 'Operating Point $W_f$ [pps]', 'Interpreter', 'latex', 'FontSize', 11);

set(ax2, 'TickLabelInterpreter', 'latex', 'FontSize', 11);
xlabel(ax2, 'Singular Value Index', 'Interpreter', 'latex', 'FontSize', 12);
ylabel(ax2, 'Singular Value $\sigma_i$ (Log Scale)', 'Interpreter', 'latex', 'FontSize', 12);
title(ax2, '\textbf{Hankel Singular Value Spectrum Across Operating Envelope}', 'Interpreter', 'latex', 'FontSize', 13);
legend(ax2, 'Location', 'northeast', 'Interpreter', 'latex');
grid(ax2, 'on');

%% ===== Figure 3: Slack Variable (\sigma) Norm over Time =====
fig3 = figure('Name', 'Slack Variable History');
ax3 = axes(fig3);
hold(ax3, 'on');

plot(ax3, time_sim, sigma_norm_hist, '-', 'Color', [0.60 0.20 0.80], 'LineWidth', 1.8);

set(ax3, 'TickLabelInterpreter', 'latex', 'FontSize', 11);
xlabel(ax3, 'Time [s]', 'Interpreter', 'latex', 'FontSize', 12);
ylabel(ax3, '$\|\sigma(t)\|_2$ (Slack Vector Norm)', 'Interpreter', 'latex', 'FontSize', 12);
title(ax3, '\textbf{Feasibility Slack Variable Norm $\|\sigma(t)\|_2$ over Time}', 'Interpreter', 'latex', 'FontSize', 13);
grid(ax3, 'on');

%% ===== Figure 4: g-Vector Weight Distribution =====
fig4 = figure('Name', 'g-Vector Analysis');
t4 = tiledlayout(fig4, 2, 1, 'TileSpacing', 'compact', 'Padding', 'compact');

% Plot 4a: Weight histogram of g across all steps
ax4a = nexttile(t4);
histogram(ax4a, g_hist_matrix(:), 'BinWidth', 0.01, 'FaceColor', [0.2 0.5 0.7], 'EdgeColor', 'w');
set(ax4a, 'TickLabelInterpreter', 'latex', 'FontSize', 11);
ylabel(ax4a, 'Frequency Count', 'Interpreter', 'latex', 'FontSize', 12);
title(ax4a, '\textbf{Global $g$-Vector Weight Distribution Histogram}', 'Interpreter', 'latex', 'FontSize', 13);
grid(ax4a, 'on');

% Plot 4b: Stem plot of g at a transient step (e.g. step 40, when step reference transitions)
ax4b = nexttile(t4);
transient_step = round(N_sim / 3);
stem(ax4b, g_hist_matrix(transient_step, :), 'Color', [0.15 0.45 0.15], 'MarkerSize', 4, 'MarkerFaceColor', [0.15 0.45 0.15]);
set(ax4b, 'TickLabelInterpreter', 'latex', 'FontSize', 11);
xlabel(ax4b, 'Hankel Column Index', 'Interpreter', 'latex', 'FontSize', 12);
ylabel(ax4b, '$g$-value', 'Interpreter', 'latex', 'FontSize', 12);
title(ax4b, sprintf('\\textbf{$g$-Vector Configuration during Transient (Step %d, $t = %.2f$\\,s)}', ...
    transient_step+T_ini, (transient_step+T_ini)*Ts), 'Interpreter', 'latex', 'FontSize', 12);
grid(ax4b, 'on');

fprintf('\nAll diagnostic figures generated successfully from .mat files!\n');

%% ===== Helpers =====
function H = block_hankel(sig, L)
    sig = sig(:); nc = numel(sig)-L+1;
    H = zeros(L, nc);
    for i = 1:L; H(i,:) = sig(i:i+nc-1).'; end
end

function out_dyn = fetch_out_dyn(data)
    if isstruct(data) && isfield(data,'out_Dyn'); out_dyn = data.out_Dyn; return; end
    if isa(data,'Simulink.SimulationOutput')
        try; out_dyn = data.get('out_Dyn'); return; catch; end
    end
    if isstruct(data) && isfield(data,'eng'); out_dyn = data; return; end
    error('Could not find out_Dyn');
end

function u = extract_input_trace(out_dyn)
    if isfield(out_dyn,'cntrl')
        if isfield(out_dyn.cntrl,'Wfact'); u = sv(out_dyn.cntrl.Wfact); return; end
        if isfield(out_dyn.cntrl,'Wf'); u = sv(out_dyn.cntrl.Wf); return; end
    end
    error('No Wf trace found');
end

function y = measure_tracking_series(out_dyn, name)
    nl = lower(string(name));
    if contains(nl,'n_h'); y = sv(out_dyn.eng.Shaft.N_HPC);
    elseif contains(nl,'n_l'); y = sv(out_dyn.eng.Shaft.N_Fan);
    elseif contains(nl,'thrust')||contains(nl,'fnet'); y = sv(out_dyn.eng.Perf.Fnet);
    else; y = sv(out_dyn.eng.Shaft.N_HPC); end
end

function x = sv(s)
    if isstruct(s)&&isfield(s,'Data'); x = double(s.Data(:));
    elseif isa(s,'timeseries'); x = double(s.Data(:));
    else; x = double(s(:)); end
end
