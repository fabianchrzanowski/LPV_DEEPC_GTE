% ========================================================================
% LPV-DeePC: Steady-State vs PRBS Data Comparison
%
% This script runs LPV-DeePC twice on the same trajectory:
%   1. Using ONLY the steady-state portion of each cloud (before PRBS)
%   2. Using ONLY the PRBS-excited portion of each cloud (after tstart)
%
% Demonstrates that DeePC fundamentally requires persistently exciting
% data and cannot function with boring steady-state recordings, even
% when plenty of data is available.
% ========================================================================
clearvars; close all; clc
import casadi.*

proj = currentProject;
cd(proj.RootFolder)
plot_options;

%% ===== Force surrogate mode and load params =====
FAST_SURROGATE_MODE = true;
run('setup_simulation_params.m');

%% ===== Settings =====
u_min = 0.35; u_max = 2.0; du_max = 0.1;
q_track = 1; r_u = 1; r_du = 25;
lambda_g = 50; lambda_sigma = 10;
n_sys = 4;
rank_limit = L + n_sys;

% Time boundary between steady-state and PRBS in each recording
tstart_prbs = 19.5;  % [s] — PRBS excitation begins here
n_ss_samples = floor(tstart_prbs / Ts);  % ~1300 samples of steady-state

%% ===== Load Training Data =====
training_folder = fullfile(proj.RootFolder,'Params','test_files_fuel_only','test_identification','AGTF30');
models = load(fullfile(proj.RootFolder,'Params','tables','matrices_DE_wf.mat'));
wf_data = load(fullfile(proj.RootFolder,'Params','test_files_fuel_only','test_identification','wf_values.mat'));
offset = models.offset;
sys_eng = models.sys_eng;
wf_values_all = wf_data.Wf_values(:)';

idx_track = find(contains(lower(string(sys_eng.OutputName)), "n_h"), 1);
if isempty(idx_track); idx_track = 2; end
track_name = sys_eng.OutputName{idx_track};

filelist_all = dir(fullfile(training_folder,'*.mat'));
[~,idx_sort] = sort({filelist_all.name});
filelist_all = filelist_all(idx_sort);
n_total_files = min(numel(filelist_all), numel(wf_values_all));

%% ===== Pre-load ALL raw data clouds =====
fprintf('Pre-loading %d raw data clouds...\n', n_total_files);
raw_u = cell(n_total_files, 1);
raw_y = cell(n_total_files, 1);
for k = 1:n_total_files
    data = load(fullfile(filelist_all(k).folder, filelist_all(k).name));
    out_dyn = fetch_out_dyn(data);
    u = extract_input_trace(out_dyn);
    y = measure_tracking_series(out_dyn, track_name);
    n_common = min(numel(u), numel(y));
    raw_u{k} = u(1:n_common);
    raw_y{k} = y(1:n_common);
end
fprintf('Done. Each cloud has ~%d samples. SS portion: first %d samples.\n', ...
    numel(raw_u{1}), n_ss_samples);

%% ===== Build Hankels for BOTH data sources =====
n_op = n_total_files;
hankel_wf = wf_values_all(1:n_op);

data_labels = {'Steady-State (before PRBS)', 'PRBS-Excited (after tstart)'};
n_modes = 2;

% Storage for both modes
all_hankel_U = cell(n_modes, 1);
all_hankel_Y = cell(n_modes, 1);
all_hankel_u_mean = cell(n_modes, 1);
all_hankel_y_mean = cell(n_modes, 1);
all_hankel_cols = cell(n_modes, 1);
all_hankel_ranks = cell(n_modes, 1);

for mode = 1:n_modes
    fprintf('\n===== Building Hankels: %s =====\n', data_labels{mode});
    
    hU = cell(n_op, 1);
    hY = cell(n_op, 1);
    h_u_mean = zeros(n_op, 1);
    h_y_mean = zeros(n_op, 1);
    h_cols = zeros(n_op, 1);
    h_ranks = zeros(n_op, 1);
    
    for k = 1:n_op
        u_full = raw_u{k};
        y_full = raw_y{k};
        
        if mode == 1
            % STEADY-STATE: use only samples BEFORE PRBS starts
            n_use = min(n_ss_samples, numel(u_full));
            u = u_full(1:n_use);
            y = y_full(1:n_use);
        else
            % PRBS: use only samples AFTER PRBS starts
            if numel(u_full) > n_ss_samples
                u = u_full(n_ss_samples+1:end);
                y = y_full(n_ss_samples+1:end);
            else
                u = u_full;
                y = y_full;
            end
        end
        
        if numel(u) < L
            warning('Wf=%.3f [%s]: not enough data (%d < %d)', ...
                hankel_wf(k), data_labels{mode}, numel(u), L);
            continue
        end
        
        h_u_mean(k) = mean(u);
        h_y_mean(k) = mean(y);
        
        du = u - h_u_mean(k);
        dy = y - h_y_mean(k);
        Hu = block_hankel(du(:), L);
        Hy = block_hankel(dy(:), L);
        
        r = rank([Hu; Hy]);
        h_ranks(k) = r;
        h_cols(k) = size(Hu, 2);
        
        hU{k} = Hu;
        hY{k} = Hy;
        
        fprintf('  Wf=%.3f: %4d cols, rank=%2d (need>=%d) %s\n', ...
            hankel_wf(k), h_cols(k), r, rank_limit, ...
            ternary(r >= rank_limit, 'OK', 'RANK DEFICIENT'));
    end
    
    all_hankel_U{mode} = hU;
    all_hankel_Y{mode} = hY;
    all_hankel_u_mean{mode} = h_u_mean;
    all_hankel_y_mean{mode} = h_y_mean;
    all_hankel_cols{mode} = h_cols;
    all_hankel_ranks{mode} = h_ranks;
    
    valid = ~cellfun(@isempty, hU);
    fprintf('  Valid: %d/%d | Rank range: [%d .. %d] | Willems: %s\n', ...
        sum(valid), n_op, min(h_ranks(valid)), max(h_ranks(valid)), ...
        ternary(min(h_ranks(valid)) >= rank_limit, 'SATISFIED', 'VIOLATED'));
end

%% ===== Build CasADi Solver (ONCE) =====
fprintf('\nCompiling parameterized CasADi solver...\n');

% Use the max columns from either mode
max_cols_either = max(cellfun(@max, all_hankel_cols));
n_cols_std = max_cols_either;

g_var     = SX.sym('g', n_cols_std, 1);
U_var     = SX.sym('U', N_pred, 1);
Y_var     = SX.sym('Y', N_pred, 1);
sigma_var = SX.sym('sigma', T_ini, 1);

p_u_ini  = SX.sym('u_ini', T_ini, 1);
p_y_ini  = SX.sym('y_ini', T_ini, 1);
p_ref    = SX.sym('ref', N_pred, 1);
p_u_prev = SX.sym('u_prev', 1, 1);
p_u0     = SX.sym('u0', 1, 1);
p_Up     = SX.sym('Up', T_ini, n_cols_std);
p_Yp     = SX.sym('Yp', T_ini, n_cols_std);
p_Uf     = SX.sym('Uf', N_pred, n_cols_std);
p_Yf     = SX.sym('Yf', N_pred, n_cols_std);

x = [g_var; U_var; Y_var; sigma_var];
p = [p_u_ini; p_y_ini; p_ref; p_u_prev; p_u0; p_Up(:); p_Yp(:); p_Uf(:); p_Yf(:)];

g_eq = [p_Up*g_var - p_u_ini;
        p_Yp*g_var - p_y_ini - sigma_var;
        p_Uf*g_var - U_var;
        p_Yf*g_var - Y_var];

g_ineq = SX.zeros(N_pred, 1);
for i = 1:N_pred
    if i == 1
        g_ineq(i) = U_var(i) - p_u_prev;
    else
        g_ineq(i) = U_var(i) - U_var(i-1);
    end
end

obj = lambda_g*(g_var'*g_var) + lambda_sigma*(sigma_var'*sigma_var);
for i = 1:N_pred
    if i == 1; du_i = U_var(i) - p_u_prev;
    else;      du_i = U_var(i) - U_var(i-1); end
    obj = obj + q_track*(Y_var(i)-p_ref(i))^2 + r_u*(U_var(i)-p_u0)^2 + r_du*du_i^2;
end

g_all = [g_eq; g_ineq];
n_eq = length(g_eq);

nlp = struct('x', x, 'f', obj, 'g', g_all, 'p', p);
opts = struct;
opts.ipopt.print_level = 0;
opts.print_time = 0;
opts.ipopt.max_iter = 200;
opts.ipopt.tol = 1e-6;

solver_casadi = nlpsol('deepc_ss_vs_prbs', 'ipopt', nlp, opts);
lbg_val = [zeros(n_eq,1); -du_max*ones(N_pred,1)];
ubg_val = [zeros(n_eq,1);  du_max*ones(N_pred,1)];
n_x_total = n_cols_std + N_pred + N_pred + T_ini;

fprintf('CasADi solver compiled.\n');

%% ===== Run Closed-Loop for BOTH Modes =====
y_results = cell(n_modes, 1);
u_results = cell(n_modes, 1);
rmse_values = zeros(n_modes, 1);

for mode = 1:n_modes
    fprintf('\n===== Running LPV-DeePC with: %s =====\n', data_labels{mode});
    
    hU = all_hankel_U{mode};
    hY = all_hankel_Y{mode};
    h_u_mean = all_hankel_u_mean{mode};
    h_y_mean = all_hankel_y_mean{mode};
    h_cols = all_hankel_cols{mode};
    
    % Remove empty entries
    valid = ~cellfun(@isempty, hU);
    hU = hU(valid);
    hY = hY(valid);
    hwf = hankel_wf(valid);
    h_u_mean = h_u_mean(valid);
    h_y_mean = h_y_mean(valid);
    h_cols = h_cols(valid);
    n_valid = sum(valid);
    
    % Pad Hankels to n_cols_std
    for k = 1:n_valid
        nc = size(hU{k}, 2);
        if nc < n_cols_std
            hU{k} = [hU{k}, zeros(L, n_cols_std - nc)];
            hY{k} = [hY{k}, zeros(L, n_cols_std - nc)];
        elseif nc > n_cols_std
            hU{k} = hU{k}(:, 1:n_cols_std);
            hY{k} = hY{k}(:, 1:n_cols_std);
        end
    end
    
    % Simulate
    y_hist = nan(N_total, 1);
    u_hist = zeros(N_total, 1);
    y_hist(1:T_ini) = z_idle;
    u_hist(1:T_ini) = u0;
    
    x_plant_fast = zeros(size(models.A_dt,1), 1);
    x0_guess = zeros(n_x_total, 1);
    u_prev = u0;
    
    % Warm-up
    for t = 1:T_ini
        [y_hist(t), x_plant_fast] = LPV_step_fast(x_plant_fast, u0, u0, models);
    end
    
    % Control loop
    wf_filtered = u0;
    for t = T_ini+1:N_total
        si = t - T_ini;
        
        wf_for_sched = max(min(u_prev, max(hwf)), min(hwf));
        if si == 1
            wf_filtered = wf_for_sched;
        else
            wf_filtered = 0.3*wf_for_sched + 0.7*wf_filtered;
        end
        [~, k_sel] = min(abs(hwf - wf_filtered));
        
        u_mean_k = h_u_mean(k_sel);
        y_mean_k = h_y_mean(k_sel);
        
        Hu_sel = hU{k_sel};
        Hy_sel = hY{k_sel};
        Up_sel = Hu_sel(1:T_ini, :);
        Uf_sel = Hu_sel(T_ini+1:end, :);
        Yp_sel = Hy_sel(1:T_ini, :);
        Yf_sel = Hy_sel(T_ini+1:end, :);
        
        du_ini = u_hist(t-T_ini:t-1) - u_mean_k;
        dy_ini = y_hist(t-T_ini:t-1) - y_mean_k;
        dref_seg = ref(t:t+N_pred-1) - y_mean_k;
        du_prev = u_prev - u_mean_k;
        du0 = u0 - u_mean_k;
        
        p_val = [du_ini; dy_ini; dref_seg; du_prev; du0; ...
                 Up_sel(:); Yp_sel(:); Uf_sel(:); Yf_sel(:)];
        
        lbx = [-inf(n_cols_std,1); (u_min-u_mean_k)*ones(N_pred,1); -inf(N_pred,1); -inf(T_ini,1)];
        ubx = [ inf(n_cols_std,1); (u_max-u_mean_k)*ones(N_pred,1);  inf(N_pred,1);  inf(T_ini,1)];
        
        try
            sol = solver_casadi('x0', x0_guess, 'p', p_val, ...
                'lbx', lbx, 'ubx', ubx, 'lbg', lbg_val, 'ubg', ubg_val);
            x_sol = full(sol.x);
            x0_guess = x_sol;
            du_now = x_sol(n_cols_std + 1);
            u_now = du_now + u_mean_k;
        catch
            u_now = u_prev;
        end
        
        if isnan(u_now); u_now = u_prev; end
        u_now = max(u_min, min(u_max, u_now));
        u_hist(t) = u_now;
        
        [y_hist(t), x_plant_fast] = LPV_step_fast(x_plant_fast, u_now, wf_filtered, models);
        u_prev = u_now;
        
        if mod(si, 50) == 0
            fprintf('  Step %d/%d | u=%.3f | wf_sched=%.3f\n', si, N_sim, u_now, hwf(k_sel));
        end
    end
    
    y_results{mode} = y_hist;
    u_results{mode} = u_hist;
    
    idx_eval = T_ini+1 : N_total;
    rmse_values(mode) = sqrt(mean((y_hist(idx_eval) - ref(idx_eval)).^2, 'omitnan'));
    fprintf('  RMSE: %.1f rpm\n', rmse_values(mode));
end

%% ===== FIGURE 1: Tracking Comparison =====
time = (0:N_total-1)'*Ts;
control_time = time(T_ini+1:end);

figure('Name','Steady-State vs PRBS Data: LPV-DeePC Tracking')
tiledlayout(3,1,'TileSpacing','compact','Padding','compact');

% Tracking
ax1 = nexttile;
plot(ax1, control_time, ref(T_ini+1:N_total), 'k--', 'LineWidth', 1.5, 'DisplayName', 'Reference');
hold(ax1, 'on');
colors_cmp = [0.85 0.33 0.10; 0.00 0.45 0.74];
for mode = 1:n_modes
    plot(ax1, control_time, y_results{mode}(T_ini+1:end), '-', ...
        'Color', colors_cmp(mode,:), 'LineWidth', 1.4, ...
        'DisplayName', sprintf('%s (RMSE=%.0f)', data_labels{mode}, rmse_values(mode)));
end
ylabel(ax1, track_name);
legend(ax1, 'Location', 'best');
title(ax1, 'LPV-DeePC: Steady-State Data vs PRBS Data');
grid(ax1, 'on');

% Control effort
ax2 = nexttile;
for mode = 1:n_modes
    stairs(ax2, control_time, u_results{mode}(T_ini+1:end), '-', ...
        'Color', colors_cmp(mode,:), 'LineWidth', 1.4, ...
        'DisplayName', data_labels{mode});
    hold(ax2, 'on');
end
ylabel(ax2, '$W_f$ [pps]');
legend(ax2, 'Location', 'best');
title(ax2, 'Applied Fuel Command');
grid(ax2, 'on');

% Rank comparison bar chart
ax3 = nexttile;
ranks_ss   = all_hankel_ranks{1};
ranks_prbs = all_hankel_ranks{2};
valid_both = (ranks_ss > 0) & (ranks_prbs > 0);
bar_data = [ranks_ss(valid_both), ranks_prbs(valid_both)];
b = bar(ax3, hankel_wf(valid_both), bar_data);
b(1).FaceColor = colors_cmp(1,:);
b(2).FaceColor = colors_cmp(2,:);
hold(ax3, 'on');
yline(ax3, rank_limit, 'r--', 'LineWidth', 2, 'DisplayName', ...
    sprintf('Willems min rank = %d', rank_limit));
xlabel(ax3, '$W_f$ [pps]');
ylabel(ax3, 'Hankel rank');
legend(ax3, [data_labels, {sprintf('Min rank needed = %d', rank_limit)}], 'Location', 'best');
title(ax3, 'Hankel Matrix Rank at Each Operating Point');
grid(ax3, 'on');

%% ===== Summary =====
fprintf('\n======================================\n');
fprintf('  STEADY-STATE vs PRBS COMPARISON\n');
fprintf('======================================\n');
fprintf('  %-30s  RMSE = %8.1f rpm\n', data_labels{1}, rmse_values(1));
fprintf('  %-30s  RMSE = %8.1f rpm\n', data_labels{2}, rmse_values(2));
fprintf('  SS rank range:   [%d .. %d]  (need >= %d)\n', ...
    min(ranks_ss(ranks_ss>0)), max(ranks_ss), rank_limit);
fprintf('  PRBS rank range: [%d .. %d]  (need >= %d)\n', ...
    min(ranks_prbs(ranks_prbs>0)), max(ranks_prbs), rank_limit);
fprintf('======================================\n');

%% ===== Save =====
save('Results/results_ss_vs_prbs.mat', ...
    'y_results', 'u_results', 'rmse_values', 'data_labels', ...
    'all_hankel_ranks', 'hankel_wf', 'rank_limit', 'time', 'ref');
fprintf('Saved to Results/results_ss_vs_prbs.mat\n');

%% ===== Helper Functions =====
function s = ternary(cond, a, b)
    if cond; s = a; else; s = b; end
end

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
