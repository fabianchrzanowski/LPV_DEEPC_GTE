% ========================================================================
% LPV-DeePC Data Sensitivity Sweep
%
% Automated sweep over:
%   1. max_cols_per_wf  — tests the Persistency of Excitation cliff
%   2. wf_decimation_factor — tests how sparse the LPV scheduling can be
%
% Uses FAST_SURROGATE_MODE for speed. Compiles one parameterized CasADi
% solver, then re-uses it for every sweep point.
%
% Output: dissertation-quality analytical figures + saved .mat results.
% ========================================================================
clearvars; close all; clc
import casadi.*

proj = currentProject;
cd(proj.RootFolder)
plot_options;

%% ===== Force surrogate mode and load params =====
FAST_SURROGATE_MODE = true;
run('setup_simulation_params.m');

%% ===== Sweep Grid =====
col_sweeps        = [15, 25, 35, 40, 44, 50, 60, 80, 100, 150, 300];
decimation_sweeps = [1, 2, 3];
n_sys = 4;
% Theoretical minimum columns for persistency of excitation
rank_limit = L + n_sys;  % L = T_ini + N_pred, n_sys = 4

%% ===== Common Controller Settings =====
u_min = 0.35; u_max = 2.0; du_max = 0.1;
q_track = 1; r_u = 1; r_du = 25;
lambda_g = 50; lambda_sigma = 10;


%% ===== Load Training Data (ONCE) =====
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
fprintf('Done.\n');

%% ===== Build CasADi Solver (ONCE) =====
fprintf('Compiling parameterized CasADi solver...\n');
n_cols_std = max(col_sweeps);  % pad all Hankels to this width

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

% Constraints
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

% Objective
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

solver_casadi = nlpsol('deepc_sweep', 'ipopt', nlp, opts);
lbg = [zeros(n_eq,1); -du_max*ones(N_pred,1)];
ubg = [zeros(n_eq,1);  du_max*ones(N_pred,1)];
n_x_total = n_cols_std + N_pred + N_pred + T_ini;

fprintf('CasADi solver compiled.\n');

%% ===== MAIN SWEEP =====
n_dec = length(decimation_sweeps);
n_col = length(col_sweeps);
rmse_results  = nan(n_dec, n_col);
solve_results = nan(n_dec, n_col);
rank_min_results = nan(n_dec, n_col);  % minimum rank across operating points
rank_max_results = nan(n_dec, n_col);  % maximum rank across operating points

total_runs = n_dec * n_col;
run_count = 0;
sweep_tic = tic;

for d_idx = 1:n_dec
    wf_dec = decimation_sweeps(d_idx);
    
    % Decimate operating points
    keep_idx = 1:wf_dec:n_total_files;
    n_op = length(keep_idx);
    hankel_wf = wf_values_all(keep_idx);
    
    % Build full Hankels for this decimation (computed once per decimation)
    hankel_U_full = cell(n_op, 1);
    hankel_Y_full = cell(n_op, 1);
    hankel_u_mean = zeros(n_op, 1);
    hankel_y_mean = zeros(n_op, 1);
    
    for k = 1:n_op
        kk = keep_idx(k);
        u = raw_u{kk};
        y = raw_y{kk};
        if numel(u) < L; continue; end
        
        hankel_u_mean(k) = mean(u);
        hankel_y_mean(k) = mean(y);
        du = u - hankel_u_mean(k);
        dy = y - hankel_y_mean(k);
        hankel_U_full{k} = block_hankel(du(:), L);
        hankel_Y_full{k} = block_hankel(dy(:), L);
    end
    
    % Remove empties
    valid = ~cellfun(@isempty, hankel_U_full);
    hankel_U_full_v = hankel_U_full(valid);
    hankel_Y_full_v = hankel_Y_full(valid);
    hankel_wf_v     = hankel_wf(valid);
    hankel_u_mean_v = hankel_u_mean(valid);
    hankel_y_mean_v = hankel_y_mean(valid);
    n_valid = sum(valid);
    
    for c_idx = 1:n_col
        max_cols = col_sweeps(c_idx);
        run_count = run_count + 1;
        
        % Chop and pad Hankels, compute rank of each
        hankel_U = cell(n_valid, 1);
        hankel_Y = cell(n_valid, 1);
        hankel_cols = zeros(n_valid, 1);
        hankel_ranks = zeros(n_valid, 1);
        for k = 1:n_valid
            Hu = hankel_U_full_v{k};
            Hy = hankel_Y_full_v{k};
            nc = size(Hu, 2);
            if nc > max_cols
                keep = unique(round(linspace(1, nc, max_cols)));
                Hu = Hu(:, keep);
                Hy = Hy(:, keep);
            end
            nc_actual = size(Hu, 2);
            % Compute rank of the CHOPPED Hankel (before padding)
            hankel_ranks(k) = rank([Hu; Hy]);
            % Zero-pad to n_cols_std
            hankel_U{k} = [Hu, zeros(L, n_cols_std - nc_actual)];
            hankel_Y{k} = [Hy, zeros(L, n_cols_std - nc_actual)];
            hankel_cols(k) = nc_actual;
        end
        
        r_min = min(hankel_ranks);
        r_max = max(hankel_ranks);
        rank_min_results(d_idx, c_idx) = r_min;
        rank_max_results(d_idx, c_idx) = r_max;
        willems_ok = ternary(r_min >= rank_limit, 'YES', 'NO');
        
        fprintf('[%d/%d] dec=%d, cols=%d, n_op=%d | rank=[%d..%d] need>=%d Willems=%s ... ', ...
            run_count, total_runs, wf_dec, max_cols, n_valid, ...
            r_min, r_max, rank_limit, willems_ok);
        
        % ===== Run closed-loop simulation =====
        y_hist = nan(N_total, 1);
        u_hist = zeros(N_total, 1);
        y_hist(1:T_ini) = z_idle;
        u_hist(1:T_ini) = u0;
        
        x_plant_fast = zeros(size(models.A_dt,1), 1);
        x0_guess = zeros(n_x_total, 1);
        u_prev = u0;
        solve_times = zeros(N_sim, 1);
        
        % Warm-up
        for t = 1:T_ini
            [y_hist(t), x_plant_fast] = LPV_step_fast(x_plant_fast, u0, u0, models);
        end
        
        % Control loop
        wf_filtered = u0;
        for t = T_ini+1:N_total
            si = t - T_ini;
            
            % Scheduler
            wf_for_sched = max(min(u_prev, max(hankel_wf_v)), min(hankel_wf_v));
            if si == 1
                wf_filtered = wf_for_sched;
            else
                wf_filtered = 0.3*wf_for_sched + 0.7*wf_filtered;
            end
            [~, k_sel] = min(abs(hankel_wf_v - wf_filtered));
            
            u_mean_k = hankel_u_mean_v(k_sel);
            y_mean_k = hankel_y_mean_v(k_sel);
            
            Hu_sel = hankel_U{k_sel};
            Hy_sel = hankel_Y{k_sel};
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
            
            tic;
            try
                sol = solver_casadi('x0', x0_guess, 'p', p_val, ...
                    'lbx', lbx, 'ubx', ubx, 'lbg', lbg, 'ubg', ubg);
                solve_times(si) = toc;
                x_sol = full(sol.x);
                x0_guess = x_sol;
                du_now = x_sol(n_cols_std + 1);
                u_now = du_now + u_mean_k;
            catch
                solve_times(si) = toc;
                u_now = u_prev;
            end
            
            if isnan(u_now); u_now = u_prev; end
            u_now = max(u_min, min(u_max, u_now));
            u_hist(t) = u_now;
            
            [y_hist(t), x_plant_fast] = LPV_step_fast(x_plant_fast, u_now, wf_filtered, models);
            u_prev = u_now;
        end
        
        % Compute RMSE
        idx_eval = T_ini+1 : N_total;
        rmse = sqrt(mean((y_hist(idx_eval) - ref(idx_eval)).^2, 'omitnan'));
        rmse_results(d_idx, c_idx) = rmse;
        solve_results(d_idx, c_idx) = mean(solve_times)*1000;
        
        fprintf('RMSE = %8.1f rpm | solve = %.1f ms\n', rmse, solve_results(d_idx, c_idx));
    end
end

elapsed = toc(sweep_tic);
fprintf('\n===== Sweep completed in %.1f s (%d runs) =====\n', elapsed, total_runs);

%% ===== Save Results =====
save('Results/results_deepc_data_sweep.mat', ...
    'col_sweeps', 'decimation_sweeps', 'rmse_results', 'solve_results', ...
    'rank_limit', 'L', 'n_sys');
fprintf('Saved to Results/results_deepc_data_sweep.mat\n');

%% ===== FIGURE 1: Persistency of Excitation Cliff =====
figure('Name','Persistency of Excitation Analysis')

colors_sweep = [
    0.00 0.45 0.74;   % blue
    0.85 0.33 0.10;   % orange
    0.20 0.70 0.30;   % green
    0.93 0.69 0.13;   % yellow
    0.49 0.18 0.56;   % violet
    0.30 0.75 0.93;   % cyan
];

for d_idx = 1:n_dec
    c = colors_sweep(mod(d_idx-1, size(colors_sweep,1))+1, :);
    plot(col_sweeps, rmse_results(d_idx, :), '-o', ...
        'Color', c, 'MarkerFaceColor', c, 'LineWidth', 1.8, 'MarkerSize', 6, ...
        'DisplayName', sprintf('Decimation = %d (%d ops)', ...
            decimation_sweeps(d_idx), ...
            length(1:decimation_sweeps(d_idx):n_total_files)));
    hold on;
end

xline(rank_limit, 'r--', 'LineWidth', 2, ...
    'DisplayName', sprintf('Theoretical min ($L + n = %d$)', rank_limit));

set(gca, 'YScale', 'log');
xlabel('Number of Hankel columns per operating point');
ylabel('Tracking RMSE [rpm]');
title('LPV-DeePC: Persistency of Excitation Sensitivity');
legend('Location', 'northeast');
grid on;

%% ===== FIGURE 2: Decimation Sensitivity (at best column count) =====
figure('Name','LPV Scheduling Density Sensitivity')

% Plot RMSE vs decimation at a few representative column counts
col_show_vals = [50, 100, 300];
col_show_idx = zeros(size(col_show_vals));
for i = 1:numel(col_show_vals)
    [~, col_show_idx(i)] = min(abs(col_sweeps - col_show_vals(i)));
end

n_ops_per_dec = arrayfun(@(d) length(1:d:n_total_files), decimation_sweeps);

for i = 1:numel(col_show_idx)
    ci = col_show_idx(i);
    c = colors_sweep(mod(i-1, size(colors_sweep,1))+1, :);
    plot(n_ops_per_dec, rmse_results(:, ci), '-s', ...
        'Color', c, 'MarkerFaceColor', c, 'LineWidth', 1.8, 'MarkerSize', 8, ...
        'DisplayName', sprintf('%d columns', col_sweeps(ci)));
    hold on;
end

xlabel('Number of operating points used');
ylabel('Tracking RMSE [rpm]');
title('LPV-DeePC: Scheduling Density Sensitivity');
legend('Location', 'northeast');
set(gca, 'XDir', 'reverse');  % More ops on left (better), fewer on right (worse)
grid on;

%% ===== FIGURE 3: Heatmap =====
figure('Name','Data Sensitivity Heatmap')

n_ops_labels = arrayfun(@(d) sprintf('%d ops (dec=%d)', ...
    length(1:d:n_total_files), d), decimation_sweeps, 'UniformOutput', false);

imagesc(rmse_results);
colormap(flipud(hot));
colorbar;
set(gca, 'XTick', 1:n_col, 'XTickLabel', col_sweeps, 'XTickLabelRotation', 45);
set(gca, 'YTick', 1:n_dec, 'YTickLabel', n_ops_labels);
xlabel('Max Hankel columns per operating point');
ylabel('Operating point density');
title('LPV-DeePC Tracking RMSE Heatmap [rpm]');

% Annotate each cell with RMSE value
for d_idx = 1:n_dec
    for c_idx = 1:n_col
        val = rmse_results(d_idx, c_idx);
        if val < median(rmse_results(:))
            txt_color = 'w';
        else
            txt_color = 'k';
        end
        text(c_idx, d_idx, sprintf('%.0f', val), ...
            'HorizontalAlignment', 'center', 'Color', txt_color, 'FontSize', 8);
    end
end

fprintf('\nAll figures generated.\n');

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
