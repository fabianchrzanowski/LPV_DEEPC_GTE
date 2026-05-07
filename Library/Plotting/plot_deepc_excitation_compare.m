% ========================================================================
% LPV-DeePC: Excitation Signal Comparison
%
% Generates synthetic training data at each operating point using the
% fast LPV surrogate model, with different excitation signals:
%   1. PRBS (original — gold standard)
%   2. Single step
%   3. Sinusoidal sweep
%   4. White noise (small amplitude)
%   5. Steady-state (no excitation at all)
%
% Then runs LPV-DeePC closed-loop with each dataset and compares the
% tracking performance + Hankel rank to show WHY PRBS is the standard
% choice for data-driven control.
% ========================================================================
clearvars; close all; clc
import casadi.*

proj = currentProject;
cd(proj.RootFolder)
plot_options;

%% ===== Load params =====
FAST_SURROGATE_MODE = true;
run('setup_simulation_params.m');

%% ===== Settings =====
u_min = 0.35; u_max = 2.0; du_max = 0.1;
q_track = 1; r_u = 1; r_du = 25;
lambda_g = 50; lambda_sigma = 10;
n_sys = 4;
rank_limit = L + n_sys;

%% ===== Load models =====
models = load(fullfile(proj.RootFolder,'Params','tables','matrices_DE_wf.mat'));
wf_data = load(fullfile(proj.RootFolder,'Params','test_files_fuel_only','test_identification','wf_values.mat'));
offset = models.offset;
sys_eng = models.sys_eng;
wf_values_all = wf_data.Wf_values(:)';

idx_track = find(contains(lower(string(sys_eng.OutputName)), "n_h"), 1);
if isempty(idx_track); idx_track = 2; end
track_name = sys_eng.OutputName{idx_track};

n_op = numel(wf_values_all);
hankel_wf = wf_values_all;

%% ===== Define Excitation Signals =====
% Each signal generates N_gen samples of fuel perturbation around the 
% equilibrium at each operating point.
N_gen = 2000;  % number of training samples per operating point
rng(42);       % reproducible

excitation_modes = {
    'PRBS',           @(wf_eq, amp) wf_eq + amp*(2*(rand(N_gen,1) > 0.5) - 1);
    'Single Step',    @(wf_eq, amp) [wf_eq*ones(N_gen/2,1); (wf_eq+amp)*ones(N_gen/2,1)];
    'Sinusoidal',     @(wf_eq, amp) wf_eq + amp*sin(2*pi*(1:N_gen)'/50);
    'White Noise',    @(wf_eq, amp) wf_eq + amp*randn(N_gen,1);
    'Steady-State',   @(wf_eq, amp) wf_eq*ones(N_gen,1);
};

n_modes = size(excitation_modes, 1);
mode_names = excitation_modes(:,1);
mode_fns = excitation_modes(:,2);

% Perturbation amplitude (fraction of operating point)
amp_fraction = 0.05;  % 5% perturbation around equilibrium

%% ===== Generate Synthetic Data + Build Hankels =====
fprintf('Generating synthetic training data for %d excitation modes...\n', n_modes);

all_hankel_U = cell(n_modes, 1);
all_hankel_Y = cell(n_modes, 1);
all_hankel_u_mean = cell(n_modes, 1);
all_hankel_y_mean = cell(n_modes, 1);
all_hankel_ranks = cell(n_modes, 1);
all_hankel_cols = cell(n_modes, 1);

for mode = 1:n_modes
    fprintf('\n--- %s ---\n', mode_names{mode});
    
    hU = cell(n_op, 1);
    hY = cell(n_op, 1);
    h_u_mean = zeros(n_op, 1);
    h_y_mean = zeros(n_op, 1);
    h_ranks = zeros(n_op, 1);
    h_cols = zeros(n_op, 1);
    
    for k = 1:n_op
        wf_eq = wf_values_all(k);
        amp = amp_fraction * wf_eq;
        
        % Generate excitation signal
        rng(42 + k);  % same seed per cloud for fairness
        u_train = mode_fns{mode}(wf_eq, amp);
        u_train = max(u_min, min(u_max, u_train));  % clamp
        
        % Simulate plant response using LPV surrogate
        x_state = zeros(size(models.A_dt, 1), 1);
        y_train = zeros(N_gen, 1);
        
        for i = 1:N_gen
            [y_train(i), x_state] = LPV_step_fast(x_state, u_train(i), wf_eq, models);
        end
        
        % Build Hankel in deviation space
        h_u_mean(k) = mean(u_train);
        h_y_mean(k) = mean(y_train);
        
        du = u_train - h_u_mean(k);
        dy = y_train - h_y_mean(k);
        Hu = block_hankel(du(:), L);
        Hy = block_hankel(dy(:), L);
        
        r = rank([Hu; Hy]);
        h_ranks(k) = r;
        h_cols(k) = size(Hu, 2);
        
        hU{k} = Hu;
        hY{k} = Hy;
        
        fprintf('  Wf=%.3f: %4d cols, rank=%2d (need>=%d) %s\n', ...
            wf_eq, h_cols(k), r, rank_limit, ...
            ternary(r >= rank_limit, 'OK', 'LOW'));
    end
    
    all_hankel_U{mode} = hU;
    all_hankel_Y{mode} = hY;
    all_hankel_u_mean{mode} = h_u_mean;
    all_hankel_y_mean{mode} = h_y_mean;
    all_hankel_ranks{mode} = h_ranks;
    all_hankel_cols{mode} = h_cols;
    
    valid = h_ranks > 0;
    fprintf('  Rank range: [%d .. %d] | Willems: %s\n', ...
        min(h_ranks(valid)), max(h_ranks(valid)), ...
        ternary(min(h_ranks(valid)) >= rank_limit, 'SATISFIED', 'VIOLATED'));
end

%% ===== Build CasADi Solver (ONCE) =====
fprintf('\nCompiling CasADi solver...\n');
n_cols_std = max(cellfun(@max, all_hankel_cols));

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
p_sym = [p_u_ini; p_y_ini; p_ref; p_u_prev; p_u0; p_Up(:); p_Yp(:); p_Uf(:); p_Yf(:)];

g_eq = [p_Up*g_var - p_u_ini;
        p_Yp*g_var - p_y_ini - sigma_var;
        p_Uf*g_var - U_var;
        p_Yf*g_var - Y_var];

g_ineq = SX.zeros(N_pred, 1);
for i = 1:N_pred
    if i == 1; g_ineq(i) = U_var(i) - p_u_prev;
    else;      g_ineq(i) = U_var(i) - U_var(i-1); end
end

obj = lambda_g*(g_var'*g_var) + lambda_sigma*(sigma_var'*sigma_var);
for i = 1:N_pred
    if i == 1; du_i = U_var(i) - p_u_prev;
    else;      du_i = U_var(i) - U_var(i-1); end
    obj = obj + q_track*(Y_var(i)-p_ref(i))^2 + r_u*(U_var(i)-p_u0)^2 + r_du*du_i^2;
end

g_all = [g_eq; g_ineq];
n_eq = length(g_eq);

nlp = struct('x', x, 'f', obj, 'g', g_all, 'p', p_sym);
opts = struct;
opts.ipopt.print_level = 0;
opts.print_time = 0;
opts.ipopt.max_iter = 200;
opts.ipopt.tol = 1e-6;

solver_casadi = nlpsol('deepc_excite', 'ipopt', nlp, opts);
lbg_val = [zeros(n_eq,1); -du_max*ones(N_pred,1)];
ubg_val = [zeros(n_eq,1);  du_max*ones(N_pred,1)];
n_x_total = n_cols_std + N_pred + N_pred + T_ini;

fprintf('CasADi solver compiled.\n');

%% ===== Run Closed-Loop for ALL Modes =====
y_results = cell(n_modes, 1);
u_results = cell(n_modes, 1);
rmse_values = zeros(n_modes, 1);

for mode = 1:n_modes
    fprintf('\n===== Running LPV-DeePC: %s =====\n', mode_names{mode});
    
    hU = all_hankel_U{mode};
    hY = all_hankel_Y{mode};
    h_u_mean = all_hankel_u_mean{mode};
    h_y_mean = all_hankel_y_mean{mode};
    h_cols = all_hankel_cols{mode};
    
    valid = ~cellfun(@isempty, hU);
    hU = hU(valid);
    hY = hY(valid);
    hwf = hankel_wf(valid);
    h_u_mean = h_u_mean(valid);
    h_y_mean = h_y_mean(valid);
    h_cols = h_cols(valid);
    n_valid = sum(valid);
    
    % Pad Hankels
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
    
    for t = 1:T_ini
        [y_hist(t), x_plant_fast] = LPV_step_fast(x_plant_fast, u0, u0, models);
    end
    
    wf_filtered = u0;
    for t = T_ini+1:N_total
        si = t - T_ini;
        
        wf_for_sched = max(min(u_prev, max(hwf)), min(hwf));
        if si == 1; wf_filtered = wf_for_sched;
        else; wf_filtered = 0.3*wf_for_sched + 0.7*wf_filtered; end
        [~, k_sel] = min(abs(hwf - wf_filtered));
        
        u_mean_k = h_u_mean(k_sel);
        y_mean_k = h_y_mean(k_sel);
        
        Hu_sel = hU{k_sel}; Hy_sel = hY{k_sel};
        Up_sel = Hu_sel(1:T_ini, :); Uf_sel = Hu_sel(T_ini+1:end, :);
        Yp_sel = Hy_sel(1:T_ini, :); Yf_sel = Hy_sel(T_ini+1:end, :);
        
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
            u_now = x_sol(n_cols_std + 1) + u_mean_k;
        catch
            u_now = u_prev;
        end
        
        if isnan(u_now); u_now = u_prev; end
        u_now = max(u_min, min(u_max, u_now));
        u_hist(t) = u_now;
        
        [y_hist(t), x_plant_fast] = LPV_step_fast(x_plant_fast, u_now, wf_filtered, models);
        u_prev = u_now;
    end
    
    y_results{mode} = y_hist;
    u_results{mode} = u_hist;
    
    idx_eval = T_ini+1:N_total;
    rmse_values(mode) = sqrt(mean((y_hist(idx_eval) - ref(idx_eval)).^2, 'omitnan'));
    fprintf('  RMSE: %.1f rpm\n', rmse_values(mode));
end

%% ===== FIGURE 1: Tracking Comparison =====
time_vec = (0:N_total-1)'*Ts;
control_time = time_vec(T_ini+1:end);

colors_cmp = [
    0.00 0.45 0.74;   % PRBS — blue
    0.85 0.33 0.10;   % Step — orange
    0.20 0.70 0.30;   % Sine — green
    0.93 0.69 0.13;   % Noise — yellow
    0.49 0.18 0.56;   % SS — violet
];

figure('Name','Excitation Signal Comparison: Tracking')
tiledlayout(2,1,'TileSpacing','compact','Padding','compact');

ax1 = nexttile;
plot(ax1, control_time, ref(T_ini+1:N_total), 'k--', 'LineWidth', 1.5, 'DisplayName', 'Reference');
hold(ax1, 'on');
for mode = 1:n_modes
    plot(ax1, control_time, y_results{mode}(T_ini+1:end), '-', ...
        'Color', colors_cmp(mode,:), 'LineWidth', 1.3, ...
        'DisplayName', sprintf('%s (RMSE=%.0f)', mode_names{mode}, rmse_values(mode)));
end
ylabel(ax1, track_name);
legend(ax1, 'Location', 'best');
title(ax1, 'LPV-DeePC: Impact of Training Excitation Signal');
grid(ax1, 'on');

ax2 = nexttile;
for mode = 1:n_modes
    stairs(ax2, control_time, u_results{mode}(T_ini+1:end), '-', ...
        'Color', colors_cmp(mode,:), 'LineWidth', 1.2);
    hold(ax2, 'on');
end
ylabel(ax2, '$W_f$ [pps]');
xlabel(ax2, 'Time [s]');
title(ax2, 'Applied Fuel Command');
grid(ax2, 'on');

%% ===== FIGURE 2: Rank Comparison Bar Chart =====
figure('Name','Hankel Rank by Excitation Signal')

% Show average rank across operating points for each mode
avg_ranks = zeros(n_modes, 1);
min_ranks = zeros(n_modes, 1);
for mode = 1:n_modes
    r = all_hankel_ranks{mode};
    avg_ranks(mode) = mean(r(r > 0));
    min_ranks(mode) = min(r(r > 0));
end

bar_data = [min_ranks, avg_ranks];
b = bar(categorical(mode_names), bar_data);
b(1).FaceColor = [0.85 0.33 0.10];
b(2).FaceColor = [0.00 0.45 0.74];
hold on;
yline(rank_limit, 'r--', 'LineWidth', 2, 'DisplayName', ...
    sprintf('Willems min rank = %d', rank_limit));
ylabel('Hankel matrix rank');
legend({'Min rank', 'Avg rank', sprintf('Required rank = %d', rank_limit)}, 'Location', 'best');
title('Hankel Rank vs Excitation Signal Type');
grid on;

%% ===== FIGURE 3: RMSE Bar Chart =====
figure('Name','RMSE by Excitation Signal')
b2 = bar(categorical(mode_names), rmse_values);
b2.FaceColor = 'flat';
for mode = 1:n_modes
    b2.CData(mode,:) = colors_cmp(mode,:);
end
ylabel('Tracking RMSE [rpm]');
title('LPV-DeePC Performance vs Training Excitation Signal');
grid on;

%% ===== Summary =====
fprintf('\n============================================\n');
fprintf('  EXCITATION SIGNAL COMPARISON SUMMARY\n');
fprintf('============================================\n');
fprintf('%-20s  %8s  %8s  %8s\n', 'Signal', 'MinRank', 'AvgRank', 'RMSE');
fprintf('%-20s  %8s  %8s  %8s\n', '------', '-------', '-------', '----');
for mode = 1:n_modes
    fprintf('%-20s  %8d  %8.1f  %8.1f\n', mode_names{mode}, min_ranks(mode), avg_ranks(mode), rmse_values(mode));
end
fprintf('Required rank (L+n): %d\n', rank_limit);
fprintf('============================================\n');

%% ===== Save =====
save('Results/results_excitation_comparison.mat', ...
    'mode_names', 'rmse_values', 'avg_ranks', 'min_ranks', ...
    'all_hankel_ranks', 'rank_limit', 'y_results', 'u_results');
fprintf('Saved to Results/results_excitation_comparison.mat\n');

%% ===== Helper Functions =====
function s = ternary(cond, a, b)
    if cond; s = a; else; s = b; end
end

function H = block_hankel(sig, L)
    sig = sig(:); nc = numel(sig)-L+1;
    H = zeros(L, nc);
    for i = 1:L; H(i,:) = sig(i:i+nc-1).'; end
end
