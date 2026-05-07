% ========================================================================
% LPV-DeePC from STEP RESPONSE data only
%
% Instead of loading PRBS clouds from sys_id, this script generates
% synthetic step-response data at each operating point using the fast
% LPV surrogate, builds Hankels from it, and runs closed-loop control.
%
% Toggle FAST_SURROGATE_MODE to test on surrogate vs real Simulink plant.
% ========================================================================
clearvars; close all; clc
import casadi.*

proj = currentProject;
cd(proj.RootFolder)
plot_options;

%% ===== Load params =====
FAST_SURROGATE_MODE = true;  % Set false to test on nonlinear Simulink
run('setup_simulation_params.m');

%% ===== Settings =====
u_min = 0.35; u_max = 2.0; du_max = 0.1;
q_track = 1; r_u = 1; r_du = 25;
lambda_g = 50; lambda_sigma = 10;
n_sys = 4;
rank_limit = L + n_sys;

%% ===== Generate step-response training data =====
models = load(fullfile(proj.RootFolder,'Params','tables','matrices_DE_wf.mat'));
wf_data = load(fullfile(proj.RootFolder,'Params','test_files_fuel_only','test_identification','wf_values.mat'));
offset = models.offset;
sys_eng = models.sys_eng;
wf_values = wf_data.Wf_values(:)';

idx_track = find(contains(lower(string(sys_eng.OutputName)), "n_h"), 1);
if isempty(idx_track); idx_track = 2; end
track_name = sys_eng.OutputName{idx_track};

n_op = numel(wf_values);
N_gen = 2000;         % samples per operating point
amp_fraction = 0.05;  % 5% step amplitude

fprintf('Generating step-response data at %d operating points...\n', n_op);

hankel_U = cell(n_op, 1);
hankel_Y = cell(n_op, 1);
hankel_wf = wf_values;
hankel_u_mean = zeros(n_op, 1);
hankel_y_mean = zeros(n_op, 1);
hankel_cols = zeros(n_op, 1);
hankel_rank = zeros(n_op, 1);

for k = 1:n_op
    wf_eq = wf_values(k);
    amp = amp_fraction * wf_eq;
    
    % Generate step input: hold at wf_eq, then step to wf_eq + amp
    u_train = [wf_eq*ones(N_gen/2, 1); (wf_eq + amp)*ones(N_gen/2, 1)];
    u_train = max(u_min, min(u_max, u_train));
    
    % Simulate using fast surrogate
    x_state = zeros(size(models.A_dt, 1), 1);
    y_train = zeros(N_gen, 1);
    for i = 1:N_gen
        [y_train(i), x_state] = LPV_step_fast(x_state, u_train(i), wf_eq, models);
    end
    
    % Build Hankel in deviation space
    hankel_u_mean(k) = mean(u_train);
    hankel_y_mean(k) = mean(y_train);
    du = u_train - hankel_u_mean(k);
    dy = y_train - hankel_y_mean(k);
    Hu = block_hankel(du(:), L);
    Hy = block_hankel(dy(:), L);
    
    r = rank([Hu; Hy]);
    hankel_rank(k) = r;
    hankel_cols(k) = size(Hu, 2);
    hankel_U{k} = Hu;
    hankel_Y{k} = Hy;
    
    fprintf('  Wf=%.3f: %d cols, rank=%d (need>=%d) %s\n', ...
        wf_eq, hankel_cols(k), r, rank_limit, ...
        ternary(r >= rank_limit, 'OK', 'LOW'));
end

%% ===== Build CasADi Solver =====
fprintf('\nCompiling CasADi solver...\n');
n_cols_std = max(hankel_cols);

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

nlp = struct('x', x, 'f', obj, 'g', g_all, 'p', p);
opts = struct; opts.ipopt.print_level = 0; opts.print_time = 0;
opts.ipopt.max_iter = 200; opts.ipopt.tol = 1e-6;

S = struct();
S.solver = nlpsol('deepc_step', 'ipopt', nlp, opts);
S.lbg = [zeros(n_eq,1); -du_max*ones(N_pred,1)];
S.ubg = [zeros(n_eq,1);  du_max*ones(N_pred,1)];
S.n_x = n_cols_std + N_pred + N_pred + T_ini;
S.n_cols = n_cols_std;

fprintf('Done.\n');

%% ===== Pad Hankels =====
for k = 1:n_op
    nc = size(hankel_U{k}, 2);
    if nc < n_cols_std
        hankel_U{k} = [hankel_U{k}, zeros(L, n_cols_std - nc)];
        hankel_Y{k} = [hankel_Y{k}, zeros(L, n_cols_std - nc)];
    end
end

%% ===== Simulink setup (if needed) =====
if ~FAST_SURROGATE_MODE
    MWS = struct();
    MWS.engName = 'AGTF30';
    MWS = setup_Controller(MWS);
    MWS = setup_AllEng(MWS);
    busVars = load(fullfile(proj.RootFolder,'Params','Eng_Bus.mat'));
    fns = fieldnames(busVars);
    for ib = 1:numel(fns); assignin('base',fns{ib},busVars.(fns{ib})); end
    MWS.In.Ts = Ts; MWS.In.Tsim = Ts;
    MWS.In.Alt  = timeseries([0;0],[0 Ts],'Name','Alt');
    MWS.In.MN   = timeseries([0;0],[0 Ts],'Name','MN');
    MWS.In.dT   = timeseries([0;0],[0 Ts],'Name','dT');
    MWS.In.dVBV = timeseries([0;0],[0 Ts],'Name','dVBV');
    MWS.In.dNz  = timeseries([0;0],[0 Ts],'Name','dNz');
    MWS.In.Wf   = timeseries([u0;u0],[0 Ts],'Name','Wf');
    MWS = AGTF30_initial_conditions(MWS);
end
model = 'AGTF30SysDyn';

%% ===== Closed-loop =====
u_hist = zeros(N_total, 1);
y_hist = nan(N_total, 1);
u_hist(1:T_ini) = u0;
u_prev = u0;
if FAST_SURROGATE_MODE
    x_plant_fast = zeros(size(models.A_dt,1), 1);
end
x0_guess = zeros(S.n_x, 1);
solve_times = zeros(N_sim, 1);
wf_scheduler_hist = zeros(N_sim, 1);

% Warm-up
fprintf('Warm-up (%d steps)...\n', T_ini);
for t = 1:T_ini
    if FAST_SURROGATE_MODE
        [y_hist(t), x_plant_fast] = LPV_step_fast(x_plant_fast, u0, u0, models);
    else
        t_step = (0:t)'*Ts; u_step = [u0; u_hist(1:t)];
        MWS.In.Tsim = t_step(end);
        MWS.In.Wf = timeseries(u_step, t_step, 'Name','Wf');
        simIn = Simulink.SimulationInput(model);
        simIn = simIn.setVariable('MWS',MWS);
        simIn = simIn.setModelParameter('StartTime','0','StopTime',num2str(t_step(end)),'ReturnWorkspaceOutputs','on');
        simOut = sim(simIn);
        y_hist(t) = measure_tracking_output(fetch_out_dyn(simOut), track_name);
    end
end

% Control loop
wf_filtered = u0;
fprintf('LPV-DeePC (step data) closed-loop (%d steps)...\n', N_sim);
for t = T_ini+1:N_total
    si = t - T_ini;
    
    wf_for_sched = max(min(u_prev, max(hankel_wf)), min(hankel_wf));
    if si == 1; wf_filtered = wf_for_sched;
    else; wf_filtered = 0.3*wf_for_sched + 0.7*wf_filtered; end
    [~, k_sel] = min(abs(hankel_wf - wf_filtered));
    wf_scheduler_hist(si) = hankel_wf(k_sel);
    
    u_mean_k = hankel_u_mean(k_sel);
    y_mean_k = hankel_y_mean(k_sel);
    
    Hu_sel = hankel_U{k_sel}; Hy_sel = hankel_Y{k_sel};
    Up_sel = Hu_sel(1:T_ini,:); Uf_sel = Hu_sel(T_ini+1:end,:);
    Yp_sel = Hy_sel(1:T_ini,:); Yf_sel = Hy_sel(T_ini+1:end,:);
    
    du_ini = u_hist(t-T_ini:t-1) - u_mean_k;
    dy_ini = y_hist(t-T_ini:t-1) - y_mean_k;
    dref_seg = ref(t:t+N_pred-1) - y_mean_k;
    du_prev = u_prev - u_mean_k;
    du0 = u0 - u_mean_k;
    
    p_val = [du_ini; dy_ini; dref_seg; du_prev; du0; Up_sel(:); Yp_sel(:); Uf_sel(:); Yf_sel(:)];
    
    lbx = [-inf(S.n_cols,1); (u_min-u_mean_k)*ones(N_pred,1); -inf(N_pred,1); -inf(T_ini,1)];
    ubx = [ inf(S.n_cols,1); (u_max-u_mean_k)*ones(N_pred,1);  inf(N_pred,1);  inf(T_ini,1)];
    
    tic;
    sol = S.solver('x0', x0_guess, 'p', p_val, 'lbx', lbx, 'ubx', ubx, 'lbg', S.lbg, 'ubg', S.ubg);
    solve_times(si) = toc;
    
    x_sol = full(sol.x);
    x0_guess = x_sol;
    u_now = x_sol(S.n_cols + 1) + u_mean_k;
    
    if isnan(u_now); u_now = u_prev; end
    u_now = max(u_min, min(u_max, u_now));
    u_hist(t) = u_now;
    
    if FAST_SURROGATE_MODE
        [y_hist(t), x_plant_fast] = LPV_step_fast(x_plant_fast, u_now, wf_filtered, models);
    else
        t_step = (0:t)'*Ts; u_step = [u0; u_hist(1:t)];
        MWS.In.Tsim = t_step(end);
        MWS.In.Wf = timeseries(u_step, t_step, 'Name','Wf');
        simIn = Simulink.SimulationInput(model);
        simIn = simIn.setVariable('MWS',MWS);
        simIn = simIn.setModelParameter('StartTime','0','StopTime',num2str(t_step(end)),'ReturnWorkspaceOutputs','on');
        simOut = sim(simIn);
        y_hist(t) = measure_tracking_output(fetch_out_dyn(simOut), track_name);
    end
    
    u_prev = u_now;
    if mod(si,50)==0
        fprintf('  Step %d/%d  Wf=%.3f -> Hankel #%d (Wf=%.3f)  solve=%.1fms\n', ...
            si, N_sim, u_prev, k_sel, hankel_wf(k_sel), solve_times(si)*1000);
    end
end

%% ===== Plot =====
time = (0:N_total-1)'*Ts;
ref_plot = ref(1:N_total);
control_time = time(T_ini+1:end);

figure('Name','LPV-DeePC from STEP DATA')
subplot(3,1,1)
plot(control_time, ref_plot(T_ini+1:N_total), 'k--', 'LineWidth', 1.5)
hold on
plot(control_time, y_hist(T_ini+1:end), 'LineWidth', 1.4)
grid on; ylabel(track_name);
legend('Reference', 'Step-trained DeePC', 'Location', 'best');
title('LPV-DeePC tracking (trained from step response only)')

subplot(3,1,2)
stairs(control_time, u_hist(T_ini+1:end), 'LineWidth', 1.4)
grid on; ylabel('$W_f$ [pps]'); title('Applied fuel command')

subplot(3,1,3)
stairs(1:N_sim, wf_scheduler_hist, 'LineWidth', 1.4)
grid on; xlabel('MPC step'); ylabel('$W_f$ scheduled');
title('LPV Scheduling variable ($\rho$)')

y_ctrl = y_hist(T_ini+1:end);
ref_ctrl = ref_plot(T_ini+1:N_total);
rmse = sqrt(mean((y_ctrl - ref_ctrl).^2, 'omitnan'));
fprintf('\n===== LPV-DeePC (Step Data) =====\n');
fprintf('RMSE: %.1f rpm\n', rmse);
fprintf('Mean solve: %.1f ms\n', mean(solve_times)*1000);

%% ===== Save =====
controller_name = 'LPV DeePC (Step)';
save('Results/results_lpv_deepc_step.mat', 'y_hist', 'u_hist', 'wf_scheduler_hist', ...
    'solve_times', 'time', 'ref_plot', 'controller_name', 'track_name');

%% ===== Helpers =====
function s = ternary(c,a,b); if c; s=a; else; s=b; end; end

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

function y = measure_tracking_output(out_dyn, name)
    nl = lower(string(name));
    if contains(nl,'n_h'); v = sv(out_dyn.eng.Shaft.N_HPC);
    elseif contains(nl,'n_l'); v = sv(out_dyn.eng.Shaft.N_Fan);
    else; v = sv(out_dyn.eng.Shaft.N_HPC); end
    y = v(end);
end

function x = sv(s)
    if isstruct(s)&&isfield(s,'Data'); x = double(s.Data(:));
    elseif isa(s,'timeseries'); x = double(s.Data(:));
    else; x = double(s(:)); end
end
