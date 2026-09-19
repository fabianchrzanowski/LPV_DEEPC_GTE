% Always load params
run('setup_simulation_params.m');
% ========================================================================
% LPV-DeePC
% Changed from IPOPT to OSQP still using CasADi, can use qpoases
% ========================================================================

try; proj = currentProject; catch; proj = openProject('GTE.prj'); end
cd(proj.RootFolder)

import casadi.*

%% ===== Load data =====
training_folder = fullfile(proj.RootFolder,'Params','test_files_fuel_only','test_identification','AGTF30');
models = load(fullfile(proj.RootFolder,'Params','tables','matrices_DE_wf.mat'));
wf_data = load(fullfile(proj.RootFolder,'Params','test_files_fuel_only','test_identification','wf_values.mat'));

offset = models.offset;
sys_eng = models.sys_eng;
wf_values = wf_data.Wf_values(:)';

output_names_lower = lower(string(sys_eng.OutputName));
idx_track = find(contains(lower(string(sys_eng.OutputName)), "n_h"), 1);
if isempty(idx_track); idx_track = 2; end
track_name = sys_eng.OutputName{idx_track};
fprintf('Tracking output #%d: %s (equilibrium = %.1f)\n', idx_track, track_name, z0(idx_track));


u_min = 0.35;
u_max = 2.0;
du_max = 0.15;

q_track      = 15;
r_u          = 1;
r_du         = 1;
lambda_g     = 300;
lambda_sigma = 1e5;

if exist('OVERRIDE_MAX_COLS', 'var')
    max_cols_per_wf = OVERRIDE_MAX_COLS;
else
    max_cols_per_wf = 100;
end
n_sys = 4;

%% ===== Build per-Wf Hankels =====
filelist = dir(fullfile(training_folder,'*.mat'));
[~,idx_sort] = sort({filelist.name});
filelist = filelist(idx_sort);

if exist('OVERRIDE_NUM_HANKELS', 'var') && OVERRIDE_NUM_HANKELS < numel(filelist)
    keep_idx = unique(round(linspace(1, numel(filelist), OVERRIDE_NUM_HANKELS)));
    filelist = filelist(keep_idx);
    hankel_wf = wf_values(keep_idx);
    n_op = numel(filelist);
else
    n_op = min(numel(filelist), numel(wf_values));
    hankel_wf = wf_values(1:n_op);
end

fprintf('Building per-Wf Hankels (%d operating points)...\n', n_op);

hankel_U      = cell(n_op, 1);
hankel_Y      = cell(n_op, 1);
hankel_rank   = zeros(n_op, 1);
hankel_cols   = zeros(n_op, 1);
hankel_u_mean = zeros(n_op, 1);
hankel_y_mean = zeros(n_op, 1);

for k = 1:n_op
    data = load(fullfile(filelist(k).folder, filelist(k).name));
    out_dyn = fetch_out_dyn(data);
    u = extract_input_trace(out_dyn);
    y = measure_tracking_series(out_dyn, track_name);

    n_common = min(numel(u), numel(y));
    u = u(1:n_common);
    y = y(1:n_common);

    if numel(u) < L
        warning('Wf=%.3f: not enough data (%d < %d)', hankel_wf(k), numel(u), L);
        continue
    end

    hankel_u_mean(k) = mean(u);
    hankel_y_mean(k) = mean(y);

    % Build Hankels in DEVIATION space
    du = u - hankel_u_mean(k);
    dy = y - hankel_y_mean(k);
    Hu = block_hankel(du(:), L);
    Hy = block_hankel(dy(:), L);

    % Rank check (Willems PE condition)
    Hu_check  = block_hankel(du(:), L + n_sys);
    r         = svd_rank(Hu_check, 1e-6);
    hankel_rank(k) = r;
    rank_needed    = min(L + n_sys, min(size(Hu_check)));

    % Downsample columns
    nc = size(Hu, 2);
    if nc > max_cols_per_wf
        keep = unique(round(linspace(1, nc, max_cols_per_wf)));
        Hu = Hu(:, keep);
        Hy = Hy(:, keep);
    end

    hankel_U{k}    = Hu;
    hankel_Y{k}    = Hy;
    hankel_cols(k) = size(Hu, 2);

    fprintf('  Wf=%.3f: u_mean=%.3f y_mean=%.0f, %d cols, rank=%d (need %d) %s\n', ...
        hankel_wf(k), hankel_u_mean(k), hankel_y_mean(k), hankel_cols(k), ...
        r, rank_needed, ternary(r >= rank_needed, 'OK', 'LOW'));
end

% Remove empty entries
valid         = ~cellfun(@isempty, hankel_U);
hankel_U      = hankel_U(valid);
hankel_Y      = hankel_Y(valid);
hankel_wf     = hankel_wf(valid);
hankel_rank   = hankel_rank(valid);
hankel_cols   = hankel_cols(valid);
hankel_u_mean = hankel_u_mean(valid);
hankel_y_mean = hankel_y_mean(valid);
n_valid       = sum(valid);

fprintf('\nValid Hankels: %d / %d\n', n_valid, n_op);

%% ===== Pre-build CasADi solver =====
fprintf('Pre-building CasADi solver...\n');
tic;

n_cols_std = max(hankel_cols);

% Decision variables
g_var     = SX.sym('g',     n_cols_std, 1);
U_var     = SX.sym('U',     N_pred,     1);
Y_var     = SX.sym('Y',     N_pred,     1);
sigma_var = SX.sym('sigma', T_ini,      1);

% Parameters
p_u_ini  = SX.sym('u_ini',  T_ini,       1);
p_y_ini  = SX.sym('y_ini',  T_ini,       1);
p_ref    = SX.sym('ref',    N_pred,      1);
p_u_prev = SX.sym('u_prev', 1,           1);
p_u0     = SX.sym('u0',     1,           1);
p_Up     = SX.sym('Up',     T_ini,       n_cols_std);
p_Yp     = SX.sym('Yp',     T_ini,       n_cols_std);
p_Uf     = SX.sym('Uf',     N_pred,      n_cols_std);
p_Yf     = SX.sym('Yf',     N_pred,      n_cols_std);

x = [g_var; U_var; Y_var; sigma_var];
p = [p_u_ini; p_y_ini; p_ref; p_u_prev; p_u0; ...
     p_Up(:); p_Yp(:); p_Uf(:); p_Yf(:)];

% Equality constraints — Hankel predictor
g_eq = [p_Up*g_var - p_u_ini;
        p_Yp*g_var - p_y_ini - sigma_var;
        p_Uf*g_var - U_var;
        p_Yf*g_var - Y_var];

% Inequality constraints — rate limits
g_ineq = SX.zeros(N_pred, 1);
for i = 1:N_pred
    if i == 1
        g_ineq(i) = U_var(i) - p_u_prev;
    else
        g_ineq(i) = U_var(i) - U_var(i-1);
    end
end

% Quadratic objective
obj = lambda_g*(g_var'*g_var) + lambda_sigma*(sigma_var'*sigma_var);
for i = 1:N_pred
    if i == 1; du_i = U_var(i) - p_u_prev;
    else;      du_i = U_var(i) - U_var(i-1); end
    obj = obj + q_track*(Y_var(i) - p_ref(i))^2 ...
              + r_u*(U_var(i) - p_u0)^2 ...
              + r_du*du_i^2;
end

g_all = [g_eq; g_ineq];
n_eq  = length(g_eq);

nlp = struct('x', x, 'f', obj, 'g', g_all, 'p', p);

% --------------------------------------
% SOLVER Settings
% ---------------------------------------
opts = struct;
opts.qpsol                          = 'osqp';
opts.qpsol_options.verbose          = false;    
opts.max_iter                       = 1;        
opts.print_time                     = 0;
opts.verbose                        = false;
opts.print_header                   = false;
opts.print_iteration                = false;

S.solver  = nlpsol('deepc_fast', 'sqpmethod', nlp, opts);
S.lbg     = [zeros(n_eq,1); -du_max*ones(N_pred,1)];
S.ubg     = [zeros(n_eq,1);  du_max*ones(N_pred,1)];
S.n_x     = n_cols_std + N_pred + N_pred + T_ini;
S.n_cols  = n_cols_std;

build_time = toc;
fprintf('Solver built in %.1f s\n', build_time);

%% ===== Setup plant =====

MWS = struct();
MWS.engName = 'AGTF30';
MWS = setup_Controller(MWS);
MWS = setup_AllEng(MWS);

busVars     = load(fullfile(proj.RootFolder,'Params','Eng_Bus.mat'));
busVarNames = fieldnames(busVars);
for ib = 1:numel(busVarNames)
    assignin('base', busVarNames{ib}, busVars.(busVarNames{ib}));
end

MWS.In.Ts   = Ts;
MWS.In.Tsim = Ts;
MWS.In.Alt  = timeseries([0;0],[0 Ts],'Name','Alt');
MWS.In.MN   = timeseries([0;0],[0 Ts],'Name','MN');
MWS.In.dT   = timeseries([0;0],[0 Ts],'Name','dT');
MWS.In.dVBV = timeseries([0;0],[0 Ts],'Name','dVBV');
MWS.In.dNz  = timeseries([0;0],[0 Ts],'Name','dNz');
MWS.In.Wf   = timeseries([u0;u0],[0 Ts],'Name','Wf');
MWS = AGTF30_initial_conditions(MWS);

model = 'AGTF30SysDyn';

%% ===== Initialise histories =====
u_hist         = zeros(N_total, 1);
y_hist         = nan(N_total, 1);
sm_hpc_hist    = nan(N_total, 1);
u_hist(1:T_ini) = u0;
u_prev         = u0;

x0_guess       = zeros(S.n_x, 1);
solve_times    = zeros(N_sim, 1);
sched_idx      = zeros(N_sim, 1);
wf_scheduler_hist = zeros(N_sim, 1);
cost_hist      = zeros(N_sim, 1);

%% ===== Warm-up =====
fprintf('Warm-up (%d steps)...\n', T_ini);
for t = 1:T_ini

    t_step = (0:t)'*Ts;
    u_step = [u0; u_hist(1:t)];
    MWS.In.Tsim = t_step(end);
    MWS.In.Alt  = timeseries([0;0],[0 t_step(end)],'Name','Alt');
    MWS.In.MN   = timeseries([0;0],[0 t_step(end)],'Name','MN');
    MWS.In.dT   = timeseries([0;0],[0 t_step(end)],'Name','dT');
    MWS.In.dVBV = timeseries([0;0],[0 t_step(end)],'Name','dVBV');
    MWS.In.dNz  = timeseries([0;0],[0 t_step(end)],'Name','dNz');
    MWS.In.Wf   = timeseries(u_step, t_step, 'Name','Wf');
    simIn = Simulink.SimulationInput(model);
    simIn = simIn.setVariable('MWS', MWS);
    simIn = simIn.setModelParameter('StartTime','0','StopTime',num2str(t_step(end)),'ReturnWorkspaceOutputs','on');
    simOut = sim(simIn);
    out_dyn_k  = fetch_out_dyn(simOut);
    y_hist(t)  = measure_tracking_output(out_dyn_k, track_name);
    sm_v = sv(out_dyn_k.eng.SM.SMHPC); sm_hpc_hist(t) = sm_v(end);

end

%% ===== Closed-loop =====
wf_filtered = u0;
fprintf('LPV-DeePC closed-loop (%d steps)...\n', N_sim);

for t = T_ini+1:N_total
    si = t - T_ini;

    % ---- STEP 1: LPV scheduling ----
    % Pick the Hankel closest to the current fuel flow.
    % A low-pass filter smooths the scheduling to avoid chattering.
    wf_for_sched = max(min(u_prev, max(hankel_wf)), min(hankel_wf));
    if si == 1
        wf_filtered = wf_for_sched;
    else
        alpha       = 0.3;
        wf_filtered = alpha * wf_for_sched + (1 - alpha) * wf_filtered;
    end
    [~, k_sel]          = min(abs(hankel_wf - wf_filtered));
    sched_idx(si)        = k_sel;
    wf_scheduler_hist(si) = hankel_wf(k_sel);

    % ---- STEP 2: Convert signals to deviation space ----
    % Subtract the operating-point means so the Hankel data is consistent.
    u_mean_k = hankel_u_mean(k_sel);
    y_mean_k = hankel_y_mean(k_sel);

    du_ini   = u_hist(t-T_ini:t-1) - u_mean_k;
    dy_ini   = y_hist(t-T_ini:t-1) - y_mean_k;
    dref_seg = ref(t:t+N_pred-1)   - y_mean_k;
    du_prev  = u_prev - u_mean_k;
    du0      = u0     - u_mean_k;

    % ---- STEP 3: Load the Hankel matrices for the selected operating point ----
    % Pad with zeros if this Hankel has fewer columns than the solver expects.
    nc_sel = hankel_cols(k_sel);
    Hu_sel = hankel_U{k_sel};
    Hy_sel = hankel_Y{k_sel};
    if nc_sel < S.n_cols
        Hu_sel = [Hu_sel, zeros(L, S.n_cols - nc_sel)];
        Hy_sel = [Hy_sel, zeros(L, S.n_cols - nc_sel)];
    end
    Up_sel = Hu_sel(1:T_ini,   :);
    Uf_sel = Hu_sel(T_ini+1:end, :);
    Yp_sel = Hy_sel(1:T_ini,   :);
    Yf_sel = Hy_sel(T_ini+1:end, :);

    % Pack all parameters into a single vector for CasADi
    p_val = [du_ini; dy_ini; dref_seg; du_prev; du0; ...
             Up_sel(:); Yp_sel(:); Uf_sel(:); Yf_sel(:)];

    % ---- STEP 4: Set variable bounds ----
    % Pin zero-padded g elements to zero so regulariser only acts on real data.
    % Input bounds are shifted into deviation space.
    g_lb = -inf(S.n_cols, 1);
    g_ub =  inf(S.n_cols, 1);
    if nc_sel < S.n_cols
        g_lb(nc_sel+1:end) = 0;
        g_ub(nc_sel+1:end) = 0;
    end
    lbx = [g_lb; (u_min - u_mean_k)*ones(N_pred,1); -inf(N_pred,1); -inf(T_ini,1)];
    ubx = [g_ub; (u_max - u_mean_k)*ones(N_pred,1);  inf(N_pred,1);  inf(T_ini,1)];

    % ---- STEP 5: Solve the DeePC QP (warm-started from previous solution) ----
    tic;
    sol = S.solver('x0',  x0_guess, ...
                   'p',   p_val,    ...
                   'lbx', lbx,      ...
                   'ubx', ubx,      ...
                   'lbg', S.lbg,    ...
                   'ubg', S.ubg);
    solve_times(si) = toc;

    x_sol      = full(sol.x);
    cost_hist(si) = full(sol.f);
    x0_guess   = x_sol;   % warm-start the next solve

    % ---- STEP 6: Extract the first optimal input and convert back to absolute ----
    du_now = x_sol(S.n_cols + 1);
    u_now  = du_now + u_mean_k;

    if isnan(u_now)
        warning('Solver returned NaN at step %d — fallback to u_prev.', t);
        u_now = u_prev;
    end

    u_now       = max(u_min, min(u_max, u_now));  % clamp to actuator limits
    u_hist(t)   = u_now;

    if mod(t, 5) == 0 || t == 1
        fprintf('Step %d / %d | u_now: %.4f pps | wf_sched: %.4f\n', ...
            t, N_total, u_now, wf_filtered);
    end

    % ---- STEP 7: Simulate the plant for one step ----

    % Run the full nonlinear Simulink model
    t_step = (0:t)'*Ts;
    u_step = [u0; u_hist(1:t)];
    MWS.In.Tsim = t_step(end);
    MWS.In.Alt  = timeseries([0;0],[0 t_step(end)],'Name','Alt');
    MWS.In.MN   = timeseries([0;0],[0 t_step(end)],'Name','MN');
    MWS.In.dT   = timeseries([0;0],[0 t_step(end)],'Name','dT');
    MWS.In.dVBV = timeseries([0;0],[0 t_step(end)],'Name','dVBV');
    MWS.In.dNz  = timeseries([0;0],[0 t_step(end)],'Name','dNz');
    MWS.In.Wf   = timeseries(u_step, t_step, 'Name','Wf');
    simIn = Simulink.SimulationInput(model);
    simIn = simIn.setVariable('MWS', MWS);
    simIn = simIn.setModelParameter('StartTime','0','StopTime',num2str(t_step(end)),'ReturnWorkspaceOutputs','on');
    simOut    = sim(simIn);
    out_dyn_k = fetch_out_dyn(simOut);
    y_hist(t) = measure_tracking_output(out_dyn_k, track_name);
    sm_v = sv(out_dyn_k.eng.SM.SMHPC); sm_hpc_hist(t) = sm_v(end);


    % Update previous input for the next iteration
    u_prev = u_now;

    if mod(si, 20) == 0
        fprintf('  Step %d/%d  Wf=%.3f -> Hankel #%d (Wf=%.3f)  solve=%.1fms\n', ...
            si, N_sim, u_prev, k_sel, hankel_wf(k_sel), solve_times(si)*1000);
    end
end

%% ===== Save results =====
time     = (0:N_total-1)'*Ts;
ref_plot = ref(1:N_total);
scen_tag = lower(SCENARIO);

if exist('OVERRIDE_MAT_FILENAME', 'var') && ~isempty(OVERRIDE_MAT_FILENAME)
    mat_filename = OVERRIDE_MAT_FILENAME;
else
    mat_filename = sprintf('Results/results_lpv_deepc_osqp_%s_%dcol_%dhank_lg%.0e_ls%.0e.mat', ...
        scen_tag, max_cols_per_wf, n_op, lambda_g, lambda_sigma);
end

controller_name = sprintf('LPV DeePC (OSQP) %s (%dcol, %dhank, \\lambda_g=%.0e, \\lambda_\\sigma=%.0e)', ...
    SCENARIO, max_cols_per_wf, n_op, lambda_g, lambda_sigma);

save(mat_filename, 'y_hist', 'u_hist', 'wf_scheduler_hist', 'cost_hist', ...
     'sm_hpc_hist', 'solve_times', 'time', 'ref_plot', 'controller_name', 'track_name');
fprintf('Saved results to %s\n', mat_filename);

%% ===== Plot results =====
control_time = time(T_ini+1:end);

figure('Name','LPV-DeePC closed-loop on nonlinear plant')
subplot(3,1,1)
plot(control_time, ref_plot(T_ini+1:N_total), '--', 'LineWidth', 1.2)
hold on
plot(control_time, y_hist(T_ini+1:end), 'LineWidth', 1.4)
xline(T_ini*Ts, 'r--', 'DeePC starts', 'LabelHorizontalAlignment', 'left')
grid on
xlabel('Time [s]')
ylabel(track_name, 'Interpreter', 'tex')
legend('Reference', 'AGTF30 output', 'Location', 'best')
title('LPV-DeePC tracking (OSQP)')

subplot(3,1,2)
stairs(control_time, u_hist(T_ini+1:end), 'LineWidth', 1.4)
xline(T_ini*Ts, 'r--')
grid on
xlabel('Time [s]')
ylabel('W_f [pps]')
title('Applied fuel command')

subplot(3,1,3)
stairs(1:N_sim, wf_scheduler_hist, 'LineWidth', 1.4)
grid on
xlabel('MPC step')
ylabel('W_f scheduled')
title('LPV Scheduling variable (\rho)')

%% ===== Print summary =====
y_ctrl   = y_hist(T_ini+1:end);
ref_ctrl = ref_plot(T_ini+1:N_total);
rmse     = sqrt(mean((y_ctrl - ref_ctrl).^2));

fprintf('\n===== LPV-DeePC OSQP Results =====\n');
fprintf('Tracking output : %s\n',   track_name);
fprintf('RMSE            : %.1f\n', rmse);
fprintf('Mean solve time : %.2f ms\n', mean(solve_times)*1000);
fprintf('Max  solve time : %.2f ms\n', max(solve_times)*1000);
fprintf('Std  solve time : %.2f ms\n', std(solve_times)*1000);

% Solve time histogram
figure('Name', 'Solve time distribution')
histogram(solve_times*1000, 30)
xlabel('Solve time (ms)')
ylabel('Count')
title('QP solve time distribution (OSQP)')
grid on
xline(mean(solve_times)*1000, 'r--', sprintf('Mean %.2f ms', mean(solve_times)*1000))
xline(max(solve_times)*1000,  'k--', sprintf('Max  %.2f ms', max(solve_times)*1000))

%% ===== Helper functions =====
function s = ternary(cond, a, b)
    if cond; s = a; else; s = b; end
end

function H = block_hankel(sig, L)
    sig = sig(:); nc = numel(sig) - L + 1;
    H = zeros(L, nc);
    for i = 1:L; H(i,:) = sig(i:i+nc-1).'; end
end

function r = svd_rank(A, tol)
    s = svd(A);
    if isempty(s) || s(1) == 0; r = 0; return; end
    r = sum(s > tol * s(1));
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
        if isfield(out_dyn.cntrl,'Wf');    u = sv(out_dyn.cntrl.Wf);    return; end
    end
    error('No Wf trace found');
end

function y = measure_tracking_series(out_dyn, name)
    nl = lower(string(name));
    if contains(nl,'n_h');                        y = sv(out_dyn.eng.Shaft.N_HPC);
    elseif contains(nl,'n_l');                    y = sv(out_dyn.eng.Shaft.N_Fan);
    elseif contains(nl,'thrust')||contains(nl,'fnet'); y = sv(out_dyn.eng.Perf.Fnet);
    else;                                         y = sv(out_dyn.eng.Shaft.N_HPC);
    end
end

function y = measure_tracking_output(out_dyn, name)
    v = measure_tracking_series(out_dyn, name); y = v(end);
end

function x = sv(s)
    if isstruct(s) && isfield(s,'Data'); x = double(s.Data(:));
    elseif isa(s,'timeseries');          x = double(s.Data(:));
    else;                                x = double(s(:));
    end
end