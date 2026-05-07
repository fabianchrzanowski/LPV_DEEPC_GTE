% Always load latest params
run('setup_simulation_params.m');
% ========================================================================
% LPV-DeePC: Gain-Scheduled Data-Enabled Predictive Control
%
% Instead of one big Hankel from all operating points, build a separate
% Hankel at each Wf. At runtime, select the Hankel closest to the
% current operating point — same scheduling idea as LPV-MPC.
%
% Uses CasADi + IPOPT. One solver pre-built per operating point.
% ========================================================================

proj = currentProject;
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

% Variables now sourced from setup_simulation_params.m

u_min = 0.35;
u_max = 2.0;
du_max = 0.2;

q_track = 1;
r_u = 1;
r_du = 25;
lambda_g = 50;        % Heavily penalize large g to prevent overfitting to DC offsets
lambda_sigma = 10;    % Allow sigma to easily absorb DC offset mismatches

max_cols_per_wf = 3000;
n_sys = 4;

% The trajectory ref is now sourced from setup_simulation_params.m

%% ===== Build per-Wf Hankels =====
filelist = dir(fullfile(training_folder,'*.mat'));
[~,idx_sort] = sort({filelist.name});
filelist = filelist(idx_sort);
n_op = min(numel(filelist), numel(wf_values));

fprintf('Building per-Wf Hankels (%d operating points)...\n', n_op);

hankel_U = cell(n_op, 1);
hankel_Y = cell(n_op, 1);
hankel_rank = zeros(n_op, 1);
hankel_cols = zeros(n_op, 1);
hankel_wf = wf_values(1:n_op);
hankel_u_mean = zeros(n_op, 1);  % mean input at each operating point
hankel_y_mean = zeros(n_op, 1);  % mean output at each operating point

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

    % Store means for deviation conversion
    hankel_u_mean(k) = mean(u);
    hankel_y_mean(k) = mean(y);

    % Build Hankels in DEVIATION space (like LPV offset handling)
    du = u - hankel_u_mean(k);
    dy = y - hankel_y_mean(k);
    Hu = block_hankel(du(:), L);
    Hy = block_hankel(dy(:), L);

    % Rank check
    r = rank([Hu; Hy]);
    hankel_rank(k) = r;
    rank_needed = L + n_sys;

    % Downsample
    nc = size(Hu, 2);
    if nc > max_cols_per_wf
        keep = unique(round(linspace(1, nc, max_cols_per_wf)));
        Hu = Hu(:, keep);
        Hy = Hy(:, keep);
    end

    hankel_U{k} = Hu;
    hankel_Y{k} = Hy;
    hankel_cols(k) = size(Hu, 2);

    fprintf('  Wf=%.3f: u_mean=%.3f y_mean=%.0f, %d cols, rank=%d (need %d) %s\n', ...
        hankel_wf(k), hankel_u_mean(k), hankel_y_mean(k), hankel_cols(k), ...
        r, rank_needed, ternary(r >= rank_needed, 'OK', 'LOW'));
end

% Remove empty entries
valid = ~cellfun(@isempty, hankel_U);
hankel_U = hankel_U(valid);
hankel_Y = hankel_Y(valid);
hankel_wf = hankel_wf(valid);
hankel_rank = hankel_rank(valid);
hankel_cols = hankel_cols(valid);
hankel_u_mean = hankel_u_mean(valid);
hankel_y_mean = hankel_y_mean(valid);
n_valid = sum(valid);

fprintf('\nValid Hankels: %d / %d\n', n_valid, n_op);

%% ===== Pre-build one CasADi solver per operating point =====
fprintf('Pre-building %d CasADi solvers...\n', n_valid);
tic;

% Standardise column count (pad smaller ones with zeros)
n_cols_std = max(hankel_cols);
solvers = cell(n_valid, 1);

for k = 1:n_valid
    nc = hankel_cols(k);
    Hu = hankel_U{k};
    Hy = hankel_Y{k};

    % Pad to standard size if needed
    if nc < n_cols_std
        Hu = [Hu, zeros(L, n_cols_std - nc)];
        Hy = [Hy, zeros(L, n_cols_std - nc)];
    end

    Up = Hu(1:T_ini, :);
    Uf = Hu(T_ini+1:end, :);
    Yp = Hy(1:T_ini, :);
    Yf = Hy(T_ini+1:end, :);

    % Decision variables
    g_var     = SX.sym('g', n_cols_std, 1);
    U_var     = SX.sym('U', N_pred, 1);
    Y_var     = SX.sym('Y', N_pred, 1);
    sigma_var = SX.sym('sigma', T_ini, 1);

    % Parameters
    p_u_ini  = SX.sym('u_ini', T_ini, 1);
    p_y_ini  = SX.sym('y_ini', T_ini, 1);
    p_ref    = SX.sym('ref', N_pred, 1);
    p_u_prev = SX.sym('u_prev', 1, 1);
    p_u0     = SX.sym('u0', 1, 1);

    x = [g_var; U_var; Y_var; sigma_var];
    p = [p_u_ini; p_y_ini; p_ref; p_u_prev; p_u0];

    % Constraints
    g_eq = [Up*g_var - p_u_ini;
            Yp*g_var - p_y_ini - sigma_var;
            Uf*g_var - U_var;
            Yf*g_var - Y_var];

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

    solvers{k} = struct();
    solvers{k}.solver = nlpsol(sprintf('deepc_%d',k), 'ipopt', nlp, opts);
    solvers{k}.lbg = [zeros(n_eq,1); -du_max*ones(N_pred,1)];
    solvers{k}.ubg = [zeros(n_eq,1);  du_max*ones(N_pred,1)];
    solvers{k}.lbx = [-inf(n_cols_std,1); (u_min - hankel_u_mean(k))*ones(N_pred,1); -inf(N_pred,1); -inf(T_ini,1)];
    solvers{k}.ubx = [ inf(n_cols_std,1); (u_max - hankel_u_mean(k))*ones(N_pred,1);  inf(N_pred,1);  inf(T_ini,1)];
    solvers{k}.n_x = n_cols_std + N_pred + N_pred + T_ini;
    solvers{k}.n_cols = n_cols_std;

    if mod(k,5)==0; fprintf('  %d / %d built\n', k, n_valid); end
end

build_time = toc;
fprintf('All %d solvers built in %.1f s\n', n_valid, build_time);

if ~FAST_SURROGATE_MODE
    % Build baseline MWS once
    MWS = struct(); % Prevent old arrays from causing dimension mismatch on 2nd run
    MWS.engName = 'AGTF30';
    MWS = setup_Controller(MWS);
    MWS = setup_AllEng(MWS);

    busVars = load(fullfile(proj.RootFolder,'Params','Eng_Bus.mat'));
    busVarNames = fieldnames(busVars);
    for ib = 1:numel(busVarNames)
        assignin('base',busVarNames{ib},busVars.(busVarNames{ib}));
    end

    MWS.In.Ts = Ts;
    MWS.In.Tsim = Ts;
    MWS.In.Alt = timeseries([0;0],[0 Ts],'Name','Alt');
    MWS.In.MN  = timeseries([0;0],[0 Ts],'Name','MN');
    MWS.In.dT  = timeseries([0;0],[0 Ts],'Name','dT');
    MWS.In.dVBV = timeseries([0;0],[0 Ts],'Name','dVBV');
    MWS.In.dNz  = timeseries([0;0],[0 Ts],'Name','dNz');
    MWS.In.Wf   = timeseries([u0;u0],[0 Ts],'Name','Wf');
    MWS = AGTF30_initial_conditions(MWS);
end

model = 'AGTF30SysDyn';

%% ===== Closed-loop =====
u_hist = zeros(N_total,1);
y_hist = nan(N_total,1);
u_hist(1:T_ini) = u0;
u_prev = u0;
if FAST_SURROGATE_MODE
    x_plant_fast = zeros(size(models.A_dt,1),1);
end
x0_guess = zeros(solvers{1}.n_x, 1);
solve_times = zeros(N_sim, 1);
sched_idx = zeros(N_sim, 1);
wf_scheduler_hist = zeros(N_sim, 1);

% Warm-up
fprintf('Warm-up (%d steps)...\n', T_ini);
for t = 1:T_ini
    if FAST_SURROGATE_MODE
        [y_hist(t), x_plant_fast] = LPV_step_fast(x_plant_fast, u0, u0, models);
    else
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
        simIn = simIn.setVariable('MWS',MWS);
        simIn = simIn.setModelParameter('StartTime','0','StopTime',num2str(t_step(end)),'ReturnWorkspaceOutputs','on');
        simOut = sim(simIn);
        y_hist(t) = measure_tracking_output(fetch_out_dyn(simOut), track_name);
    end
end

% Control loop
wf_filtered = u0;  % start scheduling from equilibrium
fprintf('LPV-DeePC closed-loop (%d steps)...\n', N_sim);
for t = T_ini+1:N_total
    si = t - T_ini;

    % Schedule: pick closest Wf operating point
    wf_for_sched = max(min(u_prev, max(hankel_wf)), min(hankel_wf));  % clamp to Hankel range
    if si == 1
        wf_filtered = wf_for_sched;
    else
        alpha = 0.3;  % smoothing factor (0 = ignore new, 1 = no smoothing)
        wf_filtered = alpha * wf_for_sched + (1 - alpha) * wf_filtered;
    end
    [~, k_sel] = min(abs(hankel_wf - wf_filtered));
    sched_idx(si) = k_sel;
    wf_scheduler_hist(si) = hankel_wf(k_sel);

    % Convert to deviation space for the selected operating point
    u_mean_k = hankel_u_mean(k_sel);
    y_mean_k = hankel_y_mean(k_sel);

    du_ini = u_hist(t-T_ini:t-1) - u_mean_k;
    dy_ini = y_hist(t-T_ini:t-1) - y_mean_k;
    dref_seg = ref(t:t+N_pred-1) - y_mean_k;
    du_prev = u_prev - u_mean_k;
    du0 = u0 - u_mean_k;

    p_val = [du_ini; dy_ini; dref_seg; du_prev; du0];

    S = solvers{k_sel};
    
    % FIX FOR MATLAB CRASH: 
    % When we switch operating points, the new Hankel matrix might have a different 
    % number of columns. If we pass the old x0_guess (which is the wrong length) 
    % into CasADi's compiled C++ MEX function, it causes a Segmentation Fault and 
    % instantly kills MATLAB. We must resize it to match S.n_x!
    if length(x0_guess) ~= S.n_x
        if length(x0_guess) > S.n_x
            x0_guess = x0_guess(1:S.n_x);
        else
            x0_guess = [x0_guess; zeros(S.n_x - length(x0_guess), 1)];
        end
    end

    tic;
    sol = S.solver('x0', x0_guess, 'p', p_val, ...
                   'lbx', S.lbx, 'ubx', S.ubx, ...
                   'lbg', S.lbg, 'ubg', S.ubg);
    solve_times(si) = toc;

    x_sol = full(sol.x);
    x0_guess = x_sol;

    % Convert back from deviation space to absolute
    du_now = x_sol(S.n_cols + 1);
    u_now = du_now + u_mean_k;
    
    if isnan(u_now)
        warning('CasADi solver returned NaN at step %d! Fallback to u_prev.', t);
        u_now = u_prev;
    end
    
    u_now = max(u_min, min(u_max, u_now));  % clamp
    u_hist(t) = u_now;
    
    if mod(t, 5) == 0 || t == 1
        fprintf('Step %d / %d | u_now: %.4f pps | wf_sched: %.4f\n', t, N_total, u_now, wf_filtered);
    end

    % Simulate
    if FAST_SURROGATE_MODE
        [y_hist(t), x_plant_fast] = LPV_step_fast(x_plant_fast, u_now, wf_filtered, models);
    else
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
        simIn = simIn.setVariable('MWS',MWS);
        simIn = simIn.setModelParameter('StartTime','0','StopTime',num2str(t_step(end)),'ReturnWorkspaceOutputs','on');
        simOut = sim(simIn);
        y_hist(t) = measure_tracking_output(fetch_out_dyn(simOut), track_name);
    end

    u_prev = u_now;
    if mod(si,20)==0
        fprintf('  Step %d/%d  Wf=%.3f -> Hankel #%d (Wf=%.3f)  solve=%.1fms\n', ...
            si, N_sim, u_prev, k_sel, hankel_wf(k_sel), solve_times(si)*1000);
    end
end

time = (0:N_total-1)'*Ts;
ref_plot = ref(1:N_total);
controller_name = 'LPV DeePC';
save('Results/results_lpv_deepc.mat', 'y_hist', 'u_hist', 'wf_scheduler_hist', 'solve_times', 'time', 'ref_plot', 'controller_name', 'track_name');
fprintf('Saved results to Results/results_lpv_deepc.mat\n');

%% ===== Results =====
if ~FAST_SURROGATE_MODE
    control_time = time(T_ini+1:end);

    figure('Name','LPV-DeePC closed-loop on nonlinear plant')
    subplot(3,1,1)
    plot(control_time,ref_plot(T_ini+1:N_total),'--','LineWidth',1.2)
    hold on
    plot(control_time,y_hist(T_ini+1:end),'LineWidth',1.4)
    xline(T_ini*Ts, 'r--', 'DeePC starts', 'LabelHorizontalAlignment','left')
    grid on
    xlabel('Time [s]')
    ylabel(track_name,'Interpreter','tex')
    legend('Reference','AGTF30 measured output','Location','best')
    title('LPV-DeePC tracking')

    subplot(3,1,2)
    plot(control_time,u_hist(T_ini+1:end),'LineWidth',1.4)
    xline(T_ini*Ts, 'r--')
    grid on
    xlabel('Time [s]')
    ylabel('W_f [pps]')
    title('Applied fuel command')

    subplot(3,1,3)
    plot(1:N_sim, wf_scheduler_hist, 'LineWidth', 1.4)
    grid on
    xlabel('MPC step')
    ylabel('W_f scheduled')
    title('LPV Scheduling variable (\rho)')

    % RMSE
    y_ctrl = y_hist(T_ini+1:end);
    ref_ctrl = ref_plot(T_ini+1:N_total);
    rmse = sqrt(mean((y_ctrl - ref_ctrl).^2));
    fprintf('\n===== LPV-DeePC CasADi Results =====\n');
    fprintf('Tracking output: %s\n', track_name);
    fprintf('RMSE: %.1f\n', rmse);
    fprintf('Mean solve time: %.1f ms\n', mean(solve_times)*1000);
    fprintf('Max solve time:  %.1f ms\n', max(solve_times)*1000);
end

%% ===== Helpers =====
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

function y = measure_tracking_output(out_dyn, name)
    v = measure_tracking_series(out_dyn, name); y = v(end);
end

function x = sv(s)
    if isstruct(s)&&isfield(s,'Data'); x = double(s.Data(:));
    elseif isa(s,'timeseries'); x = double(s.Data(:));
    else; x = double(s(:)); end
end
