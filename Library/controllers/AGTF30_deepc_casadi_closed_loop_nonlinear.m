% Always load latest params
run('setup_simulation_params.m');
% ========================================================================
% DeePC (Data-Enabled Predictive Control) — CasADi version
%
% Same as AGTF30_deepc_yalmip_closed_loop_nonlinear.m but using CasADi
% for the QP. Key speedup: the problem is built ONCE symbolically and
% the solver is compiled. Each MPC step only updates the parameter vector
% (u_ini, y_ini, ref, u_prev) — no symbolic rebuilding.
% ========================================================================

proj = currentProject;
cd(proj.RootFolder)

import casadi.*

% Variables now sourced from setup_simulation_params.m

training_folder = fullfile(proj.RootFolder,'Params','test_files_fuel_only','test_identification','AGTF30');
if ~isfolder(training_folder)
    error('DeePC training folder not found: %s',training_folder);
end

idx_track = find(contains(lower(string(models.sys_eng.OutputName)), "n_h"), 1);
if isempty(idx_track); idx_track = 2; end
track_name = models.sys_eng.OutputName{idx_track};

% Extract the specific Cruise state space model for DeePC tracking point
wf_frozen = 1.2;
[~, k_frozen] = min(abs(wf_values - wf_frozen));
u0_frozen = offset.u(:, k_frozen);
z0_frozen = offset.z(:, k_frozen);

A_all = models.A_dt;
B_all = models.B_dt;
C_all = models.C_dt;
D_all = models.D_dt;

u_min = 0.2;
u_max = 2.2;
du_max = 0.10;           % was 0.04 — too tight, controller couldn't move

q_track = 10;            % was 1  — boost tracking importance
r_u = 0.01;              % was 0.05
r_du = 2;                % was 25 — was 25x heavier than tracking → flat
lambda_g = 1e-6;         % was 1e-4 — was suppressing g → zero predictions
lambda_sigma = 1e4;

% The trajectory ref is now sourced from setup_simulation_params.m
%% ===== Build Hankel matrices (SINGLE operating point) =====
% DeePC (Willems' fundamental lemma) requires data from a SINGLE LTI
% system. Using data from all 32 Wf points mixes different dynamics,
% causing the g-vector to average everything out → flat response.
filelist_train = dir(fullfile(training_folder,'*.mat'));
[~,idx_sort] = sort({filelist_train.name});
filelist_train = filelist_train(idx_sort);
[~, k_file] = min(abs(wf_values - wf_frozen));
k_file = min(k_file, numel(filelist_train));
fprintf('Building Hankel from file #%d (Wf ~ %.3f pps) only\n', k_file, wf_values(k_file));

data_train = load(fullfile(filelist_train(k_file).folder, filelist_train(k_file).name));
out_train = fetch_out_dyn(data_train);
u_train = extract_input_trace(out_train);
y_train = measure_tracking_series(out_train, track_name);
n_common = min(numel(u_train), numel(y_train));
u_train = u_train(1:n_common) - u0_frozen;
y_train = y_train(1:n_common) - z0_frozen(idx_track);

L = T_ini + N_pred;
U_hankel_full = block_hankel(u_train(:), L);
Y_hankel_full = block_hankel(y_train(:), L);
fprintf('Raw Hankel: L=%d, columns=%d\n', size(U_hankel_full,1), size(U_hankel_full,2));

% Persistent excitation check
n_sys = 4;
H_stacked = [U_hankel_full; Y_hankel_full];
r = rank(H_stacked);
rank_needed = L + n_sys;
if r >= rank_needed
    fprintf('Rank = %d (need >= %d) OK\n', r, rank_needed);
else
    warning('Rank %d < %d: data may NOT be persistently exciting!', r, rank_needed);
end
clear H_stacked;

% --- Downsample to a manageable size ---
% g has n_cols decision variables — more columns = slower solver.
% Guidelines: 1000-3000 is a good range. More helps if solver is fast.
max_hankel_cols = 3000;
n_cols_raw = size(U_hankel_full, 2);
[U_hankel, Y_hankel] = reduce_hankel_columns(U_hankel_full, Y_hankel_full, max_hankel_cols);
clear U_hankel_full Y_hankel_full;

U_p = U_hankel(1:T_ini,:);
U_f = U_hankel(T_ini+1:end,:);
Y_p = Y_hankel(1:T_ini,:);
Y_f = Y_hankel(T_ini+1:end,:);

n_cols = size(U_hankel,2);
if n_cols < 1
    error('Not enough training data to build a DeePC Hankel matrix.');
end

fprintf('Downsampled Hankel: L=%d, columns=%d (from %d raw)\n', ...
    size(U_hankel,1), n_cols, n_cols_raw);

%% ===== Build CasADi QP (ONCE) =====
fprintf('Building CasADi QP...\n');
tic;

% Decision variables
g_var     = SX.sym('g', n_cols, 1);
U_var     = SX.sym('U', N_pred, 1);
Y_var     = SX.sym('Y', N_pred, 1);
sigma_var = SX.sym('sigma', T_ini, 1);

% Parameters (updated each MPC step)
p_du_ini  = SX.sym('du_ini', T_ini, 1);
p_dy_ini  = SX.sym('dy_ini', T_ini, 1);
p_dref    = SX.sym('dref', N_pred, 1);
p_du_prev = SX.sym('du_prev', 1, 1);
p_du0     = SX.sym('du0', 1, 1);

% Stack all decision variables
x = [g_var; U_var; Y_var; sigma_var];

% Equality constraints: Hankel data consistency
eq1 = U_p * g_var - p_du_ini;               % U_p * g == du_ini
eq2 = Y_p * g_var - p_dy_ini - sigma_var;    % Y_p * g == dy_ini + sigma
eq3 = U_f * g_var - U_var;                  % U_f * g == dU
eq4 = Y_f * g_var - Y_var;                  % Y_f * g == dY

g_eq = [eq1; eq2; eq3; eq4];

% Inequality constraints: input bounds and rate limits
g_ineq = [];
for i = 1:N_pred
    if i == 1
        du_i = U_var(i) - p_du_prev;
    else
        du_i = U_var(i) - U_var(i-1);
    end
    g_ineq = [g_ineq; du_i];  % -du_max <= du_i <= du_max
end

% Objective
obj = lambda_g * (g_var' * g_var) + lambda_sigma * (sigma_var' * sigma_var);
for i = 1:N_pred
    if i == 1
        du_i = U_var(i) - p_du_prev;
    else
        du_i = U_var(i) - U_var(i-1);
    end
    obj = obj + q_track * (Y_var(i) - p_dref(i))^2 ...
              + r_u * (U_var(i) - p_du0)^2 ...
              + r_du * du_i^2;
end

% Combine all constraints
g_all = [g_eq; g_ineq];

% Bounds
n_eq = length(g_eq);
n_ineq = N_pred;

lbg = [zeros(n_eq, 1); -du_max * ones(n_ineq, 1)];
ubg = [zeros(n_eq, 1);  du_max * ones(n_ineq, 1)];

% Variable bounds
lbx = [-inf(n_cols, 1);      % g: unbounded
       (u_min - u0_frozen) * ones(N_pred, 1);  % U: bounded
       -inf(N_pred, 1);       % Y: unbounded
       -inf(T_ini, 1)];       % sigma: unbounded

ubx = [inf(n_cols, 1);
       (u_max - u0_frozen) * ones(N_pred, 1);
       inf(N_pred, 1);
       inf(T_ini, 1)];

% Parameters vector
p = [p_du_ini; p_dy_ini; p_dref; p_du_prev; p_du0];

% Create solver
nlp = struct('x', x, 'f', obj, 'g', g_all, 'p', p);
opts = struct;
opts.ipopt.print_level = 0;
opts.print_time = 0;
opts.ipopt.max_iter = 200;
opts.ipopt.tol = 1e-6;

solver = nlpsol('deepc_solver', 'ipopt', nlp, opts);

build_time = toc;
fprintf('CasADi QP built in %.2f s (this is one-time cost)\n', build_time);

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

%% ===== Warm-up phase (T_ini steps, no control) =====
u_hist = zeros(N_total,1);
y_hist = nan(N_total,1);
u_hist(1:T_ini) = u0;
u_prev = u0;
if FAST_SURROGATE_MODE
    x_plant_fast = zeros(size(A_all,1),1);
end

% Initial guess for the solver (reused across steps)
x0_guess = zeros(n_cols + N_pred + N_pred + T_ini, 1);

fprintf('Running warm-up phase (%d steps)...\n', T_ini);
for t = 1:T_ini
    if FAST_SURROGATE_MODE
        [~, k_plant] = min(abs(wf_values - u_prev));
        A_p = A_all(:,:,k_plant); B_p = B_all(:,:,k_plant);
        C_p = C_all(:,:,k_plant); D_p = D_all(:,:,k_plant);
        u0_p = offset.u(:,k_plant); z0_p = offset.z(:,k_plant);
        
        y_meas = z0_p(idx_track) + C_p(idx_track,:)*x_plant_fast + D_p(idx_track,:)*(u_hist(t) - u0_p);
        x_plant_fast = A_p*x_plant_fast + B_p*(u_hist(t) - u0_p);
        y_hist(t) = y_meas;
    else
        t_step = (0:t)'*Ts;
        u_step = [u0; u_hist(1:t)];
        MWS.In.Tsim = t_step(end);
        MWS.In.Alt = timeseries([0;0],[0 t_step(end)],'Name','Alt');
        MWS.In.MN  = timeseries([0;0],[0 t_step(end)],'Name','MN');
        MWS.In.dT  = timeseries([0;0],[0 t_step(end)],'Name','dT');
        MWS.In.dVBV = timeseries([0;0],[0 t_step(end)],'Name','dVBV');
        MWS.In.dNz  = timeseries([0;0],[0 t_step(end)],'Name','dNz');
        MWS.In.Wf = timeseries(u_step,t_step,'Name','Wf');

        simIn = Simulink.SimulationInput(model);
        simIn = simIn.setVariable('MWS',MWS);
        simIn = simIn.setModelParameter('StartTime','0','StopTime',num2str(t_step(end)),...
            'ReturnWorkspaceOutputs','on');

        simOut = sim(simIn);
        out_step = fetch_out_dyn(simOut);
        y_hist(t) = measure_tracking_output(out_step, track_name);
    end
end

%% ===== Closed-loop DeePC =====
solve_times = zeros(N_sim, 1);

fprintf('Starting CasADi DeePC closed-loop (%d steps)...\n', N_sim);
for t = T_ini+1:N_total
    step_idx = t - T_ini;

    u_ini = u_hist(t-T_ini:t-1);
    y_ini = y_hist(t-T_ini:t-1);
    ref_seg = ref(t:t+N_pred-1);

    % Convert to deviation space
    du_ini = u_ini - u0_frozen;
    dy_ini = y_ini - z0_frozen(idx_track);
    dref_seg = ref_seg - z0_frozen(idx_track);
    du_prev = u_prev - u0_frozen;
    du0 = u0_frozen - u0_frozen; % which is 0, but keeping semantics

    % Update parameters
    p_val = [du_ini; dy_ini; dref_seg; du_prev; du0];

    % Solve
    tic;
    sol = solver('x0', x0_guess, 'p', p_val, ...
                 'lbx', lbx, 'ubx', ubx, ...
                 'lbg', lbg, 'ubg', ubg);
    solve_times(step_idx) = toc;

    x_sol = full(sol.x);

    % Warm-start next step
    x0_guess = x_sol;

    % Extract first control input
    du_now = x_sol(n_cols + 1);  % U(1) is right after g
    u_now = du_now + u0_frozen;
    
    if isnan(u_now)
        warning('CasADi solver returned NaN at step %d! Fallback to u_prev.', t);
        u_now = u_prev;
    end
    
    u_now = max(u_min, min(u_max, u_now)); % Safety clamp
    u_hist(t) = u_now;
    
    if mod(t, 5) == 0 || t == 1
        fprintf('Step %d / %d | u_now: %.4f pps\n', t, N_total, u_now);
    end

    % --- Diagnostic (first control step only) ---
    if step_idx == 1
        Y_pred = x_sol(n_cols + N_pred + 1 : n_cols + 2*N_pred);
        fprintf('\n=== DeePC Diagnostic (step 1) ===\n');
        fprintf('u_now = %.4f (u0 = %.4f, diff = %.4f)\n', u_now, u0, u_now - u0);
        fprintf('y_ini range: [%.0f, %.0f]\n', min(y_ini), max(y_ini));
        fprintf('ref_seg range: [%.0f, %.0f]\n', min(ref_seg), max(ref_seg));
        fprintf('Y_pred range: [%.0f, %.0f]\n', min(Y_pred), max(Y_pred));
        fprintf('||g|| = %.6f\n', norm(x_sol(1:n_cols)));
        if abs(u_now - u0) < 1e-4
            warning('DeePC:flat', 'Controller NOT moving - check tuning!');
        end
        fprintf('================================\n\n');
    end

    % Simulate plant
    if FAST_SURROGATE_MODE
        [~, k_plant] = min(abs(wf_values - u_prev));
        A_p = A_all(:,:,k_plant); B_p = B_all(:,:,k_plant);
        C_p = C_all(:,:,k_plant); D_p = D_all(:,:,k_plant);
        u0_p = offset.u(:,k_plant); z0_p = offset.z(:,k_plant);
        
        y_meas = z0_p(idx_track) + C_p(idx_track,:)*x_plant_fast + D_p(idx_track,:)*(u_hist(t) - u0_p);
        x_plant_fast = A_p*x_plant_fast + B_p*(u_hist(t) - u0_p);
        y_hist(t) = y_meas;
    else
        t_step = (0:t)'*Ts;
        u_step = [u0; u_hist(1:t)];
        MWS.In.Tsim = t_step(end);
        MWS.In.Alt = timeseries([0;0],[0 t_step(end)],'Name','Alt');
        MWS.In.MN  = timeseries([0;0],[0 t_step(end)],'Name','MN');
        MWS.In.dT  = timeseries([0;0],[0 t_step(end)],'Name','dT');
        MWS.In.dVBV = timeseries([0;0],[0 t_step(end)],'Name','dVBV');
        MWS.In.dNz  = timeseries([0;0],[0 t_step(end)],'Name','dNz');
        MWS.In.Wf = timeseries(u_step,t_step,'Name','Wf');

        simIn = Simulink.SimulationInput(model);
        simIn = simIn.setVariable('MWS',MWS);
        simIn = simIn.setModelParameter('StartTime','0','StopTime',num2str(t_step(end)),...
            'ReturnWorkspaceOutputs','on');

        simOut = sim(simIn);
        out_step = fetch_out_dyn(simOut);
        y_hist(t) = measure_tracking_output(out_step, track_name);
    end

    u_prev = u_now;

    if mod(step_idx, 20) == 0
        fprintf('DeePC step %d / %d  (solve %.4f s, mean %.4f s)\n', ...
            step_idx, N_sim, solve_times(step_idx), mean(solve_times(1:step_idx)));
    end
end

%% ===== Plotting =====
time = (0:N_total-1)'*Ts;
ref_plot = ref(1:N_total);
controller_name = 'Linear DeePC';
save('Results/results_deepc.mat', 'y_hist', 'u_hist', 'solve_times', 'time', 'ref_plot', 'controller_name', 'track_name');
fprintf('Saved results to Results/results_deepc.mat\n');

if ~FAST_SURROGATE_MODE
    control_time = time(T_ini+1:end);

    figure('Name','DeePC CasADi closed-loop on nonlinear plant')
    subplot(3,1,1)
    plot(control_time,ref_plot(T_ini+1:N_total),'--','LineWidth',1.2)
    hold on
    plot(control_time,y_hist(T_ini+1:end),'LineWidth',1.4)
    grid on
    xlabel('Time [s]')
    ylabel(track_name,'Interpreter','tex')
    legend('Reference','AGTF30 measured output','Location','best')
    title('DeePC tracking (CasADi)')

    subplot(3,1,2)
    plot(control_time,u_hist(T_ini+1:end),'LineWidth',1.4)
    grid on
    xlabel('Time [s]')
    ylabel('W_f [pps]')
    title('Applied fuel command')

    subplot(3,1,3)
    plot(1:N_sim, solve_times*1000, 'LineWidth', 1.0)
    grid on
    xlabel('MPC step')
    ylabel('Solve time [ms]')
    title(sprintf('Solver time (mean = %.1f ms, max = %.1f ms)', ...
        mean(solve_times)*1000, max(solve_times)*1000))

    % RMSE
    y_ctrl = y_hist(T_ini+1:end);
    ref_ctrl = ref_plot(T_ini+1:N_total);
    rmse = sqrt(mean((y_ctrl - ref_ctrl).^2));
    fprintf('\n===== DeePC CasADi Results =====\n');
    fprintf('Tracking output: %s\n', track_name);
    fprintf('RMSE: %.1f\n', rmse);
    fprintf('Mean solve time: %.1f ms\n', mean(solve_times)*1000);
    fprintf('Max solve time:  %.1f ms\n', max(solve_times)*1000);
    fprintf('Hankel columns: %d\n', n_cols);
end


%% ===== Helper functions =====
function [U_hankel, Y_hankel] = build_deepc_hankel(training_folder, track_name, L)
    filelist = dir(fullfile(training_folder,'*.mat'));
    if isempty(filelist)
        filelist = dir(fullfile(training_folder,'**','*.mat'));
    end

    u_cells = {};
    y_cells = {};

    for iFile = 1:numel(filelist)
        data = load(fullfile(filelist(iFile).folder,filelist(iFile).name));
        out_dyn = fetch_out_dyn(data);
        u = extract_input_trace(out_dyn);
        y = measure_tracking_series(out_dyn, track_name);

        if numel(u) ~= numel(y)
            n_common = min(numel(u),numel(y));
            u = u(1:n_common);
            y = y(1:n_common);
        end

        if numel(u) >= L
            u_cells{end+1} = u(:); %#ok<AGROW>
            y_cells{end+1} = y(:); %#ok<AGROW>
        end
    end

    if isempty(u_cells)
        error('No usable DeePC training trajectories were found in %s.',training_folder);
    end

    U_hankel = [];
    Y_hankel = [];
    for iTraj = 1:numel(u_cells)
        U_hankel = [U_hankel, block_hankel(u_cells{iTraj}, L)]; %#ok<AGROW>
        Y_hankel = [Y_hankel, block_hankel(y_cells{iTraj}, L)]; %#ok<AGROW>
    end
end

function H = block_hankel(signal, L)
    signal = signal(:);
    n_cols = numel(signal) - L + 1;
    if n_cols < 1
        H = zeros(L,0);
        return;
    end

    H = zeros(L,n_cols);
    for i = 1:L
        H(i,:) = signal(i:i+n_cols-1).';
    end
end

function [U_hankel, Y_hankel] = reduce_hankel_columns(U_hankel, Y_hankel, max_cols)
    n_cols = size(U_hankel,2);
    if n_cols <= max_cols
        return;
    end

    keep_idx = unique(round(linspace(1,n_cols,max_cols)));
    U_hankel = U_hankel(:,keep_idx);
    Y_hankel = Y_hankel(:,keep_idx);
end

function out_dyn = fetch_out_dyn(data)
    if isstruct(data) && isfield(data,'out_Dyn')
        out_dyn = data.out_Dyn;
        return;
    end

    if isa(data,'Simulink.SimulationOutput')
        try
            out_dyn = data.get('out_Dyn');
            return;
        catch
        end
    end

    if isstruct(data) && isfield(data,'eng') && isfield(data,'cntrl')
        out_dyn = data;
        return;
    end

    error('Could not locate out_Dyn in the supplied data file.');
end

function u = extract_input_trace(out_dyn)
    if isfield(out_dyn,'cntrl')
        ctrl = out_dyn.cntrl;
        if isfield(ctrl,'Wfact')
            u = series_to_vector(ctrl.Wfact);
            return;
        end
        if isfield(ctrl,'Wf')
            u = series_to_vector(ctrl.Wf);
            return;
        end
    end

    error('Could not find a Wf input trace in the logged data.');
end

function y = measure_tracking_series(out_dyn, track_name)
    track_name_lower = lower(string(track_name));
    if contains(track_name_lower,'n_h')
        y = series_to_vector(out_dyn.eng.Shaft.N_HPC);
    elseif contains(track_name_lower,'n_l')
        y = series_to_vector(out_dyn.eng.Shaft.N_Fan);
    elseif contains(track_name_lower,'thrust') || contains(track_name_lower,'fnet')
        y = series_to_vector(out_dyn.eng.Perf.Fnet);
    else
        y = series_to_vector(out_dyn.eng.Shaft.N_HPC);
    end
end

function y = measure_tracking_output(out_dyn, track_name)
    track_name_lower = lower(string(track_name));
    if contains(track_name_lower,'n_h')
        y = series_to_value(out_dyn.eng.Shaft.N_HPC);
    elseif contains(track_name_lower,'n_l')
        y = series_to_value(out_dyn.eng.Shaft.N_Fan);
    elseif contains(track_name_lower,'thrust') || contains(track_name_lower,'fnet')
        y = series_to_value(out_dyn.eng.Perf.Fnet);
    else
        y = series_to_value(out_dyn.eng.Shaft.N_HPC);
    end
end

function x = series_to_vector(series_like)
    if isstruct(series_like) && isfield(series_like,'Data')
        x = series_like.Data(:);
    elseif isa(series_like,'timeseries')
        x = series_like.Data(:);
    else
        x = series_like(:);
    end
    x = double(x);
end

function x = series_to_value(series_like)
    x = series_to_vector(series_like);
    x = x(end);
end
