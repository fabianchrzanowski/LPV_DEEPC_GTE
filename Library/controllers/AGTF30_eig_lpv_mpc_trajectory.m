% ========================================================================
% LPV-MPC trajectory tracking with eigenvalue-interpolated models
%
% This script runs a gain-scheduled LPV-MPC controller on the nonlinear
% AGTF30 Simulink plant. It uses the eigenvalue-interpolated LUT from
% AGTF30_sys_id_4.m for the A matrix, guaranteeing stability at every
% intermediate operating point.
%
% The user can define any reference trajectory (thrust/N_H or fuel flow)
% at the top of the script.
%
% ========================================================================
% For AGTF30 Dissertation
% Data-driven MPC with stability-guaranteed LPV interpolation
% ========================================================================

clearvars; close all; clc
yalmip('clear');

solver_name = 'mosek';

proj = currentProject;
cd(proj.RootFolder)

%% ===== LOAD EIGENVALUE-INTERPOLATED LUT =====
lpv_file = fullfile(proj.RootFolder,'Params','tables','models',...
    'AGTF30_model_tables__1_input.mat');
if ~isfile(lpv_file)
    error(['LPV table file not found: %s.\n'...
        'Run Library/system_identification/AGTF30_sys_id_4.m first.'], lpv_file);
end

LUT = load(lpv_file);
Ts = LUT.matrices.Ts;

%% ===== MPC PARAMETERS =====
N_pred   = 20;       % prediction horizon
N_sim    = 160;      % total simulation steps

u_min    = 0.2;      % fuel flow lower bound [pps]
u_max    = 2.2;      % fuel flow upper bound [pps]
du_max   = 0.04;     % max fuel flow rate change per step

q_track  = 1;        % tracking weight
r_du     = 25;       % rate-of-change penalty

%% ===== OUTPUT SELECTION =====
output_names_lower = lower(string(LUT.names.z));

% Try to track N_H (high-pressure spool speed)
idx_track = find(output_names_lower == "n_h", 1);
if isempty(idx_track)
    idx_track = find(output_names_lower == "n_l", 1);
end
if isempty(idx_track)
    idx_track = 1;
end
track_name = LUT.names.z{idx_track};
fprintf('Tracking output: %s (index %d)\n', track_name, idx_track);

%% ===== USER-DEFINED REFERENCE TRAJECTORY =====
% Choose one of the following trajectory types, or define your own.
% The reference is defined in terms of the tracked output (e.g. N_H in rpm).

% Nominal operating point
rho_eq = 1.2;  % nominal fuel flow for scheduling
z0_nom = off2mat(LUT.offset.z, rho_eq);
track0 = z0_nom(idx_track);

t_sim = (0:N_sim + N_pred - 1)' * Ts;

% --- Trajectory type selector ---
trajectory_type = 'sinusoidal';
% Options: 'sinusoidal', 'step_sequence', 'takeoff_profile', 'custom'

switch trajectory_type
    case 'sinusoidal'
        % Sinusoidal demand around nominal
        amp = 300;    % amplitude [rpm]
        period = 3;   % [seconds]
        ref = track0 + amp * sin(2*pi*t_sim / period);

    case 'step_sequence'
        % Step sequence: idle -> cruise -> max -> cruise -> idle
        z0_low  = off2mat(LUT.offset.z, 0.5);
        z0_mid  = off2mat(LUT.offset.z, 1.2);
        z0_high = off2mat(LUT.offset.z, 1.8);

        step_times  = [0, 1.5, 3.0, 4.5, 6.0];
        step_values = [z0_low(idx_track), z0_mid(idx_track), ...
                       z0_high(idx_track), z0_mid(idx_track), ...
                       z0_low(idx_track)];
        ref = interp1(step_times, step_values, t_sim, 'previous', 'extrap');

    case 'takeoff_profile'
        % Realistic take-off: idle -> spool up -> climb -> cruise
        z0_idle   = off2mat(LUT.offset.z, 0.5);
        z0_climb  = off2mat(LUT.offset.z, 1.6);
        z0_cruise = off2mat(LUT.offset.z, 1.2);

        profile_times  = [0,  0.5,  2.0,  3.5,  5.5];
        profile_values = [z0_idle(idx_track), z0_idle(idx_track), ...
                          z0_climb(idx_track), z0_climb(idx_track), ...
                          z0_cruise(idx_track)];
        ref = interp1(profile_times, profile_values, t_sim, 'pchip', 'extrap');

    case 'custom'
        % Define your own: ref must be a vector of length (N_sim + N_pred)
        ref = track0 * ones(size(t_sim));
        warning('Using constant reference. Edit the custom case to define your trajectory.');

    otherwise
        error('Unknown trajectory_type: %s', trajectory_type);
end

%% ===== SIMULINK PLANT SETUP =====
MWS.engName = 'AGTF30';
MWS = setup_Controller(MWS);
MWS = setup_AllEng(MWS);

busVars = load(fullfile(proj.RootFolder,'Params','Eng_Bus.mat'));
busVarNames = fieldnames(busVars);
for ib = 1:numel(busVarNames)
    assignin('base', busVarNames{ib}, busVars.(busVarNames{ib}));
end

% Initial scheduling point
u0_init = off2mat(LUT.offset.u, rho_eq);

MWS.In.Ts   = Ts;
MWS.In.Tsim = Ts;
MWS.In.Alt  = timeseries([0;0],[0 Ts],'Name','Alt');
MWS.In.MN   = timeseries([0;0],[0 Ts],'Name','MN');
MWS.In.dT   = timeseries([0;0],[0 Ts],'Name','dT');
MWS.In.dVBV = timeseries([0;0],[0 Ts],'Name','dVBV');
MWS.In.dNz  = timeseries([0;0],[0 Ts],'Name','dNz');
MWS.In.Wf   = timeseries([u0_init; u0_init],[0 Ts],'Name','Wf');
MWS = AGTF30_initial_conditions(MWS);

model = 'AGTF30SysDyn';

%% ===== INITIAL STATE =====
x_dev   = zeros(size(LPV2mat(LUT.matrices.A, rho_eq), 1), 1);
u_prev  = u0_init;
rho_now = rho_eq;

%% ===== STORAGE =====
u_hist   = zeros(N_sim, 1);
y_hist   = nan(N_sim, 1);
rho_hist = nan(N_sim, 1);
z_full   = nan(N_sim, numel(LUT.names.z));

ops = sdpsettings('solver', solver_name, 'verbose', 1);

%% ===== MPC LOOP =====
fprintf('Starting LPV-MPC trajectory tracking (%d steps)...\n', N_sim);
t_start = tic;

for t = 1:N_sim
    t
    % --- Gain scheduling: update linear model from current rho ---
    A = LPV2mat(LUT.matrices.A, rho_now);
    B = LPV2mat(LUT.matrices.B, rho_now);
    C = LPV2mat(LUT.matrices.C, rho_now);
    D = LPV2mat(LUT.matrices.D, rho_now);

    u0_step = off2mat(LUT.offset.u, rho_now);
    z0_step = off2mat(LUT.offset.z, rho_now);

    % --- Formulate QP ---
    U = sdpvar(1, N_pred, 'full');
    x_pred = x_dev;
    constraints = [];
    objective   = 0;

    for i = 1:N_pred
        x_pred = A * x_pred + B * (U(i) - u0_step);
        z_pred = z0_step + C * x_pred + D * (U(i) - u0_step);

        constraints = [constraints, u_min <= U(i), U(i) <= u_max];

        if i == 1
            du_i = U(i) - u_prev;
        else
            du_i = U(i) - U(i-1);
        end

        constraints = [constraints, -du_max <= du_i, du_i <= du_max];
        objective = objective ...
            + q_track * (z_pred(idx_track) - ref(t + i - 1))^2 ...
            + r_du * du_i^2;
    end

    % --- Solve ---
    diagnostics = optimize(constraints, objective, ops);
    if diagnostics.problem ~= 0
        warning('LPV-MPC optimization failed at step %d (code %d). Holding previous input.', ...
            t, diagnostics.problem);
        u_now = u_prev;
    else
        u_now = value(U(1));
    end

    u_hist(t)   = u_now;
    rho_hist(t) = rho_now;

    % --- Simulate nonlinear plant up to current step ---
    t_step = (0:t)' * Ts;
    u_step = [u0_init; u_hist(1:t)];

    MWS.In.Tsim = t_step(end);
    MWS.In.Alt  = timeseries([0;0],[0 t_step(end)],'Name','Alt');
    MWS.In.MN   = timeseries([0;0],[0 t_step(end)],'Name','MN');
    MWS.In.dT   = timeseries([0;0],[0 t_step(end)],'Name','dT');
    MWS.In.dVBV = timeseries([0;0],[0 t_step(end)],'Name','dVBV');
    MWS.In.dNz  = timeseries([0;0],[0 t_step(end)],'Name','dNz');
    MWS.In.Wf   = timeseries(u_step, t_step, 'Name', 'Wf');

    simIn = Simulink.SimulationInput(model);
    simIn = simIn.setVariable('MWS', MWS);
    simIn = simIn.setModelParameter('StartTime','0', ...
        'StopTime', num2str(t_step(end)), ...
        'ReturnWorkspaceOutputs','on');

    simOut = sim(simIn);

    % Extract plant output
    out_step = [];
    try
        out_step = simOut.get('out_Dyn');
    catch
        out_step = [];
    end
    if isempty(out_step)
        try
            out_step = evalin('base','out_Dyn');
        catch
            error('Could not retrieve out_Dyn at step %d', t);
        end
    end

    if contains(lower(track_name),'n_h')
        y_meas = out_step.eng.Shaft.N_HPC.Data(end);
    elseif contains(lower(track_name),'n_l')
        y_meas = out_step.eng.Shaft.N_Fan.Data(end);
    else
        y_meas = z0_step(idx_track) + C(idx_track,:)*x_dev ...
            + D(idx_track,:)*(u_now - u0_step);
    end

    y_hist(t) = y_meas;

    % --- Update LPV state and scheduling ---
    x_dev   = A * x_dev + B * (u_now - u0_step);
    u_prev  = u_now;
    rho_now = u_now;  % scheduling variable = fuel flow

    % --- Progress ---
    if mod(t, 50) == 0
        elapsed = toc(t_start);
        fprintf('Step %d / %d  (%.1f s elapsed, ~%.0f s remaining)\n', ...
            t, N_sim, elapsed, elapsed/t*(N_sim-t));
    end
end

elapsed = toc(t_start);
fprintf('LPV-MPC trajectory tracking complete. Total time: %.1f s\n', elapsed);

%% ===== RESULTS PLOTS =====
time     = (0:N_sim-1)' * Ts;
ref_plot = ref(1:N_sim);

figure('Name','Eigenvalue-interpolated LPV-MPC: Trajectory Tracking')

% --- Tracking performance ---
subplot(3,1,1)
plot(time, ref_plot, '--', 'LineWidth', 1.2, 'Color', [0.5 0.5 0.5])
hold on
plot(time, y_hist, 'LineWidth', 1.4, 'Color', [0 0.45 0.74])
grid on
xlabel('Time [s]')
ylabel(track_name, 'Interpreter', 'tex')
legend('Reference', 'AGTF30 plant output', 'Location', 'best')
title(sprintf('Trajectory tracking: %s', trajectory_type), 'Interpreter', 'none')

% --- Fuel flow command ---
subplot(3,1,2)
plot(time, u_hist, 'LineWidth', 1.4, 'Color', [0.85 0.33 0.1])
grid on
xlabel('Time [s]')
ylabel('W_f [pps]', 'Interpreter', 'none')
title('MPC fuel flow command')
yline(u_min, '--r', 'u_{min}');
yline(u_max, '--r', 'u_{max}');

% --- Scheduling trajectory ---
subplot(3,1,3)
plot(time, rho_hist, 'LineWidth', 1.4, 'Color', [0.47 0.67 0.19])
grid on
xlabel('Time [s]')
ylabel('\rho = W_f [pps]', 'Interpreter', 'none')
title('Scheduling variable trajectory')

% --- Tracking error ---
figure('Name','Tracking Error')
e_track = y_hist - ref_plot;
plot(time, e_track, 'LineWidth', 1.2)
grid on
xlabel('Time [s]')
ylabel('Error [rpm]')
title(sprintf('Tracking error (RMSE = %.1f rpm)', rms(e_track(~isnan(e_track)))))

fprintf('Tracking RMSE: %.2f\n', rms(e_track(~isnan(e_track))));
fprintf('Max absolute error: %.2f\n', max(abs(e_track(~isnan(e_track)))));
fprintf('Trajectory type: %s\n', trajectory_type);
