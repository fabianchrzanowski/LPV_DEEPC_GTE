clearvars;
% Always load latest params
run('setup_simulation_params.m');
% ========================================================================
% LPV-MPC: Disturbance Rejection Test
%
% Applies same disturbances at the DeePC to test the disturbance rejection
%
% Uses CasADi + IPOPT
% ========================================================================

try proj = currentProject; catch; proj = openProject('GTE.prj'); end
cd(proj.RootFolder)

import casadi.*

lpv_file = fullfile(proj.RootFolder,'Params','tables','models','AGTF30_model_tables__1_input.mat');
if ~isfile(lpv_file)
    error('LPV table file not found: %s. Run Library/system_identification/AGTF30_sys_id_4.m first.',lpv_file);
end

LUT = load(lpv_file);
wf_data = load(fullfile(proj.RootFolder,'Params','test_files_fuel_only','test_identification','wf_values.mat'));
wf_values = wf_data.Wf_values(:)';

%% ===== Controller tuning =====
u_min = 0.35;
u_max = 2.0;
du_max = 0.15;

q_track = 15;
r_du = 1;

% Scheduling variable is fuel flow Wf
rho_eq = wf_idle;

A_eq = LPV2mat(LUT.matrices.A,rho_eq);
B_eq = LPV2mat(LUT.matrices.B,rho_eq);
C_eq = LPV2mat(LUT.matrices.C,rho_eq);
D_eq = LPV2mat(LUT.matrices.D,rho_eq);

u0 = off2mat(LUT.offset.u,rho_eq);
z0 = off2mat(LUT.offset.z,rho_eq);

track_name = LUT.names.z{idx_NH};

%% ===== Disturbance Configuration =====
% EXACTLY matches the DeePC disturbance test for a fair comparison
DIST_OUTPUT = [
    round(N_sim*0.20), round(N_sim*0.35), +500;    
    round(N_sim*0.55), round(N_sim*0.70), -400;    
];

DIST_INPUT = [
    % round(N_sim*0.40), round(N_sim*0.48), +0.15;   
    0;
];

fprintf('\n--- LPV-MPC Disturbance Rejection Test ---\n');
fprintf('Output disturbances:\n');
for i = 1:size(DIST_OUTPUT,1)
    fprintf('  Steps %d-%d: %+.0f RPM\n', DIST_OUTPUT(i,1), DIST_OUTPUT(i,2), DIST_OUTPUT(i,3));
end
fprintf('Input disturbances:\n');
for i = 1:size(DIST_INPUT,1)
    fprintf('  Steps %d-%d: %+.3f pps\n', DIST_INPUT(i,1), DIST_INPUT(i,2), DIST_INPUT(i,3));
end

%% ===== Initialization =====
x_dev = zeros(size(A_eq,1),1);

u_prev = u0;
rho_now = rho_eq;

u_hist = zeros(N_sim,1);
u_applied = zeros(N_sim,1);
y_hist = nan(N_sim,1);
y_true = nan(N_sim,1);
dist_y_hist = zeros(N_sim,1);
dist_u_hist = zeros(N_sim,1);
rho_hist = nan(N_sim,1);
solve_times = nan(N_sim,1);

MWS = struct();
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

model = 'AGTF30SysDyn';

%% ===== Closed-loop with disturbances =====
for t = 1:N_sim
    % ---- STEP 1: Compute disturbances active at this step ----
    % d_y = additive output disturbance (sensor bias / load change)
    % d_u = additive input disturbance (actuator fault)
    d_y = 0;
    for di = 1:size(DIST_OUTPUT,1)
        if t >= DIST_OUTPUT(di,1) && t <= DIST_OUTPUT(di,2)
            d_y = d_y + DIST_OUTPUT(di,3);
        end
    end
    d_u = 0;
    for di = 1:size(DIST_INPUT,1)
        if t >= DIST_INPUT(di,1) && t <= DIST_INPUT(di,2)
            d_u = d_u + DIST_INPUT(di,3);
        end
    end
    dist_y_hist(t) = d_y;
    dist_u_hist(t) = d_u;

    % ---- STEP 2: Gain scheduling ----
    % Look up the A, B, C, D matrices and offsets for the current fuel flow.
    A = LPV2mat(LUT.matrices.A,rho_now);
    B = LPV2mat(LUT.matrices.B,rho_now);
    C = LPV2mat(LUT.matrices.C,rho_now);
    D = LPV2mat(LUT.matrices.D,rho_now);

    u0_step = off2mat(LUT.offset.u,rho_now);
    z0_step = off2mat(LUT.offset.z,rho_now);

    % ---- STEP 3: Build the MPC optimisation ----
    opti = casadi.Opti();
    U = opti.variable(1, N_pred);
    
    % Bias correction: difference between last measurement and model prediction.
    % This helps reject disturbances the model doesn't know about.
    if t == 1
        y_bias = 0;
    else
        y_model = z0_step(idx_NH) + C(idx_NH,:)*x_dev + D(idx_NH,:)*(u_prev - u0_step);
        y_bias = y_meas - y_model;
    end

    x_pred = x_dev;
    objective = 0;

    % ---- STEP 4: Roll out the prediction over N_pred steps ----
    for i = 1:N_pred
        % Propagate the state one step using the LPV model
        x_pred = A*x_pred + B*(U(i)-u0_step);
        z_pred = z0_step + C*x_pred + D*(U(i)-u0_step);

        % Input constraints: actuator limits and rate limits
        opti.subject_to( u_min <= U(i) <= u_max );

        if i == 1
            du_i = U(i) - u_prev;
        else
            du_i = U(i) - U(i-1);
        end

        opti.subject_to( -du_max <= du_i <= du_max );

        % Accumulate cost: tracking + rate penalty (with bias correction)
        objective = objective + q_track*(z_pred(idx_NH) + y_bias - ref(t+i-1))^2 + r_du*(du_i)^2;
    end

    % ---- STEP 5: Solve the MPC QP ----
    opti.minimize(objective);
    p_opts = struct('print_time', false);
    s_opts = struct('print_level', 0, 'sb', 'yes');
    opti.solver('ipopt', p_opts, s_opts);

    try
        tic;
        sol = opti.solve();
        solve_times(t) = toc;
        u_now = sol.value(U(1));  % apply only the first input (receding horizon)
    catch
        warning('CasADi LPV-MPC optimization failed at step %d',t);
        u_now = u_prev;
    end

    u_hist(t) = u_now;
    
    % ---- STEP 6: Apply input disturbance AFTER the controller ----
    % The controller doesn't know about d_u — it corrupts the applied input.
    u_actual = max(u_min, min(u_max, u_now + d_u));
    u_applied(t) = u_actual;
    
    % ---- STEP 7: Simulate the nonlinear plant with the disturbed input ----

    t_step = (0:t)'*Ts;
    u_step = [u0; u_applied(1:t)];

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

    % Extract the plant output from the Simulink results
    out_step = [];
    try
        out_step = simOut.get('out_Dyn');
    catch
    end
    if isempty(out_step)
        try
            out_step = evalin('base','out_Dyn');
        catch
            error('Could not retrieve out_Dyn at step %d',t);
        end
    end

    if contains(lower(track_name),'n_h')
        y_plant = out_step.eng.Shaft.N_HPC.Data(end);
    elseif contains(lower(track_name),'n_l')
        y_plant = out_step.eng.Shaft.N_Fan.Data(end);
    else
        y_plant = z0_step(idx_NH) + C(idx_NH,:)*x_dev + D(idx_NH,:)*(u_now-u0_step);
    end
    
    y_true(t) = y_plant;             % true undisturbed plant output
    y_meas = y_plant + d_y;          % controller sees this (with output disturbance)


    y_hist(t) = y_meas;
    
    % ---- STEP 8: Propagate the internal model state ----
    % NOTE: State update uses the internal model and the un-disturbed input.
    % It does NOT use y_meas to correct the state (no observer).
    x_dev = A*x_dev + B*(u_now-u0_step);
    
    % Update previous input and scheduling variable for the next iteration
    u_prev = u_now;
    rho_now = max(min(u_now, max(wf_values)), min(wf_values));
    rho_hist(t) = rho_now;
    
    if mod(t, 20) == 0 || t == 1
        dist_str = '';
        if d_y ~= 0; dist_str = [dist_str sprintf(' d_y=%+.0f', d_y)]; end
        if d_u ~= 0; dist_str = [dist_str sprintf(' d_u=%+.3f', d_u)]; end
        fprintf('  Step %d/%d | u_now: %.4f | solve: %.1fms%s\n', ...
            t, N_sim, u_now, solve_times(t)*1000, dist_str);
    end

end

%% ===== Plot Results =====
time = (0:N_sim-1)'*Ts;
ref_plot = ref(1:N_sim);

figure('Name','LPV-MPC Disturbance Rejection')

% --- Subplot 1: Tracking ---
subplot(4,1,1)
plot(time, ref_plot, '--', 'LineWidth',1.2, 'Color',[0.5 0.5 0.5])
hold on
plot(time, y_hist, 'Color', [0.85 0.33 0.10], 'LineWidth', 1.8)
plot(time, y_true, '--', 'Color', [0.85 0.45 0.45], 'LineWidth', 1.2)
% Shade disturbance regions
for di = 1:size(DIST_OUTPUT,1)
    t1 = (DIST_OUTPUT(di,1)-1)*Ts;
    t2 = (DIST_OUTPUT(di,2)-1)*Ts;
    yl = ylim;
    patch([t1 t2 t2 t1], [yl(1) yl(1) yl(2) yl(2)], ...
        [1 0.8 0.8], 'FaceAlpha',0.25, 'EdgeColor','none');
end
for di = 1:size(DIST_INPUT,1)
    t1 = (DIST_INPUT(di,1)-1)*Ts;
    t2 = (DIST_INPUT(di,2)-1)*Ts;
    yl = ylim;
    patch([t1 t2 t2 t1], [yl(1) yl(1) yl(2) yl(2)], ...
        [0.8 0.8 1], 'FaceAlpha',0.25, 'EdgeColor','none');
end
grid on
xlabel('Time [s]')
ylabel(track_name,'Interpreter','tex')
legend('Reference', 'Measured output (disturbed)', 'Controller belief (undisturbed)', 'Location', 'best')
title('Tracking with Disturbances (LPV-MPC)')

% --- Subplot 2: Control input ---
subplot(4,1,2)
stairs(time, u_hist, 'b', 'LineWidth',1.4)
hold on
stairs(time, u_applied, 'r--', 'LineWidth',1.0)
grid on
xlabel('Time [s]')
ylabel('W_f [pps]')
legend('Controller output','Actual applied (with d_u)','Location','best')
title('Fuel Command')

% --- Subplot 3: Disturbance signals ---
subplot(4,1,3)
yyaxis left
stairs(time, dist_y_hist, 'r', 'LineWidth',1.4)
ylabel('Output dist. [RPM]')
yyaxis right
stairs(time, dist_u_hist, 'b', 'LineWidth',1.4)
ylabel('Input dist. [pps]')
grid on
xlabel('Time [s]')
title('Applied Disturbances')

% --- Subplot 4: Tracking error ---
subplot(4,1,4)
err_true = y_true - ref_plot;
plot(time, err_true, 'b', 'LineWidth',1.2)
hold on
yline(0, 'k--')
grid on
xlabel('Time [s]')
ylabel('Error [RPM]')
title(sprintf('True Tracking Error (RMSE = %.1f RPM)', sqrt(mean(err_true.^2))))

sgtitle('LPV-MPC Disturbance Rejection Test', 'FontSize', 14, 'FontWeight', 'bold')

%% ===== Summary =====
rmse = sqrt(mean((y_true - ref_plot).^2));
fprintf('\n===== Disturbance Rejection Results =====\n');
fprintf('Tracking output: %s\n', track_name);
fprintf('RMSE (true output): %.1f RPM\n', rmse);
fprintf('Mean solve time: %.1f ms\n', mean(solve_times)*1000);

scen_tag = lower(SCENARIO);
mat_filename = sprintf('Results/results_lpv_mpc_disturbance_%s.mat', scen_tag);
controller_name = sprintf('LPV-MPC (Disturbance) %s', SCENARIO);
save(mat_filename, ...
    'y_hist', 'y_true', 'u_hist', 'u_applied', ...
    'dist_y_hist', 'dist_u_hist', 'rho_hist', ...
    'solve_times', 'time', 'ref_plot', 'controller_name', 'track_name', ...
    'DIST_OUTPUT', 'DIST_INPUT');
fprintf('Saved results to %s\n', mat_filename);
