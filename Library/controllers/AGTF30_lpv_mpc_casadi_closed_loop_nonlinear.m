run('setup_simulation_params.m');   

% LPV-MPC 
% uses the identified models from the PRBS data
% same idea of scheduling as LPV-DeePC
% uses Casadi and IPOPT

import casadi.*

try; proj = currentProject; catch; proj = openProject('GTE.prj'); end
cd(proj.RootFolder)

lpv_file = fullfile(proj.RootFolder,'Params','tables','models','AGTF30_model_tables__1_input.mat');
if ~isfile(lpv_file)
    error('LPV table file not found: %s. Run Library/system_identification/AGTF30_sys_id_4.m first.',lpv_file);
end

LUT = load(lpv_file);

% Variables now sourced from setup_simulation_params.m

u_min = 0.35;
u_max = 2.0;
du_max = 0.15;

q_track = 15;
r_du = 1;
r_u = 1;
% Scheduling variable is fuel flow Wf
rho_eq = wf_idle;

A_eq = LPV2mat(LUT.matrices.A,rho_eq);
B_eq = LPV2mat(LUT.matrices.B,rho_eq);
C_eq = LPV2mat(LUT.matrices.C,rho_eq);
D_eq = LPV2mat(LUT.matrices.D,rho_eq);

u0 = off2mat(LUT.offset.u,rho_eq);
z0 = off2mat(LUT.offset.z,rho_eq);

track_name = LUT.names.z{idx_NH};

% The trajectory ref is now sourced from setup_simulation_params.m

X_abs_meas = [off2mat(LUT.offset.x, rho_eq); off2mat(LUT.offset.y, rho_eq)];
x_dev = zeros(size(A_eq,1),1);

u_prev = u0;
rho_now = rho_eq;
rho_prev = rho_eq;

u_hist = zeros(N_total,1);
y_hist = nan(N_total,1);
nl_hist = nan(N_total,1);
nh_hist = nan(N_total,1);
rho_hist = nan(N_total,1);
solve_times = nan(N_total,1);
cost_hist = nan(N_total,1);
cost_track_hist = nan(N_total,1);
cost_du_hist = nan(N_total,1);

% Pad initial history
u_hist(1:T_ini) = u0;
y_hist(1:T_ini) = ref(1:T_ini);
rho_hist(1:T_ini) = rho_eq;


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

model = 'AGTF30SysDyn';

for t = T_ini+1:N_total
    % ---- STEP 1: Gain scheduling ----
    % Look up the A, B, C, D matrices and offsets for the current fuel flow.
    A = LPV2mat(LUT.matrices.A,rho_now);
    B = LPV2mat(LUT.matrices.B,rho_now);
    C = LPV2mat(LUT.matrices.C,rho_now);
    D = LPV2mat(LUT.matrices.D,rho_now);

    u0_step = off2mat(LUT.offset.u,rho_now);
    z0_step = off2mat(LUT.offset.z,rho_now);

    x0_new = off2mat(LUT.offset.x, rho_now);
    y0_new = off2mat(LUT.offset.y, rho_now);
    X0_new = [x0_new; y0_new];

    % ---- STEP 2: Full state feedback ----
    % Keep the controller state exactly aligned with the measured physics.
    x_dev = X_abs_meas - X0_new;

    % ---- STEP 3: Build the MPC optimisation ----
    opti = casadi.Opti();
    U = opti.variable(1, N_pred);
    
    % Bias correction: difference between last measurement and model prediction.
    % This eliminates steady-state offset caused by model mismatch.
    if t == T_ini+1
        y_bias = 0;
    else
        y_model = z0_step(idx_NH) + C(idx_NH,:)*x_dev + D(idx_NH,:)*(u_prev - u0_step);
        y_bias = y_meas - y_model;
    end

    x_pred = x_dev;
    objective = 0;
    C_y = 0; C_du = 0; C_u = 0;

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

        % Accumulate cost: tracking + input penalty + rate penalty
        % Bias correction added to the predicted output
        C_y = C_y + q_track*(z_pred(idx_NH) + y_bias - ref(t+i-1))^2;
        C_u = C_u + r_u*(U(i) - u0_step)^2;
        C_du = C_du + r_du*(du_i)^2;
        objective = objective + q_track*(z_pred(idx_NH) + y_bias - ref(t+i-1))^2 + r_u*(U(i) - u0_step)^2 + r_du*(du_i)^2;
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
        
        cost_hist(t) = sol.value(objective);
        cost_track_hist(t) = sol.value(C_y);
        cost_du_hist(t) = sol.value(C_du + C_u);
    catch
        error('CasADi LPV-MPC optimization failed at step %d',t);
    end

    u_hist(t) = u_now;
    rho_hist(t) = rho_now;
    
    if mod(t, 5) == 0 || t == T_ini+1
        fprintf('LPV closed-loop step %d / %d | u_now: %.4f pps | rho: %.4f | avg solve: %.1f ms\n', t, N_total, u_now, rho_now, mean(solve_times(T_ini+1:t))*1000);
    end

    % ---- STEP 6: Simulate the nonlinear plant for one step ----
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

    % Extract the plant output from the Simulink results
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
            error('Could not retrieve out_Dyn at step %d',t);
        end
    end

    if contains(lower(track_name),'n_h')
        y_meas = out_step.eng.Shaft.N_HPC.Data(end);
    elseif contains(lower(track_name),'n_l')
        y_meas = out_step.eng.Shaft.N_Fan.Data(end);
    else
        y_meas = z0_step(idx_NH) + C(idx_NH,:)*x_dev + D(idx_NH,:)*(u_now-u0_step);
    end
    
    % ---- STEP 7: Full state feedback update for next step ----
    % Read all four states from the nonlinear plant to keep the model aligned.
    NL_meas = out_step.eng.Shaft.N_Fan.Data(end);
    NH_meas = out_step.eng.Shaft.N_HPC.Data(end);
    Pt21_meas = out_step.eng.S21.Pt.Data(end);
    Pt25_meas = out_step.eng.S25.Pt.Data(end);
    X_abs_meas = [NL_meas; NH_meas; Pt21_meas; Pt25_meas];


    y_hist(t) = y_meas;
    nl_hist(t) = NL_meas;
    nh_hist(t) = NH_meas;

    % Update previous input and scheduling variable for the next iteration
    u_prev = u_now;
    rho_prev = rho_now;
    
    % Clamp scheduling variable to bounds and apply low-pass filter
    % to prevent chattering (same as LPV-DeePC).
    u_sched_target = max(min(u_now, max(wf_values)), min(wf_values));
    alpha = 0.3; % Matched to DeePC 
    rho_now = alpha * u_sched_target + (1 - alpha) * rho_now;

end

time = (0:N_total-1)'*Ts;
ref_plot = ref(1:N_total);
scen_tag = lower(SCENARIO);
mat_filename = sprintf('Results/results_lpv_mpc_casadi_%s.mat', scen_tag);
controller_name = sprintf('LPV-MPC (CasADi) %s', SCENARIO);
save(mat_filename, 'y_hist', 'u_hist', 'rho_hist', 'nl_hist', 'nh_hist', 'cost_hist', 'cost_track_hist', 'cost_du_hist', 'solve_times', 'time', 'ref_plot', 'controller_name', 'track_name');
fprintf('Saved results to %s\n', mat_filename);
fprintf('==============================================\n');
fprintf('Mean solver time: %.2f ms per step\n', mean(solve_times(T_ini+1:end))*1000);
fprintf('==============================================\n');


figure('Name','Closed-loop LPV-MPC (CasADi) on nonlinear plant')
subplot(3,1,1)
plot(time,ref_plot,'--','LineWidth',1.2)
hold on
plot(time,y_hist,'LineWidth',1.4)
grid on
xlabel('Time [s]')
ylabel(track_name,'Interpreter','tex')
legend('Reference','AGTF30 measured output','Location','best')

subplot(3,1,2)
plot(time,u_hist,'LineWidth',1.4)
grid on
xlabel('Time [s]')
ylabel('W_f [pps]','Interpreter','none')
title('Applied fuel command')

subplot(3,1,3)
plot(time,rho_hist,'LineWidth',1.4)
grid on
xlabel('Time [s]')
ylabel('W_f eq [pps]')
title('Scheduling variable \rho')

disp('Closed-loop LPV-MPC (CasADi) run complete.')
disp(['Tracking output: ', track_name])

disp(['Nominal scheduling point rho_eq = ', num2str(rho_eq)])
