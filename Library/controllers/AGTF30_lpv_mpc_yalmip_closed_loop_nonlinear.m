if ~exist('FAST_SURROGATE_MODE', 'var')
    clearvars; close all; clc
    yalmip('clear');
    run('setup_simulation_params.m');
end

solver_name = 'mosek';

proj = currentProject;
cd(proj.RootFolder)

lpv_file = fullfile(proj.RootFolder,'Params','tables','models','AGTF30_model_tables__1_input.mat');
if ~isfile(lpv_file)
    error('LPV table file not found: %s. Run Library/system_identification/AGTF30_sys_id_4.m first.',lpv_file);
end

LUT = load(lpv_file);

% Variables now sourced from setup_simulation_params.m

u_min = 0.2;
u_max = 2.2;
du_max = 0.04;

q_track = 1;
r_du = 25;

% Scheduling variable is fuel flow Wf
rho_eq = 1.2;

A_eq = LPV2mat(LUT.matrices.A,rho_eq);
B_eq = LPV2mat(LUT.matrices.B,rho_eq);
C_eq = LPV2mat(LUT.matrices.C,rho_eq);
D_eq = LPV2mat(LUT.matrices.D,rho_eq);

u0 = off2mat(LUT.offset.u,rho_eq);
z0 = off2mat(LUT.offset.z,rho_eq);

track_name = LUT.names.z{idx_NH};

% The trajectory ref is now sourced from setup_simulation_params.m

x_dev = zeros(size(A_eq,1),1);
if FAST_SURROGATE_MODE
    x_plant_fast = zeros(size(A_eq,1),1);
end
u_prev = u0;
rho_now = rho_eq;

u_hist = zeros(N_sim,1);
y_hist = nan(N_sim,1);
rho_hist = nan(N_sim,1);
solve_times = nan(N_sim,1);

ops = sdpsettings('solver',solver_name,'verbose',0);

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

for t = 1:N_sim
    % Gain scheduling: update linear model from current scheduling value
    A = LPV2mat(LUT.matrices.A,rho_now);
    B = LPV2mat(LUT.matrices.B,rho_now);
    C = LPV2mat(LUT.matrices.C,rho_now);
    D = LPV2mat(LUT.matrices.D,rho_now);

    u0_step = off2mat(LUT.offset.u,rho_now);
    z0_step = off2mat(LUT.offset.z,rho_now);

    U = sdpvar(1,N_pred,'full');
    x_pred = x_dev;
    constraints = [];
    objective = 0;

    for i = 1:N_pred
        x_pred = A*x_pred + B*(U(i)-u0_step);
        z_pred = z0_step + C*x_pred + D*(U(i)-u0_step);

        constraints = [constraints, u_min <= U(i), U(i) <= u_max];

        if i == 1
            du_i = U(i) - u_prev;
        else
            du_i = U(i) - U(i-1);
        end

        constraints = [constraints, -du_max <= du_i, du_i <= du_max];
        objective = objective + q_track*(z_pred(idx_NH) - ref(t+i-1))^2 + r_du*(du_i)^2;
    end

    tic;
    diagnostics = optimize(constraints,objective,ops);
    solve_times(t) = toc;
    if diagnostics.problem ~= 0
        error('LPV-MPC optimization failed at step %d with code %d',t,diagnostics.problem);
    end

    u_now = value(U(1));
    u_hist(t) = u_now;
    rho_hist(t) = rho_now;
    
    if mod(t, 5) == 0 || t == 1
        fprintf('Step %d / %d | u_now: %.4f pps | rho: %.4f\n', t, N_total, u_now, rho_now);
    end

    if FAST_SURROGATE_MODE
        % Fast LPV plant update using LUT models based on current Wf
        A_p = LPV2mat(LUT.matrices.A, u_prev);
        B_p = LPV2mat(LUT.matrices.B, u_prev);
        C_p = LPV2mat(LUT.matrices.C, u_prev);
        D_p = LPV2mat(LUT.matrices.D, u_prev);
        u0_p = off2mat(LUT.offset.u, u_prev);
        z0_p = off2mat(LUT.offset.z, u_prev);
        
        y_meas = z0_p(idx_NH) + C_p(idx_NH,:)*x_plant_fast + D_p(idx_NH,:)*(u_prev - u0_p);
        x_plant_fast = A_p*x_plant_fast + B_p*(u_now - u0_p);
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
    end

    y_hist(t) = y_meas;

    x_dev = A*x_dev + B*(u_now-u0_step);
    u_prev = u_now;
    rho_now = u_now;

    if mod(t,20) == 0
        fprintf('LPV closed-loop step %d / %d\n',t,N_sim);
    end
end

time = (0:N_sim-1)'*Ts;
ref_plot = ref(1:N_sim);
controller_name = 'LPV MPC';
save('Results/results_lpv_mpc.mat', 'y_hist', 'u_hist', 'rho_hist', 'solve_times', 'time', 'ref_plot', 'controller_name', 'track_name');
fprintf('Saved results to Results/results_lpv_mpc.mat\n');

if ~FAST_SURROGATE_MODE
    figure('Name','Closed-loop LPV-MPC on nonlinear plant')
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

    disp('Closed-loop LPV-MPC run complete.')
    disp(['Tracking output: ', track_name])
end
disp(['Nominal scheduling point rho_eq = ', num2str(rho_eq)])
