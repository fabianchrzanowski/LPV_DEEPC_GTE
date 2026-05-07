if ~exist('FAST_SURROGATE_MODE', 'var')
    clearvars; close all; clc
    yalmip('clear');
    run('setup_simulation_params.m');
end

solver_name = 'mosek';

proj = currentProject;
cd(proj.RootFolder)

% Variables now sourced from setup_simulation_params.m

track_name = models.sys_eng.OutputName{idx_NH};

A_all = models.A_dt;
B_all = models.B_dt;
C_all = models.C_dt;
D_all = models.D_dt;
offset = models.offset;

% Extract the specific Cruise state space model (for Frozen MPC)
wf_frozen = 1.2;
[~, k_frozen] = min(abs(wf_values - wf_frozen));
A = A_all(:,:,k_frozen);
B = B_all(:,:,k_frozen);
C = C_all(:,:,k_frozen);
D = D_all(:,:,k_frozen);
u0 = offset.u(:, k_frozen);
z0 = offset.z(:, k_frozen);



u_min = 0.2;
u_max = 2.2;
du_max = 0.04;

q_track = 1;
r_du = 25;

if ~FAST_SURROGATE_MODE
    % The trajectory ref is now sourced from setup_simulation_params.m
end

x_dev = zeros(size(A,1),1);
if FAST_SURROGATE_MODE
    x_plant_fast = zeros(size(A,1),1);
end
u_prev = u0;

u_hist = zeros(N_sim,1);
y_hist = nan(N_sim,1);

ops = sdpsettings('solver',solver_name,'verbose',1);

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
    U = sdpvar(1,N_pred,'full');
    x_pred = x_dev;
    t
    constraints = [];
    objective = 0;

    for i = 1:N_pred
        x_pred = A*x_pred + B*(U(i)-u0);
        z_pred = z0 + C*x_pred + D*(U(i)-u0);

        constraints = [constraints, u_min <= U(i), U(i) <= u_max];

        if i == 1
            du_i = U(i) - u_prev;
        else
            du_i = U(i) - U(i-1);
        end

        constraints = [constraints, -du_max <= du_i, du_i <= du_max];
        objective = objective + q_track*(z_pred(idx_track) - ref(t+i-1))^2 + r_du*(du_i)^2;
    end

    diagnostics = optimize(constraints,objective,ops);
    if diagnostics.problem ~= 0
        error('MPC optimization failed at step %d with code %d',t,diagnostics.problem);
    end

    u_now = value(U(1));
    u_hist(t) = u_now;
    
    if mod(t, 5) == 0 || t == 1
        fprintf('Step %d / %d | u_now: %.4f pps\n', t, N_total, u_now);
    end

    if FAST_SURROGATE_MODE
        [~, k_plant] = min(abs(wf_values - u_prev));
        A_p = A_all(:,:,k_plant); B_p = B_all(:,:,k_plant);
        C_p = C_all(:,:,k_plant); D_p = D_all(:,:,k_plant);
        u0_p = offset.u(:,k_plant); z0_p = offset.z(:,k_plant);
        
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
            y_meas = z0(idx_NH) + C(idx_NH,:)*x_dev + D(idx_NH,:)*(u_now-u0);
        end
    end

    y_hist(t) = y_meas;

    x_dev = A*x_dev + B*(u_now-u0);
    u_prev = u_now;

    if mod(t,20) == 0
        fprintf('closed-loop step %d / %d\n',t,N_sim);
    end
end

time = (0:N_sim-1)'*Ts;
ref_plot = ref(1:N_sim);
controller_name = 'Frozen MPC';
save('Results/results_mpc.mat', 'y_hist', 'u_hist', 'solve_times', 'time', 'ref_plot', 'controller_name', 'track_name');
fprintf('Saved results to Results/results_mpc.mat\n');

if ~FAST_SURROGATE_MODE
    figure('Name','Closed-loop MPC on nonlinear plant')
    subplot(2,1,1)
    plot(time,ref_plot,'--','LineWidth',1.2)
    hold on
    plot(time,y_hist,'LineWidth',1.4)
    grid on
    xlabel('Time [s]')
    ylabel(track_name,'Interpreter','tex')
    legend('Reference','AGTF30 measured output','Location','best')

    subplot(2,1,2)
    plot(time,u_hist,'LineWidth',1.4)
    grid on
    xlabel('Time [s]')
    ylabel('W_f [pps]','Interpreter','none')
    title('Applied fuel command')

    disp('Closed-loop run complete.')
    disp(['Tracking output: ', track_name])
end
