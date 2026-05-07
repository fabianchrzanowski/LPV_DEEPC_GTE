clearvars; close all; clc
yalmip('clear');

solver_name = 'mosek';


proj = currentProject;
cd(proj.RootFolder)

% Load the data from the System ID
models = load(fullfile(proj.RootFolder,'Params','tables','matrices_DE_wf.mat'));
wf_data = load(fullfile(proj.RootFolder,'Params','test_files_fuel_only','test_identification','wf_values.mat'));

A_all = models.A_dt;
B_all = models.B_dt;
C_all = models.C_dt;
D_all = models.D_dt;
offset = models.offset;
sys_eng = models.sys_eng;
wf_values = wf_data.Wf_values(:)';

wf_eq = 1.2;
[~,k] = min(abs(wf_values - wf_eq));

A = A_all(:,:,k);
B = B_all(:,:,k);
C = C_all(:,:,k);
D = D_all(:,:,k);

x0 = [offset.x(:,k); offset.y(:,k)];
u0 = offset.u(:,k);
z0 = offset.z(:,k);

output_names_lower = lower(string(sys_eng.OutputName));
idx_track = find(output_names_lower == "thrust",1);
if isempty(idx_track)
    idx_track = find(output_names_lower == "fnet",1);
end
if isempty(idx_track)
    idx_track = find(output_names_lower == "n_h",1);
end
if isempty(idx_track)
    idx_track = 1;
end

track_name = sys_eng.OutputName{idx_track};

Ts = 0.015;
N_pred = 20;
N_sim = 350;

u_min = 0.2;
u_max = 2.2;
du_max = 0.04;

q_track = 1;
r_du = 25;

track0 = z0(idx_track);
if contains(lower(track_name),'n_h') || contains(lower(track_name),'n_l')
    amp = 300;
else
    amp = max(0.02*abs(track0),1e-3);
end

ref = track0 + amp*sin(2*pi*(0:N_sim+N_pred-1)'/140);

x_dev = zeros(size(A,1),1);
u_prev = u0;

u_hist = zeros(N_sim,1);
z_track_hist = zeros(N_sim,1);


ops = sdpsettings('solver',solver_name,'verbose',1);

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
        error('MPC optimization failed at step %d with code %d. Try another solver in solver_name.',t,diagnostics.problem);
    end

    u_now = value(U(1));

    x_dev = A*x_dev + B*(u_now-u0);
    z_now = z0 + C*x_dev + D*(u_now-u0);

    u_hist(t) = u_now;
    z_track_hist(t) = z_now(idx_track);
    u_prev = u_now;
end

time = (0:N_sim-1)'*Ts;

figure('Name','MPC tracking')
subplot(2,1,1)
plot(time,ref(1:N_sim),'--','LineWidth',1.2)
hold on
plot(time,z_track_hist,'LineWidth',1.4)
grid on
xlabel('Time [s]')
ylabel(track_name,'Interpreter','tex')
legend('Reference','Tracked output','Location','best')

subplot(2,1,2)
plot(time,u_hist,'LineWidth',1.4)
grid on
xlabel('Time [s]')
ylabel('W_f [pps]','Interpreter','none')
title('Fuel flow command')

%% Replay on nonlinear Simulink plant (AGTF30SysDyn)
run_nonlinear_test = true;
if run_nonlinear_test
    model = 'AGTF30SysDyn';
    Tsim = time(end);

    MWS.engName = 'AGTF30';
    MWS = setup_Controller(MWS);
    MWS = setup_AllEng(MWS);

    busVars = load(fullfile(proj.RootFolder,'Params','Eng_Bus.mat'));
    busVarNames = fieldnames(busVars);
    for ib = 1:numel(busVarNames)
        assignin('base',busVarNames{ib},busVars.(busVarNames{ib}));
    end

    MWS.In.Ts = Ts;
    MWS.In.Tsim = Tsim;

    MWS.In.Alt = timeseries([0;0],[0 Tsim],'Name','Alt');
    MWS.In.MN  = timeseries([0;0],[0 Tsim],'Name','MN');
    MWS.In.dT  = timeseries([0;0],[0 Tsim],'Name','dT');

    MWS.In.Wf   = timeseries([u_hist(1);u_hist(1)],[0 Tsim],'Name','Wf');
    MWS.In.dVBV = timeseries([0;0],[0 Tsim],'Name','dVBV');
    MWS.In.dNz  = timeseries([0;0],[0 Tsim],'Name','dNz');
    MWS = AGTF30_initial_conditions(MWS);

    MWS.In.Wf = timeseries(u_hist,time,'Name','Wf');

    simOut = sim(model);
    out_Dyn_local = [];
    if isa(simOut,'Simulink.SimulationOutput')
        try
            out_Dyn_local = simOut.get('out_Dyn');
        catch
            out_Dyn_local = [];
        end
    elseif isstruct(simOut) && isfield(simOut,'eng')
        out_Dyn_local = simOut;
    end

    if isempty(out_Dyn_local)
        try
            out_Dyn_local = evalin('base','out_Dyn');
        catch
            error('Could not retrieve out_Dyn from simulation output or base workspace.');
        end
    end

    track_name_lower = lower(track_name);
    if contains(track_name_lower,'n_h')
        y_plant = out_Dyn_local.eng.Shaft.N_HPC.Data;
        t_plant = out_Dyn_local.eng.Shaft.N_HPC.Time;
    elseif contains(track_name_lower,'n_l')
        y_plant = out_Dyn_local.eng.Shaft.N_Fan.Data;
        t_plant = out_Dyn_local.eng.Shaft.N_Fan.Time;
    else
        y_plant = [];
        t_plant = [];
        warning('Nonlinear replay currently supports N_H/N_L tracking extraction. Add thrust signal path to extend.');
    end

    if ~isempty(y_plant)
        ref_plant = interp1(time,ref(1:N_sim),t_plant,'linear','extrap');
        figure('Name','Nonlinear plant replay')
        subplot(2,1,1)
        plot(t_plant,ref_plant,'--','LineWidth',1.2)
        hold on
        plot(t_plant,y_plant,'LineWidth',1.3)
        grid on
        xlabel('Time [s]')
        ylabel(track_name,'Interpreter','tex')
        legend('Reference','AGTF30 plant','Location','best')

        subplot(2,1,2)
        plot(time,u_hist,'LineWidth',1.3)
        grid on
        xlabel('Time [s]')
        ylabel('W_f [pps]','Interpreter','none')
        title('Applied fuel command')
    end
end

disp(['Tracking output: ', track_name])
disp(['Operating point index: ', num2str(k), ', Wf_eq = ', num2str(wf_values(k))])
