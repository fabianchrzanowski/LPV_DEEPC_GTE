clearvars; close all; clc
yalmip('clear');

% takes forever, doesnt work that great

solver_name = 'mosek'; % Change to 'quadprog' if you have a QP solver available

proj = currentProject;
cd(proj.RootFolder)

training_folder = fullfile(proj.RootFolder,'Params','test_files_fuel_only','test_identification','AGTF30');
if ~isfolder(training_folder)
    error('DeePC training folder not found: %s',training_folder);
end

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
    idx_track = find(output_names_lower == "n_l",1);
end
if isempty(idx_track)
    idx_track = 1;
end
track_name = sys_eng.OutputName{idx_track};

Ts = 0.015;
T_ini = 20;
N_pred = 20;
N_sim = 160;
N_total = T_ini + N_sim;

u_min = 0.2;
u_max = 2.2;
du_max = 0.04;

q_track = 1;
r_u = 0.05;
r_du = 25;
lambda_g = 1e-4;
lambda_sigma = 1e4;

track0 = z0(idx_track);
if contains(lower(track_name),'n_h') || contains(lower(track_name),'n_l')
    amp = 300;
else
    amp = max(0.02*abs(track0),1e-3);
end
ref = track0 + amp*sin(2*pi*(0:N_total+N_pred-1)'/140);

[U_hankel, Y_hankel] = build_deepc_hankel(training_folder, track_name, T_ini + N_pred);
max_hankel_cols = 800;
[U_hankel, Y_hankel] = reduce_hankel_columns(U_hankel, Y_hankel, max_hankel_cols);
U_p = U_hankel(1:T_ini,:);
U_f = U_hankel(T_ini+1:end,:);
Y_p = Y_hankel(1:T_ini,:);
Y_f = Y_hankel(T_ini+1:end,:);

if size(U_hankel,2) < 1
    error('Not enough training data to build a DeePC Hankel matrix.');
end

fprintf('DeePC Hankel size: L=%d, columns=%d\n',size(U_hankel,1),size(U_hankel,2));

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
u_hist = zeros(N_total,1);
y_hist = nan(N_total,1);
u_hist(1:T_ini) = u0;

u_prev = u0;
ops = sdpsettings('solver',solver_name,'verbose',1,'usex0',1);
if strcmpi(solver_name,'fmincon')
    ops = sdpsettings(ops,...
    'fmincon.algorithm','interior-point',...
    'fmincon.maxiter',150,...
    'fmincon.maxfunevals',20000,...
    'fmincon.tolfun',1e-6,...
    'fmincon.tolcon',1e-6,...
    'fmincon.display','off');
end

for t = 1:T_ini
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

for t = T_ini+1:N_total
    u_ini = u_hist(t-T_ini:t-1);
    y_ini = y_hist(t-T_ini:t-1);
    ref_seg = ref(t:t+N_pred-1);

    g = sdpvar(size(U_hankel,2),1,'full');
    U = sdpvar(N_pred,1,'full');
    Y = sdpvar(N_pred,1,'full');
    sigma = sdpvar(T_ini,1,'full');

    constraints = [];
    objective = lambda_g*sum_square(g) + lambda_sigma*sum_square(sigma);

    constraints = [constraints, U_p*g == u_ini];
    constraints = [constraints, Y_p*g == y_ini + sigma];
    constraints = [constraints, U_f*g == U];
    constraints = [constraints, Y_f*g == Y];

    for i = 1:N_pred
        constraints = [constraints, u_min <= U(i), U(i) <= u_max];

        if i == 1
            du_i = U(i) - u_prev;
        else
            du_i = U(i) - U(i-1);
        end

        constraints = [constraints, -du_max <= du_i, du_i <= du_max];
        objective = objective + q_track*(Y(i) - ref_seg(i))^2 + r_u*(U(i) - u0)^2 + r_du*(du_i)^2;
    end

    if strcmpi(solver_name,'fmincon')
    assign(g, zeros(size(U_hankel,2),1));
    assign(U, u_prev*ones(N_pred,1));
    assign(Y, y_ini(end)*ones(N_pred,1));
    assign(sigma, zeros(T_ini,1));
    end

    tic;
    diagnostics = optimize(constraints,objective,ops);
    solve_time = toc;
    if diagnostics.problem ~= 0
        error('DeePC optimization failed at step %d with code %d',t - T_ini,diagnostics.problem);
    end

    u_now = value(U(1));
    u_hist(t) = u_now;

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

    u_prev = u_now;

    if mod(t - T_ini,20) == 0
        fprintf('DeePC closed-loop step %d / %d (solve %.3f s)\n',t - T_ini,N_sim,solve_time);
    end
end

time = (0:N_total-1)'*Ts;
control_time = time(T_ini+1:end);
ref_plot = ref(1:N_total);

%% plotting
figure('Name','DeePC closed-loop on nonlinear plant')
subplot(2,1,1)
plot(control_time,ref_plot(T_ini+1:N_total),'--','LineWidth',1.2)
hold on
plot(control_time,y_hist(T_ini+1:end),'LineWidth',1.4)
grid on
xlabel('Time [s]')
ylabel(track_name,'Interpreter','tex')
legend('Reference','AGTF30 measured output','Location','best')

subplot(2,1,2)
plot(control_time,u_hist(T_ini+1:end),'LineWidth',1.4)
grid on
xlabel('Time [s]')
ylabel('W_f [pps]','Interpreter','none')
title('Applied fuel command')

disp('DeePC closed-loop run complete.')
disp(['Tracking output: ', track_name])
disp(['Operating point index: ', num2str(k), ', Wf_eq = ', num2str(wf_values(k))])


%% Functions
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