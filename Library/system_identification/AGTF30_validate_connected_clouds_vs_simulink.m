% ========================================================================
% Validate connected-cloud LPV model against nonlinear Simulink plant.
% Compares:
%   1) Smooth adjacent-cloud interpolation model
%   2) Hard-switch cloud model (nearest operating point)
% against AGTF30SysDyn simulation data.
% ========================================================================

clearvars; close all; clc

proj = currentProject;
cd(proj.RootFolder)

%% Settings
Ts = 0.015;
model = 'AGTF30SysDyn';

% Input profile settings
% test_mode options:
%   'nominal_ladder'  - moderate stepped up/down profile (default)
%   'aggressive_steps' - larger/faster steps + stronger jitter
%   'prbs'            - pseudo-random binary sequence style excitation
%   'edge_hammer'     - repeatedly hits near min/max Wf bounds
test_mode = 'nominal_ladder';
rng_seed = 7;

% Common constraints
wf_min = 0.2;
wf_max = 2.2;

% Nominal ladder parameters
dt_set = 0.25;
wf_levels = 0.45:0.2:1.85;
repeat_each = 4;
rand_amp = 0.03;

% Aggressive settings
dt_set_aggr = 0.12;
wf_levels_aggr = 0.25:0.25:2.1;
repeat_each_aggr = 2;
rand_amp_aggr = 0.08;

% PRBS settings
dt_set_prbs = 0.08;
n_prbs = 250;

% Edge-hammer settings
dt_set_edge = 0.10;
n_edge = 140;

% Plot/analysis settings
state_names = {'N_L','N_H','Pt_{21}','Pt_{25}'};
max_plot_abs = 1e6;

%% Load identified tables
id_file = fullfile(proj.RootFolder,'Params','tables','matrices_DE_wf.mat');
wf_file = fullfile(proj.RootFolder,'Params','test_files_fuel_only','test_identification','wf_values.mat');

if ~isfile(id_file)
    error('Missing ID file: %s',id_file);
end
if ~isfile(wf_file)
    error('Missing wf_values file: %s',wf_file);
end

models = load(id_file);
wf_data = load(wf_file);

A_dt = models.A_dt;
B_dt = models.B_dt;
offset = models.offset;
wf_knots = wf_data.Wf_values(:);

n_op = size(A_dt,3);
if numel(wf_knots) ~= n_op
    n_op = min(numel(wf_knots),n_op);
    wf_knots = wf_knots(1:n_op);
    A_dt = A_dt(:,:,1:n_op);
    B_dt = B_dt(:,:,1:n_op);
    offset.x = offset.x(:,1:n_op);
    offset.y = offset.y(:,1:n_op);
    offset.u = offset.u(:,1:n_op);
end

%% Build validation input profile
rng(rng_seed);

switch lower(test_mode)
    case 'nominal_ladder'
        wf_seq = kron([wf_levels, fliplr(wf_levels)],ones(1,repeat_each));
        wf_seq = wf_seq + rand_amp*(2*rand(size(wf_seq)) - 1);
        dt_profile = dt_set;

    case 'aggressive_steps'
        wf_seq = kron([wf_levels_aggr, fliplr(wf_levels_aggr)],ones(1,repeat_each_aggr));
        wf_seq = wf_seq + rand_amp_aggr*(2*rand(size(wf_seq)) - 1);
        dt_profile = dt_set_aggr;

    case 'prbs'
        raw = 2*(rand(1,n_prbs) > 0.5) - 1;
        wf_mid = 0.5*(wf_min + wf_max);
        wf_amp = 0.40*(wf_max - wf_min);
        wf_seq = wf_mid + wf_amp*raw;
        wf_seq = wf_seq + 0.03*(2*rand(size(wf_seq)) - 1);
        dt_profile = dt_set_prbs;

    case 'edge_hammer'
        hi = wf_max - 0.02;
        lo = wf_min + 0.02;
        wf_seq = repmat([lo, hi, lo, hi],1,ceil(n_edge/4));
        wf_seq = wf_seq(1:n_edge);
        wf_seq = wf_seq + 0.02*(2*rand(size(wf_seq)) - 1);
        dt_profile = dt_set_edge;

    otherwise
        error('Unknown test_mode: %s',test_mode);
end

wf_seq = max(min(wf_seq,wf_max),wf_min);

t_in = (0:numel(wf_seq)-1)'*dt_profile;
Tsim = t_in(end);

fprintf('Validation input mode: %s | points: %d | dt: %.4f s\n',...
    test_mode,numel(wf_seq),dt_profile);

%% Configure and run nonlinear Simulink plant
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
MWS.In.Alt  = timeseries([0;0],[0 Tsim],'Name','Alt');
MWS.In.MN   = timeseries([0;0],[0 Tsim],'Name','MN');
MWS.In.dT   = timeseries([0;0],[0 Tsim],'Name','dT');
MWS.In.dVBV = timeseries([0;0],[0 Tsim],'Name','dVBV');
MWS.In.dNz  = timeseries([0;0],[0 Tsim],'Name','dNz');
MWS.In.Wf   = timeseries([wf_seq(1);wf_seq(1)],[0 Tsim],'Name','Wf');
MWS = AGTF30_initial_conditions(MWS);
MWS.In.Wf   = timeseries(wf_seq(:),t_in,'Name','Wf');

simIn = Simulink.SimulationInput(model);
simIn = simIn.setVariable('MWS',MWS);
simIn = simIn.setModelParameter('StartTime','0','StopTime',num2str(Tsim),...
    'ReturnWorkspaceOutputs','on');

simOut = sim(simIn);
out_dyn = fetch_out_dyn(simOut);

%% Extract and resample measured signals
t_raw = signal_time(out_dyn.eng.S0.W);
t = (0:Ts:t_raw(end)).';

u_raw = signal_data(out_dyn.cntrl.Wfact);
psi_raw = [signal_data(out_dyn.eng.Shaft.N_Fan), ...
           signal_data(out_dyn.eng.Shaft.N_HPC), ...
           signal_data(out_dyn.eng.S21.Pt), ...
           signal_data(out_dyn.eng.S25.Pt)];

u = interp1(t_raw,u_raw,t,'linear','extrap');
psi = interp1(t_raw,psi_raw,t,'linear','extrap');

N = numel(t);
nx = size(psi,2);

%% Simulate models (free-run)
psi_smooth = nan(N,nx);
psi_hard = nan(N,nx);
psi_smooth(1,:) = psi(1,:);
psi_hard(1,:) = psi(1,:);

idx_hist = nan(N-1,1);
alpha_hist = nan(N-1,1);

for k = 1:N-1
    rho = u(k);

    % Smooth connected-cloud model
    [A_s,B_s,x0_s,y0_s,u0_s,i0_s,alpha_s] = eval_smooth_ab(rho,wf_knots,A_dt,B_dt,offset);
    psi0_s = [x0_s; y0_s];
    xdev_s = psi_smooth(k,:).' - psi0_s;
    psi_smooth(k+1,:) = (psi0_s + A_s*xdev_s + B_s*(rho-u0_s)).';

    % Hard-switch model (nearest cloud)
    [~,j] = min(abs(wf_knots - rho));
    A_h = A_dt(:,:,j);
    B_h = B_dt(:,:,j);
    psi0_h = [offset.x(:,j); offset.y(:,j)];
    u0_h = offset.u(:,j);
    xdev_h = psi_hard(k,:).' - psi0_h;
    psi_hard(k+1,:) = (psi0_h + A_h*xdev_h + B_h*(rho-u0_h)).';

    idx_hist(k) = i0_s;
    alpha_hist(k) = alpha_s;
end

%% One-step predictions (for discrepancy decomposition)
psi_smooth_1 = nan(N,nx);
psi_hard_1 = nan(N,nx);
psi_smooth_1(1,:) = psi(1,:);
psi_hard_1(1,:) = psi(1,:);

for k = 1:N-1
    rho = u(k);

    [A_s,B_s,x0_s,y0_s,u0_s] = eval_smooth_ab(rho,wf_knots,A_dt,B_dt,offset);
    psi0_s = [x0_s; y0_s];
    xdev_s = psi(k,:).' - psi0_s;
    psi_smooth_1(k+1,:) = (psi0_s + A_s*xdev_s + B_s*(rho-u0_s)).';

    [~,j] = min(abs(wf_knots - rho));
    A_h = A_dt(:,:,j);
    B_h = B_dt(:,:,j);
    psi0_h = [offset.x(:,j); offset.y(:,j)];
    u0_h = offset.u(:,j);
    xdev_h = psi(k,:).' - psi0_h;
    psi_hard_1(k+1,:) = (psi0_h + A_h*xdev_h + B_h*(rho-u0_h)).';
end

% Guard against extreme diverging values in plots
psi_smooth(abs(psi_smooth) > max_plot_abs) = NaN;
psi_hard(abs(psi_hard) > max_plot_abs) = NaN;

%% Metrics
rmse_free_smooth = rmse_cols(psi(:,1:2),psi_smooth(:,1:2));
rmse_free_hard   = rmse_cols(psi(:,1:2),psi_hard(:,1:2));
rmse_step_smooth = rmse_cols(psi(:,1:2),psi_smooth_1(:,1:2));
rmse_step_hard   = rmse_cols(psi(:,1:2),psi_hard_1(:,1:2));

fprintf('\n=== Validation summary: AGTF30SysDyn vs LPV surrogates ===\n');
fprintf('Free-run RMSE N_L [smooth, hard] = [%.4f, %.4f]\n',rmse_free_smooth(1),rmse_free_hard(1));
fprintf('Free-run RMSE N_H [smooth, hard] = [%.4f, %.4f]\n',rmse_free_smooth(2),rmse_free_hard(2));
fprintf('1-step RMSE  N_L [smooth, hard] = [%.4f, %.4f]\n',rmse_step_smooth(1),rmse_step_hard(1));
fprintf('1-step RMSE  N_H [smooth, hard] = [%.4f, %.4f]\n',rmse_step_smooth(2),rmse_step_hard(2));

%% Plot 1: output traces
figure('Name','Validation traces: nonlinear plant vs connected-cloud models')
subplot(3,1,1)
plot(t,psi(:,1),'k','LineWidth',1.2)
hold on
plot(t,psi_smooth(:,1),'b--','LineWidth',1.1)
plot(t,psi_hard(:,1),'r-.','LineWidth',1.1)
grid on
ylabel('N_L [rpm]')
legend('Nonlinear plant','Smooth LPV','Hard-switch LPV','Location','best')

subplot(3,1,2)
plot(t,psi(:,2),'k','LineWidth',1.2)
hold on
plot(t,psi_smooth(:,2),'b--','LineWidth',1.1)
plot(t,psi_hard(:,2),'r-.','LineWidth',1.1)
grid on
ylabel('N_H [rpm]')

subplot(3,1,3)
plot(t,u,'LineWidth',1.2)
grid on
xlabel('Time [s]')
ylabel('W_f [pps]')

%% Plot 2: error comparison
e_free_smooth = psi_smooth(:,1:2) - psi(:,1:2);
e_free_hard = psi_hard(:,1:2) - psi(:,1:2);

e_step_smooth = psi_smooth_1(:,1:2) - psi(:,1:2);
e_step_hard = psi_hard_1(:,1:2) - psi(:,1:2);

figure('Name','Validation errors: free-run vs one-step')
subplot(2,2,1)
plot(t,e_free_smooth(:,1),'b',t,e_free_hard(:,1),'r')
grid on
title('N_L free-run error')
legend('Smooth','Hard','Location','best')

subplot(2,2,2)
plot(t,e_free_smooth(:,2),'b',t,e_free_hard(:,2),'r')
grid on
title('N_H free-run error')

subplot(2,2,3)
plot(t,e_step_smooth(:,1),'b',t,e_step_hard(:,1),'r')
grid on
title('N_L one-step error')
xlabel('Time [s]')

subplot(2,2,4)
plot(t,e_step_smooth(:,2),'b',t,e_step_hard(:,2),'r')
grid on
title('N_H one-step error')
xlabel('Time [s]')

%% Plot 3: blending diagnostics
figure('Name','Connected-cloud blending diagnostics')
subplot(2,1,1)
plot(t(1:end-1),idx_hist,'LineWidth',1.2)
grid on
ylabel('Lower knot index i')
title('Active interval for smooth blending')

subplot(2,1,2)
plot(t(1:end-1),alpha_hist,'LineWidth',1.2)
grid on
xlabel('Time [s]')
ylabel('\alpha')
title('\alpha in M(rho) = (1-\alpha)M_i + \alpha M_{i+1}')

%% Plot 4: RMSE bars
figure('Name','RMSE summary')
bar_data = [rmse_free_smooth(:), rmse_free_hard(:), rmse_step_smooth(:), rmse_step_hard(:)];
bar(categorical({'N_L','N_H'}),bar_data)
grid on
ylabel('RMSE')
legend('Smooth free-run','Hard free-run','Smooth one-step','Hard one-step','Location','best')

disp('Validation complete.')

%% Helpers
function out_dyn = fetch_out_dyn(simOut)
    out_dyn = [];
    if isa(simOut,'Simulink.SimulationOutput')
        try
            out_dyn = simOut.get('out_Dyn');
        catch
            out_dyn = [];
        end
    end
    if isempty(out_dyn)
        try
            out_dyn = evalin('base','out_Dyn');
        catch
            error('Could not retrieve out_Dyn after simulation.');
        end
    end
end

function [A,B,x0,y0,u0,i0,alpha] = eval_smooth_ab(rho,wf_knots,A_all,B_all,offset)
    [i0,i1,alpha] = blend_index(rho,wf_knots);

    A = (1-alpha)*A_all(:,:,i0) + alpha*A_all(:,:,i1);
    B = (1-alpha)*B_all(:,:,i0) + alpha*B_all(:,:,i1);

    x0 = (1-alpha)*offset.x(:,i0) + alpha*offset.x(:,i1);
    y0 = (1-alpha)*offset.y(:,i0) + alpha*offset.y(:,i1);
    u0 = (1-alpha)*offset.u(:,i0) + alpha*offset.u(:,i1);
end

function [i0,i1,alpha] = blend_index(rho,wf_knots)
    if rho <= wf_knots(1)
        i0 = 1; i1 = min(2,numel(wf_knots)); alpha = 0;
        return
    end
    if rho >= wf_knots(end)
        i1 = numel(wf_knots);
        i0 = max(1,i1-1);
        alpha = 1;
        return
    end

    i0 = find(wf_knots <= rho,1,'last');
    i1 = i0 + 1;
    alpha = (rho - wf_knots(i0)) / (wf_knots(i1) - wf_knots(i0));
end

function e = rmse_cols(y, yhat)
    n = size(y,2);
    e = nan(1,n);
    for i = 1:n
        d = yhat(:,i) - y(:,i);
        d = d(isfinite(d));
        if isempty(d)
            e(i) = NaN;
        else
            e(i) = sqrt(mean(d.^2));
        end
    end
end

function d = signal_data(sig)
    if isa(sig,'timeseries')
        d = sig.Data(:);
    elseif isstruct(sig) && isfield(sig,'Data')
        d = sig.Data(:);
    else
        d = sig(:);
    end
    d = double(d);
end

function t = signal_time(sig)
    if isa(sig,'timeseries')
        t = sig.Time(:);
    elseif isstruct(sig) && isfield(sig,'Time')
        t = sig.Time(:);
    else
        error('Signal time vector not found.');
    end
    t = double(t);
end
