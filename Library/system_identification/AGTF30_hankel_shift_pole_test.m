% ========================================================================
% Shifted-Hankel pole test from input-output data (standalone experiment)
%
% Theory used:
%   Build delay-embedded lifted state from measured u,y data,
%   then identify a reduced square operator via SVD + shift relation:
%
%     X_+ ~= A_r X
%     X ~= U_r S_r V_r'
%     A_r = U_r' * X_+ * V_r * inv(S_r)
%
%   Poles are eig(A_r).
%
% This script is for testing only and does not modify the existing workflow.
% ========================================================================

clearvars; close all; clc

proj = currentProject;
cd(proj.RootFolder)

%% User settings
Ts = 0.015;
tstart = 19.5;

lift_rows = 30;      % Hankel block rows for each signal
model_order = 6;     % reduced model order for A_r
signal_name = 'N_HPC'; % 'N_HPC' | 'N_FAN' | 'FNET'

%% Paths
id_file = fullfile(proj.RootFolder,'Params','tables','matrices_DE_wf.mat');
wf_file = fullfile(proj.RootFolder,'Params','test_files_fuel_only','test_identification','wf_values.mat');
data_folder = fullfile(proj.RootFolder,'Params','test_files_fuel_only','test_identification','AGTF30');

if ~isfile(id_file)
    error('Missing ID model file: %s',id_file);
end
if ~isfile(wf_file)
    error('Missing Wf value file: %s',wf_file);
end
if ~isfolder(data_folder)
    error('Missing data folder: %s',data_folder);
end

models = load(id_file);
wf_data = load(wf_file);
wf_values = wf_data.Wf_values(:);
A_id = models.A_dt;

n_id = size(A_id,3);

filelist = dir(fullfile(data_folder,'*.mat'));
if isempty(filelist)
    error('No test files found in %s',data_folder);
end

[~,idx_sort] = sort({filelist.name});
filelist = filelist(idx_sort);

n_used = min([numel(filelist), numel(wf_values), n_id]);
wf_values = wf_values(1:n_used);
A_id = A_id(:,:,1:n_used);
filelist = filelist(1:n_used);

%% Containers
rho_id = nan(n_used,1);
rho_hankel = nan(n_used,1);
stable_id = false(n_used,1);
stable_hankel = false(n_used,1);
ev_id = nan(size(A_id,1),n_used);
ev_hankel = nan(model_order,n_used);
fit_rel_error = nan(n_used,1);

for k = 1:n_used
    % Baseline ID poles
    ev_id(:,k) = eig(A_id(:,:,k));
    rho_id(k) = max(abs(ev_id(:,k)));
    stable_id(k) = rho_id(k) < 1;

    % Load measured cloud
    sample = load(fullfile(filelist(k).folder,filelist(k).name));
    out_dyn = get_out_dyn(sample);

    t_raw = signal_time(out_dyn.eng.S0.W);
    u_raw = signal_data(out_dyn.cntrl.Wfact);
    y_raw = pick_output(out_dyn, signal_name);

    t = (tstart:Ts:t_raw(end)).';
    u = interp1(t_raw,u_raw,t,'linear','extrap');
    y = interp1(t_raw,y_raw,t,'linear','extrap');

    % Remove offsets so lifted data focuses on local dynamics
    u = u - mean(u);
    y = y - mean(y);

    % Build delay embeddings (shifted-Hankel style)
    [X, Xp] = build_shifted_lifted_data(u, y, lift_rows);

    if isempty(X)
        warning('Not enough data at index %d; skipping.',k);
        continue
    end

    r = min([model_order, size(X,1), size(X,2)]);
    [U,S,V] = svd(X,'econ');
    Ur = U(:,1:r);
    Sr = S(1:r,1:r);
    Vr = V(:,1:r);

    Ar = Ur' * Xp * Vr / Sr;

    ev = eig(Ar);
    if numel(ev) < model_order
        ev = [ev; nan(model_order-numel(ev),1)];
    end
    ev_hankel(:,k) = ev;

    rho_hankel(k) = max(abs(ev(~isnan(ev))));
    stable_hankel(k) = rho_hankel(k) < 1;

    Xp_hat = Ar * (Ur' * X);
    Xp_proj = Ur * Xp_hat;
    fit_rel_error(k) = norm(Xp - Xp_proj,'fro') / max(norm(Xp,'fro'),1e-12);
end

fprintf('\n=== Shifted-Hankel Pole Test Summary ===\n');
fprintf('ID stable points: %d / %d\n',sum(stable_id),n_used);
fprintf('Hankel-subspace stable points: %d / %d\n',sum(stable_hankel),n_used);
if any(isfinite(fit_rel_error))
    fprintf('Median lifted-shift relative fit error: %.4e\n',median(fit_rel_error(isfinite(fit_rel_error))));
end

%% Plot 1: eigenvalue plane
th = linspace(0,2*pi,400);
figure('Name','Shifted-Hankel poles vs SYS-ID poles')
plot(cos(th),sin(th),'k--','LineWidth',1.0)
hold on
axis equal
grid on
xlabel('Re(\lambda)')
ylabel('Im(\lambda)')
title(sprintf('Pole map from shifted-Hankel (%s) vs SYS-ID',signal_name))

cmap = turbo(n_used);
for k = 1:n_used
    scatter(real(ev_id(:,k)),imag(ev_id(:,k)),24,cmap(k,:),'filled','o')
    this_ev = ev_hankel(:,k);
    this_ev = this_ev(~isnan(this_ev));
    scatter(real(this_ev),imag(this_ev),24,cmap(k,:),'x')
end
legend('Unit circle','SYS-ID poles','Shifted-Hankel poles','Location','best')

%% Plot 2: spectral radius trends
figure('Name','Spectral radius: SYS-ID vs shifted-Hankel')
plot(wf_values,rho_id,'-o','LineWidth',1.3,'DisplayName','SYS-ID')
hold on
plot(wf_values,rho_hankel,'-x','LineWidth',1.3,'DisplayName','Shifted-Hankel')
yline(1,'r--','|\lambda| = 1','LineWidth',1.1)
grid on
xlabel('W_f operating point [pps]')
ylabel('Spectral radius')
legend('Location','best')
title('Stability trend comparison')

%% Plot 3: identification quality indicator
figure('Name','Shifted relation fit quality')
plot(wf_values,fit_rel_error,'-s','LineWidth',1.3)
grid on
xlabel('W_f operating point [pps]')
ylabel('Relative error ||X_+ - X_+^{hat}||_F / ||X_+||_F')
title('Lifted shift-fit quality (lower is better)')

%% Helpers
function [X, Xp] = build_shifted_lifted_data(u, y, L)
    u = u(:);
    y = y(:);
    N = numel(u);

    n_cols = N - L;
    if n_cols < 2
        X = [];
        Xp = [];
        return
    end

    Hu = block_hankel(u, L); % L x (N-L+1)
    Hy = block_hankel(y, L); % L x (N-L+1)

    Xi = [Hu; Hy];
    X = Xi(:,1:n_cols);
    Xp = Xi(:,2:n_cols+1);
end

function H = block_hankel(sig, L)
    sig = sig(:);
    n_cols = numel(sig) - L + 1;
    H = zeros(L,n_cols);
    for i = 1:L
        H(i,:) = sig(i:i+n_cols-1).';
    end
end

function out_dyn = get_out_dyn(data)
    if isstruct(data) && isfield(data,'out_Dyn')
        out_dyn = data.out_Dyn;
        return
    end
    error('Could not find out_Dyn in loaded data file.');
end

function y = pick_output(out_dyn, name)
    switch upper(name)
        case 'N_HPC'
            y = signal_data(out_dyn.eng.Shaft.N_HPC);
        case 'N_FAN'
            y = signal_data(out_dyn.eng.Shaft.N_Fan);
        case 'FNET'
            y = signal_data(out_dyn.eng.Perf.Fnet);
        otherwise
            error('Unsupported signal_name: %s',name);
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
