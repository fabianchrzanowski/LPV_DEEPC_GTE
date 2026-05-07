% ========================================================================
% Compare linearizations from system-ID tables vs direct data-cloud fitting
% and inspect pole motion with stability-constrained pole projection.
% ========================================================================

clearvars; close all; clc

proj = currentProject;
cd(proj.RootFolder)

Ts = 0.015;
tstart = 19.5;
ridge_lambda = 1e-6;
rho_margin = 0.98; % keep projected poles inside this radius
beta_lin = 0.5;    % pole interpolation factor between adjacent SYS-ID clouds (0..1)
n_show_clouds = 8;  % number of representative clouds for summary system map plot
max_plot_abs = 1e6; % reject non-physical one-step predictions in summary map

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
B_id = models.B_dt;

n_op = size(A_id,3);
nx = size(A_id,1);

if numel(wf_values) ~= n_op
    warning('wf_values length (%d) and number of A_dt slices (%d) differ. Using min length.',...
        numel(wf_values),n_op);
    n_op = min(numel(wf_values),n_op);
    wf_values = wf_values(1:n_op);
    A_id = A_id(:,:,1:n_op);
    B_id = B_id(:,:,1:n_op);
end

filelist = dir(fullfile(data_folder,'*.mat'));
if isempty(filelist)
    error('No .mat data files found in %s',data_folder);
end

% Sort files by name to match sys_id save order.
[~,idx_sort] = sort({filelist.name});
filelist = filelist(idx_sort);

if numel(filelist) < n_op
    warning('Only %d data files found for %d operating points. Truncating.',numel(filelist),n_op);
    n_op = numel(filelist);
    wf_values = wf_values(1:n_op);
    A_id = A_id(:,:,1:n_op);
    B_id = B_id(:,:,1:n_op);
end

A_data = nan(nx,nx,n_op);
B_data = nan(nx,1,n_op);
rho_id = nan(n_op,1);
rho_data = nan(n_op,1);
stable_id = false(n_op,1);
stable_data = false(n_op,1);
ev_id = nan(nx,n_op);
ev_data = nan(nx,n_op);
fit_rmse = nan(n_op,1);
ev_lin = nan(nx,n_op-1);

for k = 1:n_op
    rho_id(k) = max(abs(eig(A_id(:,:,k))));
    stable_id(k) = rho_id(k) < 1;
    ev_id(:,k) = eig(A_id(:,:,k));

    sample = load(fullfile(filelist(k).folder,filelist(k).name));
    out_dyn = get_out_dyn(sample);

    [psi, u] = extract_cloud_signals(out_dyn, Ts, tstart);
    [A_k, B_k, fit_rmse(k)] = fit_local_linear_model(psi, u, ridge_lambda);

    A_data(:,:,k) = A_k;
    B_data(:,:,k) = B_k;
    ev_data(:,k) = eig(A_k);
    rho_data(k) = max(abs(ev_data(:,k)));
    stable_data(k) = rho_data(k) < 1;
end

fprintf('\n=== Stability summary (discrete-time) ===\n');
fprintf('ID models stable points:   %d / %d\n',sum(stable_id),n_op);
fprintf('Data-fit models stable points: %d / %d\n',sum(stable_data),n_op);

pair_score = 0.5*(fit_rmse(1:end-1) + fit_rmse(2:end));
[~,k_pair] = max(pair_score);
pair_idx = [k_pair, k_pair+1];

fprintf('Most nonlinear adjacent pair (RMSE score): k = [%d, %d], Wf = [%.3f, %.3f] pps\n',...
    pair_idx(1),pair_idx(2),wf_values(pair_idx(1)),wf_values(pair_idx(2)));

for k = 1:n_op-1
    ev_next_aligned = align_poles(ev_id(:,k),ev_id(:,k+1));
    ev_lin(:,k) = (1 - beta_lin)*ev_id(:,k) + beta_lin*ev_next_aligned;
end

rho_lin = max(abs(ev_lin),[],1).';
stable_lin = rho_lin < 1;
fprintf('Eigen-linearized stable points: %d / %d\n',sum(stable_lin),n_op-1);

%% Plot 1: Pole clouds in eigenvalue plane
figure('Name','Eigenvalue plane: ID vs data-fit')
th = linspace(0,2*pi,400);
plot(cos(th),sin(th),'k--','LineWidth',1.0)
hold on
axis equal
grid on
xlabel('Re(\lambda)')
ylabel('Im(\lambda)')
title('Discrete-time poles across operating points')

cmap = parula(n_op);
for k = 1:n_op
    scatter(real(ev_id(:,k)),imag(ev_id(:,k)),26,cmap(k,:),'filled','o')
    scatter(real(ev_data(:,k)),imag(ev_data(:,k)),26,cmap(k,:),'x')
end
legend('Unit circle','ID poles','Data-fit poles','Location','best')

%% Plot 1b: all poles (ID, data-fit, eigen-linearized from ID)
figure('Name','All poles: ID vs data-fit vs eigen-linearized(ID)')
plot(cos(th),sin(th),'k--','LineWidth',1.0)
hold on
axis equal
grid on
xlabel('Re(\lambda)')
ylabel('Im(\lambda)')
title(sprintf('All poles in eigenvalue plane (ID interpolation, beta = %.2f)',beta_lin))

scatter(real(ev_id(:)),imag(ev_id(:)),20,[0 0.45 0.74],'o','DisplayName','ID poles')
scatter(real(ev_data(:)),imag(ev_data(:)),20,[0.85 0.33 0.1],'x','DisplayName','Data-fit poles')
scatter(real(ev_lin(:)),imag(ev_lin(:)),26,[0.1 0.65 0.1],'d','DisplayName','Eigen-linearized poles (from ID)')
legend('Unit circle','ID poles','Data-fit poles','Eigen-linearized poles (from ID)','Location','best')

%% Plot 2: Spectral radius vs operating point
figure('Name','Spectral radius vs Wf')
plot(wf_values,rho_id,'-o','LineWidth',1.3,'DisplayName','ID A_{dt}')
hold on
plot(wf_values,rho_data,'-x','LineWidth',1.3,'DisplayName','Data-fit A')
yline(1,'r--','|\lambda| = 1 boundary','LineWidth',1.1)
grid on
xlabel('W_f operating point [pps]')
ylabel('Spectral radius')
legend('Location','best')
title('Pole radius movement over operating points')

%% Plot 3: nonlinearity score over operating points
figure('Name','Data-cloud linearization mismatch (RMSE)')
plot(wf_values,fit_rmse,'-o','LineWidth',1.3)
hold on
plot(wf_values(pair_idx),fit_rmse(pair_idx),'rs','MarkerSize',8,'LineWidth',1.2)
grid on
xlabel('W_f operating point [pps]')
ylabel('One-step fit RMSE (state-space)')
title('Cloud nonlinearity indicator from local linear fit residual')
legend('RMSE per cloud','Selected adjacent pair','Location','best')

%% Plot 4: summary system map (similar spirit to sys_id_1)
idx_show = unique(round(linspace(1,n_op,min(n_show_clouds,n_op))));
fig_map = figure('Name','System map summary: measured vs SYS-ID vs data-fit vs eigen-linearized');
tiledlayout(fig_map,1,2,'Padding','compact','TileSpacing','compact');

ax1 = nexttile;
hold(ax1,'on'); grid(ax1,'on');
xlabel(ax1,'W_f [pps]');
ylabel(ax1,'N_L [rpm]');
title(ax1,'N_L map across representative clouds');

ax2 = nexttile;
hold(ax2,'on'); grid(ax2,'on');
xlabel(ax2,'W_f [pps]');
ylabel(ax2,'N_H [rpm]');
title(ax2,'N_H map across representative clouds');

cmap_show = turbo(numel(idx_show));
for ii_show = 1:numel(idx_show)
    k = idx_show(ii_show);

    sample_k = load(fullfile(filelist(k).folder,filelist(k).name));
    out_k = get_out_dyn(sample_k);
    [psi_k, u_k] = extract_cloud_signals(out_k, Ts, tstart);

    % Measured cloud
    plot(ax1,u_k,psi_k(:,1),'-','Color',0.55*cmap_show(ii_show,:),'LineWidth',1.2);
    plot(ax2,u_k,psi_k(:,2),'-','Color',0.55*cmap_show(ii_show,:),'LineWidth',1.2);

    % One-step predictions are used here (instead of free-run) to avoid
    % unstable local models blowing up the summary map scale.
    psi0_id = [models.offset.x(:,k); models.offset.y(:,k)].';
    u0_id = models.offset.u(:,k);
    psi_id_abs = one_step_predict_id(A_id(:,:,k), B_id(:,:,k), psi_k, u_k, psi0_id, u0_id);

    psi_df_abs = one_step_predict_absolute(A_data(:,:,k), B_data(:,:,k), psi_k, u_k);

    % Eigen-linearized model from adjacent SYS-ID poles
    if k < n_op
        k_nei = k + 1;
    else
        k_nei = k - 1;
    end
    lambda_id_k = eig(A_id(:,:,k));
    lambda_id_nei = align_poles(lambda_id_k, eig(A_id(:,:,k_nei)));
    lambda_lin_k = (1 - beta_lin)*lambda_id_k + beta_lin*lambda_id_nei;
    A_eig_lin_k = impose_eigenvalues(A_id(:,:,k), lambda_lin_k);
    B_eig_lin_k = (1 - beta_lin)*B_id(:,:,k) + beta_lin*B_id(:,:,k_nei);
    psi_el_abs = one_step_predict_id(A_eig_lin_k, B_eig_lin_k, psi_k, u_k, psi0_id, u0_id);

    psi_id_abs(abs(psi_id_abs) > max_plot_abs) = NaN;
    psi_df_abs(abs(psi_df_abs) > max_plot_abs) = NaN;
    psi_el_abs(abs(psi_el_abs) > max_plot_abs) = NaN;

    % Overlay model trajectories
    plot(ax1,u_k,psi_id_abs(:,1),'g-','LineWidth',0.9);
    plot(ax2,u_k,psi_id_abs(:,2),'g-','LineWidth',0.9);

    plot(ax1,u_k,psi_df_abs(:,1),'b--','LineWidth',1.0);
    plot(ax2,u_k,psi_df_abs(:,2),'b--','LineWidth',1.0);

    plot(ax1,u_k,psi_el_abs(:,1),'m-.','LineWidth',1.0);
    plot(ax2,u_k,psi_el_abs(:,2),'m-.','LineWidth',1.0);
end

% Legend handles (dummy lines for clean legend)
plot(ax1,nan,nan,'k-','LineWidth',1.2,'DisplayName','Measured cloud')
plot(ax1,nan,nan,'g-','LineWidth',1.0,'DisplayName','SYS-ID model')
plot(ax1,nan,nan,'b--','LineWidth',1.0,'DisplayName','Data-fit model')
plot(ax1,nan,nan,'m-.','LineWidth',1.0,'DisplayName','Eigen-linearized model (from ID)')
legend(ax1,'Location','best')

%% Pole projection on the most nonlinear adjacent pair + system-plane impact
fig_pair = figure('Name','Most nonlinear adjacent pair: pole and system-plane effect (ID interpolation)');
tiledlayout(fig_pair,2,2,'Padding','compact','TileSpacing','compact');

for ii = 1:2
    k_demo = pair_idx(ii);
    A_demo = A_data(:,:,k_demo);
    B_demo = B_data(:,:,k_demo);
    lambda_demo = eig(A_demo);
    rho_demo = max(abs(lambda_demo));

    A_id_demo = A_id(:,:,k_demo);
    B_id_demo = B_id(:,:,k_demo);

    if k_demo < n_op
        k_nei = k_demo + 1;
    else
        k_nei = k_demo - 1;
    end
    lambda_id_demo = eig(A_id_demo);
    lambda_id_nei_aligned = align_poles(lambda_id_demo, eig(A_id(:,:,k_nei)));
    lambda_eig_lin = (1 - beta_lin)*lambda_id_demo + beta_lin*lambda_id_nei_aligned;
    A_demo_eig_lin = impose_eigenvalues(A_id_demo, lambda_eig_lin);
    B_demo_eig_lin = (1 - beta_lin)*B_id_demo + beta_lin*B_id(:,:,k_nei);

    [A_demo_proj, alpha] = project_to_stable_disk(A_demo, rho_margin);
    B_demo_proj = alpha*B_demo;
    lambda_demo_proj = eig(A_demo_proj);

    sample_demo = load(fullfile(filelist(k_demo).folder,filelist(k_demo).name));
    out_demo = get_out_dyn(sample_demo);
    [psi_demo, u_demo] = extract_cloud_signals(out_demo, Ts, tstart);

    u_dev = u_demo - u_demo(1);
    psi_meas_dev = psi_demo - psi_demo(1,:);

    psi_id = simulate_deviation(A_id_demo, B_id_demo, u_dev);
    psi_lin = simulate_deviation(A_demo, B_demo, u_dev);
    psi_eig_lin = simulate_deviation(A_demo_eig_lin, B_demo_eig_lin, u_dev);
    psi_proj = simulate_deviation(A_demo_proj, B_demo_proj, u_dev);

    nexttile
    plot(real(lambda_demo),imag(lambda_demo),'bx','MarkerSize',8,'LineWidth',1.4)
    hold on
    plot(real(lambda_id_demo),imag(lambda_id_demo),'gs','MarkerSize',6,'LineWidth',1.2)
    plot(real(lambda_eig_lin),imag(lambda_eig_lin),'md','MarkerSize',6,'LineWidth',1.2)
    plot(real(lambda_demo_proj),imag(lambda_demo_proj),'ro','MarkerSize',6,'LineWidth',1.2)
    plot(cos(th),sin(th),'k--','LineWidth',1.0)
    axis equal
    grid on
    xlabel('Re(\lambda)')
    ylabel('Im(\lambda)')
    title(sprintf('Pole map, W_f = %.3f pps',wf_values(k_demo)))
    legend('Data-fit','ID','Eigen-linearized (from ID)','Projected stable (data-fit)','Unit circle','Location','best')

    nexttile
    plot(psi_meas_dev(:,1),psi_meas_dev(:,2),'k','LineWidth',1.2,'DisplayName','Measured cloud')
    hold on
    plot(psi_id(:,1),psi_id(:,2),'g-','LineWidth',1.1,'DisplayName','ID model')
    plot(psi_lin(:,1),psi_lin(:,2),'b--','LineWidth',1.2,'DisplayName','Data-fit model')
    plot(psi_eig_lin(:,1),psi_eig_lin(:,2),'m-.','LineWidth',1.2,'DisplayName','Eigen-linearized model (from ID)')
    plot(psi_proj(:,1),psi_proj(:,2),'r-.','LineWidth',1.2,'DisplayName','Projected-pole model')
    grid on
    xlabel('\Delta N_L')
    ylabel('\Delta N_H')
    title('System plane')
    legend('Location','best')

    fprintf('\nSelected cloud Wf = %.3f pps\n',wf_values(k_demo));
    fprintf('  Original spectral radius: %.5f\n',rho_demo);
    fprintf('  Eigen-linearized (from ID) spectral radius: %.5f\n',max(abs(lambda_eig_lin)));
    fprintf('  Projected spectral radius: %.5f\n',max(abs(lambda_demo_proj)));
    fprintf('  Scaling alpha: %.5f\n',alpha);
end

%% Helpers
function out_dyn = get_out_dyn(data)
    if isstruct(data) && isfield(data,'out_Dyn')
        out_dyn = data.out_Dyn;
        return
    end
    error('Could not find out_Dyn in loaded data file.');
end

function [psi, u] = extract_cloud_signals(out_dyn, Ts, tstart)
    t_raw = signal_time(out_dyn.eng.S0.W);
    t = (tstart:Ts:t_raw(end)).';

    x1 = signal_data(out_dyn.eng.Shaft.N_Fan);
    x2 = signal_data(out_dyn.eng.Shaft.N_HPC);
    y1 = signal_data(out_dyn.eng.S21.Pt);
    y2 = signal_data(out_dyn.eng.S25.Pt);
    u_raw = signal_data(out_dyn.cntrl.Wfact);

    x1_i = interp1(t_raw,x1,t,'linear','extrap');
    x2_i = interp1(t_raw,x2,t,'linear','extrap');
    y1_i = interp1(t_raw,y1,t,'linear','extrap');
    y2_i = interp1(t_raw,y2,t,'linear','extrap');
    u_i  = interp1(t_raw,u_raw,t,'linear','extrap');

    psi = [x1_i x2_i y1_i y2_i];
    u = u_i;
end

function [A, B, rmse] = fit_local_linear_model(psi, u, ridge_lambda)
    % Fit psi_{k+1} = A psi_k + B u_k from a single data cloud.
    Xk = psi(1:end-1,:).';
    Xkp = psi(2:end,:).';
    Uk = u(1:end-1).';

    reg = [Xk; Uk];
    nreg = size(reg,1);

    M = reg*reg.' + ridge_lambda*eye(nreg);
    Theta = (Xkp*reg.') / M;

    nx = size(psi,2);
    A = Theta(:,1:nx);
    B = Theta(:,nx+1:end);

    pred = (A*Xk + B*Uk).';
    err = psi(2:end,:) - pred;
    rmse = sqrt(mean(err(:).^2));
end

function [A_proj, alpha] = project_to_stable_disk(A, rho_margin)
    rho = max(abs(eig(A)));
    if rho > rho_margin
        alpha = rho_margin/rho;
    else
        alpha = 1;
    end
    A_proj = alpha*A;
end

function ev2_aligned = align_poles(ev_ref, ev2)
    n = numel(ev_ref);
    ev2_aligned = zeros(size(ev2));
    used = false(n,1);
    for i = 1:n
        d = abs(ev2 - ev_ref(i));
        d(used) = inf;
        [~,idx] = min(d);
        ev2_aligned(i) = ev2(idx);
        used(idx) = true;
    end
end

function A_new = impose_eigenvalues(A_base, lambda_target)
    [V,D] = eig(A_base);
    lambda_base = diag(D);
    lambda_t_aligned = align_poles(lambda_base, lambda_target);

    if rcond(V) < 1e-10
        A_new = A_base;
        return
    end

    A_new = real(V*diag(lambda_t_aligned)/V);
end

function psi_hat = one_step_predict_id(A, B, psi_meas, u_meas, psi0, u0)
    % One-step predictor for deviation-form ID model:
    % psi_{k+1} = psi0 + A(psi_k-psi0) + B(u_k-u0)
    N = size(psi_meas,1);
    nx = size(psi_meas,2);
    psi_hat = nan(N,nx);
    psi_hat(1,:) = psi_meas(1,:);

    for k = 1:N-1
        xk = psi_meas(k,:).' - psi0(:);
        uk = u_meas(k,:) - u0(:).';
        xkp1 = psi0(:) + A*xk + B*uk.';
        psi_hat(k+1,:) = xkp1.';
    end
end

function psi_hat = one_step_predict_absolute(A, B, psi_meas, u_meas)
    % One-step predictor for absolute-form model:
    % psi_{k+1} = A psi_k + B u_k
    N = size(psi_meas,1);
    nx = size(psi_meas,2);
    psi_hat = nan(N,nx);
    psi_hat(1,:) = psi_meas(1,:);

    for k = 1:N-1
        xkp1 = A*psi_meas(k,:).' + B*u_meas(k,:).';
        psi_hat(k+1,:) = xkp1.';
    end
end

function xdev = simulate_deviation(A, B, u_dev)
    N = numel(u_dev);
    nx = size(A,1);
    x = zeros(nx,1);
    xdev = zeros(N,nx);
    for k = 1:N
        x = A*x + B*u_dev(k);
        xdev(k,:) = x.';
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
