% ========================================================================
% COMPARISON: element-wise vs eigenvalue interpolation for LPV models
%
% Produces dissertation-quality figures comparing:
%   1. Eigenvalue locus across rho (stability comparison)
%   2. Open-loop prediction accuracy vs Simulink plant
%   3. Closed-loop MPC tracking performance comparison
%
% Run AGTF30_sys_id_4.m first to generate the LUT.
%
% ========================================================================
% For AGTF30 Dissertation
% ========================================================================

clearvars; close all; clc

proj = currentProject;
cd(proj.RootFolder)
plot_options;

%% ===== LOAD DATA =====
lpv_file = fullfile(proj.RootFolder,'Params','tables','models',...
    'AGTF30_model_tables__1_input.mat');
LUT = load(lpv_file);

raw_file = fullfile(proj.RootFolder,'Params','tables','matrices_DE_wf.mat');
raw = load(raw_file);

wf_file = fullfile(proj.RootFolder,'Params','test_files_fuel_only',...
    'test_identification','wf_values.mat');
wf_data = load(wf_file);
Wf_values = wf_data.Wf_values(:)';

Ts = LUT.matrices.Ts;

%% ===== INTERPOLATION GRIDS =====
rho = raw.offset.u(1,:)';
d_rho = 0.001;
rho_es = rho(1):d_rho:rho(end);

%% ===== SUBSAMPLE: sparse operating points to stress-test interpolation =====
% With all ~32 points, element-wise interpolation stays stable.
% With only every 5th point, the gaps are wide enough that element-wise
% can produce unstable intermediate eigenvalues — this is the realistic
% scenario where you can't afford many identification experiments.
sparse_step = 5;  % use every 5th operating point (try 3, 5, 7)
sparse_idx  = 1:sparse_step:numel(rho);

% Ensure first and last point are included
if sparse_idx(end) ~= numel(rho)
    sparse_idx = [sparse_idx, numel(rho)];
end

rho_sparse = rho(sparse_idx);
A_sparse   = raw.A_dt(:,:,sparse_idx);

fprintf('Using %d of %d operating points (every %dth)\n', ...
    numel(sparse_idx), numel(rho), sparse_step);
fprintf('  rho grid: %s\n', mat2str(rho_sparse', 3));

%% ===== BUILD BOTH A MATRICES (from sparse grid) =====
fprintf('Building element-wise A matrix (sparse)... ');
A_elem = interp_matrix(rho_sparse, A_sparse, rho_es);
fprintf('done.\n');

fprintf('Building eigenvalue-interpolated A matrix (sparse, MAC)... ');
eig_opts.method  = 'pchip';
eig_opts.nom_idx = round(numel(rho_sparse)/2);
[A_eig, eig_info] = interp_matrix_eig(rho_sparse, A_sparse, rho_es, eig_opts);
fprintf('done.\n');

% B, C, D from sparse grid too (for consistency)
B = interp_matrix(rho_sparse, raw.B_dt(:,:,sparse_idx), rho_es);
C = interp_matrix(rho_sparse, raw.C_dt(:,:,sparse_idx), rho_es);
D = interp_matrix(rho_sparse, raw.D_dt(:,:,sparse_idx), rho_es);

% Offsets (still use all points for smooth offsets)
Xss = interp1(rho', raw.offset.x', rho_es, 'pchip', 'extrap')';
Yss = interp1(rho', raw.offset.y', rho_es, 'pchip', 'extrap')';
Uss = interp1(rho', raw.offset.u', rho_es, 'pchip', 'extrap')';
if size(Uss,1) >= size(Uss,2); Uss = Uss'; end
Zss = interp1(rho', raw.offset.z', rho_es, 'pchip', 'extrap')';

nx = size(raw.offset.x, 1);
ny = size(raw.offset.y, 1);

%% ===== FIGURE 1: EIGENVALUE LOCUS COMPARISON =====
fprintf('Generating eigenvalue comparison figure...\n');

n = size(A_elem, 1);
l = size(A_elem, 3);
ev_elem = nan(n, l);
ev_eig  = nan(n, l);
ev_elem_complex = nan(n, l);
ev_eig_complex  = nan(n, l);

for i = 1:l
    evals = sort(eig(A_elem(:,:,i)));
    ev_elem(:,i)         = abs(evals);
    ev_elem_complex(:,i) = evals;

    evals = sort(eig(A_eig(:,:,i)));
    ev_eig(:,i)         = abs(evals);
    ev_eig_complex(:,i) = evals;
end

% Count unstable points
n_unstable_elem = sum(ev_elem(:) > 1);
n_unstable_eig  = sum(ev_eig(:) > 1);

figure('Name','Fig 1: Eigenvalue magnitude comparison')
tiledlayout(2,2, 'TileSpacing','compact', 'Padding','compact');

% Top-left: element-wise |lambda| vs rho
nexttile;
plot(rho_es, ev_elem, '.', 'MarkerSize', 4); hold on;
yline(1, '--r', 'LineWidth', 1.5);
ylabel('$|\lambda|$');
xlabel('$\rho$ (W_f) [pps]');
grid on; set(gca,'FontSize',12);
title(sprintf('Element-wise interpolation (%d unstable pts)', n_unstable_elem));

% Top-right: eigenvalue-interp |lambda| vs rho
nexttile;
plot(rho_es, ev_eig, '.', 'MarkerSize', 4); hold on;
yline(1, '--r', 'LineWidth', 1.5);
ylabel('$|\lambda|$');
xlabel('$\rho$ (W_f) [pps]');
grid on; set(gca,'FontSize',12);
title(sprintf('Eigenvalue interpolation (%d unstable pts)', n_unstable_eig));

% Bottom-left: element-wise complex plane
nexttile;
theta_uc = linspace(0, 2*pi, 200);
plot(cos(theta_uc), sin(theta_uc), '--r', 'LineWidth', 1); hold on;
plot(real(ev_elem_complex(:)), imag(ev_elem_complex(:)), '.', 'MarkerSize', 3);
xlabel('$Re(\lambda)$'); ylabel('$Im(\lambda)$');
grid on; set(gca,'FontSize',12);
axis equal;
title('Complex plane (element-wise)');

% Bottom-right: eigenvalue-interp complex plane
nexttile;
plot(cos(theta_uc), sin(theta_uc), '--r', 'LineWidth', 1); hold on;
plot(real(ev_eig_complex(:)), imag(ev_eig_complex(:)), '.', 'MarkerSize', 3);
xlabel('$Re(\lambda)$'); ylabel('$Im(\lambda)$');
grid on; set(gca,'FontSize',12);
axis equal;
title('Complex plane (eigenvalue interp)');

fprintf('  Element-wise: max |lambda| = %.6f (%d unstable points)\n', ...
    max(ev_elem(:)), n_unstable_elem);
fprintf('  Eigenvalue:   max |lambda| = %.6f (%d unstable points)\n', ...
    max(ev_eig(:)), n_unstable_eig);

%% ===== FIGURE 2: OPEN-LOOP PREDICTION COMPARISON =====
% Simulate LPV model (both versions) against the Simulink data
% Uses the same validation trajectory as sys_id_5

fprintf('\nSetting up open-loop validation...\n');

model = 'AGTF30SysDyn';
MWS.engName = 'AGTF30';
MWS = setup_Controller(MWS);
MWS = setup_AllEng(MWS);
evalin('base','load(''Eng_Bus.mat'')');
MWS.In.Ts = Ts;

% Validation trajectory: aggressive profile spanning full operating range
% Segments: idle hold → ramp to max → hold at max → rapid steps → 
%           sinusoidal perturbation → takeoff-cruise-descent → return idle
dt_set = 0.15;
dt_trn = 0.5;

seg1 = 0.5*ones(1,20);                          % idle hold (3s)
seg2 = linspace(0.5, 1.9, 30);                   % ramp idle→max (4.5s)
seg3 = 1.9*ones(1,15);                           % hold at max (2.25s)
seg4 = [1.9 0.8 1.6 0.5 1.4 1.0 1.8 0.6 1.5];  % rapid steps
seg5 = 1.2 + 0.6*sin(2*pi*(0:39)/20);           % sinusoidal around 1.2
seg6 = [0.5*ones(1,10), linspace(0.5,1.7,15), ...% takeoff-cruise-descent
        1.7*ones(1,20), linspace(1.7,1.2,10), ...
        1.2*ones(1,15), linspace(1.2,0.5,10), ...
        0.5*ones(1,10)];

dx_Wf = [seg1, seg2, seg3, seg4, seg5, seg6];
dx_Wf = dx_Wf + 0.03*(2*rand(1,numel(dx_Wf)) - 1);  % small noise
Tsim   = dt_set * size(dx_Wf,2);
instants = 0:dt_set:Tsim;
t_model  = 0:Ts:Tsim;

MWS.In.Tsim = Tsim;
MWS.In.Alt  = timeseries([0;0],[0 Tsim],'Name','Alt');
MWS.In.MN   = timeseries([0;0],[0 Tsim],'Name','MN');
MWS.In.dT   = timeseries([0;0],[0 Tsim],'Name','dT');
MWS.In.Wf   = timeseries([dx_Wf(1);dx_Wf(1)],[0 Tsim],'Name','Wf');
MWS.In.dVBV = timeseries([0;0],[0 Tsim],'Name','dVBV');
MWS.In.dNz  = timeseries([0;0],[0 Tsim],'Name','dNz');
MWS = AGTF30_initial_conditions(MWS);

x0     = dx_Wf(1);
values = max(min(dx_Wf,20),0.2);
[x_sig, t_sig] = signal_transition_times(values, instants, dt_trn, x0);
MWS.In.Wf = timeseries(x_sig, t_sig);

fprintf('  Running Simulink plant...\n');
tic; out = sim(model); sim_time = toc;
fprintf('  Simulink done in %.1f s\n', sim_time);

% Engine tape
x_tape(:,1) = out_Dyn.eng.Shaft.N_Fan.Data;
x_tape(:,2) = out_Dyn.eng.Shaft.N_HPC.Data;
y_tape(:,1) = out_Dyn.eng.S21.Pt.Data;
y_tape(:,2) = out_Dyn.eng.S25.Pt.Data;
u_tape      = out_Dyn.cntrl.Wfact.Data;

% Run BOTH LPV models
methods = {'Element-wise', 'Eigenvalue'};
A_set = {A_elem, A_eig};
x_models = cell(1,2);
e_models = cell(1,2);

for mm = 1:2
    A_m = A_set{mm};

    x_state = [x_tape(1,1); x_tape(1,2)];
    y_state = [y_tape(1,1); y_tape(1,2)];
    x_model_m = zeros(numel(t_model), 2);
    y_model_m = zeros(numel(t_model), 2);

    for kk = 0:numel(t_model)-1
        rho_k = u_tape(kk+1);

        % Nearest-neighbour index into the fine grid
        [~, idx] = min(abs(rho_es - rho_k));

        A_k = A_m(:,:,idx);
        B_k = B(:,:,idx);

        x0_k = Xss(:, idx);
        y0_k = Yss(:, idx);
        u0_k = Uss(:, idx);

        x_model_m(kk+1,:) = x_state';
        y_model_m(kk+1,:) = y_state';

        state = [x0_k; y0_k] + A_k * ([x_state; y_state] - [x0_k; y0_k]) ...
                + B_k * (u_tape(kk+1) - u0_k);
        x_state = state(1:nx);
        y_state = state(nx+1:nx+ny);
    end

    x_models{mm} = x_model_m;
    e_models{mm} = abs(x_model_m - x_tape) ./ abs(x_tape) * 100;
end

% Plot
figure('Name','Fig 2: Open-loop prediction comparison')
tiledlayout(2,1,'TileSpacing','compact');

for j = 1:2
    nexttile;
    plot(t_model, x_tape(:,j), 'k', 'LineWidth', 1.2); hold on;
    plot(t_model, x_models{1}(:,j), '--', 'LineWidth', 1.0, 'Color', [0.85 0.33 0.1]);
    plot(t_model, x_models{2}(:,j), '-.', 'LineWidth', 1.0, 'Color', [0 0.45 0.74]);
    grid on; set(gca,'FontSize',12);
    if j == 1; ylabel('N_L [rpm]'); else; ylabel('N_H [rpm]'); end
    xlabel('Time [s]');
    legend('Simulink plant', 'Element-wise LPV', 'Eigenvalue LPV', ...
        'Location','best');
end

% Error statistics
figure('Name','Fig 3: Prediction error histograms')
tiledlayout(1,2,'TileSpacing','compact');

titles = {'N_L error [%]', 'N_H error [%]'};
for j = 1:2
    nexttile;
    histogram(e_models{1}(:,j), 50, 'Normalization','probability', ...
        'FaceAlpha',0.5, 'DisplayName','Element-wise'); hold on;
    histogram(e_models{2}(:,j), 50, 'Normalization','probability', ...
        'FaceAlpha',0.5, 'DisplayName','Eigenvalue');
    xlabel(titles{j}); ylabel('Probability');
    legend('Location','best');
    set(gca,'FontSize',12); grid on;
    title(sprintf('Mean err: elem=%.2f%%, eig=%.2f%%', ...
        mean(e_models{1}(:,j)), mean(e_models{2}(:,j))));
end

%% ===== SUMMARY TABLE =====
fprintf('\n========== COMPARISON SUMMARY ==========\n');
fprintf('%-25s  %-15s  %-15s\n', 'Metric', 'Element-wise', 'Eigenvalue');
fprintf('%-25s  %-15d  %-15d\n', 'Unstable |lambda|>1 pts', n_unstable_elem, n_unstable_eig);
fprintf('%-25s  %-15.6f  %-15.6f\n', 'Max |lambda|', max(ev_elem(:)), max(ev_eig(:)));
fprintf('%-25s  %-15.2f  %-15.2f\n', 'Mean N_L error [%]', ...
    mean(e_models{1}(:,1)), mean(e_models{2}(:,1)));
fprintf('%-25s  %-15.2f  %-15.2f\n', 'Mean N_H error [%]', ...
    mean(e_models{1}(:,2)), mean(e_models{2}(:,2)));
fprintf('=========================================\n');
