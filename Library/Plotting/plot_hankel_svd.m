% ========================================================================
% Plot Hankel Singular Value Spectrum
%
% This script loads the raw training/identification datasets directly,
% constructs the Hankel matrices for all 32 operating points, computes
% their SVDs, and plots the singular value spectrum color-coded by fuel flow.
% ========================================================================

clearvars; close all; clc;

% ===== USER CONFIGURATION =====
DATA_SOURCE   = 'SINE';      % Options: 'PRBS' (loads from AGTF30) or 'STEP' (loads from AGTF30_Step_ALL)
SHOW_MIDPOINT = false;       % Set true to overlay and highlight the midpoint OP with black markers
COLOR_MODE    = 'turbo'; % Options: 'duo-color' (blue-to-red gradient) or 'turbo' (rainbow colormap)
% ==============================

try proj = currentProject; catch; proj = openProject('GTE.prj'); end
cd(proj.RootFolder);

% Load formatting options
if isfile(fullfile(proj.RootFolder, 'Library', 'Plotting', 'plot_options.m'))
    run('plot_options.m');
end

% Set up paths and load basic simulation parameters (for Ts, T_ini, N_pred, L)
run('setup_simulation_params.m');

addpath(fullfile(proj.RootFolder, 'Library', 'controllers'));

models = load(fullfile(proj.RootFolder,'Params','tables','matrices_DE_wf.mat'));
wf_data = load(fullfile(proj.RootFolder,'Params','test_files_fuel_only','test_identification','wf_values.mat'));
sys_eng = models.sys_eng;
wf_values = wf_data.Wf_values(:)';

idx_track = find(contains(lower(string(sys_eng.OutputName)), "n_h"), 1);
if isempty(idx_track); idx_track = 2; end
track_name = sys_eng.OutputName{idx_track};

u_min = 0.2; u_max = 2.2; % fuel constraints
n_sys = 4;
rank_limit = L + n_sys;

is_synthetic = ismember(DATA_SOURCE, {'SINE', 'WHITE_NOISE', 'STEADY_STATE'});

if is_synthetic
    N_gen = 2000;
    amp_fraction = 0.05; % 5% perturbation
    n_op = numel(wf_values);
    
    fprintf('Generating synthetic %s data and computing SVD of Hankel matrices...\n', DATA_SOURCE);
    tic;
    
    hankel_U = cell(n_op, 1);
    hankel_Y = cell(n_op, 1);
    hankel_wf = wf_values;
    
    for k = 1:n_op
        wf_eq = wf_values(k);
        amp = amp_fraction * wf_eq;
        
        rng(42 + k);
        if strcmp(DATA_SOURCE, 'SINE')
            u_train = wf_eq + amp*sin(2*pi*(1:N_gen)'/50);
        elseif strcmp(DATA_SOURCE, 'WHITE_NOISE')
            u_train = wf_eq + amp*randn(N_gen,1);
        elseif strcmp(DATA_SOURCE, 'STEADY_STATE')
            u_train = wf_eq*ones(N_gen,1);
        end
        u_train = max(u_min, min(u_max, u_train)); % clamp
        
        % Simulate using fast LPV plant
        x_state = zeros(size(models.A_dt, 1), 1);
        y_train = zeros(N_gen, 1);
        
        for i = 1:N_gen
            [y_train(i), x_state] = LPV_step_fast(x_state, u_train(i), wf_eq, models);
        end
        
        u_mean = mean(u_train);
        y_mean = mean(y_train);
        du = u_train - u_mean;
        dy = y_train - y_mean;
        
        hankel_U{k} = block_hankel(du(:), L);
        hankel_Y{k} = block_hankel(dy(:), L);
    end
else
    % File-based loading
    if strcmp(DATA_SOURCE, 'STEP')
        folder_name = 'AGTF30_Step_ALL';
    else
        folder_name = 'AGTF30';
    end
    
    training_folder = fullfile(proj.RootFolder,'Params','test_files_fuel_only','test_identification', folder_name);
    filelist = dir(fullfile(training_folder, '*.mat'));
    [~, idx_sort] = sort({filelist.name});
    filelist = filelist(idx_sort);
    n_op = min(numel(filelist), numel(wf_values));
    
    fprintf('Loading %d %s datasets and computing SVD of Hankel matrices...\n', n_op, DATA_SOURCE);
    tic;
    
    hankel_U = cell(n_op, 1);
    hankel_Y = cell(n_op, 1);
    hankel_wf = wf_values(1:n_op);
    
    for k = 1:n_op
        data = load(fullfile(filelist(k).folder, filelist(k).name));
        out_dyn = fetch_out_dyn(data);
        u = extract_input_trace(out_dyn);
        y = measure_tracking_series(out_dyn, track_name);
        n_common = min(numel(u), numel(y));
        
        u_mean = mean(u(1:n_common));
        y_mean = mean(y(1:n_common));
        du = u(1:n_common) - u_mean;
        dy = y(1:n_common) - y_mean;
        
        hankel_U{k} = block_hankel(du(:), L);
        hankel_Y{k} = block_hankel(dy(:), L);
    end
end

% Remove empty entries
valid = ~cellfun(@isempty, hankel_U);
hankel_U = hankel_U(valid);
hankel_Y = hankel_Y(valid);
hankel_wf = hankel_wf(valid);
n_valid = sum(valid);

% Create figure
fig = figure('Name', sprintf('Global Hankel SVD Spectrum (%s)', DATA_SOURCE));
ax = axes(fig);
hold(ax, 'on');
set(ax, 'YScale', 'log');

% Colormap setup based on Wf range
wf_min = min(hankel_wf);
wf_max = max(hankel_wf);

if strcmp(COLOR_MODE, 'duo-color')
    % Custom gradient from Deep Blue (idle) to Vibrant Red (max thrust)
    c_start = [0.00 0.45 0.74]; % Deep Blue
    c_end   = [0.85 0.10 0.10]; % Vibrant Red
    colors = [linspace(c_start(1), c_end(1), n_valid)', ...
              linspace(c_start(2), c_end(2), n_valid)', ...
              linspace(c_start(3), c_end(3), n_valid)'];
    colormap(ax, colors);
else
    colors = turbo(n_valid);
    colormap(ax, 'turbo');
end

% Sort operating points by Wf to ensure smooth color transition
[sorted_wf, sort_idx] = sort(hankel_wf);

for idx = 1:n_valid
    k = sort_idx(idx);
    Hu_k = hankel_U{k};
    Hy_k = hankel_Y{k};
    H_k = [Hu_k; Hy_k];
    
    s_vals = svd(H_k);
    
    % Color from colormap based on sorted index
    c = colors(idx, :);
    
    plot(ax, s_vals, '-', 'Color', c, 'LineWidth', 1.2, 'HandleVisibility', 'off');
end

% Overlay the midpoint 
if SHOW_MIDPOINT
    mid_idx = sort_idx(round(n_valid / 2));
    H_mid = [hankel_U{mid_idx}; hankel_Y{mid_idx}];
    s_vals_mid = svd(H_mid);
    semilogy(ax, s_vals_mid, 'o--', 'Color', 'k', 'LineWidth', 1.8, ...
        'MarkerFaceColor', 'k', 'MarkerSize', 4, ...
        'DisplayName', sprintf('Midpoint OP ($W_f = %.2f$ pps)', hankel_wf(mid_idx)));
end

% Plot numerical rank cutoff 
xline(ax, rank_limit, 'r--', 'LineWidth', 2.0, ...
    'DisplayName', sprintf('Willems Rank Limit ($L+n_x = %d$)', rank_limit));

% Add colorbar to explain color range
cb = colorbar(ax);
if exist('clim', 'file')
    clim(ax, [wf_min, wf_max]);
else
    caxis(ax, [wf_min, wf_max]);
end
ylabel(cb, 'Operating Point $W_f$ [pps]', 'Interpreter', 'latex', 'FontSize', 11);

set(ax, 'TickLabelInterpreter', 'latex', 'FontSize', 11);
xlabel(ax, 'Singular Value Index', 'Interpreter', 'latex', 'FontSize', 12);
ylabel(ax, 'Singular Value $\sigma_i$ (Log Scale)', 'Interpreter', 'latex', 'FontSize', 12);
title(ax, sprintf('\\textbf{Hankel Singular Value Spectrum Across Operating Envelope (%s Data)}', DATA_SOURCE), 'Interpreter', 'latex', 'FontSize', 13);
legend(ax, 'Location', 'northeast', 'Interpreter', 'latex');
grid(ax, 'on');

elapsed_time = toc;
fprintf('SVD spectrum plot generated successfully in %.3f seconds.\n', elapsed_time);

%% ===== Local Helper Functions =====
function H = block_hankel(sig, L)
    sig = sig(:); nc = numel(sig)-L+1;
    H = zeros(L, nc);
    for i = 1:L; H(i,:) = sig(i:i+nc-1).'; end
end

function out_dyn = fetch_out_dyn(data)
    if isstruct(data) && isfield(data,'out_Dyn'); out_dyn = data.out_Dyn; return; end
    if isa(data,'Simulink.SimulationOutput')
        try; out_dyn = data.get('out_Dyn'); return; catch; end
    end
    if isstruct(data) && isfield(data,'eng'); out_dyn = data; return; end
    error('Could not find out_Dyn');
end

function u = extract_input_trace(out_dyn)
    if isfield(out_dyn,'cntrl')
        if isfield(out_dyn.cntrl,'Wfact'); u = sv(out_dyn.cntrl.Wfact); return; end
        if isfield(out_dyn.cntrl,'Wf'); u = sv(out_dyn.cntrl.Wf); return; end
    end
    error('No Wf trace found');
end

function y = measure_tracking_series(out_dyn, name)
    nl = lower(string(name));
    if contains(nl,'n_h'); y = sv(out_dyn.eng.Shaft.N_HPC);
    elseif contains(nl,'n_l'); y = sv(out_dyn.eng.Shaft.N_Fan);
    elseif contains(nl,'thrust')||contains(nl,'fnet'); y = sv(out_dyn.eng.Perf.Fnet);
    else; y = sv(out_dyn.eng.Shaft.N_HPC); end
end

function x = sv(s)
    if isstruct(s)&&isfield(s,'Data'); x = double(s.Data(:));
    elseif isa(s,'timeseries'); x = double(s.Data(:));
    else; x = double(s(:)); end
end
