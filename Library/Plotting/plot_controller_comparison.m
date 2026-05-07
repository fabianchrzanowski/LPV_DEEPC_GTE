% ========================================================================
% Controller Comparison Dashboard (Modular)
%
% This script scans the Results/ folder for any generated controller
% data files. It allows the user to select which controllers to plot
% via the command window, and dynamically generates the comparison plots.
% ========================================================================
% Fabian Chrzanowski — University of Sheffield — 2026
% ========================================================================

clearvars; close all; clc
plot_options;

proj = currentProject;
cd(proj.RootFolder)

results_dir = fullfile(proj.RootFolder, 'Results');
if ~isfolder(results_dir)
    error('Results folder not found. Please run the controller scripts first.');
end

% Scan for available .mat files
files = dir(fullfile(results_dir, 'results_*.mat'));
if isempty(files)
    error('No results found in %s. Please run one or more controller scripts first.', results_dir);
end

fprintf('\n=======================================================\n');
fprintf('   CONTROLLER COMPARISON DASHBOARD\n');
fprintf('=======================================================\n');
fprintf('Found the following controller results:\n');

% Load headers to get controller names
available_results = cell(length(files), 1);
for i = 1:length(files)
    data = load(fullfile(results_dir, files(i).name));
    if isfield(data, 'controller_name')
        available_results{i} = data.controller_name;
    else
        % Fallback: derive name from filename (e.g. results_mpc.mat -> mpc)
        [~, fname] = fileparts(files(i).name);
        available_results{i} = strrep(strrep(fname, 'results_', ''), '_', ' ');
    end
    fprintf('  [%d] %s\n', i, available_results{i});
end

% Ask user which ones to plot
disp(' ');
selection = input('Enter the numbers of the controllers to plot as an array (e.g., [1, 2, 4]), or press Enter for all: ');

if isempty(selection)
    selection = 1:length(files);
end

% Validate selection
selection = selection(selection >= 1 & selection <= length(files));
if isempty(selection)
    error('No valid controllers selected.');
end

%% ===== Load Selected Data =====
% We will store data in struct arrays
loaded_data = struct();
for idx = 1:length(selection)
    i = selection(idx);
    file_path = fullfile(results_dir, files(i).name);
    loaded_data(idx).data = load(file_path);
    loaded_data(idx).name = available_results{i};
end

% Assume all controllers were run with the same setup_simulation_params
% so they share the same time vector and reference plot.
time = loaded_data(1).data.time;
ref_plot = loaded_data(1).data.ref_plot;
track_name = loaded_data(1).data.track_name;

N_total = length(time);
% In setup_simulation_params, T_ini = 10, steps_per_seg = 200
T_ini = 10;
steps_per_seg = 200;
Ts = 0.015;
seg_bounds = T_ini*Ts + (1:2)*steps_per_seg*Ts;

% Calculate RMSEs (evaluating from step T_ini+1 to end)
idx_eval = (T_ini+1) : N_total;
rmses = zeros(length(selection), 1);

disp(' ');
disp('--- Root Mean Square Error (RMSE) ---');
for idx = 1:length(selection)
    y = loaded_data(idx).data.y_hist;
    % Align arrays if one was missing warmup padding (e.g. MPC vs DeePC)
    if length(y) < N_total
        % Pad start with NaNs
        y = [nan(N_total - length(y), 1); y];
        loaded_data(idx).data.y_hist = y; % Save back
        
        u = loaded_data(idx).data.u_hist;
        u = [nan(N_total - length(u), 1); u];
        loaded_data(idx).data.u_hist = u;
    end
    
    rmse = sqrt(mean((y(idx_eval) - ref_plot(idx_eval)).^2, 'omitnan'));
    rmses(idx) = rmse;
    fprintf('%-15s: %6.1f rpm\n', loaded_data(idx).name, rmse);
end

%% ===== Plot Dashboard =====
fprintf('\nGenerating dashboard plots...\n');

% Define a robust color palette
colors = [
    0.85 0.33 0.10;   % Orange (Frozen MPC)
    0.00 0.45 0.74;   % Blue (LPV MPC)
    0.60 0.20 0.80;   % Purple (Frozen DeePC)
    0.20 0.70 0.30;   % Green (LPV DeePC)
    0.93 0.69 0.13;   % Yellow
    0.49 0.18 0.56;   % Violet
];

% --- Figure 1: Tracking Performance ---
figure('Name','Controller Comparison: Tracking', 'Units', 'pixels', 'Position', [100 100 900 600]);
tiledlayout(2,1, 'TileSpacing','compact', 'Padding','compact');

ax1 = nexttile; hold(ax1, 'on');
plot(ax1, time, ref_plot, 'k--', 'LineWidth', 1.5, 'DisplayName', 'Reference');

for idx = 1:length(selection)
    c_idx = mod(idx-1, size(colors,1)) + 1;
    lw = 1.2;
    if contains(loaded_data(idx).name, 'LPV')
        lw = 1.5; % Make LPV lines slightly thicker
    end
    plot(ax1, loaded_data(idx).data.time, loaded_data(idx).data.y_hist, '-', 'Color', colors(c_idx, :), ...
        'LineWidth', lw, 'DisplayName', sprintf('%s (RMSE: %.0f)', loaded_data(idx).name, rmses(idx)));
end

for ii=1:length(seg_bounds)
    xline(ax1, seg_bounds(ii), 'k:', 'HandleVisibility', 'off');
end
ylabel(ax1, sprintf('%s [rpm]', track_name), 'Interpreter', 'none');
title(ax1, 'Controller Comparison: Tracking Performance Across Envelope', 'Interpreter', 'none');
legend(ax1, 'Location', 'best', 'Interpreter', 'none');
grid(ax1, 'on');

ax2 = nexttile; hold(ax2, 'on');
for idx = 1:length(selection)
    c_idx = mod(idx-1, size(colors,1)) + 1;
    lw = 1.2;
    if contains(loaded_data(idx).name, 'LPV')
        lw = 1.5;
    end
    plot(ax2, loaded_data(idx).data.time, loaded_data(idx).data.u_hist, '-', 'Color', colors(c_idx, :), 'LineWidth', lw);
end

for ii=1:length(seg_bounds)
    xline(ax2, seg_bounds(ii), 'k:', 'HandleVisibility', 'off');
end
ylabel(ax2, 'W_f [pps]', 'Interpreter', 'none');
xlabel(ax2, 'Time [s]');
title(ax2, 'Control Effort (Fuel Flow)', 'Interpreter', 'none');
grid(ax2, 'on');

% --- Figure 2: Tracking Error Bar Chart ---
figure('Name', 'Controller Comparison: Error Summary', 'Units', 'pixels', 'Position', [1050 100 500 300]);
b = bar(categorical({loaded_data.name}), rmses);
b.FaceColor = 'flat';
for idx = 1:length(selection)
    c_idx = mod(idx-1, size(colors,1)) + 1;
    b.CData(idx,:) = colors(c_idx, :);
end
ylabel('RMSE [rpm]');
title('Tracking Error Comparison', 'Interpreter', 'none');
grid on;
