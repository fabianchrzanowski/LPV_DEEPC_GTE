% ========================================================================
% Controller Comparison
%
% This script scans the Results/ folder for any generated controller
% data files. It allows the user to select which controllers to plot
% via the command window, and dynamically generates the comparison plots.
% ========================================================================

clearvars; close all; clc
plot_options;

% ===== USER CONFIGURATION =====
% Toggle to plot vertical lines for database (Hankel) switches
SHOW_SWITCH_LINES = false; 
% ==============================

try proj = currentProject; catch; proj = openProject('GTE.prj'); end
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
selection = input('Enter the controllers (e.g., [1, 2, 4]), or Enter for all: ');

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

% Find the global time vector (the longest one, which includes warmup)
max_len = 0;
for i = 1:length(loaded_data)
    if length(loaded_data(i).data.time) > max_len
        time = loaded_data(i).data.time;
        ref_plot = loaded_data(i).data.ref_plot;
        track_name = loaded_data(i).data.track_name;
        max_len = length(time);
    end
end
N_total = max_len;
% In setup_simulation_params, T_ini = 20, steps_per_seg = 15 or 50
T_ini = 20;
Ts = 0.015;
N_sim = N_total - T_ini;

% OVERRIDE time vector so t=0 is the start of the active control phase
time = (-T_ini : N_sim-1)' * Ts;
% ref_plot is already loaded correctly

% Assuming standard scenario for segment boundaries (can be generalized)
steps_per_seg = (N_sim) / 3; % rough estimate for drawing vertical lines
seg_bounds = (1:2)*steps_per_seg*Ts;

% Calculate RMSEs (evaluating from step T_ini+1 to end)
idx_eval = (T_ini+1) : N_total;
rmses = zeros(length(selection), 1);

disp(' ');
disp('--- Root Mean Square Error (RMSE) ---');
for idx = 1:length(selection)
    y = loaded_data(idx).data.y_hist;
    u = loaded_data(idx).data.u_hist;
    
    % If an array is shorter than N_total, pad it at the END 
    % so its internal indices (like the jump at t=36) stay aligned with DeePC
    if length(y) < N_total
        y = [y; nan(N_total - length(y), 1)];
        u = [u; nan(N_total - length(u), 1)];
    end
    
    % Fill warmup NaNs for DeePC so its line visually starts at t=0
    if length(y) == N_total && sum(isnan(y(1:T_ini))) > 0
        y(1:T_ini) = y(T_ini+1);
        u(1:T_ini) = u(T_ini+1);
    end
    
    loaded_data(idx).data.y_hist = y; % Save back for plotting
    loaded_data(idx).data.u_hist = u;
    
    rmse = sqrt(mean((y(idx_eval) - ref_plot(idx_eval)).^2, 'omitnan'));
    rmses(idx) = rmse;
    fprintf('%-15s: %6.1f rpm\n', loaded_data(idx).name, rmse);
end

%% ===== Plot Dashboard =====
fprintf('\nGenerating dashboard plots...\n');

% Define color palette
colors = [
    0.85 0.33 0.10;   % Orange (Frozen MPC)
    0.00 0.45 0.74;   % Blue (LPV MPC)
    0.60 0.20 0.80;   % Purple (Frozen DeePC)
    0.20 0.70 0.30;   % Green (LPV DeePC)
    0.93 0.69 0.13;   % Yellow
    0.49 0.18 0.56;   % Violet
];

% --- Figure 1: Tracking Performance & Control Effort ---
fig1 = figure('Name','Controller Comparison: Tracking & Control Effort');
tiledlayout(2,1, 'TileSpacing','compact', 'Padding','compact');

ax1 = nexttile; hold(ax1, 'on');

plot(ax1, time, ref_plot, 'k--', 'LineWidth', 1.5, 'DisplayName', 'Reference');

for idx = 1:length(selection)
    c_idx = mod(idx-1, size(colors,1)) + 1;
    lw = 1.2;
    if contains(loaded_data(idx).name, 'LPV')
        lw = 1.8; % Make LPV lines thicker and clear
    end
    
    clean_name = strrep(loaded_data(idx).name, '_', '\_');
    % Use the globally aligned 'time' vector, not the internal one which might not be padded
    plot(ax1, time, loaded_data(idx).data.y_hist, '-', 'Color', colors(c_idx, :), ...
        'LineWidth', lw, 'DisplayName', sprintf('%s (RMSE: %.0f)', clean_name, rmses(idx)));
end

for ii=1:length(seg_bounds)
    xline(ax1, seg_bounds(ii), 'k:', 'HandleVisibility', 'off');
end

% Format track name for LaTeX
if contains(lower(track_name), 'n_h')
    ytarget_label = '$N_H$ [rpm]';
elseif contains(lower(track_name), 'n_l')
    ytarget_label = '$N_L$ [rpm]';
else
    ytarget_label = strrep(track_name, '_', '\_');
end

ylabel(ax1, ytarget_label, 'Interpreter', 'latex', 'FontSize', 12);
title(ax1, 'Tracking Performance Across Envelope', 'Interpreter', 'latex', 'FontSize', 13);
legend(ax1, 'Location', 'best', 'Interpreter', 'latex', 'FontSize', 10);
grid(ax1, 'on');
set(ax1, 'TickLabelInterpreter', 'latex', 'FontSize', 11);

ax2 = nexttile; hold(ax2, 'on');

for idx = 1:length(selection)
    c_idx = mod(idx-1, size(colors,1)) + 1;
    lw = 1.2;
    if contains(loaded_data(idx).name, 'LPV')
        lw = 1.8;
    end
    stairs(ax2, time, loaded_data(idx).data.u_hist, '-', 'Color', colors(c_idx, :), 'LineWidth', lw);
end

for ii=1:length(seg_bounds)
    xline(ax2, seg_bounds(ii), 'k:', 'HandleVisibility', 'off');
end
ylabel(ax2, '$W_f$ [pps]', 'Interpreter', 'latex', 'FontSize', 12);
xlabel(ax2, 'Time [s]', 'Interpreter', 'latex', 'FontSize', 12);
title(ax2, 'Control Effort (Fuel Flow)', 'Interpreter', 'latex', 'FontSize', 13);
grid(ax2, 'on');
set(ax2, 'TickLabelInterpreter', 'latex', 'FontSize', 11);

% --- Figure 2: Scheduling Variable ---
deepc_idx = find(cellfun(@(x) contains(lower(x), 'lpv') && contains(lower(x), 'deepc'), {loaded_data.name}), 1);
if ~isempty(deepc_idx) && isfield(loaded_data(deepc_idx).data, 'wf_scheduler_hist')
    fig2 = figure('Name','Controller Comparison: Scheduling Variable');
    ax3 = axes(fig2); hold(ax3, 'on');
    
    wf_sched = loaded_data(deepc_idx).data.wf_scheduler_hist;
    
    % Ensure t_sched is EXACTLY the same length as wf_sched by taking the
    % active portion of the globally defined 'time' vector.
    t_sched = time(end - length(wf_sched) + 1 : end);
    
    stairs(ax3, t_sched, wf_sched, 'LineWidth', 1.8, 'Color', colors(4, :)); % Green color for LPV DeePC
    
    for ii=1:length(seg_bounds)
        xline(ax3, seg_bounds(ii), 'k:', 'HandleVisibility', 'off');
    end
    ylabel(ax3, 'Active OP ($W_f$) [pps]', 'Interpreter', 'latex', 'FontSize', 12);
    xlabel(ax3, 'Time [s]', 'Interpreter', 'latex', 'FontSize', 12);
    title(ax3, 'Controller Comparison: Scheduling Variable Trace', 'Interpreter', 'latex', 'FontSize', 13);
    grid(ax3, 'on');
    set(ax3, 'TickLabelInterpreter', 'latex', 'FontSize', 11);
end

% --- Figure 3: Tracking Error Bar Chart ---
figure('Name', 'Controller Comparison: Error Summary');
b = bar(rmses);
xticks(1:length(selection));
clean_names = cell(length(loaded_data), 1);
for idx = 1:length(loaded_data)
    clean_names{idx} = strrep(loaded_data(idx).name, '_', '\_');
end
xticklabels(clean_names);
set(gca, 'TickLabelInterpreter', 'latex', 'FontSize', 11);
b.FaceColor = 'flat';
for idx = 1:length(selection)
    c_idx = mod(idx-1, size(colors,1)) + 1;
    b.CData(idx,:) = colors(c_idx, :);
end
ylabel('RMSE [rpm]', 'Interpreter', 'latex', 'FontSize', 12);
title('Tracking Error Comparison', 'Interpreter', 'latex', 'FontSize', 13);
grid on;
