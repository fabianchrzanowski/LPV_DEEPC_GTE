% ========================================================================
% LPV-DeePC Dissertation Plotter: Objective Cost Over Time
%
% Generates a plot of the CasADi objective function cost (sol.f) over time.
% Allows comparing multiple controller runs.
%
% ========================================================================

clearvars; close all; clc

try; proj = currentProject; catch; proj = openProject('GTE.prj'); end
cd(proj.RootFolder)

plot_options; % load latex formatting

% ========================================================================
% Leave empty '' to auto-scan the Results/ folder
FILE_1 = '';
FILE_2 = '';
FILE_3 = ''; 
% ========================================================================

results_dir = fullfile(proj.RootFolder, 'Results');

% Helper function to pick files
function fpath = pick_file(results_dir, prompt_text)
    hits = dir(fullfile(results_dir, '*.mat'));
    if isempty(hits)
        fpath = '';
        return
    end
    fprintf('\n--- %s ---\n', prompt_text);
    fprintf('  0) [Skip this file]\n');
    for i = 1:length(hits)
        fprintf('%3d) %s\n', i, hits(i).name);
    end
    
    sel = input('Select file number (or 0 to skip): ');
    if isempty(sel) || sel < 1 || sel > length(hits)
        fpath = '';
    else
        fpath = fullfile(hits(sel).folder, hits(sel).name);
    end
end

% Prompt if not set
if isempty(FILE_1); FILE_1 = pick_file(results_dir, 'Select First File to Plot'); end
if isempty(FILE_2); FILE_2 = pick_file(results_dir, 'Select Second File to Compare (optional)'); end
if isempty(FILE_3); FILE_3 = pick_file(results_dir, 'Select Third File to Compare (optional)'); end

files = {FILE_1, FILE_2, FILE_3};
colors = lines(3);

figure('Name', 'Objective Cost Comparison', 'Position', [100, 100, 700, 400]);
hold on; grid on; box on;

valid_legs = {};
for i = 1:length(files)
    if isempty(files{i}); continue; end
    if ~isfile(files{i})
        warning('File not found: %s', files{i});
        continue;
    end
    
    data = load(files{i});
    
    if ~isfield(data, 'cost_hist')
        warning('The file %s does not contain "cost_hist". This is an older run from before the cost tracking was added. Skipping.', files{i});
        continue;
    end
    
    cost_plot = data.cost_hist;
    
    % If cost_plot has N_total elements (like MPC), extract the active control region
    if length(cost_plot) == length(data.time)
        T_ini = 20; % Default T_ini
        cost_plot = cost_plot(T_ini+1:end);
    end
    
    if contains(lower(files{i}), 'mpc')
        % Shift MPC forward by 1 step to match DeePC alignment
        cost_plot = [cost_plot(2:end); cost_plot(end)];
    end
    
    % Plot against relative control time
    plot(data.time(1:length(cost_plot)), cost_plot, 'LineWidth', 1.5, 'Color', colors(i,:));
    
    % Generate label
    [~, name, ~] = fileparts(files{i});
    name = strrep(name, 'results_lpv_deepc_fast_', '');
    name = strrep(name, 'results_lpv_mpc_', 'MPC_');
    name = strrep(name, '_', '\_'); % Escape underscores for latex
    valid_legs{end+1} = name; 
end

if isempty(valid_legs)
    fprintf('\nNo valid files with cost_hist were selected. Exiting.\n');
    close(gcf);
    return;
end

xlabel('Time [s]');
ylabel('CasADi Objective Cost ($J$)');
title('Controller Objective Cost Over Time');
legend(valid_legs, 'Location', 'best', 'Interpreter', 'latex', 'FontSize', 9);

% Set y-axis to log scale (costs usually drop exponentially and have large spikes)
set(gca, 'YScale', 'log');

fprintf('\nPlot generated successfully!\n');
