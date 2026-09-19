% ========================================================================
% LPV-DeePC Solve Times Comparison
%
% Generates a plot of the CasADi solver execution times per control step.
% Calculates and prints the average and peak solve times, comparing them
% against the 15ms Simulink.
% ========================================================================

clearvars; close all; clc

try; proj = currentProject; catch; proj = openProject('GTE.prj'); end
cd(proj.RootFolder)

plot_options; % load latex formatting

% ========================================================================
% ===== USER CONFIGURATION ===============================================
% ========================================================================
% Leave empty '' to auto-scan the Results/ folder
FILE_1 = '';
FILE_2 = '';
FILE_3 = ''; 
% ========================================================================

FADEC_LIMIT_MS = 15; % 0.015s sampling time

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

figure('Name', 'Solve Time Comparison', 'Position', [100, 100, 700, 400]);
hold on; grid on; box on;

yline(FADEC_LIMIT_MS, 'r--', 'LineWidth', 2, 'DisplayName', 'FADEC Hard Real-Time Limit (15 ms)');

valid_legs = {'FADEC Limit (15 ms)'};

fprintf('\n==========================================================\n');
fprintf('                SOLVE TIME ANALYSIS SUMMARY\n');
fprintf('==========================================================\n');

for i = 1:length(files)
    if isempty(files{i}); continue; end
    if ~isfile(files{i})
        warning('File not found: %s', files{i});
        continue;
    end
    
    data = load(files{i});
    
    if ~isfield(data, 'solve_times')
        warning('The file %s does not contain "solve_times". Skipping.', files{i});
        continue;
    end
    
    % Convert to ms
    st_ms = data.solve_times * 1000;
    
    % The first step in CasADi often involves internal memory allocation / JIT overhead,
    % so we typically exclude it from the "steady-state" average/peak statistics.
    if length(st_ms) > 1
        st_steady = st_ms(2:end);
    else
        st_steady = st_ms;
    end
    
    avg_st = mean(st_steady, 'omitnan');
    peak_st = max(st_steady, [], 'omitnan');
    
    first_valid_idx = find(~isnan(st_ms), 1);
    if isempty(first_valid_idx)
        first_st = NaN;
    else
        first_st = st_ms(first_valid_idx);
    end
    
    pct_viol = sum(st_steady > FADEC_LIMIT_MS) / sum(~isnan(st_steady)) * 100;
    
    [~, name, ~] = fileparts(files{i});
    name = strrep(name, 'results_lpv_deepc_fast_', '');
    name_disp = strrep(name, 'results_lpv_mpc_', 'MPC_');
    name_disp = strrep(name_disp, '_', '\_');
    
    fprintf('\n[%s]\n', name_disp);
    fprintf('  Average Solve Time: %6.2f ms\n', avg_st);
    fprintf('  Peak Solve Time:    %6.2f ms (excluding Step 1 init)\n', peak_st);
    fprintf('  First Step Init:    %6.2f ms\n', first_st);
    if pct_viol > 0
        fprintf('  [WARNING] Exceeded FADEC limit %.1f%% of the time.\n', pct_viol);
    else
        fprintf('  [OK] 100%% hard real-time feasible.\n');
    end
    
    % Plot (Scatter dots are better for seeing density)
    scatter(1:length(st_ms), st_ms, 15, colors(i,:), 'filled', 'MarkerFaceAlpha', 0.6);
    
    valid_legs{end+1} = name_disp; %#ok<AGROW>
end

fprintf('\n==========================================================\n');

if length(valid_legs) == 1
    fprintf('\nNo valid files with solve_times were selected. Exiting.\n');
    close(gcf);
    return;
end

xlabel('Control Step');
ylabel('Solve Time [ms]');
title('Real-Time Feasibility: CasADi Optimization Time');
legend(valid_legs, 'Location', 'best', 'Interpreter', 'latex', 'FontSize', 9);
ylim([0, max(30, max(yticks))]); % Ensure minimum y-axis height to see 15ms clearly

fprintf('\nPlot generated successfully!\n');
