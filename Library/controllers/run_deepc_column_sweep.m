% ========================================================================
% run_deepc_column_sweep.m
%
% Automated sweep running the LPV-DeePC closed loop controller across
% a variety of Hankel matrix column widths.
% ========================================================================
clearvars; close all; clc;

try proj = currentProject; catch; proj = openProject('GTE.prj'); end
cd(proj.RootFolder);

run('setup_simulation_params.m');
scen_tag = lower(SCENARIO);

% Define the sweep of columns requested:
column_sweep = [1500, 1250, 1000, 750, 500, 250, 125, 100, 90, 80, 70, 60, 56, 55, 54, 53, 45];

fprintf('======================================================\n');
fprintf(' Starting LPV-DeePC Column Sweep (%d runs)\n', length(column_sweep));
fprintf('======================================================\n\n');

for i = 1:length(column_sweep)
    cols = column_sweep(i);
    fprintf('\n\n>>>>> SWEEP RUN %d/%d : %d columns <<<<<\n', i, length(column_sweep), cols);
    
    % Set overrides for the main script
    OVERRIDE_MAX_COLS = cols;
    
    % Run the controller script
    try
        run('AGTF30_lpv_deepc_casadi_fast_closed_loop_nonlinear.m');
        fprintf('\n>>> Successfully completed run with %d columns.\n', cols);
    catch ME
        fprintf('\n>>> ERROR during run with %d columns: %s\n', cols, ME.message);
    end
end

fprintf('\n\n======================================================\n');
fprintf(' Sweep complete! %d files saved to Results/ directory.\n', length(column_sweep));
fprintf('======================================================\n');

