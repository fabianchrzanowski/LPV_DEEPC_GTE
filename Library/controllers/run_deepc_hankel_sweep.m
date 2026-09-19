% ========================================================================
% run_deepc_hankel_sweep.m
%
% Automated sweep running the LPV-DeePC closed loop controller across
% a varying number of operating points (Hankel matrices).
% ========================================================================
clearvars; close all; clc;

proj = openProject('GTE.prj');
cd(proj.RootFolder);

run('setup_simulation_params.m');
scen_tag = lower(SCENARIO);

% Define the sweep of Hankel counts requested:
hankel_sweep = [32, 16, 8, 6, 5, 4, 3, 2 ,1];

% Set the number of columns
% If not set, it defaults to 3000
OVERRIDE_MAX_COLS = 500; 

fprintf('======================================================\n');
fprintf(' Starting LPV-DeePC Hankel Sweep (%d runs)\n', length(hankel_sweep));
fprintf('======================================================\n\n');

for i = 1:length(hankel_sweep)
    h_count = hankel_sweep(i);
    fprintf('\n\n>>>>> SWEEP RUN %d/%d : %d Hankel Matrices <<<<<\n', i, length(hankel_sweep), h_count);
    
    % Set override for the main script
    OVERRIDE_NUM_HANKELS = h_count;
    
    % Run the controller script
    try
        run('AGTF30_lpv_deepc_casadi_fast_closed_loop_nonlinear.m');
        fprintf('\n>>> Successfully completed run with %d Hankels.\n', h_count);
    catch ME
        fprintf('\n>>> ERROR during run with %d Hankels: %s\n', h_count, ME.message);
    end
end

fprintf('\n\n======================================================\n');
fprintf(' Sweep complete! %d files saved to Results/ directory.\n', length(hankel_sweep));
fprintf('======================================================\n');
