% ========================================================================
% Visualise the difference between all 32 PRBS data clouds and Step Responses
% Plots the full nonlinear envelope of the AGTF30 engine.
% ========================================================================
clearvars; close all; clc

%% ===== USER CONFIGURATION =====
PLOT_STEP_DATA = false;  % Set to false to hide step response traces and clouds
NUM_CLOUDS_TO_PLOT = 32; % Number of clouds to plot (e.g. 8 for evenly spaced subset, 32 for all)
% ==============================

try; proj = currentProject; catch; proj = openProject('GTE.prj'); end
cd(proj.RootFolder)

fprintf('Loading data for all 32 operating points...\n');

%% Load models and Wf values
models = load(fullfile(proj.RootFolder,'Params','tables','matrices_DE_wf.mat'));
wf_data = load(fullfile(proj.RootFolder,'Params','test_files_fuel_only','test_identification','wf_values.mat'));
wf_values = wf_data.Wf_values(:)';
n_op = numel(wf_values);

% Preallocate arrays for massive plotting
all_u_prbs = []; all_NL_prbs = []; all_NH_prbs = []; all_c_prbs = [];
all_u_step = []; all_NL_step = []; all_NH_step = []; all_c_step = [];

% Indices for outputs
idx_NL = find(contains(models.sys_eng.OutputName, 'N_L'));
idx_NH = find(contains(models.sys_eng.OutputName, 'N_H'));
if isempty(idx_NL); idx_NL = 1; end
if isempty(idx_NH); idx_NH = 2; end
Ts = 0.015;

% Determine which clouds to plot (evenly spaced)
if isempty(NUM_CLOUDS_TO_PLOT) || NUM_CLOUDS_TO_PLOT <= 0
    NUM_CLOUDS_TO_PLOT = n_op;
end
plot_indices = unique(round(linspace(1, n_op, min(n_op, NUM_CLOUDS_TO_PLOT))));
n_plot = length(plot_indices);

% Generate distinct discrete colors spanning ONLY the plotted clouds
cloud_colors = turbo(n_plot);

% Pre-create Time History figures and axes for incremental plotting inside the loop
fig5 = figure('Name', 'Time History: Fan Speed (N_L)');
ax5 = axes(fig5);
hold(ax5, 'on');

fig6 = figure('Name', 'Time History: Core Speed (N_H)');
ax6 = axes(fig6);
hold(ax6, 'on');

% Pre-create SysID-style distinct color figures
fig7 = figure('Name', 'SysID-Style: Fan Speed vs Wf');
ax7 = axes(fig7); hold(ax7, 'on');

fig8 = figure('Name', 'SysID-Style: Core Speed vs Wf');
ax8 = axes(fig8); hold(ax8, 'on');

sysid_colors = lines(7); % Default MATLAB colormap cycle

for p_idx = 1:n_plot
    index = plot_indices(p_idx);
    wf_eq = wf_values(index);
    base_color = cloud_colors(p_idx, :);
    
    % Slightly darken every alternating cloud to create visible banding where they overlap
    if mod(p_idx, 2) == 0
        base_color = base_color * 0.75;
    end
    
    dark_color = base_color * 0.4;
    
    %% 1. Load PRBS Data from Simulink
    prbs_file = fullfile(proj.RootFolder,'Params','test_files_fuel_only', ...
        'test_identification','AGTF30', sprintf('test_id__AGTF30__Wf%04d__.mat', index));
    
    if ~isfile(prbs_file)
        warning('Could not find PRBS data for index %d. Skipping.', index);
        continue;
    end
    
    data = load(prbs_file);
    out_dyn = data.out_Dyn;
    
    t_raw = out_dyn.eng.S0.W.Time(:);
    u_prbs = out_dyn.cntrl.Wfact.Data(:);
    N_L_prbs = out_dyn.eng.Shaft.N_Fan.Data(:);
    N_H_prbs = out_dyn.eng.Shaft.N_HPC.Data(:);
    
    % The PRBS signal starts at t=19.5s. only the active part.
    idx_active = find(t_raw >= 19.5);
    idx_active = idx_active(1:5:end); % Take every 5th point
    
    all_u_prbs  = [all_u_prbs; u_prbs(idx_active)];
    all_NL_prbs = [all_NL_prbs; N_L_prbs(idx_active)];
    all_NH_prbs = [all_NH_prbs; N_H_prbs(idx_active)];
    
    % Assign the exact RGB color to every point in this cloud
    all_c_prbs  = [all_c_prbs; repmat(base_color, length(idx_active), 1)];
    
    % --- Add SysID-style distinct plot (scatter instead of lines) ---
    c_sysid = sysid_colors(mod(p_idx-1, 7)+1, :);
    scatter(ax7, u_prbs, N_L_prbs, 5, c_sysid, 'filled', 'MarkerFaceAlpha', 0.6);
    scatter(ax8, u_prbs, N_H_prbs, 5, c_sysid, 'filled', 'MarkerFaceAlpha', 0.6);
    % -------------------------------------------------------
    
    %% 2. Load Step Data from Simulink
    if PLOT_STEP_DATA
        step_file = fullfile(proj.RootFolder,'Params','test_files_fuel_only', ...
            'test_identification','AGTF30_Step_ALL', sprintf('test_id__AGTF30_Step__Wf%04d__.mat', index));
        
        if isfile(step_file)
            step_data = load(step_file);
            out_step = step_data.out_Dyn;
            
            t_step_raw = out_step.eng.S0.W.Time(:);
            
            if isfield(out_step.cntrl, 'Wfact')
                u_step = out_step.cntrl.Wfact.Data(:);
            else
                u_step = out_step.cntrl.Wf.Data(:);
            end
            N_L_step = out_step.eng.Shaft.N_Fan.Data(:);
            N_H_step = out_step.eng.Shaft.N_HPC.Data(:);
            
            % We only want the active part (e.g. t >= 19.5, same as PRBS)
            idx_step_active = find(t_step_raw >= 19.5);
            
            % Downsample step data for plotting
            all_u_step  = [all_u_step; u_step(idx_step_active(1:5:end))];
            all_NL_step = [all_NL_step; N_L_step(idx_step_active(1:5:end))];
            all_NH_step = [all_NH_step; N_H_step(idx_step_active(1:5:end))];
            
            % Assign the exact DARK RGB color to every point in this step path
            all_c_step  = [all_c_step; repmat(dark_color, length(idx_step_active(1:5:end)), 1)];

            t_step_ds = t_step_raw(1:5:end);
            NL_step_ds = N_L_step(1:5:end);
            NH_step_ds = N_H_step(1:5:end);

            % Step is semi-transparent but visible
            plot(ax5, t_step_ds, NL_step_ds, 'Color', [base_color, 0.65], 'LineWidth', 1.2)
            plot(ax6, t_step_ds, NH_step_ds, 'Color', [base_color, 0.65], 'LineWidth', 1.2)
        end
    end

    %% 2b. Plot Time-Domain Signals directly (full duration, absolute values)
    t_prbs_ds = t_raw(1:5:end);
    NL_prbs_ds = N_L_prbs(1:5:end);
    NH_prbs_ds = N_H_prbs(1:5:end);

    % PRBS is normal saturation
    plot(ax5, t_prbs_ds, NL_prbs_ds, 'Color', base_color, 'LineWidth', 1.3)
    plot(ax6, t_prbs_ds, NH_prbs_ds, 'Color', base_color, 'LineWidth', 1.3)
end

%% 3. Plot the Clouds in Separate Figures
fprintf('Plotting %d total data points...\n', numel(all_u_prbs));

% -- Figure 1: 3D Scatter (Wf vs NL vs NH) --
figure('Name', '3D Global State-Space Excitation');
scatter3(all_NL_prbs, all_NH_prbs, all_u_prbs, 5, all_c_prbs, 'filled', 'MarkerFaceAlpha', 0.1)
hold on
if PLOT_STEP_DATA && ~isempty(all_NL_step)
    scatter3(all_NL_step, all_NH_step, all_u_step, 8, all_c_step, 'filled', 'MarkerFaceAlpha', 1)
end
grid on
set(gca, 'TickLabelInterpreter', 'latex')
xlabel('Fan Speed ($N_L$) [rpm]', 'Interpreter', 'latex')
ylabel('Core Speed ($N_H$) [rpm]', 'Interpreter', 'latex')
zlabel('Fuel Flow ($W_f$) [pps]', 'Interpreter', 'latex')
title('\textbf{3D Global State-Space Excitation (All 32 Clouds)}', 'Interpreter', 'latex')
view(-35, 25)

% -- Figure 2: 2D State-Space (NL vs NH) --
figure('Name', '2D Projection: Fan vs Core Speed');
scatter(all_NL_prbs, all_NH_prbs, 5, all_c_prbs, 'filled', 'MarkerFaceAlpha', 0.1)
hold on
if PLOT_STEP_DATA && ~isempty(all_NL_step)
    scatter(all_NL_step, all_NH_step, 8, all_c_step, 'filled', 'MarkerFaceAlpha', 1)
end
grid on
set(gca, 'TickLabelInterpreter', 'latex')
xlabel('Fan Speed ($N_L$) [rpm]', 'Interpreter', 'latex')
ylabel('Core Speed ($N_H$) [rpm]', 'Interpreter', 'latex')
title('\textbf{2D Projection ($N_L$ vs $N_H$)}', 'Interpreter', 'latex')

% -- Figure 3: Operating Line (Wf vs NL) --
figure('Name', 'Operating Line: Fan Speed');
scatter(all_u_prbs, all_NL_prbs, 5, all_c_prbs, 'filled', 'MarkerFaceAlpha', 0.1)
hold on
if PLOT_STEP_DATA && ~isempty(all_NL_step)
    scatter(all_u_step, all_NL_step, 8, all_c_step, 'filled', 'MarkerFaceAlpha', 1)
end
grid on
set(gca, 'TickLabelInterpreter', 'latex')
xlabel('Fuel Flow ($W_f$) [pps]', 'Interpreter', 'latex')
ylabel('Fan Speed ($N_L$) [rpm]', 'Interpreter', 'latex')
title('\textbf{Operating Line: Fan Speed ($N_L$)}', 'Interpreter', 'latex')

% -- Figure 4: Operating Line (Wf vs NH) --
figure('Name', 'Operating Line: Core Speed');
scatter(all_u_prbs, all_NH_prbs, 5, all_c_prbs, 'filled', 'MarkerFaceAlpha', 0.1)
hold on
if PLOT_STEP_DATA && ~isempty(all_NL_step)
    scatter(all_u_step, all_NH_step, 8, all_c_step, 'filled', 'MarkerFaceAlpha', 1)
end
grid on
set(gca, 'TickLabelInterpreter', 'latex')
xlabel('Fuel Flow ($W_f$) [pps]', 'Interpreter', 'latex')
ylabel('Core Speed ($N_H$) [rpm]', 'Interpreter', 'latex')
title('\textbf{Operating Line: Core Speed ($N_H$)}', 'Interpreter', 'latex')

% -- Format Figure 5: Fan Speed vs Time --
figure(fig5);
grid(ax5, 'on');
set(ax5, 'TickLabelInterpreter', 'latex')
xlabel(ax5, 'Time [s]', 'Interpreter', 'latex')
ylabel(ax5, 'Fan Speed ($N_L$) [rpm]', 'Interpreter', 'latex')
title(ax5, '\textbf{Time History: Fan Speed ($N_L$) (PRBS vs Step)}', 'Interpreter', 'latex')

% -- Format Figure 6: Core Speed vs Time --
figure(fig6);
grid(ax6, 'on');
set(ax6, 'TickLabelInterpreter', 'latex')
xlabel(ax6, 'Time [s]', 'Interpreter', 'latex')
ylabel(ax6, 'Core Speed ($N_H$) [rpm]', 'Interpreter', 'latex')
title(ax6, '\textbf{Time History: Core Speed ($N_H$) (PRBS vs Step)}', 'Interpreter', 'latex')

% -- Format Figure 7: SysID-style Fan Speed --
figure(fig7);
grid on; set(gca, 'TickLabelInterpreter', 'latex')
xlabel('Fuel Flow ($W_f$) [pps]', 'Interpreter', 'latex')
ylabel('Fan Speed ($N_L$) [rpm]', 'Interpreter', 'latex')
title('\textbf{Fan Speed ($N_L$) vs Wf (SysID Colors)}', 'Interpreter', 'latex')

% -- Format Figure 8: SysID-style Core Speed --
figure(fig8);
grid on; set(gca, 'TickLabelInterpreter', 'latex')
xlabel('Fuel Flow ($W_f$) [pps]', 'Interpreter', 'latex')
ylabel('Core Speed ($N_H$) [rpm]', 'Interpreter', 'latex')
title('\textbf{Core Speed ($N_H$) vs Wf (SysID Colors)}', 'Interpreter', 'latex')

fprintf('Finished generating separate figures for full envelope cloud comparison!\n');
