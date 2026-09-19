% ========================================================================
% LPV-DeePC Constraint Violation Analysis
%
% Checks if the controller respected the physical input constraints and 
% input rate limits throughout the simulation.
% ========================================================================

clearvars; close all; clc

try proj = currentProject; catch; proj = openProject('GTE.prj'); end
cd(proj.RootFolder)

plot_options;

% Define constraints (matches controller settings)
U_MIN = 0.35;
U_MAX = 2.0;
DU_MAX = 0.15;

results_dir = fullfile(proj.RootFolder, 'Results');
hits = dir(fullfile(results_dir, '*.mat'));
if isempty(hits)
    error('No .mat files found in Results folder.');
end

fprintf('\n--- Select a File to Check Constraints ---\n');
for i = 1:length(hits)
    fprintf('%3d) %s\n', i, hits(i).name);
end
sel = input('Select file number: ');
if isempty(sel) || sel < 1 || sel > length(hits); return; end

fpath = fullfile(hits(sel).folder, hits(sel).name);
data = load(fpath);

u_hist = data.u_hist;
du_hist = [0; diff(u_hist)]; % Rate of change

t = data.time(1:length(u_hist));

%% Plot 1: Absolute Input Limits
figure('Name', 'Constraint Analysis: Input Magnitude', 'Position', [100, 100, 600, 300]);
hold on; grid on; box on;

% Draw constraint bounds
yregion(U_MIN, U_MAX, 'FaceColor', [0.8 0.9 0.8], 'FaceAlpha', 0.5, 'DisplayName', 'Feasible Region');
yline(U_MAX, 'r--', 'LineWidth', 1.5, 'DisplayName', 'U_{MAX}');
yline(U_MIN, 'r--', 'LineWidth', 1.5, 'DisplayName', 'U_{MIN}');

% Plot actual input
plot(t, u_hist, 'b-', 'LineWidth', 1.5, 'DisplayName', 'W_f Commanded');

xlabel('Time [s]');
ylabel('Fuel Flow W_f [pps]');
title(sprintf('Input Magnitude Constraints (U_{min}=%.2f, U_{max}=%.2f)', U_MIN, U_MAX));
legend('Location', 'best');

% Check violations
u_viol_high = sum(u_hist > U_MAX + 1e-5);
u_viol_low  = sum(u_hist < U_MIN - 1e-5);
fprintf('\n--- Input Magnitude (W_f) ---\n');
if u_viol_high > 0; fprintf('  [WARNING] Exceeded U_MAX %d times!\n', u_viol_high); else; fprintf('  [OK] Never exceeded U_MAX.\n'); end
if u_viol_low > 0;  fprintf('  [WARNING] Dropped below U_MIN %d times!\n', u_viol_low); else; fprintf('  [OK] Never dropped below U_MIN.\n'); end

%% Plot 2: Input Rate Limits
figure('Name', 'Constraint Analysis: Input Rate', 'Position', [150, 150, 600, 300]);
hold on; grid on; box on;

% Draw constraint bounds
yregion(-DU_MAX, DU_MAX, 'FaceColor', [0.8 0.9 0.8], 'FaceAlpha', 0.5, 'DisplayName', 'Feasible Region');
yline(DU_MAX, 'r--', 'LineWidth', 1.5, 'DisplayName', '+\Delta U_{MAX}');
yline(-DU_MAX, 'r--', 'LineWidth', 1.5, 'DisplayName', '-\Delta U_{MAX}');

% Plot actual rate
stairs(t, du_hist, 'k-', 'LineWidth', 1.2, 'DisplayName', '\Delta W_f Rate');

xlabel('Time [s]');
ylabel('Rate of Change \Delta W_f [pps/step]');
title(sprintf('Input Rate Constraints (\\pm %.2f pps/step)', DU_MAX));
legend('Location', 'best');

% Check violations
du_viol_high = sum(du_hist > DU_MAX + 1e-5);
du_viol_low  = sum(du_hist < -DU_MAX - 1e-5);
fprintf('--- Input Rate (\\Delta W_f) ---\n');
if du_viol_high > 0; fprintf('  [WARNING] Exceeded +DU_MAX %d times!\n', du_viol_high); else; fprintf('  [OK] Never exceeded +DU_MAX.\n'); end
if du_viol_low > 0;  fprintf('  [WARNING] Exceeded -DU_MAX %d times!\n', du_viol_low); else; fprintf('  [OK] Never exceeded -DU_MAX.\n'); end
fprintf('\n');
