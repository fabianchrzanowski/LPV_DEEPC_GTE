clearvars; close all; clc;

% Plots the eigenvalues of identified models on the unit circle
% as they migrate around the unit circle -> the system is nonlinear

try proj = currentProject; catch; proj = openProject('GTE.prj'); end
cd(proj.RootFolder);

% Load the LPV state-space matrices
models = load(fullfile(proj.RootFolder, 'Params', 'tables', 'matrices_DE_wf.mat'));
wf_data = load(fullfile(proj.RootFolder, 'Params', 'test_files_fuel_only', 'test_identification', 'wf_values.mat'));
Wf_values = wf_data.Wf_values(:)';

M = models.A_dt;
[n, m, l] = size(M);
ev = nan(max(n, m), l);

for i = 1:l
    ev(:, i) = sort(eig(M(:, :, i)));
end

% Create the complex plane eigenvalue plot
figure();
plot(real(ev)', imag(ev)', 'x', 'LineWidth', 2, 'MarkerSize', 8, 'Color', [0 0.4470 0.7410]); 
grid on;
set(gca, 'FontSize', 12, 'LineWidth', 1);    
xlabel('Real Axis ($Re(\lambda)$)', 'Interpreter', 'latex', 'FontSize', 14);
ylabel('Imaginary Axis ($Im(\lambda)$)', 'Interpreter', 'latex', 'FontSize', 14);

% draw unit circle for discrete time stability reference
hold on;
th = linspace(0, 2*pi, 100);
plot(cos(th), sin(th), 'k--', 'LineWidth', 1);
axis equal; % Ensure circle looks like a circle
xlim([0.6 1.05]); % Focus on the relevant stable region

fprintf('Plot created successfully.\n');
