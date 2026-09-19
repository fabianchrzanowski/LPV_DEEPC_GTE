% ========================================================================
% Plot: Frequency Content Comparison (PRBS vs Step)
%
% Generates a log-log FFT magnitude spectrum of the input excitation
% signal (Wf) to demonstrate that PRBS provides superior persistent 
% excitation
% ========================================================================
clearvars; close all; clc;

try proj = currentProject; catch; proj = openProject('GTE.prj'); end
cd(proj.RootFolder);
addpath(genpath('Library'));
plot_options;

% ===== Configuration =====
op_idx = 16; % Middle operating point (Wf ~ 0.94)
Fs = 1/0.015; % Sampling frequency (Ts = 15ms)

% File paths
prbs_file = fullfile(proj.RootFolder, 'Params', 'test_files_fuel_only', 'test_identification', 'AGTF30', sprintf('test_id__AGTF30__Wf%04d__.mat', op_idx));
step_file = fullfile(proj.RootFolder, 'Params', 'test_files_fuel_only', 'test_identification', 'AGTF30_Step_ALL', sprintf('test_id__AGTF30_Step__Wf%04d__.mat', op_idx));

if ~isfile(prbs_file) || ~isfile(step_file)
    error('Data files missing! Please ensure both PRBS and Step data exist for index %d', op_idx);
end

% ===== Load Data =====
% PRBS
data_prbs = load(prbs_file);
out_dyn_prbs = fetch_out_dyn(data_prbs);
u_prbs = extract_input_trace(out_dyn_prbs);
u_prbs = u_prbs - mean(u_prbs); % Remove DC offset for cleaner FFT

% Step
data_step = load(step_file);
out_dyn_step = fetch_out_dyn(data_step);
u_step = extract_input_trace(out_dyn_step);
u_step = u_step - mean(u_step); % Remove DC offset

% Normalize lengths to match so the FFT resolution is identical
L = min(length(u_prbs), length(u_step));
u_prbs = u_prbs(1:L);
u_step = u_step(1:L);

% ===== Compute FFT =====
% Function for single-sided amplitude spectrum
compute_fft = @(x) get_single_sided_fft(x, L);

[f, P1_prbs] = compute_fft(u_prbs);
[~, P1_step] = compute_fft(u_step);

% Apply smoothing to PRBS spectrum to make the trend clearer on log-log
smooth_span = max(1, floor(L/200));
P1_prbs_smooth = smoothdata(P1_prbs, 'movmean', smooth_span);
P1_step_smooth = smoothdata(P1_step, 'movmean', smooth_span);

% ===== Plotting =====
fig1 = figure('Name', 'FFT Comparison: PRBS vs Step');
ax1 = axes(fig1); hold(ax1, 'on');

% Plot raw spectrums with low alpha for background
loglog(ax1, f, P1_prbs, 'Color', [0 0.4470 0.7410 0.3], 'LineWidth', 0.5, 'HandleVisibility', 'off');
loglog(ax1, f, P1_step, 'Color', [0.8500 0.3250 0.0980 0.3], 'LineWidth', 0.5, 'HandleVisibility', 'off');

% Plot smoothed trends
loglog(ax1, f, P1_prbs_smooth, 'Color', [0 0.4470 0.7410], 'LineWidth', 2.0, 'DisplayName', 'PRBS Spectrum');
loglog(ax1, f, P1_step_smooth, 'Color', [0.8500 0.3250 0.0980], 'LineWidth', 2.0, 'DisplayName', 'Step Spectrum');

% Formatting
grid(ax1, 'on');
ax1.XMinorGrid = 'on';
ax1.YMinorGrid = 'on';
xlim(ax1, [f(2), f(end)]); % Ignore exact zero Hz due to log scale

xlabel(ax1, 'Frequency [Hz]', 'Interpreter', 'latex', 'FontSize', 12);
ylabel(ax1, 'Magnitude $|U(f)|$', 'Interpreter', 'latex', 'FontSize', 12);
title(ax1, 'Frequency Content Comparison: PRBS vs. Step Input ($W_f$)', 'Interpreter', 'latex', 'FontSize', 13);
legend(ax1, 'Location', 'southwest', 'Interpreter', 'latex', 'FontSize', 11);
set(ax1, 'TickLabelInterpreter', 'latex', 'FontSize', 11);

% Optional: Annotate
text(ax1, 0.05, 0.1, 'Broadband Excitation', 'Units', 'normalized', 'Interpreter', 'latex', 'Color', [0 0.4470 0.7410], 'FontSize', 11);
text(ax1, 0.05, 0.9, 'Low-Frequency Dominant', 'Units', 'normalized', 'Interpreter', 'latex', 'Color', [0.8500 0.3250 0.0980], 'FontSize', 11);

fprintf('Plot generated successfully.\n');


% --- Helper Function ---
function [f, P1] = get_single_sided_fft(X, L)
    Fs = 1/0.015; % 15ms sampling
    Y = fft(X);
    P2 = abs(Y/L);
    P1 = P2(1:floor(L/2)+1);
    P1(2:end-1) = 2*P1(2:end-1);
    f = Fs*(0:(L/2))/L;
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

function x = sv(s)
    if isstruct(s)&&isfield(s,'Data'); x = double(s.Data(:));
    elseif isa(s,'timeseries'); x = double(s.Data(:));
    else; x = double(s(:)); end
end
