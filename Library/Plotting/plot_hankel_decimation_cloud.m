% ========================================================================
% plot_hankel_decimation_cloud.m
% Visualise the downsampled Hankel PRBS data cloud used by DeePC
% ========================================================================
clearvars; close all; clc

try proj = currentProject; catch; proj = openProject('GTE.prj'); end
cd(proj.RootFolder)

% USER CONFIGURABLE PARAMETERS
NUM_HANKELS = 8;        % Select how many operating points to keep (max 32)
NUM_COLUMNS = 3000;     % Select how many columns per OP to keep (max ~3000)
TRACK_NAME  = 'N_H';    % Output to plot

%% Load models and Wf values
wf_data = load(fullfile(proj.RootFolder,'Params','test_files_fuel_only','test_identification','wf_values.mat'));
wf_values = wf_data.Wf_values(:)';

% 1. Hankel Decimation Logic (same as DeePC)
filelist = dir(fullfile(proj.RootFolder,'Params','test_files_fuel_only','test_identification','AGTF30','*.mat'));
[~,idx_sort] = sort({filelist.name});
filelist = filelist(idx_sort);

if NUM_HANKELS < numel(filelist)
    keep_op_idx = unique(round(linspace(1, numel(filelist), NUM_HANKELS)));
else
    keep_op_idx = 1:numel(filelist);
end

filelist = filelist(keep_op_idx);
hankel_wf = wf_values(keep_op_idx);
n_op = numel(filelist);

cloud_colors = turbo(n_op);

all_u_prbs = [];
all_y_prbs = [];
all_c_prbs = [];

fprintf('Loading data for %d downsampled operating points...\n', n_op);

for index = 1:n_op
    data = load(fullfile(filelist(index).folder, filelist(index).name));
    
    out_dyn = fetch_out_dyn(data);
    u_raw = extract_input_trace(out_dyn);
    y_raw = measure_tracking_series(out_dyn, TRACK_NAME);
    
    n_common = min(numel(u_raw), numel(y_raw));
    u_raw = u_raw(1:n_common);
    y_raw = y_raw(1:n_common);
    
    % 2. Column Decimation Logic (same as DeePC)
    nc = numel(u_raw);
    if nc > NUM_COLUMNS
        keep_col_idx = unique(round(linspace(1, nc, NUM_COLUMNS)));
        u_prbs = u_raw(keep_col_idx);
        y_prbs = y_raw(keep_col_idx);
    else
        u_prbs = u_raw;
        y_prbs = y_raw;
    end
    
    all_u_prbs = [all_u_prbs; u_prbs(:)];
    all_y_prbs = [all_y_prbs; y_prbs(:)];
    
    base_color = cloud_colors(index, :);
    all_c_prbs = [all_c_prbs; repmat(base_color, length(u_prbs), 1)];
end

%% Plot the Cloud
figure('Name', sprintf('Data Cloud (%d Hankels, %d Cols)', n_op, NUM_COLUMNS));
scatter(all_u_prbs, all_y_prbs, 8, all_c_prbs, 'filled', 'MarkerFaceAlpha', 0.2)
grid on; box on;
title(sprintf('PRBS Data Cloud (Decimated to %d Hankels, %d Columns)', n_op, NUM_COLUMNS));
xlabel('Fuel Flow W_f (pps)');
ylabel(sprintf('%s', TRACK_NAME), 'Interpreter', 'none');

% Plot the equilibrium line as reference
hold on;
y_means = zeros(n_op, 1);
for i = 1:n_op
    % Retrieve means based on color blocks
    block_size = length(all_u_prbs)/n_op;
    y_means(i) = mean(all_y_prbs((i-1)*block_size + 1 : i*block_size));
end
plot(hankel_wf, y_means, 'kx-', 'LineWidth', 2, 'MarkerSize', 8)
legend('Decimated PRBS Data', 'Equilibrium Operating Points', 'Location', 'Best');

fprintf('Done! Plotted %d points.\n', numel(all_u_prbs));

%% ===== Helpers =====
function out_dyn = fetch_out_dyn(data)
    if isfield(data, 'out_Dyn')
        out_dyn = data.out_Dyn;
    else
        out_dyn = data;
    end
end

function u = extract_input_trace(out_dyn)
    if isfield(out_dyn, 'cntrl')
        u = out_dyn.cntrl.Wfact.Data(:);
    else
        u = out_dyn.u(:);
    end
end

function y = measure_tracking_series(out_dyn, track_name)
    if isfield(out_dyn, 'eng')
        switch track_name
            case 'N_L', y = out_dyn.eng.Shaft.N_Fan.Data(:);
            case 'N_H', y = out_dyn.eng.Shaft.N_HPC.Data(:);
            otherwise, error('Unknown track_name');
        end
    else
        % Fallback for direct array formats (e.g. y(:,2) = N_H)
        if strcmp(track_name, 'N_H')
            y = out_dyn.y(:,2);
        else
            y = out_dyn.y(:,1);
        end
    end
end
