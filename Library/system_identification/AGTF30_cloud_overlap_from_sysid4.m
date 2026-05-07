% ========================================================================
% Cloud overlap analysis for AGTF30 system identification results.
%
% What this script does:
%   1) Loads the interpolated LPV tables produced by sys_id_4.
%   2) Loads the saved identification clouds (out_Dyn) used to build them.
%   3) Plots representative clouds in state space / output space.
%   4) Quantifies overlap between adjacent clouds with a 2-D occupancy score.
%
% The goal is to show where neighboring operating-point clouds overlap and
% where the data-driven model needs patching or interpolation.
% ========================================================================

clearvars; close all; clc

proj = currentProject;
cd(proj.RootFolder)

%% Settings
Ts = 0.015;
tstart = 19.5;
nbins = 50;          % histogram bins for overlap score
n_show = 8;          % representative clouds to visualize
use_outputs = true;   % true: N_L vs N_H, false: N_L vs Wf
local_center = [];    % optionally choose a local cloud window center

%% Inputs / files
model_file = fullfile(proj.RootFolder,'Params','tables','models','AGTF30_model_tables__1_input.mat');
data_folder = fullfile(proj.RootFolder,'Params','test_files_fuel_only','test_identification','AGTF30');
wf_file = fullfile(proj.RootFolder,'Params','test_files_fuel_only','test_identification','wf_values.mat');

if ~isfile(model_file)
    error('Missing sys_id_4 model file: %s',model_file);
end
if ~isfolder(data_folder)
    error('Missing data folder: %s',data_folder);
end
if ~isfile(wf_file)
    error('Missing wf_values file: %s',wf_file);
end

models = load(model_file);
wf_data = load(wf_file);
wf_values = wf_data.Wf_values(:);

filelist = dir(fullfile(data_folder,'*.mat'));
if isempty(filelist)
    error('No .mat files found in %s',data_folder);
end
[~,idx_sort] = sort({filelist.name});
filelist = filelist(idx_sort);

n_used = min(numel(filelist), numel(wf_values));
filelist = filelist(1:n_used);
wf_values = wf_values(1:n_used);

%% Load clouds
clouds = cell(n_used,1);
cloud_desc = strings(n_used,1);

for k = 1:n_used
    sample = load(fullfile(filelist(k).folder,filelist(k).name));
    out_dyn = get_out_dyn(sample);

    t_raw = signal_time(out_dyn.eng.S0.W);
    t = (tstart:Ts:t_raw(end)).';

    if use_outputs
        x1 = signal_data(out_dyn.eng.Shaft.N_Fan);
        x2 = signal_data(out_dyn.eng.Shaft.N_HPC);
        x1_i = interp1(t_raw,x1,t,'linear','extrap');
        x2_i = interp1(t_raw,x2,t,'linear','extrap');
        clouds{k} = [x1_i, x2_i];
        cloud_desc(k) = sprintf('Wf = %.3f',wf_values(k));
    else
        x1 = signal_data(out_dyn.eng.Shaft.N_Fan);
        u  = signal_data(out_dyn.cntrl.Wfact);
        x1_i = interp1(t_raw,x1,t,'linear','extrap');
        u_i  = interp1(t_raw,u,t,'linear','extrap');
        clouds{k} = [u_i, x1_i];
        cloud_desc(k) = sprintf('Wf = %.3f',wf_values(k));
    end
end

%% Plot 1: representative clouds with a smooth LPV reference curve
if isempty(local_center)
    local_center = max(1,round(n_used/2));
end
half_win = floor(n_show/2);
idx_show = max(1,local_center-half_win):min(n_used,local_center+half_win);
idx_show = unique(idx_show);
figure('Name','Representative clouds from identification runs')
hold on
if use_outputs
    xlabel('N_L [rpm]')
    ylabel('N_H [rpm]')
    title('Cloud overlap in output space')
else
    xlabel('W_f [pps]')
    ylabel('N_L [rpm]')
    title('Cloud overlap in input-output space')
end
grid on

cmap = turbo(numel(idx_show));
for ii = 1:numel(idx_show)
    k = idx_show(ii);
    pts = clouds{k};
    plot(pts(:,1),pts(:,2),'-','Color',0.55*cmap(ii,:),'LineWidth',1.1)
    scatter(pts(1,1),pts(1,2),18,cmap(ii,:),'filled')
end
legend(arrayfun(@(k) sprintf('cloud %d',k),idx_show,'UniformOutput',false),'Location','best')

%% Plot 2: sys_id_4 reference curve if available (A(1,1) shown as example)
% This is not the main cloud plot; it is only a visual reference that the
% tables are smoothly varying.
if isfield(models,'matrices') && isfield(models.matrices,'A')
    A_tab = models.matrices.A;
elseif isfield(models,'A_dt')
    A_tab = models.A_dt;
else
    A_tab = [];
end

if ~isempty(A_tab)
    rho = extract_numeric_axis(models.offset.u);
    rho = rho(:);
    rho_es = linspace(min(rho),max(rho),200).';
    A11 = squeeze(A_tab(1,1,:));
    if size(A11,1) > 1
        if numel(A11) == numel(rho)
            A11_i = interp1(rho,A11,rho_es,'pchip','extrap');
        else
            A11_i = []; % skip if the axis is not compatible
        end
        figure('Name','sys_id_4 coefficient reference')
        plot(rho,A11,'o','LineWidth',1.2)
        hold on
        if ~isempty(A11_i)
            plot(rho_es,A11_i,'LineWidth',1.5)
        end
        grid on
        xlabel('W_f [pps]')
        ylabel('A(1,1)')
        title('sys_id_4 interpolated model coefficient reference')
        if isempty(A11_i)
            legend('Table points','Location','best')
        else
            legend('Table points','Interpolated','Location','best')
        end
    end
end

%% Overlap scores for adjacent clouds
adj_score = nan(n_used-1,1);
adj_common = nan(n_used-1,1);
for k = 1:n_used-1
    [adj_score(k), adj_common(k)] = cloud_overlap_score(clouds{k}, clouds{k+1}, nbins);
end

figure('Name','Adjacent cloud overlap score')
plot(wf_values(1:end-1),adj_score,'-o','LineWidth',1.3)
grid on
xlabel('W_f at lower cloud [pps]')
ylabel('Overlap coefficient')
title(sprintf('Adjacent cloud overlap (2-D histogram, %d bins)',nbins))

figure('Name','Adjacent cloud common occupancy')
plot(wf_values(1:end-1),adj_common,'-s','LineWidth',1.3)
grid on
xlabel('W_f at lower cloud [pps]')
ylabel('Common occupied-bin fraction')
title('How much the cloud areas occupy the same bins')

%% Pairwise comparison for the most overlapping and least overlapping neighbors
[~,idx_min] = min(adj_score);
[~,idx_max] = max(adj_score);

figure('Name','Most and least overlapping neighbors')
tiledlayout(1,2,'Padding','compact','TileSpacing','compact')

for jj = 1:2
    if jj == 1
        k = idx_max;
        fig_title = sprintf('Most overlap: clouds %d and %d',k,k+1);
    else
        k = idx_min;
        fig_title = sprintf('Least overlap: clouds %d and %d',k,k+1);
    end

    nexttile
    hold on; grid on
    if use_outputs
        xlabel('N_L [rpm]')
        ylabel('N_H [rpm]')
    else
        xlabel('W_f [pps]')
        ylabel('N_L [rpm]')
    end
    title(sprintf('%s (score = %.3f)',fig_title,adj_score(k)))

    plot(clouds{k}(:,1),clouds{k}(:,2),'b','LineWidth',1.1)
    plot(clouds{k+1}(:,1),clouds{k+1}(:,2),'r','LineWidth',1.1)
    legend(sprintf('cloud %d',k),sprintf('cloud %d',k+1),'Location','best')
end

%% Summary table
fprintf('\n=== Cloud overlap summary ===\n');
for k = 1:n_used-1
    fprintf('pair %2d-%2d | Wf = %.3f -> %.3f | overlap = %.4f | common bins = %.4f\n', ...
        k, k+1, wf_values(k), wf_values(k+1), adj_score(k), adj_common(k));
end
fprintf('\nHighest overlap pair: %d-%d (score %.4f)\n',idx_max,idx_max+1,adj_score(idx_max));
fprintf('Lowest overlap pair:  %d-%d (score %.4f)\n',idx_min,idx_min+1,adj_score(idx_min));

%% Helpers
function out_dyn = get_out_dyn(data)
    if isstruct(data) && isfield(data,'out_Dyn')
        out_dyn = data.out_Dyn;
        return
    end
    error('Could not find out_Dyn in loaded data file.');
end

function y = signal_data(sig)
    if isa(sig,'timeseries')
        y = sig.Data(:);
    elseif isstruct(sig) && isfield(sig,'Data')
        y = sig.Data(:);
    else
        y = sig(:);
    end
    y = double(y);
end

function t = signal_time(sig)
    if isa(sig,'timeseries')
        t = sig.Time(:);
    elseif isstruct(sig) && isfield(sig,'Time')
        t = sig.Time(:);
    else
        error('Could not extract signal time.');
    end
    t = double(t);
end

function [overlap_score, common_fraction] = cloud_overlap_score(P, Q, nbins)
    % 2-D occupancy overlap using histogram bins.
    all_pts = [P; Q];
    xmin = min(all_pts(:,1)); xmax = max(all_pts(:,1));
    ymin = min(all_pts(:,2)); ymax = max(all_pts(:,2));

    if xmin == xmax
        xmax = xmin + 1;
    end
    if ymin == ymax
        ymax = ymin + 1;
    end

    edges1 = linspace(xmin,xmax,nbins+1);
    edges2 = linspace(ymin,ymax,nbins+1);

    Hp = histcounts2(P(:,1),P(:,2),edges1,edges2,'Normalization','probability');
    Hq = histcounts2(Q(:,1),Q(:,2),edges1,edges2,'Normalization','probability');

    % Histogram intersection coefficient: sum(min(p,q))
    overlap_score = sum(min(Hp(:),Hq(:)));

    % Fraction of occupied bins shared by both clouds
    occ_p = Hp > 0;
    occ_q = Hq > 0;
    common_fraction = nnz(occ_p & occ_q) / max(nnz(occ_p | occ_q),1);
end

function x = extract_numeric_axis(v)
    if isnumeric(v)
        x = v;
        return
    end
    if isstruct(v)
        if isfield(v,'dat')
            x = v.dat;
            return
        end
        if isfield(v,'Data')
            x = v.Data;
            return
        end
    end
    if isa(v,'timeseries')
        x = v.Data;
        return
    end
    error('Unsupported axis type for sys_id_4 coefficient plot.');
end
