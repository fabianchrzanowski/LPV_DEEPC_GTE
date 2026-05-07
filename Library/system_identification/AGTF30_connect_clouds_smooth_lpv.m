% ========================================================================
% Build and inspect a smooth LPV model by explicitly connecting adjacent
% identification clouds with linear blending in scheduling variable Wf.
%
% This script does NOT run Simulink. It uses saved ID tables/data only.
% ========================================================================

clearvars; close all; clc

proj = currentProject;
cd(proj.RootFolder)

Ts = 0.015;
use_demo_data = true;

id_file = fullfile(proj.RootFolder,'Params','tables','matrices_DE_wf.mat');
wf_file = fullfile(proj.RootFolder,'Params','test_files_fuel_only','test_identification','wf_values.mat');
demo_data_folder = fullfile(proj.RootFolder,'Params','test_files_fuel_only','test_identification','AGTF30');

if ~isfile(id_file)
    error('Missing ID file: %s',id_file);
end
if ~isfile(wf_file)
    error('Missing wf_values file: %s',wf_file);
end

models = load(id_file);
wf_data = load(wf_file);

A_dt = models.A_dt;
B_dt = models.B_dt;
C_dt = models.C_dt;
D_dt = models.D_dt;
offset = models.offset;
wf_knots = wf_data.Wf_values(:);

n_op = size(A_dt,3);
if numel(wf_knots) ~= n_op
    n_op = min(numel(wf_knots),n_op);
    wf_knots = wf_knots(1:n_op);
    A_dt = A_dt(:,:,1:n_op);
    B_dt = B_dt(:,:,1:n_op);
    C_dt = C_dt(:,:,1:n_op);
    D_dt = D_dt(:,:,1:n_op);
    offset.x = offset.x(:,1:n_op);
    offset.y = offset.y(:,1:n_op);
    offset.z = offset.z(:,1:n_op);
    offset.u = offset.u(:,1:n_op);
end

% ------------------------------------------------------------------------
% Expression used to connect adjacent clouds (explicitly):
% For rho in [rho_i, rho_{i+1}], alpha = (rho-rho_i)/(rho_{i+1}-rho_i)
% M(rho) = (1-alpha) M_i + alpha M_{i+1}
% where M can be A, B, C, D and offsets.
% ------------------------------------------------------------------------

rho_grid = linspace(min(wf_knots),max(wf_knots),400).';
rho_spec = nan(size(rho_grid));

for k = 1:numel(rho_grid)
    rho = rho_grid(k);
    [A_rho,~,~,~,~,~,~,~] = eval_smooth_lpv(rho,wf_knots,A_dt,B_dt,C_dt,D_dt,offset);
    rho_spec(k) = max(abs(eig(A_rho)));
end

figure('Name','Smooth LPV from connected clouds')
subplot(2,1,1)
plot(rho_grid,rho_spec,'LineWidth',1.4)
hold on
yline(1,'r--','|\lambda|=1')
plot(wf_knots,arrayfun(@(i) max(abs(eig(A_dt(:,:,i)))),1:n_op),'ko','MarkerSize',4)
grid on
xlabel('W_f [pps]')
ylabel('Spectral radius')
title('Stability trend of smooth cloud-connected model')
legend('Interpolated model','Stability boundary','Original cloud poles','Location','best')

subplot(2,1,2)
A11 = squeeze(A_dt(1,1,:));
A11_grid = nan(size(rho_grid));
for k = 1:numel(rho_grid)
    [A_rho,~,~,~,~,~,~,~] = eval_smooth_lpv(rho_grid(k),wf_knots,A_dt,B_dt,C_dt,D_dt,offset);
    A11_grid(k) = A_rho(1,1);
end
plot(wf_knots,A11,'o','LineWidth',1.2,'DisplayName','Cloud values')
hold on
plot(rho_grid,A11_grid,'LineWidth',1.4,'DisplayName','Connected (smooth)')
grid on
xlabel('W_f [pps]')
ylabel('A(1,1)')
title('Example coefficient continuity')
legend('Location','best')

% Optional trajectory demo using recorded Wf from one dataset
if use_demo_data
    filelist = dir(fullfile(demo_data_folder,'*.mat'));
    if ~isempty(filelist)
        [~,idx_sort] = sort({filelist.name});
        filelist = filelist(idx_sort);
        data = load(fullfile(filelist(1).folder,filelist(1).name));
        out_dyn = get_out_dyn(data);

        t = signal_time(out_dyn.eng.S0.W);
        u = signal_data(out_dyn.cntrl.Wfact);

        % States used in identification: [N_Fan, N_HPC, S21.Pt, S25.Pt]
        psi_meas = [signal_data(out_dyn.eng.Shaft.N_Fan), ...
                    signal_data(out_dyn.eng.Shaft.N_HPC), ...
                    signal_data(out_dyn.eng.S21.Pt), ...
                    signal_data(out_dyn.eng.S25.Pt)];

        % Resample to Ts grid
        t_u = (t(1):Ts:t(end)).';
        u_i = interp1(t,u,t_u,'linear','extrap');
        psi_i = interp1(t,psi_meas,t_u,'linear','extrap');

        psi_hat = nan(size(psi_i));
        psi_hat(1,:) = psi_i(1,:);

        for k = 1:numel(t_u)-1
            rho = u_i(k);
            [A_rho,B_rho,~,~,x0,y0,~,u0] = eval_smooth_lpv(rho,wf_knots,A_dt,B_dt,C_dt,D_dt,offset);
            psi0 = [x0; y0];

            % Deviation-form one-step propagation
            xk = psi_hat(k,:).' - psi0;
            uk = u_i(k) - u0;
            xkp1 = psi0 + A_rho*xk + B_rho*uk;
            psi_hat(k+1,:) = xkp1.';
        end

        figure('Name','Connected-cloud model demo (no patch switching)')
        subplot(2,1,1)
        plot(t_u,psi_i(:,1),'k','LineWidth',1.2)
        hold on
        plot(t_u,psi_hat(:,1),'b--','LineWidth',1.1)
        grid on
        xlabel('Time [s]')
        ylabel('N_L [rpm]')
        legend('Measured','Connected-cloud model','Location','best')

        subplot(2,1,2)
        plot(t_u,psi_i(:,2),'k','LineWidth',1.2)
        hold on
        plot(t_u,psi_hat(:,2),'b--','LineWidth',1.1)
        grid on
        xlabel('Time [s]')
        ylabel('N_H [rpm]')

        rmse_nl = sqrt(mean((psi_hat(:,1)-psi_i(:,1)).^2));
        rmse_nh = sqrt(mean((psi_hat(:,2)-psi_i(:,2)).^2));
        fprintf('\nConnected-cloud demo RMSE: N_L = %.4f, N_H = %.4f\n',rmse_nl,rmse_nh);
    else
        warning('No demo .mat files found in %s',demo_data_folder);
    end
end

fprintf('\nDone. This script builds one continuous LPV model by explicit adjacent-cloud blending.\n');

%% Helpers
function [A,B,C,D,x0,y0,z0,u0] = eval_smooth_lpv(rho,wf_knots,A_all,B_all,C_all,D_all,offset)
    [i0,i1,alpha] = blend_index(rho,wf_knots);

    A = (1-alpha)*A_all(:,:,i0) + alpha*A_all(:,:,i1);
    B = (1-alpha)*B_all(:,:,i0) + alpha*B_all(:,:,i1);
    C = (1-alpha)*C_all(:,:,i0) + alpha*C_all(:,:,i1);
    D = (1-alpha)*D_all(:,:,i0) + alpha*D_all(:,:,i1);

    x0 = (1-alpha)*offset.x(:,i0) + alpha*offset.x(:,i1);
    y0 = (1-alpha)*offset.y(:,i0) + alpha*offset.y(:,i1);
    z0 = (1-alpha)*offset.z(:,i0) + alpha*offset.z(:,i1);
    u0 = (1-alpha)*offset.u(:,i0) + alpha*offset.u(:,i1);
end

function [i0,i1,alpha] = blend_index(rho,wf_knots)
    if rho <= wf_knots(1)
        i0 = 1; i1 = 2; alpha = 0;
        return
    end
    if rho >= wf_knots(end)
        i0 = numel(wf_knots)-1; i1 = numel(wf_knots); alpha = 1;
        return
    end

    i0 = find(wf_knots <= rho,1,'last');
    i1 = i0 + 1;
    alpha = (rho - wf_knots(i0)) / (wf_knots(i1) - wf_knots(i0));
end

function out_dyn = get_out_dyn(data)
    if isstruct(data) && isfield(data,'out_Dyn')
        out_dyn = data.out_Dyn;
        return
    end
    error('Could not locate out_Dyn in supplied data.');
end

function x = signal_data(sig)
    if isa(sig,'timeseries')
        x = sig.Data(:);
    elseif isstruct(sig) && isfield(sig,'Data')
        x = sig.Data(:);
    else
        x = sig(:);
    end
    x = double(x);
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
