% ========================================================================
% LPV system id: interpolation of the models
% 
% ========================================================================
%
%
% Pablo Baldivieso Monasterios
% The University of Sheffield
% 04/07/2024


clearvars;close all;clc

plot_options;
try; proj = currentProject; catch; proj = openProject('GTE.prj'); end % try; proj = currentProject; catch; proj = openProject('GTE.prj'); end
cd(proj.RootFolder)


% location of data generated
loc      = fullfile(proj.RootFolder,"/Params/tables/");
folder_content = what(loc);
filelist = fullfile(folder_content.path,folder_content.mat);
suffix = {'1_input'};

f = 1;

% load Wf files
load(fullfile(proj.RootFolder,"/Params/test_files_fuel_only/test_identification/wf_values.mat"))


models = load(filelist{f});

% Sampling time for identification
Ts     = 0.015;

plot_matrix_coefficients(models.A_dt,'A');
plot_matrix_coefficients(models.B_dt,'B');
plot_matrix_coefficients(models.C_dt,'C');
plot_matrix_coefficients(models.D_dt,'D');
plot_eigenvalues(models.A_dt,Wf_values)


rho = models.offset.u(1,:)';
d_rho = 0.001;
rho_es = rho(1):d_rho:rho(end);

nx = size(models.offset.x,1);
O = zeros(nx,length(rho_es));
Xss = interp1(rho',models.offset.x',rho_es,'pchip','extrap')';
offset.x = lookup_array(Xss,rho_es,'state_diff');


ny = size(models.offset.y,1);
O = zeros(ny,length(rho_es));
Yss = interp1(rho',models.offset.y',rho_es,'pchip','extrap')';
offset.y = lookup_array(Yss,rho_es,'state_alg');

nu = size(models.offset.u,1);
O = zeros(nu,length(rho_es));
Uss = interp1(rho',models.offset.u',rho_es,'pchip','extrap')';
aux = size(Uss);
if aux(1) >= aux(2)
    Uss = Uss';
end
offset.u = lookup_array(Uss,rho_es,'input');

nz = size(models.offset.z,1);
O = zeros(nz,length(rho_es));
Zss = interp1(rho',models.offset.z',rho_es,'pchip','extrap')';
offset.z = lookup_array(Zss,rho_es,'output');

LUT.offset = offset;

LUT.names.x = models.sys_eng.StateName(1:nx);
LUT.names.y = models.sys_eng.StateName(1+nx:nx+ny);
LUT.names.z = models.sys_eng.OutputName;
LUT.names.u = models.sys_eng.InputName;


% Interpolation
% A matrix: eigenvalue-based interpolation with MAC mode tracking
eig_opts.method  = 'pchip';
eig_opts.nom_idx = round(numel(rho)/2);  % midpoint eigenvectors
[A, eig_info] = interp_matrix_eig(rho, models.A_dt, rho_es, eig_opts);

% B, C, D: standard element-wise interpolation (no stability concern)
B = interp_matrix(rho,models.B_dt,rho_es);
C = interp_matrix(rho,models.C_dt,rho_es);
D = interp_matrix(rho,models.D_dt,rho_es);

% Storing the matrices
LUT.matrices.A = lookup_matrix(A,rho_es,'A');
LUT.matrices.B = lookup_matrix(B,rho_es,'B');
LUT.matrices.C = lookup_matrix(C,rho_es,'C');
LUT.matrices.D = lookup_matrix(D,rho_es,'D');
LUT.matrices.Ts = Ts;

% Eigenvalues: tracked trajectories from eigenvalue interpolation
h_eig = figure;
h_eig.Name = 'Eigenvalues (tracked from eigenvalue interp)';
n_eig = size(eig_info.eig_mag_eig, 1);
tiledlayout(2,1,"TileSpacing","compact");
nexttile;
    plot(rho_es, eig_info.eig_mag_eig, '-', 'LineWidth', 2);
    ylabel('$|\lambda|$'); grid on; set(gca,'FontSize',14);
    title('Tracked eigenvalue magnitudes');
    xlabel('$\rho$ (W_f)');
    yline(1, '--r', 'LineWidth', 1.5);
    leg = arrayfun(@(i) sprintf('\\lambda_%d', i), 1:n_eig, 'UniformOutput', false);
    legend([leg, {'|\\lambda|=1'}], 'Location','best');
nexttile;
    lam = eig_info.lambda_tracked;  % eigenvalues at identified points
    plot(rho_es, real(eig_info.lambda_interp)', '.', 'MarkerSize', 4);
    hold on;
    plot(rho_es, imag(eig_info.lambda_interp)', '.', 'MarkerSize', 4);
    xlabel('$\rho$ (W_f)');
    ylabel('$Re / Im$'); grid on; set(gca,'FontSize',14);
    title('Real and imaginary parts');

% Comparison: element-wise vs eigenvalue interpolation
A_elemwise = interp_matrix(rho,models.A_dt,rho_es);
plot_eigenvalue_comparison(A_elemwise, rho_es, eig_info);

plot_matrix_coefficients(A,'A (eig interp)');
plot_matrix_coefficients(B,'B');
plot_matrix_coefficients(C,'C');
plot_matrix_coefficients(D,'D');

savename = fullfile(proj.RootFolder,'Params','tables/models',...
    ['AGTF30_model_tables__' suffix{f} '.mat']);
save(savename,'-struct','LUT');
addFile(proj,savename);


 %% Auxiliary functions
function plot_eigenvalues(M,xdata,suffix)
    h = figure;
    h.Name = 'Eigenvalues';
    if nargin < 3; suffix = '';
    else; suffix = [ ' (' suffix ')']; end

    [n,m,l] = size(M);
    ev = nan(max(n,m),l);

    for i = 1:l
        ev(:,i) = sort(eig(M(:,:,i)));
    end
    tiledlayout(2,1,"TileSpacing","compact");
    nexttile;
        plot(xdata,abs(ev),'x','LineWidth',3); ylabel('$abs(\lambda)$')
        grid on; set(gca,'FontSize',14);
        title(['eigenvalues' suffix]);
    nexttile;
        plot(real(ev)',imag(ev)','x','LineWidth',3); 
        xlabel('$Re(\lambda)$');
        ylabel('$Im(\lambda)$')
        grid on; set(gca,'FontSize',14);    

end

function plot_matrix_coefficients(M,name)

[nr,nc,~] = size(M);

h = figure;

h.Name = name;

tiledlayout(nr,nc,"TileSpacing","compact");

for r = 1:nr
    for c = 1:nc
        nexttile;
        plot(squeeze(M(r,c,:)));
        set(gca,'FontSize',14); grid on;
        title(sprintf('(%d,%d)',r,c));
    end
end

end

function plot_eigenvalue_comparison(A_elem, xdata, eig_info)
    % Comparison plot: element-wise vs eigenvalue interpolation
    % Uses tracked eigenvalue trajectories from eig_info for the
    % eigenvalue panel (smooth), and re-computes for element-wise.
    h = figure;
    h.Name = 'Eigenvalue comparison: element-wise vs eigenvalue interp';

    [n,~,l] = size(A_elem);
    ev_elem = nan(n, l);
    for i = 1:l
        ev_elem(:,i) = sort(abs(eig(A_elem(:,:,i))), 'descend');
    end

    % Tracked eigenvalue magnitudes from interp_matrix_eig (smooth)
    ev_eig = eig_info.eig_mag_eig;

    tiledlayout(2,1,"TileSpacing","compact");

    nexttile;
    plot(xdata, ev_elem, '-', 'LineWidth', 1.5); hold on;
    yline(1, '--r', 'LineWidth', 1.5);
    ylabel('$|\lambda|$');
    grid on; set(gca,'FontSize',14);
    title('Element-wise interpolation');
    xlabel('$\rho$ (W_f)');
    leg = arrayfun(@(i) sprintf('\lambda_%d', i), 1:n, 'UniformOutput', false);
    legend([leg, {'|\lambda|=1'}], 'Location','best');

    nexttile;
    plot(xdata, ev_eig, '-', 'LineWidth', 1.5); hold on;
    yline(1, '--r', 'LineWidth', 1.5);
    ylabel('$|\lambda|$');
    grid on; set(gca,'FontSize',14);
    title('Eigenvalue interpolation (tracked trajectories)');
    xlabel('$\rho$ (W_f)');
    legend([leg, {'|\lambda|=1'}], 'Location','best');
end