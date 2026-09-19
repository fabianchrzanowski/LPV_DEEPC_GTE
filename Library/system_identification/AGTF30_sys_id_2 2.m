% ========================================================================
% LPV system id: experiment definition and execution
% 
% ========================================================================
%
%
% Pablo Baldivieso Monasterios
% The University of Sheffield
% 18/12/2023


clearvars;close all;clc

plot_options;
try; proj = currentProject; catch; proj = openProject('GTE.prj'); end % try; proj = currentProject; catch; proj = openProject('GTE.prj'); end
cd(proj.RootFolder)

% location of data generated
loc      = fullfile(proj.RootFolder,"/Params/test_files_fuel_only/test_identification/AGTF30/");
folder_content = what(loc);
filelist = fullfile(folder_content.path,folder_content.mat);

% load Wf files
load(fullfile(proj.RootFolder,"/Params/test_files_fuel_only/test_identification/wf_values.mat"))

% suffix for saving the tables
suffix = 'AGTF30';


unames = {'Wf'};

% outputs
znames = {'N_L','N_H',...
          'W_{20}','W_{25}', ...
          'Pt_{25}','Tt_{25}',...
          'Pt_{21}','Tt_{21}','Pt_{24}','Tt_{24}',...
          'Pt_{36}','Tt_{40}',...
          'FAR_4',...
          'SMF','SMH'};

% states
xnames = {'N_L','N_H'};
ynames = {'Pt_{21}','Pt_{25}'};

% 
nf = length(filelist);
nu = length(unames);
ny = length(ynames);
nx = length(xnames);
nz = length(znames);


% Memory preallocation
A_dt = zeros(nx+ny,nx+ny,nf);
B_dt = zeros(nx+ny,nu,nf);
C_dt = zeros(nz,nx+ny,nf);
D_dt = zeros(nz,nu,nf);


% Sampling time for identification
Ts     = 0.015;
tstart = 19.5;


data_z = struct();
data_x = struct();
data_y = struct();
data_u = struct();
% Identification process per file
for f = 1:nf
    data = load(filelist{f});
    tout = data.out_Dyn.eng.S0.W.Time;

    % time array (trimmed and resampled)
    t = tstart:Ts:tout(end);

    % inputs
    u = data.out_Dyn.cntrl.Wfact;

    % differential states
    x(:,1) = data.out_Dyn.eng.Shaft.N_Fan;
    x(:,2) = data.out_Dyn.eng.Shaft.N_HPC;

    % algebraic states
    y(:,1) = data.out_Dyn.eng.S21.Pt;
    y(:,2) = data.out_Dyn.eng.S25.Pt;

    % outputs
    z(:,1) = data.out_Dyn.eng.Shaft.N_Fan;
    z(:,2) = data.out_Dyn.eng.Shaft.N_HPC;

    z(:,3) = data.out_Dyn.eng.S2.W;
    z(:,4) = data.out_Dyn.eng.S25.W;
    
    z(:,5) = data.out_Dyn.eng.S25.Pt;
    z(:,6) = data.out_Dyn.eng.S25.Tt;

    z(:,7)  = data.out_Dyn.eng.S21.Pt;
    z(:,8) = data.out_Dyn.eng.S21.Tt;
    z(:,9) = data.out_Dyn.eng.S24.Pt;
    z(:,10) = data.out_Dyn.eng.S24.Tt;

    z(:,11) = data.out_Dyn.eng.S36.Pt;
    z(:,12) = data.out_Dyn.eng.S4.Tt;
    z(:,13) = data.out_Dyn.eng.S4.FAR;
    z(:,14) = data.out_Dyn.eng.SM.SMFan;
    z(:,15) = data.out_Dyn.eng.SM.SMHPC;


    % process data
    
    data_z = proc_data(z,znames,t,data_z,f);
    data_x = proc_data(x,xnames,t,data_x,f);
    data_y = proc_data(y,ynames,t,data_y,f);
    data_u = proc_data(u,unames,t,data_u,f);

    

    PSI = [data_x.proc(:,:,f)';data_y.proc(:,:,f)';data_u.proc(:,:,f)'];
    LHS = [PSI(1:end-nu,2:end);data_z.proc(1:end-1,:,f)'];

    regressors = PSI(:,1:end-1);
    N_data = size(regressors,2);

    theta = sdpvar((nx+ny+nz),(nx+ny+nu),'full');
    sol_op = sdpsettings('solver','fmincon');
%     cons = [eye(N_data*(nx+ny));-eye(N_data*(nx+ny))]*kron(regressors',eye(nx+ny))*theta <= 0.001;
%     for j = 1:nx+ny
%         cons = [cons;theta((j*nx+(j-1)*ny+1):(j*(nx+ny))) == 0];
%     end
    cons = [];
    objective = sum(sum(abs(theta))) + sum(sum((LHS - theta*regressors).^2));
    info = optimize(cons,objective,sol_op);
    sol = value(theta);
    SXY = blkdiag(diag(data_x.scale.max(:,f)),diag(data_y.scale.max(:,f)));
    SU  = diag(data_u.scale.max(:,f));
    SZ  = diag(data_z.scale.max(:,f));
    A_dt(:,:,f) = SXY*sol(1:nx+ny,1:nx+ny)*SXY^-1;
    B_dt(:,:,f) = SXY*sol(1:nx+ny,1+nx+ny:nx+ny+nu)*SU;
    C_dt(:,:,f) = SZ*sol(1+nx+ny:nx+ny+nz,1:nx+ny)*SXY^-1;
    D_dt(:,:,f) = SZ*sol(1+nx+ny:nx+ny+nz,1+nx+ny:nx+ny+nu)*SU^-1;

end

offset.x = data_x.offset;
offset.y = data_y.offset;
offset.z = data_z.offset;
offset.u = data_u.offset;

sys_eng = ss(A_dt,B_dt,C_dt,D_dt,Ts,...
    'StateName',[xnames,ynames],'OutputName',znames,'InputName',unames);

prefix   = fullfile(proj.RootFolder,"Params","tables");
save_mat = fullfile(prefix,'/matrices_DE_wf.mat');
save(save_mat,'sys_eng','A_dt','B_dt','C_dt','D_dt','offset');
addFile(proj,save_mat);

