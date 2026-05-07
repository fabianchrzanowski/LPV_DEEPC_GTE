% ========================================================================
% LPV system id: validation at each equilibrium point
% 
% ========================================================================
%
%
% Pablo Baldivieso Monasterios
% The University of Sheffield
% 04/07/2024


clearvars;close all;clc

plot_options;
proj = currentProject; % proj = currentProject;
cd(proj.RootFolder)

% location of data generated
loc      = fullfile(proj.RootFolder,"/Params/tables/");
folder_content = what(loc);
filelist = fullfile(folder_content.path,folder_content.mat);

model = 'AGTF30SysDyn';
MWS.engName = 'AGTF30';
MWS = setup_Controller(MWS);
MWS = setup_AllEng(MWS);

% load bus
evalin('base','load(''Eng_Bus.mat'')');

% sampling time
Ts = 0.015;
MWS.In.Ts = Ts;

% testing time
Tsim = 30;

% simulation time
MWS.In.Tsim = Tsim;

%% Environmental conditions. To be kept constant across the experiments
MWS.In.Alt  = timeseries([0;0],[0 Tsim],'Name','Alt');
MWS.In.MN   = timeseries([0;0],[0 Tsim],'Name','MN');
MWS.In.dT   = timeseries([0;0],[0 Tsim],'Name','dT');


% load Wf files
load(fullfile(proj.RootFolder,...
    "/Params/test_files_fuel_only/test_identification/wf_values.mat"))

models = load(filelist{1});


n_set = 10;   % number of intervals for the system to settleou
dt_set = 1;
dt_trn = 0.001;
instants = 0:dt_set:Tsim;
dx_Wf  = 0.1 * [zeros(1,n_set),ones(1,n_set),zeros(1,n_set)];

t_model = 0:Ts:Tsim;



for i = 1:numel(Wf_values)
    % initial condition computation
    MWS.In.Wf   = timeseries([Wf_values(i);Wf_values(i)],[0 Tsim],'Name','Wf');
    MWS.In.dVBV = timeseries([0;0],[0 Tsim],'Name','dVBV');
    MWS.In.dNz  = timeseries([0;0],[0 Tsim],'Name','dNz');
    % 
    MWS = AGTF30_initial_conditions(MWS);

    x0     = Wf_values(i);
    values = max(min(dx_Wf + x0,20),0.2);
    [x,t]  = signal_transition_times(values,instants,dt_trn,x0);
    Wf     = timeseries(x,t);
    MWS.In.Wf = Wf;

    tic
    out    = sim(model);
    toc;


    % engine model tape
    % differential states
    x_tape(:,1) = out_Dyn.eng.Shaft.N_Fan.Data;
    x_tape(:,2) = out_Dyn.eng.Shaft.N_HPC.Data;

    % algebraic states
    y_tape(:,1) = out_Dyn.eng.S21.Pt.Data;
    y_tape(:,2) = out_Dyn.eng.S25.Pt.Data;

    % outputs
    z_tape(:,1) = out_Dyn.eng.Shaft.N_Fan.Data;
    z_tape(:,2) = out_Dyn.eng.Shaft.N_HPC.Data;

    z_tape(:,3) = out_Dyn.eng.S2.W.Data;
    z_tape(:,4) = out_Dyn.eng.S25.W.Data;
    
    z_tape(:,5) = out_Dyn.eng.S25.Pt.Data;
    z_tape(:,6) = out_Dyn.eng.S25.Tt.Data;

    z_tape(:,7)  = out_Dyn.eng.S21.Pt.Data;
    z_tape(:,8) = out_Dyn.eng.S21.Tt.Data;
    z_tape(:,9) = out_Dyn.eng.S24.Pt.Data;
    z_tape(:,10) = out_Dyn.eng.S24.Tt.Data;

    z_tape(:,11) = out_Dyn.eng.S36.Pt.Data;
    z_tape(:,12) = out_Dyn.eng.S4.Tt.Data;
    z_tape(:,13) = out_Dyn.eng.S4.FAR.Data;
    z_tape(:,14) = out_Dyn.eng.SM.SMFan.Data;
    z_tape(:,15) = out_Dyn.eng.SM.SMHPC.Data;


    u_tape = out_Dyn.cntrl.Wfact.Data;


    clearvars x y z
    % initial conditions

    % differential states 
    x(1,1) = out_Dyn.eng.Shaft.N_Fan.Data(1);
    x(2,1) = out_Dyn.eng.Shaft.N_HPC.Data(1);

    % algebraic states
    y(1,1) = out_Dyn.eng.S21.Pt.Data(1);
    y(2,1) = out_Dyn.eng.S25.Pt.Data(1);

    
    
    u0 = models.offset.u(i);
    x0 = models.offset.x(:,i);
    y0 = models.offset.y(:,i);
    z0 = models.offset.z(:,i);

    A = models.A_dt(:,:,i);
    B = models.B_dt(:,:,i);
    C = models.C_dt(:,:,i);
    D = models.D_dt(:,:,i);


    nx = size(x0,1);
    ny = size(y0,1);
    nu = size(u0,1);
    nz = size(z0,1);

    for kk = 0:numel(t_model)-1
        
        z = z0 + C*[(x-x0);(y-y0)] + D*(u_tape(kk+1)-u0);
        
        x_model(kk+1,:) = x';
        y_model(kk+1,:) = y';
        z_model(kk+1,:) = z';
        
        
        aux = [x0;y0] + A*[(x-x0);(y-y0)] + B*(u_tape(kk+1)-u0);
        x = aux(1:nx);
        y = aux(1+nx:nx+ny);
    end
    

    % compute errors
    
    e_x = abs(x_model-x_tape)./abs(x_tape)*100;
    e_y = abs(y_model-y_tape)./abs(y_tape)*100;
    e_z = abs(z_model-z_tape)./abs(z_tape)*100;

    fig_tape = figure(1+2*(i-1));
    fig_tape.Name = ['time response at Wf = ' num2str(Wf_values(i)) 'pps'];
    subplot(311)
    plot(t_model,x_tape(:,2))
    hold on
    plot(t_model,x_model(:,2))
    subplot(311)
    xlabel('$t$')
    ylabel(['$' models.sys_eng.StateName{1} '$'])

    subplot(312)
    plot(t_model,x_tape(:,1))
    hold on
    plot(t_model,x_model(:,1))
    xlabel('$t$')
    ylabel(['$' models.sys_eng.StateName{2} '$'])

    subplot(313)
    plot(t_model,z_tape(:,3))
    hold on
    plot(t_model,z_model(:,3))
    xlabel('$t$')
    ylabel(['$' models.sys_eng.OutputName{3} '$'])

    fig_error = figure(2*i);
    fig_error.Name = ['Error histogram at Wf = ' num2str(Wf_values(i)) 'pps'];
    subplot(131)
    for j = 1:nx
    histogram(e_x(:,j),'Normalization','probability');
    hold on
    end
    xlabel('$error~\%$')
    ylabel('$\%~of~time$')
    title('State differential')
    legend(models.sys_eng.StateName(1:nx),'Interpreter','tex')
    
    subplot(132)
    for j = 1:ny
    histogram(e_y(:,j),'Normalization','probability');
    hold on
    end
    xlabel('$error~\%$')
    ylabel('$\%~of~time$')
    title('State algebraic')
    legend(models.sys_eng.StateName(1+nx:ny+nx),'Interpreter','tex')

    subplot(133)
    for j = 1:nz
    histogram(e_z(:,j),'Normalization','probability');
    hold on
    end
    xlabel('$error~\%$')
    ylabel('$\%~of~time$')
    title('Output')
    legend(models.sys_eng.OutputName,'Interpreter','tex')

end



