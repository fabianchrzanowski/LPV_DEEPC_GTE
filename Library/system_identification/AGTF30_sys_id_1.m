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
proj = currentProject; % proj = currentProject;
cd(proj.RootFolder)

%% Save folder
savefolder = 'test_identification/AGTF30';
prefix = fullfile(proj.RootFolder,"Params","test_files_fuel_only",savefolder);
if ~exist(prefix,"dir"); mkdir(prefix); end
addFile(proj,prefix);
%%

model = 'AGTF30SysDyn';
MWS.engName = 'AGTF30';
MWS = setup_Controller(MWS);
MWS = setup_AllEng(MWS);

% load bus
evalin('base','load(''Eng_Bus.mat'')');

% sampling time
Ts = 0.015;
MWS.In.Ts = Ts;



%% Design of experiments

% number of actuation variables
nx = 3;

% individual signals
dx_indiv = kron(eye(nx),[1,0]);

n_set = 200;   % number of intervals for the system to settleou
n_trn = 700;  % number of transients
% Time profile
dt_set = 0.1;  % time for each change to settle in seconds
dt_trn = 0.001; % transition time in the steps to avoid interpolation issues
Tsim = dt_set * (n_set + n_trn + length(dx_indiv) + 2); % total simulation time
instants = 0:dt_set:Tsim;

% simulation time
MWS.In.Tsim = Tsim;

%% Environmental conditions. To be kept constant across the experiments
MWS.In.Alt  = timeseries([0;0],[0 Tsim],'Name','Alt');
MWS.In.MN   = timeseries([0;0],[0 Tsim],'Name','MN');
MWS.In.dT   = timeseries([0;0],[0 Tsim],'Name','dT');



%%

Wf_values = 0.4:0.05:1.95;



for i = 1:numel(Wf_values)



    % initial condition computation
    MWS.In.Wf   = timeseries([Wf_values(i);Wf_values(i)],[0 Tsim],'Name','Wf');
    MWS.In.dVBV = timeseries([0;0],[0 Tsim],'Name','dVBV');
    MWS.In.dNz  = timeseries([0;0],[0 Tsim],'Name','dNz');
    % 
    MWS = AGTF30_initial_conditions(MWS);

    dx_Wf  = 0.1 * [zeros(1,n_set),dx_indiv(1,:),(2*rand(1,n_trn)-1),zeros(1,2)];
    dx_VBV = 0.00 * [zeros(1,n_set),dx_indiv(2,:),(2*rand(1,n_trn)-1),zeros(1,2)];
    dx_Nz  = 000  * [zeros(1,n_set),dx_indiv(3,:),(2*rand(1,n_trn)-1),zeros(1,2)];


    x0     = Wf_values(i);
    values = max(min(dx_Wf + x0,20),0.2);
    [x,t]  = signal_transition_times(values,instants,dt_trn,x0);
    Wf     = timeseries(x,t);
    MWS.In.Wf = Wf;


    x0     = 0;
    values = max(min(dx_VBV + x0,20),-inf);
    [x,t]  = signal_transition_times(values,instants,dt_trn,x0);
    dVBV   = timeseries(x,t);
    MWS.In.dVBV = dVBV;


    x0     = 0;
    values = max(min(dx_Nz + x0,8000),-inf);
    [x,t]  = signal_transition_times(values,instants,dt_trn,x0);
    dNz    = timeseries(x,t);
    MWS.In.dNz = dNz;

    tic
    out    = sim(model);
    toc;

    savename = fullfile(prefix,...
        sprintf('test_id__AGTF30__Wf%04d__.mat',i));

    save(savename,'out_Dyn');
    addFile(proj,savename);

    u_name = fieldnames(out_Dyn.cntrl);
    fig_u = figure(1);
    fig_u.Name = 'Inputs';
    for j = 1:numel(u_name)
        subplot(3,1,j)
        plot(out_Dyn.cntrl.(u_name{j}))
        hold on
    end
    
    fig_x = figure(2);
    fig_u.Name = 'Speeds';
    subplot(311)

    plot(out_Dyn.eng.Shaft.N_LPC)
    hold on
%     plot(out_Dyn.eng.Shaft.N_Fan,'--')
    
    subplot(312)
    plot(out_Dyn.eng.Shaft.N_HPC)
    hold on
    subplot(325)
    plot(out_Dyn.cntrl.Wfact.Data,out_Dyn.eng.Shaft.N_LPC.Data)
    hold on

    subplot(326)
    plot(out_Dyn.cntrl.Wfact.Data,out_Dyn.eng.Shaft.N_HPC.Data)

    hold on
end
   
    subplot(325)
    xlabel('Wf [pps]')
    ylabel('N_L [rpm]','Interpreter','tex')

    subplot(326)
    xlabel('Wf [pps]')
    ylabel('N_H [rpm]','Interpreter','tex')


saveWf = fullfile(prefix,'../wf_values.mat');
save(saveWf,'Wf_values');
addFile(proj,saveWf);

