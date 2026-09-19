% ========================================================================
% LPV system id: validation of the LPV
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
loc      = fullfile(proj.RootFolder,"/Params/tables/models");
folder_content = what(loc);
filelist = fullfile(folder_content.path,folder_content.mat);

models = load(filelist{1});


%% Configure 
model = 'AGTF30SysDyn';
MWS.engName = 'AGTF30';
MWS = setup_Controller(MWS);
MWS = setup_AllEng(MWS);

% load bus
evalin('base','load(''Eng_Bus.mat'')');

% sampling time
Ts = 0.015;
MWS.In.Ts = Ts;






%% define validation inputs

n_set  = 5;   % number of intervals for the system to settleou
dt_set = 1;
dt_trn = 0.001;
wf_values = 0.4:0.2:1.9;
dx_Wf  = kron([wf_values flip(wf_values)],ones(1,n_set));
Tsim = dt_set * (size(dx_Wf,2)); % total simulation time
instants = 0:dt_set:Tsim;


t_model = 0:Ts:Tsim;



% simulation time
MWS.In.Tsim = Tsim;

%% Environmental conditions. To be kept constant across the experiments
MWS.In.Alt  = timeseries([0;0],[0 Tsim],'Name','Alt');
MWS.In.MN   = timeseries([0;0],[0 Tsim],'Name','MN');
MWS.In.dT   = timeseries([0;0],[0 Tsim],'Name','dT');
%%

% initial condition computation
MWS.In.Wf   = timeseries([wf_values(1);wf_values(1)],[0 Tsim],'Name','Wf');
MWS.In.dVBV = timeseries([0;0],[0 Tsim],'Name','dVBV');
MWS.In.dNz  = timeseries([0;0],[0 Tsim],'Name','dNz');
%
MWS = AGTF30_initial_conditions(MWS);

%% Run engine

x0     = wf_values(1);
values = max(min(dx_Wf,20),0.2);
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
%% initial conditions

% differential states
x(1,1) = out_Dyn.eng.Shaft.N_Fan.Data(1);
x(2,1) = out_Dyn.eng.Shaft.N_HPC.Data(1);

% algebraic states
y(1,1) = out_Dyn.eng.S21.Pt.Data(1);
y(2,1) = out_Dyn.eng.S25.Pt.Data(1);

%% run the model

nx = size(models.offset.x.dat,1);
ny = size(models.offset.y.dat,1);
nu = size(models.offset.u.dat,1);
nz = size(models.offset.z.dat,1);



for kk = 0:numel(t_model)-1

    rho = u_tape(kk+1);

    A = LPV2mat(models.matrices.A,rho);
    B = LPV2mat(models.matrices.B,rho);
    C = LPV2mat(models.matrices.C,rho);
    D = LPV2mat(models.matrices.D,rho);

    x0 = off2mat(models.offset.x,rho);
    y0 = off2mat(models.offset.y,rho);
    u0 = off2mat(models.offset.u,rho);
    z0 = off2mat(models.offset.z,rho);
    
    z = z0 + C*[(x-x0);(y-y0)] + D*(u_tape(kk+1)-u0);

    x_model(kk+1,:) = x';
    y_model(kk+1,:) = y';
    z_model(kk+1,:) = z';


    aux = [x0;y0] + A*[(x-x0);(y-y0)] + B*(u_tape(kk+1)-u0);
    x = aux(1:nx);
    y = aux(1+nx:nx+ny);
end
    
%%
% compute errors

e_x = abs(x_model-x_tape)./abs(x_tape)*100;
e_y = abs(y_model-y_tape)./abs(y_tape)*100;
e_z = abs(z_model-z_tape)./abs(z_tape)*100;
i = 1;

fig_tape = figure(1+2*(i-1));
fig_tape.Name = ['time response'];
subplot(311)
plot(t_model,x_tape(:,2))
hold on
plot(t_model,x_model(:,2),'--')
subplot(311)
xlabel('$t$')
ylabel(['$' models.names.x{1} '$'])

subplot(312)
plot(t_model,x_tape(:,1))
hold on
plot(t_model,x_model(:,1),'--')
xlabel('$t$')
ylabel(['$' models.names.x{2} '$'])

subplot(313)
plot(t_model,z_tape(:,3))
hold on
plot(t_model,z_model(:,3),'--')
xlabel('$t$')
ylabel(['$' models.names.z{3} '$'])

legend({'Engine','LPV model'})

fig_error = figure(2*i);
fig_error.Name = ['Error histogram'];
subplot(131)
for j = 1:nx
    histogram(e_x(:,j),'Normalization','probability');
    hold on
end
xlabel('$error~\%$')
ylabel('$\%~of~time$')
title('State differential')
legend(models.names.x,'Interpreter','tex')

subplot(132)
for j = 1:ny
    histogram(e_y(:,j),'Normalization','probability');
    hold on
end
xlabel('$error~\%$')
ylabel('$\%~of~time$')
title('State algebraic')
legend(models.names.y,'Interpreter','tex')

subplot(133)
for j = 1:nz
    histogram(e_z(:,j),'Normalization','probability');
    hold on
end
xlabel('$error~\%$')
ylabel('$\%~of~time$')
title('Output')
legend(models.names.z,'Interpreter','tex')




