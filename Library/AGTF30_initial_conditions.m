% AGTF30_INITIAL_CONDITIONS
% ========================================================================
% 
%
% Pablo Baldivieso Monasterios
% The University of Sheffield
% 14/02/2024

function MWS = AGTF30_initial_conditions(MWS)


MWS.In.AltIC = MWS.In.Alt.Data(1,1);
MWS.In.MNIC  = MWS.In.MN.Data(1,1);
MWS.In.dTambIC = MWS.In.dT.Data(1,1);


[cond,MWS.In.AltIC,MWS.In.MNIC,MWS.In.dTambIC] = AGTF30.inEnvelope(MWS);
Alt = MWS.In.AltIC;
MN  = MWS.In.MNIC;
dT = MWS.In.dTambIC;


% Get IC points
load('SS_Mtrx.mat');

% FOR A GIVEN INPUT, WE COMPUTE THE REQUIRED VALUE

Wfreq = MWS.In.Wf.Data(1,1);


% create Alt and MN axis
AltVec = SS.Alt;
NcVec  = SS.Nc;
MNVec  = SS.MN;


% Verify table inputs are within bounds
Alt    = CheckMnMx(Alt   ,min(AltVec)  , max(AltVec),'Alt');
MN     = CheckMnMx(MN    ,min(MNVec)   , max(MNVec) ,'MN');
        

Nc_guess = 1500;
Wfguess = interpn(AltVec,NcVec,MNVec,SS.Ind{10}{1},Alt,Nc_guess,MN);
Wf_Er = Wfreq - Wfguess;
Nc_guessNew = 1501;
iter = 1;
Er_T = 0.001;

while ((abs(Wf_Er) > Er_T) && (iter < 200))
    % gather Error and iteration info
    Wf_ErOld = Wf_Er;
    Nc_guessOld = Nc_guess;
    Nc_guess = Nc_guessNew;
    % Perform table lookup with new values
    Wfguess = interpn(AltVec,NcVec,MNVec,SS.Ind{10}{1},Alt,Nc_guess,MN);
    Wf_Er = Wfreq - Wfguess;
    % if error is large use secant algorithm to guess a new Nc
    % value
    if abs(Wf_Er) > Er_T
        if abs(Wf_Er - Wf_ErOld) < 0.000001
            Wf_Er = Wf_Er*1.01;
        end
        Nc_guessNew = Nc_guess - Wf_Er * (Nc_guess - Nc_guessOld)/(Wf_Er - Wf_ErOld);
    end

    iter = iter + 1;
end

if iter >= 200
    fprintf(['Can not automatically determine corrected speed that matches with specified fuel flow, Wfref\n']);
end
Nc = Nc_guess;

% check speed meets vector requirements
Nc     = CheckMnMx(Nc,min(NcVec), max(NcVec) ,'Nc');
% check speed meets limiter requirements
NcMin = interp2(MNVec,AltVec,SS.Ncmin,MN,Alt);
NcMax = interp2(MNVec,AltVec,SS.Ncmax,MN,Alt);
Nc     = CheckMnMx(Nc,NcMin, NcMax,'Nc');


% Set final and limited corrected speed value
MWS.In.ssNc1 = Nc;

% Set VAFN throat area
MWS.In.VAFNIC = interp2(MWS.Cntrl.VAFN_MN,MWS.Cntrl.VAFN_Nc1,MWS.Cntrl.VAFN_sch, MN, Nc);

% Set VBV open position
MWS.In.VBVIC = interp2(MWS.Cntrl.VBV_MN,MWS.Cntrl.VBV_Nc1,MWS.Cntrl.VBV_sch, MN, Nc);

% Read starting points for IC generation
% dT adjustment
dT_adj = 1 + dT/1000; % rule of thumb adjustment

MWS.In.ICcom(1) = interpn(AltVec,NcVec,MNVec,SS.Ind{1}{1},Alt,Nc,MN)/ dT_adj;       % Flow (pps)
MWS.In.ICcom(2) = interpn(AltVec,NcVec,MNVec,SS.Ind{2}{1},Alt,Nc,MN)/ (1 + (dT_adj-1) * MN);               % HPT Pressure Ratio
MWS.In.ICcom(3) = interpn(AltVec,NcVec,MNVec,SS.Ind{3}{1},Alt,Nc,MN)/ (1 + (dT_adj-1) * MN);               % HPC Rline
MWS.In.ICcom(4) = interpn(AltVec,NcVec,MNVec,SS.Ind{4}{1},Alt,Nc,MN);               % Fan_Rline
MWS.In.ICcom(5) = interpn(AltVec,NcVec,MNVec,SS.Ind{5}{1},Alt,Nc,MN);               % Branch Pressure Ratio
MWS.In.ICcom(6) = interpn(AltVec,NcVec,MNVec,SS.Ind{6}{1},Alt,Nc,MN);               % LPC_Rline
MWS.In.ICcom(7) = interpn(AltVec,NcVec,MNVec,SS.Ind{7}{1},Alt,Nc,MN);               % LPT Pressure Ratio

MWS.In.ICss(1)  = interpn(AltVec,NcVec,MNVec,SS.Ind{8}{1},Alt,Nc,MN)* dT_adj;       % Low Pressure Shaft Speed (rpm)
MWS.In.ICss(2)  = interpn(AltVec,NcVec,MNVec,SS.Ind{9}{1},Alt,Nc,MN)* dT_adj;       % High Pressure Shaft Speed (rpm)
MWS.In.ICss(3)  = interpn(AltVec,NcVec,MNVec,SS.Ind{10}{1},Alt,Nc,MN)* dT_adj;      % Fuel Flow (pps)
MWS.In.ICss(4)  = interpn(AltVec,NcVec,MNVec,SS.Ind{11}{1},Alt,Nc,MN);              % VAFN area, VAFN (in2)
MWS.In.ICss(5)  = interpn(AltVec,NcVec,MNVec,SS.Ind{12}{1},Alt,Nc,MN) + 1;          % VBV position plus 1 (1-closed, 2-open)

MWS.In.N2IC = MWS.In.ICss(1);
MWS.In.N3IC = MWS.In.ICss(2);
MWS.In.WfIC = MWS.In.ICss(3);

MWS.In.VAFNManEn = 0;
MWS.In.VBVManEn  = 0;


MWS = SetSolverVars(MWS);


assignin('base', 'MWS',MWS);
fprintf(['Generating ICs with steady-state solver...\n']);
sim('AGTF30SysSS');
fprintf(['IC generation complete\n']);

SSconv = out_SS.converged.Data(end);
        
if SSconv == 1
    try
        MWS.In.ICcom = out_SS.independents.Data(end,1:7);
        MWS.In.ICss  = out_SS.independents.Data(end,8:end);
        MWS.In.N2IC = MWS.In.ICss(1);
        MWS.In.N3IC = MWS.In.ICss(2);
        MWS.In.WfIC = MWS.In.ICss(3);
    catch
        fprintf(['Sensor initial condtions were not fully generated\n']);
    end
end

try
    % Set Sensor ICs
    MWS.In.Sensor.ICN1 = MWS.In.N2IC/3.1;
    MWS.In.Sensor.ICN2 = MWS.In.N2IC;
    MWS.In.Sensor.ICN3 = MWS.In.N3IC;

    MWS.In.Sensor.ICPa  = out_SS.Pa.Data(end);
    MWS.In.Sensor.ICP2  = out_SS.S2.Pt.Data(end);
    MWS.In.Sensor.ICP21 = out_SS.S21.Pt.Data(end);
    MWS.In.Sensor.ICP25 = out_SS.S25.Pt.Data(end);
    MWS.In.Sensor.ICPs3 = out_SS.S36.Pt.Data(end);
    MWS.In.Sensor.ICP5  = out_SS.S5.Pt.Data(end);

    MWS.In.Sensor.ICT2  = out_SS.S2.Tt.Data(end);
    MWS.In.Sensor.ICT25 = out_SS.S25.Tt.Data(end);
    MWS.In.Sensor.ICT3 = out_SS.S36.Tt.Data(end);
    MWS.In.Sensor.ICT45  = out_SS.S45.Tt.Data(end);
catch
    fprintf(['Sensor initial condtions were not fully generated\n']);
end


disp('** AGTF30 ready to execute **')




end



function Vf = CheckMnMx(Vo, Mn, Mx,VarName)
    % Check input value by min and max values producing an error message and
    % limiting the value if it fails the check
    %-------------------------------------------------------------------------
    if Vo < Mn
        Vf = Mn;
        fprintf(['Initial condition ',VarName,' too low, updating to %d \n'],Mn)
    elseif Vo > Mx
        Vf = Mx;
        fprintf(['Initial condition ',VarName,' too high, updating to %d \n'],Mx)
    else
        Vf = Vo;
    end

end

function MWS = SetSolverVars(MWS)
    % Solver Input Settings (for SS and dyn)
    %---------------------------------
    %Solver inputs
    ICnumCom = length(MWS.In.ICcom);
    ICnumTot = ICnumCom +  length(MWS.In.ICss);
    %set solver perturbation sizes
    MWS.In.JPerSS = [0.001*ones(ICnumTot,1) -0.001*ones(ICnumTot,1)]';
    MWS.In.JPerDyn = [0.001*ones(ICnumCom,1) -0.001*ones(ICnumCom,1)]';
    %set Condition limit
    MWS.In.Lim = 1e-5;
    %set max number of solver iterations
    %for run (steady-state) or per time step (dyn)
    MWS.In.Max_IterDyn = numel(MWS.In.JPerDyn) + 50;
    MWS.In.Max_IterSS = 2*(numel(MWS.In.JPerSS) + 40);
    % set number of consecuative time steps for determining a dynamic run is
    % failing to converge.  Once limit is reached the dynamic simulation will
    % stop
    MWS.In.DynRunNonConv = 400;
    %set number of solver attempts before Jacobian re-calc.
    MWS.In.NRADyn = 25;
    MWS.In.NRASS = MWS.In.Max_IterSS/2;
    %set max % step change for solver
    MWS.In.dX = 1;
    %set Perturbation values for X in the linear solver
    FracStep = [0.01; 0.005;-0.005;-0.01];
    MWS.In.XPerLin = [MWS.In.N2IC*FracStep,MWS.In.N3IC*FracStep];
    % Set Perturabation steps for u in the linear solver
    MWS.In.UPerLin = MWS.In.WfIC * FracStep;
end
