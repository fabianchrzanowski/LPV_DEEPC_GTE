%% ========================================================================
% setup_simulation_params.m
% 
% This is a configuration file which affects the specification of all the
% controllers.
%
% Set trajectory, Initial steps and prediction horizon
% Choose from:
% - QUICK_TEST -> short 140 steps step reference trajectory
% - MID_TEST -> around 600-700 steps trajectory with increase and decrease
% of core speed 
% - LONG_TEST -> around 2000 steps trajectory with big transients, as well
% as gradual increase and decrease
% ========================================================================

% Choose the trajectory profile
if ~exist('SCENARIO', 'var')
    SCENARIO = 'MID_TEST'; % Options: 'QUICK_TEST', 'MID_TEST', 'LONG_TEST'
end

%% 3. Timing Parameters
Ts = 0.015;          % Sample time [s]
T_ini = 20;          % Initial steps, runs also for lpv-mpc
N_pred = 30;         % Prediction horizon

if strcmp(SCENARIO, 'QUICK_TEST')
    steps_per_seg = 70;
    n_segments = 2;
elseif strcmp(SCENARIO, 'MID_TEST')
    % Mid length trajectory: ~600 steps
    dx_Wf = [0.5*ones(1,5), ...                     % idle hold
             linspace(0.5, 1.2, 8), ...              % ramp to cruise
             1.2*ones(1,5), ...                      % cruise hold
             1.8, 1.8, 1.8, ...                      % step to high
             linspace(1.8, 0.7, 6), ...              % ramp down
             0.7, 0.7, 0.7, ...                      % low hold
             1.5, 0.9, 1.5, 0.9, ...                 % rapid toggles
             linspace(0.9, 0.5, 4), ...              % ramp back to idle
             0.5*ones(1,3)];                          % idle settle
    
    steps_per_seg = 15;  % 0.225s per segment
    n_segments = numel(dx_Wf);
elseif strcmp(SCENARIO, 'LONG_TEST')
    % Long trajectory: ~2000 steps
    seg1 = 0.5*ones(1,20);                          % idle hold
    seg2 = linspace(0.5, 1.9, 30);                   % ramp idle->max
    seg3 = 1.9*ones(1,15);                           % hold at max
    seg4 = [1.9 0.8 1.6 0.5 1.4 1.0 1.8 0.6 1.5];  % rapid steps
    seg5 = 1.2 + 0.6*sin(2*pi*(0:39)/20);           % sinusoidal
    seg6 = [0.5*ones(1,10), linspace(0.5,1.7,15), ...% takeoff-cruise-descent
            1.7*ones(1,20), linspace(1.7,1.2,10), ...
            1.2*ones(1,15), linspace(1.2,0.5,10), ...
            0.5*ones(1,10)];
    dx_Wf = [seg1, seg2, seg3, seg4, seg5, seg6];
    
    steps_per_seg = 10; % 0.15s per segment / 0.015s Ts = 10 steps
    n_segments = numel(dx_Wf);
else
    steps_per_seg = 50; % 200 * 0.015 = 3 seconds
    n_segments = 3;
end

% simulation steps
N_sim = n_segments * steps_per_seg;
N_total = T_ini + N_sim;
L = T_ini + N_pred;  % DeePC trajectory length

%% 3. Load Identification Data 
try; proj = currentProject; catch; proj = openProject('GTE.prj'); end
models = load(fullfile(proj.RootFolder,'Params','tables','matrices_DE_wf.mat'));
wf_data = load(fullfile(proj.RootFolder,'Params','test_files_fuel_only','test_identification','wf_values.mat'));
offset = models.offset;
wf_values = wf_data.Wf_values(:)';

idx_NH = find(contains(lower(string(models.sys_eng.OutputName)), "n_h"), 1);
if isempty(idx_NH); idx_NH = 2; end

%% 5. Define Transient Scenario
wf_idle = 0.5;

if strcmp(SCENARIO, 'QUICK_TEST')
    % Gentle step from Idle to a slightly higher power setting
    wf_max = 0.8;
    wf_cruise = 0.8;
else
    % Idle -> Max -> Cruise
    wf_max = 1.9;
    wf_cruise = 1.2;
end

[~, k_idle]   = min(abs(wf_values - wf_idle));
[~, k_max]    = min(abs(wf_values - wf_max));
[~, k_cruise] = min(abs(wf_values - wf_cruise));

u0 = offset.u(:, k_idle);
z0 = offset.z(:, k_idle);

z_idle   = offset.z(idx_NH, k_idle);
z_max    = offset.z(idx_NH, k_max);
z_cruise = offset.z(idx_NH, k_cruise);

% Build step reference
ref = zeros(N_total + N_pred, 1);
ref(1:T_ini) = z_idle;                            % Warmup

% quick test is a quick step
if strcmp(SCENARIO, 'QUICK_TEST')
    idx1 = T_ini + 1; idx2 = T_ini + steps_per_seg;   % Idle
    ref(idx1:idx2) = z_idle;
    idx1 = idx2 + 1; idx2 = idx2 + steps_per_seg;     % Gentle Step up
    ref(idx1:end) = z_max;

% mid test and long test 
elseif strcmp(SCENARIO, 'MID_TEST') || strcmp(SCENARIO, 'LONG_TEST')
    idx1 = T_ini + 1;
    for i = 1:numel(dx_Wf)
        idx2 = idx1 + steps_per_seg - 1;
        % map wf to the spool speed using the data from identification
        z_target = interp1(offset.u(1,:)', offset.z(idx_NH,:)', dx_Wf(i), 'pchip');
        ref(idx1:idx2) = z_target;
        idx1 = idx2 + 1;
    end
    ref(idx1:end) = ref(idx1-1); % Pad the horizon
else
    idx1 = T_ini + 1; idx2 = T_ini + steps_per_seg;   % Idle
    ref(idx1:idx2) = z_idle;
    idx1 = idx2 + 1; idx2 = idx2 + steps_per_seg;     % Max
    ref(idx1:idx2) = z_max;
    idx1 = idx2 + 1; idx2 = idx2 + steps_per_seg;     % Cruise
    ref(idx1:end) = z_cruise;
end

fprintf('--- Master Parameters Loaded ---\n');
fprintf('N_pred: %d\n', N_pred);
fprintf('Total Sim Steps: %d\n', N_total);
fprintf('Trajectory: %s\n', SCENARIO);
fprintf('--------------------------------\n');
