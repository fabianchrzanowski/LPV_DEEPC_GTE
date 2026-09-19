% ========================================================================
% LPV-DeePC Disturbance Rejection plotter
%
% loads pre-computed disturbance rejection results
% and generates figures.
%
% ========================================================================

% ===== USER CONFIGURATION =====
% DeePC disturbance result (leave '' to auto-scan *disturbance* files)
DEEPC_DIST_FILE = '';  % e.g. 'Results/results_lpv_deepc_disturbance4.mat'
% MPC disturbance result (if available)
MPC_DIST_FILE = '';    % leave '' if not yet generated
% ==============================

clearvars -except DEEPC_DIST_FILE MPC_DIST_FILE; close all; clc

try proj = currentProject; catch; proj = openProject('GTE.prj'); end
cd(proj.RootFolder)

plot_options;

% -------------------------------------------------------------------------
%  Color palette (consistent with dissertation conventions)
% -------------------------------------------------------------------------
c_deepc   = [0.00 0.45 0.74];   % Blue  — DeePC main tracking color (disturbed output)
c_mpc     = [0.85 0.33 0.10];   % Orange — MPC main tracking color (disturbed output)
c_ref     = [0.20 0.20 0.20];   % Dark charcoal — reference
c_belief  = [0.85 0.45 0.45];   % Faint Red — controller internal belief (undisturbed true output)
c_applied = [0.80 0.10 0.10];   % Red dashed — actual applied Wf (with d_u)
c_dist_y  = [0.85 0.10 0.10];   % Red — output disturbance
c_dist_u  = [0.00 0.45 0.74];   % Blue — input disturbance
c_shade_y = [1.00 0.82 0.82];   % Light red patch — output disturbance window
c_shade_u = [0.82 0.88 1.00];   % Light blue patch — input disturbance window

% =========================================================================
%  FILE SELECTION: DeePC disturbance results
% =========================================================================
results_dir = fullfile(proj.RootFolder, 'Results');

if isempty(DEEPC_DIST_FILE)
    % Auto-scan for *disturbance* files
    candidates = dir(fullfile(results_dir, '*disturbance*.mat'));
    if isempty(candidates)
        error('No disturbance result files found in %s.\nRun AGTF30_lpv_deepc_disturbance_test.m first.', results_dir);
    end
    fprintf('\nFound %d disturbance result file(s):\n', numel(candidates));
    for ii = 1:numel(candidates)
        d_info = dir(fullfile(results_dir, candidates(ii).name));
        fprintf('  [%d] %s  (%.1f kB, %s)\n', ii, candidates(ii).name, ...
            d_info.bytes/1024, datestr(d_info.datenum, 'yyyy-mm-dd HH:MM'));
    end
    choice = input('\nSelect file number for DeePC disturbance results: ');
    if ~isnumeric(choice) || choice < 1 || choice > numel(candidates)
        error('Invalid selection.');
    end
    deepc_file = fullfile(results_dir, candidates(choice).name);
else
    % User-specified path (absolute or relative to project root)
    if isAbsolutePath_local(DEEPC_DIST_FILE)
        deepc_file = DEEPC_DIST_FILE;
    else
        deepc_file = fullfile(proj.RootFolder, DEEPC_DIST_FILE);
    end
end

if ~isfile(deepc_file)
    error('DeePC disturbance results not found at:\n  %s', deepc_file);
end

fprintf('\nLoading DeePC disturbance results from:\n  %s\n', deepc_file);
D = load(deepc_file);

% If MPC file is not specified, try to guess it from the DeePC filename
if isempty(MPC_DIST_FILE)
    [~, deepc_name, ext] = fileparts(deepc_file);
    mpc_name = strrep(deepc_name, 'deepc', 'mpc');
    mpc_file_guess = fullfile(results_dir, [mpc_name, ext]);
    if isfile(mpc_file_guess)
        MPC_DIST_FILE = mpc_file_guess;
        fprintf('Auto-detected matching MPC disturbance file: %s\n', [mpc_name, ext]);
    end
end

% =========================================================================
%  FILE SELECTION: MPC disturbance results (optional)
% =========================================================================
have_mpc = false;
if ~isempty(MPC_DIST_FILE)
    if isAbsolutePath_local(MPC_DIST_FILE)
        mpc_file = MPC_DIST_FILE;
    else
        mpc_file = fullfile(proj.RootFolder, MPC_DIST_FILE);
    end
    if isfile(mpc_file)
        fprintf('Loading MPC disturbance results from:\n  %s\n', mpc_file);
        M = load(mpc_file);
        have_mpc = true;
    else
        warning('MPC disturbance file not found at:\n  %s\nProceeding with DeePC-only plots.', mpc_file);
    end
else
    fprintf('MPC_DIST_FILE not set — DeePC-only plots will be generated.\n');
end

% =========================================================================
%  PARSE DeePC DATA
% =========================================================================
deepc_time    = D.time(:);
deepc_y_hist  = D.y_hist(:);
deepc_y_true  = D.y_true(:);
deepc_u_hist  = D.u_hist(:);
deepc_u_app   = D.u_applied(:);
deepc_dist_y  = D.dist_y_hist(:);
deepc_dist_u  = D.dist_u_hist(:);
deepc_ref     = D.ref_plot(:);
deepc_stimes  = D.solve_times(:);

% track_name for axis label
track_name = '';
if isfield(D, 'track_name'); track_name = D.track_name; end
ytarget_label = make_ytarget_label(track_name);

% Controller name / scenario label
ctrl_label_deepc = 'LPV-DeePC';
if isfield(D, 'controller_name')
    ctrl_label_deepc = D.controller_name;
end

% -------------------------------------------------------------------------
%  Determine which portion of the time vector to plot.
%  DeePC files include warmup (T_ini steps at the start where solve_times
%  entries are zeros).  We plot the FULL time vector as stored — the
%  dist_y/dist_u vectors tell us exactly when disturbances were active.
% -------------------------------------------------------------------------
N_deepc = length(deepc_time);
N_st    = length(deepc_stimes);

% The control period starts where solve_times has non-zero length
% (after the warmup).  The offset between N_deepc and N_st is T_ini.
T_ini_deepc = N_deepc - N_st;   % 0 if DeePC file has no warmup padding

% Derive disturbance windows from the history vectors (robust — works even
% if DIST_OUTPUT/DIST_INPUT were not saved in older files).
[dy_wins_deepc] = derive_dist_windows(deepc_dist_y, deepc_time);
[du_wins_deepc] = derive_dist_windows(deepc_dist_u, deepc_time);
n_dy_events = size(dy_wins_deepc, 1);
n_du_events = size(du_wins_deepc, 1);
n_dist_events = n_dy_events + n_du_events;

% RMSE (over full control window, exclude warmup)
idx_ctrl_deepc = (T_ini_deepc + 1):N_deepc;
err_true_deepc = deepc_y_true(idx_ctrl_deepc) - deepc_ref(idx_ctrl_deepc);
rmse_deepc = sqrt(mean(err_true_deepc.^2, 'omitnan'));

% =========================================================================
%  PARSE MPC DATA (if available)
% =========================================================================
if have_mpc
    mpc_time   = M.time(:);
    mpc_y_hist = M.y_hist(:);
    mpc_y_true = M.y_true(:);
    mpc_u_hist = M.u_hist(:);
    mpc_u_app  = M.u_applied(:);
    mpc_dist_y = M.dist_y_hist(:);
    mpc_dist_u = M.dist_u_hist(:);
    mpc_ref    = M.ref_plot(:);
    mpc_stimes = M.solve_times(:);

    ctrl_label_mpc = 'LPV-MPC';
    if isfield(M, 'controller_name'); ctrl_label_mpc = M.controller_name; end

    [dy_wins_mpc] = derive_dist_windows(mpc_dist_y, mpc_time);
    [du_wins_mpc] = derive_dist_windows(mpc_dist_u, mpc_time);

    err_true_mpc = mpc_y_true - mpc_ref;
    rmse_mpc = sqrt(mean(err_true_mpc.^2, 'omitnan'));
end

% =========================================================================
%  CONSOLE SUMMARY
% =========================================================================
fprintf('\n============================================================\n');
fprintf('  DISTURBANCE REJECTION SUMMARY\n');
fprintf('============================================================\n');
fprintf('Controller: %s\n', ctrl_label_deepc);
fprintf('File:       %s\n', deepc_file);
fprintf('Sim steps:  %d  (warmup = %d, control = %d)\n', N_deepc, T_ini_deepc, N_st);
fprintf('Output disturbance events: %d\n', n_dy_events);
fprintf('Input  disturbance events: %d\n', n_du_events);
fprintf('RMSE (true output):        %.1f RPM\n', rmse_deepc);
fprintf('Solve time — mean: %.2f ms | max: %.2f ms\n', ...
    mean(deepc_stimes,'omitnan')*1000, max(deepc_stimes,[],'omitnan')*1000);

if have_mpc
    fprintf('\n--- MPC ---\n');
    fprintf('Controller: %s\n', ctrl_label_mpc);
    fprintf('File:       %s\n', mpc_file);
    fprintf('RMSE (true output):   %.1f RPM\n', rmse_mpc);
    fprintf('Solve time — mean: %.2f ms | max: %.2f ms\n', ...
        mean(mpc_stimes,'omitnan')*1000, max(mpc_stimes,[],'omitnan')*1000);
    if rmse_deepc < rmse_mpc
        pct = (rmse_mpc - rmse_deepc) / rmse_mpc * 100;
        fprintf('DeePC improvement over MPC: %.1f%%\n', pct);
    else
        pct = (rmse_deepc - rmse_mpc) / rmse_deepc * 100;
        fprintf('MPC improvement over DeePC: %.1f%%\n', pct);
    end
end
fprintf('============================================================\n\n');

% =========================================================================
%  FIGURE 1: DeePC Disturbance Rejection (Tracking & Fuel Flow)
% =========================================================================
fig1 = figure('Name', 'DeePC: Tracking & Fuel Flow');
tiledlayout(2,1, 'TileSpacing','compact', 'Padding','compact');

% ---- Panel 1: NH tracking ----
ax1 = nexttile;
hold(ax1, 'on');

plot(ax1, deepc_time, deepc_ref,    '--', 'Color', c_ref,   'LineWidth', 1.5, ...
    'DisplayName', 'Reference');
plot(ax1, deepc_time, deepc_y_hist, '-',  'Color', c_deepc, 'LineWidth', 1.8, ...
    'DisplayName', 'Measured output (disturbed)');

% Shade disturbance windows
shade_windows(ax1, dy_wins_deepc, c_shade_y);
shade_windows(ax1, du_wins_deepc, c_shade_u);

% Add invisible patches for legend entries
if ~isempty(dy_wins_deepc)
    patch(ax1, NaN, NaN, c_shade_y, 'FaceAlpha', 0.5, 'EdgeColor', 'none', ...
        'DisplayName', 'Output dist. window');
end
if ~isempty(du_wins_deepc)
    patch(ax1, NaN, NaN, c_shade_u, 'FaceAlpha', 0.5, 'EdgeColor', 'none', ...
        'DisplayName', 'Input dist. window');
end

ylabel(ax1, ytarget_label, 'Interpreter', 'latex', 'FontSize', 11);
title(ax1, sprintf('\\textbf{$N_H$ Tracking} (RMSE$_{\\mathrm{true}}$ = %.1f rpm)', rmse_deepc), ...
    'Interpreter', 'latex', 'FontSize', 11);
legend(ax1, 'Location', 'best', 'Interpreter', 'latex', 'FontSize', 9);
grid(ax1, 'on');
set(ax1, 'XTickLabel', {}, 'TickLabelInterpreter', 'latex');

% ---- Panel 2: Wf command ----
ax2 = nexttile;
hold(ax2, 'on');

stairs(ax2, deepc_time, deepc_u_hist, '-',  'Color', c_deepc,   'LineWidth', 1.8, ...
    'DisplayName', 'Controller output ($u$)');
stairs(ax2, deepc_time, deepc_u_app,  '--', 'Color', c_applied, 'LineWidth', 1.2, ...
    'DisplayName', 'Applied ($u + d_u$)');

% Shade disturbance windows
shade_windows(ax2, dy_wins_deepc, c_shade_y);
shade_windows(ax2, du_wins_deepc, c_shade_u);

ylabel(ax2, '$W_f$ [pps]', 'Interpreter', 'latex', 'FontSize', 11);
title(ax2, '\textbf{Fuel Flow Command}', 'Interpreter', 'latex', 'FontSize', 11);
legend(ax2, 'Location', 'best', 'Interpreter', 'latex', 'FontSize', 9);
grid(ax2, 'on');
xlabel(ax2, 'Time [s]', 'Interpreter', 'latex', 'FontSize', 11);
sgtitle(fig1, '\textbf{LPV-DeePC: Tracking \& Fuel Flow}', 'Interpreter', 'latex', 'FontSize', 13);
% linkaxes([ax1 ax2], 'x');

% =========================================================================
%  FIGURE 2: DeePC Applied Disturbances
% =========================================================================
fig1_dist = figure('Name', 'DeePC: Applied Disturbances');
ax3 = axes(fig1_dist);
hold(ax3, 'on');

yyaxis(ax3, 'left');
stairs(ax3, deepc_time, deepc_dist_y, '-', 'Color', c_dist_y, 'LineWidth', 1.8);
ylabel(ax3, 'Output dist. [rpm]', 'Interpreter', 'latex', 'FontSize', 11);
ax3.YAxis(1).Color = c_dist_y;

yyaxis(ax3, 'right');
stairs(ax3, deepc_time, deepc_dist_u, '-', 'Color', c_dist_u, 'LineWidth', 1.8);
ylabel(ax3, 'Input dist. [pps]', 'Interpreter', 'latex', 'FontSize', 11);
ax3.YAxis(2).Color = c_dist_u;

title(ax3, '\textbf{Applied Disturbances}', 'Interpreter', 'latex', 'FontSize', 11);
grid(ax3, 'on');
xlabel(ax3, 'Time [s]', 'Interpreter', 'latex', 'FontSize', 11);

% =========================================================================
%  FIGURE 3: DeePC True Tracking Error
% =========================================================================
fig1_err = figure('Name', 'DeePC: True Tracking Error');
ax4 = axes(fig1_err);
hold(ax4, 'on');

plot(ax4, deepc_time, deepc_y_true - deepc_ref, '-', ...
    'Color', c_deepc, 'LineWidth', 1.8, 'DisplayName', 'True error');

% Shade disturbance windows
shade_windows(ax4, dy_wins_deepc, c_shade_y);
shade_windows(ax4, du_wins_deepc, c_shade_u);
yline(ax4, 0, 'k--', 'LineWidth', 1.2, 'HandleVisibility', 'off');

xlabel(ax4, 'Time [s]', 'Interpreter', 'latex', 'FontSize', 11);
ylabel(ax4, 'Error [rpm]', 'Interpreter', 'latex', 'FontSize', 11);
title(ax4, sprintf('\\textbf{True Tracking Error} (RMSE = %.1f rpm)', rmse_deepc), ...
    'Interpreter', 'latex', 'FontSize', 11);
legend(ax4, 'Location', 'best', 'Interpreter', 'latex', 'FontSize', 9);
grid(ax4, 'on');
set(ax4, 'TickLabelInterpreter', 'latex');

% =========================================================================
%  FIGURE 4: MPC Disturbance Rejection (Tracking & Fuel Flow)
% =========================================================================
if have_mpc
    fig2 = figure('Name', 'MPC: Tracking & Fuel Flow');
    tiledlayout(2,1, 'TileSpacing','compact', 'Padding','compact');

    % ---- Panel 1: NH tracking ----
    bx1 = nexttile;
    hold(bx1, 'on');

    plot(bx1, mpc_time, mpc_ref,    '--', 'Color', c_ref,   'LineWidth', 1.5, ...
        'DisplayName', 'Reference');
    plot(bx1, mpc_time, mpc_y_hist, '-',  'Color', c_mpc,   'LineWidth', 1.8, ...
        'DisplayName', 'Measured output (disturbed)');

    % Shade disturbance windows
    shade_windows(bx1, dy_wins_mpc, c_shade_y);
    shade_windows(bx1, du_wins_mpc, c_shade_u);

    if ~isempty(dy_wins_mpc)
        patch(bx1, NaN, NaN, c_shade_y, 'FaceAlpha', 0.5, 'EdgeColor', 'none', ...
            'DisplayName', 'Output dist. window');
    end
    if ~isempty(du_wins_mpc)
        patch(bx1, NaN, NaN, c_shade_u, 'FaceAlpha', 0.5, 'EdgeColor', 'none', ...
            'DisplayName', 'Input dist. window');
    end

    ylabel(bx1, ytarget_label, 'Interpreter', 'latex', 'FontSize', 11);
    title(bx1, sprintf('\\textbf{$N_H$ Tracking} (RMSE$_{\\mathrm{true}}$ = %.1f rpm)', rmse_mpc), ...
        'Interpreter', 'latex', 'FontSize', 11);
    legend(bx1, 'Location', 'best', 'Interpreter', 'latex', 'FontSize', 9);
    grid(bx1, 'on');
    set(bx1, 'XTickLabel', {}, 'TickLabelInterpreter', 'latex');

    % ---- Panel 2: Wf command ----
    bx2 = nexttile;
    hold(bx2, 'on');

    stairs(bx2, mpc_time, mpc_u_hist, '-',  'Color', c_mpc,     'LineWidth', 1.8, ...
        'DisplayName', 'Controller output ($u$)');
    stairs(bx2, mpc_time, mpc_u_app,  '--', 'Color', c_applied, 'LineWidth', 1.2, ...
        'DisplayName', 'Applied ($u + d_u$)');

    % Shade disturbance windows
    shade_windows(bx2, dy_wins_mpc, c_shade_y);
    shade_windows(bx2, du_wins_mpc, c_shade_u);

    ylabel(bx2, '$W_f$ [pps]', 'Interpreter', 'latex', 'FontSize', 11);
    title(bx2, '\textbf{Fuel Flow Command}', 'Interpreter', 'latex', 'FontSize', 11);
    legend(bx2, 'Location', 'best', 'Interpreter', 'latex', 'FontSize', 9);
    grid(bx2, 'on');
    xlabel(bx2, 'Time [s]', 'Interpreter', 'latex', 'FontSize', 11);
    sgtitle(fig2, '\textbf{LPV-MPC: Tracking \& Fuel Flow}', 'Interpreter', 'latex', 'FontSize', 13);
    % linkaxes([bx1 bx2], 'x');

    % =========================================================================
    %  FIGURE 5: MPC Applied Disturbances
    % =========================================================================
    fig2_dist = figure('Name', 'MPC: Applied Disturbances');
    bx3 = axes(fig2_dist);
    hold(bx3, 'on');

    yyaxis(bx3, 'left');
    stairs(bx3, mpc_time, mpc_dist_y, '-', 'Color', c_dist_y, 'LineWidth', 1.8);
    ylabel(bx3, 'Output dist. [rpm]', 'Interpreter', 'latex', 'FontSize', 11);
    bx3.YAxis(1).Color = c_dist_y;

    yyaxis(bx3, 'right');
    stairs(bx3, mpc_time, mpc_dist_u, '-', 'Color', c_dist_u, 'LineWidth', 1.8);
    ylabel(bx3, 'Input dist. [pps]', 'Interpreter', 'latex', 'FontSize', 11);
    bx3.YAxis(2).Color = c_dist_u;

    title(bx3, '\textbf{Applied Disturbances}', 'Interpreter', 'latex', 'FontSize', 11);
    grid(bx3, 'on');
    xlabel(bx3, 'Time [s]', 'Interpreter', 'latex', 'FontSize', 11);

    % =========================================================================
    %  FIGURE 6: MPC True Tracking Error
    % =========================================================================
    fig2_err = figure('Name', 'MPC: True Tracking Error');
    bx4 = axes(fig2_err);
    hold(bx4, 'on');

    plot(bx4, mpc_time, mpc_y_true - mpc_ref, '-', ...
        'Color', c_mpc, 'LineWidth', 1.8, 'DisplayName', 'True error');

    % Shade disturbance windows
    shade_windows(bx4, dy_wins_mpc, c_shade_y);
    shade_windows(bx4, du_wins_mpc, c_shade_u);
    yline(bx4, 0, 'k--', 'LineWidth', 1.2, 'HandleVisibility', 'off');

    xlabel(bx4, 'Time [s]', 'Interpreter', 'latex', 'FontSize', 11);
    ylabel(bx4, 'Error [rpm]', 'Interpreter', 'latex', 'FontSize', 11);
    title(bx4, sprintf('\\textbf{True Tracking Error} (RMSE = %.1f rpm)', rmse_mpc), ...
        'Interpreter', 'latex', 'FontSize', 11);
    legend(bx4, 'Location', 'best', 'Interpreter', 'latex', 'FontSize', 9);
    grid(bx4, 'on');
    set(bx4, 'TickLabelInterpreter', 'latex');

end

% =========================================================================
%  FIGURE 3: Side-by-side RMSE comparison (only if both files exist)
% =========================================================================
if have_mpc
    fig3 = figure('Name', 'RMSE Comparison');
    ax_bar = axes(fig3);
    hold(ax_bar, 'on');

    bar_vals   = [rmse_deepc, rmse_mpc];
    bar_colors = [c_deepc; c_mpc];
    bar_labels = {'LPV-DeePC', 'LPV-MPC'};

    for ii = 1:2
        b = bar(ax_bar, ii, bar_vals(ii), 0.55, 'FaceColor', bar_colors(ii,:), ...
            'EdgeColor', 'none');
        % Value annotation above each bar
        text(ax_bar, ii, bar_vals(ii) + max(bar_vals)*0.02, ...
            sprintf('%.1f', bar_vals(ii)), ...
            'HorizontalAlignment', 'center', 'VerticalAlignment', 'bottom', ...
            'Interpreter', 'latex', 'FontSize', 11, 'FontWeight', 'bold');
    end

    % Improvement annotation
    if rmse_deepc < rmse_mpc
        pct_imp = (rmse_mpc - rmse_deepc) / rmse_mpc * 100;
        ann_str = sprintf('DeePC: $-%.1f\\%%$ vs MPC', pct_imp);
        ann_color = c_deepc;
    else
        pct_imp = (rmse_deepc - rmse_mpc) / rmse_deepc * 100;
        ann_str = sprintf('MPC: $-%.1f\\%%$ vs DeePC', pct_imp);
        ann_color = c_mpc;
    end

    text(ax_bar, 1.5, max(bar_vals)*0.88, ann_str, ...
        'HorizontalAlignment', 'center', 'Interpreter', 'latex', ...
        'FontSize', 10, 'Color', ann_color, 'FontWeight', 'bold');

    set(ax_bar, 'XTick', 1:2, 'XTickLabel', bar_labels, ...
        'TickLabelInterpreter', 'latex', 'FontSize', 11);
    xlim(ax_bar, [0.3 2.7]);
    ylim(ax_bar, [0, max(bar_vals) * 1.20]);
    ylabel(ax_bar, 'RMSE [rpm]', 'Interpreter', 'latex', 'FontSize', 12);
    title(ax_bar, '\textbf{Disturbance Rejection: RMSE Comparison}', ...
        'Interpreter', 'latex', 'FontSize', 12);
    grid(ax_bar, 'on');
    ax_bar.GridAlpha = 0.3;
    set(ax_bar, 'TickLabelInterpreter', 'latex');
end

fprintf('Done. Generated %d figure(s).\n', 1 + have_mpc + have_mpc);

% =========================================================================
%  LOCAL HELPER FUNCTIONS
% =========================================================================

function lbl = make_ytarget_label(track_name)
    if contains(lower(track_name), 'n_h') || isempty(track_name)
        lbl = 'Core Speed ($N_H$) [rpm]';
    elseif contains(lower(track_name), 'n_l')
        lbl = 'Fan Speed ($N_L$) [rpm]';
    else
        lbl = strrep(track_name, '_', '\_');
    end
end

function wins = derive_dist_windows(dist_vec, time_vec)
%DERIVE_DIST_WINDOWS  Find contiguous non-zero segments in dist_vec.
%   Returns an Mx2 matrix where each row is [t_start, t_end] in seconds.
    dist_vec = dist_vec(:);
    time_vec = time_vec(:);
    active = (dist_vec ~= 0);
    wins = [];
    if ~any(active); return; end
    % find rising/falling edges
    d_active = diff([0; active; 0]);
    starts = find(d_active ==  1);
    ends   = find(d_active == -1) - 1;
    for ii = 1:length(starts)
        wins(end+1, :) = [time_vec(starts(ii)), time_vec(ends(ii))]; %#ok<AGROW>
    end
end

function shade_windows(ax, wins, rgb)
%SHADE_WINDOWS  Draw semi-transparent shaded patches over disturbance windows.
    if isempty(wins); return; end
    yl = ylim(ax);
    for ii = 1:size(wins, 1)
        t1 = wins(ii, 1);
        t2 = wins(ii, 2);
        h = patch(ax, [t1 t2 t2 t1], [yl(1) yl(1) yl(2) yl(2)], rgb, ...
            'FaceAlpha', 0.30, 'EdgeColor', 'none', 'HandleVisibility', 'off');
        uistack(h, 'bottom');
    end
    ylim(ax, yl);
end

function tf = isAbsolutePath_local(p)
%ISABSOLUTEPATH_LOCAL  Returns true if p starts with / or a drive letter.
    tf = ~isempty(p) && (p(1) == '/' || p(1) == '\' || ...
         (length(p) >= 2 && p(2) == ':'));
end
