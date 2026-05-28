% Before running this code, run Compute_Save_TrajInfo_forPlotting.m first
% This code loads payload and plots for selected num_traj and show_id
tic
clear
close all

%% cd to current directory
scriptFullName = matlab.desktop.editor.getActiveFilename();
scriptDir = fileparts(scriptFullName);
if ~isempty(scriptDir), cd(scriptDir); end
cur_dir = pwd;

%% -------- Choose payload file --------
payload_dir = fullfile(cur_dir,'TrajectoryComparison_payloads');
mat_list = dir(fullfile(payload_dir,'TrajComparePayload_*.mat'));
if isempty(mat_list)
    error('No payload files found in:\n%s', payload_dir);
end
[~, idxNewest] = max([mat_list.datenum]);
payload_path = fullfile(payload_dir, mat_list(idxNewest).name);

fprintf('Loading payload:\n%s\n', payload_path);
S = load(payload_path,'payload');        % load into struct (no workspace overwrite)
payload = S.payload;

%% -------- User choices --------
save_plots = 1;
num_traj = 10;   % which "num trajectories" model to use (1..num_trajs_upto)
show_id  = 6;    % which trajectory within that bundle (1..num_traj)

%% Validate
if num_traj < 1 || num_traj > payload.meta.num_trajs_upto
    error('num_traj must be in [1, %d].', payload.meta.num_trajs_upto);
end
if show_id < 1 || show_id > num_traj
    error('show_id must be in [1, num_traj].');
end

t_grid = payload.t_grid;
time_grid = payload.time_grid;
Ssel = payload.S;

J = payload.meta.J;

%% Map to the selected IC index in X6_cell and mean-field row
idx_sel = J(Ssel{num_traj});
id = idx_sel(show_id);

%% Pull true mean + envelope
data = payload.X6_cell{id};
U_mean = data(:,1);
E_mean = data(:,2);
P_mean = data(:,5);

min_data = payload.X6_mins{id};
max_data = payload.X6_maxes{id};

U_min = min_data(:,1)';  U_max = max_data(:,1)';
E_min = min_data(:,2)';  E_max = max_data(:,2)';
P_min = min_data(:,5)';  P_max = max_data(:,5)';

%% Pull inferred trajectory for this num_traj + show_id
U_inf = payload.U_inferred_trajs{num_traj}(show_id,:);
E_inf = payload.E_inferred_trajs{num_traj}(show_id,:);
P_inf = payload.P_inferred_trajs{num_traj}(show_id,:);

%% Pull mean-field solution row corresponding to IC index "id"
U_mf = payload.Usol_mf(id,:);
E_mf = payload.Esol_mf(id,:);
P_mf = payload.Psol_mf(id,:);

%% Output directory
out_dir = fullfile(cur_dir, sprintf(['TrajectoryComparisonPlots_fromPayload_' ...
    'N=%u_k=%u_numTraj=%u_show=%u_noise=%.5g_NetID=%u'], ...
    payload.meta.N, payload.meta.k, num_traj, show_id, payload.meta.noise_ratio, payload.meta.network_id));
if ~exist(out_dir,'dir'), mkdir(out_dir); end

%% ===== Settings =====
% Uses Okabe–Ito palette:
%   True      = Black
%   Inferred  = Vermillion
%   Mean-field= Blue

cb_true = [0 0 0];               % black
cb_inf  = [0.835 0.369 0.000];   % vermillion (#D55E00)
cb_mf   = [0.000 0.447 0.698];   % blue       (#0072B2)

% Line widths
lw_true = 4.0;     % truth dominant
lw_inf  = 2.5;
lw_mf   = 2.1;

% Line styles
ls_true = '-';
ls_inf  = '--';
ls_mf   = '-.';

% Markers (only for inferred & mean-field by default)
mk_inf = 'o';
mk_mf  = 's';

% Marker sizes (can differ per line)
ms_inf = 11;
ms_mf  = 8;

% Marker edge "length/thickness" in MATLAB is controlled by LineWidth
% (so mf marker edges will be slightly thicker than inferred).
lw_inf_marker = lw_inf;
lw_mf_marker  = lw_mf;

% Marker face: turn OFF (hollow markers)
mfc_off = 'none';

% How many markers to show (avoid clutter)
mk_every_true = max(1, floor(numel(t_grid)/15));
mk_every_inf  = max(1, floor(numel(time_grid)/12));
mk_every_mf   = max(1, floor(numel(time_grid)/12));

% Envelope fill (subtle, doesn’t fight parula-like backgrounds)
env_alpha = 0.18;

%% ======================= U plot =======================
fU = figure('Visible','off','Units','inches','Position',[1 1 8 6]);
axU = axes('Parent',fU); hold(axU,'on');

% Envelope first
fill([t_grid, fliplr(t_grid)], [U_max, fliplr(U_min)], ...
     cb_true, 'FaceAlpha', env_alpha, 'EdgeColor','none');

% True: thick solid black, optional sparse markers (comment out Marker lines if prefer none)
p1 = plot(t_grid, U_mean, ...
    'LineStyle', ls_true, 'Color', cb_true, 'LineWidth', lw_true, ...
    'Marker', 'none');  % set to 'o' if you want markers on truth too

% Inferred: dashed vermillion + hollow circles
p2 = plot(time_grid, U_inf, ...
    'LineStyle', ls_inf, 'Color', cb_inf, 'LineWidth', lw_inf_marker, ...
    'Marker', mk_inf, 'MarkerSize', ms_inf, ...
    'MarkerFaceColor', mfc_off, 'MarkerEdgeColor', cb_inf, ...
    'MarkerIndices', 1:mk_every_inf:numel(time_grid));

% Mean-field: dash-dot blue + hollow squares
p3 = plot(time_grid, U_mf, ...
    'LineStyle', ls_mf, 'Color', cb_mf, 'LineWidth', lw_mf_marker, ...
    'Marker', mk_mf, 'MarkerSize', ms_mf, ...
    'MarkerFaceColor', mfc_off, 'MarkerEdgeColor', cb_mf, ...
    'MarkerIndices', 1:mk_every_mf:numel(time_grid));
set(axU,'FontSize',25);
legend([p1 p2 p3], 'True', 'Inferred', 'Mean-field', 'Location','best','FontSize',26);
xlabel('Time','FontSize',32);
ylabel('$U$ Proportion','Interpreter','latex','FontSize',32);
grid on; box on;

fnameU = sprintf('U_numTraj=%u_show=%u.pdf', num_traj, show_id);
try
    exportgraphics(axU, fullfile(out_dir,fnameU), 'ContentType','vector','BackgroundColor','none');
catch
    set(fU,'PaperPositionMode','auto');
    print(fU, fullfile(out_dir,fnameU), '-dpdf', '-painters');
end
close(fU);

%% ======================= E plot =======================
fE = figure('Visible','off','Units','inches','Position',[1 1 8 6]);
axE = axes('Parent',fE); hold(axE,'on');

fill([t_grid, fliplr(t_grid)], [E_max, fliplr(E_min)], ...
     cb_true, 'FaceAlpha', env_alpha, 'EdgeColor','none');

p4 = plot(t_grid, E_mean, ...
    'LineStyle', ls_true, 'Color', cb_true, 'LineWidth', lw_true, ...
    'Marker', 'none');

p5 = plot(time_grid, E_inf, ...
    'LineStyle', ls_inf, 'Color', cb_inf, 'LineWidth', lw_inf_marker, ...
    'Marker', mk_inf, 'MarkerSize', ms_inf, ...
    'MarkerFaceColor', mfc_off, 'MarkerEdgeColor', cb_inf, ...
    'MarkerIndices', 1:mk_every_inf:numel(time_grid));

p6 = plot(time_grid, E_mf, ...
    'LineStyle', ls_mf, 'Color', cb_mf, 'LineWidth', lw_mf_marker, ...
    'Marker', mk_mf, 'MarkerSize', ms_mf, ...
    'MarkerFaceColor', mfc_off, 'MarkerEdgeColor', cb_mf, ...
    'MarkerIndices', 1:mk_every_mf:numel(time_grid));

set(axE,'FontSize',25);
xlabel('Time','FontSize',30);
ylabel('$E$ Proportion','Interpreter','latex','FontSize',30);
grid on; box on;

fnameE = sprintf('E_numTraj=%u_show=%u.pdf', num_traj, show_id);
try
    exportgraphics(axE, fullfile(out_dir,fnameE), 'ContentType','vector','BackgroundColor','none');
catch
    set(fE,'PaperPositionMode','auto');
    print(fE, fullfile(out_dir,fnameE), '-dpdf', '-painters');
end
close(fE);

%% ======================= P plot =======================
fP = figure('Visible','off','Units','inches','Position',[1 1 8 6]);
axP = axes('Parent',fP); hold(axP,'on');

fill([t_grid, fliplr(t_grid)], [P_max, fliplr(P_min)], ...
     cb_true, 'FaceAlpha', env_alpha, 'EdgeColor','none');

p7 = plot(t_grid, P_mean, ...
    'LineStyle', ls_true, 'Color', cb_true, 'LineWidth', lw_true, ...
    'Marker', 'none');

p8 = plot(time_grid, P_inf, ...
    'LineStyle', ls_inf, 'Color', cb_inf, 'LineWidth', lw_inf_marker, ...
    'Marker', mk_inf, 'MarkerSize', ms_inf, ...
    'MarkerFaceColor', mfc_off, 'MarkerEdgeColor', cb_inf, ...
    'MarkerIndices', 1:mk_every_inf:numel(time_grid));

p9 = plot(time_grid, P_mf, ...
    'LineStyle', ls_mf, 'Color', cb_mf, 'LineWidth', lw_mf_marker, ...
    'Marker', mk_mf, 'MarkerSize', ms_mf, ...
    'MarkerFaceColor', mfc_off, 'MarkerEdgeColor', cb_mf, ...
    'MarkerIndices', 1:mk_every_mf:numel(time_grid));

set(axP,'FontSize',25);
xlabel('Time','FontSize',30);
ylabel('$P$ Proportion','Interpreter','latex','FontSize',30);
grid on; box on;

fnameP = sprintf('P_numTraj=%u_show=%u.pdf', num_traj, show_id);
try
    exportgraphics(axP, fullfile(out_dir,fnameP), 'ContentType','vector','BackgroundColor','none');
catch
    set(fP,'PaperPositionMode','auto');
    print(fP, fullfile(out_dir,fnameP), '-dpdf', '-painters');
end
close(fP);


toc
