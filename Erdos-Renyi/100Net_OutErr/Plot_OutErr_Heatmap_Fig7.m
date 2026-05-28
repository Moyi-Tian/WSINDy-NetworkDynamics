%% Plot output-error heatmaps from preprocessed 100-network summary matrices
% Before running this script, run Compute_Store_AvgOutErr.m first to get
% summary matrices.

% This script does not recompute averages. It only loads saved summary
% matrices and plots heatmaps for U, E, and P.
%
% cap_vals is state-specific:
%   cap_vals = [cap_U, cap_E, cap_P]
%   Use NaN for a state if no color saturation/capping should be applied.

clear; close all; tic

%% Set Directory
scriptFullName = matlab.desktop.editor.getActiveFilename();
scriptDir = fileparts(scriptFullName);
if ~isempty(scriptDir), cd(scriptDir); end
cur_dir = pwd;

%% Load preprocessed summary file
in_dir = fullfile(cur_dir, 'TrajectoryErr_matrices_100Networks');

mat_list = dir(fullfile(in_dir, 'TrajectoryErrSummary_100Networks_FromOneK_*.mat'));
if isempty(mat_list)
    mat_list = dir(fullfile(in_dir, 'TrajectoryErrSummary_100Networks_*.mat'));
end

if isempty(mat_list)
    error('No saved trajectory-error summary file found in:\n%s', in_dir);
end

[~, idxNewest] = max([mat_list.datenum]);
in_fullpath = fullfile(mat_list(idxNewest).folder, mat_list(idxNewest).name);

fprintf('Loading preprocessed summary file:\n%s\n', in_fullpath);

S = load(in_fullpath, ...
    'TrajectoryErr_mean_all', ...
    'TrajectoryErr_median_all', ...
    'TrajectoryErr_numfinite_all', ...
    'TrajectoryErr_numslots_all', ...
    'state_labels', ...
    'summary_meta');

state_labels = S.state_labels;
meta = S.summary_meta;

%% Choose summary type to plot
% Options: 'mean' or 'median'
summary_to_plot = 'mean';

switch lower(summary_to_plot)
    case 'mean'
        TrajectoryErr_plot_all = S.TrajectoryErr_mean_all;
        summary_label = 'Mean';
    case 'median'
        TrajectoryErr_plot_all = S.TrajectoryErr_median_all;
        summary_label = 'Median';
    otherwise
        error('summary_to_plot must be either ''mean'' or ''median''.');
end

%% User-specified plot ranges
k_min = 4;
k_max = 10;

traj_min = 1;
traj_max = 5;

% State-specific cap values for log10 color saturation.
% Order is U, E, P.
% Example:
%   cap_vals = [NaN, 2, 1];  % no cap for U, cap E at 10^2, cap P at 10^1
%   cap_vals = [2, 2, 2];    % same cap for all states
%   cap_vals = [NaN, NaN, NaN]; % no cap for any state
cap_vals = [0.25, 3, 0.2];

if isscalar(cap_vals)
    cap_vals = repmat(cap_vals, 1, 3);
end

if numel(cap_vals) ~= 3
    error('cap_vals must be either a scalar or a 1-by-3 vector for U, E, P.');
end

%% Validate and slice indices
k_values_all = meta.k_values(:);

k_mask = (k_values_all >= k_min) & (k_values_all <= k_max);
k_values_plot = k_values_all(k_mask);

if isempty(k_values_plot)
    error('Requested k range [%d,%d] has no overlap with saved k_values.', k_min, k_max);
end

k_idx = find(k_mask);

if isfield(meta, 'num_trajs_upto')
    num_trajs_saved = meta.num_trajs_upto;
else
    num_trajs_saved = size(TrajectoryErr_plot_all{1}, 2);
end

traj_min = max(1, traj_min);
traj_max = min(num_trajs_saved, traj_max);

traj_idx = traj_min:traj_max;
x_values = traj_idx;

%% Plot output directory
cap_tag_parts = cell(1,3);
for s = 1:3
    if isnan(cap_vals(s))
        cap_tag_parts{s} = sprintf('%sNoCap', state_labels{s});
    else
        tmp = sprintf('%sCap%.3g', state_labels{s}, cap_vals(s));
        tmp = strrep(tmp, '.', 'p');
        cap_tag_parts{s} = tmp;
    end
end
cap_tag = strjoin(cap_tag_parts, '_');

plot_path = fullfile(cur_dir, ...
    sprintf('Heatmaps_100Networks_%s_k=%d-%d_traj=%d-%d_%s', ...
    lower(summary_to_plot), k_min, k_max, traj_min, traj_max, cap_tag));

if ~exist(plot_path, 'dir')
    mkdir(plot_path);
end

%% Make heatmaps for U/E/P
for s = 1:3

    cap_val_s = cap_vals(s);

    figure('Units', 'inches', 'Position', [6 2 6.2 8]);

    M_full = TrajectoryErr_plot_all{s};
    M = M_full(k_idx, traj_idx);

    Z = log10(M);

    nonfinite_mask = ~isfinite(Z);
    finite_vals = Z(~nonfinite_mask);

    if isempty(finite_vals)
        finite_vals = 0;
    end

    Z_show = Z;
    Z_show(nonfinite_mask) = max(finite_vals);

    %% Apply optional state-specific cap
    if isnan(cap_val_s)
        overflow_mask = false(size(Z_show));
    else
        overflow_mask = Z_show > cap_val_s;
        Z_show(overflow_mask) = cap_val_s;
    end

    imagesc(x_values, k_values_plot, Z_show);
    set(gca, 'YDir', 'normal');

    colormap("parula");
    box on; grid off;

    %% Vertical separator between 2 and 3 trajectories
    hold on

    yl = ylim;
    delta = 0.2;

    h = plot([2.5 2.5], ...
             [yl(1)-delta, yl(2)+delta], ...
             '-', ...
             'Color', [0.70 0.15 0.15], ...
             'LineWidth', 4);

    h.Clipping = 'off';

    hold off

    %% Coordinates for annotations
    [Xgrid, Ygrid] = meshgrid(x_values, k_values_plot);

    %% Annotate nonfinite entries
    [rN, cN] = find(nonfinite_mask);
    for ii = 1:numel(rN)
        text(Xgrid(rN(ii), cN(ii)), Ygrid(rN(ii), cN(ii)), 'NaN', ...
            'Color', 'r', ...
            'FontWeight', 'bold', ...
            'FontSize', 14, ...
            'HorizontalAlignment', 'center', ...
            'VerticalAlignment', 'middle');
    end

    %% Annotate overflow values only when cap is active for this state
    if ~isnan(cap_val_s)
        valid_overflow = overflow_mask & ~nonfinite_mask;
        [rO, cO] = find(valid_overflow);

        for ii = 1:numel(rO)
            raw = Z(rO(ii), cO(ii));
            txt_val = sprintf('%.1e', 10^raw);

            text(Xgrid(rO(ii), cO(ii)), Ygrid(rO(ii), cO(ii)), txt_val, ...
                'Color', 'r', ...
                'FontWeight', 'bold', ...
                'FontSize', 14, ...
                'HorizontalAlignment', 'center', ...
                'VerticalAlignment', 'middle');
        end
    end

    %% Colorbar
    c = colorbar;
    c.FontSize = 18;

    % Only the first heatmap has colorbar label
    if s == 1
        c.Label.String = sprintf('%s Output Error (log)', summary_label);
        c.Label.FontSize = 22;
    else
        c.Label.String = '';
    end

    %% Axes
    ax = gca;
    set(ax, 'FontSize', 20);

    xlabel('Number of Trajectories', 'FontSize', 25);

    % Only the first heatmap has y-axis label
    if s == 1
        ylabel('$k$', 'Interpreter', 'latex', 'FontSize', 28);
    else
        ylabel('');
    end

    %% Color limits
    cmin = min(finite_vals);

    if isnan(cap_val_s)
        cmax = max(finite_vals);
    else
        cmax = min(max(finite_vals), cap_val_s);
    end

    if cmin == cmax
        cmin = cmin - 0.5;
        cmax = cmax + 0.5;
    end

    caxis([cmin cmax]);

    %% Save
    if isnan(cap_val_s)
        state_cap_tag = 'noCap';
    else
        state_cap_tag = sprintf('cap%.3g', cap_val_s);
        state_cap_tag = strrep(state_cap_tag, '.', 'p');
    end

    fname = sprintf('Heatmap_%s_%sOutputErr_k=%d-%d_traj=%d-%d_100Networks_%s.pdf', ...
        state_labels{s}, summary_label, k_min, k_max, traj_min, traj_max, state_cap_tag);

    exportgraphics(gcf, fullfile(plot_path, fname), ...
        'ContentType', 'vector', ...
        'BackgroundColor', 'none');
end

toc