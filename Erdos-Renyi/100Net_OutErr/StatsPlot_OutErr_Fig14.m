%% Appendix boxplots: stochastic output error cross-section at one k

% Before running this script, run Compute_Store_AvgOutErr.m first to get
% summary matrices.

% Loads preprocessed raw output-error values from
% TrajectoryErrSummary_100Networks_FromOneK_*.mat.
% Plots U, E, and P separately. Only U has the y-axis label.

clear; close all; tic

%% Set Directory
scriptFullName = matlab.desktop.editor.getActiveFilename();
scriptDir = fileparts(scriptFullName);
if ~isempty(scriptDir), cd(scriptDir); end
cur_dir = pwd;

%% User choices
selected_k = 5;
plot_trajNum_to = 5;
use_log_scale = true;

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
    'TrajectoryErr_values_all', ...
    'state_labels', ...
    'summary_meta');

TrajectoryErr_values_all = S.TrajectoryErr_values_all;
state_labels = S.state_labels;
meta = S.summary_meta;

%% Find selected k
k_values_all = meta.k_values(:);
[~, k_pos] = min(abs(k_values_all - selected_k));
k_actual = k_values_all(k_pos);

fprintf('Requested k = %g\n', selected_k);
fprintf('Using k = %g\n', k_actual);

if isfield(meta, 'num_trajs_upto')
    num_trajs_saved = meta.num_trajs_upto;
else
    num_trajs_saved = size(TrajectoryErr_values_all{1}, 2);
end

plot_trajNum_to = min(plot_trajNum_to, num_trajs_saved);
traj_idx = 1:plot_trajNum_to;

%% Output folder
plot_path = fullfile(cur_dir, 'plots_variability_appendix');
if ~exist(plot_path, 'dir')
    mkdir(plot_path);
end

%% Determine common y-axis limits across U/E/P
allVals = [];

for s = 1:3
    for jj = 1:plot_trajNum_to
        vals = TrajectoryErr_values_all{s}{k_pos, jj};
        vals = vals(:);
        vals = vals(isfinite(vals));
        allVals = [allVals; vals]; %#ok<AGROW>
    end
end

if isempty(allVals)
    error('No finite output-error values found for k = %g.', k_actual);
end

if use_log_scale && all(allVals > 0)
    y_min = min(allVals);
    y_max = max(allVals);

    y_min_pow = floor(log10(y_min));
    y_max_pow = ceil(log10(y_max));

    common_ylim = [10^y_min_pow, 10^y_max_pow];
    use_log_scale_final = true;
else
    y_min = min(allVals);
    y_max = max(allVals);

    if y_min == y_max
        common_ylim = [y_min - 0.1*abs(y_min+eps), y_max + 0.1*abs(y_max+eps)];
    else
        pad = 0.08 * (y_max - y_min);
        common_ylim = [max(0, y_min - pad), y_max + pad];
    end

    use_log_scale_final = false;
end

%% Generate boxplots
SummaryRows = {};

for s = 1:3

    dataVec = [];
    groupVec = [];
    medVals = nan(plot_trajNum_to, 1);

    for jj = 1:plot_trajNum_to
        vals = TrajectoryErr_values_all{s}{k_pos, jj};
        vals = vals(:);
        vals = vals(isfinite(vals));

        dataVec = [dataVec; vals]; %#ok<AGROW>
        groupVec = [groupVec; jj * ones(numel(vals), 1)]; %#ok<AGROW>

        medVals(jj) = median(vals, 'omitnan');

        if isempty(vals)
            q25 = NaN;
            q75 = NaN;
            iqrVal = NaN;
            nVals = 0;
        else
            q25 = prctile(vals, 25);
            q75 = prctile(vals, 75);
            iqrVal = q75 - q25;
            nVals = numel(vals);
        end

        SummaryRows(end+1, :) = {state_labels{s}, k_actual, jj, ...
            medVals(jj), q25, q75, iqrVal, nVals}; %#ok<SAGROW>
    end

    keep = isfinite(dataVec);
    dataVec_plot = dataVec(keep);
    groupVec_plot = groupVec(keep);

    figure(300 + s); clf;
    set(gcf, 'Units', 'inches', 'Position', [3 3 6.2 5.4]);

    boxplot(dataVec_plot, groupVec_plot, ...
        'Labels', arrayfun(@num2str, traj_idx, 'UniformOutput', false), ...
        'Symbol', '.', ...
        'Widths', 0.55);

    hold on

    plot(1:plot_trajNum_to, medVals, 'ko-', ...
        'LineWidth', 1.5, ...
        'MarkerFaceColor', 'k', ...
        'MarkerSize', 5);

    ax = gca;
    box on
    ax.Position = [0.15 0.15 0.80 0.80];

    set(ax, 'FontSize', 22);

    xlabel('Number of Trajectories', ...
        'FontSize', 28);

    if s == 1
        ylabel('Output Error', ...
            'FontSize', 28);
    else
        ylabel('');
    end

    ylim(common_ylim);
    xlim([0.5, plot_trajNum_to + 0.5]);

    if use_log_scale_final
        set(ax, 'YScale', 'log');

        ytick_powers = floor(log10(common_ylim(1))) : ceil(log10(common_ylim(2)));
        ytick_vals = 10.^ytick_powers;

        keep_ticks = ytick_vals >= common_ylim(1) & ytick_vals <= common_ylim(2);
        ytick_vals = ytick_vals(keep_ticks);
        ytick_powers = ytick_powers(keep_ticks);

        ax.YTick = ytick_vals;

        ax.TickLabelInterpreter = 'tex';
        ax.YTickLabel = arrayfun(@(p) sprintf('10^{%d}', p), ...
            ytick_powers, 'UniformOutput', false);
    else
        ax.TickLabelInterpreter = 'latex';
        ax.YTickMode = 'auto';
        ax.YTickLabelMode = 'auto';
    end

    drawnow;

    %% Save figure
    label = state_labels{s};
    label_clean = regexprep(label, '[^a-zA-Z0-9]', '');

    filename = sprintf('ER_OutputErr_%s_CrossSection_k=%u.pdf', ...
        label_clean, k_actual);

    exportgraphics(gcf, fullfile(plot_path, filename), ...
        'ContentType', 'vector', ...
        'BackgroundColor', 'none', ...
        'Resolution', 300);
end

%% Save summary table
SummaryTable = cell2table(SummaryRows, ...
    'VariableNames', {'State', 'k', 'NumTrajectories', ...
    'Median', 'Q25', 'Q75', 'IQR', 'NumFiniteValues'});

disp(SummaryTable)

writetable(SummaryTable, fullfile(plot_path, ...
    sprintf('ER_OutputErr_CrossSection_SummaryTable_k=%u.csv', k_actual)));

toc