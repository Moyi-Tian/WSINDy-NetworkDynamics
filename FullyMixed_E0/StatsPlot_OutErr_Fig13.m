%% Appendix boxplots: output error cross-section at one noise level
% Plot E and P only, no U. Show y-axis label only for E.

clear
close all

%% cd to current directory
scriptFullName = matlab.desktop.editor.getActiveFilename();
scriptDir = fileparts(scriptFullName);
if ~isempty(scriptDir)
    cd(scriptDir);
end

%% Load Data
cur_dir = pwd;

beta = 0.5;
theta = 0.4;
eta = 0.2;
gamma_i = 0.1;
gamma_p = 0.3;

time_range = [0 50];
t_num = 1000;

start_val = 0.01;
end_val = 0.3;
num_trajs_upto = 10;

noise_ratios = linspace(0,0.1,21);
numNoise = length(noise_ratios);

rept = 100;
save_path = fullfile(cur_dir, 'WSINDy_data_3StateInput');

fname = sprintf(['WSINDy Results for Fully-Mixed Network Dynamics E0 Batch Data ', ...
    'numTrajUpTo=%u, E0_=%.5g-%.5g, beta=%.3g, theta=%.3g, eta=%.3g, ', ...
    'gamma_i=%.3g, gamma_p=%.3g, tmax=%u, t_grid_num=%u, NumNoise=%u, ', ...
    'repeat=%u_noise_%.5g-%.5g_MSTLS.mat'], ...
    num_trajs_upto, start_val, end_val, beta, theta, eta, gamma_i, gamma_p, ...
    time_range(end), t_num, numNoise, rept, noise_ratios(1), noise_ratios(end));

load(fullfile(save_path, fname))

%% User choices
plot_trajNum_to = 5;
selected_noise_level = 0.05;    % change this as needed
use_log_scale = true;           % set false if you want linear scale

eq_letters = {'U','E','P'};
nEq = 3;

% Plot E and P only
eq_to_plot = [2, 3];

%% Find nearest available noise level
[~, noise_idx] = min(abs(noise_ratios - selected_noise_level));
alpha_val = noise_ratios(noise_idx);

fprintf('Requested noise level: %.5g\n', selected_noise_level);
fprintf('Using nearest available noise level: %.5g\n', alpha_val);

alpha_str = sprintf('%g', alpha_val);
alpha_str = strrep(alpha_str, '.', 'p');

%% Output folder
plot_path = fullfile(cur_dir, 'plots_variability_appendix');
if ~exist(plot_path, 'dir')
    mkdir(plot_path);
end

%% Extract output error values
% OutErrVals has dimensions:
% equation/state x number of trajectories x repeat
OutErrVals = NaN(nEq, plot_trajNum_to, rept);

for eq = eq_to_plot
    for nt = 1:plot_trajNum_to
        for r = 1:rept
            curVal = trajectory_errors_cell{nt, noise_idx, r};
            if ~isempty(curVal) && numel(curVal) >= eq
                OutErrVals(eq, nt, r) = curVal(eq);
            end
        end
    end
end

%% Determine common y-axis limits using E and P only
allVals = OutErrVals(eq_to_plot, :, :);
allVals = allVals(:);
allVals = allVals(isfinite(allVals));

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

%% Generate one boxplot per state, E and P only
SummaryRows = {};

for eq = eq_to_plot

    letter = eq_letters{eq};
    vals = squeeze(OutErrVals(eq, :, :));
    % vals dimensions: number of trajectories x repeat

    dataVec = [];
    groupVec = [];

    for nt = 1:plot_trajNum_to
        curVals = vals(nt, :).';
        dataVec = [dataVec; curVals];
        groupVec = [groupVec; nt * ones(length(curVals), 1)];
    end

    keep = isfinite(dataVec);
    dataVec_plot = dataVec(keep);
    groupVec_plot = groupVec(keep);

    figure(300 + eq); clf;
    set(gcf, 'Units', 'inches', 'Position', [3 3 6.2 5.4]);

    boxplot(dataVec_plot, groupVec_plot, ...
        'Labels', arrayfun(@num2str, 1:plot_trajNum_to, 'UniformOutput', false), ...
        'Symbol', '.', ...
        'Widths', 0.55);

    hold on

    % Median line
    medVals = median(vals, 2, 'omitnan');
    plot(1:plot_trajNum_to, medVals, 'ko-', ...
        'LineWidth', 1.5, ...
        'MarkerFaceColor', 'k', ...
        'MarkerSize', 5);

    %% Axis formatting
    ax = gca;
    box on

    % Use enough left/bottom margin so tick labels are not clipped
    ax.Position = [0.15 0.15 0.80 0.80];

    set(ax, 'FontSize', 18);

    xlabel('Number of Trajectories', ...
        'FontSize', 22, ...
        'Interpreter', 'latex');

    % Show y-axis label only for E, not P
    if eq == 2
        ylabel('Output Error', ...
            'FontSize', 22, ...
            'Interpreter', 'latex');
    else
        ylabel('');
    end

    ylim(common_ylim);
    xlim([0.5, plot_trajNum_to + 0.5]);

    if use_log_scale_final
        set(ax, 'YScale', 'log');

        % Force clean powers-of-ten ticks
        ytick_powers = floor(log10(common_ylim(1))) : ceil(log10(common_ylim(2)));
        ytick_vals = 10.^ytick_powers;

        keep_ticks = ytick_vals >= common_ylim(1) & ytick_vals <= common_ylim(2);
        ytick_vals = ytick_vals(keep_ticks);
        ytick_powers = ytick_powers(keep_ticks);

        ax.YTick = ytick_vals;

        % Use tex for robust superscripts in exported MATLAB PDFs
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
    filename = sprintf('%s_OutputErr_CrossSection_Boxplot_noise_%s.pdf', ...
        letter, alpha_str);

    exportgraphics(gcf, fullfile(plot_path, filename), ...
        'ContentType', 'vector', ...
        'BackgroundColor', 'none', ...
        'Resolution', 300);

    %% Summary table information
    for nt = 1:plot_trajNum_to
        curVals = vals(nt, :);
        curVals = curVals(isfinite(curVals));

        medVal = median(curVals, 'omitnan');
        q25 = prctile(curVals, 25);
        q75 = prctile(curVals, 75);
        iqrVal = q75 - q25;

        SummaryRows(end+1, :) = {letter, alpha_val, nt, medVal, q25, q75, iqrVal};
    end
end

%% Save summary table
SummaryTable = cell2table(SummaryRows, ...
    'VariableNames', {'State', 'NoiseLevel', 'NumTrajectories', ...
    'Median', 'Q25', 'Q75', 'IQR'});

disp(SummaryTable)

writetable(SummaryTable, fullfile(plot_path, ...
    sprintf('OutputErr_CrossSection_SummaryTable_EP_only_noise_%s.csv', alpha_str)));