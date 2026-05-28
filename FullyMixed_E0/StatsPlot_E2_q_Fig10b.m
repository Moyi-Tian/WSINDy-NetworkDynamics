%% Appendix boxplot: q cross-section at one noise level
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
max_traj_plot = 5;
selected_noise_level = 0.02;   % change this if needed

metric_label = '$q$';
metric_file_tag = 'q';

%% Find nearest available noise level
[~, noise_idx] = min(abs(noise_ratios - selected_noise_level));
alpha_val = noise_ratios(noise_idx);

fprintf('Requested noise level: %.5g\n', selected_noise_level);
fprintf('Using nearest available noise level: %.5g\n', alpha_val);

%% Compute q values at the selected noise level
% q_vals dimensions: max_traj_plot x rept
q_vals = NaN(max_traj_plot, rept);

den = norm(weights_true);

for nt = 1:max_traj_plot
    W = cell2mat(Weights_learned_cell(nt, noise_idx, 1:rept));
    W = reshape(W, [], rept);

    if den > 0
        diffs = W - weights_true;
        q_r = vecnorm(diffs, 2, 1) / den;
    else
        q_r = vecnorm(W, 2, 1);
    end

    q_vals(nt, :) = q_r;
end

vals = q_vals;

%% Prepare data for boxplot
dataVec = [];
groupVec = [];

for nt = 1:max_traj_plot
    curVals = vals(nt, :).';
    dataVec = [dataVec; curVals];
    groupVec = [groupVec; nt * ones(length(curVals), 1)];
end

% Remove NaN and Inf values from boxplot
keep = isfinite(dataVec);
dataVec_plot = dataVec(keep);
groupVec_plot = groupVec(keep);

%% Output folder
plot_path = fullfile(cur_dir, 'plots_variability_appendix');
if ~exist(plot_path, 'dir')
    mkdir(plot_path);
end

%% Plot
figure;
set(gcf, 'Units', 'inches', 'Position', [3 3 7 7]);

boxplot(dataVec_plot, groupVec_plot, ...
    'Labels', arrayfun(@num2str, 1:max_traj_plot, 'UniformOutput', false), ...
    'Symbol', '.', ...
    'Widths', 0.55);

hold on

% Median line
medVals = median(vals, 2, 'omitnan');
plot(1:max_traj_plot, medVals, 'ko-', ...
    'LineWidth', 1.5, ...
    'MarkerFaceColor', 'k', ...
    'MarkerSize', 5);

set(gca, 'FontSize', 18);
box on

% q is often skewed. Use log scale if all finite plotted values are positive.
if all(dataVec_plot > 0)
    set(gca, 'YScale', 'log');

    ax = gca;
    ax.TickLabelInterpreter = 'latex';

    % Choose log ticks manually from the data range
    y_min = min(dataVec_plot);
    y_max = max(dataVec_plot);
    ytick_powers = floor(log10(y_min)) : ceil(log10(y_max));

    ax.YTick = 10.^ytick_powers;
    ax.YTickLabel = arrayfun(@(p) sprintf('$10^{%d}$', p), ...
        ytick_powers, 'UniformOutput', false);
end

set(gca, 'FontSize', 18);
xlabel('Number of Trajectories', 'FontSize', 20);
ylabel(metric_label, 'FontSize', 24, 'Interpreter', 'latex');
drawnow;

%% Save figure
alpha_str = sprintf('%g', alpha_val);
alpha_str = strrep(alpha_str, '.', 'p');

figure_name = sprintf('%s_CrossSection_Boxplot_noise_%s.pdf', ...
    metric_file_tag, alpha_str);

exportgraphics(gcf, fullfile(plot_path, figure_name), ...
    'ContentType', 'vector', ...
    'BackgroundColor', 'none', ...
    'Resolution', 300);

%% Summary table: medians, IQRs, and changes
SummaryRows = {};
ChangeRows = {};

for nt = 1:max_traj_plot
    curVals = vals(nt, :);
    curVals = curVals(isfinite(curVals));

    medVal = median(curVals, 'omitnan');
    q25 = prctile(curVals, 25);
    q75 = prctile(curVals, 75);
    iqrVal = q75 - q25;

    SummaryRows(end+1, :) = {alpha_val, nt, medVal, q25, q75, iqrVal};
end

for nt = 1:(max_traj_plot-1)
    beforeVals = vals(nt, :);
    afterVals = vals(nt+1, :);

    validPair = isfinite(beforeVals) & isfinite(afterVals);
    beforeVals = beforeVals(validPair);
    afterVals = afterVals(validPair);

    % For q, improvement means q decreases.
    pairedReduction = beforeVals - afterVals;

    medBefore = median(beforeVals, 'omitnan');
    medAfter = median(afterVals, 'omitnan');

    medReduction = median(pairedReduction, 'omitnan');
    q25Reduction = prctile(pairedReduction, 25);
    q75Reduction = prctile(pairedReduction, 75);

    medianBasedPercentReduction = 100 * (medBefore - medAfter) / max(abs(medBefore), eps);

    ChangeRows(end+1, :) = {alpha_val, nt, nt+1, ...
        medBefore, medAfter, medReduction, q25Reduction, q75Reduction, ...
        medianBasedPercentReduction};
end

SummaryTable = cell2table(SummaryRows, ...
    'VariableNames', {'NoiseLevel', 'NumTrajectories', 'Median', 'Q25', 'Q75', 'IQR'});

ChangeTable = cell2table(ChangeRows, ...
    'VariableNames', {'NoiseLevel', 'TrajBefore', 'TrajAfter', ...
    'MedianBefore', 'MedianAfter', 'MedianPairedReduction', ...
    'Q25PairedReduction', 'Q75PairedReduction', ...
    'MedianBasedPercentReduction'});

disp(SummaryTable)
disp(ChangeTable)

writetable(SummaryTable, fullfile(plot_path, ...
    sprintf('%s_CrossSection_SummaryTable_noise_%s.csv', ...
    metric_file_tag, alpha_str)));

writetable(ChangeTable, fullfile(plot_path, ...
    sprintf('%s_CrossSection_ChangeTable_noise_%s.csv', ...
    metric_file_tag, alpha_str)));