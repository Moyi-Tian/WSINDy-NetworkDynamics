%% Cross-sectional variability plot for TP score heatmap
% This script loads the same data as Figure 2(a), extracts selected
% noise-level cross sections, and plots distributions over 100 repeats.

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
E0 = linspace(0.01,0.99,99);

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

%% Data array
% A has dimensions:
% number of trajectories x number of noise levels x number of repeats
A = reshape(cell2mat(Tps_cell), [num_trajs_upto, numNoise, rept]);

%% User choices
max_traj_plot = 5;
selected_noise_levels = [0.02];

metric_label = 'TPR';
metric_file_tag = 'TPR';

% Use "higher" for TP score.
% For error metrics, change this to "lower".
metric_direction = "higher";

%% Find nearest available noise indices
selected_noise_ids = zeros(size(selected_noise_levels));
for ii = 1:length(selected_noise_levels)
    [~, selected_noise_ids(ii)] = min(abs(noise_ratios - selected_noise_levels(ii)));
end

selected_noise_levels_actual = noise_ratios(selected_noise_ids);

%% Output folder
plot_path = fullfile(cur_dir, 'plots_variability_appendix');
if ~exist(plot_path, 'dir')
    mkdir(plot_path);
end

%% Make cross-sectional boxplot figure
figure;
set(gcf, 'Units', 'inches', 'Position', [3 3 8 6.5]);

SummaryRows = {};
ChangeRows = {};

for ii = 1:length(selected_noise_ids)

    noise_idx = selected_noise_ids(ii);
    alpha_val = selected_noise_levels_actual(ii);

    vals = squeeze(A(1:max_traj_plot, noise_idx, :));
    % vals is max_traj_plot x rept

    subplot(1, length(selected_noise_ids), ii);
    hold on

    dataVec = [];
    groupVec = [];

    for nt = 1:max_traj_plot
        curVals = vals(nt, :).';
        dataVec = [dataVec; curVals];
        groupVec = [groupVec; nt * ones(length(curVals), 1)];
    end

    keep = isfinite(dataVec);
    boxplot(dataVec(keep), groupVec(keep), ...
        'Labels', arrayfun(@num2str, 1:max_traj_plot, 'UniformOutput', false), ...
        'Symbol', '.', ...
        'Widths', 0.55);

    medVals = median(vals, 2, 'omitnan');
    plot(1:max_traj_plot, medVals, 'ko-', ...
        'LineWidth', 1.5, ...
        'MarkerFaceColor', 'k', ...
        'MarkerSize', 5);

    set(gca, 'FontSize', 18);

    xlabel('Number of Trajectories', 'FontSize', 22);
    ylabel(metric_label, 'FontSize', 22);

    box on

    if strcmp(metric_file_tag, 'TPR')
        ylim([0 1.05]);
    end

    %% Summary table for medians and IQRs
    for nt = 1:max_traj_plot
        curVals = vals(nt, :);
        curVals = curVals(isfinite(curVals));

        medVal = median(curVals, 'omitnan');
        q25 = prctile(curVals, 25);
        q75 = prctile(curVals, 75);
        iqrVal = q75 - q25;

        SummaryRows(end+1, :) = {alpha_val, nt, medVal, q25, q75, iqrVal};
    end

    %% Paired changes between consecutive trajectory counts
    for nt = 1:(max_traj_plot-1)
        beforeVals = vals(nt, :);
        afterVals = vals(nt+1, :);

        validPair = isfinite(beforeVals) & isfinite(afterVals);
        beforeVals = beforeVals(validPair);
        afterVals = afterVals(validPair);

        if strcmp(metric_direction, "higher")
            pairedImprovement = afterVals - beforeVals;
            medianBasedPercentChange = 100 * ...
                (median(afterVals, 'omitnan') - median(beforeVals, 'omitnan')) / ...
                max(abs(median(beforeVals, 'omitnan')), eps);
        else
            pairedImprovement = beforeVals - afterVals;
            medianBasedPercentChange = 100 * ...
                (median(beforeVals, 'omitnan') - median(afterVals, 'omitnan')) / ...
                max(abs(median(beforeVals, 'omitnan')), eps);
        end

        medImprove = median(pairedImprovement, 'omitnan');
        q25Improve = prctile(pairedImprovement, 25);
        q75Improve = prctile(pairedImprovement, 75);

        ChangeRows(end+1, :) = {alpha_val, nt, nt+1, ...
            median(beforeVals, 'omitnan'), median(afterVals, 'omitnan'), ...
            medImprove, q25Improve, q75Improve, medianBasedPercentChange};
    end
end

%% Save figure
figure_name = sprintf('%s_CrossSection_Boxplots_SelectedNoise.pdf', metric_file_tag);
exportgraphics(gcf, fullfile(plot_path, figure_name), ...
    'ContentType', 'vector', ...
    'BackgroundColor', 'none', ...
    'Resolution', 300);

%% Save tables
SummaryTable = cell2table(SummaryRows, ...
    'VariableNames', {'NoiseLevel', 'NumTrajectories', 'Median', 'Q25', 'Q75', 'IQR'});

ChangeTable = cell2table(ChangeRows, ...
    'VariableNames', {'NoiseLevel', 'TrajBefore', 'TrajAfter', ...
    'MedianBefore', 'MedianAfter', 'MedianPairedImprovement', ...
    'Q25PairedImprovement', 'Q75PairedImprovement', ...
    'MedianBasedPercentChange'});

disp(SummaryTable)
disp(ChangeTable)

writetable(SummaryTable, fullfile(plot_path, ...
    sprintf('%s_CrossSection_SummaryTable.csv', metric_file_tag)));

writetable(ChangeTable, fullfile(plot_path, ...
    sprintf('%s_CrossSection_ChangeTable.csv', metric_file_tag)));