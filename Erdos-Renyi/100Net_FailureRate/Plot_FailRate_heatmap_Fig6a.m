%% Failure Rate Heatmap across k and # of trajectories
clear; close all; tic

%% --- cd to current directory ---
scriptFullName = matlab.desktop.editor.getActiveFilename();
scriptDir = fileparts(scriptFullName);
if ~isempty(scriptDir)
    cd(scriptDir);
end

cur_dir = pwd;
data_path = fullfile(cur_dir, 'FailRate_data');

%% --- parameters ---
N = 1000;
beta = 0.2; theta = 0.5; eta = 0.4;
gamma_i = 0.1; gamma_p = 0.2;
t_max = 200; t_end = 40;
NumNetworkRepeat = 100;
NumIteration = 100;
NumGrids = 1000;

noise_ratio = 0;
num_trajs_upto = 4;
k_values = 4:10;  
num_k = length(k_values);

Initial_E_percents = linspace(0.01,0.3,30);

%% --- Preallocate 3 matrices ---
FailureRate_all_1 = nan(num_k, num_trajs_upto);
FailureRate_all_2 = nan(num_k, num_trajs_upto);
FailureRate_all_3 = nan(num_k, num_trajs_upto);

%% --- Loop over k values ---
for kk = 1:num_k
    k = k_values(kk);

    fname_pattern = sprintf(['FailureRate Gillepsie batched E0 on %u Fixed Erdos-Renyi N=%u, k=%u'], ...
        NumNetworkRepeat, N, k);

    file_list = dir(fullfile(data_path, [fname_pattern, '*.mat']));

    if isempty(file_list)
        warning('No file found for k = %d', k);
        continue;
    end

    fname_full = fullfile(data_path, file_list(1).name);
    fprintf('Loading file: %s\n', file_list(1).name);
    S = load(fname_full, 'fail_rate1','fail_rate2','fail_rate3');

    if isfield(S,'fail_rate1')
        FailureRate_all_1(kk, 1:min(num_trajs_upto, numel(S.fail_rate1))) = S.fail_rate1(1:min(num_trajs_upto, numel(S.fail_rate1))).';
    end
    if isfield(S,'fail_rate2')
        FailureRate_all_2(kk, 1:min(num_trajs_upto, numel(S.fail_rate2))) = S.fail_rate2(1:min(num_trajs_upto, numel(S.fail_rate2))).';
    end
    if isfield(S,'fail_rate3')
        FailureRate_all_3(kk, 1:min(num_trajs_upto, numel(S.fail_rate3))) = S.fail_rate3(1:min(num_trajs_upto, numel(S.fail_rate3))).';
    end

end

%% --- Plotting setup ---
x_values = 1:num_trajs_upto;
plot_path = fullfile(cur_dir, 'Plots_FailureRateHeatmap_traj=1-4');

if ~exist(plot_path, 'dir')
    mkdir(plot_path);
end

%% Helper function
% -----------------------------------------------------
% Helper function for plotting heatmaps
% -----------------------------------------------------
function plot_heatmap(data_matrix, title_text, filename, x_values, k_values, plot_path)
    figure('Units','inches','Position',[6 2 6 8.5]);

    imagesc(x_values, k_values, data_matrix);
    set(gca,'YDir','normal');
    set(gca,'FontSize',23);

    colormap("parula"); 
    c = colorbar;
    c.Label.String = 'Mean Failure Rate';
    c.Label.Interpreter = 'none';
    c.Label.FontSize = 24;


    xlabel('Number of Trajectories','FontSize',24);
    ylabel('$k$','Interpreter','latex','FontSize',28);

    exportgraphics(gcf, fullfile(plot_path, filename), ...
        'ContentType','vector','BackgroundColor','none','Resolution',300);
end

%% --- Plot the 3 heatmaps ---
plot_heatmap(FailureRate_all_1, 'Failure Rate 1', ...
    sprintf('FailureRate1Heatmap_k=%u-%u.pdf', k_values(1), k_values(end)), ...
    x_values, k_values, plot_path);

plot_heatmap(FailureRate_all_2, 'Failure Rate 2', ...
    sprintf('FailureRate2Heatmap_k=%u-%u.pdf', k_values(1), k_values(end)), ...
    x_values, k_values, plot_path);

plot_heatmap(FailureRate_all_3, 'Failure Rate 3', ...
    sprintf('FailureRate3Heatmap_k=%u-%u.pdf', k_values(1), k_values(end)), ...
    x_values, k_values, plot_path);

toc
