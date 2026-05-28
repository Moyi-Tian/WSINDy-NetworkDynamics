% Find the q for each truncation percentage and noise level
% Empirical-reference version:
% replace the true coefficient vector by the empirical mean learned vector

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
save_path = [cur_dir,'/WSINDy_data_3StateInput'];

fname = sprintf(['WSINDy Results for Fully-Mixed Network Dynamics E0 Batch Data ', ...
    'numTrajUpTo=%u, E0_=%.5g-%.5g, beta=%.3g, theta=%.3g, eta=%.3g, ', ...
    'gamma_i=%.3g, gamma_p=%.3g, tmax=%u, t_grid_num=%u, NumNoise=%u, ', ...
    'repeat=%u_noise_%.5g-%.5g_MSTLS.mat'], ...
    num_trajs_upto, start_val, end_val, beta, theta, eta, gamma_i, gamma_p, ...
    time_range(end), t_num, numNoise, rept, noise_ratios(1), noise_ratios(end));

load(fullfile(save_path, fname))

plot_trajNum_to = 5;

%% Get q and plot heatmap using empirical mean as reference
q_mat = zeros(plot_trajNum_to, numNoise);

for i = 1:plot_trajNum_to
    for j = 1:numNoise

        % Each column is one learned coefficient vector from one realization
        W = cell2mat(Weights_learned_cell(i,j,1:rept));
        W = reshape(W, [], rept);

        % Replace true W by empirical mean learned W
        W_ref = mean(W, 2, 'omitnan');
        den = norm(W_ref);

        if den > 0
            diffs = W - W_ref;
            q_r = vecnorm(diffs, 2, 1) / den;
        else
            q_r = vecnorm(W, 2, 1);
        end

        q_mat(i,j) = mean(q_r, 'omitnan');
    end
end

%% Plot heatmap
figure;
set(gcf, 'Units', 'inches', 'Position', [6 2 6.5 10]); 

imagesc(q_mat.');
set(gca, 'YDir', 'normal');

% Create and label colorbar
cb = colorbar;
cb.Label.String = '';
cb.Label.Interpreter = 'latex';
cb.Label.FontSize = 30;
cb.Label.FontWeight = 'normal';
cb.Label.Rotation = 0;
cb.Label.VerticalAlignment = 'bottom';
cb.Label.HorizontalAlignment = 'center';
cb.Label.Position(1) = cb.Label.Position(1) - 1.5;
cb.Label.Position(2) = cb.Label.Position(2) + 12.5;

% Ticks and tick labels
xticks(1:plot_trajNum_to);
xticklabels(arrayfun(@(x) sprintf('%u', x), ...
    1:plot_trajNum_to, 'UniformOutput', false));

yticks(1:numNoise);
yticklabels(arrayfun(@(x) sprintf('%.3f', x), ...
    noise_ratios, 'UniformOutput', false));

set(gca, 'FontSize', 20);

% Axis labels
xlabel('Number of Trajectories', 'FontSize', 30);

axis tight;

%% Save plot
plot_path = fullfile(cur_dir, 'plots');
if ~exist(plot_path, 'dir')
    mkdir(plot_path);
end

filename = sprintf('E2_q_empMean_heatmap_%uNoises_%uRepeat_plotUpTo%u_vectorLevel.pdf', ...
    numNoise, rept, plot_trajNum_to);

exportgraphics(gcf, fullfile(plot_path, filename), ...
    'ContentType', 'vector', 'BackgroundColor', 'none');