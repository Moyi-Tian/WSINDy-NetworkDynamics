clear
close all

%% cd to current directory
scriptFullName = matlab.desktop.editor.getActiveFilename();
scriptDir = fileparts(scriptFullName);
% change to current folder
if ~isempty(scriptDir)
    cd(scriptDir);
end

%% Load Data
cur_dir=pwd;

beta = 0.5; % beta: transmission rate due to online network NotEngaged/Uninterested(U) -> Engaged(E)
theta = 0.4; % theta: transmission rate due to offline protesting ratio NotEngaged/Uninterested(U) -> Engaged(E)
eta = 0.2; % eta: fraction of Engaged(E) online from UnProtesting(UP) -> Protesting(P) offline
gamma_i = 0.1; % gamma_i: recovery rate online Engaged(E) -> DoneEngaging/DisEngaged(D)
gamma_p = 0.3; % gamma_p: recovery rate offline Protesting(P) -> DoneProtesting(R)

time_range = [0 50];
t_num = 1000;
E0 = linspace(0.01,0.99,99);

start_val = 0.01;
end_val = 0.3;
num_trajs_upto = 10;

noise_ratios = linspace(0,0.1,21);
numNoise = length(noise_ratios);

rept = 100;
save_path=[cur_dir,'/WSINDy_data_3StateInput'];
fname = sprintf('WSINDy Results for Fully-Mixed Network Dynamics E0 Batch Data numTrajUpTo=%u, E0_=%.5g-%.5g, beta=%.3g, theta=%.3g, eta=%.3g, gamma_i=%.3g, gamma_p=%.3g, tmax=%u, t_grid_num=%u, NumNoise=%u, repeat=%u_noise_%.5g-%.5g_MSTLS.mat',num_trajs_upto,start_val,end_val,beta,theta,eta,gamma_i,gamma_p,time_range(end),t_num,numNoise,rept,noise_ratios(1),noise_ratios(end));
load(fullfile(save_path, fname))

%% data
% NaN-robust mean across the 3rd dim (k)
A = reshape(cell2mat(Tps_cell), [num_trajs_upto, numNoise, rept]);   % to numeric
TPS_median = median(A, 3, 'omitnan');    % slice by noise level double

TPS_median_plot = TPS_median(1:5,:);

%% TPS Heatmap
% Thresholds
low_thresh  = 0;   % values <= this show bold light green
high_thresh = 1;    % values >= this show bold red

% Heatmap with (1,1) at lower-left and slice on X, noise level on Y
figure; 
set(gcf, 'Units', 'inches', 'Position', [6 2 5.5 10]);
imagesc(TPS_median_plot.');     % transpose: rows=noise level (Y), cols=slice (X)
colormap(flipud(parula));
set(gca, 'YDir', 'normal');   % row 1 at bottom
c = colorbar;
c.Label.String = 'Median TPR';
c.Label.FontSize = 22;

% Ticks and tick labels
xticks(1:num_trajs_upto);
xticklabels(arrayfun(@(x) sprintf('%u', x), 1:num_trajs_upto, 'UniformOutput', false));
yticks(1:numNoise);
yticklabels(arrayfun(@(x) sprintf('%.3f', x), noise_ratios, 'UniformOutput', false));
set(gca, 'FontSize', 15);

% Axis labels and title
xlabel('Number of Trajectories','FontSize',25);
ylabel('Noise Level','FontSize',25);

axis tight;

% Annotate values
CoeffErr_disp = TPS_median_plot.';  % matrix actually plotted (rows=noise, cols=slice)
[nrows, ncols] = size(CoeffErr_disp);

for j = 1:nrows          % row = noise index (Y)
    for i = 1:ncols      % col = slice index (X)
        val = CoeffErr_disp(j,i);
        if val < low_thresh
            text(i, j, sprintf('%.2f', val), ...
                'Color',[0.7 1 0.7], ...   % very light green
                'FontWeight','bold', ...
                'HorizontalAlignment','center', ...
                'VerticalAlignment','middle');
        elseif val > high_thresh
            text(i, j, sprintf('%.2f', val), ...
                'Color','r', ...
                'FontWeight','bold', ...
                'HorizontalAlignment','center', ...
                'VerticalAlignment','middle');
        end
    end
end

%% Force identical heatmap and colorbar geometry across paired figures
ax = gca;

ax.Units = 'normalized';
ax.Position = [0.18 0.12 0.58 0.82];

c.Units = 'normalized';
c.Position = [0.82 0.12 0.04 0.82];

%% Save
% ensure plot path exists
plot_path = fullfile(cur_dir, 'plots');
if ~exist(plot_path, 'dir'), mkdir(plot_path); end

figure_name = "TPS_Heatmap_100repeats_median.pdf";
exportgraphics(gcf, fullfile(plot_path, figure_name), 'ContentType','vector', 'BackgroundColor','none');


