clear
close all

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

% NaN-robust median across repeats
A = reshape(cell2mat(CoeffErr_cell), [num_trajs_upto, numNoise, rept]);
CoeffErr_median = median(A, 3, 'omitnan');

CoeffErr_median_plot = CoeffErr_median(1:5,:);

%% E2 Heatmap with optional thresholds

% Thresholds (set to NaN if unused)
low_thresh  = NaN;   % example: NaN = ignore, or 0 = apply
high_thresh = NaN;     % example: NaN = ignore, or 3 = apply

CoeffErr_disp = CoeffErr_median_plot.';  % rows=noise, cols=trajectories

% Apply thresholding for display
Z_plot = CoeffErr_disp;
if ~isnan(low_thresh)
    Z_plot(Z_plot < low_thresh) = low_thresh;
end
if ~isnan(high_thresh)
    Z_plot(Z_plot > high_thresh) = high_thresh;
end

% Plot
figure;
set(gcf, 'Units', 'inches', 'Position', [6 2 5.5 10]); 
imagesc(Z_plot);
colormap parula;
set(gca, 'YDir', 'normal');
c = colorbar;
c.Label.String = 'Median Parameter Error';
c.Label.FontSize = 22;

% Adjust color limits if thresholds are active
if ~isnan(low_thresh) || ~isnan(high_thresh)
    lo = low_thresh; if isnan(lo), lo = min(CoeffErr_disp(:)); end
    hi = high_thresh; if isnan(hi), hi = max(CoeffErr_disp(:)); end
    caxis([lo hi]);
end

% Ticks and tick labels
xticks(1:num_trajs_upto);
xticklabels(arrayfun(@(x) sprintf('%u', x), 1:num_trajs_upto, 'UniformOutput', false));
yticks(1:numNoise);
yticklabels(arrayfun(@(x) sprintf('%.3f', x), noise_ratios, 'UniformOutput', false));
set(gca, 'FontSize', 15);

% Axis labels
xlabel('Number of Trajectories','FontSize',25);

axis tight;

% Annotate threshold-exceeding values
[nrows, ncols] = size(CoeffErr_disp);
for j = 1:nrows
    for i = 1:ncols
        val = CoeffErr_disp(j,i);
        if (~isnan(low_thresh) && val < low_thresh) || ...
           (~isnan(high_thresh) && val > high_thresh)
            if isinf(val)
                txt = 'Inf';
            else
                txt = sprintf('%.1f', val);
            end
            text(i, j, txt, ...
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

figure_name = "E2_Heatmap_100repeats_median_withoutCutoffs.pdf";
exportgraphics(gcf, fullfile(plot_path, figure_name), 'ContentType','vector', 'BackgroundColor','none');

