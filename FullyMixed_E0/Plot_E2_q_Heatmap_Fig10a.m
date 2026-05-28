% Find the q for each truncation percentage and noise level

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

plot_trajNum_to = 5;

%% Get q and plot heatmap - vector level
q_mat = zeros(plot_trajNum_to,numNoise);
den = norm(weights_true);
for i=1:plot_trajNum_to
    for j=1:numNoise

        W = cell2mat(Weights_learned_cell(i,j,1:rept));
        W = reshape(W, [], rept);

        if den > 0
            diffs = W - weights_true;    % implicit expansion (R2016b+)
            q_r = vecnorm(diffs, 2, 1) / den;
        else
            % If w_true is the zero vector, use absolute L2 norm:
            % q_r(r) = ‖w_hat(r)‖
            q_r = vecnorm(W, 2, 1);
        end

        % Average over repetitions (ignore any NaNs)
        q_mat(i,j) = mean(q_r, 'omitnan');
    end
end

figure;
set(gcf, 'Units', 'inches', 'Position', [6 2 7 10]); 
imagesc(q_mat.');     % transpose: rows=noise level (Y), cols=slice (X)
set(gca, 'YDir', 'normal');   % row 1 at bottom

% Create and label colorbar
cb = colorbar;
cb.Label.String = '$q$';             % LaTeX-style math label
cb.Label.Interpreter = 'latex';      % enable LaTeX rendering
cb.Label.FontSize = 30;              % adjust size
cb.Label.FontWeight = 'normal';
cb.Label.Rotation = 0;               % horizontal orientation
cb.Label.VerticalAlignment = 'bottom';
cb.Label.HorizontalAlignment = 'center';
cb.Label.Position(1) = cb.Label.Position(1) - 1.5;  % optional: lift label upward
cb.Label.Position(2) = cb.Label.Position(2) + 12.5;  % optional: lift label upward


% Ticks and tick labels
xticks(1:plot_trajNum_to);
xticklabels(arrayfun(@(x) sprintf('%u', x), 1:plot_trajNum_to, 'UniformOutput', false));
yticks(1:numNoise);
yticklabels(arrayfun(@(x) sprintf('%.3f', x), noise_ratios, 'UniformOutput', false));
set(gca, 'FontSize', 20);

% Axis labels
xlabel('Number of Trajectories','FontSize',30);
ylabel('Noise Level','FontSize',30);

axis tight;


%% save plot
plot_path = fullfile(cur_dir, 'plots');
if ~exist(plot_path, 'dir'), mkdir(plot_path); end
filename = sprintf('E2_q_heatmap_%uNoises_%uRepeat_plotUpTo%u_vectorLevel.pdf', ...
                    numNoise, rept,plot_trajNum_to);
exportgraphics(gcf, fullfile(plot_path, filename), ...
    'ContentType','vector','BackgroundColor','none');

%% Plot q surface in 3D
% Build meshgrid from x and y axis
[X, Y] = meshgrid(1:plot_trajNum_to, noise_ratios);

% Plot surface
figure;
surf(X, Y, q_mat');

set(gca, 'FontSize', 12);

xlabel('Number of Trajectories','FontSize',15,'Rotation',335);
ylabel('Noise Level','FontSize',15,'Rotation',20);
zlabel('q','FontSize',15,'Rotation',0);

ax = gca;
% x label location
posX = ax.XLabel.Position;
posX(3) = posX(3) + 0.2*range(zlim);   % lift along z
posX(2) = posX(2) - 0.05*range(ylim);  % pull closer in y, if needed
ax.XLabel.Position = posX;
% z label position
posZ = ax.ZLabel.Position;
posZ(1) = posZ(1) + 0.05*range(xlim);
posZ(2) = posZ(2) - 1.2*range(ylim);
ax.ZLabel.Position = posZ;

shading interp;         % smooth shading (optional)
colorbar;               % keep a colorbar for reference
colormap(parula);       % choose a colormap
view(45,30);            % set 3D viewing angle


%% Plot contour map of q surface
figure;
set(gcf, 'Units', 'inches', 'Position', [6 2 9 7.5]); 

% Build meshgrid consistent with the surface plot
[X, Y] = meshgrid(1:plot_trajNum_to, noise_ratios);

% Specify contour levels to display
contour_levels = [0.01 0.1 0.2 0.5 1 2 5 10];  

% Filled contour plot (smooth)
contourf(X, Y, q_mat', contour_levels, 'LineStyle', 'none');
hold on;

% Draw contour lines
[C, h] = contour(X, Y, q_mat', contour_levels, ...
                 'LineColor', 'k', 'LineWidth', 1.2);

% Label contour lines and store handles
% clabel(C, h, 'FontSize', 13, 'FontWeight', 'bold', 'LabelSpacing', 800);
clabel(C, h, 'FontSize', 15, 'FontWeight', 'bold', 'LabelSpacing', 800, 'Color', '#FF5F1F');

% Axes and appearance
xlabel('Number of Trajectories', 'FontSize', 15);
ylabel('Noise level', 'FontSize', 15);
title('Contour map of q', 'FontSize', 16);
set(gca, 'FontSize', 12, 'YDir', 'normal');
colormap(parula);
colorbar;
axis tight;
grid off;
box on;

%% save as PDF
plot_path = fullfile(cur_dir, 'plots');
if ~exist(plot_path, 'dir'), mkdir(plot_path); end
filename = sprintf('E2_q_contour_%uNoises_%uRepeat.pdf', numNoise, rept);
exportgraphics(gcf, fullfile(plot_path, filename), ...
    'ContentType','vector','BackgroundColor','none');







