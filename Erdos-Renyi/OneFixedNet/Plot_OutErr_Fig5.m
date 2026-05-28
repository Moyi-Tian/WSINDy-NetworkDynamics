tic
clear
close all

%% cd to current directory
scriptFullName = matlab.desktop.editor.getActiveFilename();
scriptDir = fileparts(scriptFullName);
% change to current folder
if ~isempty(scriptDir)
    cd(scriptDir);
end

%% load data
cur_dir=pwd;

data_path=[cur_dir,'/WSINDy_data_3StateInput'];

N = 1000;
k = 5;

network_id = 1;

beta = 0.2; % beta: transmission rate due to online network NotEngaged/Uninterested(U) -> Engaged(E)
theta = 0.5; % theta: transmission rate due to offline protesting ratio NotEngaged/Uninterested(U) -> Engaged(E)
eta = 0.4; % eta: fraction of Engaged(E) online from UnProtesting(UP) -> Protesting(P) offline
gamma_i = 0.1; % gamma_i: recovery rate online Engaged(E) -> DoneEngaging/DisEngaged(D)
gamma_p = 0.2; % gamma_p: recovery rate offline Protesting(P) -> DoneProtesting(R)

t_max = 200;
t_end = 40;

Initial_E_percents = linspace(0.01,0.3,30);
len_E_percents = length(Initial_E_percents);

NumIteration = 100;
NumGrids = 1000;

start_val = 0.01;
end_val = 0.3;
J = find(Initial_E_percents >= start_val & Initial_E_percents <= end_val);
num_ICs = numel(J);
num_trajs_upto = 10;

noise_ratio = 0;

rept = 1;

fname = sprintf('WSINDy Results for Gillepsie batched E0 on Fixed Erdos-Renyi N=%u, k=%u, numTrajUpTo=%u, E0_=%.5g-%.5g, beta=%.3g, theta=%.3g, eta=%.3g, gamma_i=%.3g, gamma_p=%.3g, t_max=%u, t_end=%.3g, t_grid_num=%u, noise=%.5g, repeat=%u, iter=%u_MSTLS_3rdLib - NetID-%u.mat',N,k,num_trajs_upto,start_val,end_val,beta,theta,eta,gamma_i,gamma_p,t_max,t_end,NumGrids,noise_ratio,rept,NumIteration,network_id);

load(fullfile(data_path, fname))

%% E0 values
x_values = 1:num_trajs_upto;
plot_path = fullfile(cur_dir, sprintf('OutputError_Plots'));
if ~exist(plot_path, 'dir')
    mkdir(plot_path);
end

%% OutErr_avgs
[nRows, ~] = size(trajectory_errors_cell);
nStates = 3;
state_labels = {'$U$', '$E$', '$P$'};

% Preallocate median error values for each state
TrajectoryErr_mean = zeros(nRows, nStates);

% Compute median for each row and each state
for i = 1:nRows
    trajectory_errors_cur = trajectory_errors_cell{i};
    [num_trajs, ~] = size(trajectory_errors_cur);
    state_vals = nan(num_trajs, nStates);
    for j = 1:num_trajs
        err_vec = trajectory_errors_cur(j,:);  % 1 x 3
        if ~isempty(err_vec) && all(isnumeric(err_vec))
            state_vals(j, :) = err_vec;
        end
    end
    % Compute median across the runs for each state, ignoring NaNs
    TrajectoryErr_mean(i, :) = mean(state_vals, 'omitnan');
end

% Plot
markers = {'o','s','^'};   % reuse or redefine
colors  = lines(nStates);

figure; hold on;
ax = gca; hold(ax,'on'); box(ax,'on'); grid(ax,'on');
for s = 1:nStates
    plot(x_values, TrajectoryErr_mean(:,s), ...
        'LineWidth', 2, ...
        'Marker', markers{s}, ...
        'MarkerSize', 12,...
        'Color', colors(s,:), ...
        'DisplayName', state_labels{s});
end
hold off;

set(ax,'FontSize',18,'TickLabelInterpreter','latex');
legend('Interpreter','latex','Location','best','FontSize',22);
xlabel('Number of Trajectories','FontSize',22);
ylabel('Output Error (log)','FontSize',22);
xlim([1 num_trajs_upto]);
set(gca,'YScale','log');

% Save plot
filename = sprintf('MeanOutErr_N=%u_k=%u_noise=%.5g_rept=%u_iter=%u_logScale_new.pdf',N,k,noise_ratio,rept,NumIteration);
exportgraphics(gcf, fullfile(plot_path,filename), 'ContentType', 'vector', 'BackgroundColor', 'none', 'Resolution', 300);


toc