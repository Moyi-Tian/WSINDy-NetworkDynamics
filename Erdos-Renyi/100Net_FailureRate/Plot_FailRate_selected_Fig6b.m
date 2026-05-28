tic
clear; close all

%% cd to current directory
scriptFullName = matlab.desktop.editor.getActiveFilename();
scriptDir = fileparts(scriptFullName);
if ~isempty(scriptDir)
    cd(scriptDir);
end

%% load data
cur_dir  = pwd;
data_path = fullfile(cur_dir,'FailRate_data');

%% --- parameters used in filename ---
N = 1000;

beta = 0.2; theta = 0.5; eta = 0.4;
gamma_i = 0.1; gamma_p = 0.2;

t_max = 200;
t_end = 40;

Initial_E_percents = linspace(0.01,0.3,30);
start_val = 0.01;
end_val   = 0.3;
J = find(Initial_E_percents >= start_val & Initial_E_percents <= end_val); %#ok<NASGU>
num_ICs = numel(J); %#ok<NASGU>

NumIteration     = 100;
NumGrids         = 1000;
NumNetworkRepeat = 100;

noise_ratio = 0;
rept = 1;

num_trajs_upto_load = 10;   % what's stored in file
num_trajs_upto      = 4;    % what you want to plot
x_values = 1:num_trajs_upto;

%% --- choose which k values to overlay ---
k_list = [4 6 8 10];   % <<< change this to whatever you want

%% --- plot folder ---
plot_path = fullfile(cur_dir, sprintf('Line_Plots_FailureRate_traj=1-%u_overlay', num_trajs_upto));
if ~exist(plot_path,'dir'); mkdir(plot_path); end

%% --- styles for overlay (reused cyclically if k_list is long) ---
colors     = lines(max(7, numel(k_list)));         % distinct colors
lineStyles = {'-','--',':','-.'};
markers    = {'o','s','^','d','v','>','<','p','h','x','+'};

%% --- helper: load one k safely (only needed vars, and truncate) ---
load_failrates = @(k) load( fullfile(data_path, sprintf( ...
    'FailureRate Gillepsie batched E0 on %u Fixed Erdos-Renyi N=%u, k=%u, numTrajUpTo=%u, E0_=%.5g-%.5g, beta=%.3g, theta=%.3g, eta=%.3g, gamma_i=%.3g, gamma_p=%.3g, t_max=%u, t_end=%.3g, t_grid_num=%u, noise=%.5g, repeat=%u, iter=%u_MSTLS_3rdLib.mat', ...
    NumNetworkRepeat,N,k,num_trajs_upto_load,start_val,end_val,beta,theta,eta,gamma_i,gamma_p,t_max,t_end,NumGrids,noise_ratio,rept,NumIteration)), ...
    'fail_rate1','fail_rate2','fail_rate3');

%% ---------------------------
% Plot Failure Rate 1 overlay
%% ---------------------------
fig = figure('Units','inches','Position',[6 2 7 5]);
ax = gca; hold(ax,'on'); box(ax,'on'); grid(ax,'on');

for ii = 1:numel(k_list)
    k = k_list(ii);

    S = load_failrates(k);

    y = S.fail_rate1(:);
    L = min(num_trajs_upto, numel(y));
    y = y(1:L);

    c  = colors(ii,:);
    ls = lineStyles{mod(ii-1, numel(lineStyles))+1};
    mk = markers{mod(ii-1, numel(markers))+1};

    plot(x_values(1:L), y, ...
        'Color', c, 'LineStyle', ls, 'Marker', mk, ...
        'LineWidth', 2, 'MarkerSize', 7, ...
        'DisplayName', sprintf('$k=%d$', k));
end

set(ax,'FontSize',14,'TickLabelInterpreter','latex');
xlabel('Number of Trajectories','FontSize',18,'Interpreter','latex');
ylabel('Failure Rate 1','FontSize',18,'Interpreter','latex');
xlim([1 num_trajs_upto]);
xticks(1:num_trajs_upto);              % <<< ADD THIS
xticklabels(string(1:num_trajs_upto)); % <<< optional but robust
ylim([0 1]);
legend('Location','best','Interpreter','latex');

filename1 = sprintf('FailureRate1_overlay_k=%s_%uFixedER_N=%u_noise=%.5g_rept=%u_iter=%u.pdf', ...
    strrep(num2str(k_list),'  ','_'), NumNetworkRepeat,N,noise_ratio,rept,NumIteration);
exportgraphics(fig, fullfile(plot_path,filename1), 'ContentType','vector','BackgroundColor','none','Resolution',300);

%% ---------------------------
% Plot Failure Rate 2 overlay
%% ---------------------------
fig = figure('Units','inches','Position',[6 2 7 5]);
ax = gca; hold(ax,'on'); box(ax,'on'); grid(ax,'on');

for ii = 1:numel(k_list)
    k = k_list(ii);

    S = load_failrates(k);

    y = S.fail_rate2(:);
    L = min(num_trajs_upto, numel(y));
    y = y(1:L);

    c  = colors(ii,:);
    ls = lineStyles{mod(ii-1, numel(lineStyles))+1};
    mk = markers{mod(ii-1, numel(markers))+1};

    plot(x_values(1:L), y, ...
        'Color', c, 'LineStyle', ls, 'Marker', mk, ...
        'LineWidth', 2, 'MarkerSize', 7, ...
        'DisplayName', sprintf('$k=%d$', k));
end

set(ax,'FontSize',14,'TickLabelInterpreter','latex');
xlabel('Number of Trajectories','FontSize',18,'Interpreter','latex');
ylabel('Failure Rate 2','FontSize',18,'Interpreter','latex');
xlim([1 num_trajs_upto]);
xticks(1:num_trajs_upto);              % <<< ADD THIS
xticklabels(string(1:num_trajs_upto)); % <<< optional but robust

ylim([0 1]);
legend('Location','best','Interpreter','latex');

filename2 = sprintf('FailureRate2_overlay_k=%s_%uFixedER_N=%u_noise=%.5g_rept=%u_iter=%u.pdf', ...
    strrep(num2str(k_list),'  ','_'), NumNetworkRepeat,N,noise_ratio,rept,NumIteration);
exportgraphics(fig, fullfile(plot_path,filename2), 'ContentType','vector','BackgroundColor','none','Resolution',300);

%% ---------------------------
% Plot Failure Rate 3 overlay
%% ---------------------------
fig = figure('Units','inches','Position',[6 2 7 8]);
ax = gca; hold(ax,'on'); box(ax,'on'); grid(ax,'on');

for ii = 1:numel(k_list)
    k = k_list(ii);

    S = load_failrates(k);

    y = S.fail_rate3(:);
    L = min(num_trajs_upto, numel(y));
    y = y(1:L);

    c  = colors(ii,:);
    ls = lineStyles{mod(ii-1, numel(lineStyles))+1};
    mk = markers{mod(ii-1, numel(markers))+1};

    plot(x_values(1:L), y, ...
        'Color', c, 'LineStyle', ls, 'Marker', mk, ...
        'LineWidth', 2, 'MarkerSize', 10, ...
        'DisplayName', sprintf('$k=%d$', k));
end

set(ax,'FontSize',18);
xlabel('Number of Trajectories','FontSize',23);
ylabel('Failure Rate','FontSize',23);
xlim([1 num_trajs_upto]);
xticks(1:num_trajs_upto);              % <<< ADD THIS
xticklabels(string(1:num_trajs_upto)); % <<< optional but robust

legend('Location','best','Interpreter','latex','FontSize',23);

filename3 = sprintf('FailureRate3_overlay_k=%s_%uFixedER_N=%u_noise=%.5g_rept=%u_iter=%u.pdf', ...
    strrep(num2str(k_list),'  ','_'), NumNetworkRepeat,N,noise_ratio,rept,NumIteration);
exportgraphics(fig, fullfile(plot_path,filename3), 'ContentType','vector','BackgroundColor','none','Resolution',300);

toc
