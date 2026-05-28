% Compute and save comparison payload: true (mean+envelope), inferred, mean-field
tic
clear
close all

%% cd to current directory
scriptFullName = matlab.desktop.editor.getActiveFilename();
scriptDir = fileparts(scriptFullName);
if ~isempty(scriptDir), cd(scriptDir); end
cur_dir = pwd;

%% Paths
addpath(genpath(cur_dir), '../../../Generate_Synthetic_Data/Gillespie_on_ErdosRenyi/data');
data_path = fullfile(cur_dir, 'WSINDy_data_3StateInput');

%% Parameters (same as your script)
N = 1000;
k = 5;
network_id = 1;

beta = 0.2;
theta = 0.5;
eta = 0.4;
gamma_i = 0.1;
gamma_p = 0.2;

t_max = 200;
t_end = 40;

Initial_E_percents = linspace(0.01,0.3,30);
len_E_percents = length(Initial_E_percents);

start_val = 0.01;
end_val   = 0.3;
J = find(Initial_E_percents >= start_val & Initial_E_percents <= end_val);
ICs = Initial_E_percents(J);

NumIteration = 100;
NumGrids = 1000;

num_ICs = numel(J);
num_trajs_upto = 10;

noise_ratio = 0;
rept = 1;

%% Load original dynamics
fname_data = sprintf(['Dynamics from Gillepsie batched E0 on Fixed Erdos-Renyi N=%u, k=%u, ' ...
    'init_E=%.5g-%.5g, num_IC=%u, beta=%.3g, theta=%.3g, eta=%.3g, gamma_i=%.3g, gamma_p=%.3g, ' ...
    'tmax=%u, tend=%.3g, num_t=%u, iter=%u_ThreshT - NetID-%u.mat'], ...
    N,k,Initial_E_percents(1),Initial_E_percents(end),len_E_percents, ...
    beta,theta,eta,gamma_i,gamma_p,t_max,t_end,NumGrids,NumIteration,network_id);

load(fname_data);  % expects: t_grid, X6_cell, X6_mins, X6_maxes, G, etc.
t_grid_saved = t_grid;
time_grid = t_grid_saved(:);

%% Load learned WSINDy results
fname_learned = sprintf(['WSINDy Results for Gillepsie batched E0 on Fixed Erdos-Renyi N=%u, k=%u, ' ...
    'numTrajUpTo=%u, E0_=%.5g-%.5g, beta=%.3g, theta=%.3g, eta=%.3g, gamma_i=%.3g, gamma_p=%.3g, ' ...
    't_max=%u, t_end=%.3g, t_grid_num=%u, noise=%.5g, repeat=%u, iter=%u_MSTLS_3rdLib - NetID-%u.mat'], ...
    N,k,num_trajs_upto,start_val,end_val,beta,theta,eta,gamma_i,gamma_p,t_max,t_end,NumGrids,noise_ratio,rept,NumIteration,network_id);

load(fullfile(data_path, fname_learned)); % expects: WS_cell, Uobj_cell, t_grid (same grid)
% keep time_grid from loaded dynamics for consistency
% (but also save the learned t_grid if you want)
t_grid_learned = t_grid;

%% Progressive selection (needed later for mapping)
S = progressive_ic_selection(ICs, num_trajs_upto);

%% -------- Mean-field single-level dynamics (same as your code) --------
Usol_mf = zeros(num_ICs,length(time_grid));
Esol_mf = zeros(num_ICs,length(time_grid));
Dsol_mf = zeros(num_ICs,length(time_grid));
Psol_mf = zeros(num_ICs,length(time_grid));
Rsol_mf = zeros(num_ICs,length(time_grid));

pars(1) = beta;
pars(2) = theta;
pars(3) = eta;
pars(4) = gamma_i;
pars(5) = gamma_p;
pars(6) = N;

degs = sum(G);
maxk = max(degs);
pars(7) = maxk;
counts = histcounts(degs,[0:maxk+1]);
pars(8) = sum([0:maxk].*counts);

for i=1:num_ICs
    Usol_mf_all_repts = zeros(NumIteration,length(time_grid));
    Esol_mf_all_repts = zeros(NumIteration,length(time_grid));
    Dsol_mf_all_repts = zeros(NumIteration,length(time_grid));
    Psol_mf_all_repts = zeros(NumIteration,length(time_grid));
    Rsol_mf_all_repts = zeros(NumIteration,length(time_grid));

    Initial_E_percent = Initial_E_percents(J(i));

    E0 = (Initial_E_percent*counts)';
    U0 = counts'-E0;
    D0 = zeros(maxk+1,1);
    P0 = zeros(maxk+1,1);
    R0 = zeros(maxk+1,1);

    for r=1:NumIteration
        y0 = [U0;E0;D0;P0;R0];
        opts = odeset('AbsTol',1e-15,'RelTol',1e-12,'Stats','off');
        odeFunc = @(t,y) HeterogeneousSingleLevelApprox_ODE(t,y,pars);
        [~,ysol] = ode45(odeFunc, time_grid, y0, opts);

        Usol_mf_all_repts(r,:) = sum(ysol(:,1:(maxk+1)),2)/N;
        Esol_mf_all_repts(r,:) = sum(ysol(:,(maxk+2):(2*(maxk+1))),2)/N;
        Dsol_mf_all_repts(r,:) = sum(ysol(:,(2*(maxk+1)+1):(3*(maxk+1))),2)/N;
        Psol_mf_all_repts(r,:) = sum(ysol(:,(3*(maxk+1)+1):(4*(maxk+1))),2)/N;
        Rsol_mf_all_repts(r,:) = sum(ysol(:,(4*(maxk+1)+1):(5*(maxk+1))),2)/N;
    end

    Usol_mf(i,:) = mean(Usol_mf_all_repts,1);
    Esol_mf(i,:) = mean(Esol_mf_all_repts,1);
    Dsol_mf(i,:) = mean(Dsol_mf_all_repts,1);
    Psol_mf(i,:) = mean(Psol_mf_all_repts,1);
    Rsol_mf(i,:) = mean(Rsol_mf_all_repts,1);
end

%% -------- Inferred trajectories (same as your code) --------
nstates = 3;
U_inferred_trajs = cell(num_trajs_upto,1);
E_inferred_trajs = cell(num_trajs_upto,1);
P_inferred_trajs = cell(num_trajs_upto,1);

for i=1:num_trajs_upto
    WS = WS_cell{i};

    rhs_learned = WS.get_rhs('w',cell2mat(WS.reshape_w));
    tol_dd = 1e-12;
    rhs_bounded = @(x) rhs_learned(max(0, min(1, x)));

    U_inferred_trajs{i} = [];
    E_inferred_trajs{i} = [];
    P_inferred_trajs{i} = [];

    Uobj = Uobj_cell{i};
    for q=1:i
        x0_reduced = Uobj(q).get_x0([]);
        options_ode_sim = odeset('RelTol',tol_dd,'AbsTol',tol_dd*ones(1,nstates),'NonNegative', 1:length(x0_reduced));
        [~,xH0_learned]=ode15s(@(t,x)rhs_bounded(x),time_grid,x0_reduced,options_ode_sim);

        U_inferred_trajs{i} = [U_inferred_trajs{i}; xH0_learned(:,1)'];
        E_inferred_trajs{i} = [E_inferred_trajs{i}; xH0_learned(:,2)'];
        P_inferred_trajs{i} = [P_inferred_trajs{i}; xH0_learned(:,3)'];
    end
end

%% -------- Pack & save payload --------
payload = struct();
payload.meta = struct( ...
    'N',N,'k',k,'network_id',network_id, ...
    'beta',beta,'theta',theta,'eta',eta,'gamma_i',gamma_i,'gamma_p',gamma_p, ...
    't_max',t_max,'t_end',t_end,'NumIteration',NumIteration,'NumGrids',NumGrids, ...
    'noise_ratio',noise_ratio,'rept',rept, ...
    'Initial_E_percents',Initial_E_percents,'start_val',start_val,'end_val',end_val, ...
    'J',J,'ICs',ICs,'num_ICs',num_ICs,'num_trajs_upto',num_trajs_upto );

payload.t_grid = t_grid_saved;
payload.time_grid = time_grid;
payload.S = S;

% True dynamics containers needed for selecting U/E/P mean+envelope later
payload.X6_cell  = X6_cell;
payload.X6_mins  = X6_mins;
payload.X6_maxes = X6_maxes;

% Inferred trajectories (cell by num_traj)
payload.U_inferred_trajs = U_inferred_trajs;
payload.E_inferred_trajs = E_inferred_trajs;
payload.P_inferred_trajs = P_inferred_trajs;

% Mean-field solutions (indexed by IC index i=1:num_ICs)
payload.Usol_mf = Usol_mf;
payload.Esol_mf = Esol_mf;
payload.Psol_mf = Psol_mf;

% Also keep learned file / grid for reference
payload.t_grid_learned = t_grid_learned;

% Save location (DO NOT overwrite)
save_dir = fullfile(cur_dir,'TrajectoryComparison_payloads');
if ~exist(save_dir,'dir'), mkdir(save_dir); end

save_name = sprintf('TrajComparePayload_N=%u_k=%u_tend=%.3g_numTrajUpTo=%u_noise=%.5g_NetID=%u.mat', ...
    N,k,t_end,num_trajs_upto,noise_ratio,network_id);
save_full = fullfile(save_dir, save_name);

if exist(save_full,'file')
    error('Will not overwrite existing payload file:\n%s\nRename save_name or delete manually.', save_full);
end

save(save_full,'payload','-v7.3');
fprintf('\nSaved payload:\n%s\n', save_full);

toc
