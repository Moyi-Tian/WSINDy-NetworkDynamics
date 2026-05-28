% initialization of nodes in E ensures none is isolated node in each iteration 

clc
clear
close all

tic

%% cd to current directory
scriptFullName = matlab.desktop.editor.getActiveFilename();
scriptDir = fileparts(scriptFullName);
% change to current folder
if ~isempty(scriptDir)
    cd(scriptDir);
end


%% Set parameters
N = 1000;
k = 8;

% create network if not exist, otherwise load the network
tempdir=pwd;

beta = 0.2;
theta = 0.5;
eta = 0.4;
gamma_i = 0.1;
gamma_p = 0.2;

pars = [beta, theta, eta, gamma_i, gamma_p];

t_max = 200;
t_end = 40;

Initial_E_percents = linspace(0.01,0.3,30);
num_ICs = length(Initial_E_percents);

NumIteration = 100;
NumGrids = 1000;

NumNetworkRepeat = 100;

save_data = 1; % = 1 save data; = 0 not save data

%% Preallocate result arrays

t_grid_cell = cell(NumNetworkRepeat,1);
X6_cell_cell = cell(NumNetworkRepeat,1);
X6_mins_cell = cell(NumNetworkRepeat,1);
X6_maxes_cell = cell(NumNetworkRepeat,1);

StartClock_A=tic;
for net_r = 1:NumNetworkRepeat

    if mod(net_r,1) == 0
        fprintf("Running i = %i out of %i ... \n", net_r, NumNetworkRepeat)
        EndClock_A=toc(StartClock_A);
        fprintf('Incremental Runtime up to Now: %f  seconds \n', EndClock_A)
        StartClock_A=tic;
    end

    times_all = cell(num_ICs,NumIteration);
    U_all = cell(num_ICs,NumIteration);
    E_all = cell(num_ICs,NumIteration);
    D_all = cell(num_ICs,NumIteration);
    UP_all = cell(num_ICs,NumIteration);
    P_all = cell(num_ICs,NumIteration);
    R_all = cell(num_ICs,NumIteration);
    
    G = sparse(O2O_DesignNetwork_NodesAveDegK(N,k));
    
    %% Run Gillespie
    for i = 1:num_ICs
    
        Initial_E_percent = Initial_E_percents(i);
    
        iter = 1;
        while iter <= NumIteration
            [times,U,E,D,UP,P,R] = func_run_Gillespie(G,pars,t_max,Initial_E_percent);
            if times(end) >= t_end
                times_all{i,iter} = times;
                U_all{i,iter} = U;
                E_all{i,iter} = E;
                D_all{i,iter} = D;
                UP_all{i,iter} = UP;
                P_all{i,iter} = P;
                R_all{i,iter} = R;
                
                iter = iter + 1;
            end
        end
    end
    
    %% Interpolate all runs onto the same uniform time grid and find the mean
    t_grid = linspace(0,t_end,NumGrids);
    X6_cell = cell(num_ICs, 1);
    X6_mins = cell(num_ICs, 1);
    X6_maxes = cell(num_ICs, 1);
    
    for IC_id = 1:num_ICs
        U_mat = zeros(NumIteration,NumGrids);
        E_mat = zeros(NumIteration,NumGrids);
        D_mat = zeros(NumIteration,NumGrids);
        UP_mat = zeros(NumIteration,NumGrids);
        P_mat = zeros(NumIteration,NumGrids);
        R_mat = zeros(NumIteration,NumGrids);
        
        for iter = 1:NumIteration
            U_mat(iter,:) = interp1(times_all{IC_id,iter},U_all{IC_id,iter},t_grid);
            E_mat(iter,:) = interp1(times_all{IC_id,iter},E_all{IC_id,iter},t_grid);
            D_mat(iter,:) = interp1(times_all{IC_id,iter},D_all{IC_id,iter},t_grid);
            UP_mat(iter,:) = interp1(times_all{IC_id,iter},UP_all{IC_id,iter},t_grid);
            P_mat(iter,:) = interp1(times_all{IC_id,iter},P_all{IC_id,iter},t_grid);
            R_mat(iter,:) = interp1(times_all{IC_id,iter},R_all{IC_id,iter},t_grid);
        end
        
        U_mean = mean(U_mat,1);
        U_max = max(U_mat);
        U_min = min(U_mat);
        
        E_mean = mean(E_mat,1);
        E_max = max(E_mat);
        E_min = min(E_mat);
        
        D_mean = mean(D_mat,1);
        D_max = max(D_mat);
        D_min = min(D_mat);
        
        UP_mean = mean(UP_mat,1);
        UP_max = max(UP_mat);
        UP_min = min(UP_mat);
        
        P_mean = mean(P_mat,1);
        P_max = max(P_mat);
        P_min = min(P_mat);
        
        R_mean = mean(R_mat,1);
        R_max = max(R_mat);
        R_min = min(R_mat);
        
        X6_cell{IC_id} = [U_mean', E_mean', D_mean', UP_mean', P_mean', R_mean'];
    
        X6_mins{IC_id} = [U_min', E_min', D_min', UP_min', P_min', R_min'];
        X6_maxes{IC_id} = [U_max', E_max', D_max', UP_max', P_max', R_max']; 
    end
    
    t_grid_cell{net_r} = t_grid;
    X6_cell_cell{net_r} = X6_cell;
    X6_mins_cell{net_r} = X6_mins;
    X6_maxes_cell{net_r} = X6_maxes;
end

%% Save Data
if save_data == 1

    save_path=[tempdir,'/data'];

    parameters.N = N;
    parameters.k = k;
    parameters.pars = pars;
    parameters.t_max = t_max;
    parameters.t_end = t_end;
    parameters.initial_E_percents = Initial_E_percents;
    parameters.number_iternation = NumIteration;
    
    fname = sprintf('Dynamics from Gillepsie batched E0 on %u Fixed Erdos-Renyi N=%u, k=%u, init_E=%.5g-%.5g, num_IC=%u, beta=%.3g, theta=%.3g, eta=%.3g, gamma_i=%.3g, gamma_p=%.3g, tmax=%u, tend=%.3g, num_t=%u, iter=%u_ThreshT, numNetworkRepeat=%u.mat',NumNetworkRepeat,N,k,Initial_E_percents(1),Initial_E_percents(end),num_ICs,beta,theta,eta,gamma_i,gamma_p,t_max,t_end,NumGrids,NumIteration,NumNetworkRepeat);
    save(fullfile(save_path, fname), 'parameters',...
                                     't_grid_cell', 'X6_cell_cell',...
                                     'X6_mins_cell', 'X6_maxes_cell');

end

toc