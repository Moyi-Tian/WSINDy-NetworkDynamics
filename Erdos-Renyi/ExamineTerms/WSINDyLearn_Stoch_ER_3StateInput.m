% choose number of trajectories with different ICs 

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

%% Add path and load data
cur_dir=pwd;
addpath(genpath(cur_dir),'provide_path_to_WSINDy',...
    './Generate_Data/data')

N = 1000;
k = 5;

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
NumNetworkRepeat = 100;

start_val = 0.01;
end_val = 0.3;
J = find(Initial_E_percents >= start_val & Initial_E_percents <= end_val);
ICs = Initial_E_percents(J);
num_ICs = numel(J);
num_trajs = 10;

noise_ratio = 0;

% Select unseen trajectories
select_id_num = 10;     % number of unseen trajectories to test
S = progressive_ic_selection(Initial_E_percents(J), num_trajs);
Unseen_ids = setdiff(1:30, S{end,1});
select_ids_raw = progressive_ic_selection(Unseen_ids, select_id_num);
select_ids = Unseen_ids(select_ids_raw{end,1});

% Candidate k values to loop over

fname = sprintf('Dynamics from Gillepsie batched E0 on %u Fixed Erdos-Renyi N=%u, k=%u, init_E=%.5g-%.5g, num_IC=%u, beta=%.3g, theta=%.3g, eta=%.3g, gamma_i=%.3g, gamma_p=%.3g, tmax=%u, tend=%.3g, num_t=%u, iter=%u_ThreshT, numNetworkRepeat=%u.mat',NumNetworkRepeat,N,k,Initial_E_percents(1),Initial_E_percents(end),num_ICs,beta,theta,eta,gamma_i,gamma_p,t_max,t_end,NumGrids,NumIteration,NumNetworkRepeat);
load(fname);

save_result = 1; %=1 save; =0 not save

%% Pre-allocate results

S = progressive_ic_selection(ICs, num_trajs);

% map back to indices in E0 and values
idx_sel  = J(S{num_trajs});

% Storage: one long table-like cell buffer ===
% Each row: [net_r, r, eq, E0, term_idx, w_scaled, w_unscaled, term_str]
LearnedRows = {};     % will convert to table at end
row_ctr = 0;

StartClock_A=tic;
for net_r = 1:NumNetworkRepeat
    if mod(net_r,10) == 0
        fprintf("Running i = %i out of %i ... \n", net_r, NumNetworkRepeat)
        EndClock_A=toc(StartClock_A);
        fprintf('Incremental Runtime up to Now: %f  seconds \n', EndClock_A)
        StartClock_A=tic;
    end

    X_data_6 = X6_cell_cell{net_r};

    X_data_6_select = {X_data_6{idx_sel}};

    X_data = cellfun(@(X) X(:, [1 2 5]), X_data_6_select, 'UniformOutput', false);
    t_data = t_grid_cell{net_r};

    %% get wsindy_data object
    xcell = X_data;
    
    nx = ones(1,size(xcell{1},2));
    xred = cellfun(@(x)x.*nx,xcell,'uni',0);
    ntraj = length(xred);
    tred = t_data;
    Uobj = arrayfun(@(i)wsindy_data(xred{i},tred(:)),(1:ntraj)');
    
    % if many time-series
    arrayfun(@(U) U.coarsen(ceil(U.dims/500)), Uobj);
    
    nstates = Uobj.nstates;
    
    rng('shuffle')
    rng_seed = rng().Seed; rng(rng_seed);
    
    % if many time-series
    arrayfun(@(U) U.addnoise(noise_ratio,'seed',rng_seed), Uobj);
    
    % if many time-series
    % Uobj(1).set_scales([],'nrm',inf,'val',2); % if want to scale
    Uobj(1).set_scales(1); % if not scaling
    scales = Uobj(1).scales;
    arrayfun(@(U) U.set_scales(scales), Uobj(2:end));

    %% get lib tags
    polys = 1:3;
    tags = get_tags(polys,[],nstates);
    lib = library('tags',tags);
                
    %% get test function

    tf_meth = 'FFT'; tf_param = 1;
    p = 12;
    phifun = @(t)(1-t.^2).^p;     
    tf = arrayfun(@(U) arrayfun(@(j) testfcn(U,'phifuns',phifun,'meth',tf_meth,'param',tf_param,'subinds',-4,'stateind',j),(1:nstates)','uni',0),Uobj(:),'uni',0);
   
    %% build WSINDy linear system
    
    WS = wsindy_model(Uobj,lib,tf);
    
    % rescale the coefficents
    Mscale = arrayfun(@(L)L.get_scales(scales),WS.lib(:),'un',0);
    lhs_scales = cellfun(@(t)t.get_scale(scales),WS.lhsterms(:),'un',0);
    Mscale = cellfun(@(M,L)M/L,Mscale,lhs_scales,'un',0);
    
    %% solve
    
    optm = WS_opt();
    
    % only run weak SINDy
    toggle_wendy = 0;
    
    if toggle_wendy==0
        lambdas = 10.^linspace(-4,0,40);
        toggle_jointthreshold = 2;
        [WS,~,~,~,~] = optm.MSTLS(WS,'lambdas',lambdas,'toggle_jointthresh',toggle_jointthreshold, 'M_diag',Mscale);
    elseif toggle_wendy==1
        [WS,~,~,~,~] = optm.wendy(WS,'maxits',100,'regmeth','MSTLS');
    elseif toggle_wendy==2
        lambdas = 10.^linspace(-4,0,40);
        toggle_jointthreshold = 2;
        [WS,~,~,~,~] = optm.MSTLS(WS,'lambdas',lambdas,'toggle_jointthresh',toggle_jointthreshold, 'M_diag',Mscale);
        [WS,~,~,~,~] = optm.wendy(WS);
    elseif toggle_wendy==3
        [WS,~,~,w_its,~,~,~] = optm.MSTLS_WENDy(WS,'maxits_wendy',2,'lambda',10.^linspace(-4,-1,50),'verbose',1);
        disp(['wendy its at optimal lambda=',num2str(size(w_its,2))])
    end
    
    weights = WS.weights;
    
    N_cols = size(tags, 2);  % get number of columns in tags_mat
    if mod(length(weights), N_cols) ~= 0
        error('Length of weights (%d) is not divisible by N_cols (%d)', length(weights), N_cols);
    end

    %% Extract learned terms + weights (scaled and unscaled) and store
    
    % WS.get_supp and WS.reshape_w are the clean access points (as in your plotting script)
    n_eq = numel(WS.get_supp);
    
    for eq = 1:n_eq
        supp_eq = WS.get_supp{eq};    % active term indices for this equation
        if isempty(supp_eq), continue; end
    
        w_scaled_all = WS.reshape_w{eq};     % scaled vector over all candidate terms
        mdiag = Mscale{eq};                 % scaling diag for this equation
    
        for kk = 1:numel(supp_eq)
            term_idx = supp_eq(kk);
    
            if term_idx > numel(w_scaled_all), continue; end
            w_scaled = w_scaled_all(term_idx);
    
            term_str = tag_to_str_tex(tags(term_idx,:));  % same helper as your plotting code
    
            row_ctr = row_ctr + 1;
            LearnedRows(row_ctr,:) = { ...
                net_r, eq, term_idx, w_scaled, term_str ...
            };
        end
    end
end


%% Save data
if save_result == 1
    save_path=[cur_dir,'/LearnedTerms'];
    if ~exist(save_path, 'dir'), mkdir(save_path); end

    % Convert rows to table and save
    LearnedTermsTable = cell2table(LearnedRows, ...
        'VariableNames', {'network_repeat_id','eq','term_idx','w','term_str'});

    fname = sprintf('LearnedTermsTable_Gillespie_ER_N=%u_k=%u_numTrajUpTo=%u_E0=%.5g-%.5g_beta=%.3g_theta=%.3g_eta=%.3g_gi=%.3g_gp=%.3g_tmax=%u_tend=%.3g_numGrid=%u_noise=%.5g_iter=%u_numNetRepeat=%u.mat', ...
        N,k,num_trajs,start_val,end_val,beta,theta,eta,gamma_i,gamma_p,t_max,t_end,NumGrids,noise_ratio,NumIteration,NumNetworkRepeat);

    save(fullfile(save_path, fname), ...
        'LearnedTermsTable','tags', ...
        'N','k','beta','theta','eta','gamma_i','gamma_p', ...
        't_max','t_end','NumGrids','NumIteration','NumNetworkRepeat','num_trajs','noise_ratio');

end

toc



function out = tag_to_str_tex(v)
    syms_tex = {'U','E','P'};
    parts = {};
    for k = 1:numel(v)
        e = v(k);
        if e == 1
            parts{end+1} = syms_tex{k}; %#ok<AGROW>
        elseif e > 1
            parts{end+1} = sprintf('%s^{%d}', syms_tex{k}, e); %#ok<AGROW>
        end
    end
    if isempty(parts)
        out = '1';
    else
        out = strjoin(parts, ' \\cdot ');
    end
end
