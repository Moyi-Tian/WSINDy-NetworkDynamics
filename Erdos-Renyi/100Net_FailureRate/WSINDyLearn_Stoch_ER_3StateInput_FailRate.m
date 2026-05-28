% choose varying number of trajectories with different ICs 
% and run with different levels of noise
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
k = 5; % change k

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
num_trajs_upto = 10;

noise_ratio = 0;

rept = 1;

% Select unseen trajectories
select_id_num = 10;     % number of unseen trajectories to test
S = progressive_ic_selection(Initial_E_percents(J), num_trajs_upto);
Unseen_ids = setdiff(1:30, S{end,1});
select_ids_raw = progressive_ic_selection(Unseen_ids, select_id_num);
select_ids = Unseen_ids(select_ids_raw{end,1});

% Candidate k values to loop over

fname = sprintf('Dynamics from Gillepsie batched E0 on %u Fixed Erdos-Renyi N=%u, k=%u, init_E=%.5g-%.5g, num_IC=%u, beta=%.3g, theta=%.3g, eta=%.3g, gamma_i=%.3g, gamma_p=%.3g, tmax=%u, tend=%.3g, num_t=%u, iter=%u_ThreshT, numNetworkRepeat=%u.mat',NumNetworkRepeat,N,k,Initial_E_percents(1),Initial_E_percents(end),num_ICs,beta,theta,eta,gamma_i,gamma_p,t_max,t_end,NumGrids,NumIteration,NumNetworkRepeat);
load(fname);

save_result = 1; %=1 save; =0 not save

%% Pre-allocate results

fail_rate1 = zeros(num_trajs_upto,1);
fail_rate2 = zeros(num_trajs_upto,1);
fail_rate3 = zeros(num_trajs_upto,1);

S = progressive_ic_selection(ICs, num_trajs_upto);

StartClock_A=tic;
for i=1:num_trajs_upto
    fail1_ct = 0;
    fail2_ct = 0;
    fail3_ct = 0;

    total_ct = 0;

    if mod(i,2) == 0
        fprintf("Running i = %i out of %i ... \n", i, num_trajs_upto)
        EndClock_A=toc(StartClock_A);
        fprintf('Incremental Runtime up to Now: %f  seconds \n', EndClock_A)
        StartClock_A=tic;
    end

    % map back to indices in E0 and values
    idx_sel  = J(S{i});

    for net_r = 1:NumNetworkRepeat

        X_data_6 = X6_cell_cell{net_r};

        X_data_6_select = {X_data_6{idx_sel}};

        X_data = cellfun(@(X) X(:, [1 2 5]), X_data_6_select, 'UniformOutput', false);
        t_data = t_grid_cell{net_r};
    
        for r=1:rept
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
            
            %% simulate learned and true reduced systems
            rhs_learned = WS.get_rhs('w', cell2mat(WS.reshape_w));
            tol_dd = 1e-12;
            rhs_unbounded = @(x) rhs_learned(x);

            for test_r = 1:select_id_num
                X6_test_mat = X_data_6{select_ids(test_r)};
                    
                x0_reduced = [X6_test_mat(1,1), ...
                              X6_test_mat(1,2), ...
                              X6_test_mat(1,5)];

                had_exception = false;
                try
                    [~, xH0_learned] = ode15s(@(t,x) rhs_unbounded(x), ...
                        t_data, x0_reduced, ...
                        odeset('RelTol', tol_dd, 'AbsTol', tol_dd*ones(1,3), 'NonNegative', 1:3));     
                catch
                    had_exception = true;
                    xH0_learned = nan(length(t_data), 3);
                end
        
                % failure diagnostics
                T_true = size(X6_test_mat,1);
                T_learned = size(xH0_learned,1);

                has_nan   = any(~isfinite(xH0_learned(:)));
                bad_len   = (T_learned ~= T_true);
                out_range = any(xH0_learned(:) < -1e-12) || any(xH0_learned(:) > 1+1e-12);
    
                % flat trajectory check — if any of the three states is constant
                state_var = max(xH0_learned) - min(xH0_learned);   % 1×3 vector
                flat_sol  = any(state_var < 1e-12);

                is_fail1 = had_exception || has_nan || bad_len;
                if is_fail1
                    fail1_ct = fail1_ct + 1;
                end

                is_fail2 = is_fail1 || out_range;
                if is_fail2
                    fail2_ct = fail2_ct + 1;
                end

                is_fail3 = is_fail2 || flat_sol;
                if is_fail3
                    fail3_ct = fail3_ct + 1;
                end

                total_ct = total_ct + 1;
            end
            
        end
    end

    fail_rate1(i) = fail1_ct/total_ct;
    fail_rate2(i) = fail2_ct/total_ct;
    fail_rate3(i) = fail3_ct/total_ct;
end

%% Save data
if save_result == 1
    save_path=[cur_dir,'/FailRate_data'];
    if ~exist(save_path, 'dir')
        mkdir(save_path);
    end
    
    fname = sprintf('FailureRate Gillepsie batched E0 on %u Fixed Erdos-Renyi N=%u, k=%u, numTrajUpTo=%u, E0_=%.5g-%.5g, beta=%.3g, theta=%.3g, eta=%.3g, gamma_i=%.3g, gamma_p=%.3g, t_max=%u, t_end=%.3g, t_grid_num=%u, noise=%.5g, repeat=%u, iter=%u_MSTLS_3rdLib.mat',NumNetworkRepeat,N,k,num_trajs_upto,start_val,end_val,beta,theta,eta,gamma_i,gamma_p,t_max,t_end,NumGrids,noise_ratio,rept,NumIteration);
    save(fullfile(save_path, fname),'fail_rate1','fail_rate2','fail_rate3');
end

toc