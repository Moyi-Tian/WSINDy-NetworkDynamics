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

fname = sprintf('Dynamics from Gillepsie batched E0 on Fixed Erdos-Renyi N=%u, k=%u, init_E=%.5g-%.5g, num_IC=%u, beta=%.3g, theta=%.3g, eta=%.3g, gamma_i=%.3g, gamma_p=%.3g, tmax=%u, tend=%.3g, num_t=%u, iter=%u_ThreshT - NetID-%u.mat',N,k,Initial_E_percents(1),Initial_E_percents(end),len_E_percents,beta,theta,eta,gamma_i,gamma_p,t_max,t_end,NumGrids,NumIteration,network_id);
load(fname);

start_val = 0.01;
end_val = 0.3;
J = find(Initial_E_percents >= start_val & Initial_E_percents <= end_val);
ICs = Initial_E_percents(J);
num_ICs = numel(J);
num_trajs_upto = 10;

noise_ratio = 0;

rept = 1;

save_result = 0; %=1 save; =0 not save
save_figs = 0;
save_log = 0;

%% Pre-allocate results

EquationErr_cell = cell(num_trajs_upto,rept);
radii_cell = cell(num_trajs_upto,rept);
G_CondNum_cell = cell(num_trajs_upto,rept);
trajectory_errors_cell = cell(num_trajs_upto,rept);
Weights_learned_cell = cell(num_trajs_upto,rept);

Uobj_cell = cell(num_trajs_upto,rept);
WS_cell = cell(num_trajs_upto,rept);
Str_mod_cell = cell(num_trajs_upto,rept);

S = progressive_ic_selection(ICs, num_trajs_upto);

StartClock_A=tic;
for i=1:num_trajs_upto
    if mod(i,2) == 0
        fprintf("Running i = %i out of %i ... \n", i, num_trajs_upto)
        EndClock_A=toc(StartClock_A);
        fprintf('Incremental Runtime up to Now: %f  seconds \n', EndClock_A)
        StartClock_A=tic;
    end

    % map back to indices in E0 and values
    idx_sel  = J(S{i});

    X_data = {X3_reduced_cell{idx_sel}};
    t_data = t_grid;

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
        M = Uobj.dims;
        
        rng('shuffle')
        rng_seed = rng().Seed; rng(rng_seed);
        
        % if many time-series
        arrayfun(@(U) U.addnoise(noise_ratio,'seed',rng_seed), Uobj);
        
        % if many time-series
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
        Mscale_W = cell2mat(Mscale);
        
       
        %% solve
        
        optm = WS_opt();
        
        % only run weak SINDy
        toggle_wendy = 0;
        
        if toggle_wendy==0
            lambdas = 10.^linspace(-4,0,40);
            toggle_jointthreshold = 2;
            [WS,loss_wsindy,its,G,b] = optm.MSTLS(WS,'lambdas',lambdas,'toggle_jointthresh',toggle_jointthreshold, 'M_diag',Mscale);
        elseif toggle_wendy==1
            [WS,w_its,res,res_0,CovW] = optm.wendy(WS,'maxits',100,'regmeth','MSTLS');
        elseif toggle_wendy==2
            lambdas = 10.^linspace(-4,0,40);
            toggle_jointthreshold = 2;
            [WS,loss_wsindy,its,G,b] = optm.MSTLS(WS,'lambdas',lambdas,'toggle_jointthresh',toggle_jointthreshold, 'M_diag',Mscale);
            [WS,w_its,res,res_0,CovW] = optm.wendy(WS);
        elseif toggle_wendy==3
            [WS,loss_wsindy,lambda,w_its,res,res_0,CovW] = optm.MSTLS_WENDy(WS,'maxits_wendy',2,'lambda',10.^linspace(-4,-1,50),'verbose',1);
            disp(['wendy its at optimal lambda=',num2str(size(w_its,2))])
        end
        
        weights = WS.weights;
        
        N_cols = size(tags, 2);  % get number of columns in tags_mat
        if mod(length(weights), N_cols) ~= 0
            error('Length of weights (%d) is not divisible by N_cols (%d)', length(weights), N_cols);
        end
        weights = reshape(weights, [], N_cols);  % reshape
        
        W_nd = cellfun(@(w,m)w./m,WS.reshape_w,Mscale,'un',0);
        
        %% simulate learned and true reduced systems
        
        toggle_compare = 1:ntraj;
        if ~isempty(toggle_compare)
            w_plot = WS.weights;
            rhs_learned = WS.get_rhs('w',cell2mat(W_nd)); % learned dynamics in the scaled coefficients
            tol_dd = 10^-12;
            rhs_bounded = @(x) rhs_learned(max(0, min(1, x)));
            
            trajectory_errors = zeros(ntraj, nstates);
            for q=toggle_compare
                t_train = Uobj(q).grid{1};
                x0_reduced = Uobj(q).get_x0([]);
                options_ode_sim = odeset('RelTol',tol_dd,'AbsTol',tol_dd*ones(1,nstates),'NonNegative', 1:length(x0_reduced));
                [t_learned,xH0_learned]=ode15s(@(t,x)rhs_bounded(x),t_train,x0_reduced,options_ode_sim);
                
                if save_figs == 1
                    figure(q);clf
                end
                for l=1:nstates
                    if size(xH0_learned,1) == size(Uobj(q).Uobs{l},1)
                        trajectory_errors(q,l) = norm(xH0_learned(:,l) - Uobj(q).Uobs{l}) / norm(Uobj(q).Uobs{l});
                    else
                        trajectory_errors(q,l) = NaN;
                        warning('Skipped error computation for traj %d, state %d due to length mismatch.', q, l);
                    end

                    if save_figs == 1
                        subplot(nstates,1,l)
                        plot(Uobj(q).grid{1},Uobj(q).Uobs{l},'b-o',t_learned,xH0_learned(:,l),'r-.','linewidth',2)
                        ylim([0 1])
                        try
                            title(['rel err=',num2str(trajectory_errors(q,l))])
                        catch
                        end
                        legend({'data','learned'})

                        % Save figure
                        fig_folder = sprintf('figures_3StateInput_N=%u_k=%u_numTrajUpTo_%u_IC_E0_=%.5g-%.5g_rept=%u_noise_%.5g_MSTLS_tau=%.5g, theta=%.5g, eta=%.5g, gamma_i=%.5g, gamma_p=%.5g, tend=%.3g, iter=%u_3rdLib - NetID-%u', N,k,num_trajs_upto,start_val,end_val,rept,noise_ratio,beta,theta,eta,gamma_i,gamma_p,t_end,NumIteration,network_id);
                        if ~exist(fullfile(cur_dir,fig_folder), 'dir')
                            mkdir(fullfile(cur_dir,fig_folder));
                        end
                    
                        % Customize filename
                        fig_name = sprintf('traj_compare_numTraj=%u_noise=%.5g_rept%u_traj%u',i, noise_ratio,r,q);
                        full_path_jpg = fullfile(fig_folder, [fig_name, '.jpg']);        
                        % Save
                        saveas(gcf, full_path_jpg);
                    end
                end
            end
        
        end
    
        Str_mod = WS.disp_mod('w',WS.weights);
        resids = cellfun(@(G, w, b) norm(G*w - b)/norm(b), WS.Gs{1}, W_nd, WS.bs{1});
        
        %% Write to file
        if save_log == 1
            % Generate file name using sprintf
            filename = sprintf('Results_3StateInput_N=%u_k=%u_numTrajUpTo_%u_IC_E0_=%.5g-%.5g_rept=%u_noise_%.5g_MSTLS_tau=%.5g, theta=%.5g, eta=%.5g, gamma_i=%.5g, gamma_p=%.5g, tend=%.3g, iter=%u_3rdLib - NetID-%u.txt', N,k,num_trajs_upto,start_val,end_val,rept,noise_ratio,beta,theta,eta,gamma_i,gamma_p,t_end,NumIteration,network_id);
            
            % Full file path
            folder = 'txt_logs';
            if ~exist(fullfile(cur_dir,folder), 'dir')
                mkdir(fullfile(cur_dir,folder));
            end
            filepath = fullfile(folder, filename);
            
            % Open file in append mode (creates if not exist)
            fid = fopen(filepath, 'a');
            if fid == -1
                error('Failed to open log file for writing.');
            end
            
            % write to the file
            fprintf(fid, '\nWith num Trajs = %u, noise=%.5g, rept=%u', i, noise_ratio, r);
            fprintf(fid, '\n------------------------------------------------------------------------\n');
            fprintf(fid, '\ndata dims ='); fprintf(fid, '%u ', Uobj.dims); fprintf(fid, '\n');
            
            for c = 1:WS.numeq
                fprintf(fid, '----------Eq %i----------\n', c);
                cellfun(@(s) fprintf(fid, '%s\n', s), Str_mod{c});
            end
            
            fprintf(fid, '\n');
            arrayfun(@(r) fprintf(fid, 'rel resid=%g\n', r), resids);
            cellfun(@(s) fprintf(fid, 'sparsity (number of terms)=%d\n', length(s)), WS.get_supp);
            fprintf(fid, '\ntf rads='); fprintf(fid, '%u ', WS.tf{1}{1}.rads);
            fprintf(fid, '\nsize G =%i', size(WS.Gs{1}{2}));
            fprintf(fid, '\ncond G =%i \n', cond(WS.Gs{1}{2}));
            fprintf(fid, '\n------------------------------------------------------------------------\n');
            
            fclose(fid);  % close the file
        end

        EquationErr_cell{i,r} = resids;
        radii_cell{i,r} = WS.tf{1}{1}.rads;
        G_CondNum_cell{i,r} = cond(WS.Gs{1}{2});
        trajectory_errors_cell{i,r} = trajectory_errors;
        Weights_learned_cell{i,r} = WS.weights;

        Uobj_cell{i,r} = Uobj;
        WS_cell{i,r} = WS;
        Str_mod_cell{i,r} = Str_mod;
    
    end
end

%% Save data
if save_result == 1
    save_path=[cur_dir,'/WSINDy_data_3StateInput'];
    if ~exist(save_path, 'dir')
        mkdir(save_path);
    end

    fname = sprintf('WSINDy Results for Gillepsie batched E0 on Fixed Erdos-Renyi N=%u, k=%u, numTrajUpTo=%u, E0_=%.5g-%.5g, beta=%.3g, theta=%.3g, eta=%.3g, gamma_i=%.3g, gamma_p=%.3g, t_max=%u, t_end=%.3g, t_grid_num=%u, noise=%.5g, repeat=%u, iter=%u_MSTLS_3rdLib - NetID-%u.mat',N,k,num_trajs_upto,start_val,end_val,beta,theta,eta,gamma_i,gamma_p,t_max,t_end,NumGrids,noise_ratio,rept,NumIteration,network_id);
    save(fullfile(save_path, fname),...
                'Uobj_cell', 'WS_cell', 'Str_mod_cell', ...
                'EquationErr_cell', ...
                'radii_cell','G_CondNum_cell','trajectory_errors_cell',...
                'Weights_learned_cell');
end

toc