%% WSINDy on 100 fixed ER network repeats for one selected k
% For one k:
%   For each network repeat net_r = 1,...,100:
%       For each number of training trajectories i:
%           Learn one WSINDy model from that network-repeat mean trajectory data.
%           Test on 10 unseen initial conditions from the same network repeat.
%
% Saved:
%   TrajectoryErr_values_by_traj{i}: 1000 x 3 matrix
%       columns are U, E, P output errors
%
%   TrajectoryErr_values_all{s}{1,i}: 1000 x 1 vector
%       s = 1: U, s = 2: E, s = 3: P
%
%   TrajectoryErr_mean_all{s}(1,i): mean over 1000 slots, omitting NaNs
%   TrajectoryErr_median_all{s}(1,i): median over 1000 slots, omitting NaNs
%   TrajectoryErr_numfinite_all{s}(1,i): number of finite values

tic
clear
close all

%% cd to current directory
scriptFullName = matlab.desktop.editor.getActiveFilename();
scriptDir = fileparts(scriptFullName);
if ~isempty(scriptDir)
    cd(scriptDir);
end

%% Add path and load data
cur_dir = pwd;

data_path = fullfile(cur_dir, 'Gillespie_data_networkr');
if ~exist(data_path, 'dir')
    error('Data folder not found:\n%s', data_path);
end

addpath(genpath(cur_dir), ...
    '../../../wsindy_obj_base-main_06_25_25/', ...
    data_path)

%% ---------------- User settings ----------------
N = 1000;
k = 9;   % choose one k here

beta = 0.2;
theta = 0.5;
eta = 0.4;
gamma_i = 0.1;
gamma_p = 0.2;

t_max = 200;
t_end = 40;

Initial_E_percents = linspace(0.01, 0.3, 30);
num_ICs = length(Initial_E_percents);

NumIteration = 100;
NumGrids = 1000;
NumNetworkRepeat = 100;

start_val = 0.01;
end_val = 0.3;

J = find(Initial_E_percents >= start_val & Initial_E_percents <= end_val);
ICs = Initial_E_percents(J);

num_trajs_upto = 5;
select_id_num = 10;

noise_ratio = 0;
rept = 1;

save_result = 1;
save_models = 0;   % set to 1 only if you want to save all WS objects

checkpoint_save = 1;
checkpoint_every = 10;

state_labels = {'U','E','P'};

%% Suppress warnings during batch run
old_warning_state = warning('off','all');

%% Load 100-network data file for selected k
fname_pattern = sprintf(['Dynamics from Gillepsie batched E0 on 100 Fixed Erdos-Renyi ', ...
    'N=%u, k=%u, init_E=%.5g-%.5g, num_IC=%u, beta=%.3g, theta=%.3g, ', ...
    'eta=%.3g, gamma_i=%.3g, gamma_p=%.3g, tmax=%u, tend=%.3g, ', ...
    'num_t=%u, iter=%u_ThreshT, numNetworkRepeat=%u.mat'], ...
    N, k, Initial_E_percents(1), Initial_E_percents(end), num_ICs, ...
    beta, theta, eta, gamma_i, gamma_p, t_max, t_end, NumGrids, ...
    NumIteration, NumNetworkRepeat);

fname_full = fullfile(data_path, fname_pattern);

if ~exist(fname_full, 'file')
    file_list = dir(fullfile(data_path, sprintf('*N=%u, k=%u,*numNetworkRepeat=%u.mat', ...
        N, k, NumNetworkRepeat)));

    if isempty(file_list)
        warning(old_warning_state);
        error('No data file found for k=%u in:\n%s', k, data_path);
    end

    [~, idxNewest] = max([file_list.datenum]);
    fname_full = fullfile(file_list(idxNewest).folder, file_list(idxNewest).name);
end

D = load(fname_full, 'parameters', 't_grid_cell', 'X6_cell_cell', ...
    'X6_mins_cell', 'X6_maxes_cell');

if ~isfield(D, 'X6_cell_cell') || ~isfield(D, 't_grid_cell')
    warning(old_warning_state);
    error('Loaded file does not contain X6_cell_cell and t_grid_cell.');
end

X6_cell_cell = D.X6_cell_cell;
t_grid_cell = D.t_grid_cell;

num_net_use = min(NumNetworkRepeat, numel(X6_cell_cell));

%% Select training and test trajectories
S = progressive_ic_selection(ICs, num_trajs_upto);

max_train_ids = J(S{end});
Unseen_ids = setdiff(J, max_train_ids);

S_test = progressive_ic_selection(Initial_E_percents(Unseen_ids), select_id_num);
select_ids = Unseen_ids(S_test{end});

%% Preallocate outputs
total_slots = num_net_use * select_id_num;

TrajectoryErr_values_by_traj = cell(num_trajs_upto, 1);
for i = 1:num_trajs_upto
    TrajectoryErr_values_by_traj{i} = nan(total_slots, 3);
end

TrajectoryErr_values_all = cell(1, 3);
TrajectoryErr_mean_all = cell(1, 3);
TrajectoryErr_median_all = cell(1, 3);
TrajectoryErr_numfinite_all = cell(1, 3);

for s = 1:3
    TrajectoryErr_values_all{s} = cell(1, num_trajs_upto);
    TrajectoryErr_mean_all{s} = nan(1, num_trajs_upto);
    TrajectoryErr_median_all{s} = nan(1, num_trajs_upto);
    TrajectoryErr_numfinite_all{s} = nan(1, num_trajs_upto);
end

if save_models == 1
    WS_cell_all = cell(num_net_use, num_trajs_upto);
    Str_mod_cell_all = cell(num_net_use, num_trajs_upto);
    Weights_learned_cell_all = cell(num_net_use, num_trajs_upto);
end

%% Save path and metadata
save_path = fullfile(cur_dir, 'WSINDy_OutputErrors_100Networks');
if ~exist(save_path, 'dir')
    mkdir(save_path);
end

fname_out = sprintf(['TrajectoryErrValues_100Networks_OneK_N=%u_k=%u_', ...
    'traj=1-%u_test=%u.mat'], ...
    N, k, num_trajs_upto, select_id_num);

out_fullpath = fullfile(save_path, fname_out);

meta = struct();
meta.N = N;
meta.k = k;
meta.beta = beta;
meta.theta = theta;
meta.eta = eta;
meta.gamma_i = gamma_i;
meta.gamma_p = gamma_p;
meta.t_max = t_max;
meta.t_end = t_end;
meta.Initial_E_percents = Initial_E_percents;
meta.start_val = start_val;
meta.end_val = end_val;
meta.NumIteration = NumIteration;
meta.NumGrids = NumGrids;
meta.NumNetworkRepeat = NumNetworkRepeat;
meta.num_net_use = num_net_use;
meta.num_trajs_upto = num_trajs_upto;
meta.select_id_num = select_id_num;
meta.noise_ratio = noise_ratio;
meta.train_selection = S;
meta.max_train_ids = max_train_ids;
meta.select_ids = select_ids;
meta.state_labels = state_labels;
meta.source_data_file = fname_full;
meta.note = ['TrajectoryErr_values_by_traj{i} is a 1000 x 3 matrix for ', ...
             'the selected k, with rows corresponding to network-repeat/test-IC pairs.'];

%% Main loop over network repeats
for net_r = 1:num_net_use

    if mod(net_r, 5) == 0
        fprintf('Running network %u out of %u\n', net_r, num_net_use);
    end

    X6_cell = X6_cell_cell{net_r};
    t_data = t_grid_cell{net_r};

    if isempty(X6_cell) || isempty(t_data)
        continue;
    end

    t_data = t_data(:);

    %% Build reduced 3-state cell from X6 data
    % X6 columns: U, E, D, UP, P, R
    X3_reduced_cell = cell(size(X6_cell));

    for IC_id = 1:numel(X6_cell)
        X = X6_cell{IC_id};
        X3_reduced_cell{IC_id} = [X(:,1), X(:,2), X(:,5)];
    end

    %% Loop over number of training trajectories
    for i = 1:num_trajs_upto

        idx_sel = J(S{i});
        X_data = {X3_reduced_cell{idx_sel}};

        %% get wsindy_data object
        xcell = X_data;

        try
            nx = ones(1, size(xcell{1}, 2));
            xred = cellfun(@(x) x .* nx, xcell, 'uni', 0);
            ntraj = length(xred);
            tred = t_data;

            Uobj = arrayfun(@(q) wsindy_data(xred{q}, tred(:)), (1:ntraj)');
        catch
            continue;
        end

        %% Coarsen data as in old working code
        try
            arrayfun(@(U) U.coarsen(ceil(U.dims/500)), Uobj);
        catch
            continue;
        end

        nstates = Uobj(1).nstates;

        rng('shuffle')
        rng_seed = rng().Seed;
        rng(rng_seed);

        %% Add noise, usually zero here
        try
            arrayfun(@(U) U.addnoise(noise_ratio, 'seed', rng_seed), Uobj);
        catch
            continue;
        end

        %% Scaling
        try
            Uobj(1).set_scales(1);
            scales = Uobj(1).scales;
            arrayfun(@(U) U.set_scales(scales), Uobj(2:end));
        catch
            continue;
        end

        %% Get library tags
        polys = 1:3;
        tags = get_tags(polys, [], nstates);
        lib = library('tags', tags);

        %% Get test functions
        tf_meth = 'FFT';
        tf_param = 1;
        p = 12;
        phifun = @(t) (1 - t.^2).^p;

        try
            tf = arrayfun(@(U) ...
                arrayfun(@(jstate) ...
                    testfcn(U, ...
                    'phifuns', phifun, ...
                    'meth', tf_meth, ...
                    'param', tf_param, ...
                    'subinds', -4, ...
                    'stateind', jstate), ...
                (1:nstates)', 'uni', 0), ...
                Uobj(:), 'uni', 0);
        catch
            continue;
        end

        %% Build WSINDy linear system
        try
            WS = wsindy_model(Uobj, lib, tf);
        catch
            continue;
        end

        %% Rescale coefficients
        try
            Mscale = arrayfun(@(L) L.get_scales(scales), WS.lib(:), 'un', 0);
            lhs_scales = cellfun(@(t) t.get_scale(scales), WS.lhsterms(:), 'un', 0);
            Mscale = cellfun(@(M, L) M/L, Mscale, lhs_scales, 'un', 0);
        catch
            continue;
        end

        %% Solve with MSTLS
        optm = WS_opt();

        try
            lambdas = 10.^linspace(-4, 0, 40);
            toggle_jointthreshold = 2;
            [WS, ~, ~, ~, ~] = optm.MSTLS(WS, ...
                'lambdas', lambdas, ...
                'toggle_jointthresh', toggle_jointthreshold, ...
                'M_diag', Mscale);
        catch
            continue;
        end

        %% Learned RHS, same structure as old working code
        try
            weights = WS.weights;

            N_cols = size(tags, 2);
            if mod(length(weights), N_cols) ~= 0
                continue;
            end
            weights = reshape(weights, [], N_cols); %#ok<NASGU>

            W_nd = cellfun(@(w, m) w ./ m, WS.reshape_w, Mscale, 'un', 0);

            rhs_learned = WS.get_rhs('w', cell2mat(W_nd));
            rhs_unbounded = @(x) rhs_learned(x);
        catch
            continue;
        end

        %% Test on unseen trajectories
        trajectory_errors_test = nan(select_id_num, nstates);
        tol_dd = 1e-12;

        for q = 1:select_id_num

            test_id = select_ids(q);
            X_true = X3_reduced_cell{test_id};

            t_test = t_data;
            x0_reduced = X_true(1, :);

            options_ode_sim = odeset( ...
                'RelTol', tol_dd, ...
                'AbsTol', tol_dd * ones(1, nstates), ...
                'NonNegative', 1:length(x0_reduced));

            try
                [~, xH0_learned] = ode15s(@(t, x) rhs_unbounded(x), ...
                    t_test, x0_reduced, options_ode_sim);
            catch
                continue;
            end

            for l = 1:nstates
                if size(xH0_learned, 1) == size(X_true, 1)
                    denom = norm(X_true(:, l));
                    if denom > 0
                        trajectory_errors_test(q, l) = ...
                            norm(xH0_learned(:, l) - X_true(:, l)) / denom;
                    else
                        trajectory_errors_test(q, l) = ...
                            norm(xH0_learned(:, l) - X_true(:, l));
                    end
                else
                    trajectory_errors_test(q, l) = NaN;
                end
            end
        end

        %% Store test errors into 1000-slot matrix
        row_idx = (net_r - 1) * select_id_num + (1:select_id_num);
        TrajectoryErr_values_by_traj{i}(row_idx, :) = trajectory_errors_test;

        if save_models == 1
            WS_cell_all{net_r, i} = WS;
            Str_mod_cell_all{net_r, i} = WS.disp_mod('w', WS.weights);
            Weights_learned_cell_all{net_r, i} = WS.weights;
        end
    end

    %% Update summaries
    for s = 1:3
        for i = 1:num_trajs_upto
            vals = TrajectoryErr_values_by_traj{i}(:, s);
            TrajectoryErr_values_all{s}{1, i} = vals;
            TrajectoryErr_mean_all{s}(1, i) = mean(vals, 'omitnan');
            TrajectoryErr_median_all{s}(1, i) = median(vals, 'omitnan');
            TrajectoryErr_numfinite_all{s}(1, i) = sum(isfinite(vals));
        end
    end

    %% Checkpoint save
    if save_result == 1 && checkpoint_save == 1 && mod(net_r, checkpoint_every) == 0
        if save_models == 1
            save(out_fullpath, ...
                'TrajectoryErr_values_by_traj', ...
                'TrajectoryErr_values_all', ...
                'TrajectoryErr_mean_all', ...
                'TrajectoryErr_median_all', ...
                'TrajectoryErr_numfinite_all', ...
                'state_labels', ...
                'meta', ...
                'WS_cell_all', ...
                'Str_mod_cell_all', ...
                'Weights_learned_cell_all', ...
                '-v7.3');
        else
            save(out_fullpath, ...
                'TrajectoryErr_values_by_traj', ...
                'TrajectoryErr_values_all', ...
                'TrajectoryErr_mean_all', ...
                'TrajectoryErr_median_all', ...
                'TrajectoryErr_numfinite_all', ...
                'state_labels', ...
                'meta', ...
                '-v7.3');
        end
    end
end

%% Final save
if save_result == 1
    if save_models == 1
        save(out_fullpath, ...
            'TrajectoryErr_values_by_traj', ...
            'TrajectoryErr_values_all', ...
            'TrajectoryErr_mean_all', ...
            'TrajectoryErr_median_all', ...
            'TrajectoryErr_numfinite_all', ...
            'state_labels', ...
            'meta', ...
            'WS_cell_all', ...
            'Str_mod_cell_all', ...
            'Weights_learned_cell_all', ...
            '-v7.3');
    else
        save(out_fullpath, ...
            'TrajectoryErr_values_by_traj', ...
            'TrajectoryErr_values_all', ...
            'TrajectoryErr_mean_all', ...
            'TrajectoryErr_median_all', ...
            'TrajectoryErr_numfinite_all', ...
            'state_labels', ...
            'meta', ...
            '-v7.3');
    end
end

%% Restore warnings
warning(old_warning_state);

toc