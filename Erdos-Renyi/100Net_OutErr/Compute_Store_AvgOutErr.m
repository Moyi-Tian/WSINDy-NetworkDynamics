%% Preprocess one-k raw 100-network output-error files into reusable matrices
% This script loads files such as
% TrajectoryErrValues_100Networks_OneK_N=1000_k=4_traj=1-5_test=10.mat
% for k = 4,...,10, combines them into k-by-trajectory matrices,
% and saves summary matrices plus raw values for later plotting.

clear; close all; tic

%% Set Directory
scriptFullName = matlab.desktop.editor.getActiveFilename();
scriptDir = fileparts(scriptFullName);
if ~isempty(scriptDir), cd(scriptDir); end
cur_dir = pwd;

%% User settings
N = 1000;
k_values = 4:10;
num_trajs_upto = 5;
select_id_num = 10;

in_dir = fullfile(cur_dir, 'WSINDy_OutputErrors_100Networks');
if ~exist(in_dir, 'dir')
    error('Input folder not found:\n%s', in_dir);
end

out_dir = fullfile(cur_dir, 'TrajectoryErr_matrices_100Networks');
if ~exist(out_dir, 'dir')
    mkdir(out_dir);
end

%% Preallocate
num_k = numel(k_values);
state_labels = {'U','E','P'};
num_states = numel(state_labels);

TrajectoryErr_values_all = cell(1, num_states);
TrajectoryErr_mean_all = cell(1, num_states);
TrajectoryErr_median_all = cell(1, num_states);
TrajectoryErr_q25_all = cell(1, num_states);
TrajectoryErr_q75_all = cell(1, num_states);
TrajectoryErr_iqr_all = cell(1, num_states);
TrajectoryErr_numfinite_all = cell(1, num_states);
TrajectoryErr_numslots_all = cell(1, num_states);
TrajectoryErr_nanfrac_all = cell(1, num_states);

for s = 1:num_states
    TrajectoryErr_values_all{s} = cell(num_k, num_trajs_upto);
    TrajectoryErr_mean_all{s} = nan(num_k, num_trajs_upto);
    TrajectoryErr_median_all{s} = nan(num_k, num_trajs_upto);
    TrajectoryErr_q25_all{s} = nan(num_k, num_trajs_upto);
    TrajectoryErr_q75_all{s} = nan(num_k, num_trajs_upto);
    TrajectoryErr_iqr_all{s} = nan(num_k, num_trajs_upto);
    TrajectoryErr_numfinite_all{s} = nan(num_k, num_trajs_upto);
    TrajectoryErr_numslots_all{s} = nan(num_k, num_trajs_upto);
    TrajectoryErr_nanfrac_all{s} = nan(num_k, num_trajs_upto);
end

source_files = cell(num_k, 1);
source_meta = cell(num_k, 1);

%% Load and process each k file
for kk = 1:num_k
    k = k_values(kk);

    pattern = sprintf('TrajectoryErrValues_100Networks_OneK_N=%u_k=%u_traj=1-%u_test=%u.mat', ...
        N, k, num_trajs_upto, select_id_num);

    file_list = dir(fullfile(in_dir, pattern));

    if isempty(file_list)
        broad_pattern = sprintf('TrajectoryErrValues_100Networks_OneK_N=%u_k=%u_*.mat', N, k);
        file_list = dir(fullfile(in_dir, broad_pattern));
    end

    if isempty(file_list)
        warning('No file found for k = %u. Leaving this row as NaN.', k);
        continue;
    end

    [~, idxNewest] = max([file_list.datenum]);
    fname_full = fullfile(file_list(idxNewest).folder, file_list(idxNewest).name);

    fprintf('Processing k = %u from file:\n%s\n', k, file_list(idxNewest).name);

    S = load(fname_full);

    source_files{kk} = fname_full;

    if isfield(S, 'meta')
        source_meta{kk} = S.meta;
    end

    if isfield(S, 'state_labels')
        state_labels = S.state_labels;
    end

    for s = 1:num_states
        for j = 1:num_trajs_upto

            vals = [];

            if isfield(S, 'TrajectoryErr_values_all')
                try
                    vals = S.TrajectoryErr_values_all{s}{1,j};
                catch
                    vals = [];
                end
            end

            if isempty(vals) && isfield(S, 'TrajectoryErr_values_by_traj')
                try
                    vals = S.TrajectoryErr_values_by_traj{j}(:,s);
                catch
                    vals = [];
                end
            end

            if isempty(vals)
                TrajectoryErr_values_all{s}{kk,j} = [];
                TrajectoryErr_numslots_all{s}(kk,j) = 0;
                TrajectoryErr_numfinite_all{s}(kk,j) = 0;
                TrajectoryErr_nanfrac_all{s}(kk,j) = NaN;
                continue;
            end

            vals = vals(:);
            TrajectoryErr_values_all{s}{kk,j} = vals;

            num_slots = numel(vals);
            finite_vals = vals(isfinite(vals));
            num_finite = numel(finite_vals);

            TrajectoryErr_numslots_all{s}(kk,j) = num_slots;
            TrajectoryErr_numfinite_all{s}(kk,j) = num_finite;
            TrajectoryErr_nanfrac_all{s}(kk,j) = 1 - num_finite / num_slots;

            if num_finite == 0
                continue;
            end

            TrajectoryErr_mean_all{s}(kk,j) = mean(finite_vals, 'omitnan');
            TrajectoryErr_median_all{s}(kk,j) = median(finite_vals, 'omitnan');

            q25 = prctile(finite_vals, 25);
            q75 = prctile(finite_vals, 75);

            TrajectoryErr_q25_all{s}(kk,j) = q25;
            TrajectoryErr_q75_all{s}(kk,j) = q75;
            TrajectoryErr_iqr_all{s}(kk,j) = q75 - q25;
        end
    end
end

%% Metadata
summary_meta = struct();
summary_meta.N = N;
summary_meta.k_values = k_values;
summary_meta.num_trajs_upto = num_trajs_upto;
summary_meta.select_id_num = select_id_num;
summary_meta.state_labels = state_labels;
summary_meta.source_files = source_files;
summary_meta.source_meta = source_meta;
summary_meta.created_by = mfilename;
summary_meta.summary_note = ['Summary matrices computed from one-k raw files. ', ...
    'TrajectoryErr_values_all stores raw output-error values for boxplots. ', ...
    'Summary matrices use omitnan.'];

%% Save
out_fname = sprintf('TrajectoryErrSummary_100Networks_N=%u_k=%u-%u_traj=1-%u.mat', ...
    N, k_values(1), k_values(end), num_trajs_upto);

out_fullpath = fullfile(out_dir, out_fname);

save(out_fullpath, ...
    'TrajectoryErr_values_all', ...
    'TrajectoryErr_mean_all', ...
    'TrajectoryErr_median_all', ...
    'TrajectoryErr_q25_all', ...
    'TrajectoryErr_q75_all', ...
    'TrajectoryErr_iqr_all', ...
    'TrajectoryErr_numfinite_all', ...
    'TrajectoryErr_numslots_all', ...
    'TrajectoryErr_nanfrac_all', ...
    'state_labels', ...
    'summary_meta', ...
    '-v7.3');

fprintf('\nSaved preprocessed output-error summary to:\n%s\n', out_fullpath);

toc