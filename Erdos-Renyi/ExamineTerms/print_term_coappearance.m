% Pair co-appearance report:
% For each equation separately, report:
% 1) pairs that ALWAYS co-appear (whenever either appears, the other also appears)
% 2) pairs that NEVER co-appear (they never appear together)
%
% Uses LearnedTermsTable with columns: network_repeat_id, eq, term_idx, w, term_str

clear
close all

%% cd to current directory
scriptFullName = matlab.desktop.editor.getActiveFilename();
scriptDir = fileparts(scriptFullName);
if ~isempty(scriptDir), cd(scriptDir); end

%% Load
cur_dir   = pwd;
save_path = fullfile(cur_dir, 'LearnedTerms');

files = dir(fullfile(save_path, 'LearnedTermsTable_*.mat'));
if isempty(files), error('No LearnedTermsTable_*.mat found in %s', save_path); end
[~, idx_new] = max([files.datenum]);
S = load(fullfile(save_path, files(idx_new).name));
T = S.LearnedTermsTable;

%% Basic cleaning / canonical term key
T.term_key = cell(height(T),1);
for i = 1:height(T)
    ts = T.term_str{i};
    if isstring(ts), ts = char(ts); end
    T.term_key{i} = canonical_term(ts);
end

eq_list = sort(unique(T.eq(:)).');

for eq = eq_list
    Teq = T(T.eq == eq, :);
    if isempty(Teq), continue; end

    % Define what counts as a "learned system instance":
    % If you have multiple trajectories per network_repeat_id, and you want
    % co-appearance at the *network repeat* level, group by network_repeat_id only.
    % If you later add a trajectory id, include it here too.
    sys_ids = unique(Teq.network_repeat_id(:)).';
    n_sys   = numel(sys_ids);

    % Term universe for this equation
    terms = unique(Teq.term_key);
    n_terms = numel(terms);

    % Build presence matrix A (n_sys x n_terms), A(s,t)=1 if term t appears in system s
    A = false(n_sys, n_terms);

    % Build a map from term_key to column index
    term2col = containers.Map(terms, num2cell(1:n_terms));

    % Fill presence
    for s = 1:n_sys
        sid = sys_ids(s);
        rows = Teq.network_repeat_id == sid;

        present_keys = unique(Teq.term_key(rows));
        for kidx = 1:numel(present_keys)
            key = present_keys{kidx};
            if isKey(term2col, key)
                A(s, term2col(key)) = true;
            end
        end
    end

    % Co-appearance counts
    % C(i,j) = number of systems where i and j both present
    C = A' * A;                    % n_terms x n_terms (double)
    N = double(n_sys);

    % Individual presence counts
    pres = sum(A, 1);              % 1 x n_terms

    % Pair always co-appear:
    % interpret as: support(i) == support(j) and >0
    % i.e., they appear in exactly the same set of learned systems.
    always_pairs = [];
    for i = 1:n_terms
        if pres(i) == 0, continue; end
        for j = i+1:n_terms
            if pres(j) == 0, continue; end
            % exact same presence pattern
            if all(A(:,i) == A(:,j))
                always_pairs(end+1,:) = [i j]; %#ok<AGROW>
            end
        end
    end

    % Pair never co-appear:
    % they never appear together in any system, but both appear at least once.
    never_pairs = [];
    for i = 1:n_terms
        if pres(i) == 0, continue; end
        for j = i+1:n_terms
            if pres(j) == 0, continue; end
            if C(i,j) == 0
                never_pairs(end+1,:) = [i j]; %#ok<AGROW>
            end
        end
    end

    % Optional: rank never-pairs by "strength" (both common but never co-occur)
    % score = min(pres(i), pres(j)) or pres(i)*pres(j)
    if ~isempty(never_pairs)
        score = zeros(size(never_pairs,1),1);
        for r = 1:size(never_pairs,1)
            i = never_pairs(r,1); j = never_pairs(r,2);
            score(r) = min(pres(i), pres(j));
        end
        [~,ord] = sort(score, 'descend');
        never_pairs = never_pairs(ord,:);
        score = score(ord);
    end

    % Print
    fprintf('\n============================================================\n');
    fprintf('Equation %d: co-appearance pairs over %d learned systems (network_repeat_id)\n', eq, n_sys);
    fprintf('------------------------------------------------------------\n');

    % Show term prevalence (helpful context)
    fprintf('Term prevalence (top 15):\n');
    [pres_sorted,ordp] = sort(pres,'descend');
    topk = min(15, n_terms);
    for kkk = 1:topk
        tcol = ordp(kkk);
        fprintf('  %-18s  %4d / %d\n', terms{tcol}, pres_sorted(kkk), n_sys);
    end

    fprintf('\n--- ALWAYS co-appear (identical support sets, both present >= 1 time) ---\n');
    if isempty(always_pairs)
        fprintf('  (none)\n');
    else
        for r = 1:size(always_pairs,1)
            i = always_pairs(r,1); j = always_pairs(r,2);
            fprintf('  %-18s <-> %-18s   (appeared together %d / %d systems)\n', ...
                terms{i}, terms{j}, C(i,j), n_sys);
        end
    end

    fprintf('\n--- NEVER co-appear (C(i,j)=0, both present >= 1 time) ---\n');
    if isempty(never_pairs)
        fprintf('  (none)\n');
    else
        max_print = min(50, size(never_pairs,1)); % avoid huge print
        fprintf('  printing top %d by min(prevalence)\n', max_print);
        for r = 1:max_print
            i = never_pairs(r,1); j = never_pairs(r,2);
            fprintf('  %-18s X  %-18s   pres=(%d,%d)  min=%d\n', ...
                terms{i}, terms{j}, pres(i), pres(j), min(pres(i), pres(j)));
        end
        if size(never_pairs,1) > max_print
            fprintf('  ... (%d more never-pairs not shown)\n', size(never_pairs,1)-max_print);
        end
    end

    % Optional: show a compact summary
    fprintf('\nSummary: %d terms, %d always-pairs, %d never-pairs\n', ...
        n_terms, size(always_pairs,1), size(never_pairs,1));
end

%% ---------------- local helpers ----------------
function key = canonical_term(ts)
    ts = strrep(ts, '\\cdot', '\cdot');
    ts = strrep(ts, '\cdot', '*');
    ts = strrep(ts, ' ', '');
    ts = strrep(ts, '$', '');
    key = char(ts);
end