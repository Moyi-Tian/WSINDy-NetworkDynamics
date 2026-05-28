% Print per-equation term weight summaries (ordered like boxed boxplot)
% No plots. Just console output.

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

%% Term boxes (same as your plotting script)
box_terms = cell(3,1);
box_terms{1} = {'U*E','U*P'};          % equation 1
box_terms{2} = {'E','U*E','U*P'};      % equation 2
box_terms{3} = {'E','P'};              % equation 3

%% Loop equations
eq_list = sort(unique(T.eq(:)).');

for eq = eq_list
    Teq = T(T.eq == eq, :);
    if isempty(Teq), continue; end

    % Ensure weight column is numeric vector
    w_all = Teq.w;
    if iscell(w_all), w_all = cell2mat(w_all); end
    w_all = double(w_all);

    term_ids = unique(Teq.term_idx);
    n_terms  = numel(term_ids);

    counts = zeros(n_terms,1);
    labels_plain = cell(n_terms,1);   % print-friendly
    keys   = cell(n_terms,1);         % canonical keys for boxing + sorting
    feats  = nan(n_terms, 5);         % [group, deg, u, e, p]

    % Build counts + labels + keys + features
    for i = 1:n_terms
        tid = term_ids(i);
        mask_tid = (Teq.term_idx == tid);

        counts(i) = sum(mask_tid);

        ts = Teq.term_str(mask_tid);
        if iscell(ts), ts = ts{1}; else, ts = ts(1); end
        ts = char(ts);

        ts = strrep(ts, '\\cdot', '\cdot');  % LaTeX-friendly display
        labels_plain{i} = ts;

        keys{i} = canonical_term(ts);

        feat = term_features(keys{i});  % [group, deg, u, e, p, 0, 0]
        feats(i,:) = feat(1:5);
    end

    % -------- ORDERING (polynomial-type + U/E/P priority) --------
    group = feats(:,1);
    uexp  = feats(:,3);
    eexp  = feats(:,4);
    pexp  = feats(:,5);

    [~, ord] = sortrows([group, -uexp, -eexp, -pexp, (1:n_terms)']); %#ok<ASGLU>

    term_ids     = term_ids(ord);
    labels_plain = labels_plain(ord);
    keys         = keys(ord);
    counts       = counts(ord);

    % boxed set for this equation
    if eq >= 1 && eq <= 3
        box_set = box_terms{eq};
    else
        box_set = {};
    end
    box_set = cellfun(@(s) canonical_term(s), box_set, 'UniformOutput', false);

    %% ---- Print header ----
    fprintf('\n============================================================\n');
    fprintf('Equation %d: term weight summaries (ordered)\n', eq);
    fprintf('------------------------------------------------------------\n');
    fprintf('%-4s %-25s %-6s %-6s %-10s %-10s %-10s %-10s %-10s %-10s %-10s %-8s\n', ...
        '#','term','n','boxed','mean','median','std','min','q25','q75','max','nOut');
    fprintf('%s\n', repmat('-',1,120));

    %% ---- Per-term stats in this order ----
    for i = 1:n_terms
        tid = term_ids(i);

        mask_tid = (Teq.term_idx == tid);
        w = w_all(mask_tid);
        w = w(isfinite(w));   % safety

        n = numel(w);
        if n == 0
            continue;
        end

        mu  = mean(w);
        med = median(w);
        sd  = std(w);
        mn  = min(w);
        mx  = max(w);
        q25 = prctile(w,25);
        q75 = prctile(w,75);

        % boxplot-style outliers (1.5*IQR fences)
        iqrw = q75 - q25;
        lo_fence = q25 - 1.5*iqrw;
        hi_fence = q75 + 1.5*iqrw;
        nOut = sum(w < lo_fence | w > hi_fence);

        is_boxed = any(strcmp(keys{i}, box_set));
        boxed_str = "";
        if is_boxed, boxed_str = "YES"; end

        fprintf('%-4d %-25s %-6d %-6s %-10.3e %-10.3e %-10.3e %-10.3e %-10.3e %-10.3e %-10.3e %-8d\n', ...
            i, labels_plain{i}, n, boxed_str, mu, med, sd, mn, q25, q75, mx, nOut);
    end

    fprintf('%s\n', repmat('-',1,120));
    fprintf('Total rows in Teq: %d\n', height(Teq));
    fprintf('Sum of per-term n: %d\n', sum(counts));  % should match height(Teq)
end

%% ---------------- local helpers  ----------------
function key = canonical_term(ts)
    ts = strrep(ts, '\\cdot', '\cdot');
    ts = strrep(ts, '\cdot', '*');
    ts = strrep(ts, ' ', '');
    ts = strrep(ts, '$', '');
    key = ts;
end

function feat = term_features(key)
    [u,e,p] = parse_exponents(key);
    deg = u + e + p;

    if deg == 1
        group = 1; % linear
    elseif deg == 2
        if (u==2 && e==0 && p==0) || (u==0 && e==2 && p==0) || (u==0 && e==0 && p==2)
            group = 3; % quadratic (square)
        elseif (u<=1 && e<=1 && p<=1) && (u+e+p==2)
            group = 2; % bilinear
        else
            group = 5; % other
        end
    elseif deg == 3
        group = 4; % cubic
    else
        group = 5; % other
    end

    feat = [group, deg, u, e, p, 0, 0];
end

function [u,e,p] = parse_exponents(key)
    u = extract_exp(key, 'U');
    e = extract_exp(key, 'E');
    p = extract_exp(key, 'P');
end

function expn = extract_exp(s, symb)
    expn = 0;

    tok = regexp(s, [symb '\^(\d+)'], 'tokens');
    if ~isempty(tok)
        for i = 1:numel(tok)
            expn = expn + str2double(tok{i}{1});
        end
        s = regexprep(s, [symb '\^\d+'], '');
    end

    expn = expn + numel(regexp(s, symb, 'match'));
end
