% Boxplot per equation: weights grouped by term_idx, with LaTeX tick labels

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

%% LaTeX
set(groot,'defaultTextInterpreter','latex');
set(groot,'defaultAxesTickLabelInterpreter','latex');

%% Optional save
save_figs = 1;
out_dir = fullfile(cur_dir,'figs_for_paper');
if save_figs && ~exist(out_dir,'dir'), mkdir(out_dir); end

%% Term boxes to draw (canonical form, see canonical_term() below)
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
    labels = cell(n_terms,1);
    keys   = cell(n_terms,1);   % canonical keys for boxing + sorting
    feats  = nan(n_terms, 5);   % [group, deg, u, e, p]

    % Build counts + labels + keys + features
    for i = 1:n_terms
        tid = term_ids(i);
        mask_tid = (Teq.term_idx == tid);

        counts(i) = sum(mask_tid);

        % representative term_str for this tid
        ts = Teq.term_str(mask_tid);
        if iscell(ts), ts = ts{1}; else, ts = ts(1); end
        ts = char(ts);

        % LaTeX-friendly
        ts = strrep(ts, '\\cdot', '\cdot');
        labels{i} = ['$' ts '$'];

        % canonical key for matching/ordering
        keys{i} = canonical_term(ts);

        % features for ordering
        feat = term_features(keys{i});     % [group, deg, u, e, p, 0, 0]
        feats(i,:) = feat(1:5);
    end

    % -------- ORDERING (polynomial-type + U/E/P priority) --------
    group = feats(:,1);
    uexp  = feats(:,3);
    eexp  = feats(:,4);
    pexp  = feats(:,5);

    % Sort:
    % 1) group (linear -> bilinear -> quadratic -> cubic -> other)
    % 2) within group, prioritize U then E then P (higher exponent first)
    % 3) stable tie-breaker: canonical key
    [~, ord] = sortrows([group, -uexp, -eexp, -pexp, (1:n_terms)']); %#ok<ASGLU>

    term_ids = term_ids(ord);
    labels   = labels(ord);
    keys     = keys(ord);
    counts   = counts(ord); %#ok<NASGU>

    % Grouping vector following this NEW term order
    [is_in, loc] = ismember(Teq.term_idx, term_ids);
    w_plot   = w_all(is_in);
    group_id = loc(is_in);   % numeric groups 1..n_terms

    %% ---- Box plot ----
    f = figure('Color','w','Units','inches','Position',[1 1 20 5]);
    ax = gca;

    boxplot(w_plot, group_id, 'Symbol','', 'Whisker', Inf, 'BoxStyle','outline');

    % Make all boxplot lines consistent (removes grey “shadow” look)
    L = findall(ax, 'Type', 'Line');
    set(L, 'LineStyle','-', 'Color','k', 'LineWidth',1.5);

    % Emphasis
    set(findall(ax,'Tag','Median'), 'Color', 'r', 'LineWidth', 2.5);
    set(findall(ax,'Tag','Box'),    'Color', 'b', 'LineWidth', 1.8);

    % Slightly thinner whiskers
    set(findall(ax,'Tag','Whisker'),               'LineWidth',1.2);
    set(findall(ax,'Tag','Upper Adjacent Value'),  'LineWidth',1.2);
    set(findall(ax,'Tag','Lower Adjacent Value'),  'LineWidth',1.2);

    % Axis look
    grid on; box on
    ax.FontSize = 20;
    ylabel('Weight', 'FontSize',25);

    % ---- Custom tick labels so we can box them ----
    ax.XTick = 1:n_terms;
    ax.XTickLabel = repmat({''}, size(1:n_terms)); % hide default labels
    ax.TickLabelInterpreter = 'latex';
    ax.XAxis.FontSize = 10;
    ax.XLim = [0.5, n_terms+0.5];

    % % Make room for rotated labels
    % ax.Position(2) = ax.Position(2) + 0.08;
    % ax.Position(4) = ax.Position(4) - 0.08;
    % [left bottom width height]
    ax.Position = [ ...
        0.08 ...   % left   (increase if y-label is big)
        0.28 ...   % bottom (increase for rotated x-labels)
        0.90 ...   % width
        0.62 ...   % height (reduced but safe)
    ];

    % Place labels slightly below bottom
    yl = ylim(ax);
    y_text = yl(1) - 0.06*(yl(2)-yl(1));

    % which terms to box for this equation
    if eq >= 1 && eq <= 3
        box_set = box_terms{eq};
    else
        box_set = {};
    end
    box_set = cellfun(@(s) canonical_term(s), box_set, 'UniformOutput', false);

    edge_col = [0.55 0.10 0.15];  % dark reddish-brown
    lw = 2;

    for ii = 1:n_terms
        do_box = any(strcmp(keys{ii}, box_set));

        if do_box
            text(ax, ii, y_text, labels{ii}, ...
                'Interpreter','latex', ...
                'Rotation',45, ...
                'HorizontalAlignment','right', ...
                'VerticalAlignment','top', ...
                'FontSize',25, ...
                'Margin',4, ...
                'EdgeColor', edge_col, ...
                'LineWidth', lw, ...
                'BackgroundColor','none');
        else
            text(ax, ii, y_text, labels{ii}, ...
                'Interpreter','latex', ...
                'Rotation',45, ...
                'HorizontalAlignment','right', ...
                'VerticalAlignment','top', ...
                'FontSize',25);
        end
    end

    % keep plot elements visible even with below-axis text
    ax.Clipping = 'off';

    if save_figs
        pdf_path = fullfile(out_dir, sprintf('Eq%d_term_weight_boxplot_ordered_boxed.pdf', eq));
        try
            exportgraphics(f, pdf_path, 'ContentType','vector','BackgroundColor','none');
        catch
            set(f,'PaperPositionMode','auto');
            print(f, pdf_path, '-dpdf','-painters');
        end
        fprintf('Saved %s\n', pdf_path);
    end
end

%% ---------------- local helpers ----------------
function key = canonical_term(ts)
    % Canonical form for matching + sorting:
    % - remove spaces
    % - convert \cdot to *
    % - strip $ if any
    ts = strrep(ts, '\\cdot', '\cdot');
    ts = strrep(ts, '\cdot', '*');
    ts = strrep(ts, ' ', '');
    ts = strrep(ts, '$', '');
    key = ts;
end

function feat = term_features(key)
    % Returns feature vector:
    % [group, deg, u, e, p, 0, 0]
    %
    % group:
    % 1 linear, 2 bilinear, 3 quadratic, 4 cubic, 5 other
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
            group = 5;
        end
    elseif deg == 3
        group = 4; % cubic
    else
        group = 5; % other
    end

    feat = [group, deg, u, e, p, 0, 0];
end

function [u,e,p] = parse_exponents(key)
    u = 0; e = 0; p = 0;
    u = u + extract_exp(key, 'U');
    e = e + extract_exp(key, 'E');
    p = p + extract_exp(key, 'P');
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
