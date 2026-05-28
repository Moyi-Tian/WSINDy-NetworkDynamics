% Histogram per equation: count by term_idx, label using term_str from table
% NEW:
% 1) Order terms by polynomial type: linear -> bilinear -> quadratic -> cubic
% 2) Within each type, order by variable priority U -> E -> P (deterministic)
% 3) Draw a box around selected terms on the x-axis (per equation) via custom tick labels

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

    term_ids = unique(Teq.term_idx);
    n_terms  = numel(term_ids);

    counts = zeros(n_terms,1);
    labels = cell(n_terms,1);
    keys   = cell(n_terms,1);   % canonical keys for boxing + sorting
    feats  = nan(n_terms, 7);   % [group, deg, u, e, p, withinKey1, withinKey2]

    for i = 1:n_terms
        tid = term_ids(i);
        mask_tid = (Teq.term_idx == tid);

        % count occurrences (rows)
        counts(i) = sum(mask_tid);

        % representative term_str
        ts = Teq.term_str(mask_tid);
        if iscell(ts), ts = ts{1}; else, ts = ts(1); end
        ts = char(ts);

        % LaTeX-friendly cleanup
        ts = strrep(ts, '\\cdot', '\cdot'); % if stored with double slashes

        % store label for display
        labels{i} = ['$' ts '$'];

        % canonical key for matching/ordering
        keys{i} = canonical_term(ts);

        % features for ordering
        feats(i,:) = term_features(keys{i});
    end

    % -------- ORDERING (NO LONGER by counts/frequency) --------
    % Primary: group (linear, bilinear, quadratic, cubic, other)
    % Secondary: variable priority U->E->P in a deterministic way via exponents
    % Tertiary: canonical key as a stable tie-breaker
    %
    % feats columns:
    % 1 group, 2 deg, 3 u, 4 e, 5 p, 6/7 reserved (not used here)
    group = feats(:,1);
    uexp  = feats(:,3);
    eexp  = feats(:,4);
    pexp  = feats(:,5);

    % Within-group ordering rule (U then E then P):
    % sort by -uexp, then -eexp, then -pexp, then key
    [~, ord] = sortrows([group, -uexp, -eexp, -pexp, (1:n_terms)']); %#ok<ASGLU>

    counts = counts(ord);
    labels = labels(ord);
    keys   = keys(ord);
    feats  = feats(ord,:);

    %% plot
    f = figure('Color','w','Units','inches','Position',[1 1 20 5]);
    ax = gca;

    x = 1:numel(counts);
    bar(x, counts);
    ylim([0 100])

    % ---- Replace tick labels with custom text objects so we can "box" them ----
    ax.XTick = x;
    ax.XTickLabel = repmat({''}, size(x)); % hide default labels
    ax.TickLabelInterpreter = 'latex';
    ax.XAxis.FontSize = 10;

    grid on; box on
    ax.FontSize = 18;
    ax.XLim = [0.5, numel(x)+0.5];
    ylabel('Count','FontSize',22);

    % % Make room for rotated tick labels
    % ax.Position(2) = ax.Position(2) + 0.08;
    % ax.Position(4) = ax.Position(4) - 0.08;

    % [left bottom width height]
    ax.Position = [ ...
        0.08 ...   % left   (increase if y-label is big)
        0.28 ...   % bottom (increase for rotated x-labels)
        0.90 ...   % width
        0.62 ...   % height (reduced but safe)
    ];


    % y-position for labels (a bit below 0)
    yl = ylim(ax);
    y_text = yl(1) - 0.06*(yl(2)-yl(1));

    % which terms to box for this equation
    if eq >= 1 && eq <= 3
        box_set = box_terms{eq};
    else
        box_set = {};
    end
    box_set = cellfun(@(s) canonical_term(s), box_set, 'UniformOutput', false);

    for ii = 1:numel(x)
        this_key = keys{ii};
        do_box = any(strcmp(this_key, box_set));

        % Choose a subtle but visible box color
        edge_col = [0.55 0.10 0.15];  % dark reddish-brown
        lw = 2;

        if do_box
            text(ax, x(ii), y_text, labels{ii}, ...
                'Interpreter','latex', ...
                'Rotation',45, ...
                'HorizontalAlignment','right', ...
                'VerticalAlignment','top', ...
                'FontSize',22, ...
                'Margin',4, ...
                'EdgeColor', edge_col, ...
                'LineWidth', lw, ...
                'BackgroundColor','none');
        else
            text(ax, x(ii), y_text, labels{ii}, ...
                'Interpreter','latex', ...
                'Rotation',45, ...
                'HorizontalAlignment','right', ...
                'VerticalAlignment','top', ...
                'FontSize',22);
        end
    end

    % keep bars visible even with below-axis text
    ax.Clipping = 'off';

    if save_figs
        pdf_path = fullfile(out_dir, sprintf('Eq%d_term_counts_ordered_boxed.pdf', eq));
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
    % - strip leading/trailing $ if any
    ts = strrep(ts, '\\cdot', '\cdot');
    ts = strrep(ts, '\cdot', '*');
    ts = strrep(ts, ' ', '');
    ts = strrep(ts, '$', '');

    % normalize common LaTeX power formatting if any
    % (keep ^ as-is, e.g., U^2)
    key = ts;
end

function feat = term_features(key)
    % Returns feature vector:
    % [group, deg, u, e, p, 0, 0]
    %
    % group:
    % 1 linear, 2 bilinear, 3 quadratic, 4 cubic, 5 other
    %
    % Parse exponents of U/E/P from canonical key like:
    % U, E, P, U*E, U^2, U^2*E, etc.

    [u,e,p] = parse_exponents(key);
    deg = u + e + p;

    % classify
    if deg == 1
        group = 1; % linear
    elseif deg == 2
        if (u==2 && e==0 && p==0) || (u==0 && e==2 && p==0) || (u==0 && e==0 && p==2)
            group = 3; % quadratic (square term)
        elseif (u<=1 && e<=1 && p<=1) && (u+e+p==2)
            group = 2; % bilinear (two distinct vars)
        else
            group = 5;
        end
    elseif deg == 3
        group = 4; % cubic (anything degree 3)
    else
        group = 5; % other (constants, higher degree, etc.)
    end

    feat = [group, deg, u, e, p, 0, 0];
end

function [u,e,p] = parse_exponents(key)
    % Extract U/E/P exponents from canonical string.
    % Supports patterns like:
    % U, E, P
    % U*E, U*P, E*P
    % U^2, E^3, U^2*E, U*E^2, etc.
    u = 0; e = 0; p = 0;

    % helper: exponent for a given symbol
    u = u + extract_exp(key, 'U');
    e = e + extract_exp(key, 'E');
    p = p + extract_exp(key, 'P');
end

function expn = extract_exp(s, symb)
    % Sum exponents for occurrences of symb in s.
    % If symb appears as "U^k", count k.
    % If appears as "U" (not followed by ^), count 1.
    expn = 0;

    % count powered occurrences first: symb^(\d+)
    tok = regexp(s, [symb '\^(\d+)'], 'tokens');
    if ~isempty(tok)
        for i = 1:numel(tok)
            expn = expn + str2double(tok{i}{1});
        end
        % remove those so we don't double count the plain symbol
        s = regexprep(s, [symb '\^\d+'], '');
    end

    % now count remaining plain occurrences
    expn = expn + numel(regexp(s, symb, 'match'));
end
