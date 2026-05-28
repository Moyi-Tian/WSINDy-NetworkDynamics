function S = progressive_ic_selection(ICs, Kmax)
% ICs: sorted vector of initial conditions (e.g., 0.01:0.01:0.30)
% Kmax: largest set size you want (e.g., 10)
% S: cell array; S{k} gives the indices of the k selected ICs (nested)

    n = numel(ICs);
    S = cell(Kmax,1);

    % 1) seed with left-middle index
    mid = floor((n+1)/2);
    selected = mid;

    for k = 1:Kmax
        if k > 1
            % compute distance to nearest selected point for each candidate
            d = inf(n,1);
            for i = 1:n
                d(i) = min(abs(ICs(i) - ICs(selected)));
            end
            d(selected) = -inf;                        % don't reselect chosen ones

            % pick farthest; tie-break by smallest index for determinism
            maxd = max(d);
            cand = find(d == maxd, 1, 'first');       % leftmost among ties
            selected = sort([selected; cand]);        % keep indices sorted
        end
        S{k} = selected;
    end
end
