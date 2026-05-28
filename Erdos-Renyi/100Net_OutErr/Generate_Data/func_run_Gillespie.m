function [times,U,E,D,UP,P,R] = func_run_Gillespie(G,pars,t_max,percent_initial_S2)
% run one Gillespie stochastic simulation
%   inputs: 
%           G - sparse adjacency matrix
%           pars - list of parameters [beta, theta, eta, gamma_i, gamma_p]
%           t_max - threshold for max time the simulation runs up to
%           percent_initial_S2 - initial percent of population in state S2
%   outputs:        
%           times - list of times
%           U,E,D,UP,P,R - each is a list of ratios in the corresponding
%                   state

    %% parameters

    [N,~] = size(G);
    beta = pars(1);
    theta = pars(2);
    eta = pars(3);
    gamma_i = pars(4);
    gamma_p = pars(5);
    numInitial_S2 = round(N*percent_initial_S2);
    
    
    %% Initialization 

    max_stepnum_guess = 10*N;
    ct_step = 1;
    
    times = zeros(1,max_stepnum_guess);
    
    S1 = [N-numInitial_S2 zeros(1,max_stepnum_guess-1)];
    S2 = [numInitial_S2 zeros(1,max_stepnum_guess-1)];
    S3 = zeros(1,max_stepnum_guess);
    S4 = zeros(1,max_stepnum_guess);
    S5 = zeros(1,max_stepnum_guess);

    % initialize list of nodes in S2 and ensure none is isolated node
    degs = sum(G);
    isolated_nodes_ids = find(degs == 0);
    flag = 0;
    while flag == 0 
        S2_0_try = randsample(N,numInitial_S2);
        if ~any(ismember(S2_0_try,isolated_nodes_ids))
            flag = 1;
        end
    end
    
    ls_S2 = [S2_0_try;zeros(N-numInitial_S2,1)];
    ls_S1 = setdiff(1:N,ls_S2)'; % an always shrinking set
    
    % create a matrix recording S2 neighbors
    % entry (i,j)=1 if node i has neighbor j in S2
    S2_neighbor_mat = sparse(N,N);
    S2_neighbor_mat(:,ls_S2(1:S2(ct_step))) = G(:,ls_S2(1:S2(ct_step)));
    num_S2_neighbors = sum(S2_neighbor_mat,2);
    
    %% Calculate rates
    
    % each S1 node's transmission rate from S1 to S2
    r1_list = (theta*S3(ct_step)/N)*ones(1,S1(ct_step));
    for i = 1:S1(ct_step)
        u = ls_S1(i);
        if num_S2_neighbors(u) ~= 0
            r1_list(i) = r1_list(i) + beta*num_S2_neighbors(u);
        end
    end
    
    r1 = sum(r1_list);
    r2 = eta*S2(ct_step);
    r3 = gamma_i*S2(ct_step);
    r4 = gamma_p*S3(ct_step);
    r_total = r1+r2+r3+r4;
    
    time = exprnd(1/r_total);
    
    
    %% Update through time   
    
    while (time < t_max) && (r_total > 1e-08)
    
        rand_num = rand*r_total;
    
        if rand_num < r2 % S2 -> S3
            if S2(ct_step) < 1
                fprintf('Intended S2 -> S3 but no nodes in S2, skipping...');
            else
                ct_step = ct_step + 1;

                S2(ct_step) = S2(ct_step-1) - 1;
                S3(ct_step) = S3(ct_step-1) + 1;
                
                if S2(ct_step) == 0
                    u = ls_S2(1);
                    ls_S2 = [];
                else
                    u = randsample(ls_S2(1:S2(ct_step-1)), 1);
                    ls_S2(ls_S2 == u) = [];
                end
                
                S1(ct_step) = S1(ct_step-1);
                S4(ct_step) = S4(ct_step-1);
                S5(ct_step) = S5(ct_step-1);

                % update the matrix recording S2 neighbors
                change_r1_rates_atNodeID = find(S2_neighbor_mat(:, u) == 1);
                S2_neighbor_mat(:,u) = zeros(1,N);
        
                % update the transmission rate S1 -> S2 for S1 nodes neighboring u
                [~,S1_IDs] = ismember(change_r1_rates_atNodeID,ls_S1);
                S1_IDs = nonzeros(S1_IDs);
                r1_list(S1_IDs) = r1_list(S1_IDs)-beta;

            end
    
        elseif (rand_num >= r2) && (rand_num < r2 + r3) % S2 -> S4
            if S2(ct_step) < 1
                fprintf('Intended S2 -> S4 but no nodes in S2, skipping...');
            else
                ct_step = ct_step + 1;

                S2(ct_step) = S2(ct_step-1) - 1;
                S4(ct_step) = S4(ct_step-1) + 1;

                if S2(ct_step) == 0
                    ls_S2 = [];
                else
                    u = randsample(ls_S2(1:S2(ct_step-1)), 1);
                    ls_S2(ls_S2 == u) = [];
                end

                S1(ct_step) = S1(ct_step-1);
                S3(ct_step) = S3(ct_step-1);
                S5(ct_step) = S5(ct_step-1);

                % update the matrix recording S2 neighbors
                change_r1_rates_atNodeID = find(S2_neighbor_mat(:, u) == 1);
                S2_neighbor_mat(:,u) = zeros(1,N);
        
                % update the transmission rate S1 -> S2 for S1 nodes neighboring u
                [~,S1_IDs] = ismember(change_r1_rates_atNodeID,ls_S1);
                S1_IDs = nonzeros(S1_IDs);
                r1_list(S1_IDs) = r1_list(S1_IDs)-beta;
            end
    
        elseif (rand_num >= r2 + r3) && (rand_num < r2 + r3 + r4) % S3 -> S5
            if S3(ct_step) < 1
                fprintf('Intended S3 -> S5 but no nodes in S3, skipping...');
            else
                ct_step = ct_step + 1;

                S3(ct_step) = S3(ct_step-1) - 1;
                S5(ct_step) = S5(ct_step-1) + 1;

                S1(ct_step) = S1(ct_step-1);
                S2(ct_step) = S2(ct_step-1);
                S4(ct_step) = S4(ct_step-1);
            end
    
        else % S1 -> S2
            if S1(ct_step) < 1
                fprintf('Intended S1 -> S2 but no nodes in S1, skipping...');
            else
                ct_step = ct_step + 1;

                S1(ct_step) = S1(ct_step-1) - 1;
                S2(ct_step) = S2(ct_step-1) + 1;
                
                if S1(ct_step) == 0
                    u = ls_S1(1);
                    ls_S1 = [];
                    r1_list = [];
                else
                    p = r1*rand;
                    cumulativeSum = cumsum(r1_list+S3(ct_step-1)*theta/N);
                    FirstID = find(cumulativeSum >= p, 1);
                    u = ls_S1(FirstID);
                    ls_S1(FirstID) = [];
                    r1_list(FirstID) = [];
                end
        
                ls_S2(S2(ct_step)) = u;

                S3(ct_step) = S3(ct_step-1);
                S4(ct_step) = S4(ct_step-1);
                S5(ct_step) = S5(ct_step-1);

                % update the matrix recording S2 neighbors
                change_r1_rates_atNodeID = find(G(:, u) == 1);
                S2_neighbor_mat(:,u) = G(:,u);
                
                % update the transmission rate S1 -> S2 for S1 nodes neighboring u
                if ~isempty(r1_list)
                    [~,S1_IDs] = ismember(change_r1_rates_atNodeID,ls_S1);
                    S1_IDs = nonzeros(S1_IDs);
                    r1_list(S1_IDs) = r1_list(S1_IDs)+beta;
                end
            end
    
        end

        % updates the rates
        r1 = sum(r1_list)+S1(ct_step)*theta*S3(ct_step)/N;
        r2 = eta*S2(ct_step);
        r3 = gamma_i*S2(ct_step);
        r4 = gamma_p*S3(ct_step);
        r_total = r1+r2+r3+r4;

        times(ct_step) = time;
        time = time + exprnd(1/r_total);
    
    end

    S1 = S1(1:ct_step);
    S2 = S2(1:ct_step);
    S3 = S3(1:ct_step);
    S4 = S4(1:ct_step);
    S5 = S5(1:ct_step);
    times = times(1:ct_step);

    U = S1/N;
    E = S2/N;
    D = (S3+S4+S5)/N;
    UP = (S1+S2+S4)/N;
    P = S3/N;
    R = S5/N;

end