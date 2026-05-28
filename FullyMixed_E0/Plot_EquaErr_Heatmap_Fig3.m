clear
close all

%% Load Data
cur_dir = pwd;

beta = 0.5;
theta = 0.4;
eta = 0.2;
gamma_i = 0.1;
gamma_p = 0.3;

time_range = [0 50];
t_num = 1000;
E0 = linspace(0.01,0.99,99);

start_val = 0.01;
end_val = 0.3;
num_trajs_upto = 10;

noise_ratios = linspace(0,0.1,21);
numNoise = length(noise_ratios);

rept = 100;
save_path = [cur_dir,'/WSINDy_data_3StateInput'];

fname = sprintf(['WSINDy Results for Fully-Mixed Network Dynamics E0 Batch Data ', ...
    'numTrajUpTo=%u, E0_=%.5g-%.5g, beta=%.3g, theta=%.3g, eta=%.3g, ', ...
    'gamma_i=%.3g, gamma_p=%.3g, tmax=%u, t_grid_num=%u, NumNoise=%u, ', ...
    'repeat=%u_noise_%.5g-%.5g_MSTLS.mat'], ...
    num_trajs_upto, start_val, end_val, beta, theta, eta, gamma_i, gamma_p, ...
    time_range(end), t_num, numNoise, rept, noise_ratios(1), noise_ratios(end));

load(fullfile(save_path, fname))

plot_trajNum_to = 5;

%% Equation Error Heatmap for E and P only
num_trajs_upto_data = size(EquationErr_cell, 1);
numNoise = size(EquationErr_cell, 2);
Repeat = size(EquationErr_cell, 3);

eq_letters = {'U','E','P'};

% Only plot E and P, no U
eq_to_plot = [2, 3];

% Per-equation cutoff; set any to NaN to disable cutoff
cutoffs = [NaN, NaN, NaN];

%% Ensure plot path exists
plot_path = fullfile(cur_dir, 'plots_new');
if ~exist(plot_path, 'dir')
    mkdir(plot_path);
end

for eq = eq_to_plot

    %% Build median matrix M: size trajectories x noise levels
    M = NaN(plot_trajNum_to, numNoise);

    for i = 1:plot_trajNum_to
        for j = 1:numNoise
            vals_ij = arrayfun(@(r) EquationErr_cell{i,j,r}(eq), 1:Repeat);
            M(i,j) = median(vals_ij, 'omitnan');
        end
    end

    % Plot MT so rows = noise levels, columns = number of trajectories
    MT = M.';

    %% Cutoff logic
    cutoff_val = cutoffs(eq);
    vmin = min(MT(:), [], 'omitnan');
    vmax = max(MT(:), [], 'omitnan');

    if isnan(cutoff_val)
        MT_disp = MT;
        clim = [vmin, vmax];
        apply_cutoff_labels = false;
    else
        MT_disp = min(MT, cutoff_val);
        clim = [vmin, cutoff_val];
        apply_cutoff_labels = true;
    end

    %% Plot
    figure(100 + eq); clf;
    set(gcf, 'Units','inches','Position',[6 2 5.5 10]);

    imagesc(MT_disp);
    set(gca, 'YDir','normal');

    %% Ticks and labels
    xticks(1:plot_trajNum_to);
    xticklabels(arrayfun(@(x) sprintf('%u', x), ...
        1:plot_trajNum_to, 'UniformOutput', false));

    yticks(1:numNoise);
    yticklabels(arrayfun(@(x) sprintf('%.3f', x), ...
        noise_ratios, 'UniformOutput', false));

    set(gca, 'FontSize', 20);
    axis tight;

    colormap parula;
    caxis(clim);

    cb = colorbar;
    cb.TickLabelInterpreter = 'latex';
    cb.FontSize = 20;
    cb.Label.FontSize = 22;

    % Show colorbar label only for E, not for P
    if eq == 2
        cb.Label.String = 'Median Equation Error';
    else
        cb.Label.String = '';
    end

    xlabel('Number of Trajectories', 'FontSize', 25);

    % Show y-axis label only for E, not for P
    if eq == 2
        ylabel('Noise Level', 'FontSize', 25);
    else
        ylabel('');
    end

    %% Optional red overlay for values above cutoff
    if apply_cutoff_labels
        [nrows, ncols] = size(MT);
        for j = 1:nrows
            for i = 1:min(ncols, plot_trajNum_to)
                val = MT(j,i);
                if ~isnan(val) && val > cutoff_val
                    text(i, j, sprintf('%.3f', val), ...
                        'Color','r', ...
                        'FontWeight','bold', ...
                        'HorizontalAlignment','center', ...
                        'VerticalAlignment','middle', ...
                        'FontSize',8);
                end
            end
        end
    end

    %% Force identical canvas, heatmap axes, and colorbar geometry
    fig_w = 6.4;
    fig_h = 10;
    
    set(gcf, 'Units', 'inches', 'Position', [6 2 fig_w fig_h]);
    set(gcf, 'PaperUnits', 'inches');
    set(gcf, 'PaperPosition', [0 0 fig_w fig_h]);
    set(gcf, 'PaperSize', [fig_w fig_h]);
    set(gcf, 'InvertHardcopy', 'off');
    
    ax = gca;
    ax.Units = 'normalized';
    
    % Enough left margin for y-axis label, enough right margin for colorbar label
    ax.Position = [0.23 0.12 0.50 0.82];
    
    cb.Units = 'normalized';
    cb.Position = [0.78 0.12 0.035 0.82];

    %% Save fixed-size PDF
    letter = eq_letters{eq};

    filename = sprintf('%s_EquationErr_heatmap_median_%uNoises_%uRepeat_no_cutoff.pdf', ...
        letter, numNoise, rept);

    print(gcf, fullfile(plot_path, filename), '-dpdf', '-painters');
end