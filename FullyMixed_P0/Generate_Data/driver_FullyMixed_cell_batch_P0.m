% Simulate dynamics for all compartments

close all
clear
tic

%% Set Parameters
pars(1) = 0.5; % beta: transmission rate due to online network NotEngaged/Uninterested(U) -> Engaged(E)
pars(2) = 0.4; % theta: transmission rate due to offline protesting ratio NotEngaged/Uninterested(U) -> Engaged(E)
pars(3) = 0.2; % eta: fraction of Engaged(E) online from UnProtesting(UP) -> Protesting(P) offline
pars(4) = 0.1; % gamma_i: recovery rate online Engaged(E) -> DoneEngaging/DisEngaged(D)
pars(5) = 0.3; % gamma_p: recovery rate offline Protesting(P) -> DoneProtesting(R)

time_range = [0 50];
t_num = 1000;

save_data = 1; % = 1 save simulation data; = 0 not save data

%% Initialization
P0 = linspace(0.01,0.99,99);
R0 = zeros(1,length(P0));

E0 = zeros(1,length(P0));
D0 = P0;
U0 = 1 - D0;

y0 = [U0; E0; D0; P0; R0];

X5_cell = cell(size(y0,2), 1);
X4_reduced_cell = cell(size(y0,2), 1);
X3_reduced_cell = cell(size(y0,2), 1);

%% Solve
opts = odeset('AbsTol',1e-15,'RelTol',1e-12,'Stats','on','OutputFcn',@odeplot);
% opts = odeset('AbsTol',1e-15,'RelTol',1e-12);
odeFunc = @(t,y) FullyMixedModel_ODE(t,y,pars);

time_array = linspace(time_range(1),time_range(2),t_num);
for i=1:size(y0,2)
    [~,ysol] = ode45(odeFunc, time_array, y0(:,i), opts);
    X5_cell{i} = [ysol(:,1),ysol(:,2),ysol(:,3),ysol(:,4),ysol(:,5)];
    X4_reduced_cell{i} = [ysol(:,1),ysol(:,2),ysol(:,4),ysol(:,5)];
    X3_reduced_cell{i} = [ysol(:,1),ysol(:,2),ysol(:,4)];
end

%% Modify Plot
legend('U','E','D','P','R');
xlabel('Time');
ylabel('Ratio');

%% Save Data
if save_data == 1
    tempdir=pwd;
    save_path=[tempdir,'/data'];
    if ~exist(save_path, 'dir')
        mkdir(save_path);
    end

    parameters.beta = pars(1);
    parameters.theta = pars(2); 
    parameters.eta = pars(3);
    parameters.gamma_i = pars(4);
    parameters.gamma_p = pars(5);
    parameters.t_max = time_range(end);

    fname = sprintf('Fully-Mixed Network Dynamics batch %u P0=%.5g-%.5g, beta=%.3g, theta=%.3g, eta=%.3g, gamma_i=%.3g, gamma_p=%.3g, tmax=%u, t_grid_num=%u.mat',length(P0),P0(1),P0(end),pars(1),pars(2),pars(3),pars(4),pars(5),time_range(end),t_num);
    save(fullfile(save_path, fname), 'parameters', 'time_array', 'X5_cell', 'X4_reduced_cell', 'X3_reduced_cell', 'y0');
end

toc