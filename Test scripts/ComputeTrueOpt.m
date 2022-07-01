%% Initialization
% Clear workspace, close all figures and clear command window
clear; close all; clc

% Add path to experimental plan function
addpath(strcat(erase(pwd, "Test scripts"), "\Functions"))
addpath(strcat(erase(pwd, "Test scripts"), "\easyGSA"))

% Run custom figures script
figurer;
rng(0, 'twister');

%% User-defined decisions
bool_SaveFigures = false;

%% Latin Hypercube Sampling Design
% Set experimental parameters
pars.design    = 'lhs';
pars.setpoints = [20 11 20 2]; % [T pH Co lambda0 tdose]
pars.lows      = [5 10 5 1];   % Low values
pars.highs     = [40 12 50 3]; % High values
pars.reps      = 1; 
pars.n         = 4;           % Number of factors
pars.seed      = 0;
pars.samples   = 250000;

%% Set Instant Lab parameters
% Set coefficient of variation for noise input and output
CV_input = [inf inf inf inf];

%
sigma_input.T      = 0;
sigma_input.pH     = 0;
sigma_input.Co     = 0;
sigma_input.CSdose = 0;

col = ["r"; "g"; "b"];

% Generate experimental plan
plan = exp_plan(pars);
plan = [plan 30*ones(length(plan), 1)];

% Generate synthetic experimental data
data = instantlab(plan, sigma_input, 0);

% Save output
YD = zeros(length(plan),1);
for i = 1:length(plan)
    YD(i,1) = data.out{i}(6); % Yie     
    YH(i,1) = data.out{i}(5); % Impurity

    y(i,:) = [data.out{i}(1) data.out{i}(2) data.out{i}(3) ...
                  data.out{i}(4) data.out{i}(5) data.out{i}(6)];
end

% Design matrix (tdose not included for DOE, GPR and ANN)
X = plan(:,1:4);

%% Compute initial feasible points for the six cases
reacstruc = reacstruccreate();
price = [0; reacstruc.optim.Price.SCrel; 2];

% Subcases:
% price  | ub[tri]
% 0      | 1 
% 0      | 0.01
% 0.5370 | 1
% 0.5370 | 0.01
% 2      | 1
% 2      | 0.01

k = find(YH <= 0.005);
k = find(YH <= 0.01);

for i = 1:length(price)
    % Compute optimum with no upper bound on triacylated
    [cost(2*i-1,1), IdxOpt(2*i-1,1)] = min((1 + price(i)*X(:,4)) ./ YD);

    % Compute optimum with upper bound on triacylated
    [cost(2*i,1), IdxTemp]  = min((1 + price(i)*X(k,4)) ./ YD(k));

    % Get true index 
    IdxOpt(2*i,1) = k(IdxTemp);
end


% idx_MM = [3 2 4 1];
idx_MM = [4 2 1 3];
% disp(plan(k(I),:))
% z0 = plan(k(I),:);
% z0 = z0([4 2 1 3]);
z0 = X(IdxOpt, idx_MM);

%% Compute true optimum
reacstruc = reacstruccreate();

% Set bounds 
for i = 1:length(reacstruc.optim.var)
    % Define helper variables
    LB(i) = reacstruc.optim.LB(i);
    UB(i) = reacstruc.optim.UB(i);

    for j = 1:6
        % Compute normalized initial condition
        x0(j,i) = (z0(j,i)-LB(i)) / (UB(i)-LB(i));
    end

    % Set normalized bounds
    con(1,i) = 0;
    con(2,i) = 1;
end

% Compute optima for each price and with/without nonlinear constraints
for i = 1:length(price)
    % Set price of current iteration
    reacstruc.optim.Price.SCrel = price(i);

    % Set bounds 
    pars.UB = reacstruc.optim.UB;
    pars.LB = reacstruc.optim.LB;
    fun = @(x) reacoptim(x,reacstruc);
    options = optimset('TolX', 1e-6, 'TolFun', 1e-6);
    
    % Compute solutions with and without nonlinear constraints
    [sol(2*i-1,:), fval(2*i-1,1)] = fminsearchcon(fun, x0(2*i-1,:), ...
                                con(1,:), con(2,:), [],[],[], options);
    [sol(2*i,:),   fval(2*i,1)] = fminsearchcon(fun, x0(2*i,:),   ...
                                con(1,:), con(2,:), [], [], @(x) reaccon(x,pars,reacstruc), options);

    options_fmincon = optimset('TolX', 1e-9, 'TolFun', 1e-9);
    % Compute solutions with fmincon
    [sol_fmincon(i,:), fval_fmincon(i,1)] = fmincon(fun, x0(2*i,:), [], [], [], [], ...
                    con(1,:), con(2,:), @(x) reaccon(x,pars,reacstruc), options);


    % Compute minima
    opt(2*i-1,:) = sol(2*i-1,:) .* (UB-LB) + LB;
    opt(2*i,:)   = sol(2*i,:) .* (UB-LB) + LB;
    opt_fmincon(i,:) = sol_fmincon(i,:) .* (UB-LB) + LB;

    % Compute yield and impurity concentrations for global minima
    [yDopt(2*i-1,1), yHopt(2*i-1,1)] = reacsim_wrapper([opt(2*i-1,:) 30], reacstruc);
    [yDopt(2*i,1), yHopt(2*i,1)]     = reacsim_wrapper([opt(2*i,:) 30], reacstruc);
    [yDopt_fmincon(i,1), yHopt_fmincon(i,1)] = reacsim_wrapper([opt_fmincon(i,:) 30], reacstruc);
end



% sol
[opt(2:2:6,:) yDopt(2:2:6,:) yHopt(2:2:6,:) fval(2:2:6,:)]

[opt_fmincon yDopt_fmincon yHopt_fmincon fval_fmincon]
% 
% % Notes:
% [opt(2:2:6,:) yDopt(2:2:6,:) yHopt(2:2:6,:) fval(2:2:6,:)]
% 
% [opt_fmincon yDopt_fmincon yHopt_fmincon fval_fmincon]
% 
% ans =
% 
%     2.9998   11.9054   39.9995   29.9989    0.9850    0.0100    1.0152
%     2.0773   11.0715   40.0000   50.0000    0.9671    0.0100    2.1870
%     2.0504   10.9662   40.0000   49.9999    0.9588    0.0100    5.3184
% 
% 
% ans =
% 
%     2.7592   11.8461   39.2548   32.8659    0.9828    0.0097    1.0174
%     2.0452   11.1251   33.7999   41.2993    0.9248    0.0091    2.2673
%     2.0675   11.0390   39.9402   49.8552    0.9640    0.0099    5.3253

% Old starting at x0 with 0.005 as lower bound
% ans =
% 
%     2.8324   11.9920   39.9983   44.0064    0.9853    0.0100    1.0149
%     2.0782   11.0742   39.9817   49.6572    0.9671    0.0100    2.1876
%     2.0519   10.9681   39.9457   49.1103    0.9587    0.0100    5.3218
% 
% 
% ans =
% 
%     2.7366   11.7283   39.4057   25.9515    0.9822    0.0097    1.0181
%     2.0646   11.0367   39.8917   49.5353    0.9624    0.0098    2.1906
%     1.9347   10.5775   38.3022   35.0131    0.8721    0.0099    5.5753

%% Notes
% fmincon results:

% sol_fmincon =
% 
%     0.8683    0.8641    0.9830    0.4656
%     0.5323    0.5183    0.9969    0.9897
%     0.4674    0.2887    0.9515    0.6670
% 
% opt_fmincon =
% 
%     2.7366   11.7283   39.4057   25.9515
%     2.0646   11.0367   39.8917   49.5353
%     1.9347   10.5775   38.3022   35.0131
% 
% fval_fmincon =
% 
%     1.0181
%     2.1906
%     5.5753
% 
% yDopt_fmincon =
% 
%     0.9866
%     0.9853
%     0.9746
% 
% yHopt_fmincon =
% 
%     0.0117
%     0.0100
%     0.0158

%%
% x0 =
% 
%     0.8909    0.9344    0.9962    0.8661
%     0.7856    0.9518    0.9939    0.4308
%     0.5484    0.4114    0.9998    0.9665
%     0.5173    0.7396    0.9781    0.8641
%     0.5381    0.4158    0.9902    0.9635
%     0.4707    0.6089    0.9722    0.8844

% z0 =
% 
%     2.7819   11.8689   39.8675   43.9748
%     2.5713   11.9036   39.7866   24.3864
%     2.0967   10.8227   39.9939   48.4912
%     2.0346   11.4792   39.2336   43.8840
%     2.0761   10.8316   39.6578   48.3564
%     1.9413   11.2177   39.0285   44.8002

% cost =
% 
%     1.0139
%     1.0951
%     2.1742
%     2.3567
%     5.2968
%     5.6887

% YD(IdxOpt)
% 
% ans =
% 
%     0.9863
%     0.9132
%     0.9778
%     0.8879
%     0.9727
%     0.8583

% YH(IdxOpt)
% 
% ans =
% 
%     0.0124
%     0.0050
%     0.0165
%     0.0050
%     0.0141
%     0.0050


% Optimum found with fminsearchcon
% sol
%     0.9416    1.0000    1.0000    0.9918
%     0.8606    0.9652    0.9998    0.8694
%     0.5386    0.3845    1.0000    1.0000
%     0.5408    0.5284    0.9905    0.8858
%     0.5300    0.3474    1.0000    1.0000
%     0.5250    0.4823    1.0000    1.0000

% opt
%     2.8833   12.0000   40.0000   49.6326
%     2.7213   11.9304   39.9945   44.1231
%     2.0773   10.7691   40.0000   50.0000
%     2.0816   11.0569   39.6663   44.8617
%     2.0600   10.6948   40.0000   50.0000
%     2.0501   10.9647   40.0000   49.9996

% fval
%     1.0136
%     1.0151
%     2.1705
%     2.1946
%     5.2820
%     5.3184

% yDopt
%     0.9866
%     0.9851
%     0.9746
%     0.9648
%     0.9691
%     0.9586

% yHopt
%     0.0117
%     0.0100
%     0.0158
%     0.0100
%     0.0155
%     0.0100






%% Multistart with fmincon
% rng default % For reproducibility
% opts = optimoptions(@fmincon);
% problem = createOptimProblem('fmincon','objective',...
%     fun,'x0',x0,'lb',con(1,:),'ub',con(2,:), 'nonlcon', ...
%     @(x) reaccon(x,pars,reacstruc), 'options', opts);
% ms = MultiStart;
% [x,f] = run(ms,problem,200)
% 
% disp(x .* (UB-LB) + LB)

% 179 out of 200 local solver runs converged with a positive local solver exit flag.
% 
% x =
% 
%     0.5459    0.5548    0.9914    0.8865    0.3120
% 
% 
% f =
% 
%     2.1971
% 
%     2.0918   11.1096   39.6977   44.8938   25.5976

