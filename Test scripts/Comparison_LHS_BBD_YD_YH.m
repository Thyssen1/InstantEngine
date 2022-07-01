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
SimulationIndex = 2;        % {1,2} = {No noise, normal noise}

%% Latin Hypercube Sampling Design
% Set experimental parameters
pars.design    = 'lhs';
pars.setpoints = [20 11 20 2]; % [T pH Co lambda0 tdose]
pars.lows      = [5 10 5 1];   % Low values
pars.highs     = [40 12 50 3]; % High values
pars.reps      = 1; 
pars.n         = 4;           % Number of factors
pars.seed      = 2; % Original was 2
pars.samples   = 27;

%% Compute with noise
% Set coefficient of variation for noise input and output
switch SimulationIndex
    case 1
        CV_input = [inf inf inf inf];
    case 2
        CV_input = [10 250 25 25];
    case 3
        CV_input = [10 100 10 10];
end
CV_output = [100 100 100 100];

% 
for i = 1:size(CV_input,1)
    sigma_input(i).T      = pars.setpoints(1) / CV_input(i,1);
    sigma_input(i).pH     = pars.setpoints(2) / CV_input(i,2);
    sigma_input(i).Co     = pars.setpoints(3) / CV_input(i,3);
    sigma_input(i).CSdose = 333 / CV_input(i,4);
end

col = ["r"; "g"; "b"];

% Generate experimental plan
plan = exp_plan(pars);
plan = [plan 30*ones(length(plan), 1)];

% Generate synthetic experimental data
data = instantlab(plan, sigma_input, 0);

% Save output
YD_LHS_noisy = zeros(length(plan),1);
YH_LHS_noisy = zeros(length(plan),1);
for i = 1:length(plan)
    YD_LHS_noisy(i,1) = data.out{i}(6); % Yield  
    YH_LHS_noisy(i,1) = data.out{i}(5); % Impurity
end

% Design matrix (tdose not included for DOE, GPR and ANN)
X_LHS = plan(:,1:4);

%% Compute with no noise
SimulationIndex = 1;

%
switch SimulationIndex
    case 1
        CV_input = [inf inf inf inf];
    case 2
        CV_input = [10 250 25 25];
    case 3
        CV_input = [10 100 10 10];
end
CV_output = [100 100 100 100];

% 
for i = 1:size(CV_input,1)
    sigma_input(i).T      = pars.setpoints(1) / CV_input(i,1);
    sigma_input(i).pH     = pars.setpoints(2) / CV_input(i,2);
    sigma_input(i).Co     = pars.setpoints(3) / CV_input(i,3);
    sigma_input(i).CSdose = 333 / CV_input(i,4);
end

col = ["r"; "g"; "b"];

% Generate experimental plan
plan = exp_plan(pars);
plan = [plan 30*ones(length(plan), 1)];

% Generate synthetic experimental data
data = instantlab(plan, sigma_input, 0);

% Save output
YD_LHS_noisefree = zeros(length(plan),1);
YH_LHS_noisefree = zeros(length(plan),1);
for i = 1:length(plan)
    YD_LHS_noisefree(i,1) = data.out{i}(6); % Yield  
    YH_LHS_noisefree(i,1) = data.out{i}(5); % Impurity
end

% Design matrix (tdose not included for DOE, GPR and ANN)
X_LHS_noisefree = plan(:,1:4);

%% Box-Behnken Design
% Set experimental parameters
pars.design    = 'bbd';
pars.setpoints = [20 11 20 2]; % [T pH Co lambda0 tdose]
pars.lows      = [5 10 5 1];   % Low values
pars.highs     = [40 12 50 3]; % High values
pars.reps      = 1; 
pars.n         = 4;           % Number of factors

%% Compute with noise
SimulationIndex = 2;

% Set coefficient of variation for noise input and output
switch SimulationIndex
    case 1
        CV_input = [inf inf inf inf];
    case 2
        CV_input = [10 250 25 25];
    case 3
        CV_input = [10 100 10 10];
end
CV_output = [100 100 100 100];

% 
for i = 1:size(CV_input,1)
    sigma_input(i).T      = pars.setpoints(1) / CV_input(i,1);
    sigma_input(i).pH     = pars.setpoints(2) / CV_input(i,2);
    sigma_input(i).Co     = pars.setpoints(3) / CV_input(i,3);
    sigma_input(i).CSdose = 333 / CV_input(i,4);
end

col = ["r"; "g"; "b"];

% Generate experimental plan
plan = exp_plan(pars);
plan = [plan 30*ones(length(plan), 1)];

% Generate synthetic experimental data
data = instantlab(plan, sigma_input, 0);

% Save output
YD_BBD_noisy = zeros(length(plan),1);
YH_BBD_noisy = zeros(length(plan),1);
for i = 1:length(plan)
    YD_BBD_noisy(i,1) = data.out{i}(6); % Yield   
    YH_BBD_noisy(i,1) = data.out{i}(5); % Impurity
end

% Design matrix (tdose not included for DOE, GPR and ANN)
X_BBD = plan(:,1:4);

%% Compute BBD with no noise
SimulationIndex = 1;

% Set coefficient of variation for noise input and output
switch SimulationIndex
    case 1
        CV_input = [inf inf inf inf];
    case 2
        CV_input = [10 250 25 25];
    case 3
        CV_input = [10 100 10 10];
end
CV_output = [100 100 100 100];

% 
for i = 1:size(CV_input,1)
    sigma_input(i).T      = pars.setpoints(1) / CV_input(i,1);
    sigma_input(i).pH     = pars.setpoints(2) / CV_input(i,2);
    sigma_input(i).Co     = pars.setpoints(3) / CV_input(i,3);
    sigma_input(i).CSdose = 333 / CV_input(i,4);
end

col = ["r"; "g"; "b"];

% Generate experimental plan
plan = exp_plan(pars);
plan = [plan 30*ones(length(plan), 1)];

% Generate synthetic experimental data
data = instantlab(plan, sigma_input, 0);

% Save output
YD_BBD_noisefree = zeros(length(plan),1);
YH_BBD_noisefree = zeros(length(plan),1);
for i = 1:length(plan)
    YD_BBD_noisefree(i,1) = data.out{i}(6); % Yield   
    YH_BBD_noisefree(i,1) = data.out{i}(5); % Impurity
end

% Design matrix (tdose not included for DOE, GPR and ANN)
X_BBD_noisefree = plan(:,1:4);

%%
figure
boxplot([YD_LHS_noisy YD_LHS_noisefree], {'YD noisy', 'YD'})
title('LHS')
ylim([0 1])

disp(norm(YD_LHS_noisy-YD_LHS_noisefree))

figure
boxplot([YD_BBD_noisy YD_BBD_noisefree], {'YD noisy', 'YD'})
title('BBD')
ylim([0 1])

figure
boxplot([YD_LHS_noisefree YD_BBD_noisefree], {'YD LHS (noisefree)', 'YD BBD (noisefree)'})
ylim([0 1])

figure
boxplot([YD_LHS_noisy YD_BBD_noisy], {'YD LHS (noise)', 'YD BBD (noise)'})
ylim([0 1])

disp(norm(YD_BBD_noisy - YD_BBD_noisefree))

figure
histogram(YD_LHS_noisy, 10)
ylim([0 1])

figure
histogram(YD_BBD_noisy, 10)
ylim([0 1])

%%
figure
boxplot([YH_LHS_noisefree YH_LHS_noisy], {'YH LHS (noisefree)', 'YH LHS (noisy)'})
ylim([0 0.5])

disp(norm(YH_LHS_noisy - YH_LHS_noisefree))

figure
boxplot([YH_BBD_noisefree YH_BBD_noisy], {'YH BBD (noisefree)', 'YH BBD (noisy)'})
ylim([0 0.5])

disp(norm(YH_BBD_noisy - YH_BBD_noisefree))

figure
boxplot([YH_LHS_noisefree YH_BBD_noisefree], {'YH LHS (noisefree)', 'YH BBD (noisefree)'})
ylim([0 0.5])

figure
boxplot([YH_LHS_noisy YH_BBD_noisy], {'YH LHS (noise)', 'YH BBD (noise)'})
ylim([0 0.5])





