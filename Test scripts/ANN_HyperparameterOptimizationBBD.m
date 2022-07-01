%% Initialization
% Clear workspace, close all figures and clear command window
clear; close all; clc

% Add path to experimental plan MATLAB function
addpath(strcat(erase(pwd, "Test scripts"), "\Functions"))

%
figurer;
rng(0, 'twister')

%% Description


%% Main 
bool_SaveFigures = true;

SimulationIndex = 1; % {1 = no noise, 2 = noise}

% Set experimental design
params.design    = 'bbd';           % Box-Behnken design
params.setpoints = [20 11 20 2];    % [T pH Co lambda0 tdose]
params.lows      = [5 10 5 1];      % Low values
params.highs     = [40 12 50 3];    % High values
params.n         = 4;               % Number of factors
params.samples   = 27; 
params.seed      = 2;

% Generate experimental plan
plan = exp_plan(params);
plan = [plan 30*ones(length(plan),1)]

% Set parameters for running experiment
switch SimulationIndex
    case 1              
        CV_input = [inf inf inf inf];
    case 2
        CV_input = [10 250 25 25];
end
        
% Set standard deviations for InstantLab simulation
sigma_input.T       = params.setpoints(1) / CV_input(1);
sigma_input.pH      = params.setpoints(2) / CV_input(2);
sigma_input.Co      = params.setpoints(3) / CV_input(3);
sigma_input.CSdose  = 333 / CV_input(4);

% Generate synthetic experimental data
data = instantlab(plan, sigma_input, 0);

% Unwrap input and data from structure array
for i = 1:size(data.out,2)
    % Inputs
    T(i,1)        = data.nom_input(i, 1);   % Temperature
    pH(i,1)       = data.nom_input(i, 2);   % pH
    Co(i,1)       = data.nom_input(i, 3);   % Initial conc of A
    lambda0(i,1)  = data.nom_input(i, 4);   % Molar ratio of S and A
    tdose(i,1)    = data.nom_input(i, 5);   % Sidechain S dosing time

    % Responses
    yD(i,1)       = data.out{i}(6); % Yield
    yH(i,1)       = data.out{i}(5); % Impurity

    % Save response data in matrix
    y = [data.out{i}(1) data.out{i}(2) data.out{i}(3) ...
              data.out{i}(4) data.out{i}(5)];
end

% Set up design matrix
X = [T pH Co lambda0];

% Set hyperparameter objects
HyperParams_Yield = hyperparameters("fitrnet", X, yD);
HyperParams_Tri   = hyperparameters("fitrnet", X, yH);

% First object
HyperParams_Yield(1).Optimize = false;  % Number of layers
HyperParams_Yield(2).Optimize = false;  % Activation function
HyperParams_Yield(3).Optimize = false;  % Standardize
HyperParams_Yield(4).Optimize = true;   % Regularization strength
HyperParams_Yield(5).Optimize = false;   % Layer Weights Initializer
HyperParams_Yield(6).Optimize = false;   % Layer Biases Initializer
HyperParams_Yield(7).Optimize = true;  % Number of neurons layer 1
HyperParams_Yield(8).Optimize = true;  % Number of neurons layer 2
HyperParams_Yield(9).Optimize = false;  % Number of neurons layer 3
HyperParams_Yield(10).Optimize = false; % Number of neurons layer 4
HyperParams_Yield(11).Optimize = false; % Number of neurons layer 5

HyperParams_Yield(7).Range = [5 50];
HyperParams_Yield(8).Range = [5 50];

% Second object
HyperParams_Tri(1).Optimize = false;  % Number of layers
HyperParams_Tri(2).Optimize = false;  % Activation function
HyperParams_Tri(3).Optimize = false;  % Standardize
HyperParams_Tri(4).Optimize = true;   % Regularization strength
HyperParams_Tri(5).Optimize = false;  % Layer Weights Initializer
HyperParams_Tri(6).Optimize = false;  % Layer Biases Initializer
HyperParams_Tri(7).Optimize = true;   % Number of neurons layer 1
HyperParams_Tri(8).Optimize = true;   % Number of neurons layer 2
HyperParams_Tri(9).Optimize = false;  % Number of neurons layer 3
HyperParams_Tri(10).Optimize = false; % Number of neurons layer 4
HyperParams_Tri(11).Optimize = false; % Number of neurons layer 5

HyperParams_Tri(7).Range = [5 50];
HyperParams_Tri(8).Range = [5 50];

% Set number of neurons in layers
NeuronsYield = [40 40];
NeuronsTri   = [40 40];

% Compare ANN models with different activation functions
ActivationFunction = {'relu', 'tanh', 'sigmoid'};
for i = 1:length(ActivationFunction)
    % Fit yield
    rng(11, 'twister')
    ANN_Yield{i} = fitrnet(X, yD,                              ...
                    'LayerSizes', NeuronsYield,                     ...
                    'Activations', ActivationFunction{i},           ...
                    'Standardize', true,                            ...
                    'OptimizeHyperparameters', HyperParams_Yield,   ...
                    'HyperParameterOptimizationOptions',            ...
                    struct("AcquisitionFunctionName",...
                          "expected-improvement-plus", ...
                          "MaxObjectiveEvaluations", 100,           ...
                    'ShowPlots', false, 'KFold', 5));

    % Fit impurity
    rng(12, 'twister')
    ANN_Tri{i} = fitrnet(X, yH,                                     ...
                    'LayerSizes', NeuronsTri,                               ...
                    'Activations', ActivationFunction{i},           ...
                    'Standardize', true,                            ...
                    'OptimizeHyperparameters', HyperParams_Tri,   ...
                    'HyperParameterOptimizationOptions',            ...
                    struct("MaxObjectiveEvaluations", 100,           ...
                            "AcquisitionFunctionName",...
                            "expected-improvement-plus", ...
                    'ShowPlots', false, 'KFold', 5));

    % Compute crossvalidation object 
    rng(14, 'twister')
    cvANN_Yield{i} = fitrnet(X, yD,                                     ...
                    'LayerSizes', ANN_Yield{i}.LayerSizes,                         ...
                    'Activations', ActivationFunction{i},               ...
                    'Standardize', true,                                ...
                    'Lambda', ANN_Yield{i}.ModelParameters.Lambda, ...
                    'KFold', 5);

    rng(15, 'twister')
    cvANN_Tri{i} = fitrnet(X, yH,                                     ...
                    'LayerSizes', ANN_Tri{i}.LayerSizes,                          ...
                    'Activations', ActivationFunction{i},           ...
                    'Standardize', true,                            ...
                    'Lambda', ANN_Tri{i}.ModelParameters.Lambda, ...
                    'KFold', 5);

    % Step 3: Compute loss (mean squared error) of crossvalidation
    Loss_Yield(i) = kfoldLoss(cvANN_Yield{i});
    Loss_Tri(i)   = kfoldLoss(cvANN_Tri{i});
end

%% Display results
[~,YieldIdx] = min(Loss_Yield);
[~,ImpurityIdx] = min(Loss_Tri);

disp(strcat("Yield kernel: ", ActivationFunction{YieldIdx}))
disp(strcat("Impurity kernel: ", ActivationFunction{ImpurityIdx}))

disp('Number of neurons in each layer (Yield):')
for i = 1:length(ActivationFunction)
    disp(ANN_Yield{i}.LayerSizes)
end

disp('Number of neurons in each layer (Impurity):')
for i = 1:length(ActivationFunction)
    disp(ANN_Tri{i}.LayerSizes)
end

disp('Regularization strength (Yield)')
disp(ANN_Yield{YieldIdx}.ModelParameters.Lambda)

disp('Regularization strength (Impurity)')
disp(ANN_Tri{ImpurityIdx}.ModelParameters.Lambda)

disp('Loss_Yield Loss_Tri')
disp([Loss_Yield' Loss_Tri'])

%% Evaluate results on test sets of n=10 observations with LHS
% Create N different datasets each consisting of n=10 observations by
% creating test plans with different seeds

% % Set experimental design
% params.design    = 'lhs';           % Box-Behnken design
% params.setpoints = [20 11 20 2];    % [T pH Co lambda0 tdose]
% params.lows      = [5 10 5 1];      % Low values
% params.highs     = [40 12 50 3];    % High values
% params.n         = 4;               % Number of factors
% params.samples   = 15;
% params.seed      = 20;
% 
% % Set standard deviations for InstantLab simulation
% sigma_input.T       = 0;
% sigma_input.pH      = 0;
% sigma_input.Co      = 0;
% sigma_input.CSdose  = 0;
% 
% % Generate experimental plan
% test_plan(:,1:4) = exp_plan(params);
% test_plan(:,5) = 30*ones(params.samples,1);        
% 
% % Generate synthetic experimental data with no noise
% data = instantlab(test_plan, sigma_input, 0);

testplan.setpoints = [20 11 20 2]; % [T pH Co lambda0 tdose]
testplan.lows      = [5 10 5 1];   % Low values
testplan.highs     = [40 12 50 3]; % High values
Ntest              = 1000; % 15

% Generate test data
rng(20, 'twister')
TestPlan = zeros(Ntest, 5);
for i = 1:length(testplan.setpoints)
    TestPlan(:,i) = icdf('Uniform', lhsdesign(Ntest,1), ...
                         testplan.lows(i), testplan.highs(i));
end
% Add tdose as a constant to the test plan
TestPlan(:,5) = 30*ones(Ntest,1);

% Generate test data
test_sigma_input.T = 0;
test_sigma_input.pH = 0;
test_sigma_input.Co = 0;
test_sigma_input.CSdose = 0;

TestData = instantlab(TestPlan, test_sigma_input, 1);
XTest = TestData.nom_input(:,1:4);
YTest_Yield    = zeros(Ntest, 1);
YTest_Impurity = zeros(Ntest, 1);
for i = 1:Ntest
    YTest_Yield(i)    = TestData.out{i}(6);
    YTest_Impurity(i) = TestData.out{i}(5);
end

for i = 1:length(ActivationFunction)
    % Compute statistics
    stats_Yield(i) = rs_stats(YTest_Yield, predict(ANN_Yield{i}, TestPlan(:,1:4)));
    stats_Impurity(i) = rs_stats(YTest_Impurity, predict(ANN_Tri{i}, TestPlan(:,1:4)));

    R2_Yield(i) = stats_Yield(i).R2;
    SSE_Yield(i) = stats_Yield(i).SSE;
    MSE_Yield(i) = stats_Yield(i).MSE;
    RMSE_Yield(i) = stats_Yield(i).RMSE;
    AAD_Yield(i) = stats_Yield(i).AAD;

    R2_Impurity(i) = stats_Impurity(i).R2;
    SSE_Impurity(i) = stats_Impurity(i).SSE;
    MSE_Impurity(i) = stats_Impurity(i).MSE;
    RMSE_Impurity(i) = stats_Impurity(i).RMSE;
    AAD_Impurity(i) = stats_Impurity(i).AAD;
end

disp('R2 Yield:')
disp(R2_Yield)
disp('SSE Yield:')
disp(SSE_Yield)
disp('MSE Yield:')
disp(MSE_Yield)
disp('RMSE Yield:')
disp(RMSE_Yield)
disp('AAD Yield:')
disp(AAD_Yield)

disp('R2 Impurity:')
disp(R2_Impurity)
disp('SSE Impurity:')
disp(SSE_Impurity)
disp('MSE Impurity:')
disp(MSE_Impurity)
disp('RMSE Impurity:')
disp(RMSE_Impurity)
disp('AAD Impurity:')
disp(AAD_Impurity)

%% Figures
switch SimulationIndex
    case 1
        str = " ANN BBD (noisefree)";
    case 2
        str = " ANN BBD (noisy)";
end

figure
cats = categorical(ActivationFunction);
cats = reordercats(cats, ActivationFunction);
bar(cats, [Loss_Yield' MSE_Yield'])
legend("CV estimate", "MSE measure", 'location', 'north')
ylabel('MSE')
title(strcat('Yield ', str))
FigureTitle(1) = "ANN_BBD_HyperparameterOptCurves_Yield";

figure
cats = categorical(ActivationFunction);
cats = reordercats(cats, ActivationFunction);
bar(cats, [Loss_Tri' MSE_Yield'])
legend("CV estimate", "MSE measure", 'location', 'northwest')
ylabel('MSE')
title(strcat('Impurity ', str))
FigureTitle(2) = "ANN_BBD_HyperparameterOptCurves_Impurity";

%% No noise
% Yield kernel: tanh
% Impurity kernel: sigmoid
% Number of neurons in each layer (Yield):
%     49    43
% 
%     46    14
% 
%     36    47
% 
% Number of neurons in each layer (Impurity):
%      7     9
% 
%     37    49
% 
%     44     5
% 
% Regularization strength (Yield)
%     0.0029
% 
% Regularization strength (Impurity)
%    5.5847e-07
% 
% Loss_Yield Loss_Tri
%     0.0114    0.0053
%     0.0110    0.0087
%     0.0135    0.0047
% 
% R2 Yield:
%     0.9018    0.9388    0.9424
% 
% SSE Yield:
%     7.1895    4.2996    3.5950
% 
% RMSE Yield:
%     0.0848    0.0656    0.0600
% 
% AAD Yield:
%     0.0627    0.0491    0.0471
% 
% R2 Impurity:
%     0.9295    0.9121    0.9677
% 
% SSE Impurity:
%     0.3488    0.5903    0.3225
% 
% RMSE Impurity:
%     0.0187    0.0243    0.0180
% 
% AAD Impurity:
%     0.0093    0.0115    0.0082

%% Noise
% Yield kernel: tanh
% Impurity kernel: sigmoid
% Number of neurons in each layer (Yield):
%     50    40
% 
%     18    50
% 
%      5    43
% 
% Number of neurons in each layer (Impurity):
%      7    39
% 
%     31     5
% 
%     46     5
% 
% Regularization strength (Yield)
%     0.0036
% 
% Regularization strength (Impurity)
%    3.9852e-07
% 
% Loss_Yield Loss_Tri
%     0.0162    0.0035
%     0.0111    0.0052
%     0.0123    0.0021
% 
% R2 Yield:
%     0.9396    0.9407    0.9345
% 
% SSE Yield:
%     3.9324    3.8079    4.1283
% 
% RMSE Yield:
%     0.0627    0.0617    0.0643
% 
% AAD Yield:
%     0.0460    0.0452    0.0467
% 
% R2 Impurity:
%     0.9130    0.8991    0.9577
% 
% SSE Impurity:
%     0.4366    0.5336    0.2660
% 
% RMSE Impurity:
%     0.0209    0.0231    0.0163
% 
% AAD Impurity:
%     0.0114    0.0124    0.0089

%% Save figures
if bool_SaveFigures
    for i = 1:length(FigureTitle) 
        saveas(figure(i), ...
            strcat(pwd, "\Images\ANN_HyperparameterOptimization\", ...
                        FigureTitle(i), "_NoiseLevel_", num2str(SimulationIndex)), 'png')
    end
end



%% No noise
% Yield kernel: relu
% Impurity kernel: relu
% Regularization strength (Yield)
%    4.9937e-05
% 
% Regularization strength (Impurity)
%    4.2694e-06


%% Noise
% Yield kernel: tanh
% Impurity kernel: relu
% Regularization strength (Yield)
%     0.0029
% 
% Regularization strength (Impurity)
%    2.0221e-06
% 
% Loss_Yield Loss_Tri
%     0.0144    0.0030
%     0.0110    0.0047
%     0.0144    0.0031




