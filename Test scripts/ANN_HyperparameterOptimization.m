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
params.design    = 'lhs';           % Box-Behnken design
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
    rng(10, 'twister')
    cvANN_Yield{i} = fitrnet(X, yD,                                     ...
                    'LayerSizes', ANN_Yield{i}.LayerSizes,                         ...
                    'Activations', ActivationFunction{i},               ...
                    'Standardize', true,                                ...
                    'Lambda', ANN_Yield{i}.ModelParameters.Lambda, ...
                    'KFold', 5);

    rng(11, 'twister')
    cvANN_Tri{i} = fitrnet(X, yH,                                     ...
                    'LayerSizes', ANN_Tri{i}.LayerSizes,                          ...
                    'Activations', ActivationFunction{i},           ...
                    'Standardize', true,                            ...
                    'Lambda', ANN_Tri{i}.ModelParameters.Lambda, ...
                    'KFold', 5);

    % Step 3: Compute loss (mean squared error) of crossvalidation
    Loss_Yield(i) = kfoldLoss(cvANN_Yield{i});
    Loss_Tri(i)   = kfoldLoss(cvANN_Tri{i});

%     % Fit models
%     rng(10, 'twister')
%     ANN_Yield{i} = fitrnet(X, yD,                                     ...
%                     'LayerSizes', NeuronsYield,                         ...
%                     'Activations', ActivationFunction{i},               ...
%                     'Standardize', true,                                ...
%                     'Lambda', ANN_Yield_Init{i}.ModelParameters.Lambda);
% 
%     rng(11, 'twister')
%     ANN_Tri{i} = fitrnet(X, yD,                                     ...
%                     'LayerSizes', NeuronsYield,                         ...
%                     'Activations', ActivationFunction{i},               ...
%                     'Standardize', true,                                ...
%                     'Lambda', ANN_Yield_Init{i}.ModelParameters.Lambda);
end

%% Display results
[~,YieldIdx] = min(Loss_Yield);
[~,ImpurityIdx] = min(Loss_Tri);

disp(strcat("Best yield kernel by CV: ", ActivationFunction{YieldIdx}))
disp(strcat("Best impurity kernel by CV: ", ActivationFunction{ImpurityIdx}))

disp('Number of neurons in each layer (Yield)')
for i = 1:length(ActivationFunction)
    disp(ANN_Yield{i}.LayerSizes)
end

disp('Number of neurons in each layer (Impurity)')
for i = 1:length(ActivationFunction)
    disp(ANN_Tri{i}.LayerSizes)
end

disp('Regularization strengths (Yield)')
for i = 1:length(ActivationFunction)
    disp(ANN_Yield{i}.ModelParameters.Lambda)
end

disp('Regularization strengths (Impurity)')
for i = 1:length(ActivationFunction)
    disp(ANN_Tri{i}.ModelParameters.Lambda)
end

%% Evaluate results on test sets of n=10 observations with LHS
% Create N different datasets each consisting of n=10 observations by
% creating test plans with different seeds

% Set experimental design
% params.design    = 'lhs';           % Box-Behnken design
% params.setpoints = [20 11 20 2];    % [T pH Co lambda0 tdose]
% params.lows      = [5 10 5 1];      % Low values
% params.highs     = [40 12 50 3];    % High values
% params.n         = 4;               % Number of factors
% params.samples   = 15;
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
% 
% yD = zeros(params.samples, length(ActivationFunction));
% yH = zeros(params.samples, length(ActivationFunction));

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
% %     params.seed = seeds(i);
% 
%     % Unwrap input and data from structure array
%     for j = 1:size(data.out,2)    
%         % Responses
%         yD(j,i) = data.out{j}(6); % Yield
%         yH(j,i) = data.out{j}(5); % Impurity
%     end

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

%% Figure
switch SimulationIndex
    case 1
        str = " ANN LHS (noisefree)";
    case 2
        str = " ANN LHS (noisy)";
end

figure
cats = categorical(ActivationFunction);
cats = reordercats(cats, ActivationFunction);
bar(cats, [Loss_Yield' MSE_Yield'])
legend("CV estimate", "MSE measure", 'location', 'north')
ylabel('MSE')
title(strcat('Yield ', str))
FigureTitle(1) = "ANN_LHS_HyperparameterOptCurves_Yield";

figure
cats = categorical(ActivationFunction);
cats = reordercats(cats, ActivationFunction);
bar(cats, [Loss_Tri' MSE_Impurity'])
legend("CV estimate", "MSE measure", 'location', 'north')
ylabel('MSE')
title(strcat('Impurity', str))
FigureTitle(2) = "ANN_LHS_HyperparameterOptCurves_Impurity";


%% No noise
% Best yield kernel by CV: relu
% Best impurity kernel by CV: sigmoid
% Number of neurons in each layer (Yield)
%     50    47
% 
%     49     6
% 
%     49     5
% 
% Number of neurons in each layer (Impurity)
%     49    31
% 
%     49    43
% 
%     50    10
% 
% Regularization strengths (Yield)
%    1.8538e-06
% 
%    8.3094e-04
% 
%    1.8421e-05
% 
% Regularization strengths (Impurity)
%    2.0985e-06
% 
%    2.5667e-06
% 
%    2.4194e-06
% 
% R2 Yield:
%     0.9306    0.9516    0.9665
% 
% SSE Yield:
%     4.2100    2.9261    2.0774
% 
% RMSE Yield:
%     0.0649    0.0541    0.0456
% 
% AAD Yield:
%     0.0425    0.0373    0.0293
% 
% R2 Impurity:
%     0.8764    0.8479    0.9634
% 
% SSE Impurity:
%     0.6474    0.7289    0.1694
% 
% RMSE Impurity:
%     0.0254    0.0270    0.0130
% 
% AAD Impurity:
%     0.0113    0.0185    0.0075 

%% Noise
% Best yield kernel by CV: tanh
% Best impurity kernel by CV: sigmoid
% Number of neurons in each layer (Yield)
%     50    12
% 
%      5    16
% 
%      5    13
% 
% Number of neurons in each layer (Impurity)
%      5    46
% 
%     26    21
% 
%     33     7
% 
% Regularization strengths (Yield)
%     0.0013
% 
%     0.0110
% 
%     0.0017
% 
% Regularization strengths (Impurity)
%    1.9071e-04
% 
%    2.2430e-04
% 
%    5.7982e-05
% 
% R2 Yield:
%     0.8866    0.8210    0.8167
% 
% SSE Yield:
%     7.1649   10.8840   11.0725
% 
% RMSE Yield:
%     0.0846    0.1043    0.1052
% 
% AAD Yield:
%     0.0656    0.0805    0.0819
% 
% R2 Impurity:
%     0.8824    0.5925    0.8190
% 
% SSE Impurity:
%     0.5829    1.9996    0.8248
% 
% RMSE Impurity:
%     0.0241    0.0447    0.0287
% 
% AAD Impurity:
%     0.0128    0.0296    0.0170


%% Save figures
if bool_SaveFigures
    for i = 1:length(FigureTitle) 
        saveas(figure(i), ...
            strcat(pwd, "\Images\ANN_HyperparameterOptimization\", ...
                        FigureTitle(i), "_NoiseLevel_", num2str(SimulationIndex)), 'png')
    end
end




%% No noise
% Yield kernel: tanh
% Impurity kernel: sigmoid
% Regularization strength (Yield)
%    3.7166e-07
% 
% Regularization strength (Impurity)
%    3.7950e-06

%% Noise
% Yield kernel: relu
% Impurity kernel: tanh
% Regularization strength (Yield)
%     0.0048
% 
% Regularization strength (Impurity)
%    6.8101e-06
%
%
