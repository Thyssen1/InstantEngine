%% Initialization
% Clear workspace, close all figures and clear command window
clear; close all; clc

% Add path to experimental plan function
addpath(strcat(erase(pwd, "Test scripts"), "\Functions"))

% Run custom figures script
figurer;

% Set seed for random number generator
rng(0, 'twister')

%% User-defined decisions
bool_SaveFigures = false;

%% Fit artificial neural network( ANN)
% Set experimental parameters and generate experimental plan
lhs.design  = 'lhs';            % Latin hypercube sampling
lhs.lows    = [5 10 5 1 10];    % Low values
lhs.highs   = [40 12 50 3 60];  % High values
lhs.setpoints = [20 11 20 2 30];% Set point
lhs.samples = 1000;              % Number of samples to draw from distribution
lhs.seed    = 0;                % Seed for random number generator
lhs.plan    = exp_plan(lhs);    % Experimental plan

% Set seeds for cvpartition
seeds = 1:10;

% Set noise parameters for running experiment
CV_input  = [inf inf inf inf; ...
             10 250 25 25];     % T pH Co CSdose
for i = 1:size(CV_input,1)
    sigma_input(i).T      = lhs.setpoints(1) / CV_input(i,1);
    sigma_input(i).pH     = lhs.setpoints(2) / CV_input(i,2);
    sigma_input(i).Co     = lhs.setpoints(3) / CV_input(i,3);
    sigma_input(i).CSdose = 333 / CV_input(i,4);
end
CV_output = [10 10 10 10 10];   % yA ymono ydi ytri yD

% Run experiments
lhs.data       = instantlab(lhs.plan, sigma_input(1), 0);
lhs.data_noise = instantlab(lhs.plan, sigma_input(2), 0);

% Pre-allocate
lhs.T       = zeros(size(lhs.data.out,2), 1);
lhs.pH      = zeros(size(lhs.data.out,2), 1);
lhs.Co      = zeros(size(lhs.data.out,2), 1);
lhs.lambda0 = zeros(size(lhs.data.out,2), 1);
lhs.tdose   = zeros(size(lhs.data.out,2), 1);
lhs.D       = zeros(size(lhs.data.out,2), 7);
lhs.yD      = zeros(size(lhs.data.out,2), 5);

for i = 1:size(lhs.data.out,2)
    lhs.T(i, 1)       = lhs.data.nom_input(i, 1);
    lhs.pH(i, 1)      = lhs.data.nom_input(i, 2);
    lhs.Co(i, 1)      = lhs.data.nom_input(i, 3);
    lhs.lambda0(i, 1) = lhs.data.nom_input(i, 4);
    lhs.tdose(i, 1)   = lhs.data.nom_input(i, 5);
    lhs.yD(i, 1)      = lhs.data.out{i}(6);
    lhs.ytri(i, 1)    = lhs.data.out{i}(5);

    lhs.yD_noise(i, 1) = lhs.data_noise.out{i}(6);
    lhs.ytri_noise(i, 1) = lhs.data_noise.out{i}(5);
    
    lhs.D(i,:) = [lhs.yD(i) lhs.ytri(i) lhs.T(i) lhs.pH(i) ...
                  lhs.Co(i) lhs.lambda0(i) lhs.tdose(i)];
    lhs.D_noise(i,:) = [lhs.yD_noise(i) lhs.ytri_noise(i) lhs.T(i) ...
                        lhs.pH(i) lhs.Co(i) lhs.lambda0(i) lhs.tdose(i)];
end

% Compute response and design matrix
lhs.Y       = lhs.D(:,1:2);         % Measured response (no noise)
lhs.Y_noise = lhs.D_noise(:,1:2);   % Measured response (noisy)
lhs.X       = lhs.D(:,3:7);         % Predictors

% Number of observations being held in the training set.
nHoldout(:,1) = [10 25 50 100 150 200 250 300 350 400 450 500 550 600];

% Fraction of observations held in validation set
p(:,1) = (lhs.samples - nHoldout) ./ lhs.samples;

lhs.samples = 100;
HypOptPlan = exp_plan(lhs);    % Experimental plan
HypOptPlan(:,5) = 30;
dataHypOpt       = instantlab(HypOptPlan, sigma_input(1), 11);
dataHypOpt_noise = instantlab(HypOptPlan, sigma_input(2), 11);

yDHypOpt       = zeros(lhs.samples, 1);
yDHypOpt_noise = zeros(lhs.samples, 1);
yHHypOpt       = zeros(lhs.samples, 1);
yHHypOpt_noise = zeros(lhs.samples, 1);
for i = 1:lhs.samples
    yDHypOpt(i,1)       = dataHypOpt.out{i}(6);
    yDHypOpt_noise(i,1) = dataHypOpt_noise.out{i}(6);
    yHHypOpt(i,1)       = dataHypOpt.out{i}(5);
    yHHypOpt_noise(i,1) = dataHypOpt_noise.out{i}(5);
end
XHypOpt = dataHypOpt.nom_input(:,1:4);

% Set activation functions
Yield_Activation          = 'relu';     % tanh
Yield_Activation_noise    = 'sigmoid';     % relu
Impurity_Activation       = 'tanh'; % sigmoid 
Impurity_Activation_noise = 'sigmoid'; % tanh


% Neurons_Yield = [40 40];
% Neurons_Impurity = [40 40];

Neurons_Yield = [50 47];
Neurons_Impurity = [50 10];
Neurons_Yield_Noise = [5 16];
Neurons_Impurity_Noise = [33 7];

rng(10, 'twister')
ANN_Yield_Init = fitrnet(XHypOpt, yDHypOpt,                              ...
                'LayerSizes', Neurons_Yield,                     ...
                'Activations', Yield_Activation,           ...
                'Standardize', true,                            ...
                'OptimizeHyperparameters', 'Lambda',   ...
                'HyperParameterOptimizationOptions',            ...
                struct("AcquisitionFunctionName",...
                      "expected-improvement-plus", ...
                      "MaxObjectiveEvaluations", 50,           ...
                      'ShowPlots', false, 'KFold', 5));

rng(11, 'twister')
ANN_Yield_Init_noise = fitrnet(XHypOpt, yDHypOpt_noise,                ...
                'LayerSizes', Neurons_Yield,                     ...
                'Activations', Yield_Activation_noise,                ...
                'Standardize', true,                            ...
                'OptimizeHyperparameters', 'Lambda',   ...
                'HyperParameterOptimizationOptions',            ...
                struct("AcquisitionFunctionName",               ...
                      "expected-improvement-plus",              ...
                      "MaxObjectiveEvaluations", 50,           ...
                      'ShowPlots', false, 'KFold', 5));

rng(12, 'twister')
ANN_Impurity_Init = fitrnet(XHypOpt, yHHypOpt,                ...
                'LayerSizes', Neurons_Impurity,                     ...
                'Activations', Impurity_Activation,                ...
                'Standardize', true,                            ...
                'OptimizeHyperparameters', 'Lambda',   ...
                'HyperParameterOptimizationOptions',            ...
                struct("AcquisitionFunctionName",               ...
                      "expected-improvement-plus",              ...
                      "MaxObjectiveEvaluations", 50,           ...
                      'ShowPlots', false, 'KFold', 5));

rng(13, 'twister')
ANN_Impurity_Init_noise = fitrnet(XHypOpt, yHHypOpt_noise,          ...
                'LayerSizes', Neurons_Impurity,                     ...
                'Activations', Impurity_Activation_noise,          ...
                'Standardize', true,                            ...
                'OptimizeHyperparameters', 'Lambda',                ...
                'HyperParameterOptimizationOptions',            ...
                struct("AcquisitionFunctionName",               ...
                      "expected-improvement-plus",              ...
                      "MaxObjectiveEvaluations", 50,           ...
                      'ShowPlots', false, 'KFold', 5));

Yield_Lambda            = ANN_Yield_Init.ModelParameters.Lambda;
Yield_Lambda_noise      = ANN_Yield_Init_noise.ModelParameters.Lambda;
Impurity_Lambda         = ANN_Impurity_Init.ModelParameters.Lambda;
Impurity_Lambda_noise   = ANN_Impurity_Init_noise.ModelParameters.Lambda;


% Yield_Lambda            = 3.7166e-07;
% Yield_Lambda_noise      = 0.0050;
% Impurity_Lambda         = 3.7950e-06;
% Impurity_Lambda_noise   = 2.1192e-05;

% Pre-allocate
YieldTrainLoss      = zeros(length(p), length(seeds));
YieldTestLoss       = zeros(length(p), length(seeds));
YieldTrainLossNoise = zeros(length(p), length(seeds));
YieldTestLossNoise  = zeros(length(p), length(seeds));

TriTrainLoss        = zeros(length(p), length(seeds));
TriTestLoss         = zeros(length(p), length(seeds));
TriTrainLossNoise   = zeros(length(p), length(seeds));
TriTestLossNoise    = zeros(length(p), length(seeds));

AverageYieldTrainLoss      = zeros(length(p), 1);
AverageYieldTestLoss       = zeros(length(p), 1);
AverageYieldTrainLossNoise = zeros(length(p), 1);
AverageYieldTestLossNoise  = zeros(length(p), 1);

AverageTriTrainLoss      = zeros(length(p), 1);
AverageTriTestLoss       = zeros(length(p), 1);
AverageTriTrainLossNoise = zeros(length(p), 1);
AverageTriTestLossNoise  = zeros(length(p), 1);

% Compute learning curve for ANN model
for i = 1:length(p)
    for j = 1:length(seeds)
        % Cross-validation indices
        rng(j) % For reproducibility
        n = length(lhs.Y);
        hpartition = cvpartition(n,'Holdout', p(i)); % Nonstratified partition
        idxTrain = training(hpartition);
        idxTest = test(hpartition);
        
        mdl_yield = fitrnet(lhs.X(idxTrain,:), lhs.Y(idxTrain,1), ...
            "LayerSizes", Neurons_Yield, ...
            'Activations', Yield_Activation, ...
            'Standardize', true, ...
            "Verbose", 0, ...
            'Lambda', Yield_Lambda);
        %1e-6 tanh

        mdl_tri = fitrnet(lhs.X(idxTrain,:), lhs.Y(idxTrain,2), ...
            "LayerSizes", Neurons_Yield, ...
            'Activations', Impurity_Activation, ...
            'Standardize', true, ...
            "Verbose", 0, ...
            'Lambda', Impurity_Lambda);
        %1e-6 sigmoid

        mdl_yield_noise = fitrnet(lhs.X(idxTrain,:), lhs.Y_noise(idxTrain,1), ...
            "LayerSizes", Neurons_Yield, ...
            'Activations', Yield_Activation_noise, ...
            'Standardize', true, ...
            "Verbose", 0, ...
            'Lambda', Yield_Lambda_noise);
        
        mdl_tri_noise = fitrnet(lhs.X(idxTrain,:), lhs.Y_noise(idxTrain,2), ...
            "LayerSizes", Neurons_Yield, ...
            'Activations', Impurity_Activation_noise, ...
            'Standardize', true, ...
            "Verbose", 0, ...
            'Lambda', Impurity_Lambda_noise);
    
        % Compute training and test loss
        YieldTrainLoss(i,j)      = sqrt(loss(mdl_yield, lhs.X(idxTrain,:), lhs.Y(idxTrain,1)));
        YieldTestLoss(i,j)       = sqrt(loss(mdl_yield, lhs.X(idxTest,:), lhs.Y(idxTest,1)));
        YieldTrainLossNoise(i,j) = sqrt(loss(mdl_yield_noise, lhs.X(idxTrain,:), ...
                                                         lhs.Y_noise(idxTrain,1)));
        YieldTestLossNoise(i,j)  = sqrt(loss(mdl_yield_noise, lhs.X(idxTest,:), ...
                                                         lhs.Y(idxTest,1)));

        TriTrainLoss(i,j)       = sqrt(loss(mdl_tri, lhs.X(idxTrain,:), lhs.Y(idxTrain,2)));
        TriTestLoss(i,j)        = sqrt(loss(mdl_tri, lhs.X(idxTest,:), lhs.Y(idxTest,2)));
        TriTrainLossNoise(i,j)  = sqrt(loss(mdl_tri_noise, lhs.X(idxTrain,:), ...
                                                     lhs.Y_noise(idxTrain,2)));
        TriTestLossNoise(i,j)   = sqrt(loss(mdl_tri_noise, lhs.X(idxTest,:), ...
                                                    lhs.Y(idxTest,2)));
    end

    % Average over seeds
    AverageYieldTrainLoss(i,1)      = mean(YieldTrainLoss(i,:));
    AverageYieldTestLoss(i,1)       = mean(YieldTestLoss(i,:));
    AverageYieldTrainLossNoise(i,1) = mean(YieldTrainLossNoise(i,:));
    AverageYieldTestLossNoise(i,1)  = mean(YieldTestLossNoise(i,:));

    AverageTriTrainLoss(i,1)      = mean(TriTrainLoss(i,:));
    AverageTriTestLoss(i,1)       = mean(TriTestLoss(i,:));
    AverageTriTrainLossNoise(i,1) = mean(TriTrainLossNoise(i,:));
    AverageTriTestLossNoise(i,1)  = mean(TriTestLossNoise(i,:));
end


figure; hold all
plot(nHoldout, AverageYieldTrainLoss, 'bo-')
plot(nHoldout, AverageYieldTrainLossNoise, 'ro-')
plot(nHoldout, AverageYieldTestLoss, 'bx-')
plot(nHoldout, AverageYieldTestLossNoise, 'rx-')
set(gca, 'YScale', 'log')
xlabel('$N_{train}$, Number of training observations', 'FontSize', 14)
ylabel('Average Loss (Root mean squared error)', 'FontSize', 14)
title('Yield (ANN)')
legend("Training", "Training (noise)", ...
       "Validation", "Validation (noise)", 'location', 'southeast', 'FontSize', 12)
FigureTitle(1) = "ANN_YieldLearningCurve_both";
ylim([10^(-5) 10^(0)])
box on

figure; hold all
plot(nHoldout, AverageTriTrainLoss, 'bo-')
plot(nHoldout, AverageTriTrainLossNoise, 'ro-')
plot(nHoldout, AverageTriTestLoss, 'bx-')
plot(nHoldout, AverageTriTestLossNoise, 'rx-')
set(gca, 'YScale', 'log')
xlabel('$N_{train}$, Number of training observations', 'FontSize', 14)
ylabel('Average Loss (Root mean squared error)', 'FontSize', 14)
title('Impurity (ANN)')
legend("Training", "Training (noise)", ...
       "Validation", "Validation (noise)", 'location', 'southeast', 'FontSize', 12)
FigureTitle(2) = "ANN_TriLearningCurve_both";
ylim([10^(-5) 10^(-0)])
box on

%%
if bool_SaveFigures
    for i = 1:length(FigureTitle)
        saveas(figure(i), strcat(pwd,"\Images\LearningCurves\",FigureTitle(i)), 'png')
    end
end