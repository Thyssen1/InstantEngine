%% Initialization
% Clear workspace, close all figures and clear command window
clear; close all; clc

% Add path to experimental plan function
addpath(strcat(erase(pwd, "Test scripts"), "\Functions"))

%
figurer;

%% User-defined decisions
bool_SaveFigures = false;

col = [1 0 0; 1 0 1; 0 0 0; 0 0.4470 0.7410; 0.8500 0.3250 0.0980;       ...
       0.9290 0.6940 0.1250; 0.4940 0.1840 0.5560; 0.4660 0.6740 0.1880; ...
       0.3010 0.7450 0.9330; 0.6350 0.0780 0.1840];

% Set experimental parameters
params.setpoints = [20 11 20 2 30]; % [T pH Co lambda0 tdose]
params.lows      = [5 10 5 1 10];   % Low values
params.highs     = [40 12 50 3 60]; % High values

% Set Instant Lab parameters
CV_input = [inf inf inf inf];

sigma_input.T      = params.setpoints(1) / CV_input(1);
sigma_input.pH     = params.setpoints(2) / CV_input(2);
sigma_input.Co     = params.setpoints(3) / CV_input(3);
sigma_input.CSdose = 333 / CV_input(4);

CV_output = [inf inf inf inf];
seed = 1:10;
Ntrain = [10 25 50 100 150 200 250 300 350 400];
Ntest = 500;

rng(0, 'twister')
for i = 1:4
    TestPlan(:,i) = icdf('Uniform', lhsdesign(Ntest,1), ...
                         params.lows(i), params.highs(i));
end
TestPlan(:,5) = 30;

TestData = instantlab(TestPlan, sigma_input, 0);
XTest = TestData.nom_input;
YTest_Yield = zeros(Ntest, 1);
YTest_Tri = zeros(Ntest, 1);
for i = 1:Ntest
    YTest_Yield(i) = TestData.out{i}(6);
    YTest_Tri(i)   = TestData.out{i}(5);
end

% Set kernel and activation functions
Kernel_Yield    = 'ardmatern32';
Kernel_Impurity = 'ardmatern32';
Sigma_Yield     = 2.4*1e-3;
Sigma_Impurity  = 1.0*1e-3;
Activation_Yield = 'relu';
Activation_Impurity = 'sigmoid';
Lambda_Yield = 1.8538*1e-6;
Lambda_Impurity = 2.4194*1e-6;
Neurons_Yield = [50 47];
Neurons_Tri   = [50 10];


for i = 1:length(Ntrain)
    for j = 1:length(seed)
        % Generate experimental plan for current 
        % number of observations and seed
        rng(seed(j), 'twister')
        for k = 1:4
            ExperimentalPlan{i,j}(:,k) = icdf('Uniform', ...
                                              lhsdesign(Ntrain(i),1), ...
                                              params.lows(k), params.highs(k));
        end
        ExperimentalPlan{i,j}(:,5) = 30;
        
        % Run InstantLab experiments with no noise
        Data{i,j} = instantlab(ExperimentalPlan{i,j}, sigma_input, 0);
        
        % Design and output matrices
        X = Data{i,j}.nom_input;
        Y_Yield = zeros(Ntrain(i), 1);
        Y_Tri   = zeros(Ntrain(i), 1);
        for p = 1:Ntrain(i)
            Y_Yield(p) = Data{i,j}.out{p}(6);
            Y_Tri(p)   = Data{i,j}.out{p}(5);
        end
       
        % Fit DOE 
        ModelDOE_Yield{i,j} = fitlm(X, Y_Yield, 'quadratic'); 
        ModelDOE_Tri{i,j}   = fitlm(X, Y_Tri, 'quadratic');
        
        % Compute estimate of kernel parameters for fitting GPR
        rng(10, 'twister')
        ModelGPR_Yield{i,j} = fitrgp(X, Y_Yield,                ...
                            'KernelFunction', Kernel_Yield,     ...
                            'Sigma', Sigma_Yield,               ...
                            'ConstantSigma', false,             ...
                            'Standardize', true);

        rng(11, 'twister')
        ModelGPR_Tri{i,j} = fitrgp(X, Y_Tri,                    ...
                            'KernelFunction', Kernel_Impurity,  ...
                            'Sigma', Sigma_Impurity,            ...
                            'ConstantSigma', false,             ...
                            'Standardize', true);

%         rng(13, 'twister')
%         ANNTri_Init = fitrnet(X, Y_Tri,                    ...
%                             'Activations', Activation_Impurity,     ...
%                             'OptimizeHyperparameters', 'Lambda',  ...
%                             'Standardize', true,                 ...
%                             'HyperparameterOptimizationOptions', ...
%                             struct("AcquisitionFunctionName",    ...
%                                    "expected-improvement-plus", ...
%                                    'ShowPlots', false,           ...
%                                    'KFold', 5,                   ...
%                                    'MaxObjectiveEvaluations', 30));

        rng(13, 'twister')
        ModelANN_Yield{i,j} = fitrnet(X, Y_Yield,       ...
                    'LayerSizes', Neurons_Yield,        ...
                    'Activations', Activation_Yield,    ...
                    'Lambda', Lambda_Yield,             ...
                    'Standardize', true);

        rng(14, 'twister')
        ModelANN_Tri{i,j} = fitrnet(X, Y_Tri,           ...
                    'LayerSizes', Neurons_Tri,          ...
                    'Activations', Activation_Impurity, ...
                    'Lambda', Lambda_Impurity,          ... 
                    'Standardize', true);

        % Set initial parameters
%         Sigma0_Yield            = gprMdlYield_Init.Sigma;
%         Sigma0_Tri              = gprMdlTri_Init.Sigma;
%         Sigma0_lb_Yield         = min(1e-2*std(Y_Yield), gprMdlYield_Init.Sigma-1e-3);
%         Sigma0_lb_Tri           = min(1e-2*std(Y_Tri), gprMdlTri_Init.Sigma-1e-3);
%         KernelParameters0_Yield = gprMdlYield_Init.KernelInformation.KernelParameters;
%         KernelParameters0_Tri   = gprMdlTri_Init.KernelInformation.KernelParameters;

        % Step 3: Fit models
%         ModelGPR_Yield{i,j} = fitrgp(X, Y_Yield,                     ...
%                         'KernelFunction', Kernel_Yield,              ...
%                         'KernelParameters', KernelParameters0_Yield, ...
%                         'Sigma', Sigma0_Yield,                       ...
%                         'SigmaLowerBound', Sigma0_lb_Yield,          ...
%                         'Standardize', true,                         ...
%                         'ConstantSigma', false);
% 
%         ModelGPR_Tri{i,j} = fitrgp(X, Y_Tri,                         ...
%                         'KernelFunction', Kernel_Impurity,           ...
%                         'KernelParameters', KernelParameters0_Yield, ...
%                         'Sigma', Sigma0_Yield,                       ...
%                         'SigmaLowerBound', Sigma0_lb_Yield,          ...
%                         'Standardize', true,                         ...
%                         'ConstantSigma', false);

        % Fit ANN
%         ModelANN_Yield{i,j} = fitrnet(X, Y_Yield, ...  
%             "LayerSizes", Neurons_Yield,    ...
%             "Standardize", true,      ...
%             "Verbose", 0,             ...
%             'Activations', Activation_Yield, ...
%             'Lambda', ANN_Yield_Init.ModelParameters.Lambda);
        %1e-6

%         ModelANN_Tri{i,j} = fitrnet(X, Y_Tri, ...  
%             "LayerSizes", Neurons_Tri,    ...
%             "Standardize", true,      ...
%             "Verbose", 0,             ...
%             'Activations', Activation_Impurity, ...
%             'Lambda', ANN_Tri_Init.ModelParameters.Lambda);
                           
%         % Compute loss
%         YieldTrainLossGPR(i,j) = loss(ModelGPR{i,j}, X(idxTrain,:), Y(idxTrain));
%         YieldTestLossGPR(i,j) = loss(ModelGPR{i,j}, X(idxVal,:), Y(idxVal));
%         YieldTrainLossANN(i,j) = loss(ModelANN{i,j}, X(idxTrain,:), Y(idxTrain));
%         YieldTestLossANN(i,j) = loss(ModelANN{i,j}, X(idxVal,:), Y(idxVal));

        % Predict on training set
        PredictionsTrainingDOE_Yield{i,j} = predict(ModelDOE_Yield{i,j}, X);
        PredictionsTrainingDOE_Tri{i,j}   = predict(ModelDOE_Tri{i,j}, X);
        PredictionsTrainingGPR_Yield{i,j} = predict(ModelGPR_Yield{i,j}, X);
        PredictionsTrainingGPR_Tri{i,j}   = predict(ModelGPR_Tri{i,j}, X);
        PredictionsTrainingANN_Yield{i,j} = predict(ModelANN_Yield{i,j}, X);    
        PredictionsTrainingANN_Tri{i,j}   = predict(ModelANN_Tri{i,j}, X); 
        
        % Predict on test set
        PredictionsTestDOE_Yield{i,j} = predict(ModelDOE_Yield{i,j}, XTest);
        PredictionsTestDOE_Tri{i,j}   = predict(ModelDOE_Tri{i,j}, XTest);
        PredictionsTestGPR_Yield{i,j} = predict(ModelGPR_Yield{i,j}, XTest);
        PredictionsTestGPR_Tri{i,j}   = predict(ModelGPR_Tri{i,j}, XTest);
        PredictionsTestANN_Yield{i,j} = predict(ModelANN_Yield{i,j}, XTest);   
        PredictionsTestANN_Tri{i,j}   = predict(ModelANN_Tri{i,j}, XTest);  
        
        % Compute statistics
        StatisticsTrainingDOE_Yield{i,j} = rs_stats(Y_Yield, PredictionsTrainingDOE_Yield{i,j});
        StatisticsTrainingDOE_Tri{i,j}   = rs_stats(Y_Tri, PredictionsTrainingDOE_Tri{i,j});
        StatisticsTestDOE_Yield{i,j}     = rs_stats(YTest_Yield, PredictionsTestDOE_Yield{i,j});
        StatisticsTestDOE_Tri{i,j}       = rs_stats(YTest_Tri, PredictionsTestDOE_Tri{i,j});
        
        StatisticsTrainingGPR_Yield{i,j} = rs_stats(Y_Yield, PredictionsTrainingGPR_Yield{i,j});
        StatisticsTrainingGPR_Tri{i,j}   = rs_stats(Y_Tri, PredictionsTrainingGPR_Tri{i,j});
        StatisticsTestGPR_Yield{i,j}     = rs_stats(YTest_Yield, PredictionsTestGPR_Yield{i,j});
        StatisticsTestGPR_Tri{i,j}       = rs_stats(YTest_Tri, PredictionsTestGPR_Tri{i,j});

        StatisticsTrainingANN_Yield{i,j} = rs_stats(Y_Yield, PredictionsTrainingANN_Yield{i,j});
        StatisticsTrainingANN_Tri{i,j}   = rs_stats(Y_Tri, PredictionsTrainingANN_Tri{i,j});
        StatisticsTestANN_Yield{i,j}     = rs_stats(YTest_Yield, PredictionsTestANN_Yield{i,j});
        StatisticsTestANN_Tri{i,j}       = rs_stats(YTest_Tri, PredictionsTestANN_Tri{i,j});
    end
end


%% R2 - Yield
DOE_Yield_R2 = zeros(length(Ntrain), length(seed));
GPR_Yield_R2 = zeros(length(Ntrain), length(seed));
ANN_Yield_R2 = zeros(length(Ntrain), length(seed));
% MEAN_DOE_Yield_R2 = zeros(length(Ntrain),1);
% MEAN_GPR_Yield_R2 = zeros(length(Ntrain),1);
% MEAN_ANN_Yield_R2 = zeros(length(Ntrain),1);
for i = 1:length(Ntrain)
    for j = 1:length(seed)
        DOE_Yield_R2(i,j) = StatisticsTestDOE_Yield{i,j}.R2;
        GPR_Yield_R2(i,j) = StatisticsTestGPR_Yield{i,j}.R2;
        ANN_Yield_R2(i,j) = StatisticsTestANN_Yield{i,j}.R2;

%         MEAN_DOE_Yield_R2(i) =  MEAN_DOE_Yield_R2(i) + ...
%                                 StatisticsTestDOE_Yield{i,j}.R2;
%         MEAN_GPR_Yield_R2(i) =  MEAN_GPR_Yield_R2(i) + ...
%                                 StatisticsTestGPR_Yield{i,j}.R2;
%         MEAN_ANN_Yield_R2(i) =  MEAN_ANN_Yield_R2(i) + ...
%                                 StatisticsTestANN_Yield{i,j}.R2;
    end
%     MEAN_DOE_Yield_R2(i) = MEAN_DOE_Yield_R2(i) / length(seed);
%     MEAN_GPR_Yield_R2(i) = MEAN_GPR_Yield_R2(i) / length(seed);
%     MEAN_ANN_Yield_R2(i) = MEAN_ANN_Yield_R2(i) / length(seed);
end

MEAN_DOE_Yield_R2 = mean(DOE_Yield_R2,2);
MEAN_GPR_Yield_R2 = mean(GPR_Yield_R2,2);
MEAN_ANN_Yield_R2 = mean(ANN_Yield_R2,2);

figure; hold all
boxplot(reshape(transpose(DOE_Yield_R2), numel(DOE_Yield_R2), 1), ...
reshape(repmat(Ntrain, length(seed), 1), length(Ntrain)*length(seed), 1))
ylim([0 1])
xlabel('$N_{obs}$, Number of observations', 'FontSize', 14)
ylabel('$R^{2}$, coefficient of determination', 'FontSize', 14)
title('MLR - Yield')
FigureTitle(1) = "Bias_Yield_DOE_R2";

figure; hold all
boxplot(reshape(transpose(GPR_Yield_R2), numel(GPR_Yield_R2), 1), ...
reshape(repmat(Ntrain, length(seed), 1), length(Ntrain)*length(seed), 1))
ylim([0 1])
xlabel('$N_{obs}$, Number of observations', 'FontSize', 14)
ylabel('$R^{2}$, coefficient of determination', 'FontSize', 14)
title('GPR - Yield')
FigureTitle(2) = "Bias_Yield_GPR_R2";

figure; hold all
boxplot(reshape(transpose(ANN_Yield_R2), numel(ANN_Yield_R2), 1), ...
reshape(repmat(Ntrain, length(seed), 1), length(Ntrain)*length(seed), 1))
ylim([0 1])
xlabel('$N_{obs}$, Number of observations', 'FontSize', 14)
ylabel('$R^{2}$, coefficient of determination', 'FontSize', 14)
title('ANN - Yield')
FigureTitle(3) = "Bias_Yield_ANN_R2";

%% Mean R2 Yield - DOE, GPR and ANN
figure; hold all
plot(Ntrain, MEAN_DOE_Yield_R2, 'x', 'MarkerSize', 8)
plot(Ntrain, MEAN_GPR_Yield_R2, 'o', 'MarkerSize', 8)
plot(Ntrain, MEAN_ANN_Yield_R2, 'd', 'MarkerSize', 8)
xlabel('$N_{obs}$, Number of observations', 'FontSize', 14)
ylabel('Average $R^{2}$', 'FontSize', 14)
title('Yield')
ylim([0 1])
legend("MLR", "GPR", "ANN", 'FontSize', 14, 'location', 'southeast', ...
    'Box', 'on')
FigureTitle(4) = "Bias_Yield_DOE_GPR_ANN_MEANR2";

%% R2 - Impurity
DOE_Impurity_R2 = zeros(length(Ntrain), length(seed));
GPR_Impurity_R2 = zeros(length(Ntrain), length(seed));
ANN_Impurity_R2 = zeros(length(Ntrain), length(seed));
% MEAN_DOE_Impurity_R2 = zeros(length(Ntrain),1);
% MEAN_GPR_Impurity_R2 = zeros(length(Ntrain),1);
% MEAN_ANN_Impurity_R2 = zeros(length(Ntrain),1);
for i = 1:length(Ntrain)
    for j = 1:length(seed)
        DOE_Impurity_R2(i,j) = StatisticsTestDOE_Tri{i,j}.R2;
        GPR_Impurity_R2(i,j) = StatisticsTestGPR_Tri{i,j}.R2;
        ANN_Impurity_R2(i,j) = StatisticsTestANN_Tri{i,j}.R2;

%         MEAN_DOE_Impurity_R2(i) =  MEAN_DOE_Impurity_R2(i) + ...
%                                     StatisticsTestDOE_Tri{i,j}.R2;
%         MEAN_GPR_Impurity_R2(i) =  MEAN_GPR_Impurity_R2(i) + ...
%                                     StatisticsTestGPR_Tri{i,j}.R2;
%         MEAN_ANN_Impurity_R2(i) =  MEAN_ANN_Impurity_R2(i) + ...
%                                     StatisticsTestANN_Tri{i,j}.R2;
    end
%     MEAN_DOE_Impurity_R2(i) = MEAN_DOE_Impurity_R2(i) / length(seed);
%     MEAN_GPR_Impurity_R2(i) = MEAN_GPR_Impurity_R2(i) / length(seed);
%     MEAN_ANN_Impurity_R2(i) = MEAN_ANN_Impurity_R2(i) / length(seed);
end

MEAN_DOE_Impurity_R2 = mean(DOE_Impurity_R2,2);
MEAN_GPR_Impurity_R2 = mean(GPR_Impurity_R2,2);
MEAN_ANN_Impurity_R2 = mean(ANN_Impurity_R2,2);

figure; hold all
boxplot(reshape(transpose(DOE_Impurity_R2), numel(DOE_Impurity_R2), 1), ...
reshape(repmat(Ntrain, length(seed), 1), length(Ntrain)*length(seed), 1))
ylim([0 1])
xlabel('$N_{obs}$, Number of observations', 'FontSize', 14)
ylabel('$R^{2}$, coefficient of determination', 'FontSize', 14)
title('MLR - Impurity')
FigureTitle(5) = "Bias_Impurity_DOE_R2";

figure; hold all
boxplot(reshape(transpose(GPR_Impurity_R2), numel(GPR_Impurity_R2), 1), ...
reshape(repmat(Ntrain, length(seed), 1), length(Ntrain)*length(seed), 1))
ylim([0 1])
xlabel('$N_{obs}$, Number of observations', 'FontSize', 14)
ylabel('$R^{2}$, coefficient of determination', 'FontSize', 14)
title('GPR - Impurity')
FigureTitle(6) = "Bias_Impurity_GPR_R2";

figure; hold all
boxplot(reshape(transpose(ANN_Impurity_R2), numel(ANN_Impurity_R2), 1), ...
reshape(repmat(Ntrain, length(seed), 1), length(Ntrain)*length(seed), 1))
ylim([0 1])
xlabel('$N_{obs}$, Number of observations', 'FontSize', 14)
ylabel('$R^{2}$, coefficient of determination', 'FontSize', 14)
title('ANN - Impurity')
FigureTitle(7) = "Bias_Impurity_ANN_R2";

%% Mean R2 Impurity - DOE, GPR and ANN
figure; hold all
plot(Ntrain, MEAN_DOE_Impurity_R2, 'x', 'MarkerSize', 8)
plot(Ntrain, MEAN_GPR_Impurity_R2, 'o', 'MarkerSize', 8)
plot(Ntrain, MEAN_ANN_Impurity_R2, 'd', 'MarkerSize', 8)
xlabel('$N_{obs}$, Number of observations', 'FontSize', 14)
ylabel('Average $R^{2}$', 'FontSize', 14)
title('Impurity')
ylim([0 1])
legend("MLR", "GPR", "ANN", 'FontSize', 14, 'location', 'southeast', ...
    'Box', 'on')
FigureTitle(8) = "Bias_Impurity_DOE_GPR_ANN_MEANR2";

%% RMSE - Yield
DOE_Yield_RMSE      = zeros(length(Ntrain), length(seed));
GPR_Yield_RMSE      = zeros(length(Ntrain), length(seed));
ANN_Yield_RMSE      = zeros(length(Ntrain), length(seed));
% MEAN_DOE_Yield_MSE = zeros(length(Ntrain), 1);
% MEAN_GPR_Yield_MSE = zeros(length(Ntrain), 1);
% MEAN_ANN_Yield_MSE = zeros(length(Ntrain), 1);
for i = 1:length(Ntrain)
    for j = 1:length(seed)
        DOE_Yield_RMSE(i,j) = StatisticsTestDOE_Yield{i,j}.RMSE;
        GPR_Yield_RMSE(i,j) = StatisticsTestGPR_Yield{i,j}.RMSE;
        ANN_Yield_RMSE(i,j) = StatisticsTestANN_Yield{i,j}.RMSE;

%         MEAN_DOE_Yield_MSE(i) =  MEAN_DOE_Yield_MSE(i) + ...
%                                 StatisticsTestDOE_Yield{i,j}.MSE;
%         MEAN_GPR_Yield_MSE(i) =  MEAN_GPR_Yield_MSE(i) + ...
%                                 StatisticsTestGPR_Yield{i,j}.MSE;
%         MEAN_ANN_Yield_MSE(i) =  MEAN_ANN_Yield_MSE(i) + ...
%                                 StatisticsTestANN_Yield{i,j}.MSE;
    end
%     MEAN_DOE_Yield_MSE(i) = MEAN_DOE_Yield_MSE(i) / length(seed);
%     MEAN_GPR_Yield_MSE(i) = MEAN_GPR_Yield_MSE(i) / length(seed);
%     MEAN_ANN_Yield_MSE(i) = MEAN_ANN_Yield_MSE(i) / length(seed);
end

MEAN_DOE_Yield_RMSE = mean(DOE_Yield_RMSE,2);
MEAN_GPR_Yield_RMSE = mean(GPR_Yield_RMSE,2);
MEAN_ANN_Yield_RMSE = mean(ANN_Yield_RMSE,2);

figure; hold all
boxplot(reshape(transpose(DOE_Yield_RMSE), numel(DOE_Yield_RMSE), 1), ...
reshape(repmat(Ntrain, length(seed), 1), length(Ntrain)*length(seed), 1))
xlabel('$N_{obs}$, Number of observations', 'FontSize', 14)
ylabel('RMSE, Root mean squared error', 'FontSize', 14)
title('MLR - Yield')
set(gca, 'YScale', 'log')
ylim([1e-3 1e1])
FigureTitle(9) = "Bias_Yield_DOE_RMSE";

figure; hold all
boxplot(reshape(transpose(GPR_Yield_RMSE), numel(GPR_Yield_RMSE), 1), ...
reshape(repmat(Ntrain, length(seed), 1), length(Ntrain)*length(seed), 1))
xlabel('$N_{obs}$, Number of observations', 'FontSize', 14)
ylabel('RMSE, Root mean squared error', 'FontSize', 14)
title('GPR - Yield')
set(gca, 'YScale', 'log')
ylim([1e-3 1e1])
FigureTitle(10) = "Bias_Yield_GPR_RMSE";

figure; hold all
boxplot(reshape(transpose(ANN_Yield_RMSE), numel(ANN_Yield_RMSE), 1), ...
reshape(repmat(Ntrain, length(seed), 1), length(Ntrain)*length(seed), 1))
xlabel('$N_{obs}$, Number of observations', 'FontSize', 14)
ylabel('RMSE, Root mean squared error', 'FontSize', 14)
title('ANN - Yield')
set(gca, 'YScale', 'log')
ylim([1e-3 1e1])
FigureTitle(11) = "Bias_Yield_ANN_RMSE";

%% Mean MSE Yield - DOE, GPR and ANN
figure; hold all
plot(Ntrain, MEAN_DOE_Yield_RMSE, 'x', 'MarkerSize', 8)
plot(Ntrain, MEAN_GPR_Yield_RMSE, 'o', 'MarkerSize', 8)
plot(Ntrain, MEAN_ANN_Yield_RMSE, 'd', 'MarkerSize', 8)
xlabel('$N_{obs}$, Number of observations', 'FontSize', 14)
ylabel('Average RMSE (Root mean squared error)', 'FontSize', 14)
title('Yield')
set(gca, 'YScale', 'log')
legend("MLR", "GPR", "ANN", 'FontSize', 14, 'location', 'northeast', ...
    'Box', 'on')
ylim([1e-3 1e1])
FigureTitle(12) = "Bias_Yield_DOE_GPR_ANN_MEANRMSE";

%% SSE - Impurity
% Pre-allocate
DOE_Impurity_RMSE = zeros(length(Ntrain), length(seed));
GPR_Impurity_RMSE = zeros(length(Ntrain), length(seed));
ANN_Impurity_RMSE = zeros(length(Ntrain), length(seed));
% MEAN_DOE_Impurity_MSE = zeros(length(Ntrain), 1);
% MEAN_GPR_Impurity_MSE = zeros(length(Ntrain), 1);
% MEAN_ANN_Impurity_MSE = zeros(length(Ntrain), 1);

for i = 1:length(Ntrain)
    for j = 1:length(seed)
        DOE_Impurity_RMSE(i,j) = StatisticsTestDOE_Tri{i,j}.RMSE;
        GPR_Impurity_RMSE(i,j) = StatisticsTestGPR_Tri{i,j}.RMSE;
        ANN_Impurity_RMSE(i,j) = StatisticsTestANN_Tri{i,j}.RMSE;
%         MEAN_DOE_Impurity_MSE(i) = MEAN_DOE_Impurity_MSE(i) + ...
%                                 StatisticsTestDOE_Tri{i,j}.SSE;
%         MEAN_GPR_Impurity_MSE(i) = MEAN_GPR_Impurity_MSE(i) + ...
%                                 StatisticsTestGPR_Tri{i,j}.SSE;
%         MEAN_ANN_Impurity_MSE(i) = MEAN_ANN_Impurity_MSE(i) + ...
%                                 StatisticsTestANN_Tri{i,j}.SSE;
    end
%     MEAN_DOE_Impurity_MSE(i) = MEAN_DOE_Impurity_MSE(i) / length(seed);
%     MEAN_GPR_Impurity_MSE(i) = MEAN_GPR_Impurity_MSE(i) / length(seed);
%     MEAN_ANN_Impurity_MSE(i) = MEAN_ANN_Impurity_MSE(i) / length(seed);
end

MEAN_DOE_Impurity_RMSE = mean(DOE_Impurity_RMSE,2);
MEAN_GPR_Impurity_RMSE = mean(GPR_Impurity_RMSE,2);
MEAN_ANN_Impurity_RMSE = mean(ANN_Impurity_RMSE,2);

figure; hold all
boxplot(reshape(transpose(DOE_Impurity_RMSE), numel(DOE_Impurity_RMSE), 1), ...
reshape(repmat(Ntrain, length(seed), 1), length(Ntrain)*length(seed), 1))
xlabel('$N_{obs}$, Number of observations', 'FontSize', 14)
ylabel('RMSE, Root mean squared error', 'FontSize', 14)
title('MLR - Impurity')
set(gca, 'YScale', 'log')
ylim([1e-3 1e1])
FigureTitle(13) = "Bias_Impurity_DOE_RMSE";

figure; hold all
boxplot(reshape(transpose(GPR_Impurity_RMSE), numel(GPR_Impurity_RMSE), 1), ...
reshape(repmat(Ntrain, length(seed), 1), length(Ntrain)*length(seed), 1))
xlabel('$N_{obs}$, Number of observations', 'FontSize', 14)
ylabel('RMSE, Root mean squared error', 'FontSize', 14)
title('GPR - Impurity')
set(gca, 'YScale', 'log')
ylim([1e-3 1e1])
FigureTitle(14) = "Bias_Impurity_GPR_RMSE";

figure; hold all
boxplot(reshape(transpose(ANN_Impurity_RMSE), numel(ANN_Impurity_RMSE), 1), ...
reshape(repmat(Ntrain, length(seed), 1), length(Ntrain)*length(seed), 1))
xlabel('$N_{obs}$, Number of observations', 'FontSize', 14)
ylabel('RMSE, Root mean squared error', 'FontSize', 14)
title('ANN - Impurity')
set(gca, 'YScale', 'log')
ylim([1e-3 1e1])
FigureTitle(15) = "Bias_Impurity_ANN_RMSE";

%% Average MSE Impurity - DOE, GPR and ANN
figure; hold all
plot(Ntrain, MEAN_DOE_Impurity_RMSE, 'x', 'MarkerSize', 8)
plot(Ntrain, MEAN_GPR_Impurity_RMSE, 'o', 'MarkerSize', 8)
plot(Ntrain, MEAN_ANN_Impurity_RMSE, 'd', 'MarkerSize', 8)
xlabel('$N_{obs}$, Number of observations', 'FontSize', 14)
ylabel('Average RMSE', 'FontSize', 14)
title('Impurity')
set(gca, 'YScale', 'log')
legend("MLR", "GPR", "ANN", 'FontSize', 14, 'location', 'northeast', ...
    'Box', 'on')
ylim([1e-3 1e1])
FigureTitle(16) = "Bias_Impurity_DOE_GPR_ANN_MEANRMSE";

%% Save figures
if bool_SaveFigures
    for i = 1:length(FigureTitle)
        saveas(figure(i), strcat(pwd, "\Images\ComputeBias\", FigureTitle(i)), 'png')
    end
end