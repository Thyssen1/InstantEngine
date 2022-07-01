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

%% Fit Gaussian Process Regression (GPR) model
% Set experimental parameters and generate experimental plan
lhs.design  = 'lhs';            % Latin hypercube sampling
lhs.lows    = [5 10 5 1];    % Low values
lhs.highs   = [40 12 50 3];  % High values
lhs.setpoints = [20 11 20 2];% Set point
lhs.samples = 1000;              % Number of samples to draw from distribution
lhs.seed    = 0;                % Seed for random number generator
lhs.plan    = exp_plan(lhs);    % Experimental plan

% Set seed for random number generator
seeds(:,1) = 1:10;
% Set parameters for running experiment
CV_input  = [inf inf inf inf; ...
             10 250 25 25];     % T pH Co CSdose

% Set standard deviations
for i = 1:size(CV_input,1)
    sigma_input(i).T      = lhs.setpoints(1) / CV_input(i,1);
    sigma_input(i).pH     = lhs.setpoints(2) / CV_input(i,2);
    sigma_input(i).Co     = lhs.setpoints(3) / CV_input(i,3);
    sigma_input(i).CSdose = 333 / CV_input(i,4);
end
CV_output = [10 10 10 10 10];   % yA ymono ydi ytri yD

lhs.plan(:,5) = 30;

% Run experiments
lhs.data       = instantlab(lhs.plan, sigma_input(1), lhs.seed);
lhs.data_noise = instantlab(lhs.plan, sigma_input(2), lhs.seed);

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
lhs.Y       = lhs.D(:,1:2);       % Measured response (no noise)
lhs.Y_noise = lhs.D_noise(:,1:2); % Measured response (noisy)
lhs.X       = lhs.D(:,3:6);       % Predictors
    
nHoldout(:,1) = [10 25 50 100 150 200 250 300 350 400 450 500 600];
p(:,1) = (lhs.samples - nHoldout) ./ lhs.samples;

sigma0 = 1000 * ones(length(nHoldout),1);
sigma0(nHoldout > 25) = 10;

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

% Set kernel for yield and impurity
Kernel_Yield          = 'ardmatern32'; % ardmatern32                    
Kernel_Yield_noise    = 'ardmatern52'; % ardmatern32
Kernel_Impurity       = 'ardmatern32';
Kernel_Impurity_noise = 'ardmatern32';

rng(10, 'twister')
% Step 1: 
% Optimize hyperparameters and compute initial estimate of kernel parameters
gprMdlYield_Init = fitrgp(XHypOpt, yDHypOpt,              ...
        'KernelFunction', Kernel_Yield,      ...
        'OptimizeHyperparameters', 'Sigma',  ...
        'Standardize', true,                 ...
        'HyperparameterOptimizationOptions', ...
        struct('ShowPlots', false,           ...
               'KFold', 5,                   ...
               'MaxObjectiveEvaluations', 50));

rng(11, 'twister')
gprMdlYield_Init_noise = fitrgp(XHypOpt, yDHypOpt_noise,              ...
        'KernelFunction', Kernel_Yield_noise,      ...
        'OptimizeHyperparameters', 'Sigma',  ...
        'Standardize', true,                 ...
        'HyperparameterOptimizationOptions', ...
        struct('ShowPlots', false,           ...
               'KFold', 5,                   ...
               'MaxObjectiveEvaluations', 50));

rng(10, 'twister')
gprMdlTri_Init = fitrgp(XHypOpt, yHHypOpt,                ...    
            'KernelFunction', Kernel_Impurity,      ...
            'OptimizeHyperparameters', 'Sigma',     ...   
            'Standardize', true,                    ...
            'HyperparameterOptimizationOptions',    ...
            struct('ShowPlots', false,              ...
                   'KFold', 5,                      ...
                   'MaxObjectiveEvaluations', 50));

rng(11, 'twister')
gprMdlTri_Init_noise = fitrgp(XHypOpt, yHHypOpt_noise,              ...    
            'KernelFunction', Kernel_Impurity_noise,    ...   
            'OptimizeHyperparameters', 'Sigma',         ...   
            'Standardize', true,                        ...
            'HyperparameterOptimizationOptions',        ...
            struct('ShowPlots', false,                  ...
                   'KFold', 5,                          ...
                   'MaxObjectiveEvaluations', 50));

% Step 2: 
% Set initial parameters
Sigma0_Yield            = gprMdlYield_Init.Sigma;
% Sigma0_lb_Yield         = min(1e-2*std(yDHypOpt), ...
%                                     1);
Sigma0_lb_Yield         = gprMdlYield_Init.ModelParameters.InitialSigmaLowerBoundTolerance;
KernelParameters0_Yield = gprMdlYield_Init.KernelInformation.KernelParameters;

Sigma0_Yield_noise            = gprMdlYield_Init_noise.Sigma;
% Sigma0_lb_Yield_noise         = min(1e-2*std(yDHypOpt_noise), ...
%                                     gprMdlYield_Init_noise.Sigma-1e-3);
Sigma0_lb_Yield_noise         = gprMdlYield_Init_noise.ModelParameters.InitialSigmaLowerBoundTolerance;
KernelParameters0_Yield_noise = gprMdlYield_Init_noise.KernelInformation.KernelParameters;

Sigma0_Tri            = gprMdlTri_Init.Sigma;
% Sigma0_lb_Tri         = min(1e-2*std(yHHypOpt), gprMdlTri_Init.Sigma-1e-3);
Sigma0_lb_Tri         = gprMdlTri_Init.ModelParameters.InitialSigmaLowerBoundTolerance;
KernelParameters0_Tri = gprMdlTri_Init.KernelInformation.KernelParameters;

Sigma0_Tri_noise            = gprMdlTri_Init_noise.Sigma;
% Sigma0_lb_Tri_noise         = min(1e-2*std(yHHypOpt_noise), ...
%                                     gprMdlTri_Init_noise.Sigma-1e-3);
Sigma0_lb_Tri_noise         = gprMdlTri_Init_noise.ModelParameters.InitialSigmaLowerBoundTolerance;
KernelParameters0_Tri_noise = gprMdlTri_Init_noise.KernelInformation.KernelParameters;

% Set initial kernel parameters found from hyperparameter optimization code
% Kernel_Yield_Params_Init = [5.9793e+04; 4.6965; 17.6527; 4.055; 0.4556];
% Sigma_Yield_Init = 0.0024;
% Sigma_Yield_LowerBound_Init = 0.0014;
% 
% Kernel_Yield_Params_Init_noise = [9.2033; 2.0635; 1.1345e+07; 1.3714; 6.6483e+03; 0.3333];
% Sigma_Yield_Init_noise = 0.0488;
% Sigma_Yield_LowerBound_Init_noise = 0.0024;
% 
% Kernel_Impurity_Params_Init = [4.1178e+04; 6.5216; 18.2146; 6.7013; 0.4025];
% Sigma_Impurity_Init = 0.0010;
% Sigma_Impurity_LowerBound_Init = 1e-06;
% 
% Kernel_Impurity_Params_Init_noise = [1.6675e+04; 4.2150; 8.3821; 3.9889; 0.2420];
% Sigma_Impurity_Init_noise = 0.0089;
% Sigma_Impurity_LowerBound_Init_noise = 0.0024;

% Pre-allocate
YieldTrainLoss        = zeros(length(p), length(seeds));
YieldValLoss          = zeros(length(p), length(seeds));
YieldTrainLoss_noise  = zeros(length(p), length(seeds));
YieldValLoss_noise    = zeros(length(p), length(seeds));

TriTrainLoss          = zeros(length(p), length(seeds));
TriValLoss            = zeros(length(p), length(seeds));
TriTrainLoss_noise    = zeros(length(p), length(seeds));
TriValLoss_noise      = zeros(length(p), length(seeds));

AverageYieldTrainLoss       = zeros(length(p), 1);
AverageYieldValLoss         = zeros(length(p), 1);
AverageYieldTrainLoss_noise = zeros(length(p), 1);
AverageYieldValLoss_noise   = zeros(length(p), 1);

AverageTriTrainLoss         = zeros(length(p), 1);
AverageTriValLoss           = zeros(length(p), 1);
AverageTriTrainLoss_noise   = zeros(length(p), 1);
AverageTriValLoss_noise     = zeros(length(p), 1);

gprMdl_Yield          = cell(length(p), length(seeds));
gprMdl_Tri            = cell(length(p), length(seeds));
gprMdl_Yield_noise    = cell(length(p), length(seeds));
gprMdl_Tri_noise      = cell(length(p), length(seeds));

n = length(lhs.Y);


for i = 1:length(p)    
    for j = 1:length(seeds)
%         disp(strcat('(i,j)=(',num2str(i),',',num2str(j),')'))

        % Cross-validation indices
        rng(j) % For reproducibility
        hpartition = cvpartition(n,'Holdout', p(i)); % Nonstratified partition
        idxTrain = training(hpartition);
        idxVal = test(hpartition);

        % Setup design and output matrices
        % Step 3: Fit models
        gprMdl_Yield{i,j} = fitrgp(lhs.X(idxTrain,:), lhs.Y(idxTrain,1), ...
                        'KernelFunction', Kernel_Yield,                  ...
                        'KernelParameters', KernelParameters0_Yield,     ...
                        'Sigma', Sigma0_Yield,                           ...
                        'SigmaLowerBound', Sigma0_lb_Yield,              ...
                        'Standardize', true,                             ...
                        'ConstantSigma', true);
        
        % Step 3: Fit models
        gprMdl_Yield_noise{i,j} = fitrgp(lhs.X(idxTrain,:), lhs.Y_noise(idxTrain,1),  ...
                        'KernelFunction', Kernel_Yield_noise,                  ...
                        'KernelParameters', KernelParameters0_Yield_noise, ...
                        'Sigma', Sigma0_Yield_noise,                           ...
                        'Standardize', true, ...
                        'ConstantSigma', true);
      

        gprMdl_Tri{i,j} = fitrgp(lhs.X(idxTrain,:), lhs.Y(idxTrain,2),  ...
                        'KernelFunction', Kernel_Impurity,              ...
                        'KernelParameters', KernelParameters0_Tri,      ...
                        'Sigma', Sigma0_Tri,                           ...
                        'Standardize', true,                            ...
                        'ConstantSigma', true);

        gprMdl_Tri_noise{i,j} = fitrgp(lhs.X(idxTrain,:), lhs.Y_noise(idxTrain,2),    ...
                        'KernelFunction', Kernel_Impurity_noise,              ...
                        'KernelParameters', KernelParameters0_Tri_noise,      ...
                        'Sigma', Sigma0_Tri_noise,                           ...
                        'Standardize', true, ...
                        'ConstantSigma', true);


        % Step 4: Compute training and test loss
        stats_Yield{i,j}        = rs_stats(lhs.Y(idxTrain,1), ...
                                    predict(gprMdl_Yield{i,j}, lhs.X(idxTrain,:)));
        stats_Yield_noise{i,j}  = rs_stats(lhs.Y_noise(idxTrain,1), ...
                                    predict(gprMdl_Yield_noise{i,j}, lhs.X(idxTrain,:)));

        stats_Tri{i,j} = rs_stats(lhs.Y(idxTrain,2), ...
                                    predict(gprMdl_Tri{i,j}, lhs.X(idxTrain,:)));
        stats_Tri_noise{i,j} = rs_stats(lhs.Y_noise(idxTrain,2), ...
                                    predict(gprMdl_Tri_noise{i,j}, lhs.X(idxTrain,:)));

        stats_YieldVal{i,j} = rs_stats(lhs.Y(idxVal,1), ...
                                    predict(gprMdl_Yield{i,j}, lhs.X(idxVal,:)));
        stats_YieldVal_noise{i,j} = rs_stats(lhs.Y_noise(idxVal,1), ...
                                    predict(gprMdl_Yield_noise{i,j}, lhs.X(idxVal,:)));

        stats_TriVal{i,j} = rs_stats(lhs.Y(idxVal,2), ...
                                    predict(gprMdl_Tri{i,j}, lhs.X(idxVal,:)));
        stats_TriVal_noise{i,j} = rs_stats(lhs.Y_noise(idxVal,2), ...
                                    predict(gprMdl_Tri_noise{i,j}, lhs.X(idxVal,:)));


        RMSE_Yield(i,j)          = stats_Yield{i,j}.RMSE;
        RMSE_Yield_noise(i,j)    = stats_Yield_noise{i,j}.RMSE;
        RMSE_Tri(i,j)            = stats_Tri{i,j}.RMSE;
        RMSE_Tri_noise(i,j)      = stats_Tri_noise{i,j}.RMSE;
        RMSE_YieldVal(i,j)       = stats_YieldVal{i,j}.RMSE;
        RMSE_YieldVal_noise(i,j) = stats_YieldVal_noise{i,j}.RMSE;
        RMSE_TriVal(i,j)         = stats_TriVal{i,j}.RMSE;
        RMSE_TriVal_noise(i,j)   = stats_TriVal_noise{i,j}.RMSE;
    end

    % Average over the seeds
    AverageYieldTrainLoss(i,1)       = mean(RMSE_Yield(i,:));
    AverageYieldTrainLoss_noise(i,1) = mean(RMSE_Yield_noise(i,:));
    AverageYieldValLoss(i,1)         = mean(RMSE_YieldVal(i,:));
    AverageYieldValLoss_noise(i,1)   = mean(RMSE_YieldVal_noise(i,:));

    AverageTriTrainLoss(i,1)       = mean(RMSE_Tri(i,:));
    AverageTriTrainLoss_noise(i,1) = mean(RMSE_Tri_noise(i,:));
    AverageTriValLoss(i,1)         = mean(RMSE_TriVal(i,:));
    AverageTriValLoss_noise(i,1)   = mean(RMSE_TriVal_noise(i,:));

%     AverageYieldTrainLoss(i,1) = mean(YieldTrainLoss(i,:));
%     AverageYieldValLoss(i,1)   = mean(YieldValLoss(i,:));
% 
%     AverageTriTrainLoss(i,1) = mean(TriTrainLoss(i,:));
%     AverageTriValLoss(i,1)   = mean(TriValLoss(i,:));
end

%%
figure; hold all
plot(nHoldout, AverageYieldTrainLoss, 'bo-')
plot(nHoldout, AverageYieldTrainLoss_noise, 'ro-')
plot(nHoldout, AverageYieldValLoss, 'bx-')
plot(nHoldout, AverageYieldValLoss_noise, 'rx-')
set(gca, 'YScale', 'log')
xlabel('$N_{obs}$, Number of observations used for training', 'FontSize', 14)
ylabel('Average loss (Root mean squared error)', 'FontSize', 14)
legend("Training", "Training (noise)", "Validation", "Validation (noise)", ...
    'location', 'southeast', 'FontSize', 10)
title('Yield (GPR)')
FigureTitle(1) = "GPR_LearningCurve_Yield_both";
ylim([10^(-5) 10^(0)])
box on

figure; hold all
plot(nHoldout, AverageTriTrainLoss, 'bo-')
plot(nHoldout, AverageTriTrainLoss_noise, 'ro-')
plot(nHoldout, AverageTriValLoss, 'bx-')
plot(nHoldout, AverageTriValLoss_noise, 'rx-')
set(gca, 'YScale', 'log')
xlabel('$N_{obs}$, Number of observations used for training', 'FontSize', 14)
ylabel('Average loss (Root nean squared error)', 'FontSize', 14)
legend("Training", "Training (noise)", "Validation", "Validation (noise)", ...
    'location', 'northeast', 'FontSize', 10)
title('Impurity (GPR)')
FigureTitle(2) = "GPR_LearningCurve_Tri_both";
ylim([10^(-5) 10^(-0)])
box on

%% Save figures
if bool_SaveFigures
    for i = 1:length(FigureTitle)
        saveas(figure(i), strcat(pwd,"\Images\LearningCurves\", ...
                                        FigureTitle(i)), 'png')
    end
end
