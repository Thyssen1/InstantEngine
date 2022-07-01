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
bool_SaveFigures = true;

%% Fit linear model to data
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
    
% nHoldout(:,1) = [25 50 75 100 125 150];
nHoldout(:,1) = [25 50 75 100 150 200 250 300 350 400 450 500 550 600];
p(:,1) = (lhs.samples - nHoldout) ./ lhs.samples;

% Compute hyperparameter Sigma and estimate of initial values 
% of kernel parameters on a subset of the data
% n = length(lhs.Y);
% idxEst = randsample(lhs.samples, 125);

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

MLR_Yield          = cell(length(p), length(seeds));
MLR_Tri            = cell(length(p), length(seeds));
MLR_Yield_noise    = cell(length(p), length(seeds));
MLR_Tri_noise      = cell(length(p), length(seeds));

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
        MLR_Yield{i,j} = fitlm(lhs.X(idxTrain,:), lhs.Y(idxTrain,1),  ...
                        'quadratic');
        
        MLR_Yield_noise{i,j} = fitlm(lhs.X(idxTrain,:), lhs.Y_noise(idxTrain,1),  ...
                        'quadratic');
      
        MLR_Tri{i,j} = fitlm(lhs.X(idxTrain,:), lhs.Y(idxTrain,2),    ...
                        'quadratic');

        MLR_Tri_noise{i,j} = fitlm(lhs.X(idxTrain,:), lhs.Y_noise(idxTrain,2),    ...
                        'quadratic');

        % Step 4: Compute training and test loss
        stats_Yield{i,j}        = rs_stats(lhs.Y(idxTrain,1), ...
                                    predict(MLR_Yield{i,j}, lhs.X(idxTrain,:)));
        stats_Yield_noise{i,j}  = rs_stats(lhs.Y_noise(idxTrain,1), ...
                                    predict(MLR_Yield_noise{i,j}, lhs.X(idxTrain,:)));

        stats_Tri{i,j} = rs_stats(lhs.Y(idxTrain,2), ...
                                    predict(MLR_Tri{i,j}, lhs.X(idxTrain,:)));
        stats_Tri_noise{i,j} = rs_stats(lhs.Y_noise(idxTrain,2), ...
                                    predict(MLR_Tri_noise{i,j}, lhs.X(idxTrain,:)));

        stats_YieldVal{i,j} = rs_stats(lhs.Y(idxVal,1), ...
                                    predict(MLR_Yield{i,j}, lhs.X(idxVal,:)));
        stats_YieldVal_noise{i,j} = rs_stats(lhs.Y_noise(idxVal,1), ...
                                    predict(MLR_Yield_noise{i,j}, lhs.X(idxVal,:)));

        stats_TriVal{i,j} = rs_stats(lhs.Y(idxVal,2), ...
                                    predict(MLR_Tri{i,j}, lhs.X(idxVal,:)));
        stats_TriVal_noise{i,j} = rs_stats(lhs.Y_noise(idxVal,2), ...
                                    predict(MLR_Tri_noise{i,j}, lhs.X(idxVal,:)));


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


figure; hold all
plot(nHoldout, AverageYieldTrainLoss, 'bo-')
plot(nHoldout, AverageYieldTrainLoss_noise, 'ro-')
plot(nHoldout, AverageYieldValLoss, 'bx-')
plot(nHoldout, AverageYieldValLoss_noise, 'rx-')
set(gca, 'YScale', 'log')
xlabel('$N_{obs}$, Number of observations used for training', 'FontSize', 14)
ylabel('Loss (Root mean squared error)', 'FontSize', 14)
legend("Training", "Training (noise)", "Test", "Test (noise)", ...
    'location', 'southeast', 'FontSize', 14)
title('Yield (MLR)')
FigureTitle(1) = "MLR_LearningCurve_Yield_both";
ylim([10^(-2) 10^(0)])
box on

figure; hold all
plot(nHoldout, AverageTriTrainLoss, 'bo-')
plot(nHoldout, AverageTriTrainLoss_noise, 'ro-')
plot(nHoldout, AverageTriValLoss, 'bx-')
plot(nHoldout, AverageTriValLoss_noise, 'rx-')
set(gca, 'YScale', 'log')
xlabel('$N_{obs}$, Number of observations used for training', 'FontSize', 14)
ylabel('Loss (Root mean squared error)', 'FontSize', 14)
legend("Training", "Training (noise)", "Test", "Test (noise)", ...
    'location', 'northeast', 'FontSize', 14)
title('Impurity (MLR)')
FigureTitle(2) = "MLR_LearningCurve_Tri_both";
ylim([10^(-2) 10^(-0)])
box on

%% Save figures
if bool_SaveFigures
    for i = 1:length(FigureTitle)
        saveas(figure(i), strcat(pwd,"\Images\LearningCurves\", ...
                                        FigureTitle(i)), 'png')
    end
end