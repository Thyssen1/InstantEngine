%% Initialization
% Clear workspace, close all figures and clear command window
clear; close all; clc

% Add path to experimental plan MATLAB function
addpath(strcat(erase(pwd, "Test scripts"), "\Functions"))

%
figurer;

%% Description


%% Main 
bool_SaveFigures = true;

SimulationIndex = 2; % {1 = no noise, 2 = noise}

% Set experimental design
params.design    = 'lhs';           % Box-Behnken design
params.setpoints = [20 11 20 2]; % [T pH Co lambda0 tdose]
params.lows      = [5 10 5 1];   % Low values
params.highs     = [40 12 50 3]; % High values
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

Kernels = {'exponential'; 'squaredexponential'; 'matern32';         ...
           'matern52'; 'rationalquadratic'; 'ardexponential';       ...
           'ardsquaredexponential'; 'ardmatern32'; 'ardmatern52';   ...
           'ardrationalquadratic'};

% Compare GPR models with different kernel functions
for i = 1:length(Kernels)
    % Step 1:
    % Compute hyperparameter Sigma and estimate of initial values 
    % of kernel parameters
    rng(10, 'twister')
    gprMdl_Yield{i} = fitrgp(X, yD,                                      ...
                    'KernelFunction', Kernels{i},  ...
                    'OptimizeHyperparameters', 'Sigma',         ...
                    'Standardize', true,                        ...
                    'HyperparameterOptimizationOptions',        ...
                    struct('ShowPlots', false,                  ...
                           'KFold', 5,                          ...
                           'MaxObjectiveEvaluations', 50));
    rng(11, 'twister')
    gprMdl_Tri{i} = fitrgp(X, yH,                           ...
                    'KernelFunction', Kernels{i},               ...
                    'OptimizeHyperparameters', 'Sigma',         ...
                    'Standardize', true,                        ...
                    'HyperparameterOptimizationOptions',        ...
                    struct('ShowPlots', false,                  ...
                           'KFold', 5,                          ...
                           'MaxObjectiveEvaluations', 50));
    
    % Step 2: 
    % Set initial parameters
    Sigma0_Yield(i,1)          = gprMdl_Yield{i}.Sigma;
    Sigma0_lb_Yield(i,1)       = min(1e-2*std(yD), gprMdl_Yield{i}.Sigma-1e-3);
    KernelParameters0_Yield{i} = gprMdl_Yield{i}.KernelInformation.KernelParameters;

    Sigma0_Tri(i,1)          = gprMdl_Tri{i}.Sigma;
    Sigma0_lb_Tri(i,1)       = min(1e-2*std(yD), gprMdl_Tri{i}.Sigma-1e-3);
    KernelParameters0_Tri{i} = gprMdl_Tri{i}.KernelInformation.KernelParameters;
    
    % Compute crossvalidation object 
%     rng(14, 'twister')
%     cvgprMdl_Yield{i} = fitrgp(X, yD,                               ...
%                     'KernelFunction', Kernels{i},                   ...
%                     'KernelParameters', KernelParameters0_Yield{i}, ...
%                     'Sigma', Sigma0_Yield(i),                       ...
%                     'SigmaLowerBound', Sigma0_lb_Yield(i),          ...
%                     'KFold', 5);
% 
%     rng(15, 'twister')
%     cvgprMdl_Tri{i} = fitrgp(X, yH,                               ...
%                     'KernelFunction', Kernels{i},                   ...
%                     'KernelParameters', KernelParameters0_Tri{i}, ...
%                     'Sigma', Sigma0_Tri(i),                       ...
%                     'SigmaLowerBound', Sigma0_lb_Tri(i),          ...
%                     'KFold', 5);
    
    % Step 3: Compute loss (mean squared error) of crossvalidation
%     Loss_Yield(i) = kfoldLoss(cvgprMdl_Yield{i});
%     Loss_Yield(i) = kfoldLoss(gprMdl_Yield{i}.crossval);
%     Loss_Tri(i)   = kfoldLoss(cvgprMdl_Tri{i});
%     Loss_Tri(i)   = kfoldLoss(gprMdl_Tri{i}.crossval);
    
    % Step 4: Fit models
    rng(12, 'twister')
    gprMdl_Yield_CV{i} = fitrgp(X, yD,                                   ...
                    'KernelFunction', Kernels{i},               ...
                    'KernelParameters', KernelParameters0_Yield{i},   ...
                    'Sigma', Sigma0_Yield(i),                         ...
                    'SigmaLowerBound', Sigma0_lb_Yield(i), 'KFold', 5);

    rng(13, 'twister')
    gprMdl_Tri_CV{i} = fitrgp(X, yH,                                   ...
                    'KernelFunction', Kernels{i},               ...
                    'KernelParameters', KernelParameters0_Tri{i},   ...
                    'Sigma', Sigma0_Tri(i),                         ...
                    'SigmaLowerBound', Sigma0_lb_Tri(i), 'KFold', 5);

    Loss_Yield(i) = kfoldLoss(gprMdl_Yield_CV{i});
% %     Loss_Tri(i)   = kfoldLoss(cvgprMdl_Tri{i});
    Loss_Tri(i)   = kfoldLoss(gprMdl_Tri_CV{i});

    % Fit models
%     rng(12, 'twister')
%     gprMdl_Yield{i} = fitrgp(X, yD,                                   ...
%                     'KernelFunction', Kernels{i},               ...
%                     'KernelParameters', KernelParameters0_Yield{i},   ...
%                     'Sigma', Sigma0_Yield(i),                         ...
%                     'SigmaLowerBound', Sigma0_lb_Yield(i));
% 
%     rng(13, 'twister')
%     gprMdl_Tri{i} = fitrgp(X, yH,                                   ...
%                     'KernelFunction', Kernels{i},               ...
%                     'KernelParameters', KernelParameters0_Tri{i},   ...
%                     'Sigma', Sigma0_Tri(i),                         ...
%                     'SigmaLowerBound', Sigma0_lb_Tri(i));

end

%% Display results
[~,YieldIdx] = min(Loss_Yield);
[~,ImpurityIdx] = min(Loss_Tri);

disp(strcat("Yield kernel: ", Kernels{YieldIdx}))
disp(strcat("Impurity kernel: ", Kernels{ImpurityIdx}))

disp('Loss Yield | Loss Tri')
disp([Loss_Yield' Loss_Tri'])

disp('Yield Kernel parameters:')
for i = 1:length(KernelParameters0_Yield{YieldIdx})
    disp(KernelParameters0_Yield{YieldIdx}(i))
end

disp('Impurity Kernel parameters:')
for i = 1:length(KernelParameters0_Tri{ImpurityIdx})
    disp(KernelParameters0_Tri{ImpurityIdx}(i))
end

disp('Sigma0_Yield:')
disp(Sigma0_Yield(YieldIdx))

disp('Sigma0_Tri:')
disp(Sigma0_Tri(YieldIdx))

disp('SigmaLowerBound_Yield:')
disp(Sigma0_lb_Yield(YieldIdx))

disp('SigmaLowerBound_Tri:')
disp(Sigma0_lb_Tri(ImpurityIdx))

%% Evaluate results on test sets of n=10 observations with LHS
% Create N different datasets each consisting of n=10 observations by
% creating test plans with different seeds

% Set experimental design
params.design    = 'lhs';           % Box-Behnken design
params.setpoints = [20 11 20 2];    % [T pH Co lambda0 tdose]
params.lows      = [5 10 5 1];      % Low values
params.highs     = [40 12 50 3];    % High values
params.n         = 4;               % Number of factors
params.samples   = 1000;

% Set standard deviations for InstantLab simulation
sigma_input.T       = 0;
sigma_input.pH      = 0;
sigma_input.Co      = 0;
sigma_input.CSdose  = 0;

plan = zeros(params.samples, 5);
yD = zeros(params.samples, length(Kernels));
yH = zeros(params.samples, length(Kernels));
for i = 1:length(Kernels)
%     params.seed = seeds(i);

    % Generate experimental plan
    plan(:,1:4) = exp_plan(params);
    plan(:,5) = 30*ones(params.samples,1);        

    % Generate synthetic experimental data with no noise
    data = instantlab(plan, sigma_input, 0);

    % Unwrap input and data from structure array
    for j = 1:size(data.out,2)    
        % Responses
        yD(j,i) = data.out{j}(6); % Yield
        yH(j,i) = data.out{j}(5); % Impurity
    end

    % Compute statistics
    stats_Yield(i) = rs_stats(yD(:,i), predict(gprMdl_Yield{i}, plan(:,1:4)));
    stats_Impurity(i) = rs_stats(yH(:,i), predict(gprMdl_Tri{i}, plan(:,1:4)));

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

%%
disp('R2_Yield')
disp(R2_Yield)
disp('SSE Yield')
disp(SSE_Yield)
disp('RMSE Yield')
disp(RMSE_Yield)
disp('AAD Yield:')
disp(AAD_Yield)

disp('R2 Impurity:')
disp(R2_Impurity)
disp('SSE Impurity:')
disp(SSE_Impurity)
disp('RMSE Impurity:')
disp(RMSE_Impurity)
disp('AAD Impurity:')
disp(AAD_Impurity)

% %% Ntest = 1000
% % Set experimental design
% params.design    = 'lhs';           % Box-Behnken design
% params.setpoints = [20 11 20 2];    % [T pH Co lambda0 tdose]
% params.lows      = [5 10 5 1];      % Low values
% params.highs     = [40 12 50 3];    % High values
% params.n         = 4;               % Number of factors
% params.samples   = 1000;
% 
% % Set standard deviations for InstantLab simulation
% sigma_input.T       = 0;
% sigma_input.pH      = 0;
% sigma_input.Co      = 0;
% sigma_input.CSdose  = 0;
% 
% plan_test2 = zeros(params.samples, 5);
% yD_test2 = zeros(params.samples, length(Kernels));
% yH_test2 = zeros(params.samples, length(Kernels));
% for i = 1:length(Kernels)
% %     params.seed = seeds(i);
% 
%     % Generate experimental plan
%     plan_test2(:,1:4) = exp_plan(params);
%     plan_test2(:,5) = 30*ones(params.samples,1);        
% 
%     % Generate synthetic experimental data with no noise
%     data = instantlab(plan_test2, sigma_input, 0);
% 
%     % Unwrap input and data from structure array
%     for j = 1:size(data.out,2)    
%         % Responses
%         yD_test2(j,i) = data.out{j}(6); % Yield
%         yH_test2(j,i) = data.out{j}(5); % Impurity
%     end
% 
%     % Compute statistics
%     stats_Yield_test2(i) = rs_stats(yD_test2(:,i), predict(gprMdl_Yield{i}, plan_test2(:,1:4)));
%     stats_Impurity_test2(i) = rs_stats(yH_test2(:,i), predict(gprMdl_Tri{i}, plan_test2(:,1:4)));
% 
%     R2_Yield_test2(i) = stats_Yield_test2(i).R2;
%     SSE_Yield_test2(i) = stats_Yield_test2(i).SSE;
%     MSE_Yield_test2(i) = stats_Yield_test2(i).MSE;
%     AAD_Yield_test2(i) = stats_Yield_test2(i).AAD;
% 
%     R2_Impurity_test2(i) = stats_Impurity_test2(i).R2;
%     SSE_Impurity_test2(i) = stats_Impurity_test2(i).SSE;
%     MSE_Impurity_test2(i) = stats_Impurity_test2(i).MSE;
%     AAD_Impurity_test2(i) = stats_Impurity_test2(i).AAD;
% end


%% Figure
switch SimulationIndex
    case 1
        str = " GPR LHS (noisefree)";
    case 2
        str = " GPR LHS (noisy)";
end

figure
cats = categorical(Kernels);
cats = reordercats(cats, Kernels);
bar(cats, [Loss_Yield' MSE_Yield'])
legend("CV estimate", "MSE measure", 'location', 'northeast')
ylabel('MSE')
title(strcat('Yield ', str))
FigureTitle(1) = "GPR_LHS_HyperparameterOptCurves_Yield";

figure
cats = categorical(Kernels);
cats = reordercats(cats, Kernels);
bar(cats, [Loss_Tri' MSE_Impurity'])
legend("CV estimate", "MSE measure", 'location', 'northeast')
ylabel('MSE')
title(strcat('Impurity', str))
FigureTitle(2) = "GPR_LHS_HyperparameterOptCurves_Impurity";

%% No noise
% Yield kernel: ardmatern32
% Impurity kernel: ardmatern32
% Loss Yield | Loss Tri
%     0.0588    0.0092
%     0.0584    0.0109
%     0.0590    0.0104
%     0.0592    0.0107
%     0.0584    0.0110
%     0.0130    0.0066
%     0.0088    0.0049
%     0.0046    0.0046
%     0.0087    0.0046
%     0.0087    0.0049
% 
% Yield Kernel parameters:
%    5.9793e+04
% 
%     4.6965
% 
%    17.6527
% 
%     4.0355
% 
%     0.4556
% 
% Impurity Kernel parameters:
%    4.1178e+04
% 
%     6.5216
% 
%    18.2146
% 
%     6.7013
% 
%     0.4025
% 
% Sigma0_Yield:
%     0.0024
% 
% Sigma0_Tri:
%     0.0010
% 
% SigmaLowerBound_Yield:
%     0.0014
% 
% SigmaLowerBound_Tri:
%    1.0000e-06
% 
% R2_Yield
%     0.8832    0.9059    0.9059    0.9092    0.9070    0.8952    0.9451    0.9411    0.9545    0.9526
% 
% SSE Yield
%     7.8632    5.8019    5.6487    5.4906    5.6927    6.1514    3.3632    3.6295    2.8258    2.9353
% 
% RMSE Yield
%     0.0887    0.0762    0.0752    0.0741    0.0755    0.0784    0.0580    0.0602    0.0532    0.0542
% 
% AAD Yield:
%     0.0714    0.0573    0.0570    0.0562    0.0570    0.0543    0.0419    0.0428    0.0363    0.0380
% 
% R2 Impurity:
%     0.7535    0.7783    0.8224    0.8259    0.7804    0.9268    0.9627    0.9418    0.9563    0.9627
% 
% SSE Impurity:
%     1.2423    1.1653    0.8941    0.8911    1.1536    0.4942    0.4358    0.6881    0.5319    0.4358
% 
% RMSE Impurity:
%     0.0352    0.0341    0.0299    0.0299    0.0340    0.0222    0.0209    0.0262    0.0231    0.0209
% 
% AAD Impurity:
%     0.0181    0.0199    0.0173    0.0174    0.0198    0.0122    0.0131    0.0144    0.0131    0.0131

%% Noise
% Yield kernel: ardmatern52
% Impurity kernel: ardmatern32
% Loss Yield | Loss Tri
%     0.0623    0.0091
%     0.0568    0.0108
%     0.0622    0.0104
%     0.0622    0.0107
%     0.0630    0.0109
%     0.0203    0.0063
%     0.0183    0.0094
%     0.0175    0.0061
%     0.0164    0.0080
%     0.0183    0.0080
% 
% Yield Kernel parameters:
%    16.0568
% 
%     3.0485
% 
%    3.2043e+04
% 
%     2.1668
% 
%     0.3910
% 
% Impurity Kernel parameters:
%    1.6269e+04
% 
%     4.1859
% 
%     8.3033
% 
%     3.9607
% 
%     0.2407
% 
% Sigma0_Yield:
%     0.0478
% 
% Sigma0_Tri:
%     0.0077
% 
% SigmaLowerBound_Yield:
%     0.0024
% 
% SigmaLowerBound_Tri:
%     0.0024
% 
% R2_Yield
%     0.8755    0.8173    0.8811    0.8668    0.8503    0.9161    0.8809    0.8877    0.8854    0.8861
% 
% SSE Yield
%     8.4416   10.7378    7.1145    7.8635    8.7770    5.4206    7.8456    7.0234    7.2606    7.4364
% 
% RMSE Yield
%     0.0919    0.1036    0.0843    0.0887    0.0937    0.0736    0.0886    0.0838    0.0852    0.0862
% 
% AAD Yield:
%     0.0750    0.0833    0.0670    0.0709    0.0753    0.0531    0.0658    0.0603    0.0618    0.0635
% 
% R2 Impurity:
%     0.7434    0.7154    0.7857    0.7670    0.7494    0.9222    0.9215    0.9438    0.9471    0.9257
% 
% SSE Impurity:
%     1.3014    1.4471    1.0736    1.1727    1.2704    0.4668    0.3872    0.5259    0.2653    0.3660
% 
% RMSE Impurity:
%     0.0361    0.0380    0.0328    0.0342    0.0356    0.0216    0.0197    0.0229    0.0163    0.0191
% 
% AAD Impurity:
%     0.0183    0.0223    0.0190    0.0198    0.0207    0.0115    0.0142    0.0129    0.0107    0.0133

%% Save figures
if bool_SaveFigures
    for i = 1:length(FigureTitle) 
        saveas(figure(i), ...
            strcat(pwd, "\Images\GPR_HyperparameterOptimization\", ...
                        FigureTitle(i), "_NoiseLevel_", num2str(SimulationIndex)), 'png')
    end
end


%% old notes

%% No noise
% Yield kernel: ardmatern32
% Impurity kernel: ardmatern32
% Yield Kernel parameters:
%    5.9793e+04
% 
%     4.6965
% 
%    17.6527
% 
%     4.0355
% 
%     0.4556
% 
% Impurity Kernel parameters:
%    4.1178e+04
% 
%     6.5216
% 
%    18.2146
% 
%     6.7013
% 
%     0.4025
% 
% Sigma0_Yield:
%     0.0024
% 
% Sigma0_Tri:
%     0.0010
% 
% SigmaLowerBound_Yield:
%     0.0014
% 
% SigmaLowerBound_Tri:
%    1.0000e-06

%% Noise
% Yield kernel: ardmatern52
% Impurity kernel: ardmatern32
% Yield Kernel parameters:
%    15.9405
% 
%     3.0409
% 
%    2.5854e+04
% 
%     2.1634
% 
%     0.3906
% 
% Impurity Kernel parameters:
%    1.6675e+04
% 
%     4.2150
% 
%     8.3821
% 
%     3.9889
% 
%     0.2420
% 
% Sigma0_Yield:
%     0.0473
% 
% Sigma0_Tri:
%     0.0077
% 
% SigmaLowerBound_Yield:
%     0.0024
% 
% SigmaLowerBound_Tri:
%     0.0024
% 
% R2_Yield
%     0.0465    0.0376    0.0438    0.0413    0.0376    0.8923    0.8738    0.8771    0.8771    0.8738
% 
% SSE Yield
%     0.5439    0.5433    0.5409    0.5416    0.5433    0.0742    0.0792    0.0771    0.0771    0.0792
% 
% RMSE Yield
%     0.2332    0.2331    0.2326    0.2327    0.2331    0.0861    0.0890    0.0878    0.0878    0.0890
% 
% AAD Yield:
% R2 Impurity:
%     0.0247    0.0016    0.0135    0.0091    0.0158    0.9182    0.8279    0.8627    0.8568    0.8279
% 
% SSE Impurity:
%     0.0479    0.0297    0.0371    0.0336    0.0354    0.0081    0.0011    0.0092    0.0012    0.0011
% 
% RMSE Impurity:
%     0.0692    0.0545    0.0609    0.0580    0.0595    0.0285    0.0104    0.0303    0.0109    0.0104
% 
% AAD Impurity:
%     0.0436    0.0302    0.0377    0.0352    0.0377    0.0156    0.0088    0.0170    0.0087    0.0088

