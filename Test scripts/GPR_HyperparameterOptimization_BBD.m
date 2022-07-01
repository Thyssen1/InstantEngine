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
params.design    = 'bbd';           % Box-Behnken design
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

%     rng(11, 'twister')
%     gprMdl_Yield{i} = fitrgp(X, yD,                                   ...
%                     'KernelFunction', Kernels{i},               ...
%                     'KernelParameters', KernelParameters0_Yield{i},   ...
%                     'Sigma', Sigma0_Yield(i),                         ...
%                     'SigmaLowerBound', Sigma0_lb_Yield(i));
% 
%     gprMdl_Tri{i} = fitrgp(X, yH,                                   ...
%                     'KernelFunction', Kernels{i},               ...
%                     'KernelParameters', KernelParameters0_Tri{i},   ...
%                     'Sigma', Sigma0_Tri(i),                         ...
%                     'SigmaLowerBound', Sigma0_lb_Tri(i));

    Loss_Yield(i) = kfoldLoss(gprMdl_Yield_CV{i});
% %     Loss_Tri(i)   = kfoldLoss(cvgprMdl_Tri{i});
    Loss_Tri(i)   = kfoldLoss(gprMdl_Tri_CV{i});

end

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
yD = zeros(params.samples, 1);
yH = zeros(params.samples, 1);
% Generate experimental plan
plan(:,1:4) = exp_plan(params);
plan(:,5) = 30*ones(params.samples,1);        

% Generate synthetic experimental data with no noise
data = instantlab(plan, sigma_input, 0);

% Unwrap input and data from structure array
for j = 1:size(data.out,2)    
    % Responses
    yD(j) = data.out{j}(6); % Yield
    yH(j) = data.out{j}(5); % Impurity
end

for i = 1:length(Kernels)
    % Compute statistics
    stats_Yield(i) = rs_stats(yD, predict(gprMdl_Yield{i}, plan(:,1:4)));
    stats_Impurity(i) = rs_stats(yH, predict(gprMdl_Tri{i}, plan(:,1:4)));

    R2_Yield(i) = stats_Yield(i).R2;
    SSE_Yield(i) = stats_Yield(i).SSE;
    MSE_Yield(i) = stats_Yield(i).MSE;
    AAD_Yield(i) = stats_Yield(i).AAD;

    R2_Impurity(i) = stats_Impurity(i).R2;
    SSE_Impurity(i) = stats_Impurity(i).SSE;
    MSE_Impurity(i) = stats_Impurity(i).MSE;
    AAD_Impurity(i) = stats_Impurity(i).AAD;
end


%% Display results
[~,YieldIdx] = min(Loss_Yield);
[~,ImpurityIdx] = min(Loss_Tri);

disp(strcat("Yield kernel: ", Kernels{YieldIdx}))
disp(strcat("Impurity kernel: ", Kernels{ImpurityIdx}))

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

disp('Loss_Yield Loss_Tri')
disp([Loss_Yield' Loss_Tri'])

disp('R2_Yield')
disp(R2_Yield)
disp('SSE Yield')
disp(SSE_Yield)
disp('MSE Yield')
disp(MSE_Yield)
disp('AAD Yield:')

disp('R2 Impurity:')
disp(R2_Impurity)
disp('SSE Impurity:')
disp(SSE_Impurity)
disp('MSE Impurity:')
disp(MSE_Impurity)
disp('AAD Impurity:')
disp(AAD_Impurity)

%% Figure
switch SimulationIndex
    case 1
        str = " GPR BBD (noisefree)";
    case 2
        str = " GPR BBD (noisy)";
end

figure
cats = categorical(Kernels);
cats = reordercats(cats, Kernels);
bar(cats, [Loss_Yield' MSE_Yield'])
legend("CV estimate", "MSE measure", 'location', 'northeast')
ylabel('MSE')
title(strcat('Yield ', str))
FigureTitle(1) = "GPR_BBD_HyperparameterOptCurves_Yield";

figure
cats = categorical(Kernels);
cats = reordercats(cats, Kernels);
bar(cats, [Loss_Tri' MSE_Impurity'])
legend("CV estimate", "MSE measure", 'location', 'northwest')
ylabel('MSE')
title(strcat('Impurity ', str))
FigureTitle(2) = "GPR_BBD_HyperparameterOptCurves_Impurity";

%% Save figures
if bool_SaveFigures
    for i = 1:length(FigureTitle) 
        saveas(figure(i), ...
            strcat(pwd, "\Images\GPR_HyperparameterOptimization\", ...
                        FigureTitle(i), "_NoiseLevel_", num2str(SimulationIndex)), 'png')
    end
end


%% No noise
% Yield kernel: ardmatern52
% Impurity kernel: ardmatern32
% Yield Kernel parameters:
%    14.7062
% 
%     3.5371
% 
%     4.6234
% 
%     1.1533
% 
%     0.3731
% 
% Impurity Kernel parameters:
%    20.7363
% 
%     2.3129
% 
%    13.2451
% 
%     1.4600
% 
%     0.1708
% 
% Sigma0_Yield:
%     0.0081
% 
% Sigma0_Tri:
%     0.0010
% 
% SigmaLowerBound_Yield:
%     0.0028
% 
% SigmaLowerBound_Tri:
%    1.0000e-06
% 
% Loss_Yield Loss_Tri
%     0.0859    0.0076
%     0.0859    0.0076
%     0.0859    0.0076
%     0.0859    0.0076
%     0.0859    0.0076
%     0.0112    0.0057
%     0.0474    0.0059
%     0.0076    0.0048
%     0.0073    0.0059
%     0.0451    0.0124
% 
% R2_Yield
%     0.9068    0.9344    0.9329    0.9357    0.9342    0.1400    0.9660    0.3401    0.9650    0.9656
% 
% SSE Yield
%     6.8740    3.8494    4.0144    3.7775    3.8558   71.2897    2.1793   78.9198    2.4288    2.1947
% 
% MSE Yield
%     0.0069    0.0038    0.0040    0.0038    0.0039    0.0713    0.0022    0.0789    0.0024    0.0022
% 
% AAD Yield:
% R2 Impurity:
%     0.8083    0.8300    0.8316    0.8337    0.8300    0.8468    0.9219    0.9182    0.9238    0.0353
% 
% SSE Impurity:
%     1.0292    0.9545    0.8348    0.8499    0.9545    0.8611    0.7510    0.6034    0.6403    4.8773
% 
% MSE Impurity:
%     0.0010    0.0010    0.0008    0.0008    0.0010    0.0009    0.0008    0.0006    0.0006    0.0049
% 
% AAD Impurity:
%     0.0170    0.0210    0.0170    0.0181    0.0210    0.0176    0.0178    0.0143    0.0150    0.0502

%% Noise
% Yield kernel: ardmatern32
% Impurity kernel: ardmatern32
% Yield Kernel parameters:
%    14.6785
% 
%     3.1564
% 
%     4.3146
% 
%     0.1596
% 
%     0.2975
% 
% Impurity Kernel parameters:
%    17.2477
% 
%     1.8816
% 
%    10.2581
% 
%     1.6491
% 
%     0.1574
% 
% Sigma0_Yield:
%     0.0632
% 
% Sigma0_Tri:
%     0.0010
% 
% SigmaLowerBound_Yield:
%     0.0029
% 
% SigmaLowerBound_Tri:
%    1.0000e-06
% 
% Loss_Yield Loss_Tri
%     0.0920    0.0078
%     0.0920    0.0078
%     0.0920    0.0078
%     0.0920    0.0078
%     0.0920    0.0078
%     0.0133    0.0060
%     0.0345    0.0055
%     0.0102    0.0053
%     0.0111    0.0057
%     0.0345    0.0110
% 
% R2_Yield
%     0.8561    0.8460    0.8764    0.8783    0.8454    0.1769    0.9220    0.3985    0.4286    0.9152
% 
% SSE Yield
%    10.5519   10.3193    8.0120    7.8377   10.3732   64.2948    6.7814   52.1527   50.2025    7.5956
% 
% MSE Yield
%     0.0106    0.0103    0.0080    0.0078    0.0104    0.0643    0.0068    0.0522    0.0502    0.0076
% 
% AAD Yield:
% R2 Impurity:
%     0.8025    0.8256    0.8282    0.8303    0.8256    0.8558    0.9078    0.9075    0.9114    0.2469
% 
% SSE Impurity:
%     1.0267    1.0135    0.8697    0.8947    1.0135    0.8684    0.9320    0.7531    0.7781    4.8754
% 
% MSE Impurity:
%     0.0010    0.0010    0.0009    0.0009    0.0010    0.0009    0.0009    0.0008    0.0008    0.0049
% 
% AAD Impurity:
%     0.0185    0.0224    0.0186    0.0197    0.0224    0.0165    0.0196    0.0157    0.0169    0.0582

