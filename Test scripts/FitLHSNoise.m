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
pars.setpoints = [20 11 20 2]; % [T pH Co lambda0]
pars.lows      = [5 10 5 1];   % Low values
pars.highs     = [40 12 50 3]; % High values
pars.reps      = 1;                             
pars.n         = 4;            % Number of factors
pars.seed      = 2;            % Original was 2
pars.samples   = 27;

%% Set Instant Lab parameters
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
YD = zeros(length(plan),1);
for i = 1:length(plan)
    YD(i,1) = data.out{i}(6); % Yield  
    YH(i,1) = data.out{i}(5); % Impurity

    y(i,:) = [data.out{i}(1) data.out{i}(2) data.out{i}(3) ...
                  data.out{i}(4) data.out{i}(5) data.out{i}(6) ...
                  data.out{i}(7) data.out{i}(8) data.out{i}(9) ...
                  data.out{i}(10) data.out{i}(11) data.out{i}(12) ...
                  data.out{i}(13)];
end

% Design matrix (tdose not included for DOE, GPR and ANN)
X = plan(:,1:4);

%% Fit mechanistic model (tdose included)
% Retrieve true parameters
reacstruc = reacstruccreate();
params = [reacstruc.model.kAref;    ...
          reacstruc.model.sref(3);  ...
          reacstruc.model.sdref;    ...
          reacstruc.model.pkA1;     ...
          reacstruc.model.pkA3;     ...
          reacstruc.model.EA1;      ...
          reacstruc.model.EA3;      ...
          reacstruc.model.EAd];
params = params(:);
clear reacstruc

% Fit model
reacstruc = reacfit_all(plan, y);
ParameterMinimum = [reacstruc.model.kAref;    ...
                    reacstruc.model.sref(3);  ...
                    reacstruc.model.sdref;     ...
                    reacstruc.model.pkA1;      ...
                    reacstruc.model.pkA3;      ...
                    reacstruc.model.EA1;       ...
                    reacstruc.model.EA3;       ...
                    reacstruc.model.EAd];
disp(abs(params - ParameterMinimum) ./ params)

%% Multiple linear regression (MLR) model
% Fit quadratic (2nd order) model
MLR_Yield    = fitlm(X, YD, 'quadratic');
MLR_Impurity = fitlm(X, YH, 'quadratic');

%% Fit Gaussian Process Regression (GPR) model
% Set kernel function
switch SimulationIndex
    case 1
        Kernel_Yield = 'ardmatern32';
        Kernel_Impurity = 'ardmatern32';
    case 2
        Kernel_Yield = 'ardmatern52'; % ardmatern52
%         Kernel_Yield = 'ardsquaredexponential';
        Kernel_Impurity = 'ardmatern32'; %armatern32
%         Kernel_Impuirty = 'ardsquaredexponential';
end

rng(10, 'twister')
% Step 1: 
% Optimize hyperparameters and compute initial estimate of kernel parameters
gprMdl_Yield = fitrgp(X, YD,              ...
        'KernelFunction', Kernel_Yield,      ...
        'OptimizeHyperparameters', 'Sigma',  ...
        'Standardize', true,                 ...
        'HyperparameterOptimizationOptions', ...
        struct('ShowPlots', false,           ...
               'KFold', 5,                   ...
               'MaxObjectiveEvaluations', 50));

rng(11, 'twister')
gprMdl_Tri = fitrgp(X, YH,                      ...    
            'KernelFunction', Kernel_Impurity,      ...
            'OptimizeHyperparameters', 'Sigma',     ...   
            'Standardize', true,                    ...
            'HyperparameterOptimizationOptions',    ...
            struct('ShowPlots', false,              ...
                   'KFold', 5,                      ...
                   'MaxObjectiveEvaluations', 50));

% % Step 2: 
% % Set initial parameters
% Sigma0_Yield            = gprMdlYield_Init.Sigma;
% Sigma0_lb_Yield         = min(1e-2*std(YD), gprMdlYield_Init.Sigma-1e-3);
% KernelParameters0_Yield = gprMdlYield_Init.KernelInformation.KernelParameters;
% 
% Sigma0_Tri            = gprMdlTri_Init.Sigma;
% Sigma0_lb_Tri         = min(1e-2*std(YH), gprMdlTri_Init.Sigma-1e-3);
% KernelParameters0_Tri = gprMdlTri_Init.KernelInformation.KernelParameters;
% 
% 
% % Step 3: Fit models 
% rng(12, 'twister')
% gprMdl_Yield = fitrgp(X, YD,                                        ...
%                     'KernelFunction', Kernel_Yield,                ...
%                     'KernelParameters', KernelParameters0_Yield,    ...
%                     'Sigma', Sigma0_Yield,                          ...
%                     'SigmaLowerBound', Sigma0_lb_Yield,             ...
%                     'Standardize', true,                            ...
%                     'ConstantSigma', false);
% rng(13, 'twister')
% gprMdl_Tri = fitrgp(X, YH,                                      ...
%                     'KernelFunction', Kernel_Impurity,            ...
%                     'KernelParameters', KernelParameters0_Tri,  ...
%                     'Sigma', Sigma0_Tri,                        ...
%                     'SigmaLowerBound', Sigma0_lb_Tri,           ...
%                     'Standardize', true,                        ...
%                     'ConstantSigma', false);

%% Fit Artificial Neural Network (ANN) model
% Set hyperparameter objects
HyperParams_Yield = hyperparameters("fitrnet", X, YD);
HyperParams_Tri   = hyperparameters("fitrnet", X, YH);

% First object
HyperParams_Yield(1).Optimize = false;  % Number of layers
HyperParams_Yield(2).Optimize = false;  % Activation function
HyperParams_Yield(3).Optimize = false;  % Standardize
HyperParams_Yield(4).Optimize = true;   % Regularization strength
HyperParams_Yield(5).Optimize = false;   % Layer Weights Initializer
HyperParams_Yield(6).Optimize = false;   % Layer Biases Initializer
HyperParams_Yield(7).Optimize = false;  % Number of neurons layer 1
HyperParams_Yield(8).Optimize = false;  % Number of neurons layer 2
HyperParams_Yield(9).Optimize = false;  % Number of neurons layer 3
HyperParams_Yield(10).Optimize = false; % Number of neurons layer 4
HyperParams_Yield(11).Optimize = false; % Number of neurons layer 5

% Second object
HyperParams_Tri(1).Optimize = false;  % Number of layers
HyperParams_Tri(2).Optimize = false;  % Activation function
HyperParams_Tri(3).Optimize = false;  % Standardize
HyperParams_Tri(4).Optimize = true;   % Regularization strength
HyperParams_Tri(5).Optimize = false;  % Layer Weights Initializer
HyperParams_Tri(6).Optimize = false;  % Layer Biases Initializer
HyperParams_Tri(7).Optimize = false;   % Number of neurons layer 1
HyperParams_Tri(8).Optimize = false;   % Number of neurons layer 2
HyperParams_Tri(9).Optimize = false;  % Number of neurons layer 3
HyperParams_Tri(10).Optimize = false; % Number of neurons layer 4
HyperParams_Tri(11).Optimize = false; % Number of neurons layer 5

Neurons_Yield = [40 40];
Neurons_Impurity = [40 40];

% Set activation functions
switch SimulationIndex
    case 1
        Activation_Yield = 'tanh';
        Activation_Impurity = 'sigmoid';

        Activation_Yield = 'relu';
        Activation_Impurity = 'sigmoid';

        Neurons_Yield = [50 47];
        Neurons_Impurity = [50 10];

%         Regularization_Yield = 1.8538e-06;
%         Regularization_Impurity = 2.4194e-06;

    case 2
        Activation_Yield = 'tanh'; %relu
        Activation_Impurity = 'sigmoid'; %relu

        Activation_Yield = 'tanh';
        Activation_Sigmoid = 'sigmoid';

        Neurons_Yield = [5 16];
        Neurons_Impurity = [33 7];

%         Regularization_Yield = 0.0110;
%         Regularization_Impurity = 5.7982e-05;
end

rng(11, 'twister')
ANN_Yield = fitrnet(X, YD,                         ...
                'LayerSizes', Neurons_Yield,            ...  
                'Activations', Activation_Yield,        ...  
                'Standardize', true,                    ...
                'OptimizeHyperparameters', HyperParams_Yield,     ...
                'HyperParameterOptimizationOptions',    ...
                struct("AcquisitionFunctionName",       ...
                      "expected-improvement-plus",      ...
                      "MaxObjectiveEvaluations", 100,    ...
                      'ShowPlots', false, 'KFold', 5));

rng(12, 'twister')
ANN_Impurity = fitrnet(X, YH,                      ...
                'LayerSizes', Neurons_Impurity,         ...   
                'Activations', Activation_Impurity,     ...
                'Standardize', true,                    ...
                'OptimizeHyperparameters', HyperParams_Tri,     ...
                'HyperParameterOptimizationOptions',    ...
                struct("AcquisitionFunctionName",       ...
                      "expected-improvement-plus",      ...
                      "MaxObjectiveEvaluations", 100,    ...
                'ShowPlots', false, 'KFold', 5));

% rng(13, 'twister')
% ANN_Yield = fitrnet(X, YD,   ...  
%     "LayerSizes", Neurons_Yield,   ...
%     "Standardize", true,     ...
%     'Activations', Activation_Yield,   ...
%     'Lambda', ANN_Yield_Init.ModelParameters.Lambda);
% rng(13, 'twister')
% ANN_Yield = fitrnet(X, YD,   ...  
%     "LayerSizes", Neurons_Yield,   ...
%     "Standardize", true,     ...
%     'Activations', Activation_Yield,   ...
%     'Lambda', Regularization_Yield);
% 
% % rng(14, 'twister')
% % ANN_Impurity = fitrnet(X, YH,   ...  
% %     "LayerSizes", Neurons_Impurity,      ...
% %     "Standardize", true,        ...
% %     'Activations', Activation_Impurity,   ...
% %     'Lambda', ANN_Impurity_Init.ModelParameters.Lambda);
% rng(14, 'twister')
% ANN_Impurity = fitrnet(X, YH,   ...  
%     "LayerSizes", Neurons_Impurity,      ...
%     "Standardize", true,        ...
%     'Activations', Activation_Impurity,   ...
%     'Lambda', Regularization_Impurity);

%% Predictions (Mechanistic model)
% Run simulations for mechanistic model fit
for k = 1:size(plan,1)
    reacstruc.process.T       = plan(k,1); % Celcius scale
    reacstruc.process.pH      = plan(k,2); % Logarithmic pH scale
    reacstruc.process.Co      = plan(k,3); % g/L
    reacstruc.process.lambda0 = plan(k,4); % mol /mol
    reacstruc.process.tdose   = plan(k,5); % g/L

    % Run sim with new parameters
    reacstruc = reacsim(reacstruc);

    % Extract yield and triacylated components
    yD(k,1) = reacstruc.out.y(end,7);
    yH(k,1) = reacstruc.out.y(end,10);
%     yD_act(k,1) = data.out{k}(6);
%     yH_act(k,1) = data.out{k}(5);
end

% Set optim options
options = optimset('display', 'iter', 'tolfun', 1.0e-06,   ...
                   'tolx', 1.0e-5, 'maxfunevals', 0); % No iteration
               
% Compute estimate of Jacobian
% Update reacstruc structure array with current parameter estimates
for i = 1:length(reacstruc.parfit.par)
    LB = reacstruc.parfit.LB(i);
    UB = reacstruc.parfit.UB(i);
    
    % Scale back the initial guess
    ParMin0(i) = (ParameterMinimum(i) - LB) / (UB-LB);
end

[~,~,residualD,~,~,~,JacobiD] = lsqnonlin(@(x) reacobj_yd(x,reacstruc), ...
                                        ParameterMinimum, [], [], options);

% The Jacobian is reported in a column vector format, 
% which is reformated below to a matrix format
JacobianD = []; JacobianD(:,:) = JacobiD;

nP = length(residualD);
pP = length(ParameterMinimum);
DegreesFreedomD = nP - pP;

SquaredSumErrorD         = residualD'*residualD;
VarianceD                = SquaredSumErrorD/DegreesFreedomD; % Variance of errors. 
CovarianceParametersD    = VarianceD*inv(JacobianD'*JacobianD);  % Covariance of parameter estimators
StandardDeviationD       = sqrt(diag(CovarianceParametersD))'; % Standard deviation of parameter estimators
CorrelationParametersD   = CovarianceParametersD ./ [StandardDeviationD'*StandardDeviationD]; % Correlation matrix for the parameters


alpha = 0.01;  % Significance level 
tcr   = tinv((1-alpha),DegreesFreedomD); % Critical t-dist value at alpha 

%+-95% confidence intervals
ConfidenceIntervalP95D = [ParameterMinimum' - StandardDeviationD*tcr; ...
                          ParameterMinimum' + StandardDeviationD*tcr];

% Calculate confidence intervals on the model output
n = pars.samples;
m = 8;

% Compute output covariance
ModelCovariance                 = JacobianD * CovarianceParametersD * JacobianD'; % Calculate model covariance # Maybe some statistical reference here?
ModelStandardDeviation          = sqrt(diag(ModelCovariance));          % Standard deviation of model outputs, vector with all data points
ModelStandardDeviationReshaped  = reshape(ModelStandardDeviation,n,m);  % Reshape to matrix with columns corresponding to each set of data points

% 95% confidence intervals, calculated only for the measured species (data) Two boundaries, each with a set of m columns 
ModelConfidenceInterval95D = [yD - ModelStandardDeviationReshaped(:,4)*tcr,...  % Lower boundary
                              yD + ModelStandardDeviationReshaped(:,4)*tcr];      % Upper boundary 
ModelConfidenceInterval95Tri = [yH - ModelStandardDeviationReshaped(:,8)*tcr,...  % Lower boundary
                              yH + ModelStandardDeviationReshaped(:,8)*tcr];      % Upper boundary 

disp(ModelConfidenceInterval95D)
disp(ModelConfidenceInterval95Tri)
%% Old

% 
% 
% % Compute estimate of Jacobian
% [~,~,residualTri,~,~,~,JacobiTri] = lsqnonlin(@(x) reacobj_yh(x,reacstruc), ...
%                                         ParameterMinimum, [], [], options);
% 
% % The Jacobian is reported in a column vector format, 
% % which is reformated below to a matrix format
% JacobianD = []; JacobianD(:,:) = JacobiD;
% JacobianTri = []; JacobianTri(:,:) = JacobiTri;
% 
% % Degrees of freedom calculation
% nD = length(residualD);
% pD = length(ParameterMinimum);
% DegreesFreedomD = nD - pD;
% nTri = length(residualTri);
% pTri = pD;
% DegreesFreedomTri = nTri - pTri;
% 
% % Statistics
% SquaredSumErrorD         = residualD'*residualD;
% VarianceD                = SquaredSumErrorD/DegreesFreedomD; % Variance of errors. 
% CovarianceParametersD    = VarianceD*inv(JacobianD'*JacobianD);  % Covariance of parameter estimators
% StandardDeviationD       = sqrt(diag(CovarianceParametersD))'; % Standard deviation of parameter estimators
% CorrelationParametersD   = CovarianceParametersD ./ [StandardDeviationD'*StandardDeviationD]; % Correlation matrix for the parameters
% SquaredSumErrorTri       = residualTri'*residualTri;
% VarianceTri              = SquaredSumErrorTri / DegreesFreedomTri; % Variance of errors. 
% CovarianceParametersTri  = VarianceTri*inv(JacobianTri'*JacobianTri);  % Covariance of parameter estimators
% StandardDeviationTri     = sqrt(diag(CovarianceParametersTri))'; % Standard deviation of parameter estimators
% CorrelationParametersTri = CovarianceParametersTri ...
%                         ./ [StandardDeviationTri'*StandardDeviationTri]; % Correlation matrix for the parameters
% 
% 
% alpha = 0.025;  % Significance level 
% tcrD   = tinv((1-alpha),DegreesFreedomD); % Critical t-dist value at alpha 
% tcrTri   = tinv((1-alpha),DegreesFreedomTri); % Critical t-dist value at alpha 
% 
% %+-95% confidence intervals
% ConfidenceIntervalP95D = [ParameterMinimum' - StandardDeviationD*tcrD; ...
%                           ParameterMinimum' + StandardDeviationD*tcrD];
% ConfidenceIntervalP95Tri = [ParameterMinimum' - StandardDeviationTri*tcrTri; ...
%                             ParameterMinimum' + StandardDeviationTri*tcrTri];
% 
% % Calculate confidence intervals on the model output
% [nD, mD] = size(yD);
% [nTri, mTri] = size(yH);
% 
% % Compute output covariance
% ModelCovarianceD                 = JacobianD * CovarianceParametersD * JacobianD'; % Calculate model covariance # Maybe some statistical reference here?
% ModelStandardDeviationD          = sqrt(diag(ModelCovarianceD));          % Standard deviation of model outputs, vector with all data points
% ModelStandardDeviationReshapedD  = reshape(ModelStandardDeviationD,nD,mD);  % Reshape to matrix with columns corresponding to each set of data points
% ModelCovarianceTri                 = JacobianTri * CovarianceParametersTri * JacobianTri'; % Calculate model covariance # Maybe some statistical reference here?
% ModelStandardDeviationTri          = sqrt(diag(ModelCovarianceTri));          % Standard deviation of model outputs, vector with all data points
% ModelStandardDeviationReshapedTri  = reshape(ModelStandardDeviationTri,nTri,mTri);  % Reshape to matrix with columns corresponding to each set of data points
% 
% % 95% confidence intervals, calculated only for the measured species (data) Two boundaries, each with a set of m columns 
% ModelConfidenceInterval95D = [yD - ModelStandardDeviationReshapedD*tcrD,...  % Lower boundary
%                              yD + ModelStandardDeviationReshapedD*tcrD];      % Upper boundary 
% ModelConfidenceInterval95Tri = [yH - ModelStandardDeviationReshapedTri*tcrTri,...  % Lower boundary
%                              yH + ModelStandardDeviationReshapedTri*tcrTri];      % Upper boundary 

%% Compute predicted response and 95% confidence interval for MLR model
[MLR_yfit_Yield, MLR_yint_Yield] =  predict(MLR_Yield, X, 'alpha', alpha);
[MLR_yfit_Impurity, MLR_yint_Impurity] =  predict(MLR_Impurity, X, 'alpha', alpha);

%% Compute predicted response and 95% confidence intervals for GPR model
[GPR_yfit_Yield, ~, GPR_yint_Yield] = predict(gprMdl_Yield, X, 'alpha', alpha);
[GPR_yfit_Impurity, ~, GPR_yint_Impurity] = predict(gprMdl_Tri, X, 'alpha', alpha);
% [yfit,~,yint(:,:,j)] = resubPredict(gprMdl,'Alpha',0.01);

%% Compute predicted response for ANN model 
ANN_yfit_Yield = predict(ANN_Yield, X);  
ANN_yfit_Impurity = predict(ANN_Impurity, X);  

%% Compute model statistics
% Compute (mechanistic)
stats_MM_Yield = rs_stats(YD, yD);
stats_MLR_Yield = rs_stats(YD, MLR_yfit_Yield);
stats_GPR_Yield = rs_stats(YD, GPR_yfit_Yield);
stats_ANN_Yield = rs_stats(YD, ANN_yfit_Yield);

stats_MM_Impurity = rs_stats(YH, yH);
stats_MLR_Impurity = rs_stats(YH, MLR_yfit_Impurity);
stats_GPR_Impurity = rs_stats(YH, GPR_yfit_Impurity);
stats_ANN_Impurity = rs_stats(YH, ANN_yfit_Impurity);

R2_Yield = [stats_MM_Yield.R2; stats_MLR_Yield.R2;       ...
            stats_GPR_Yield.R2; stats_ANN_Yield.R2];
SSE_Yield = [stats_MM_Yield.SSE; stats_MLR_Yield.SSE;    ...
            stats_GPR_Yield.SSE; stats_ANN_Yield.SSE];  
RMSE_Yield = [stats_MM_Yield.RMSE; stats_MLR_Yield.RMSE; ...
            stats_GPR_Yield.RMSE; stats_ANN_Yield.RMSE];

R2_Impurity = [stats_MM_Impurity.R2; stats_MLR_Impurity.R2;       ...
            stats_GPR_Impurity.R2; stats_ANN_Impurity.R2];
SSE_Impurity = [stats_MM_Impurity.SSE; stats_MLR_Impurity.SSE;    ...
            stats_GPR_Impurity.SSE; stats_ANN_Impurity.SSE];  
RMSE_Impurity = [stats_MM_Impurity.RMSE; stats_MLR_Impurity.RMSE; ...
            stats_GPR_Impurity.RMSE; stats_ANN_Impurity.RMSE];

%% Figures of fit to training data (Yield)
[~, IY]  = sort(YD);
[~, IYr] = sort(YD, 'descend');

figure
subplot(2, 2, 1); hold all % Mechanistic 
plot(-1:0.01:2, -1:0.01:2, 'k-')
plot(YD, yD, 'x')
xlim([-0.30 1.35])
ylim([-0.30 1.35])
xlabel('Measured yield ($\%$)', 'FontSize', 10)
ylabel('Yield ($\%$)', 'FontSize', 10)
title('Mechanistic', 'FontSize', 10)
patch([YD(IY)' YD(IYr)'],                     ...
      [ModelConfidenceInterval95D(IY,1)'      ...
       ModelConfidenceInterval95D(IYr,2)'],   ...
      'blue', 'FaceAlpha', .125)
set(gca, 'FontSize', 11)
legend("", "Predicted", "95\% CI", ...
        'location', 'northwest', 'FontSize', 8, 'box', 'on')
text(0.5,-0.05, strcat("$R^{2}$: ", num2str(round(R2_Yield(1),4)*100), "\%"), 'FontSize', 10)
text(0.5,-0.15, strcat("RMSE: ", num2str(round(RMSE_Yield(1),4))), 'FontSize', 10)

subplot(2, 2, 2); hold all  % DOE (MLR)
plot(-1:0.01:2, -1:0.01:2, 'k-')
plot(YD, MLR_yfit_Yield, 'x')
patch([YD(IY)' YD(IYr)'], [MLR_yint_Yield(IY,1)' MLR_yint_Yield(IYr,2)'], ...
      'blue', 'FaceAlpha', .125) 
xlabel('Measured yield ($\%$)', 'FontSize', 10)
ylabel('Yield ($\%$)', 'FontSize', 10)
xlim([-0.30 1.35])
ylim([-0.30 1.35])
title('MLR', 'FontSize', 10)
set(gca, 'FontSize', 11)
legend("", "Predicted", "95 \% CI", ...
       'FontSize', 8, 'location', 'northwest', 'box', 'on')
text(0.5,-0.05, strcat("$R^{2}$: ", num2str(round(R2_Yield(2),4)*100), "\%"), 'FontSize', 10)
text(0.5,-0.15, strcat("RMSE: ", num2str(round(RMSE_Yield(2),4))), 'FontSize', 10)

subplot(2, 2, 3); hold all % GPR 
plot(-1:0.01:2, -1:0.01:2, 'k-')
plot(YD, GPR_yfit_Yield, 'x')
patch([YD(IY)' YD(IYr)'], [GPR_yint_Yield(IY,1)' GPR_yint_Yield(IYr,2)'], ...
      'blue', 'FaceAlpha', .125) 
xlabel('Measured yield ($\%$)', 'FontSize', 10)
ylabel('Yield ($\%$)', 'FontSize', 10)
xlim([-0.30 1.35])
ylim([-0.30 1.35])
title('GPR', 'FontSize', 10)
set(gca, 'FontSize', 11)
legend("", "Predicted", "95 \% CI", ...
       'location', 'northwest', 'FontSize', 8, 'box', 'on')
text(0.5,-0.05, strcat("$R^{2}$: ", num2str(round(R2_Yield(3),4)*100), "\%"), 'FontSize', 10)
text(0.5,-0.15, strcat("RMSE: ", num2str(round(RMSE_Yield(3),4))), 'FontSize', 10)   

subplot(2, 2, 4); hold all % ANN
plot(-1:0.01:2, -1:0.01:2, 'k-')
plot(YD, ANN_yfit_Yield, 'x')
xlabel('Measured yield ($\%$)', 'FontSize', 10)
ylabel('Yield ($\%$)', 'FontSize', 10)
xlim([-0.30 1.35])
ylim([-0.30 1.35])
title('ANN', 'FontSize', 10)
set(gca, 'FontSize', 11)
legend("", "Predicted", ...
       'location', 'northwest', 'FontSize', 8, 'box', 'on')
FigureTitle(1) = "FitLHS_TrainingFit_Yield";
sgtitle('LHS - Training')
text(0.5,-0.05, strcat("$R^{2}$: ", num2str(round(R2_Yield(4),4)*100), "\%"), 'FontSize', 10)
text(0.5,-0.15, strcat("RMSE: ", num2str(round(RMSE_Yield(4),4))), 'FontSize', 10)

%% Figures of fit to training data (Impurity)
[~, IH]  = sort(YH);
[~, IHr] = sort(YH, 'descend');

figure
subplot(2, 2, 1); hold all % Mechanistic 
plot(-1:0.01:2, -1:0.01:2, 'k-')
plot(YH, yH, 'x')
ylim([-0.15 0.50])
xlim([-0.15 0.50])
xlabel('Measured impurity ($\%$)', 'FontSize', 10)
ylabel('Impurity ($\%$)', 'FontSize', 10)
title('Mechanistic', 'FontSize', 10)
patch([YH(IH)' YH(IHr)'],                     ...
      [ModelConfidenceInterval95Tri(IH,1)'      ...
       ModelConfidenceInterval95Tri(IHr,2)'],   ...
      'blue', 'FaceAlpha', .125)
set(gca, 'FontSize', 11)
legend("", "Predicted", "95\% CI", ...
        'location', 'northwest', 'FontSize', 8, 'box', 'on')
text(0.2,-0.05, strcat("$R^{2}$: ", num2str(round(R2_Impurity(1),4)*100), "\%"), 'FontSize', 10)
text(0.2,-0.10, strcat("RMSE: ", num2str(round(RMSE_Impurity(1),4))), 'FontSize', 10)

subplot(2, 2, 2); hold all  % DOE (MLR)
plot(-1:0.01:2, -1:0.01:2, 'k-')
plot(YH, MLR_yfit_Impurity, 'x')
patch([YH(IH)' YH(IHr)'], [MLR_yint_Impurity(IH,1)' MLR_yint_Impurity(IHr,2)'], ...
      'blue', 'FaceAlpha', .125) 
ylim([-0.15 0.50])
xlim([-0.15 0.50])
title('MLR', 'FontSize', 10)
xlabel('Measured impurity ($\%$)', 'FontSize', 14)
ylabel('Impurity ($\%$)', 'FontSize', 14)
set(gca, 'FontSize', 11)
legend("", "Predicted", "95 \% CI", ...
       'FontSize', 8, 'location', 'northwest', 'box', 'on')
text(0.2,-0.05, strcat("$R^{2}$: ", num2str(round(R2_Impurity(2),4)*100), "\%"), 'FontSize', 10)
text(0.2,-0.10, strcat("RMSE: ", num2str(round(RMSE_Impurity(2),4))), 'FontSize', 10)

subplot(2, 2, 3); hold all % GPR 
plot(-1:0.01:2, -1:0.01:2, 'k-')
plot(YH, GPR_yfit_Impurity, 'x')
patch([YH(IH)' YH(IHr)'], [GPR_yint_Impurity(IH,1)' GPR_yint_Impurity(IHr,2)'], ...
      'blue', 'FaceAlpha', .125) 
xlim([min(min(min([YH GPR_yfit_Impurity GPR_yint_Impurity])),0) ...
      min(max(max([YH GPR_yfit_Impurity GPR_yint_Impurity])),1)])
ylim([min(min(min([YH GPR_yfit_Impurity GPR_yint_Impurity])),0) ...
      min(max(max([YH GPR_yfit_Impurity GPR_yint_Impurity])),1)])
ylim([-0.15 0.50])
xlim([-0.15 0.50])
title('GPR', 'FontSize', 10)
xlabel('Measured impurity ($\%$)', 'FontSize', 10)
ylabel('Impurity ($\%$)', 'FontSize', 10)
set(gca, 'FontSize', 11)
legend("", "Predicted", "95 \% CI", ...
       'location', 'northwest', 'FontSize', 8, 'box', 'on')
text(0.2,-0.05, strcat("$R^{2}$: ", num2str(round(R2_Impurity(3),4)*100), "\%"), 'FontSize', 10)
text(0.2,-0.10, strcat("RMSE: ", num2str(round(RMSE_Impurity(3),4))), 'FontSize', 10)
   
subplot(2, 2, 4); hold all % ANN
plot(-1:0.01:2, -1:0.01:2, 'k-')
plot(YH, ANN_yfit_Impurity, 'x')
xlim([min(min(min([YH ANN_yfit_Impurity])),0) ...
      min(max(max([YH ANN_yfit_Impurity])),1)])
ylim([min(min(min([YH ANN_yfit_Impurity])),0) ...
      min(max(max([YH ANN_yfit_Impurity])),1)])
ylim([-0.15 0.50])
xlim([-0.15 0.50])
xlabel('Measured impurity ($\%$)', 'FontSize', 10)
ylabel('Impurity ($\%$)', 'FontSize', 10)
title('ANN', 'FontSize', 10)
set(gca, 'FontSize', 11)
legend("", "Predicted", ...
       'location', 'northwest', 'FontSize', 8, 'box', 'on')
FigureTitle(2) = "FitLHS_TrainingFit_Tri";
text(0.2,-0.05, strcat("$R^{2}$: ", num2str(round(R2_Impurity(4),4)*100), "\%"), 'FontSize', 10)
text(0.2,-0.10, strcat("RMSE: ", num2str(round(RMSE_Impurity(4),4))), 'FontSize', 10)
sgtitle('LHS - Training')

%% Model evaluation (MLR)
MLR_Residual_YD = YD - MLR_yfit_Yield;
MLR_Residual_YH = YH - MLR_yfit_Impurity;

figure
subplot(2,2,1)
plot(MLR_yfit_Yield, MLR_Residual_YD,'o')
set(gca, 'FontSize', 8)
xlabel('Predicted yield, $y_{D}$', 'interpreter', 'latex')
ylabel('Residual')

subplot(2,2,2); hold all
histogram(MLR_Residual_YD, 10, 'normalization', 'pdf')
plot(min(MLR_Residual_YD):0.001:max(MLR_Residual_YD), ...
    normpdf(min(MLR_Residual_YD):0.001:max(MLR_Residual_YD), 0, ...
    std(MLR_Residual_YD)), 'LineWidth', 1)
set(gca, 'FontSize', 8)
xlabel('Residual')
ylabel('Frequency')

subplot(2,2,3)
qqplot(MLR_Residual_YD)
set(gca, 'FontSize', 8)

subplot(2,2,4)
plot(MLR_yfit_Yield, abs(MLR_Residual_YD),'o')
set(gca, 'FontSize', 8)
xlabel('Predicted yield, $y_{D}$', 'interpreter', 'latex')
ylabel('Abs(Residual)')
FigureTitle(3) = "FitLHS_TrainingResiduals_Yield";

figure
subplot(2,2,1)
plot(MLR_yfit_Impurity, MLR_Residual_YH,'o')
set(gca, 'FontSize', 8)
xlabel('Predicted impurity, $y_{H}$', 'interpreter', 'latex')
ylabel('Residual')

subplot(2,2,2); hold all
histogram(MLR_Residual_YH, 10, 'normalization', 'pdf')
plot(min(MLR_Residual_YH):0.001:max(MLR_Residual_YH), ...
    normpdf(min(MLR_Residual_YH):0.001:max(MLR_Residual_YH), 0, ...
    std(MLR_Residual_YH)), 'LineWidth', 1)
set(gca, 'FontSize', 11)
xlabel('Residual')
ylabel('Frequency')

subplot(2,2,3)
qqplot(MLR_Residual_YH)
set(gca, 'FontSize', 8)

subplot(2,2,4)
plot(MLR_yfit_Impurity, abs(MLR_Residual_YH),'o')
set(gca, 'FontSize', 8)
xlabel('Predicted impurity, $y_{H}$', 'interpreter', 'latex')
ylabel('Abs(Residual)')
FigureTitle(4) = "FitLHS_TrainingResiduals_Impurity";

%% Evaluation of models on test data
testplan.setpoints = [20 11 20 2]; % [T pH Co lambda0 tdose]
testplan.lows      = [5 10 5 1];   % Low values
testplan.highs     = [40 12 50 3]; % High values
Ntest              = 15; %15

% Generate test data
rng(20, 'twister')
for i = 1:length(testplan.setpoints)
    TestPlan(:,i) = icdf('Uniform', lhsdesign(Ntest,1), ...
                         testplan.lows(i), testplan.highs(i));
end
% Add tdose as a constant to the test plan
TestPlan = [TestPlan 30*ones(Ntest,1)];

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

%% Predictions with Mechanistic model on test set
% Run simulations for mechanistic model fit
for k = 1:size(TestPlan,1)
    reacstruc.process.T       = TestPlan(k,1); % Celcius
    reacstruc.process.pH      = TestPlan(k,2); % Logarithmic pH scale
    reacstruc.process.Co      = TestPlan(k,3); % g/L
    reacstruc.process.lambda0 = TestPlan(k,4); % mol /mol
    reacstruc.process.tdose   = TestPlan(k,5); % g/L

    % Run sim with new parameters
    reacstruc = reacsim(reacstruc);

    % Extract yield and triacylated components
    yD_test(k,1) = reacstruc.out.y(end,7);
    yH_test(k,1) = reacstruc.out.y(end,10);
end

%% Compute predicted response and 95% confidence interval for MLR model on test set
[MLR_yfittest_Yield, MLR_yinttest_Yield] =  predict(MLR_Yield, XTest, 'alpha', alpha);
[MLR_yfittest_Impurity, MLR_yinttest_Impurity] =  predict(MLR_Impurity, XTest, 'alpha', alpha);

%% Compute predicted response and 95% confidence intervals for GPR model on test set
[GPR_yfittest_Yield, ~, GPR_yinttest_Yield] = predict(gprMdl_Yield, XTest, 'alpha', alpha);
[GPR_yfittest_Impurity, ~, GPR_yinttest_Impurity] = predict(gprMdl_Tri, XTest, 'alpha', alpha);

%% Compute predicted response for ANN model on test set
ANN_yfittest_Yield = predict(ANN_Yield, XTest);  
ANN_yfittest_Impurity = predict(ANN_Impurity, XTest);  

%% Compute statistics on test data
% Compute (mechanistic)
stats_test_MM_Yield = rs_stats(YTest_Yield, yD_test);
stats_test_MLR_Yield = rs_stats(YTest_Yield, MLR_yfittest_Yield);
stats_test_GPR_Yield = rs_stats(YTest_Yield, GPR_yfittest_Yield);
stats_test_ANN_Yield = rs_stats(YTest_Yield, ANN_yfittest_Yield);

stats_test_MM_Impurity = rs_stats(YTest_Impurity, yH_test);
stats_test_MLR_Impurity = rs_stats(YTest_Impurity, MLR_yfittest_Impurity);
stats_test_GPR_Impurity = rs_stats(YTest_Impurity, GPR_yfittest_Impurity);
stats_test_ANN_Impurity = rs_stats(YTest_Impurity, ANN_yfittest_Impurity);

R2_test_Yield = [stats_test_MM_Yield.R2; stats_test_MLR_Yield.R2; ...
        stats_test_GPR_Yield.R2; stats_test_ANN_Yield.R2];
SSE_test_Yield = [stats_test_MM_Yield.SSE; stats_test_MLR_Yield.SSE; ...
        stats_test_GPR_Yield.SSE; stats_test_ANN_Yield.SSE];
MSE_test_Yield = [stats_test_MM_Yield.MSE; stats_test_MLR_Yield.MSE; ...
        stats_test_GPR_Yield.MSE; stats_test_ANN_Yield.MSE];
RMSE_test_Yield = [stats_test_MM_Yield.RMSE; stats_test_MLR_Yield.RMSE; ...
        stats_test_GPR_Yield.RMSE; stats_test_ANN_Yield.RMSE];
AAD_test_Yield = [stats_test_MM_Yield.AAD; stats_test_MLR_Yield.AAD; ...
        stats_test_GPR_Yield.AAD; stats_test_ANN_Yield.AAD];
RE_test_Yield(:,1) = stats_test_MM_Yield.RE; 
RE_test_Yield(:,2) = stats_test_MLR_Yield.RE;
RE_test_Yield(:,3) = stats_test_GPR_Yield.RE; 
RE_test_Yield(:,4) = stats_test_ANN_Yield.RE;

R2_test_Impurity = [stats_test_MM_Impurity.R2; stats_test_MLR_Impurity.R2; ...
        stats_test_GPR_Impurity.R2; stats_test_ANN_Impurity.R2];
SSE_test_Impurity = [stats_test_MM_Impurity.SSE; stats_test_MLR_Impurity.SSE; ...
        stats_test_GPR_Impurity.SSE; stats_test_ANN_Impurity.SSE];
MSE_test_Impurity = [stats_test_MM_Impurity.MSE; stats_test_MLR_Impurity.MSE; ...
        stats_test_GPR_Impurity.MSE; stats_test_ANN_Impurity.MSE];
RMSE_test_Impurity = [stats_test_MM_Impurity.RMSE; stats_test_MLR_Impurity.RMSE; ...
        stats_test_GPR_Impurity.RMSE; stats_test_ANN_Impurity.RMSE];
AAD_test_Impurity = [stats_test_MM_Impurity.AAD; stats_test_MLR_Impurity.AAD; ...
        stats_test_GPR_Impurity.AAD; stats_test_ANN_Impurity.AAD];

[R2_test_Yield SSE_test_Yield RMSE_test_Yield AAD_test_Yield]
[R2_test_Impurity SSE_test_Impurity RMSE_test_Impurity AAD_test_Impurity]

%% Compute 5-fold CV errors
rng(0, 'twister')
% CV_MLR_Yield    = kfoldLoss(crossval(MLR_Yield, 'KFold', 10));
% CV_MLR_Impurity = kfoldLoss(crossval(gprMdl_Yield, 'KFold', 10));
CV_GPR_Yield    = kfoldLoss(crossval(gprMdl_Yield, 'KFold', 10));
CV_GPR_Impurity = kfoldLoss(crossval(gprMdl_Tri, 'KFold', 10));
CV_ANN_Yield    = kfoldLoss(crossval(ANN_Yield, 'KFold', 10));
CV_ANN_Impurity = kfoldLoss(crossval(ANN_Impurity, 'KFold', 10));

disp(sqrt([CV_GPR_Yield CV_GPR_Impurity; ...
      CV_ANN_Yield CV_ANN_Impurity]))

%% Figures of fit to test data (Yield)
[~, IY]  = sort(YTest_Yield);
[~, IYr] = sort(YTest_Yield, 'descend');

figure
subplot(2, 2, 1); hold all % Mechanistic 
% plot(YTest_Yield, YTest_Yield, '-')
plot(-1:0.01:2, -1:0.01:2, 'k-')
plot(YTest_Yield, yD_test, 'x')
xlim([min(min(min(YTest_Yield,yD_test)),0) max(max(max(YTest_Yield,yD_test)),1)])
ylim([min(min(min(YTest_Yield,yD_test)),0) max(max(max(YTest_Yield,yD_test)),1)])
xlabel('Measured yield ($\%$)', 'FontSize', 10)
ylabel('Yield ($\%$)', 'FontSize', 10)
xlim([-0.2 1.35])
ylim([-0.2 1.35])
title('Mechanistic', 'FontSize', 10)
set(gca, 'FontSize', 11)
legend("","Predicted", ...
        'location', 'northwest', 'FontSize', 8, 'box', 'on')
text(0.6,-0.00, strcat("$R^{2}$: ", num2str(round(R2_test_Yield(1),4)*100), "\%"), 'FontSize', 10)
text(0.6,-0.10, strcat("RMSE: ", num2str(round(RMSE_test_Yield(1),4))), 'FontSize', 10)

subplot(2, 2, 2); hold all  % DOE (MLR)
% plot(YTest_Yield, YTest_Yield, '-')
plot(-1:0.01:2, -1:0.01:2, 'k-')
plot(YTest_Yield, MLR_yfittest_Yield, 'x')
patch([YTest_Yield(IY)' YTest_Yield(IYr)'], [MLR_yinttest_Yield(IY,1)' MLR_yinttest_Yield(IYr,2)'], ...
      'blue', 'FaceAlpha', .125) 
xlabel('Measured yield ($\%$)', 'FontSize', 10)
ylabel('Yield ($\%$)', 'FontSize', 10)
% xlim([min(min(min([YTest_Yield MLR_yfittest_Yield MLR_yinttest_Yield])),0) ...
%       max(max(max([YTest_Yield MLR_yfittest_Yield MLR_yinttest_Yield])),1)])
% ylim([min(min(min([YTest_Yield MLR_yfittest_Yield MLR_yinttest_Yield])),0) ...
%       max(max(max([YTest_Yield MLR_yfittest_Yield MLR_yinttest_Yield])),1)])
xlim([-0.2 1.35])
ylim([-0.2 1.35])
title('MLR', 'FontSize', 10)
set(gca, 'FontSize', 11)
legend("", "Predicted", "95 \% PI", ...
       'FontSize', 8, 'location', 'northwest', 'box', 'on')
text(0.6,-0.00, strcat("$R^{2}$: ", num2str(round(R2_test_Yield(2),4)*100), "\%"), 'FontSize', 10)
text(0.6,-0.10, strcat("RMSE: ", num2str(round(RMSE_test_Yield(2),4))), 'FontSize', 10)

subplot(2, 2, 3); hold all % GPR 
% plot(YTest_Yield, YTest_Yield, '-')
plot(-1:0.01:2, -1:0.01:2, 'k-')
plot(YTest_Yield, GPR_yfittest_Yield, 'x')
patch([YTest_Yield(IY)' YTest_Yield(IYr)'], ...
    [GPR_yinttest_Yield(IY,1)' GPR_yinttest_Yield(IYr,2)'], ...
      'blue', 'FaceAlpha', .125) 
xlabel('Measured yield ($\%$)', 'FontSize', 10)
ylabel('Yield ($\%$)', 'FontSize', 10)
xlim([-0.2 1.35])
ylim([-0.2 1.35])
title('GPR', 'FontSize', 10)
set(gca, 'FontSize', 11)
legend("", "Predicted", "95 \% PI", ...
       'location', 'northwest', 'FontSize', 8, 'box', 'on')
text(0.6,-0.00, strcat("$R^{2}$: ", num2str(round(R2_test_Yield(3),4)*100), "\%"), 'FontSize', 10)
text(0.6,-0.10, strcat("RMSE: ", num2str(round(RMSE_test_Yield(3),4))), 'FontSize', 10)
   
subplot(2, 2, 4); hold all % ANN
% plot(YTest_Yield, YTest_Yield, '-')
plot(-1:0.01:2, -1:0.01:2, 'k-')
plot(YTest_Yield, ANN_yfittest_Yield, 'x')
xlabel('Measured yield ($\%$)', 'FontSize', 10)
ylabel('Yield ($\%$)', 'FontSize', 10)
xlim([-0.2 1.35])
ylim([-0.2 1.35])
title('ANN', 'FontSize', 10)
set(gca, 'FontSize', 11)
legend("", "Predicted", ...
       'location', 'northwest', 'FontSize', 8, 'box', 'on')
FigureTitle(5) = "FitLHS_TestFit_Yield";
text(0.6,-0.00, strcat("$R^{2}$: ", num2str(round(R2_test_Yield(4),4)*100), "\%"), 'FontSize', 10)
text(0.6,-0.10, strcat("RMSE: ", num2str(round(RMSE_test_Yield(4),4))), 'FontSize', 10)
sgtitle('LHS - Test set')

%% Figures of fit to test data (Impurity)
[~, IH]  = sort(YTest_Impurity);
[~, IHr] = sort(YTest_Impurity, 'descend');

figure
subplot(2, 2, 1); hold all % Mechanistic 
% plot(YTest_Impurity, YTest_Impurity, '-')
plot(-1:0.01:2, -1:0.01:2, 'k-')
plot(YTest_Impurity, yH_test, 'x')
% xlim([min(min(min(YTest_Impurity, yD_test)),0) min(max(max(YTest_Impurity,yH_test)),1)])
% ylim([min(min(min(YTest_Impurity, yD_test)),0) min(max(max(YTest_Impurity,yH_test)),1)])
xlabel('Measured impurity ($\%$)', 'FontSize', 10)
ylabel('Impurity ($\%$)', 'FontSize', 10)
xlim([-0.2 0.4])
ylim([-0.2 0.4])
title('Mechanistic', 'FontSize', 10)
% patch([yD_test(I)' yD_test(Ir)'],                     ...
%       [ModelConfidenceInterval95(I,1)'      ...
%        ModelConfidenceInterval95(Ir,2)'],   ...
%       'blue', 'FaceAlpha', .125)
set(gca, 'FontSize', 11)
legend("", "Predicted", ...
        'location', 'northwest', 'FontSize', 8, 'box', 'on')
text(0.15, -0.135, strcat("$R^{2}$: ", num2str(round(R2_test_Impurity(1),4)*100), "\%"), 'FontSize', 10)
text(0.15, -0.17, strcat("RMSE: ", num2str(round(RMSE_test_Impurity(1),4))), 'FontSize', 10)

subplot(2, 2, 2); hold all  % DOE (MLR)
% plot(YTest_Impurity, YTest_Impurity, '-')
plot(-1:0.01:2, -1:0.01:2, 'k-')
plot(YTest_Impurity, MLR_yfittest_Impurity, 'x')
patch([YTest_Impurity(IH)' YTest_Impurity(IHr)'], ...
    [MLR_yinttest_Impurity(IH,1)' MLR_yinttest_Impurity(IHr,2)'], ...
      'blue', 'FaceAlpha', .125) 
xlabel('Measured impurity ($\%$)', 'FontSize', 10)
ylabel('Impurity ($\%$)', 'FontSize', 10)
xlim([-0.2 0.4])
ylim([-0.2 0.4])
% xlim([min(min(min([YTest_Impurity MLR_yfittest_Impurity MLR_yinttest_Impurity])),0) ...
%       min(max(max([YTest_Impurity MLR_yfittest_Impurity MLR_yinttest_Impurity])),1)])
% ylim([min(min(min([YTest_Impurity MLR_yfittest_Impurity MLR_yinttest_Impurity])),0) ...
%       min(max(max([YTest_Impurity MLR_yfittest_Impurity MLR_yinttest_Impurity])),1)])
title('MLR', 'FontSize', 10)
set(gca, 'FontSize', 11)
legend("", "Predicted", "95 \% PI", ...
       'FontSize', 8, 'location', 'northwest', 'box', 'on')
text(0.15, -0.135, strcat("$R^{2}$: ", num2str(round(R2_test_Impurity(2),4)*100), "\%"), 'FontSize', 10)
text(0.15, -0.17, strcat("RMSE: ", num2str(round(RMSE_test_Impurity(2),4))), 'FontSize', 10)

subplot(2, 2, 3); hold all % GPR 
% plot(YTest_Impurity, YTest_Impurity, '-')
plot(-1:0.01:2, -1:0.01:2, 'k-')
plot(YTest_Impurity, GPR_yfittest_Impurity, 'x')
patch([YTest_Impurity(IH)' YTest_Impurity(IHr)'], ...
    [GPR_yinttest_Impurity(IH,1)' GPR_yinttest_Impurity(IHr,2)'], ...
      'blue', 'FaceAlpha', .125) 
xlabel('Measured impurity ($\%$)', 'FontSize', 10)
ylabel('Impurity ($\%$)', 'FontSize', 10)
xlim([-0.2 0.4])
ylim([-0.2 0.4])
title('GPR', 'FontSize', 10)
set(gca, 'FontSize', 11)
legend("", "Predicted", "95 \% PI", ...
       'location', 'northwest', 'FontSize', 8, 'box', 'on')
text(0.15, -0.135, strcat("$R^{2}$: ", num2str(round(R2_test_Impurity(3),4)*100), "\%"), 'FontSize', 10)
text(0.15, -0.17, strcat("RMSE: ", num2str(round(RMSE_test_Impurity(3),4))), 'FontSize', 10)
   
subplot(2, 2, 4); hold all % ANN
% plot(YTest_Impurity, YTest_Impurity, '-')
plot(-1:0.01:2, -1:0.01:2, 'k-')
plot(YTest_Impurity, ANN_yfittest_Impurity, 'x')
xlabel('Measured impurity ($\%$)', 'FontSize', 10)
ylabel('Impurity ($\%$)', 'FontSize', 10)
xlim([-0.2 0.4])
ylim([-0.2 0.4])
title('ANN', 'FontSize', 10)
set(gca, 'FontSize', 11)
legend("", "Predicted", ...
       'location', 'northwest', 'FontSize', 8, 'box', 'on')
FigureTitle(6) = "FitLHS_TestFit_Tri";
text(0.15, -0.135, strcat("$R^{2}$: ", num2str(round(R2_test_Impurity(4),4)*100), "\%"), 'FontSize', 10)
text(0.15, -0.17, strcat("RMSE: ", num2str(round(RMSE_test_Impurity(4),4))), 'FontSize', 10)
sgtitle('LHS - Test set')

%% Global sensitivity analysis (GSA)
%% Sobol's method
f_MLR_Yield     = @(x) predict(MLR_Yield, x);
f_MLR_Impurity  = @(x) predict(MLR_Impurity, x);
f_GPR_Yield     = @(x) predict(gprMdl_Yield, x);
f_GPR_Impurity  = @(x) predict(gprMdl_Tri, x);
f_ANN_Yield     = @(x) predict(ANN_Yield, x);
f_ANN_Impurity  = @(x) predict(ANN_Impurity, x);
f_MM_Yield      = @(x) reacsim_yDwrapper(x, reacstruc);
f_MM_Impurity   = @(x) reacsim_yHwrapper(x, reacstruc);

N = 1e5;
pars = {'T', 'pH', 'Co', 'lambda0'};
lbs = [5 10 5 1];
ubs = [40 12 50 3];

pars_MM = {'T', 'pH', 'Co', 'lambda0', 'tdose'};
lbs_MM = [lbs 10];
ubs_MM = [ubs 60];
N_MM = 1e4;
InputSpace = {'ParNames', pars, 'LowerBounds', lbs, 'UpperBounds', ubs};
InputSpace_MM = {'ParNames', pars_MM, 'LowerBounds', lbs_MM, ...
                                      'UpperBounds', ubs_MM};

% Compute Sobol indices
[Si_MLR_Yield, STi_MLR_Yield]       = easyGSA(f_MLR_Yield, N, InputSpace{:},'UseParallel',true);
[Si_MLR_Impurity, STi_MLR_Impurity] = easyGSA(f_MLR_Impurity, N, InputSpace{:},'UseParallel',true);
[Si_GPR_Yield, STi_GPR_Yield]       = easyGSA(f_GPR_Yield, N, InputSpace{:},'UseParallel',true);
[Si_GPR_Impurity, STi_GPR_Impurity] = easyGSA(f_GPR_Impurity, N, InputSpace{:},'UseParallel',true);
[Si_ANN_Yield, STi_ANN_Yield]       = easyGSA(f_ANN_Yield, N, InputSpace{:},'UseParallel',true);
[Si_ANN_Impurity, STi_ANN_Impurity] = easyGSA(f_ANN_Impurity, N, InputSpace{:},'UseParallel',true);
% [Si_MM_Yield, STi_MM_Yield]         = easyGSA(f_MM_Yield, N_MM, InputSpace_MM{:},'UseParallel',true);
% [Si_MM_Impurity, STi_MM_Impurity]   = easyGSA(f_MM_Impurity, N_MM, InputSpace_MM{:},'UseParallel',true);

disp([Si_MLR_Yield, STi_MLR_Yield])
disp([Si_MLR_Impurity STi_MLR_Impurity])
disp([Si_GPR_Yield, STi_GPR_Yield])
disp([Si_GPR_Impurity, STi_GPR_Impurity])
disp([Si_ANN_Yield, STi_ANN_Yield])
disp([Si_ANN_Impurity, STi_ANN_Impurity])
% disp([Si_MM_Yield, STi_MM_Yield])
% disp([Si_MM_Impurity, STi_MM_Impurity])

% Notes:
% No noise:
%     0.0037    0.0519
%     0.0503    0.1140
%     0.0129    0.0714
%     0.8189    0.8774
% 
%    -0.0006    0.1262
%     0.1893    0.4874
%     0.0318    0.1965
%     0.2916    0.6730
% 
%     0.0053    0.0007
%     0.0427    0.1390
%     0.0078    0.0275
%     0.8512    0.9359
% 
%     0.0056    0.0000
%     0.2658    0.5865
%     0.0463    0.1157
%     0.3398    0.6603
% 
%     0.0105    0.0272
%     0.0287    0.0808
%     0.0052    0.0470
%     0.8919    0.9205
% 
%     0.0134    0.0661
%     0.2307    0.5205
%     0.0560    0.1597
%     0.3265    0.6493
% 
%     0.0108    0.0137
%     0.0427    0.1137
%     0.0074    0.0711
%     0.8446    0.8966
%    -0.0007    0.0000
% 
%     0.0021    0.0292
%     0.2776    0.5306
%     0.0314    0.1228
%     0.3738    0.7195
%    -0.0092    0.0000

% With noise
%     0.0184    0.0472
%     0.0688    0.1493
%     0.0050    0.0427
%     0.8017    0.8682
% 
%     0.0004    0.1346
%     0.1802    0.4723
%     0.0425    0.2160
%     0.2884    0.6613
% 
%     0.0135    0.0103
%     0.0686    0.1809
%     0.0045    0.0000
%     0.8088    0.9212
% 
%     0.0050    0.0000
%     0.2515    0.5644
%     0.0561    0.1543
%     0.3353    0.6611
% 
%     0.0633    0.0661
%     0.0049    0.0171
%     0.0048    0.0052
%     0.9147    0.9271
% 
%     0.0221    0.1389
%     0.2039    0.4776
%     0.0913    0.2425
%     0.2868    0.5509
% 
%     0.0341    0.0458
%     0.0388    0.0972
%     0.0142    0.0738
%     0.8264    0.8754
%     0.0014    0.0000
% 
%    -0.0008    0.0179
%     0.2989    0.5466
%     0.0325    0.1230
%     0.3671    0.6958
%    -0.0091    0.0000

%% Propagation of uncertainty with Monte Carlo simulation
% if SimulationIndex == 2
    % Sample parameter values from multivariate normal distribution induced by
    % the linear error propagation method
    Theta_MVN = abs(mvnrnd(ParameterMinimum, CovarianceParametersD, 1000));
    
    % Setpoint
    reacstruc.process.T       = 20; % Celcius scale
    reacstruc.process.pH      = 11; % Logarithmic pH scale
    reacstruc.process.Co      = 20; % g/L
    reacstruc.process.lambda0 = 2;  % mol /mol
    reacstruc.process.tdose   = 30; % min
    
    %
    yD_theta = zeros(size(Theta_MVN,1), 1);
    yH_theta = zeros(size(Theta_MVN,1), 1);
    for k = 1:size(Theta_MVN,1)
        % Update parameters
        reacstruc.model.kAref   = Theta_MVN(k,1);
        reacstruc.model.sref(3) = Theta_MVN(k,2);
        reacstruc.model.sdref   = Theta_MVN(k,3);
        reacstruc.model.pkA1    = Theta_MVN(k,4);   
        reacstruc.model.pkA3    = Theta_MVN(k,5);   
        reacstruc.model.EA1     = Theta_MVN(k,6);    
        reacstruc.model.EA3     = Theta_MVN(k,7);   
        reacstruc.model.EAd     = Theta_MVN(k,8);
    
        % Run sim with new parameters
        reacstruc = reacsim(reacstruc);
    
        % Extract yield and triacylated components
        yD_theta(k,1) = reacstruc.out.y(end,7);
        yH_theta(k,1) = reacstruc.out.y(end,10);
    end
    
    figure
    [~,ax] = plotmatrix(Theta_MVN);
    ax(1,1).YLabel.String = "$k_{A,ref}$";
    ax(2,1).YLabel.String = "$s_{ref,3}$";
    ax(3,1).YLabel.String = "$s_{d,ref}$";
    ax(4,1).YLabel.String = "pk$_{a,1}$";
    ax(5,1).YLabel.String = "pk$_{a,3}$";
    ax(6,1).YLabel.String = "E$_{A,1}$";
    ax(7,1).YLabel.String = "E$_{A,3}$";
    ax(8,1).YLabel.String = "E$_{A,d}$";
    
    ax(8,1).XLabel.String = "$k_{A,ref}$";
    ax(8,2).XLabel.String = "$s_{ref,3}$";
    ax(8,3).XLabel.String = "$s_{d,ref}$";
    ax(8,4).XLabel.String = "pk$_{a,1}$";
    ax(8,5).XLabel.String = "pk$_{a,3}$";
    ax(8,6).XLabel.String = "E$_{A,1}$";
    ax(8,7).XLabel.String = "E$_{A,3}$";
    ax(8,8).XLabel.String = "E$_{A,d}$";
    for i = 1:8
        ax(8,i).XAxis.FontSize = 7;
        ax(i,1).YAxis.FontSize = 7;
    end
    FigureTitle(7) = "FitLHS_UA_MM_PlotMatrix";
    
    figure
    histogram(yD_theta, 'normalization', 'probability')
    ylabel('Frequency (probability)')
    xlabel('Yield, $y_{D}$')
    xline(mean(yD_theta), 'r-', 'LineWidth', 2)
    xline(0.8150, 'k-', 'LineWidth', 2)
    FigureTitle(8) = "FitLHS_UA_MM_Yield";
    
    figure
    histogram(yH_theta, 'normalization', 'probability')
    ylabel('Frequency (probability)')
    xlabel('Impurity, $y_{H}$')
    xline(mean(yH_theta), 'r-', 'LineWidth', 2)
    xline(0.0114, 'k-', 'LineWidth', 2)
    FigureTitle(9) = "FitLHS_UA_MM_Tri";
% end

%% Simulation-based optimization
% Set vector of prices
price = [0; reacstruc.optim.Price.SCrel; 2];

% 
lows      = [5 10 5 1];     % Low values
highs     = [40 12 50 3];   % High values
N_Samples = 1000000;         % Number of sample points

% Generate test data
rng(20, 'twister')
OptPlan = zeros(N_Samples, 4);
for i = 1:4
    OptPlan(:,i) = icdf('Uniform', lhsdesign(N_Samples,1), lows(i), highs(i));
end

%% MLR
% Compute predictions over entire input space with N samples
YD_MLR = predict(MLR_Yield, OptPlan);
YH_MLR = predict(MLR_Impurity, OptPlan);

% Compute optimum without and with upper bound of 0.01 impurity
k = cell(2,1);

% Find indices for feasible initial points
k{1} = find(0 < YH_MLR & YH_MLR <= 1 & 0 < YD_MLR & YD_MLR <= 1);

% Find indices for feasible initial points below for which the impurity 
% component is below the constraint
k{2} = find(0 < YH_MLR & YH_MLR <= 0.005 & 0 < YD_MLR & YD_MLR <= 1);
    
% Set structure array object for solving optimization problem with multiple
% linear regression model
optim_MLR.LB            = lows;
optim_MLR.UB            = highs;
optim_MLR.ModelYield    = MLR_Yield;
optim_MLR.ModelImpurity = MLR_Impurity;
optim_MLR.price         = 0;
optim_MLR.z0            = [20 11 20 2];
optim_MLR.glb           = [0 0];
optim_MLR.gub           = [1 1];

% Pre-allocate
MLR_sol         = zeros(6, 4);  %
MLR_sol_fmincon = zeros(6, 4);
MLR_opt         = zeros(6, 4);  %
MLR_opt_fmincon = zeros(6,4);
IdxOptGuess_MLR = zeros(6, 1);   % 
costGuess_MLR   = zeros(6, 1);
for i = 1:length(price)
    % Compute best initial guess with default constraints
    [costGuess_MLR(2*i-1), IdxOptGuess_MLR(2*i-1)] = ...
        min((1 + price(i)*OptPlan(k{1},4)) ./ YD_MLR(k{1}));
    
    % Compute best initial guess with impurity below the constraint of 0.01
    [cost(2*i), IdxOptGuess_MLR(2*i)] = ...
        min((1 + price(i)*OptPlan(k{2},4)) ./ YD_MLR(k{2}));

    % Compute initial point
%     [cost, IdxOpt] = min((1 + 0.5370*OptPlan(k,4)) ./ YD_MLR(k));
    
    % Compute minimizer with default constraints
    optim_MLR.z0     = OptPlan(k{1}(IdxOptGuess_MLR(2*i-1)),:);  % Set initial guess
    optim_MLR.glb    = [0 0];                           % Lower constraints
    optim_MLR.gub    = [1 1];                           % Upper constraints
    optim_MLR.price  = price(i);
    MLR_sol(2*i-1,:) = runoptim_lm(optim_MLR, 'fminsearchcon');
    MLR_sol_fmincon(2*i-1,:) = runoptim_lm(optim_MLR, 'fmincon');
    MLR_opt(2*i-1,:) = MLR_sol(2*i-1,:) .* (highs - lows) + lows;
    MLR_opt_fmincon(2*i-1,:) = MLR_sol_fmincon(2*i-1,:) .* (highs - lows) + lows;

    % Compute minimizer with impurity below the constraint of 0.01
    optim_MLR.z0   = OptPlan(k{2}(IdxOptGuess_MLR(2*i)),:);  % Set initial guess
    optim_MLR.glb  = [0 0];                           % Lower constraints
    optim_MLR.gub  = [1 0.01];                        % Upper constraints
    MLR_sol(2*i,:) = runoptim_lm(optim_MLR, 'fminsearchcon');
    MLR_sol_fmincon(2*i,:) = runoptim_lm(optim_MLR, 'fmincon');
    MLR_opt(2*i,:) = MLR_sol(2*i,:) .* (highs - lows) + lows;
    MLR_opt_fmincon(2*i,:) = MLR_sol_fmincon(2*i,:) .* (highs - lows) + lows;
end

%% GPR
% Compute predictions over entire input space with N samples
YD_GPR = predict(gprMdl_Yield, OptPlan);
YH_GPR = predict(gprMdl_Tri, OptPlan);

% Compute optimum without and with upper bound of 0.01 impurity
k = cell(2,1);

% Find indices for feasible initial points
k{1} = find(0 < YH_GPR & YH_GPR <= 1 & 0 < YD_GPR & YD_GPR <= 1);

% Find indices for feasible initial points below for which the impurity 
% component is below the constraint
k{2} = find(0 < YH_GPR & YH_GPR <= 0.005 & 0 < YD_GPR & YD_GPR <= 1);
    
% Set structure array object for solving optimization problem with multiple
% linear regression model
optim_GPR.LB            = lows;
optim_GPR.UB            = highs;
optim_GPR.ModelYield    = gprMdl_Yield;
optim_GPR.ModelImpurity = gprMdl_Tri;
optim_GPR.price         = 0;
optim_GPR.z0            = [20 11 20 2];
optim_GPR.glb           = [0 0];
optim_GPR.gub           = [1 1];

% Pre-allocate
GPR_sol          = zeros(6, 4);  %
GPR_sol_fmincon  = zeros(6, 4);
GPR_opt          = zeros(6, 4);  %
GPR_opt_fmincon  = zeros(6, 4);
IdxOptGuess_GPR  = zeros(6, 1);   % 
costGuess_GPR    = zeros(6, 1);
for i = 1:length(price)
    optim_GPR.price = price(i);

    % Compute best initial guess with default constraints
    [costGuess_GPR(2*i-1), IdxOptGuess_GPR(2*i-1)] = ...
        min((1 + price(i)*OptPlan(k{1},4)) ./ YD_GPR(k{1}));
    
    % Compute best initial guess with impurity below the constraint of 0.01
    [costGuess_GPR(2*i), IdxOptGuess_GPR(2*i)] = ...
        min((1 + price(i)*OptPlan(k{2},4)) ./ YD_GPR(k{2}));

    % Compute initial point
%     [cost, IdxOpt] = min((1 + 0.5370*OptPlan(k,4)) ./ YD_MLR(k));
    
    % Compute minimizer with default constraints
    optim_GPR.z0     = OptPlan(k{1}(IdxOptGuess_GPR(2*i-1)),:);  % Set initial guess
%     optim_GPR.z0     = TrueOpt(2*i-1,:);
    optim_GPR.glb    = [0 0];                           % Lower constraints
    optim_GPR.gub    = [1 1];                           % Upper constraints
    optim_GPR.price  = price(i);
    GPR_sol(2*i-1,:) = runoptim_lm(optim_GPR, 'fminsearchcon');
    GPR_sol_fmincon(2*i-1,:) = runoptim_lm(optim_GPR, 'fmincon');
    GPR_opt(2*i-1,:) = GPR_sol(2*i-1,:) .* (highs - lows) + lows;
    GPR_opt_fmincon(2*i-1,:) = GPR_sol_fmincon(2*i-1,:) .* (highs - lows) + lows;

    % Compute minimizer with impurity below the constraint of 0.01
    optim_GPR.z0   = OptPlan(k{2}(IdxOptGuess_GPR(2*i)),:);  % Set initial guess
%     optim_GPR.z0   = TrueOpt(2*i,:);
    optim_GPR.glb  = [0 0];                           % Lower constraints
    optim_GPR.gub  = [1 0.01];                        % Upper constraints
    GPR_sol(2*i,:) = runoptim_lm(optim_GPR, 'fminsearchcon');
    GPR_sol_fmincon(2*i,:) = runoptim_lm(optim_GPR, 'fmincon');
    GPR_opt(2*i,:) = GPR_sol(2*i,:) .* (highs - lows) + lows;
    GPR_opt_fmincon(2*i,:) = GPR_sol_fmincon(2*i,:) .* (highs - lows) + lows;
end

%% ANN
% Compute predictions over entire input space with N samples
YD_ANN = predict(ANN_Yield, OptPlan);
YH_ANN = predict(ANN_Impurity, OptPlan);

% Compute optimum without and with upper bound of 0.01 impurity
k = cell(2,1);

% Find indices for feasible initial points
k{1} = find(0 < YH_ANN & YH_ANN <= 1 & 0 < YD_ANN & YD_ANN <= 1);

% Find indices for feasible initial points below for which the impurity 
% component is below the constraint
k{2} = find(0 < YH_ANN & YH_ANN <= 0.005 & 0 < YD_ANN & YD_ANN <= 1);
    
% Set structure array object for solving optimization problem with multiple
% linear regression model
optim_ANN.LB            = lows;
optim_ANN.UB            = highs;
optim_ANN.ModelYield    = ANN_Yield;
optim_ANN.ModelImpurity = ANN_Impurity;
optim_ANN.price         = 0;
optim_ANN.z0            = [20 11 20 2];
optim_ANN.glb           = [0 0];
optim_ANN.gub           = [1 1];

% Pre-allocate
ANN_sol          = zeros(6, 4);  %
ANN_sol_fmincon  = zeros(6, 4);
ANN_opt          = zeros(6, 4);  %
ANN_opt_fmincon  = zeros(6, 4);
IdxOptGuess_ANN  = zeros(6, 1);   % 
costGuess_ANN    = zeros(6, 1);
for i = 1:length(price)
    optim_ANN.price = price(i);

    % Compute best initial guess with default constraints
    [costGuess_ANN(2*i-1), IdxOptGuess_ANN(2*i-1)] = ...
        min((1 + price(i)*OptPlan(k{1},4)) ./ YD_ANN(k{1}));
    
    % Compute best initial guess with impurity below the constraint of 0.01
    [costGuess_ANN(2*i), IdxOptGuess_ANN(2*i)] = ...
        min((1 + price(i)*OptPlan(k{2},4)) ./ YD_ANN(k{2}));

    % Compute initial point
%     [cost, IdxOpt] = min((1 + 0.5370*OptPlan(k,4)) ./ YD_MLR(k));
    
    % Compute minimizer with default constraints
    optim_ANN.z0     = OptPlan(k{1}(IdxOptGuess_ANN(2*i-1)),:);  % Set initial guess
    optim_ANN.glb    = [0 0];                           % Lower constraints
    optim_ANN.gub    = [1 1];                           % Upper constraints
    optim_ANN.price  = price(i);
    ANN_sol(2*i-1,:) = runoptim_lm(optim_ANN, 'fminsearchcon');
    ANN_sol_fmincon(2*i-1,:) = runoptim_lm(optim_ANN, 'fmincon');
    ANN_opt(2*i-1,:) = ANN_sol(2*i-1,:) .* (highs - lows) + lows;
    ANN_opt_fmincon(2*i-1,:) = ANN_sol_fmincon(2*i-1,:) .* (highs - lows) + lows;

    % Compute minimizer with impurity below the constraint of 0.01
    optim_ANN.z0   = OptPlan(k{2}(IdxOptGuess_ANN(2*i)),:);  % Set initial guess
    optim_ANN.glb  = [0 0];                           % Lower constraints
    optim_ANN.gub  = [1 0.01];                        % Upper constraints
    ANN_sol(2*i,:) = runoptim_lm(optim_ANN, 'fminsearchcon');
    ANN_sol_fmincon(2*i,:) = runoptim_lm(optim_ANN, 'fmincon');
    ANN_opt(2*i,:) = ANN_sol(2*i,:) .* (highs - lows) + lows;
    ANN_opt_fmincon(2*i,:) = ANN_sol_fmincon(2*i,:) .* (highs - lows) + lows;
end

MLR_opt
GPR_opt
ANN_opt
MLR_opt_fmincon
GPR_opt_fmincon
ANN_opt_fmincon

%% Run estimated optima fmincon
MLR_opt_fmincon(:,5) = 30;
GPR_opt_fmincon(:,5) = 30;
ANN_opt_fmincon(:,5) = 30;

MLR_opt_data_fmincon = instantlab(MLR_opt_fmincon, test_sigma_input, 0);
GPR_opt_data_fmincon = instantlab(GPR_opt_fmincon, test_sigma_input, 0);
ANN_opt_data_fmincon = instantlab(ANN_opt_fmincon, test_sigma_input, 0);

MLR_opt_YTest_Yield_fmincon    = zeros(6, 1);
MLR_opt_YTest_Impurity_fmincon = zeros(6, 1);
GPR_opt_YTest_Yield_fmincon    = zeros(6, 1);
GPR_opt_YTest_Impurity_fmincon = zeros(6, 1);
ANN_opt_YTest_Yield_fmincon    = zeros(6, 1);
ANN_opt_YTest_Impurity_fmincon = zeros(6, 1);
for i = 1:size(MLR_opt,1)
    MLR_opt_YTest_Yield_fmincon(i)    = MLR_opt_data_fmincon.out{i}(6);
    MLR_opt_YTest_Impurity_fmincon(i) = MLR_opt_data_fmincon.out{i}(5);

    GPR_opt_YTest_Yield_fmincon(i)    = GPR_opt_data_fmincon.out{i}(6);
    GPR_opt_YTest_Impurity_fmincon(i) = GPR_opt_data_fmincon.out{i}(5);

    ANN_opt_YTest_Yield_fmincon(i)    = ANN_opt_data_fmincon.out{i}(6);
    ANN_opt_YTest_Impurity_fmincon(i) = ANN_opt_data_fmincon.out{i}(5);
end

MLR_opt_cost_fmincon          = zeros(6, 1);
MLR_opt_YTest_ActCost_fmincon = zeros(6, 1);
GPR_opt_cost_fmincon          = zeros(6, 1);
GPR_opt_YTest_ActCost_fmincon = zeros(6, 1);
ANN_opt_cost_fmincon          = zeros(6, 1);
ANN_opt_YTest_ActCost_fmincon = zeros(6, 1);
for i = 1:length(price)
    MLR_opt_cost_fmincon(2*i-1)          = (1 + price(i)*MLR_opt_fmincon(2*i-1,4)) / predict(MLR_Yield, MLR_opt_fmincon(2*i-1,1:4));
    MLR_opt_cost_fmincon(2*i)            = (1 + price(i)*MLR_opt_fmincon(2*i,4)) / predict(MLR_Yield, MLR_opt_fmincon(2*i,1:4));
    MLR_opt_YTest_ActCost_fmincon(2*i-1) = (1 + price(i)*MLR_opt_fmincon(2*i-1,4)) / MLR_opt_YTest_Yield_fmincon(2*i-1);
    MLR_opt_YTest_ActCost_fmincon(2*i)   = (1 + price(i)*MLR_opt_fmincon(2*i,4))   / MLR_opt_YTest_Yield_fmincon(2*i);

    GPR_opt_cost_fmincon(2*i-1)          = (1 + price(i)*GPR_opt_fmincon(2*i-1,4)) / predict(gprMdl_Yield, GPR_opt_fmincon(2*i-1,1:4));
    GPR_opt_cost_fmincon(2*i)            = (1 + price(i)*GPR_opt_fmincon(2*i,4)) / predict(gprMdl_Yield, GPR_opt_fmincon(2*i,1:4));
    GPR_opt_YTest_ActCost_fmincon(2*i-1) = (1 + price(i)*GPR_opt_fmincon(2*i-1,4)) / GPR_opt_YTest_Yield_fmincon(2*i-1);
    GPR_opt_YTest_ActCost_fmincon(2*i)   = (1 + price(i)*GPR_opt_fmincon(2*i,4))   / GPR_opt_YTest_Yield_fmincon(2*i);

    ANN_opt_cost_fmincon(2*i-1)          = (1 + price(i)*ANN_opt_fmincon(2*i-1,4)) / predict(ANN_Yield, ANN_opt_fmincon(2*i-1,1:4));
    ANN_opt_cost_fmincon(2*i)            = (1 + price(i)*ANN_opt_fmincon(2*i,4)) / predict(ANN_Yield, ANN_opt_fmincon(2*i,1:4));
    ANN_opt_YTest_ActCost_fmincon(2*i-1) = (1 + price(i)*ANN_opt_fmincon(2*i-1,4)) / ANN_opt_YTest_Yield_fmincon(2*i-1);
    ANN_opt_YTest_ActCost_fmincon(2*i)   = (1 + price(i)*ANN_opt_fmincon(2*i,4))   / ANN_opt_YTest_Yield_fmincon(2*i);
end

disp('True yield with estimated opt vs predicted opt yield (MLR) fmincon')
disp([MLR_opt_YTest_Yield_fmincon predict(MLR_Yield, MLR_opt_fmincon(:,1:4))])

disp('True impurity with estimated opt vs predicted opt impurity (MLR) fmincon')
disp([MLR_opt_YTest_Impurity_fmincon predict(MLR_Impurity, MLR_opt_fmincon(:,1:4))])

disp('True yield with estimated opt vs predicted opt yield (GPR) fmincon')
disp([GPR_opt_YTest_Yield_fmincon predict(gprMdl_Yield, GPR_opt_fmincon(:,1:4))])

disp('True impurity with estimated opt vs predicted opt impurity (GPR) fmincon')
disp([GPR_opt_YTest_Impurity_fmincon predict(gprMdl_Tri, GPR_opt_fmincon(:,1:4))])

disp('True yield with estimated opt vs predicted opt yield (ANN) fmincon')
disp([ANN_opt_YTest_Yield_fmincon predict(ANN_Yield, ANN_opt_fmincon(:,1:4))])

disp('True impurity with estimated opt vs predicted opt impurity (ANN) fmincon')
disp([ANN_opt_YTest_Impurity_fmincon predict(ANN_Impurity, ANN_opt_fmincon(:,1:4))])

disp('True cost with estimated opt vs predicted opt cost (MLR) fmincon')
disp([MLR_opt_YTest_ActCost_fmincon MLR_opt_cost_fmincon])

disp('True cost with estimated opt vs predicted opt cost (GPR) fmincon')
disp([GPR_opt_YTest_ActCost_fmincon GPR_opt_cost_fmincon])

disp('True cost with estimated opt vs predicted opt cost (ANN) fmincon')
disp([ANN_opt_YTest_ActCost_fmincon ANN_opt_cost_fmincon])

%% Run estimated optima fminsearchcon
MLR_opt(:,5) = 30;
GPR_opt(:,5) = 30;
ANN_opt(:,5) = 30;

MLR_opt_data = instantlab(MLR_opt, test_sigma_input, 0);
GPR_opt_data = instantlab(GPR_opt, test_sigma_input, 0);
ANN_opt_data = instantlab(ANN_opt, test_sigma_input, 0);

MLR_opt_YTest_Yield    = zeros(6, 1);
MLR_opt_YTest_Impurity = zeros(6, 1);
GPR_opt_YTest_Yield    = zeros(6, 1);
GPR_opt_YTest_Impurity = zeros(6, 1);
ANN_opt_YTest_Yield    = zeros(6, 1);
ANN_opt_YTest_Impurity = zeros(6, 1);
for i = 1:size(MLR_opt,1)
    MLR_opt_YTest_Yield(i)    = MLR_opt_data.out{i}(6);
    MLR_opt_YTest_Impurity(i) = MLR_opt_data.out{i}(5);

    GPR_opt_YTest_Yield(i)    = GPR_opt_data.out{i}(6);
    GPR_opt_YTest_Impurity(i) = GPR_opt_data.out{i}(5);

    ANN_opt_YTest_Yield(i)    = ANN_opt_data.out{i}(6);
    ANN_opt_YTest_Impurity(i) = ANN_opt_data.out{i}(5);
end

MLR_opt_cost          = zeros(6, 1);
MLR_opt_YTest_ActCost = zeros(6, 1);
GPR_opt_cost          = zeros(6, 1);
GPR_opt_YTest_ActCost = zeros(6, 1);
ANN_opt_cost          = zeros(6, 1);
ANN_opt_YTest_ActCost = zeros(6, 1);
for i = 1:length(price)
    MLR_opt_cost(2*i-1)          = (1 + price(i)*MLR_opt(2*i-1,4)) / predict(MLR_Yield, MLR_opt(2*i-1,1:4));
    MLR_opt_cost(2*i)            = (1 + price(i)*MLR_opt(2*i,4)) / predict(MLR_Yield, MLR_opt(2*i,1:4));
    MLR_opt_YTest_ActCost(2*i-1) = (1 + price(i)*MLR_opt(2*i-1,4)) / MLR_opt_YTest_Yield(2*i-1);
    MLR_opt_YTest_ActCost(2*i)   = (1 + price(i)*MLR_opt(2*i,4))   / MLR_opt_YTest_Yield(2*i);

    GPR_opt_cost(2*i-1)          = (1 + price(i)*GPR_opt(2*i-1,4)) / predict(gprMdl_Yield, GPR_opt(2*i-1,1:4));
    GPR_opt_cost(2*i)            = (1 + price(i)*GPR_opt(2*i,4)) / predict(gprMdl_Yield, GPR_opt(2*i,1:4));
    GPR_opt_YTest_ActCost(2*i-1) = (1 + price(i)*GPR_opt(2*i-1,4)) / GPR_opt_YTest_Yield(2*i-1);
    GPR_opt_YTest_ActCost(2*i)   = (1 + price(i)*GPR_opt(2*i,4))   / GPR_opt_YTest_Yield(2*i);

    ANN_opt_cost(2*i-1)          = (1 + price(i)*ANN_opt(2*i-1,4)) / predict(ANN_Yield, ANN_opt(2*i-1,1:4));
    ANN_opt_cost(2*i)            = (1 + price(i)*ANN_opt(2*i,4)) / predict(ANN_Yield, ANN_opt(2*i,1:4));
    ANN_opt_YTest_ActCost(2*i-1) = (1 + price(i)*ANN_opt(2*i-1,4)) / ANN_opt_YTest_Yield(2*i-1);
    ANN_opt_YTest_ActCost(2*i)   = (1 + price(i)*ANN_opt(2*i,4))   / ANN_opt_YTest_Yield(2*i);
end

disp('True yield with estimated opt vs predicted opt yield (MLR) fminsearchcon')
disp([MLR_opt_YTest_Yield predict(MLR_Yield, MLR_opt(:,1:4))])

disp('True impurity with estimated opt vs predicted opt impurity (MLR) fminsearchcon')
disp([MLR_opt_YTest_Impurity predict(MLR_Impurity, MLR_opt(:,1:4))])

disp('True yield with estimated opt vs predicted opt yield (GPR) fminsearchcon')
disp([GPR_opt_YTest_Yield predict(gprMdl_Yield, GPR_opt(:,1:4))])

disp('True impurity with estimated opt vs predicted opt impurity (GPR) fminsearchcon')
disp([GPR_opt_YTest_Impurity predict(gprMdl_Tri, GPR_opt(:,1:4))])

disp('True yield with estimated opt vs predicted opt yield (ANN) fminsearchcon')
disp([ANN_opt_YTest_Yield predict(ANN_Yield, ANN_opt(:,1:4))])

disp('True impurity with estimated opt vs predicted opt impurity (ANN) fminsearchcon')
disp([ANN_opt_YTest_Impurity predict(ANN_Impurity, ANN_opt(:,1:4))])

disp('True cost with estimated opt vs predicted opt cost (MLR) fminsearchcon')
disp([MLR_opt_YTest_ActCost MLR_opt_cost])

disp('True cost with estimated opt vs predicted opt cost (GPR) fminsearchcon')
disp([GPR_opt_YTest_ActCost GPR_opt_cost])

disp('True cost with estimated opt vs predicted opt cost (ANN) fminsearchcon')
disp([ANN_opt_YTest_ActCost ANN_opt_cost])

%% Display info in a smarter way
% Define best optimum computed with fminsearch.con for case I, II and III
% BestOpt = [39.9983 11.9920 44.0064 2.8324; ...
%            39.9817 11.0742 49.6572 2.0782; ...
%            39.9457 10.9681 49.1103 2.0519];
% 
BestOpt = [39.9995 11.9054 29.9989 2.9998; ...
           40.0000 11.0715 50.0000 2.0773; ...
           40.0000 10.9662 49.9999 2.0504];

%% Display info for tables in Overleaf - fminsearchcon
% Predicted optimums in each case
for i = 1:3
    disp(strcat('Optimum predicted (MLR | GPR | ANN) case ', num2str(i)))

    disp([MLR_opt(2*i,1:4)' GPR_opt(2*i,1:4)' ANN_opt(2*i,1:4)'; ...
     predict(MLR_Yield, MLR_opt(2*i,1:4)) predict(gprMdl_Yield, GPR_opt(2*i,1:4)) ...
     predict(ANN_Yield, ANN_opt(2*i,1:4)); ...
     predict(MLR_Impurity, MLR_opt(2*i,1:4)) predict(gprMdl_Tri, GPR_opt(2*i,1:4)) ...
     predict(ANN_Impurity, ANN_opt(2*i,1:4)); ...
     MLR_opt_cost(2*i) GPR_opt_cost(2*i) ANN_opt_cost(2*i)])
end


% Difference between predicted and true minima
for i = 1:3
    disp(strcat('Difference between predicted and actual minima (MLR | GPR | ANN) case ', num2str(i)))
    disp([MLR_opt(2*i,1:4)' GPR_opt(2*i,1:4)' ANN_opt(2*i,1:4)'] - BestOpt(i,:)')
    disp([norm(MLR_opt(2*i,1:4)' - BestOpt(i,:)') ...
          norm(GPR_opt(2*i,1:4)' - BestOpt(i,:)') ...
          norm(ANN_opt(2*i,1:4)' - BestOpt(i,:)')])
end


% Difference between predicted and true minima for yD, yH and cost
for i = 1:3
    disp(strcat('Difference between optima for yD, yH and cost (MLR | GPR | ANN) case ', num2str(i)))
    
    disp([MLR_opt_YTest_Yield(2*i) MLR_opt_YTest_Impurity(2*i)  MLR_opt_YTest_ActCost(2*i) ...
        MLR_opt_YTest_Yield(2*i)-predict(MLR_Yield, MLR_opt(2*i,1:4)) ...
        MLR_opt_YTest_Impurity(2*i)-predict(MLR_Impurity, MLR_opt(2*i,1:4)) ...
        MLR_opt_YTest_ActCost(2*i)-MLR_opt_cost(2*i); ...
        GPR_opt_YTest_Yield(2*i) GPR_opt_YTest_Impurity(2*i)  GPR_opt_YTest_ActCost(2*i) ...
        GPR_opt_YTest_Yield(2*i)-predict(gprMdl_Yield, GPR_opt(2*i,1:4)) ...
        GPR_opt_YTest_Impurity(2*i)-predict(gprMdl_Tri, GPR_opt(2*i,1:4)) ...
        GPR_opt_YTest_ActCost(2*i)-GPR_opt_cost(2*i); ...
        ANN_opt_YTest_Yield(2*i) ANN_opt_YTest_Impurity(2*i)  ANN_opt_YTest_ActCost(2*i) ...
        ANN_opt_YTest_Yield(2*i)-predict(ANN_Yield, ANN_opt(2*i,1:4)) ...
        ANN_opt_YTest_Impurity(2*i)-predict(ANN_Impurity, ANN_opt(2*i,1:4)) ...
        ANN_opt_YTest_ActCost(2*i)-ANN_opt_cost(2*i)])

%      GPR_opt_YTest_Yield(2*i)-predict(gprMdl_Yield, GPR_opt(2*i,1:4)) ...
%      ANN_opt_YTest_Yield(2*i)-predict(ANN_Yield, ANN_opt(2*i,1:4))])
end

%% Display info for tables in Overleaf - fmincon.m
% Predicted optimums in each case
for i = 1:3
    disp(strcat('Optimum predicted (MLR | GPR | ANN) case ', num2str(i)))

    disp([MLR_opt_fmincon(2*i,1:4)' GPR_opt_fmincon(2*i,1:4)' ANN_opt_fmincon(2*i,1:4)'; ...
     predict(MLR_Yield, MLR_opt_fmincon(2*i,1:4)) predict(gprMdl_Yield, GPR_opt_fmincon(2*i,1:4)) ...
     predict(ANN_Yield, ANN_opt_fmincon(2*i,1:4)); ...
     predict(MLR_Impurity, MLR_opt_fmincon(2*i,1:4)) predict(gprMdl_Tri, GPR_opt_fmincon(2*i,1:4)) ...
     predict(ANN_Impurity, ANN_opt_fmincon(2*i,1:4)); ...
     MLR_opt_cost_fmincon(2*i) GPR_opt_cost_fmincon(2*i) ANN_opt_cost_fmincon(2*i)])
end


% Difference between predicted and true minima
for i = 1:3
    disp(strcat('Difference between predicted and actual minima (MLR | GPR | ANN) case ', num2str(i)))
    disp([MLR_opt_fmincon(2*i,1:4)' GPR_opt_fmincon(2*i,1:4)' ANN_opt_fmincon(2*i,1:4)'] - BestOpt(i,:)')
    disp([norm(MLR_opt_fmincon(2*i,1:4)' - BestOpt(i,:)') ...
          norm(GPR_opt_fmincon(2*i,1:4)' - BestOpt(i,:)') ...
          norm(ANN_opt_fmincon(2*i,1:4)' - BestOpt(i,:)')])
end

% Difference between predicted and true minima for yD, yH and cost
for i = 1:3
    disp(strcat('Difference between optima for yD, yH and cost (MLR | GPR | ANN) case ', num2str(i)))
    
    disp([MLR_opt_YTest_Yield_fmincon(2*i) MLR_opt_YTest_Impurity_fmincon(2*i)  MLR_opt_YTest_ActCost_fmincon(2*i) ...
        MLR_opt_YTest_Yield_fmincon(2*i)-predict(MLR_Yield, MLR_opt_fmincon(2*i,1:4)) ...
        MLR_opt_YTest_Impurity_fmincon(2*i)-predict(MLR_Impurity, MLR_opt_fmincon(2*i,1:4)) ...
        MLR_opt_YTest_ActCost_fmincon(2*i)-MLR_opt_cost_fmincon(2*i); ...
        GPR_opt_YTest_Yield_fmincon(2*i) GPR_opt_YTest_Impurity_fmincon(2*i)  GPR_opt_YTest_ActCost_fmincon(2*i) ...
        GPR_opt_YTest_Yield_fmincon(2*i)-predict(gprMdl_Yield, GPR_opt_fmincon(2*i,1:4)) ...
        GPR_opt_YTest_Impurity_fmincon(2*i)-predict(gprMdl_Tri, GPR_opt_fmincon(2*i,1:4)) ...
        GPR_opt_YTest_ActCost_fmincon(2*i)-GPR_opt_cost_fmincon(2*i); ...
        ANN_opt_YTest_Yield_fmincon(2*i) ANN_opt_YTest_Impurity_fmincon(2*i)  ANN_opt_YTest_ActCost_fmincon(2*i) ...
        ANN_opt_YTest_Yield_fmincon(2*i)-predict(ANN_Yield, ANN_opt_fmincon(2*i,1:4)) ...
        ANN_opt_YTest_Impurity_fmincon(2*i)-predict(ANN_Impurity, ANN_opt_fmincon(2*i,1:4)) ...
        ANN_opt_YTest_ActCost_fmincon(2*i)-ANN_opt_cost_fmincon(2*i)])

end


%% Compute Standardized Regression Coefficients (SRCs)
% Use Monte Carlo sampling for computing SRCs
lows      = [5 10 5 1];     % Low values
highs     = [40 12 50 3];   % High values
N_Samples = 1e4;            % Number of sample points

% Generate test data
rng(20, 'twister')
SRCsPlan = zeros(N_Samples, 4);
for i = 1:4
    SRCsPlan(:,i) = icdf('Uniform', lhsdesign(N_Samples,1), lows(i), highs(i));
end

% 
Data_MM.X  = SRCsPlan;
Data_MLR.X = SRCsPlan;
Data_GPR.X = SRCsPlan;
Data_ANN.X = SRCsPlan;

% Compute SRCs for MM
YD_SRCs_MM = zeros(size(SRCsPlan,1),1);
YH_SRCs_MM = zeros(size(SRCsPlan,1),1);
for i = 1:size(Data_MM.X,1)
    %
    reacstruc.process.T       = SRCsPlan(i,1);  % Celcius
    reacstruc.process.pH      = SRCsPlan(i,2);  % Logarithmic pH scale
    reacstruc.process.Co      = SRCsPlan(i,3);  % g/L
    reacstruc.process.lambda0 = SRCsPlan(i,4);  % mol /mol
    reacstruc.process.tdose   = 30;              % min

    % Run sim with new parameters
    reacstruc = reacsim(reacstruc);

    % Extract yield and triacylated components
    YD_SRCs_MM(i) = reacstruc.out.y(end,7);
    YH_SRCs_MM(i) = reacstruc.out.y(end,10);
end

Data_MM.Y = YD_SRCs_MM;
[SRCs_Yield_MM, results] = easyGSA('UserData', Data_MM, 'Method', 'SRC');
Data_MM.Y = YH_SRCs_MM;
[SRCs_Impurity_MM, results] = easyGSA('UserData', Data_MM, 'Method', 'SRC');

% Compute SRCs for MLR
% Data.Y = YD_MLR;
Data_MLR.Y = predict(MLR_Yield, Data_MLR.X);
[SRCs_Yield_MLR, results] = easyGSA('UserData', Data_MLR, 'Method', 'SRC');
% Data.Y = YH_MLR;
Data_MLR.Y = predict(MLR_Impurity, Data_MLR.X);
[SRCs_Impurity_MLR, results] = easyGSA('UserData', Data_MLR, 'Method', 'SRC');

% Compute SRCs for GPR
% Data.Y = YD_GPR;
Data_GPR.Y = predict(gprMdl_Yield, Data_GPR.X);
[SRCs_Yield_GPR, results] = easyGSA('UserData', Data_GPR, 'Method', 'SRC');
% Data.Y = YH_GPR;
Data_GPR.Y = predict(gprMdl_Tri, Data_GPR.X);
[SRCs_Impurity_GPR, results] = easyGSA('UserData', Data_GPR, 'Method', 'SRC');

% Compute SRCs for ANN
% Data.Y = YD_ANN;
Data_ANN.Y = predict(ANN_Yield, Data_ANN.X);
[SRCs_Yield_ANN, results] = easyGSA('UserData', Data_ANN, 'Method', 'SRC');
% Data.Y = YH_ANN;
Data_ANN.Y = predict(ANN_Impurity, Data_ANN.X);
[SRCs_Impurity_ANN, results] = easyGSA('UserData', Data_ANN, 'Method', 'SRC');

disp('SRCs:')
disp([SRCs_Yield_MM  SRCs_Impurity_MM])
disp([SRCs_Yield_MLR SRCs_Impurity_MLR])
disp([SRCs_Yield_GPR SRCs_Impurity_GPR])
disp([SRCs_Yield_ANN SRCs_Impurity_ANN])

disp('Squared SRCs:')
disp([SRCs_Yield_MLR.^2 SRCs_Impurity_MLR.^2])
disp([SRCs_Yield_GPR.^2 SRCs_Impurity_GPR.^2])
disp([SRCs_Yield_ANN.^2 SRCs_Impurity_ANN.^2])

disp('Sum of squared SRCs:')
disp([sum(SRCs_Yield_MLR.^2) sum(SRCs_Impurity_MLR.^2)])
disp([sum(SRCs_Yield_GPR.^2) sum(SRCs_Impurity_GPR.^2)])
disp([sum(SRCs_Yield_ANN.^2) sum(SRCs_Impurity_ANN.^2)])


% No noise
% SRCs:
%     0.0968   -0.0989
%    -0.1053   -0.4456
%     0.1294    0.1835
%     0.8547    0.5499
% 
%     0.0399   -0.0445
%    -0.1067   -0.4130
%     0.1194    0.1746
%     0.8908    0.5288
% 
%     0.0180    0.0049
%    -0.0506   -0.4728
%     0.0759    0.2052
%     0.8661    0.5429
% 
%     0.0836   -0.1257
%    -0.0826   -0.4650
%     0.0373    0.2330
%     0.9211    0.5440
% 
% Squared SRCs:
%     0.0016    0.0020
%     0.0114    0.1706
%     0.0142    0.0305
%     0.7936    0.2797
% 
%     0.0003    0.0000
%     0.0026    0.2236
%     0.0058    0.0421
%     0.7502    0.2947
% 
%     0.0070    0.0158
%     0.0068    0.2162
%     0.0014    0.0543
%     0.8484    0.2959
% 
% Sum of squared SRCs:
%     0.8208    0.4827
% 
%     0.7588    0.5604
% 
%     0.8636    0.5822


% Noise
% SRCs:
%     0.1631   -0.0907
%    -0.0890   -0.4638
%     0.1411    0.1846
%     0.8539    0.5388
% 
%     0.1231   -0.0376
%    -0.0578   -0.4060
%     0.0411    0.1938
%     0.8789    0.5246
% 
%     0.0903    0.0049
%    -0.0261   -0.4699
%     0.0089    0.2250
%     0.8553    0.5386
% 
%     0.2473   -0.1472
%    -0.0700   -0.4291
%     0.0562    0.2922
%     0.9508    0.5185

%% OFAT and plot
reacstruc = reacstruccreate();

OFAT_setpoint = [20 11 20 2]; % [T pH Co lambda0 tdose]

OFAT_T       = 5:1.25:40;
OFAT_pH      = 10:0.05:12;
OFAT_Co      = 5:1.25:50;
OFAT_lambda0 = 1:0.05:3;

% Create experimental plans 
OFAT_plan_T       = repmat(OFAT_setpoint, length(OFAT_T), 1);
OFAT_plan_pH      = repmat(OFAT_setpoint, length(OFAT_pH), 1);
OFAT_plan_Co      = repmat(OFAT_setpoint, length(OFAT_Co), 1);
OFAT_plan_lambda0 = repmat(OFAT_setpoint, length(OFAT_lambda0), 1);

OFAT_plan_T(:,1)       = OFAT_T;
OFAT_plan_pH(:,2)      = OFAT_pH;
OFAT_plan_Co(:,3)      = OFAT_Co;
OFAT_plan_lambda0(:,4) = OFAT_lambda0;

% Reset at setpoint
reacstruc.process.T       = 20; % Celcius
reacstruc.process.pH      = 11; % Logarithmic pH scale
reacstruc.process.Co      = 20; % g/L
reacstruc.process.lambda0 = 2;  % mol /mol
reacstruc.process.tdose   = 30; % g/L

% Pre-allocate
yD_OFAT_T       = zeros(length(OFAT_T),1);
yD_OFAT_pH      = zeros(length(OFAT_pH),1);
yD_OFAT_Co      = zeros(length(OFAT_Co),1);
yD_OFAT_lambda0 = zeros(length(OFAT_lambda0),1);
yH_OFAT_T       = zeros(length(OFAT_T),1);
yH_OFAT_pH      = zeros(length(OFAT_pH),1);
yH_OFAT_Co      = zeros(length(OFAT_Co),1);
yH_OFAT_lambda0 = zeros(length(OFAT_lambda0),1);

% Run simulations and change only temperature
for i = 1:length(OFAT_T)
    reacstruc.process.T       = OFAT_T(i); % Celcius

    % Run sim with new parameters
    reacstruc = reacsim(reacstruc);

    % Extract yield and triacylated components
    yD_OFAT_T(i,1) = reacstruc.out.y(end,7);
    yH_OFAT_T(i,1) = reacstruc.out.y(end,10);
end

% Run simulations and change only pH
reacstruc.process.T = 20;   % Reset temperature
for i = 1:length(OFAT_pH)
    reacstruc.process.pH       = OFAT_pH(i); % Celcius

    % Run sim with new parameters
    reacstruc = reacsim(reacstruc);

    % Extract yield and triacylated components
    yD_OFAT_pH(i,1) = reacstruc.out.y(end,7);
    yH_OFAT_pH(i,1) = reacstruc.out.y(end,10);
end

% Run simulations and change only Co
reacstruc.process.pH = 11;   % Reset temperature
for i = 1:length(OFAT_Co)
    reacstruc.process.Co       = OFAT_Co(i); % Celcius

    % Run sim with new parameters
    reacstruc = reacsim(reacstruc);

    % Extract yield and triacylated components
    yD_OFAT_Co(i,1) = reacstruc.out.y(end,7);
    yH_OFAT_Co(i,1) = reacstruc.out.y(end,10);
end

% Run simulations and change only Co
reacstruc.process.pH = 11;   % Reset temperature
for i = 1:length(OFAT_lambda0)
    reacstruc.process.lambda0 = OFAT_lambda0(i); % Celcius

    % Run sim with new parameters
    reacstruc = reacsim(reacstruc);

    % Extract yield and triacylated components
    yD_OFAT_lambda0(i,1) = reacstruc.out.y(end,7);
    yH_OFAT_lambda0(i,1) = reacstruc.out.y(end,10);
end

% Predictions
[OFAT_T_MLR_Yield, OFAT_MLR_T_Yield_PI] = predict(MLR_Yield, OFAT_plan_T);
[OFAT_pH_MLR_Yield, OFAT_MLR_pH_Yield_PI] = predict(MLR_Yield, OFAT_plan_pH);
[OFAT_Co_MLR_Yield, OFAT_MLR_Co_Yield_PI] = predict(MLR_Yield, OFAT_plan_Co);
[OFAT_lambda0_MLR_Yield, OFAT_MLR_lambda0_Yield_PI] = predict(MLR_Yield, OFAT_plan_lambda0);

[OFAT_T_GPR_Yield, ~,  OFAT_GPR_T_Yield_PI] = predict(gprMdl_Yield, OFAT_plan_T);
[OFAT_pH_GPR_Yield, ~, OFAT_GPR_pH_Yield_PI] = predict(gprMdl_Yield, OFAT_plan_pH);
[OFAT_Co_GPR_Yield, ~, OFAT_GPR_Co_Yield_PI] = predict(gprMdl_Yield, OFAT_plan_Co);
[OFAT_lambda0_GPR_Yield, ~, OFAT_GPR_lambda0_Yield_PI] = predict(gprMdl_Yield, OFAT_plan_lambda0);

OFAT_T_ANN_Yield       = predict(ANN_Yield, OFAT_plan_T);
OFAT_pH_ANN_Yield      = predict(ANN_Yield, OFAT_plan_pH);
OFAT_Co_ANN_Yield      = predict(ANN_Yield, OFAT_plan_Co);
OFAT_lambda0_ANN_Yield = predict(ANN_Yield, OFAT_plan_lambda0);

% Compute statistical measures
OFAT_T_stats_MLR       = rs_stats(yD_OFAT_T, OFAT_T_MLR_Yield);
OFAT_pH_stats_MLR      = rs_stats(yD_OFAT_pH, OFAT_pH_MLR_Yield);
OFAT_Co_stats_MLR      = rs_stats(yD_OFAT_Co, OFAT_Co_MLR_Yield);
OFAT_lambda0_stats_MLR = rs_stats(yD_OFAT_lambda0, OFAT_lambda0_MLR_Yield);

OFAT_T_stats_GPR       = rs_stats(yD_OFAT_T, OFAT_T_GPR_Yield);
OFAT_pH_stats_GPR      = rs_stats(yD_OFAT_pH, OFAT_pH_GPR_Yield);
OFAT_Co_stats_GPR      = rs_stats(yD_OFAT_Co, OFAT_Co_GPR_Yield);
OFAT_lambda0_stats_GPR = rs_stats(yD_OFAT_lambda0, OFAT_lambda0_GPR_Yield);

OFAT_T_stats_ANN       = rs_stats(yD_OFAT_T, OFAT_T_ANN_Yield);
OFAT_pH_stats_ANN      = rs_stats(yD_OFAT_pH, OFAT_pH_ANN_Yield);
OFAT_Co_stats_ANN      = rs_stats(yD_OFAT_Co, OFAT_Co_ANN_Yield);
OFAT_lambda0_stats_ANN = rs_stats(yD_OFAT_lambda0, OFAT_lambda0_ANN_Yield);

% Plots
OFAT_Ylim_MAX = 1.25;

figure
subplot(2,2,1); hold all
plot(OFAT_T, yD_OFAT_T, '-')
plot(OFAT_T, OFAT_T_MLR_Yield, 'x')
patch([OFAT_T flip(OFAT_T)], ...
    [OFAT_MLR_T_Yield_PI(:,1)' flip(OFAT_MLR_T_Yield_PI(:,2))'], ...
      'blue', 'FaceAlpha', .125) 
xlabel('Temperature', 'FontSize', 10)
ylabel('Yield', 'FontSize', 10)
ylim([0 OFAT_Ylim_MAX])
set(gca, 'FontSize', 8)
legend("Ground truth", "Predictions", "95 \% CI", 'location', 'southeast')
% text(5, 1.15, strcat("$R^{2}$: ", num2str(round(OFAT_T_stats_MLR.R2,4))), 'FontSize', 10)
% text(5, 1.05, strcat("RMSE: ", num2str(round(OFAT_T_stats_MLR.RMSE,4))), 'FontSize', 10)

subplot(2,2,2); hold all
plot(OFAT_pH, yD_OFAT_pH, '-')
plot(OFAT_pH, OFAT_pH_MLR_Yield, 'x')
patch([OFAT_pH flip(OFAT_pH)], ...
    [OFAT_MLR_pH_Yield_PI(:,1)' flip(OFAT_MLR_pH_Yield_PI(:,2))'], ...
      'blue', 'FaceAlpha', .125) 
xlabel('pH', 'FontSize', 10)
ylabel('Yield', 'FontSize', 10)
ylim([0 OFAT_Ylim_MAX])
set(gca, 'FontSize', 8)
legend("Ground truth", "Predictions", "95 \% CI", 'location', 'southeast')
% text(10.25, 1.15, strcat("$R^{2}$: ", num2str(round(OFAT_pH_stats_MLR.R2,4))), 'FontSize', 10)
% text(10.25, 1.05, strcat("RMSE: ", num2str(round(OFAT_pH_stats_MLR.RMSE,4))), 'FontSize', 10)

subplot(2,2,3); hold all
plot(OFAT_Co, yD_OFAT_Co, '-')
plot(OFAT_Co, OFAT_Co_MLR_Yield, 'x')
patch([OFAT_Co flip(OFAT_Co)], ...
    [OFAT_MLR_Co_Yield_PI(:,1)' flip(OFAT_MLR_Co_Yield_PI(:,2))'], ...
      'blue', 'FaceAlpha', .125) 
xlabel('Initial concentration of A, $C_{A,0}$', 'FontSize', 10)
ylabel('Yield', 'FontSize', 10)
ylim([0 OFAT_Ylim_MAX])
set(gca, 'FontSize', 8)
legend("Ground truth", "Predictions", "95 \% CI", 'location', 'southeast')
% text(5, 1.15, strcat("$R^{2}$: ", num2str(round(OFAT_Co_stats_MLR.R2,4))), 'FontSize', 10)
% text(5, 1.05, strcat("RMSE: ", num2str(round(OFAT_Co_stats_MLR.RMSE,4))), 'FontSize', 10)

subplot(2,2,4); hold all
plot(OFAT_lambda0, yD_OFAT_lambda0, '-')
plot(OFAT_lambda0, OFAT_lambda0_MLR_Yield, 'x')
patch([OFAT_lambda0 flip(OFAT_lambda0)], ...
    [OFAT_MLR_lambda0_Yield_PI(:,1)' flip(OFAT_MLR_lambda0_Yield_PI(:,2))'], ...
      'blue', 'FaceAlpha', .125) 
xlabel('Ratio of S and A, $\lambda_{0}$', 'FontSize', 10)
ylabel('Yield', 'FontSize', 10)
ylim([0 OFAT_Ylim_MAX])
set(gca, 'FontSize', 8)
legend("Ground truth", "Predictions", "95 \% CI", 'location', 'southeast')
% text(1.25, 1.15, strcat("$R^{2}$: ", num2str(round(OFAT_lambda0_stats_MLR.R2,4))), 'FontSize', 10)
% text(1.25, 1.05, strcat("RMSE: ", num2str(round(OFAT_lambda0_stats_MLR.RMSE,4))), 'FontSize', 10)
sgtitle('LHS - MLR predictions')
FigureTitle(10) = "FitLHS_OFAT_MLR";


figure
subplot(2,2,1); hold all
plot(OFAT_T, yD_OFAT_T, '-')
plot(OFAT_T, predict(gprMdl_Yield, OFAT_plan_T), 'x')
patch([OFAT_T flip(OFAT_T)], ...
    [OFAT_GPR_T_Yield_PI(:,1)' flip(OFAT_GPR_T_Yield_PI(:,2))'], ...
      'blue', 'FaceAlpha', .125) 
xlabel('Temperature', 'FontSize', 10)
ylabel('Yield', 'FontSize', 10)
ylim([0 OFAT_Ylim_MAX])
set(gca, 'FontSize', 8)
legend("Ground truth", "Predictions", "95 \% CI", 'location', 'southeast')
% text(5, 1.15, strcat("$R^{2}$: ", num2str(round(OFAT_T_stats_GPR.R2,4))), 'FontSize', 10)
% text(5, 1.05, strcat("RMSE: ", num2str(round(OFAT_T_stats_GPR.RMSE,4))), 'FontSize', 10)

subplot(2,2,2); hold all
plot(OFAT_pH, yD_OFAT_pH, '-')
plot(OFAT_pH, predict(gprMdl_Yield, OFAT_plan_pH), 'x')
patch([OFAT_pH flip(OFAT_pH)], ...
    [OFAT_GPR_pH_Yield_PI(:,1)' flip(OFAT_GPR_pH_Yield_PI(:,2))'], ...
      'blue', 'FaceAlpha', .125) 
xlabel('pH', 'FontSize', 10)
ylabel('Yield', 'FontSize', 10)
ylim([0 OFAT_Ylim_MAX])
set(gca, 'FontSize', 8)
legend("Ground truth", "Predictions", "95 \% CI", 'location', 'southeast')
% text(10.25, 1.15, strcat("$R^{2}$: ", num2str(round(OFAT_pH_stats_GPR.R2,4))), 'FontSize', 10)
% text(10.25, 1.05, strcat("RMSE: ", num2str(round(OFAT_pH_stats_GPR.RMSE,4))), 'FontSize', 10)

subplot(2,2,3); hold all
plot(OFAT_Co, yD_OFAT_Co, '-')
plot(OFAT_Co, predict(gprMdl_Yield, OFAT_plan_Co), 'x')
patch([OFAT_Co flip(OFAT_Co)], ...
    [OFAT_GPR_Co_Yield_PI(:,1)' flip(OFAT_GPR_Co_Yield_PI(:,2))'], ...
      'blue', 'FaceAlpha', .125) 
xlabel('Initial concentration of A, $C_{A,0}$', 'FontSize', 10)
ylabel('Yield', 'FontSize', 10)
ylim([0 OFAT_Ylim_MAX])
set(gca, 'FontSize', 8)
legend("Ground truth", "Predictions", "95 \% CI", 'location', 'southeast')
% text(5, 1.15, strcat("$R^{2}$: ", num2str(round(OFAT_Co_stats_GPR.R2,4))), 'FontSize', 10)
% text(5, 1.05, strcat("RMSE: ", num2str(round(OFAT_Co_stats_GPR.RMSE,4))), 'FontSize', 10)

subplot(2,2,4); hold all
plot(OFAT_lambda0, yD_OFAT_lambda0, '-')
plot(OFAT_lambda0, predict(gprMdl_Yield, OFAT_plan_lambda0), 'x')
patch([OFAT_lambda0 flip(OFAT_lambda0)], ...
    [OFAT_GPR_lambda0_Yield_PI(:,1)' flip(OFAT_GPR_lambda0_Yield_PI(:,2))'], ...
      'blue', 'FaceAlpha', .125) 
xlabel('Ratio of S and A, $\lambda_{0}$', 'FontSize', 10)
ylabel('Yield', 'FontSize', 10)
ylim([0 OFAT_Ylim_MAX])
set(gca, 'FontSize', 8)
sgtitle('LHS - GPR predictions')
legend("Ground truth", "Predictions", "95 \% CI", 'location', 'southeast')
% text(1.25, 1.15, strcat("$R^{2}$: ", num2str(round(OFAT_lambda0_stats_GPR.R2,4))), 'FontSize', 10)
% text(1.25, 1.05, strcat("RMSE: ", num2str(round(OFAT_lambda0_stats_GPR.RMSE,4))), 'FontSize', 10)
FigureTitle(11) = "FitLHS_OFAT_GPR";

figure
subplot(2,2,1); hold all
plot(OFAT_T, yD_OFAT_T, '-')
plot(OFAT_T, predict(ANN_Yield, OFAT_plan_T), 'x')
xlabel('Temperature', 'FontSize', 10)
ylabel('Yield', 'FontSize', 10)
ylim([0 OFAT_Ylim_MAX])
set(gca, 'FontSize', 8)
legend("Ground truth", "Predictions", 'location', 'southeast')
% text(5, 1.15, strcat("$R^{2}$: ", num2str(round(OFAT_T_stats_ANN.R2,4))), 'FontSize', 10)
% text(5, 1.05, strcat("RMSE: ", num2str(round(OFAT_T_stats_ANN.RMSE,4))), 'FontSize', 10)

subplot(2,2,2); hold all
plot(OFAT_pH, yD_OFAT_pH, '-')
plot(OFAT_pH, predict(ANN_Yield, OFAT_plan_pH), 'x')
xlabel('pH', 'FontSize', 10)
ylabel('Yield', 'FontSize', 10)
ylim([0 OFAT_Ylim_MAX])
set(gca, 'FontSize', 8)
legend("Ground truth", "Predictions", 'location', 'southeast')
% text(10.25, 1.15, strcat("$R^{2}$: ", num2str(round(OFAT_pH_stats_ANN.R2,4))), 'FontSize', 10)
% text(10.25, 1.05, strcat("RMSE: ", num2str(round(OFAT_pH_stats_ANN.RMSE,4))), 'FontSize', 10)

subplot(2,2,3); hold all
plot(OFAT_Co, yD_OFAT_Co, '-')
plot(OFAT_Co, predict(ANN_Yield, OFAT_plan_Co), 'x')
xlabel('Initial concentration of A, $C_{A,0}$', 'FontSize', 10)
ylabel('Yield', 'FontSize', 10)
ylim([0 OFAT_Ylim_MAX])
set(gca, 'FontSize', 8)
legend("Ground truth", "Predictions", 'location', 'southeast')
% text(10.25, 1.15, strcat("$R^{2}$: ", num2str(round(OFAT_Co_stats_ANN.R2,4))), 'FontSize', 10)
% text(10.25, 1.05, strcat("RMSE: ", num2str(round(OFAT_Co_stats_ANN.RMSE,4))), 'FontSize', 10)

subplot(2,2,4); hold all
plot(OFAT_lambda0, yD_OFAT_lambda0, '-')
plot(OFAT_lambda0, predict(ANN_Yield, OFAT_plan_lambda0), 'x')
xlabel('Ratio of S and A, $\lambda_{0}$', 'FontSize', 10)
ylabel('Yield', 'FontSize', 10)
ylim([0 OFAT_Ylim_MAX])
set(gca, 'FontSize', 8)
legend("Ground truth", "Predictions", 'location', 'southeast')
% text(1.25, 1.15, strcat("$R^{2}$: ", num2str(round(OFAT_lambda0_stats_ANN.R2,4))), 'FontSize', 10)
% text(1.25, 1.05, strcat("RMSE: ", num2str(round(OFAT_lambda0_stats_ANN.RMSE,4))), 'FontSize', 10)
sgtitle('LHS - ANN predictions')
FigureTitle(12) = "FitLHS_OFAT_ANN";


%% MLR
% Compute predictions over entire input space with N samples
% YD_MLR = predict(MLR_Yield, OptPlan);
% YH_MLR = predict(MLR_Impurity, OptPlan);
% Data.X = Opt_Plan; 
% Data.Y = ym(:,j);
% [SRCs(:,j), results] = easyGSA('UserData',Data,'Method','SRC');
% 
% figure
% b = categorical({'T','pH','Co','$\lambda_{0}$'});
% b = reordercats(b,{'T','pH','Co','$\lambda_{0}$'});
% sgtitle('Sobol indices', 'FontSize', 16)
% subplot(2,1,1)
% bar(b, Si)
% ylabel('First Order Sobol indices', 'FontSize', 16)
% legend("No noise", "Low noise", 'location', 'northwest')
% 
% subplot(2,1,2)
% bar(b, STi)
% ylabel('Total Order Sobol indices', 'FontSize', 16)
% legend("No noise", "Low noise", "High noise", 'location', 'northwest')
% FigureTitle(17) = "DOE_SobolIndices";
% 
% figure
% bar(b, SRCs.^2)
% ylabel('Standardized regression coefficients')
% xlabel('Input factors')
% legend("No noise", "Normal noise", "High noise", 'location', 'northwest')
% FigureTitle(18) = "DOE_StandardizedRegressionCoefficients";


%% GPR
% Compute 
% YD_GPR = predict(gprMdl_Yield, OptPlan);
% YH_GPR = predict(gprMdl_Tri, OptPlan);
% 
% % Find impurity component below the constraint
% k = find(0 < YH_GPR & YH_GPR <= 0.005 & 0 < YD_GPR & YD_GPR <= 1);
% 
% % Compute initial point
% [cost, IdxOpt] = min((1 + 0.5370*OptPlan(k,4)) ./ YD_GPR(k));
% 
% % Compute initial point with latin hypercube sampling
% optim_GPR.LB            = lows;
% optim_GPR.UB            = highs;
% optim_GPR.ModelYield    = gprMdl_Yield;
% optim_GPR.ModelImpurity = gprMdl_Tri;
% optim_GPR.price         = 0.5370;
% optim_GPR.z0            = OptPlan(k(IdxOpt),:);
% optim_GPR.glb           = [0; 0];
% optim_GPR.gub           = [1; 0.01];
% GPR_sol                 = runoptim_lm(optim_GPR);
% GPR_opt                 = GPR_sol(:)' .* (highs - lows) + lows;


%%
% reacstruc = reacstruccreate();
% 
% % Compute true optimum
% price = reacstruc.optim.Price.SCrel;
% 
% z0 = pars.setpoints;
% for i = 1:length(reacstruc.optim.var)
%     % Define helper variables
%     LB(i) = reacstruc.optim.LB(i);
%     UB(i) = reacstruc.optim.UB(i);

%     % Compute normalized initial condition
%     z0(i) = eval([reacstruc.optim.var{i}]);
%     x0(i) = (z0(i)-LB(i)) / (UB(i)-LB(i));
% 
%     % Set normalized constraints
%     con(1,i) = 0;
%     con(2,i) = 1;
% end
% fun = @(x) reacoptim(x,reacstruc);
% options = optimset('TolX',1e-3);
% sol = fminsearchcon(fun, x0, con(1,:), con(2,:), [],[],[],options);
% 
% % Define constraints and initial point
% optim.LB = pars.lows;
% optim.UB = pars.highs;
% optim.z0 = pars.setpoints;
% optim.price = reacstruc.optim.Price.SCrel;
% 
% % MLR
% optim.ModelYield = MLR_Yield;
% optim.ModelImpurity = MLR_Impurity;
% MLR_sol = runoptim_lm(optim);
% 
% % GPR
% optim.ModelYield = gprMdl_Yield;
% optim.ModelImpurity = gprMdl_Tri;
% GPR_sol = runoptim_lm(optim);
% 
% % ANN
% optim.ModelYield = ANN_Yield;
% optim.ModelImpurity = ANN_Impurity;
% ANN_sol = runoptim_lm(optim);

%% Display results
% True sol vs MLR vs GPR vs ANN
% disp(sol(:)' .* (reacstruc.optim.UB'-reacstruc.optim.LB') + reacstruc.optim.LB')
% disp(MLR_sol(:)' .* (optim.UB-optim.LB) + optim.LB)
% disp(GPR_sol(:)' .* (optim.UB-optim.LB) + optim.LB)
% disp(ANN_sol(:)' .* (optim.UB-optim.LB) + optim.LB)


%% Save figures
if bool_SaveFigures
    for i = 1:length(FigureTitle) 
        saveas(figure(i), ...
            strcat(pwd, "\Images\FitLHSNoise\", ...
                        FigureTitle(i), "_NoiseLevel_", num2str(SimulationIndex)), 'png')
    end
end

%% Notes
%% No noise
% R2 SSE RMSE AAD
% ans =
% 
%     1.0000    0.0000    0.0000    0.0000
%     0.9204    0.1111    0.0861    0.0673
%     0.9657    0.0386    0.0507    0.0395
%     0.9602    0.0460    0.0554    0.0428
% 
% 
% ans =
% 
%     1.0000    0.0000    0.0000    0.0000
%     0.7834    0.0141    0.0307    0.0219
%     0.9474    0.0036    0.0154    0.0111
%     0.9709    0.0008    0.0073    0.0055

%% Noise
% R2 SSE RMSE AAD
% ans =
% 
%     0.9928    0.0110    0.0270    0.0186
%     0.9030    0.1205    0.0896    0.0742
%     0.9216    0.0834    0.0746    0.0586
%     0.8999    0.1009    0.0820    0.0615
% 
% 
% ans =
% 
%     0.9933    0.0002    0.0037    0.0022
%     0.8225    0.0116    0.0278    0.0191
%     0.9471    0.0051    0.0184    0.0116
%     0.8642    0.0039    0.0161    0.0123


%% Notes - optimization results
%% No noise
% MLR_opt =
% 
%    31.9175   10.7149   28.5127    2.5556   30.0000
%    31.0059   11.4570   44.4661    2.4789   30.0000
%    39.7650   11.2334   49.8959    2.3378   30.0000
%    39.9985   11.2287   49.9190    2.3363   30.0000
%    39.9393   11.2296   49.9707    2.3364   30.0000
%    39.9393   11.2296   49.9707    2.3364   30.0000
% 
% GPR_opt =
% 
%    29.5640   11.2862   46.2529    2.6513   30.0000
%    13.1673   11.2825    8.2987    2.9024   30.0000
%    40.0000   11.0916   50.0000    2.2081   30.0000
%    37.1146   10.8971   23.1313    2.1508   30.0000
%    40.0000   11.0762   50.0000    2.1166   30.0000
%    13.5839   11.1225   48.0528    1.9756   30.0000
% 
% ANN_opt =
% 
%    33.7757   11.1914   43.6507    2.4793   30.0000
%    39.8872   10.9140    5.7091    2.9994   30.0000
%    34.6925   11.1614   45.7675    2.3137   30.0000
%    31.3707   11.3772   49.8199    2.1135   30.0000
%    25.5833   11.0234   37.3533    2.0626   30.0000
%    30.1664   11.2792   49.4746    2.0606   30.0000

%% Noise
% MLR_opt =
% 
%    24.7700   11.2050   22.9735    2.8055   30.0000
%    21.5757   11.8624   42.3806    2.8738   30.0000
%    39.8210   11.2077   39.0726    2.3159   30.0000
%    39.9943   11.2491   41.6135    2.3078   30.0000
%    40.0000   11.2255   42.3327    2.3107   30.0000
%    40.0000   11.2255   42.3327    2.3107   30.0000
% 
% GPR_opt =
% 
%    36.3193   11.5778   32.0269    2.8549   30.0000
%    34.9905   11.5292   18.4632    2.8463   30.0000
%    40.0000   10.8964   50.0000    2.2347   30.0000
%    40.0000   10.8964   10.1666    2.2347   30.0000
%    40.0000   10.7775   50.0000    2.0287   30.0000
%    40.0000   10.7775   13.3884    2.0287   30.0000
% 
% ANN_opt =
% 
%    36.4304   10.2487   47.3283    2.7338   30.0000
%    19.5549   11.9982    8.8674    2.9982   30.0000
%    40.0000   10.9809   49.9957    2.4950   30.0000
%    32.1948   11.3311    5.2463    2.7542   30.0000
%    40.0000   10.0000   50.0000    1.8906   30.0000
%    40.0000   10.4461    5.0144    2.1824   30.0000

%% More optimization results
% No noise
% True yield with estimated opt vs predicted opt yield (MLR)
%     0.9211    1.0000
%     0.9766    1.0000
%     0.9755    1.0000
%     0.9754    1.0000
%     0.9754    1.0000
%     0.9754    1.0000
% 
% True impurity with estimated opt vs predicted opt impurity (MLR)
%     0.0789    0.0735
%     0.0229    0.0017
%     0.0245    0.0000
%     0.0245    0.0000
%     0.0245    0.0000
%     0.0245    0.0000
% 
% True yield with estimated opt vs predicted opt yield (GPR)
%     0.9524    1.0000
%     0.8922    0.9515
%     0.9795    0.9664
%     0.9642    0.9066
%     0.9795    0.9392
%     0.8402    0.8712
% 
% True impurity with estimated opt vs predicted opt impurity (GPR)
%     0.0476    0.0542
%     0.0187    0.0100
%     0.0203    0.0252
%     0.0130    0.0100
%     0.0125    0.0187
%     0.0159    0.0100
% 
% True yield with estimated opt vs predicted opt yield (ANN)
%     0.9626    1.0000
%     0.9800    0.9909
%     0.9730    0.9779
%     0.9418    0.9128
%     0.9194    0.9160
%     0.9224    0.8965
% 
% True impurity with estimated opt vs predicted opt impurity (ANN)
%     0.0374    0.0372
%     0.0182    0.0100
%     0.0269    0.0274
%     0.0098    0.0100
%     0.0139    0.0194
%     0.0096    0.0100
% 
% True cost with estimated opt vs predicted opt cost (MLR)
%     1.0856    1.0000
%     1.0239    1.0000
%     2.3122    2.2555
%     2.3114    2.2547
%     5.8156    5.6728
%     5.8156    5.6728
% 
% True cost with estimated opt vs predicted opt cost (GPR)
%     1.0499    1.0000
%     1.1209    1.0509
%     2.2317    2.2620
%     2.2351    2.3772
%     5.3429    5.5718
%     5.8925    5.6834
% 
% True cost with estimated opt vs predicted opt cost (ANN)
%     1.0389    1.0000
%     1.0204    1.0092
%     2.3048    2.2933
%     2.2671    2.3391
%     5.5744    5.5952
%     5.5518    5.7123

% Noise
% True yield with estimated opt vs predicted opt yield (MLR)
%     0.9597    1.0000
%     0.9709    1.0000
%     0.9805    1.0000
%     0.9815    1.0000
%     0.9802    1.0000
%     0.9802    1.0000
% 
% True impurity with estimated opt vs predicted opt impurity (MLR)
%     0.0402    0.0403
%     0.0200    0.0027
%     0.0193    0.0095
%     0.0182    0.0020
%     0.0196    0.0031
%     0.0196    0.0031
% 
% True yield with estimated opt vs predicted opt yield (GPR)
%     0.9657    1.0000
%     0.9822    1.0000
%     0.9679    0.9258
%     0.9438    0.9258
%     0.9439    0.8695
%     0.8757    0.8695
% 
% True impurity with estimated opt vs predicted opt impurity (GPR)
%     0.0343    0.0494
%     0.0169    0.0094
%     0.0321    0.0252
%     0.0097    0.0100
%     0.0109    0.0203
%     0.0075    0.0100
% 
% True yield with estimated opt vs predicted opt yield (ANN)
%     0.9689    1.0000
%     0.5337    1.0000
%     0.9276    1.0000
%     0.9854    1.0000
%     0.8387    0.8132
%     0.9200    0.8187
% 
% True impurity with estimated opt vs predicted opt impurity (ANN)
%     0.0311    0.0279
%     0.0028    0.0082
%     0.0724    0.0628
%     0.0104    0.0100
%     0.0200    0.0872
%     0.0163    0.0100

%% Smarter optimization notes
% Noisefree
% Optimum predicted (MLR | GPR | ANN) case1
%    31.0059   13.1673   39.8872
%    11.4570   11.2825   10.9140
%    44.4661    8.2987    5.7091
%     2.4789    2.9024    2.9994
%     1.0000    0.9515    0.9909
%     0.0017    0.0100    0.0100
%     1.0000    1.0509    1.0092
% 
% Optimum predicted (MLR | GPR | ANN) case2
%    39.9985   37.1146   31.3707
%    11.2287   10.8971   11.3772
%    49.9190   23.1313   49.8199
%     2.3363    2.1508    2.1135
%     1.0000    0.9066    0.9128
%     0.0000    0.0100    0.0100
%     2.2547    2.3772    2.3391
% 
% Optimum predicted (MLR | GPR | ANN) case3
%    39.9393   13.5839   30.1664
%    11.2296   11.1225   11.2792
%    49.9707   48.0528   49.4746
%     2.3364    1.9756    2.0606
%     1.0000    0.8712    0.8965
%     0.0000    0.0100    0.0100
%     5.6728    5.6834    5.7123
% 
% Difference between predicted and actual minima (MLR | GPR | ANN) case1
%    -8.9924  -26.8310   -0.1111
%    -0.5350   -0.7095   -1.0780
%     0.4597  -35.7077  -38.2973
%    -0.3535    0.0700    0.1670
% 
%     9.0269   44.6705   38.3130
% 
% Difference between predicted and actual minima (MLR | GPR | ANN) case2
%     0.0168   -2.8671   -8.6110
%     0.1545   -0.1771    0.3030
%     0.2618  -26.5259    0.1627
%     0.2581    0.0726    0.0353
% 
%     0.3991   26.6811    8.6179
% 
% Difference between predicted and actual minima (MLR | GPR | ANN) case3
%    -0.0064  -26.3618   -9.7793
%     0.2615    0.1544    0.3111
%     0.8604   -1.0575    0.3643
%     0.2845   -0.0763    0.0087
% 
%     0.9432   26.3836    9.7910
% 
% Difference between optima for yD, yH and cost (MLR | GPR | ANN) case1
%     0.9766    0.0229    1.0239   -0.0234    0.0213    0.0239
%     0.8922    0.0187    1.1209   -0.0594    0.0087    0.0699
%     0.9800    0.0182    1.0204   -0.0109    0.0082    0.0112
% 
% Difference between optima for yD, yH and cost (MLR | GPR | ANN) case2
%     0.9754    0.0245    2.3114   -0.0246    0.0245    0.0567
%     0.9642    0.0130    2.2351    0.0576    0.0030   -0.1421
%     0.9418    0.0098    2.2671    0.0290   -0.0002   -0.0720
% 
% Difference between optima for yD, yH and cost (MLR | GPR | ANN) case3
%     0.9754    0.0245    5.8156   -0.0246    0.0245    0.1428
%     0.8402    0.0159    5.8925   -0.0309    0.0059    0.2091
%     0.9224    0.0096    5.5518    0.0259   -0.0004   -0.1605

% Noise


%% GSA notes
%% SRCs
% No noise
% SRCs:
%     0.0439   -0.0451
%    -0.1015   -0.4144
%     0.1112    0.1742
%     0.8894    0.5345
% 
%     0.0265   -0.0007
%    -0.0460   -0.4736
%     0.0674    0.2049
%     0.8645    0.5438
% 
%     0.1116   -0.1234
%    -0.0855   -0.4603
%     0.0953    0.2476
%     0.9079    0.5309
% 
% Squared SRCs:
%     0.0019    0.0020
%     0.0103    0.1717
%     0.0124    0.0303
%     0.7910    0.2857
% 
%     0.0007    0.0000
%     0.0021    0.2243
%     0.0045    0.0420
%     0.7473    0.2958
% 
%     0.0125    0.0152
%     0.0073    0.2119
%     0.0091    0.0613
%     0.8243    0.2818
% 
% Sum of squared SRCs:
%     0.8156    0.4898
% 
%     0.7546    0.5621
% 
%     0.8531    0.5702


% Noise
% SRCs:
%     0.1263   -0.0379
%    -0.0520   -0.4070
%     0.0311    0.1921
%     0.8784    0.5305
% 
%     0.0884   -0.0006
%    -0.0275   -0.4698
%     0.0001    0.2215
%     0.8586    0.5408
% 
%     0.2549   -0.1186
%    -0.0659   -0.4571
%     0.0553    0.2726
%     0.9503    0.5157
% 
% Squared SRCs:
%     0.0159    0.0014
%     0.0027    0.1656
%     0.0010    0.0369
%     0.7716    0.2815
% 
%     0.0078    0.0000
%     0.0008    0.2207
%     0.0000    0.0491
%     0.7371    0.2925
% 
%     0.0650    0.0141
%     0.0043    0.2089
%     0.0031    0.0743
%     0.9031    0.2659
% 
% Sum of squared SRCs:
%     0.7912    0.4854
% 
%     0.7457    0.5622
% 
%     0.9755    0.5632



















%% Notes
%% NO NOISE
% R2 SSE MSE AAD
% Yield
%     1.0000    0.0000    0.0000    0.0000
%     0.8521    0.4969    0.0099    0.0791
%     0.9276    0.1988    0.0040    0.0449
%     0.9483    0.1494    0.0030    0.0406
% 
%     1.0000    0.0000    0.0000    0.0000
%     0.8379    0.0657    0.0013    0.0285
%     0.9637    0.0270    0.0005    0.0127
%     0.9799    0.0078    0.0002    0.0084

%% NOISE
% Yield
%     0.9978    0.0068    0.0001    0.0079
%     0.8796    0.3421    0.0068    0.0687
%     0.8583    0.4070    0.0081    0.0644
%     0.8271    0.5066    0.0101    0.0789
% 
% Impurity
%     0.9965    0.0037    0.0001    0.0033
%     0.8422    0.0635    0.0013    0.0281
%     0.9684    0.0181    0.0004    0.0119
%     0.9486    0.0188    0.0004    0.0104







































