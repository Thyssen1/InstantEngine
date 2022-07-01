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

%% Box-Behnken Design
% Set experimental parameters
pars.design    = 'bbd';
pars.setpoints = [20 11 20 2]; % [T pH Co lambda0 tdose]
pars.lows      = [5 10 5 1];   % Low values
pars.highs     = [40 12 50 3]; % High values
pars.reps      = 1; 
pars.n         = 4;           % Number of factors
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

%% Perform stepwise multiple linear regression 
% Fit quadratic (2nd order) model
MLR_Yield    = fitlm(X, YD, 'quadratic');
MLR_Impurity = fitlm(X, YH, 'quadratic');
% MLR_Yield    = fitlm([X(:,1) X(:,2) X(:,3) X(:,4).^2], YD, 'quadratic');
% MLR_Impurity = fitlm([X(:,1) X(:,2) X(:,3) X(:,4).^2], YH, 'quadratic');

%% Fit Gaussian Process Regression (GPR) model
% Set kernel function
switch SimulationIndex
    case 1
        Kernel_Yield = 'ardmatern52';
        Kernel_Impurity = 'ardmatern32';
    case 2
        Kernel_Yield = 'ardsquaredexponential';
%         Kernel_Yield = 'ardsquaredexponential';
        Kernel_Impurity = 'ardmatern32'; % ardmatern32
%         Kernel_Impurity = 'ardsquaredexponential';
%         Kernel_Impurity = 'ardmatern52';
end

rng(10, 'twister')
% Step 1: 
% Optimize hyperparameters and compute initial estimate of kernel parameters
gprMdl_Yield = fitrgp(X, YD,             ...
        'KernelFunction', Kernel_Yield, ...
        'OptimizeHyperparameters', 'Sigma',  ...
        'Standardize', true,                 ...
        'HyperparameterOptimizationOptions', ...
        struct('ShowPlots', false,           ...
               'KFold', 5,                   ...
               'MaxObjectiveEvaluations', 50));

rng(11, 'twister')
gprMdl_Tri = fitrgp(X, YH,                     ...    
            'KernelFunction', Kernel_Impurity,  ...                                  
            'OptimizeHyperparameters', 'Sigma',     ...   
            'Standardize', true,                    ...
            'HyperparameterOptimizationOptions',    ...
            struct('ShowPlots', false,              ...
                   'KFold', 5,                      ...
                   'MaxObjectiveEvaluations', 50));

% Step 2: 
% Set initial parameters
% Sigma0_Yield            = gprMdlYield_Init.Sigma;
% Sigma0_lb_Yield         = min(1e-2*std(YD), gprMdlYield_Init.Sigma-1e-3);
% KernelParameters0_Yield = gprMdlYield_Init.KernelInformation.KernelParameters;
% 
% Sigma0_Tri            = gprMdlTri_Init.Sigma;
% Sigma0_lb_Tri         = min(1e-2*std(YH), gprMdlTri_Init.Sigma-1e-3);
% KernelParameters0_Tri = gprMdlTri_Init.KernelInformation.KernelParameters;
% 
% % Step 3: Fit models 
% rng(11, 'twister')
% gprMdl_Yield = fitrgp(X, YD,                                        ...
%                     'KernelFunction', Kernel_Yield,       ...
%                     'KernelParameters', KernelParameters0_Yield,    ...
%                     'Sigma', Sigma0_Yield,                          ...
%                     'SigmaLowerBound', Sigma0_lb_Yield,             ...
%                     'Standardize', true,                            ...
%                     'ConstantSigma', true);
% 
% gprMdl_Tri = fitrgp(X, YH,       ...
%                     'KernelFunction', Kernel_Impurity,   ...
%                     'KernelParameters', KernelParameters0_Tri,  ...
%                     'Sigma', Sigma0_Tri,                        ...
%                     'SigmaLowerBound', Sigma0_lb_Tri,           ...
%                     'Standardize', true,                        ...
%                     'ConstantSigma', true);

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

Neurons_Yield = [40 40];
Neurons_Impurity = [40 40];

% Set activation functions
switch SimulationIndex
    case 1
        Activation_Yield = 'tanh'; % relu
        Activation_Impurity = 'sigmoid'; %relu

        Activation_Yield = 'tanh';
        Activation_Impurity = 'sigmoid';
    case 2
        Activation_Yield = 'tanh';
        Activation_Impurity = 'relu';

        Activation_Yield = 'tanh';
        Activation_Impurity = 'sigmoid';
end

rng(10, 'twister')
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

rng(11, 'twister')
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

% Fit ANN
% rng(13, 'twister')
% ANN_Yield = fitrnet(X, YD,      ...  
%     "LayerSizes", Neurons_Yield,      ...
%     "Standardize", true,        ...
%     "Verbose", 0,               ...
%     'Activations', Activation_Yield,   ...
%     'Lambda', ANN_Yield_Init.ModelParameters.Lambda);
% 
% rng(14, 'twister')
% ANN_Impurity = fitrnet(X, YH,   ...  
%     "LayerSizes", Neurons_Yield,      ...
%     "Standardize", true,        ...
%     "Verbose", 0,               ...
%     'Activations', Activation_Impurity,   ...
%     'Lambda', ANN_Impurity_Init.ModelParameters.Lambda);

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
    yD_act(k,1) = data.out{k}(6);
    yH_act(k,1) = data.out{k}(5);
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


alpha = 0.025;  % Significance level 
tcr   = tinv((1-alpha),DegreesFreedomD); % Critical t-dist value at alpha 

%+-95% confidence intervals
ConfidenceIntervalP95D = [ParameterMinimum' - StandardDeviationD*tcr; ...
                          ParameterMinimum' + StandardDeviationD*tcr];

% Calculate confidence intervals on the model output
n = pars.samples;
% n = 25;
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

%% Linear error propagation method 
% % Set optim options
% options = optimset('display', 'iter', 'tolfun', 1.0e-09,   ...
%                    'tolx', 1.0e-9, 'maxfunevals', 0); % No iteration
%                
% % Compute estimate of Jacobian
% [~,~,residualD,~,~,~,JacobiD] = lsqnonlin(@(x) reacobj_yd(x,reacstruc), ...
%                                         ParameterMinimum, [], [], options);
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
% ConfidenceIntervalP95D = [ParameterMinimum - StandardDeviationD*tcrD; ...
%                           ParameterMinimum + StandardDeviationD*tcrD];
% ConfidenceIntervalP95Tri = [ParameterMinimum - StandardDeviationTri*tcrTri; ...
%                             ParameterMinimum + StandardDeviationTri*tcrTri];
% 
% % Calculate confidence intervals on the model output
% [nD, mD] = size(yD);
% [nTri, mTri] = size(yH);
% 
% % Compute output covariance
% ModelCovarianceD                 = JacobianD * CovarianceParametersD * JacobianD'; % Calculate model covariance # Maybe some statistical reference here?
% ModelStandardDeviationD          = sqrt(diag(ModelCovarianceD));          % Standard deviation of model outputs, vector with all data points
% ModelStandardDeviationReshapedD  = reshape(ModelStandardDeviationD,nD,mD);  % Reshape to matrix with columns corresponding to each set of data points
% ModelCovarianceTri                = JacobianTri * CovarianceParametersTri * JacobianTri'; % Calculate model covariance # Maybe some statistical reference here?
% ModelStandardDeviationTri         = sqrt(diag(ModelCovarianceTri));          % Standard deviation of model outputs, vector with all data points
% ModelStandardDeviationReshapedTri = reshape(ModelStandardDeviationTri,nTri,mTri);  % Reshape to matrix with columns corresponding to each set of data points
% 
% % 95% confidence intervals, calculated only for the measured species (data) Two boundaries, each with a set of m columns 
% ModelConfidenceInterval95D = [yD - ModelStandardDeviationReshapedD*tcrD,...  % Lower boundary
%                              yD + ModelStandardDeviationReshapedD*tcrD];      % Upper boundary 
% ModelConfidenceInterval95Tri = [yH - ModelStandardDeviationReshapedTri*tcrTri,...  % Lower boundary
%                              yH + ModelStandardDeviationReshapedTri*tcrTri];      % Upper boundary 

%% Compute predicted response and 95% confidence interval for MLR model
[MLR_yfit_Yield, MLR_yint_Yield] =  predict(MLR_Yield, X, 'alpha', alpha);
[MLR_yfit_Impurity, MLR_yint_Impurity] =  predict(MLR_Impurity, X, 'alpha', alpha);
%     stats(j) = rs_stats(bbd.Y(:,j), yfit(:,j));

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
MSE_Yield = [stats_MM_Yield.MSE; stats_MLR_Yield.MSE;    ...
            stats_GPR_Yield.MSE; stats_ANN_Yield.MSE];  
RMSE_Yield = [stats_MM_Yield.RMSE; stats_MLR_Yield.RMSE; ...
            stats_GPR_Yield.RMSE; stats_ANN_Yield.RMSE];
AAD_Yield = [stats_MM_Yield.AAD; stats_MLR_Yield.AAD;    ...
            stats_GPR_Yield.AAD; stats_ANN_Yield.AAD]; 

R2_Impurity = [stats_MM_Impurity.R2; stats_MLR_Impurity.R2;       ...
            stats_GPR_Impurity.R2; stats_ANN_Impurity.R2];
SSE_Impurity = [stats_MM_Impurity.SSE; stats_MLR_Impurity.SSE;    ...
            stats_GPR_Impurity.SSE; stats_ANN_Impurity.SSE];  
MSE_Impurity = [stats_MM_Impurity.MSE; stats_MLR_Impurity.MSE;    ...
            stats_GPR_Impurity.MSE; stats_ANN_Impurity.MSE]; 
RMSE_Impurity = [stats_MM_Impurity.RMSE; stats_MLR_Impurity.RMSE; ...
            stats_GPR_Impurity.RMSE; stats_ANN_Impurity.RMSE];
AAD_Impurity = [stats_MM_Impurity.AAD; stats_MLR_Impurity.AAD;    ...
            stats_GPR_Impurity.AAD; stats_ANN_Impurity.AAD]; 

%% Figures of fit to training data (Yield)
[~, IY]  = sort(YD);
[~, IYr] = sort(YD, 'descend');

figure
subplot(2, 2, 1); hold all % Mechanistic 
% plot(YD, YD, 'o')
plot(-1:0.01:2, -1:0.01:2, 'k-')
plot(YD, yD, 'x')
xlim([-0.30 1.25])
ylim([-0.30 1.25])
xlabel('Measured yield ($\%$)', 'FontSize', 10)
ylabel('Yield ($\%$)', 'FontSize', 10)
title('Mechanistic', 'FontSize', 10)
patch([YD(IY)' YD(IYr)'],                     ...
      [ModelConfidenceInterval95D(IY,1)'      ...
       ModelConfidenceInterval95D(IYr,2)'],   ...
      'blue', 'FaceAlpha', .125)
% errorbar(YD, yD, yD-ModelConfidenceInterval95D(:,1), ...
%                  ModelConfidenceInterval95D(:,2)-yD, ...
%                  'LineStyle', 'none', 'Color', 'b')
set(gca, 'FontSize', 11)
legend("", "Predicted", "95\% CI", ...
        'location', 'northwest', 'FontSize', 8, 'box', 'on')
text(0.5,-0.05, strcat("$R^{2}$: ", num2str(round(R2_Yield(1),4)*100), "\%"), 'FontSize', 10)
text(0.5,-0.15, strcat("RMSE: ", num2str(round(RMSE_Yield(1),4))), 'FontSize', 10)

subplot(2, 2, 2); hold all  % DOE (MLR)
% plot(YD, YD, 'o')
plot(-1:0.01:2, -1:0.01:2, 'k-')
plot(YD, MLR_yfit_Yield, 'x')
patch([YD(IY)' YD(IYr)'], [MLR_yint_Yield(IY,1)' MLR_yint_Yield(IYr,2)'], ...
      'blue', 'FaceAlpha', .125) 
xlabel('Measured yield ($\%$)', 'FontSize', 10)
ylabel('Yield ($\%$)', 'FontSize', 10)
xlim([-0.30 1.25])
ylim([-0.30 1.25])
title('MLR', 'FontSize', 10)
set(gca, 'FontSize', 11)
legend("", "Predicted", "95 \% CI", ...
       'FontSize', 8, 'location', 'northwest', 'box', 'on')
text(0.5,-0.05, strcat("$R^{2}$: ", num2str(round(R2_Yield(2),4)*100), "\%"), 'FontSize', 10)
text(0.5,-0.15, strcat("RMSE: ", num2str(round(RMSE_Yield(2),4))), 'FontSize', 10)

subplot(2, 2, 3); hold all % GPR 
% plot(YD, YD, 'o')
plot(-1:0.01:2, -1:0.01:2, 'k-')
plot(YD, GPR_yfit_Yield, 'x')
patch([YD(IY)' YD(IYr)'], [GPR_yint_Yield(IY,1)' GPR_yint_Yield(IYr,2)'], ...
      'blue', 'FaceAlpha', .125) 
xlabel('Measured yield ($\%$)', 'FontSize', 10)
ylabel('Yield ($\%$)', 'FontSize', 10)
title('GPR', 'FontSize', 10)
set(gca, 'FontSize', 11)
xlim([-0.30 1.25])
ylim([-0.30 1.25])
legend("", "GPR", "95 \% CI", ...
       'location', 'northwest', 'FontSize', 8, 'box', 'on')
text(0.5,-0.05, strcat("$R^{2}$: ", num2str(round(R2_Yield(3),4)*100), "\%"), 'FontSize', 10)
text(0.5,-0.15, strcat("RMSE: ", num2str(round(RMSE_Yield(3),4))), 'FontSize', 10)
   
subplot(2, 2, 4); hold all % ANN
% plot(YD, YD, 'o')
plot(-1:0.01:2, -1:0.01:2, 'k-')
plot(YD, ANN_yfit_Yield, 'x')
xlabel('Measured yield ($\%$)', 'FontSize', 10)
ylabel('Yield ($\%$)', 'FontSize', 10)
title('ANN', 'FontSize', 10)
set(gca, 'FontSize', 11)
xlim([-0.30 1.35])
ylim([-0.30 1.35])
legend("", "Predicted", ...
       'location', 'northwest', 'FontSize', 8, 'box', 'on')
FigureTitle(1) = "FitBoxBehnken_TrainingFit_Yield";
sgtitle('BBD - Training')
text(0.5,-0.05, strcat("$R^{2}$: ", num2str(round(R2_Yield(4),4)*100), "\%"), 'FontSize', 10)
text(0.5,-0.15, strcat("RMSE: ", num2str(round(RMSE_Yield(4),4))), 'FontSize', 10)

%% Figures of fit to training data (Impurity)
[~, IH]  = sort(YH);
[~, IHr] = sort(YH, 'descend');

figure
subplot(2, 2, 1); hold all % Mechanistic 
% plot(YH, YH, 'o')
plot(-1:0.01:2, -1:0.01:2, 'k-')
plot(YH, yH, 'x')
% xlim([min(min(min(YH,yH)),0) min(max(max(YH,yH)),1)])
% ylim([min(min(min(YH,yH)),0) min(max(max(YH,yH)),1)])
xlim([-0.15 0.5])
ylim([-0.15 0.5])
xlabel('Measured impurity ($\%$)', 'FontSize', 10)
ylabel('Impurity ($\%$)', 'FontSize', 10)
title('Mechanistic', 'FontSize', 10)
patch([YH(IH)' YH(IHr)'],                     ...
      [ModelConfidenceInterval95Tri(IH,1)'      ...
       ModelConfidenceInterval95Tri(IHr,2)'],   ...
      'blue', 'FaceAlpha', .125)
set(gca, 'FontSize', 11)
legend("", "Predicted", "95\% CI", ...
        'location', 'northwest',  'FontSize', 8, 'box', 'on')
text(0.2,-0.05, strcat("$R^{2}$: ", num2str(round(R2_Impurity(1),4)*100), "\%"), 'FontSize', 10)
text(0.2,-0.10, strcat("RMSE: ", num2str(round(RMSE_Impurity(1),4))), 'FontSize', 10)

subplot(2, 2, 2); hold all  % DOE (MLR)
% plot(YH, YH, 'o')
plot(-1:0.01:2, -1:0.01:2, 'k-')
plot(YH, MLR_yfit_Impurity, 'x')
patch([YH(IH)' YH(IHr)'], [MLR_yint_Impurity(IH,1)' MLR_yint_Impurity(IHr,2)'], ...
      'blue', 'FaceAlpha', .125) 
xlabel('Measured impurity ($\%$)', 'FontSize', 10)
ylabel('Impurity ($\%$)', 'FontSize', 10)
% xlim([min(min(min([YH MLR_yfit_Impurity MLR_yint_Impurity])),0) ...
%       min(max(max([YH MLR_yfit_Impurity MLR_yint_Impurity])),1)])
% ylim([min(min(min([YH MLR_yfit_Impurity MLR_yint_Impurity])),0) ...
%       min(max(max([YH MLR_yfit_Impurity MLR_yint_Impurity])),1)])
xlim([-0.15 0.5])
ylim([-0.15 0.5])
title('MLR', 'FontSize', 10)
set(gca, 'FontSize', 11)
legend("", "Predicted", "95 \% CI", ...
        'FontSize', 8, 'box', 'on', 'location', 'northwest')
text(0.2,-0.05, strcat("$R^{2}$: ", num2str(round(R2_Impurity(2),4)*100), "\%"), 'FontSize', 10)
text(0.2,-0.10, strcat("RMSE: ", num2str(round(RMSE_Impurity(2),4))), 'FontSize', 10)

subplot(2, 2, 3); hold all % GPR 
% plot(YH, YH, 'o')
plot(-1:0.01:2, -1:0.01:2, 'k-')
plot(YH, GPR_yfit_Impurity, 'x')
patch([YH(IH)' YH(IHr)'], [GPR_yint_Impurity(IH,1)' GPR_yint_Impurity(IHr,2)'], ...
      'blue', 'FaceAlpha', .125) 
% xlim([min(min(min([YH GPR_yfit_Impurity GPR_yint_Impurity])),0) ...
%       min(max(max([YH GPR_yfit_Impurity GPR_yint_Impurity])),1)])
% ylim([min(min(min([YH GPR_yfit_Impurity GPR_yint_Impurity])),0) ...
%       min(max(max([YH GPR_yfit_Impurity GPR_yint_Impurity])),1)])
xlim([-0.15 0.5])
ylim([-0.15 0.5])
xlabel('Measured impurity ($\%$)', 'FontSize', 10)
ylabel('Impurity ($\%$)', 'FontSize', 10)
title('GPR', 'FontSize', 10)
set(gca, 'FontSize', 11)
legend("", "Predicted", "95 \% CI", ...
       'location', 'northwest', 'FontSize', 8, 'box', 'on')
text(0.2,-0.05, strcat("$R^{2}$: ", num2str(round(R2_Impurity(3),4)*100), "\%"), 'FontSize', 10)
text(0.2,-0.10, strcat("RMSE: ", num2str(round(RMSE_Impurity(3),4))), 'FontSize', 10)
   
subplot(2, 2, 4); hold all % ANN
% plot(YH, YH, 'o')
plot(-1:0.01:2, -1:0.01:2, 'k-')
plot(YH, ANN_yfit_Impurity, 'x')
% xlim([min(min(min([YH ANN_yfit_Impurity])),0) ...
%       min(max(max([YH ANN_yfit_Impurity])),1)])
% ylim([min(min(min([YH ANN_yfit_Impurity])),0) ...
%       min(max(max([YH ANN_yfit_Impurity])),1)])
xlim([-0.15 0.5])
ylim([-0.15 0.5])
xlabel('Measured impurity ($\%$)', 'FontSize', 10)
ylabel('Impurity ($\%$)', 'FontSize', 10)
title('ANN', 'FontSize', 10)
set(gca, 'FontSize', 11)
legend("", "Predicted", ...
       'location', 'northwest', 'FontSize', 8, 'box', 'on')
FigureTitle(2) = "FitBoxBehnken_TrainingFit_Tri";
text(0.2,-0.05, strcat("$R^{2}$: ", num2str(round(R2_Impurity(4),4)*100), "\%"), 'FontSize', 10)
text(0.2,-0.10, strcat("RMSE: ", num2str(round(RMSE_Impurity(4),4))), 'FontSize', 10)
sgtitle('BBD - Training')

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
set(gca, 'FontSize', 11)
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

FigureTitle(3) = "FitBoxBehnken_TrainingResiduals_Yield";

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
set(gca, 'FontSize', 8)
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

FigureTitle(4) = "FitBoxBehnken_TrainingResiduals_Impurity";

%% Evaluation of models on test data
testplan.setpoints = [20 11 20 2]; % [T pH Co lambda0 tdose]
testplan.lows      = [5 10 5 1];   % Low values
testplan.highs     = [40 12 50 3]; % High values
Ntest              = 15; % 15

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

reacstruc.data.input = TestData.nom_input;
reacstruc.data.num = [];
for i = 1:Ntest
    reacstruc.data.num{i}(:,1) = TestData.out{i};
end

%% Linear error propagation method 
% Set optim options
% options = optimset('display', 'iter', 'tolfun', 1e-09,   ...
%                    'tolx', 1e-09, 'maxfunevals', 0); % No iteration
%                
% % Compute estimate of Jacobian
% [~,~,residualD,~,~,~,JacobiD] = lsqnonlin(@(x) reacobj_yd(x,reacstruc), ...
%                                         ParameterMinimum, [], [], options);
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
% ConfidenceIntervalP95D = [ParameterMinimum - StandardDeviationD*tcrD; ...
%                           ParameterMinimum + StandardDeviationD*tcrD];
% ConfidenceIntervalP95Tri = [ParameterMinimum - StandardDeviationTri*tcrTri; ...
%                             ParameterMinimum + StandardDeviationTri*tcrTri];
% 
% % Calculate confidence intervals on the model output
% [nD, mD] = size(yD_test);
% [nTri, mTri] = size(yH_test);
% 
% % Compute output covariance
% ModelCovarianceD                 = JacobianD * CovarianceParametersD * JacobianD'; % Calculate model covariance # Maybe some statistical reference here?
% ModelStandardDeviationD          = sqrt(diag(ModelCovarianceD));          % Standard deviation of model outputs, vector with all data points
% ModelStandardDeviationReshapedD  = reshape(ModelStandardDeviationD,nD,mD);  % Reshape to matrix with columns corresponding to each set of data points
% ModelCovarianceTri                 = JacobianTri*CovarianceParametersTri * JacobianTri'; % Calculate model covariance # Maybe some statistical reference here?
% ModelStandardDeviationTri          = sqrt(diag(ModelCovarianceTri));          % Standard deviation of model outputs, vector with all data points
% ModelStandardDeviationReshapedTri  = reshape(ModelStandardDeviationTri,nTri,mTri);  % Reshape to matrix with columns corresponding to each set of data points
% 
% % 95% confidence intervals, calculated only for the measured species (data) Two boundaries, each with a set of m columns 
% ModelConfidenceInterval95D = [yD_test - ModelStandardDeviationReshapedD*tcrD,...  % Lower boundary
%                              yD_test + ModelStandardDeviationReshapedD*tcrD];      % Upper boundary 
% ModelConfidenceInterval95Tri = [yH_test - ModelStandardDeviationReshapedTri*tcrTri,...  % Lower boundary
%                              yH_test + ModelStandardDeviationReshapedTri*tcrTri];      % Upper boundary 


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

% Display statistical measures
disp('Training statistical measures:')
disp([R2_Yield SSE_Yield RMSE_Yield AAD_Yield])
disp([R2_Impurity SSE_Impurity RMSE_Impurity AAD_test_Impurity])

disp('Test statistical measures:')
disp([R2_test_Yield SSE_test_Yield RMSE_test_Yield AAD_test_Yield])
disp([R2_test_Impurity SSE_test_Impurity RMSE_test_Impurity AAD_test_Impurity])

%% Figures of fit to test data (Yield)
[~, IY]  = sort(YTest_Yield);
[~, IYr] = sort(YTest_Yield, 'descend');

figure
subplot(2, 2, 1); hold all % Mechanistic 
% plot(YTest_Yield, YTest_Yield, 'o')
plot(-1:0.01:2, -1:0.01:2, 'k-')
plot(YTest_Yield, yD_test, 'x')
% xlim([min(min(min(YTest_Yield,yD_test)),0) ...
%     max(max(max(YTest_Yield,yD_test+ModelStandardDeviationReshapedD*tcrD)),1)])
% ylim([min(min(min(YTest_Yield,yD_test)),0) ...
%     max(max(max(YTest_Yield,yD_test+ModelStandardDeviationReshapedD*tcrD)),1)])
xlim([-0.2 1.35])
ylim([-0.2 1.35])
xlabel('Measured yield ($\%$)', 'FontSize', 10)
ylabel('Yield ($\%$)', 'FontSize', 10)
title('Mechanistic', 'FontSize', 10)
% patch([YTest_Yield(IY)' YTest_Yield(IYr)'],                     ...
%       [ModelConfidenceInterval95D(IY,1)'      ...
%        ModelConfidenceInterval95D(IYr,2)'],   ...
%       'blue', 'FaceAlpha', .125)
set(gca, 'FontSize', 11)
legend("", "Predicted", ...
        'location', 'northwest', 'FontSize', 8,  'box', 'on')
text(0.6,-0.00, strcat("$R^{2}$: ", num2str(round(R2_test_Yield(1),4)*100), "\%"), 'FontSize', 10)
text(0.6,-0.10, strcat("RMSE: ", num2str(round(RMSE_test_Yield(1),4))), 'FontSize', 10)

subplot(2, 2, 2); hold all  % DOE (MLR)
% plot(YTest_Yield, YTest_Yield, 'o')
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
% plot(YTest_Yield, YTest_Yield, 'o')
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
       'location', 'northwest', 'FontSize', 8,  'box', 'on')
text(0.6,-0.00, strcat("$R^{2}$: ", num2str(round(R2_test_Yield(3),4)*100), "\%"), 'FontSize', 10)
text(0.6,-0.10, strcat("RMSE: ", num2str(round(RMSE_test_Yield(3),4))), 'FontSize', 10)
   
subplot(2, 2, 4); hold all % ANN
% plot(YTest_Yield, YTest_Yield, 'o')
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
FigureTitle(5) = "FitBoxBehnken_TestFit_Yield";
text(0.6,-0.00, strcat("$R^{2}$: ", num2str(round(R2_test_Yield(4),4)*100), "\%"), 'FontSize', 10)
text(0.6,-0.10, strcat("RMSE: ", num2str(round(RMSE_test_Yield(4),4))), 'FontSize', 10)
sgtitle('BBD - Test set')

%% Figures of fit to test data (Impurity)
[~, IH]  = sort(YTest_Impurity);
[~, IHr] = sort(YTest_Impurity, 'descend');

figure
subplot(2, 2, 1); hold all % Mechanistic 
% plot(YTest_Impurity, YTest_Impurity, 'o')
plot(-1:0.01:2, -1:0.01:2, 'k-')
plot(YTest_Impurity, yH_test, 'x')
% xlim([min(min(min(YTest_Impurity, yD_test)),0) ...
%     min(max(max(YTest_Impurity,yH_test)),1)])
% ylim([min(min(min(YTest_Impurity, yD_test)),0) ...
%     min(max(max(YTest_Impurity,yH_test)),1)])
xlim([-0.2 0.4])
ylim([-0.2 0.4])
xlabel('Measured impurity ($\%$)', 'FontSize', 10)
ylabel('Impurity ($\%$)', 'FontSize', 10)
title('Mechanistic', 'FontSize', 10)
% patch([YTest_Impurity(IH)' YTest_Impurity(IHr)'],                     ...
%       [ModelConfidenceInterval95Tri(IH,1)'      ...
%        ModelConfidenceInterval95Tri(IHr,2)'],   ...
%       'blue', 'FaceAlpha', .125)
set(gca, 'FontSize', 11)
legend("", "Predicted", ...
        'location', 'northwest', 'FontSize', 8, 'box', 'on')
text(0.15, -0.135, strcat("$R^{2}$: ", num2str(round(R2_test_Impurity(1),4)*100), "\%"), 'FontSize', 10)
text(0.15, -0.17, strcat("RMSE: ", num2str(round(RMSE_test_Impurity(1),4))), 'FontSize', 10)

subplot(2, 2, 2); hold all  % DOE (MLR)
% plot(YTest_Impurity, YTest_Impurity, 'o')
plot(-1:0.01:2, -1:0.01:2, 'k-')
plot(YTest_Impurity, MLR_yfittest_Impurity, 'x')
patch([YTest_Impurity(IH)' YTest_Impurity(IHr)'], ...
    [MLR_yinttest_Impurity(IH,1)' MLR_yinttest_Impurity(IHr,2)'], ...
      'blue', 'FaceAlpha', .125) 
xlabel('Measured impurity ($\%$)', 'FontSize', 10)
ylabel('Impurity ($\%$)', 'FontSize', 10)
% xlim([min(min(min([YTest_Impurity MLR_yfittest_Impurity MLR_yinttest_Impurity])),0) ...
%       min(max(max([YTest_Impurity MLR_yfittest_Impurity MLR_yinttest_Impurity])),1)])
% ylim([min(min(min([YTest_Impurity MLR_yfittest_Impurity MLR_yinttest_Impurity])),0) ...
%       min(max(max([YTest_Impurity MLR_yfittest_Impurity MLR_yinttest_Impurity])),1)])
xlim([-0.2 0.4])
ylim([-0.2 0.4])
title('MLR', 'FontSize', 10)
set(gca, 'FontSize', 11)
legend("", "Predicted", "95 \% PI", ...
       'FontSize', 8, 'location', 'northwest', 'box', 'on')
text(0.15,-0.135, strcat("$R^{2}$: ", num2str(round(R2_test_Impurity(2),4)*100), "\%"), 'FontSize', 10)
text(0.15,-0.17, strcat("RMSE: ", num2str(round(RMSE_test_Impurity(2),4))), 'FontSize', 10)

subplot(2, 2, 3); hold all % GPR 
% plot(YTest_Impurity, YTest_Impurity, 'o')
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
% plot(YTest_Impurity, YTest_Impurity, 'o')
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
FigureTitle(6) = "FitBoxBehnken_TestFit_Tri";
text(0.15,-0.135, strcat("$R^{2}$: ", num2str(round(R2_test_Impurity(4),4)*100), "\%"), 'FontSize', 10)
text(0.15,-0.17, strcat("RMSE: ", num2str(round(RMSE_test_Impurity(4),4))), 'FontSize', 10)
sgtitle('BBD - Test set')

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

N_MM = 1e4;
pars_MM = {'T', 'pH', 'Co', 'lambda0', 'tdose'};
lbs_MM = [lbs 10];
ubs_MM = [ubs 60];
InputSpace = {'ParNames', pars, 'LowerBounds', lbs, 'UpperBounds', ubs};
InputSpace_MM = {'ParNames', pars_MM, 'LowerBounds', lbs_MM, ...
                                      'UpperBounds', ubs_MM};

% Compute Sobol indices
[Si_MLR_Yield, STi_MLR_Yield]       = easyGSA(f_MLR_Yield, N, InputSpace{:});
[Si_MLR_Impurity, STi_MLR_Impurity] = easyGSA(f_MLR_Impurity, N, InputSpace{:});
[Si_GPR_Yield, STi_GPR_Yield]       = easyGSA(f_GPR_Yield, N, InputSpace{:});
[Si_GPR_Impurity, STi_GPR_Impurity] = easyGSA(f_GPR_Impurity, N, InputSpace{:});
[Si_ANN_Yield, STi_ANN_Yield]       = easyGSA(f_ANN_Yield, N, InputSpace{:});
[Si_ANN_Impurity, STi_ANN_Impurity] = easyGSA(f_ANN_Impurity, N, InputSpace{:});
% [Si_MM_Yield, STi_MM_Yield]         = easyGSA(f_MM_Yield, N_MM, InputSpace_MM{:});
% [Si_MM_Impurity, STi_MM_Impurity]   = easyGSA(f_MM_Impurity, N_MM, InputSpace_MM{:});

disp([Si_MLR_Yield, STi_MLR_Yield])
disp([Si_MLR_Impurity STi_MLR_Impurity])
disp([Si_GPR_Yield, STi_GPR_Yield])
disp([Si_GPR_Impurity, STi_GPR_Impurity])
disp([Si_ANN_Yield, STi_ANN_Yield])
disp([Si_ANN_Impurity, STi_ANN_Impurity])
% disp([Si_MM_Yield, STi_MM_Yield])
% disp([Si_MM_Impurity, STi_MM_Impurity])

% No noise
%     0.0212    0.0249
%     0.0699    0.1153
%     0.0679    0.0934
%     0.7869    0.8155
% 
%     0.0113    0.0268
%     0.2193    0.5067
%     0.0225    0.0539
%     0.4190    0.7377
% 
%     0.0144    0.0157
%     0.0339    0.0995
%     0.0377    0.0654
%     0.8272    0.9041
% 
%     0.0077    0.0100
%     0.2154    0.5665
%     0.0146    0.0247
%     0.3916    0.7485
% 
%     0.0186    0.0220
%     0.0573    0.1307
%     0.0469    0.0827
%     0.7860    0.8631
% 
%     0.0175    0.0420
%     0.2138    0.5577
%     0.0183    0.0712
%     0.3591    0.7356
% 
%     0.0024    0.0145
%     0.0322    0.1099
%     0.0159    0.0674
%     0.8428    0.9201
%    -0.0062    0.0000
% 
%     0.0162    0.0274
%     0.2356    0.5526
%     0.0476    0.1253
%     0.3429    0.7030
%    -0.0009    0.0000
% 
%
% With noise:
% 0.0016    0.0046
%     0.0849    0.1421
%     0.0410    0.0825
%     0.8075    0.8321
% 
%     0.0041    0.0088
%     0.3521    0.6043
%     0.0251    0.0529
%     0.3384    0.6113
% 
%     0.0044    0.0062
%     0.0544    0.1326
%     0.0263    0.0529
%     0.8257    0.8964
% 
%     0.0057    0.0039
%     0.3226    0.6325
%     0.0216    0.0349
%     0.3232    0.6396
% 
%     0.0007    0.0027
%     0.0849    0.1036
%     0.0339    0.0458
%     0.8490    0.8795
% 
%     0.0073    0.0113
%     0.3131    0.6034
%     0.0304    0.1100
%     0.3203    0.6228
% 
%    -0.0054    0.0031
%     0.0339    0.1177
%     0.0089    0.0615
%     0.8518    0.9314
%    -0.0067    0.0000
% 
%     0.0086    0.0100
%     0.2145    0.5401
%     0.0517    0.1317
%     0.3635    0.7312
%     0.0010    0.0000

%% Propagation of uncertainty with Monte Carlo simulation
% if SimulationIndex == 2
    Beta_Yield = table2array(MLR_Yield.Coefficients(:,1));
    Beta_Tri   = table2array(MLR_Impurity.Coefficients(:,1));

    % Compute design matrix
    XMat = lm_mat(plan(:,1:4));

    % Compute sample standard deviation
    sd_Yield = (YD-XMat*Beta_Yield)'*(YD-XMat*Beta_Yield) / (size(XMat,1)-length(Beta_Yield));
    sd_Tri = (YH-XMat*Beta_Tri)'*(YH-XMat*Beta_Tri) / (size(XMat,1)-length(Beta_Tri));

    % Compute parameter covariances
    Cov_Yield = sd_Yield^2*inv(XMat'*XMat);
    Cov_Tri = sd_Tri^2*inv(XMat'*XMat);

    % Sample parameter values from multivariate normal distribution
    Beta_Yield_Sample = mvnrnd(Beta_Yield, Cov_Yield, 1000);
    Beta_Impurity_Sample = mvnrnd(Beta_Tri, Cov_Tri, 1000);
    
    % Setpoint
    xset = [20 11 20 2];

    %
    yD_beta = zeros(size(Beta_Yield_Sample,1), 1);
    yH_beta = zeros(size(Beta_Impurity_Sample,1), 1);
    for k = 1:size(Beta_Yield_Sample,1)
        yD_beta(k,1) = predict_lm(xset,Beta_Yield_Sample(k,:));
        yH_beta(k,1) = predict_lm(xset,Beta_Impurity_Sample(k,:));
    end
    
%     figure
%     [~,ax] = plotmatrix(Beta_Yield_Sample);
%     ax(1,1).YLabel.String = "$k_{A,ref}$";
%     ax(2,1).YLabel.String = "$s_{ref,3}$";
%     ax(3,1).YLabel.String = "$s_{d,ref}$";
%     ax(4,1).YLabel.String = "pk$_{a,1}$";
%     ax(5,1).YLabel.String = "pk$_{a,3}$";
%     ax(6,1).YLabel.String = "E$_{A,1}$";
%     ax(7,1).YLabel.String = "E$_{A,3}$";
%     ax(8,1).YLabel.String = "E$_{A,d}$";
%     
%     ax(8,1).XLabel.String = "$k_{A,ref}$";
%     ax(8,2).XLabel.String = "$s_{ref,3}$";
%     ax(8,3).XLabel.String = "$s_{d,ref}$";
%     ax(8,4).XLabel.String = "pk$_{a,1}$";
%     ax(8,5).XLabel.String = "pk$_{a,3}$";
%     ax(8,6).XLabel.String = "E$_{A,1}$";
%     ax(8,7).XLabel.String = "E$_{A,3}$";
%     ax(8,8).XLabel.String = "E$_{A,d}$";
%     for i = 1:8
%         ax(8,i).XAxis.FontSize = 7;
%         ax(i,1).YAxis.FontSize = 7;
%     end
%     FigureTitle(7) = "FitBoxBehnken_UA_PlotMatrix";
    
    figure; hold all
    histogram(yD_beta, 'normalization', 'probability')
    ylabel('Frequency (probability)')
    xlabel('Yield, $y_{D}$')
    xline(mean(yD_beta), 'r-', 'LineWidth', 2)
    xline(0.8150, 'k-', 'LineWidth', 2)
    FigureTitle(7) = "FitBoxBehnken_UA_MLR_Yield";

    figure
    histogram(yH_beta, 'normalization', 'probability')
    ylabel('Frequency (probability)')
    xlabel('Impurity, $y_{H}$')
    xline(mean(yH_beta), 'r-', 'LineWidth', 2)
    xline(0.0114, 'k-', 'LineWidth', 2)
    FigureTitle(8) = "FitBoxBehnken_UA_MLR_Tri";
% end




%% Residuals of training
% Residuals_YD = MLR_yfit_Yield-YD;
% Residuals_YH = MLR_yfit_Impurity-YH;
% 
% figure; hold all
% histogram(Residuals_YD,'BinWidth', 0.025)
% plot(min(Residuals_YD):0.01:max(Residuals_YD), ...
%     normpdf(min(Residuals_YD):0.01:max(Residuals_YD), 0, std(MLR_yfit_Yield-YD)))
% 
% figure
% qqplot(Residuals_YD)
% 
% figure
% plot(MLR_yfit_Yield, Residuals_YD, 'o')
% 
% % Compute residuals for test data
% Residuals_YD_test = MLR_yfittest_Yield - YTest_Yield;
% Residuals_YH_test = MLR_yfittest_Impurity - YTest_Impurity;
% 
% figure; hold all
% histogram(Residuals_YD_test,'BinWidth', 0.025)
% plot(min(Residuals_YD_test):0.01:max(Residuals_YD_test), ...
%     normpdf(min(Residuals_YD_test):0.01:max(Residuals_YD_test), 0, ...
%     std(Residuals_YD_test)))
% 
% figure
% qqplot(Residuals_YD_test)
% 
% figure
% plot(MLR_yfittest_Yield, Residuals_YD_test, 'o')

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
MLR_sol     = zeros(6, 4);  %
MLR_sol_fmincon = zeros(6, 4);
MLR_opt     = zeros(6, 4);  %
MLR_opt_fmincon = zeros(6, 4);
IdxOptGuess_MLR = zeros(6, 1);   % 
costGuess_MLR   = zeros(6, 1);
for i = 1:length(price)
    % Compute best initial guess with default constraints
    [costGuess_MLR(2*i-1), IdxOptGuess_MLR(2*i-1)] = ...
        min((1 + price(i)*OptPlan(k{1},4)) ./ YD_MLR(k{1}));
    
    % Compute best initial guess with impurity below the constraint of 0.01
    [costGuess_MLR(2*i), IdxOptGuess_MLR(2*i)] = ...
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
    optim_MLR.gub = [1 0.01];
    optim_MLR.gub = [1 0.005];
    optim_MLR.price  = price(i);
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
    GPR_opt_fmincon(2*i-1,:) = GPR_sol_fmincon(2*i-1,:) .* (highs -lows) + lows;

    % Compute minimizer with impurity below the constraint of 0.01
    optim_GPR.z0   = OptPlan(k{2}(IdxOptGuess_GPR(2*i)),:);  % Set initial guess
%     optim_GPR.z0   = TrueOpt(2*i,:);
    optim_GPR.glb  = [0 0];                           % Lower constraints
    optim_GPR.gub = [1 0.01];
    optim_GPR.gub = [1 0.005];
    optim_GPR.price  = price(i);
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
    optim_ANN.gub = [1 0.01];
    optim_ANN.gub = [1 0.005];
    optim_ANN.price  = price(i);
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

%% Run estimated optima for fmincon
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
for i = 1:size(MLR_opt_fmincon,1)
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

%% Run estimated optima for fminsearchcon
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
    reacstruc.process.tdose   = 30;             % min

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
%     0.1425   -0.1174
%    -0.2030   -0.4508
%     0.2550    0.1349
%     0.8237    0.5985
% 
%     0.1069   -0.0771
%    -0.1308   -0.4580
%     0.1925    0.1030
%     0.8558    0.5826
% 
%     0.1290   -0.1211
%    -0.1805   -0.4373
%     0.2155    0.1573
%     0.8304    0.5398
% 
% Squared SRCs:
%     0.0203    0.0138
%     0.0412    0.2032
%     0.0650    0.0182
%     0.6785    0.3582
% 
%     0.0114    0.0059
%     0.0171    0.2098
%     0.0371    0.0106
%     0.7324    0.3394
% 
%     0.0166    0.0147
%     0.0326    0.1913
%     0.0465    0.0247
%     0.6896    0.2914
% 
% Sum of squared SRCs:
%     0.8050    0.5934
% 
%     0.7981    0.5658
% 
%     0.7852    0.5221

% Noise
% SRCs:
%     0.0404   -0.0544
%    -0.1242   -0.4246
%     0.1118    0.1864
%     0.8509    0.5597
% 
%     0.0540   -0.0763
%    -0.2654   -0.5608
%     0.2085    0.1475
%     0.8346    0.5492
% 
%     0.0352   -0.0539
%    -0.2095   -0.5480
%     0.1641    0.1272
%     0.8542    0.5440
% 
%     0.0431   -0.0772
%    -0.2925   -0.5211
%     0.1951    0.1879
%     0.9062    0.5362
% 
% Squared SRCs:
%     0.0029    0.0058
%     0.0704    0.3145
%     0.0435    0.0218
%     0.6966    0.3016
% 
%     0.0012    0.0029
%     0.0439    0.3003
%     0.0269    0.0162
%     0.7297    0.2960
% 
%     0.0019    0.0060
%     0.0855    0.2715
%     0.0381    0.0353
%     0.8213    0.2875
% 
% Sum of squared SRCs:
%     0.8134    0.6437
% 
%     0.8018    0.6153
% 
%     0.9467    0.6003


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
sgtitle('BBD - MLR predictions')
FigureTitle(9) = "FitBoxBehnken_OFAT_MLR";


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
sgtitle('BBD - GPR predictions')
legend("Ground truth", "Predictions", "95 \% CI", 'location', 'southeast')
% text(1.25, 1.15, strcat("$R^{2}$: ", num2str(round(OFAT_lambda0_stats_GPR.R2,4))), 'FontSize', 10)
% text(1.25, 1.05, strcat("RMSE: ", num2str(round(OFAT_lambda0_stats_GPR.RMSE,4))), 'FontSize', 10)
FigureTitle(10) = "FitBoxBehnken_OFAT_GPR";

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
sgtitle('BBD - ANN predictions')
FigureTitle(11) = "FitBoxBehnken_OFAT_ANN";

%% Statistical measures results
% No noise
% Test statistical measures:
%     1.0000    0.0000    0.0000    0.0000
%     0.9522    0.0760    0.0712    0.0612
%     0.9800    0.0338    0.0474    0.0386
%     0.9532    0.0518    0.0588    0.0459
% 
%     1.0000    0.0000    0.0000    0.0000
%     0.9005    0.0046    0.0174    0.0154
%     0.9497    0.0038    0.0160    0.0114
%     0.9723    0.0025    0.0130    0.0084

% Noise
% Test statistical measures:
%     0.9959    0.0070    0.0216    0.0148
%     0.9561    0.1080    0.0849    0.0704
%     0.9548    0.0465    0.0557    0.0483
%     0.8917    0.1275    0.0922    0.0826
% 
%     0.9889    0.0003    0.0046    0.0029
%     0.8717    0.0041    0.0165    0.0134
%     0.9172    0.0044    0.0171    0.0096
%     0.9550    0.0015    0.0099    0.0069

%% Optimization results
% No noise 
% MLR_opt =
% 
%    38.0360   11.1111   36.3068    2.6771
%    39.4517   11.4684   37.6504    2.5192
%    35.6789   10.9319   37.3592    2.1441
%    38.6093   10.9672   37.6342    2.1222
%    34.3941   10.8342   37.9329    1.8536
%    34.3941   10.8342   37.9329    1.8536
% 
% 
% GPR_opt =
% 
%    39.6478   10.9223   41.7331    2.1947
%    40.0000   11.3130   25.3575    2.8645
%    40.0000   10.7012   37.2140    2.0653
%    39.9975   10.7067   36.7869    1.9976
%    40.0000   10.6610   37.4224    1.9981
%    40.0000   10.6614   37.3980    1.9925
% 
% 
% ANN_opt =
% 
%    39.6150   11.0088   40.7808    2.3943
%    39.9458   11.5268   15.9110    2.9616
%    40.0000   10.2544   50.0000    1.9193
%    40.0000   10.2967   50.0000    1.9044
%    40.0000   10.0142   50.0000    1.6996
%    40.0000   10.0142   50.0000    1.6996
% 
% 
% MLR_opt_fmincon =
% 
%    38.0358   11.1111   36.3068    2.6771
%    40.0000   11.4298   37.2455    2.5000
%    35.6789   10.9319   37.3592    2.1441
%    37.3798   10.9776   37.6278    2.1081
%    34.3940   10.8342   37.9328    1.8536
%    34.3941   10.8342   37.9328    1.8536
% 
% 
% GPR_opt_fmincon =
% 
%    37.7633   10.7832   37.3506    2.2313
%    40.0000   11.3130   25.3576    2.8645
%    40.0000   10.7012   37.2140    2.0653
%    40.0000   10.7084   36.8156    1.9977
%    40.0000   10.6610   37.4224    1.9981
%    40.0000   10.6614   37.3980    1.9925
% 
% 
% ANN_opt_fmincon =
% 
%    39.3339   10.9847   37.5395    2.3951
%    37.7446   11.4560   18.0780    2.9492
%    40.0000   10.2544   50.0000    1.9193
%    40.0000   10.2966   50.0000    1.9043
%    40.0000   10.0142   50.0000    1.6996
%    40.0000   10.0142   50.0000    1.6996


% Noise
% MLR_opt =
% 
%    40.0000   11.4225   50.0000    2.6871
%    39.9988   11.6856   49.9902    2.6035
%    40.0000   11.2153   50.0000    2.1659
%    40.0000   11.2153   50.0000    2.1659
%    39.9997   11.0430   48.3783    1.8915
%    40.0000   11.0463   48.4642    1.8945
% 
% 
% GPR_opt =
% 
%    35.3421   11.0782   46.8879    2.4102
%    39.3993   11.2975   28.3233    2.6337
%    36.9980   10.6493   38.6836    2.1256
%    39.9943   11.0554   44.6309    2.0888
%    37.3534   10.5127   37.9652    2.0101
%    39.9991   10.9363   45.0205    2.0005
% 
% 
% ANN_opt =
% 
%    40.0000   10.0000   50.0000    3.0000
%    39.3771   11.1589    8.0971    2.9999
%    40.0000   10.0000   50.0000    1.7147
%    39.7044   10.0753   49.9459    1.4843
%    40.0000   10.0000   50.0000    1.3761
%    40.0000   10.0000   50.0000    1.3761
% 
% 
% MLR_opt_fmincon =
% 
%    40.0000   11.4225   49.9999    2.6871
%    40.0000   11.6826   50.0000    2.6004
%    40.0000   11.2153   50.0000    2.1659
%    39.9999   11.2153   49.9992    2.1659
%    40.0000   11.0463   48.4641    1.8945
%    40.0000   11.0463   48.4641    1.8945
% 
% 
% GPR_opt_fmincon =
% 
%    29.1851   10.8459   35.2618    2.4240
%    39.9999   11.3217   30.4139    2.6008
%    36.9980   10.6493   38.6836    2.1256
%    40.0000   11.0660   44.6317    2.0960
%    37.3533   10.5127   37.9652    2.0101
%    40.0000   10.9362   45.0465    2.0004
% 
% 
% ANN_opt_fmincon =
% 
%    40.0000   10.0000   50.0000    3.0000
%    40.0000   11.0429    5.0000    3.0000
%    40.0000   10.0000   50.0000    1.7147
%    40.0000   10.0000   50.0000    1.4402
%    40.0000   10.0000   50.0000    1.3761
%    40.0000   10.0000   50.0000    1.3761


%% More optimization results
% No noise
% True yield with estimated opt vs predicted opt yield (MLR) fmincon
%     0.9516    0.9984
%     0.9819    0.9719
%     0.9756    0.9361
%     0.9685    0.9267
%     0.8004    0.8497
%     0.8004    0.8497
% 
% True impurity with estimated opt vs predicted opt impurity (MLR) fmincon
%     0.0484    0.0558
%     0.0179    0.0100
%     0.0169    0.0179
%     0.0125    0.0100
%     0.0065    0.0007
%     0.0065    0.0007
% 
% True yield with estimated opt vs predicted opt yield (GPR) fmincon
%     0.9666    1.0000
%     0.9724    0.9922
%     0.9630    0.9887
%     0.9215    0.9667
%     0.9234    0.9673
%     0.9193    0.9651
% 
% True impurity with estimated opt vs predicted opt impurity (GPR) fmincon
%     0.0333    0.0426
%     0.0276    0.0000
%     0.0141    0.0193
%     0.0100    0.0100
%     0.0107    0.0108
%     0.0105    0.0100
% 
% True yield with estimated opt vs predicted opt yield (ANN) fmincon
%     0.9628    1.0000
%     0.9828    1.0000
%     0.8675    0.9524
%     0.8569    0.9482
%     0.6861    0.8806
%     0.6861    0.8806
% 
% True impurity with estimated opt vs predicted opt impurity (ANN) fmincon
%     0.0372    0.0183
%     0.0166    0.0066
%     0.0147    0.0119
%     0.0131    0.0100
%     0.0110    0.0073
%     0.0110    0.0073
% 
% True cost with estimated opt vs predicted opt cost (MLR) fmincon
%     1.0508    1.0016
%     1.0184    1.0290
%     2.2052    2.2983
%     2.2014    2.3009
%     5.8807    5.5395
%     5.8807    5.5395
% 
% True cost with estimated opt vs predicted opt cost (GPR) fmincon
%     1.0346    1.0000
%     1.0283    1.0079
%     2.1902    2.1333
%     2.2495    2.1443
%     5.4110    5.1654
%     5.4224    5.1655
% 
% True cost with estimated opt vs predicted opt cost (ANN) fmincon
%     1.0386    1.0000
%     1.0175    1.0000
%     2.3409    2.1323
%     2.3605    2.1332
%     6.4114    4.9959
%     6.4114    4.9959
% 
% True yield with estimated opt vs predicted opt yield (MLR) fminsearchcon
%     0.9516    0.9984
%     0.9823    0.9706
%     0.9756    0.9361
%     0.9748    0.9297
%     0.8004    0.8497
%     0.8004    0.8497
% 
% True impurity with estimated opt vs predicted opt impurity (MLR) fminsearchcon
%     0.0484    0.0558
%     0.0175    0.0100
%     0.0169    0.0179
%     0.0133    0.0100
%     0.0065    0.0007
%     0.0065    0.0007
% 
% True yield with estimated opt vs predicted opt yield (GPR) fminsearchcon
%     0.9768    1.0000
%     0.9724    0.9922
%     0.9630    0.9887
%     0.9214    0.9667
%     0.9234    0.9673
%     0.9193    0.9651
% 
% True impurity with estimated opt vs predicted opt impurity (GPR) fminsearchcon
%     0.0229    0.0295
%     0.0276    0.0000
%     0.0141    0.0193
%     0.0100    0.0100
%     0.0107    0.0108
%     0.0105    0.0100
% 
% True yield with estimated opt vs predicted opt yield (ANN) fminsearchcon
%     0.9622    1.0000
%     0.9853    1.0000
%     0.8675    0.9524
%     0.8569    0.9482
%     0.6861    0.8806
%     0.6861    0.8806
% 
% True impurity with estimated opt vs predicted opt impurity (ANN) fminsearchcon
%     0.0378    0.0181
%     0.0121    0.0007
%     0.0147    0.0119
%     0.0131    0.0100
%     0.0110    0.0073
%     0.0110    0.0073
% 
% True cost with estimated opt vs predicted opt cost (MLR) fminsearchcon
%     1.0508    1.0016
%     1.0180    1.0302
%     2.2052    2.2983
%     2.1951    2.3015
%     5.8807    5.5395
%     5.8807    5.5395
% 
% True cost with estimated opt vs predicted opt cost (GPR) fminsearchcon
%     1.0237    1.0000
%     1.0283    1.0079
%     2.1902    2.1333
%     2.2497    2.1443
%     5.4110    5.1654
%     5.4224    5.1655
% 
% True cost with estimated opt vs predicted opt cost (ANN) fminsearchcon
%     1.0393    1.0000
%     1.0149    1.0000
%     2.3409    2.1323
%     2.3604    2.1332
%     6.4114    4.9959
%     6.4114    4.9959


% Noise
% True yield with estimated opt vs predicted opt yield (MLR) fmincon
%     0.9670    0.9912
%     0.9835    0.9802
%     0.9830    0.9310
%     0.9830    0.9310
%     0.8398    0.8516
%     0.8398    0.8516
% 
% True impurity with estimated opt vs predicted opt impurity (MLR) fmincon
%     0.0330    0.0374
%     0.0164    0.0100
%     0.0128    0.0094
%     0.0128    0.0094
%     0.0050    0.0000
%     0.0050    0.0000
% 
% True yield with estimated opt vs predicted opt yield (GPR) fmincon
%     0.9404    1.0000
%     0.9775    0.9894
%     0.9720    0.9719
%     0.9706    0.9491
%     0.9286    0.9364
%     0.9235    0.9201
% 
% True impurity with estimated opt vs predicted opt impurity (GPR) fmincon
%     0.0596    0.0711
%     0.0224    0.0100
%     0.0249    0.0524
%     0.0105    0.0100
%     0.0153    0.0493
%     0.0080    0.0100
% 
% True yield with estimated opt vs predicted opt yield (ANN) fmincon
%     0.5335    0.9355
%     0.9768    0.8969
%     0.6972    0.8243
%     0.4982    0.7446
%     0.4558    0.7210
%     0.4558    0.7210
% 
% True impurity with estimated opt vs predicted opt impurity (ANN) fmincon
%     0.4665    0.3476
%     0.0120    0.0100
%     0.0117    0.0342
%     0.0056    0.0100
%     0.0047    0.0059
%     0.0047    0.0059
% 
% True cost with estimated opt vs predicted opt cost (MLR) fmincon
%     1.0341    1.0089
%     1.0168    1.0202
%     2.2007    2.3236
%     2.2007    2.3236
%     5.7025    5.6233
%     5.7025    5.6233
% 
% True cost with estimated opt vs predicted opt cost (GPR) fmincon
%     1.0634    1.0000
%     1.0230    1.0107
%     2.2032    2.2036
%     2.1901    2.2397
%     5.4065    5.3615
%     5.4150    5.4351
% 
% True cost with estimated opt vs predicted opt cost (ANN) fmincon
%     1.8742    1.0690
%     1.0238    1.1149
%     2.7551    2.3302
%     3.5598    2.3816
%     8.2321    5.2042
%     8.2321    5.2042
% 
% True yield with estimated opt vs predicted opt yield (MLR) fminsearchcon
%     0.9670    0.9912
%     0.9835    0.9802
%     0.9830    0.9310
%     0.9830    0.9310
%     0.8375    0.8506
%     0.8398    0.8516
% 
% True impurity with estimated opt vs predicted opt impurity (MLR) fminsearchcon
%     0.0330    0.0374
%     0.0163    0.0100
%     0.0128    0.0094
%     0.0128    0.0094
%     0.0049    0.0000
%     0.0050    0.0000
% 
% True yield with estimated opt vs predicted opt yield (GPR) fminsearchcon
%     0.9581    1.0000
%     0.9764    0.9889
%     0.9720    0.9719
%     0.9683    0.9473
%     0.9286    0.9364
%     0.9236    0.9201
% 
% True impurity with estimated opt vs predicted opt impurity (GPR) fminsearchcon
%     0.0419    0.0445
%     0.0236    0.0100
%     0.0249    0.0524
%     0.0103    0.0100
%     0.0153    0.0493
%     0.0080    0.0100
% 
% True yield with estimated opt vs predicted opt yield (ANN) fminsearchcon
%     0.5335    0.9355
%     0.9828    0.8959
%     0.6972    0.8243
%     0.5296    0.7505
%     0.4558    0.7210
%     0.4558    0.7210
% 
% True impurity with estimated opt vs predicted opt impurity (ANN) fminsearchcon
%     0.4665    0.3476
%     0.0150    0.0100
%     0.0117    0.0342
%     0.0057    0.0100
%     0.0047    0.0059
%     0.0047    0.0059
% 
% True cost with estimated opt vs predicted opt cost (MLR) fminsearchcon
%     1.0341    1.0089
%     1.0168    1.0202
%     2.2007    2.3236
%     2.2007    2.3236
%     5.7107    5.6233
%     5.7025    5.6233
% 
% True cost with estimated opt vs predicted opt cost (GPR) fminsearchcon
%     1.0437    1.0000
%     1.0242    1.0112
%     2.2032    2.2036
%     2.1913    2.2398
%     5.4065    5.3615
%     5.4149    5.4351
% 
% True cost with estimated opt vs predicted opt cost (ANN) fminsearchcon
%     1.8742    1.0690
%     1.0175    1.1162
%     2.7551    2.3302
%     3.3933    2.3945
%     8.2321    5.2042
%     8.2321    5.2042

%% Smarter optimization notes (fminsearchcon)
% No noise
% Optimum predicted (MLR | GPR | ANN) case1
%    39.4517   40.0000   39.9458
%    11.4684   11.3130   11.5268
%    37.6504   25.3589   15.9110
%     2.5192    2.8645    2.9616
%     0.9706    0.9922    1.0000
%     0.0100    0.0000    0.0007
%     1.0302    1.0079    1.0000
% 
% Optimum predicted (MLR | GPR | ANN) case2
%    38.6093   39.9999   40.0000
%    10.9672   10.7100   10.2967
%    37.6342   36.8339   50.0000
%     2.1222    1.9978    1.9044
%     0.9297    0.9667    0.9482
%     0.0100    0.0100    0.0100
%     2.3015    2.1442    2.1332
% 
% Optimum predicted (MLR | GPR | ANN) case3
%    34.3941   40.0000   40.0000
%    10.8342   10.6615   10.0142
%    37.9329   37.4018   50.0000
%     1.8536    1.9925    1.6996
%     0.8497    0.9651    0.8806
%     0.0007    0.0100    0.0073
%     5.5395    5.1655    4.9959
% 
% Difference between predicted and actual minima (MLR | GPR | ANN) case1
%    -0.5466    0.0017   -0.0525
%    -0.5236   -0.6790   -0.4652
%    -6.3560  -18.6475  -28.0954
%    -0.3132    0.0321    0.1292
% 
%     6.4085   18.6599   28.0996
% 
% Difference between predicted and actual minima (MLR | GPR | ANN) case2
%    -1.3724    0.0182    0.0183
%    -0.1070   -0.3642   -0.7775
%   -12.0230  -12.8233    0.3428
%     0.0440   -0.0804   -0.1738
% 
%    12.1016   12.8287    0.8676
% 
% Difference between predicted and actual minima (MLR | GPR | ANN) case3
%    -5.5516    0.0543    0.0543
%    -0.1339   -0.3066   -0.9539
%   -11.1774  -11.7085    0.8897
%    -0.1983   -0.0594   -0.3523
% 
%    12.4825   11.7128    1.3523
% 
% Difference between optima for yD, yH and cost (MLR | GPR | ANN) case1
%     0.9823    0.0175    1.0180    0.0117    0.0075   -0.0122
%     0.9724    0.0276    1.0283   -0.0198    0.0276    0.0205
%     0.9853    0.0121    1.0149   -0.0147    0.0114    0.0149
% 
% Difference between optima for yD, yH and cost (MLR | GPR | ANN) case2
%     0.9748    0.0133    2.1951    0.0450    0.0033   -0.1064
%     0.9215    0.0100    2.2495   -0.0452   -0.0000    0.1052
%     0.8569    0.0131    2.3604   -0.0913    0.0031    0.2272
% 
% Difference between optima for yD, yH and cost (MLR | GPR | ANN) case3
%     0.8004    0.0065    5.8807   -0.0493    0.0059    0.3412
%     0.9193    0.0105    5.4224   -0.0457    0.0005    0.2569
%     0.6861    0.0110    6.4114   -0.1944    0.0037    1.4155


% Noise
% Optimum predicted (MLR | GPR | ANN) case1
%    39.9988   39.9027   39.3771
%    11.6856   11.3162   11.1589
%    49.9902   29.8567    8.0971
%     2.6035    2.6122    2.9999
%     0.9802    0.9893    0.8959
%     0.0100    0.0100    0.0100
%     1.0202    1.0108    1.1162
% 
% Optimum predicted (MLR | GPR | ANN) case2
%    40.0000   39.9065   39.7044
%    11.2153   11.0796   10.0753
%    50.0000   46.1601   49.9459
%     2.1659    2.0952    1.4843
%     0.9310    0.9483    0.7505
%     0.0094    0.0100    0.0100
%     2.3236    2.2411    2.3945
% 
% Optimum predicted (MLR | GPR | ANN) case3
%    40.0000   40.0000   40.0000
%    11.0463   10.9357   10.0000
%    48.4642   44.9249   50.0000
%     1.8945    2.0004    1.3761
%     0.8516    0.9199    0.7210
%     0.0000    0.0100    0.0059
%     5.6233    5.4362    5.2042
% 
% Difference between predicted and actual minima (MLR | GPR | ANN) case1
%     0.0005   -0.0956   -0.6212
%    -0.3064   -0.6758   -0.8331
%     5.9838  -14.1497  -35.9093
%    -0.2289   -0.2202    0.1675
% 
%     5.9961   14.1679   35.9247
% 
% Difference between predicted and actual minima (MLR | GPR | ANN) case2
%     0.0183   -0.0752   -0.2773
%     0.1411    0.0054   -0.9989
%     0.3428   -3.4971    0.2887
%     0.0877    0.0170   -0.5939
% 
%     0.3814    3.4979    1.2291
% 
% Difference between predicted and actual minima (MLR | GPR | ANN) case3
%     0.0543    0.0543    0.0543
%     0.0782   -0.0324   -0.9681
%    -0.6461   -4.1854    0.8897
%    -0.1574   -0.0515   -0.6758
% 
%     0.6718    4.1862    1.4793
% 
% Difference between of optima and actual values of yD, yH and cost (MLR | GPR | ANN) case1
%     0.9835    0.0163    1.0168    0.0033    0.0063   -0.0035
%     0.9772    0.0228    1.0233   -0.0121    0.0128    0.0125
%     0.9828    0.0150    1.0175    0.0868    0.0050   -0.0986
% 
% Difference between of optima and actual values for yD, yH and cost (MLR | GPR | ANN) case2
%     0.9830    0.0128    2.2007    0.0520    0.0034   -0.1229
%     0.9708    0.0105    2.1893    0.0224    0.0005   -0.0518
%     0.5296    0.0057    3.3933   -0.2209   -0.0043    0.9987
% 
% Difference between optima and actual values for yD, yH and cost (MLR | GPR | ANN) case3
%     0.8398    0.0050    5.7025   -0.0118    0.0050    0.0792
%     0.9234    0.0080    5.4155    0.0035   -0.0020   -0.0207
%     0.4558    0.0047    8.2321   -0.2652   -0.0012    3.0280


%% Save figures
if bool_SaveFigures
    for i = 1:length(FigureTitle) 
        saveas(figure(i), ...
            strcat(pwd, "\Images\FitBoxBehnkenNoise\", ...
                        FigureTitle(i), "_NoiseLevel_", num2str(SimulationIndex)), 'png')
    end
end


