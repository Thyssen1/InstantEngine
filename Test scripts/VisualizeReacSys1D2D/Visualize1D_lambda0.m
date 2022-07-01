%% Initialization
% Clear workspace, close all figures and clear command window
clear; close all; clc

% Add path to experimental plan function
addpath(strcat(erase(pwd, "Test scripts\VisualizeReacSys1D2D"), "\Functions"))
% addpath(strcat(erase(pwd, "Test scripts"), "\easyGSA"))

% Run custom figures script
figurer;
rng(0, 'twister');

%% User-defined decisions
bool_SaveFigures = true;

%%
reacstruc = reacstruccreate();

% Set point
reacstruc.process.T         = 20;
reacstruc.process.pH        = 11;
reacstruc.process.Co        = 20;
reacstruc.process.lambda0   = 2;
reacstruc.process.tdose     = 30;

% 
lambda0 = (1:0.1:3)';
% lambda0 = lhsdesign(10,1)*3;
% lambda0 = [1:0.05:2.5]';
% lambda0 = [1:0.05:2.5]';
% lambda0 = rand(10,1).* (2.5-1.0) + 1.0;
epsilon = normrnd(0, 2/25, length(lambda0), 1);
% epsilon = normrnd(0, 25/100, length(lambda0), 1);

yD_lambda0 = zeros(length(lambda0),1);
yH_lambda0 = zeros(length(lambda0),1);
for i = 1:length(lambda0)
    % Change molar ratio of sidechain S and component A
    reacstruc.process.lambda0 = lambda0(i)+epsilon(i);

    % Run simulation and save solutions in structure array
    reacstruc = reacsim(reacstruc);  

    % Record solution
    yD_lambda0(i) = reacstruc.out.y(end,7);        % Diacylated (product)
    yH_lambda0(i) = reacstruc.out.y(end,10);       % Triacylated (impurity)
end

%% Fit GPR
% Fit GPR model for different initial values of kernel parameters
GPR_Yield_1 = fitrgp(lambda0, yD_lambda0, 'KernelParameters', [0.01 0.1]);
GPR_Yield_2 = fitrgp(lambda0, yD_lambda0, 'KernelParameters', [0.5 0.1]);
GPR_Yield_3 = fitrgp(lambda0, yD_lambda0, 'KernelParameters', [0.01 1]);
GPR_Yield_4 = fitrgp(lambda0, yD_lambda0, 'KernelParameters', [0.5 1]);

% Generate test data
lambda0_test = (1:0.01:3)';
yD_lambda0_test = zeros(length(lambda0_test),1);
yH_lambda0_test = zeros(length(lambda0_test),1);
for i = 1:length(lambda0_test)
    % Change molar ratio of sidechain S and component A
    reacstruc.process.lambda0 = lambda0_test(i);

    % Run instant lab simulation and save solutions in structure array
    reacstruc = reacsim(reacstruc);  

    % Record solution
    yD_lambda0_test(i) = reacstruc.out.y(end,7);    % Diacylated (product)
    yH_lambda0_test(i) = reacstruc.out.y(end,10);   % Triacylated (impurity)
end

[GPR_Yield_pred1, ~, GPR_Yield_ci1] = predict(GPR_Yield_1, lambda0_test);
[GPR_Yield_pred2, ~, GPR_Yield_ci2] = predict(GPR_Yield_2, lambda0_test);
[GPR_Yield_pred3, ~, GPR_Yield_ci3] = predict(GPR_Yield_3, lambda0_test);
[GPR_Yield_pred4, ~, GPR_Yield_ci4] = predict(GPR_Yield_4, lambda0_test);

figure; hold all
plot(lambda0, yD_lambda0*100, 'o')
plot(lambda0_test, yD_lambda0_test*100, '-', 'LineWidth', 2)
plot(lambda0_test, GPR_Yield_pred1*100, '--', 'LineWidth', 2)
plot(lambda0_test, GPR_Yield_pred2*100, '--', 'LineWidth', 2)
plot(lambda0_test, GPR_Yield_pred3*100, '--', 'LineWidth', 2)
plot(lambda0_test, GPR_Yield_pred4*100, '--', 'LineWidth', 2)
% patch([lambda0_test; flipud(lambda0_test)], ...
%     [GPR_Yield_ci(:,1)*100; flipud(GPR_Yield_ci(:,2))*100], ....
%     'blue', 'FaceAlpha', .125) 
ylim([0 120])
xlabel('Molar ratio of S and A, $\lambda_{0}$', 'interpreter', 'latex')
ylabel('Yield [$\%$]')
legend("Data", "Ground truth", ...
               "GPR ($\ell = 0.01, \sigma_{f} = 0.1$)", ...
               "GPR ($\ell = 0.5, \sigma_{f} = 0.1$)", ...
               "GPR ($\ell = 0.01, \sigma_{f} = 1$)", ...
               "GPR ($\ell = 0.5, \sigma_{f} = 1$)", ...
               'interpreter', 'latex', 'location', 'southeast', 'FontSize', 12)
FigureTitle(1) = "1D_Lambda0_GPR";

% Optimize hyperparameters and compute initial estimate of kernel parameters
GPR_Yield = fitrgp(lambda0, yD_lambda0,             ...
        'KernelFunction', 'ardsquaredexponential',     ...
        'OptimizeHyperparameters', 'Sigma',  ...
        'Standardize', true,                 ...
        'HyperparameterOptimizationOptions', ...
        struct('ShowPlots', false,           ...
               'KFold', 5,                   ...
               'MaxObjectiveEvaluations', 50));

% % Set initial parameters
% Sigma0_Yield            = GPR_Yield_Init.Sigma;
% Sigma0_lb_Yield         = min(1e-2*std(yD_lambda0), GPR_Yield_Init.Sigma-1e-3);
% KernelParameters0_Yield = GPR_Yield_Init.KernelInformation.KernelParameters;
% 
% % Fit model 
% GPR_Yield = fitrgp(lambda0, yD_lambda0,                             ...
%                     'KernelFunction', 'ardsquaredexponential',      ...
%                     'KernelParameters', KernelParameters0_Yield,    ...
%                     'Sigma', Sigma0_Yield,                          ...
%                     'SigmaLowerBound', Sigma0_lb_Yield,             ...
%                     'Standardize', true,                            ...
%                     'ConstantSigma', true);

% Predictions
[GPR_Yield_pred, GPR_Yield_sds, GPR_Yield_ci] = predict(GPR_Yield, lambda0_test);

% Unwrap optimal kernel parameters
LengthScale      = GPR_Yield.KernelInformation.KernelParameters(1);
NoiseStandardDev = GPR_Yield.KernelInformation.KernelParameters(2);

figure; hold all
plot(lambda0, yD_lambda0*100, 'o')
plot(lambda0_test, yD_lambda0_test*100, '-', 'LineWidth', 2)
plot(lambda0_test, GPR_Yield_pred*100, '--', 'LineWidth', 2)
patch([lambda0_test; flipud(lambda0_test)], ...
    [GPR_Yield_ci(:,1)*100; flipud(GPR_Yield_ci(:,2))*100], ....
    'blue', 'FaceAlpha', .125) 
ylim([0 120])
xlabel('Molar ratio of S and A, $\lambda_{0}$', 'interpreter', 'latex')
ylabel('Yield [$\%$]')
legend("Data", "Ground truth", ...
               strcat("GPR ($\ell =", num2str(LengthScale), ", ", ...
                      "\sigma_{f} =", num2str(NoiseStandardDev),"$)"), ...
               "95 \% prediction interval", ...
               'interpreter', 'latex', 'location', 'southeast', 'FontSize', 12)
FigureTitle(2) = "1D_Lambda0_GPR_best";

%% Fit ANN
% Fit ANNs with different architectures
ANN_Yield_1 = fitrnet(lambda0, yD_lambda0);
ANN_Yield_2 = fitrnet(lambda0, yD_lambda0, 'LayerSizes', 100);
ANN_Yield_3 = fitrnet(lambda0, yD_lambda0, 'LayerSizes', [50 50], ...
                                           'Lambda', 1e-4);

ANN_Yield_pred1 = predict(ANN_Yield_1, lambda0_test);
ANN_Yield_pred2 = predict(ANN_Yield_2, lambda0_test);
ANN_Yield_pred3 = predict(ANN_Yield_3, lambda0_test);

figure; hold all
plot(lambda0, yD_lambda0*100, 'o')
plot(lambda0_test, yD_lambda0_test*100, '-', 'LineWidth', 2)
plot(lambda0_test, ANN_Yield_pred1*100, '-', 'LineWidth', 2)
plot(lambda0_test, ANN_Yield_pred2*100, '--', 'LineWidth', 2)
plot(lambda0_test, ANN_Yield_pred3*100, '--', 'LineWidth', 2)
xlabel('Molar ratio of S and A, $\lambda_{0}$', 'interpreter', 'latex')
ylabel('Yield [$\%$]')
ylim([0 120])
legend("Data", "Ground truth", ...
               "ANN ($n_{H_{1}} = 10, n_{H_{2}} = 0, \lambda = 0$)", ...
               "ANN ($n_{H_{1}} = 100, n_{H_{2}} = 0, \lambda = 0$)", ...
               "ANN ($n_{H_{1}} = 50, n_{H_{2}} = 50, \lambda = 1e-4$)", ...
               'interpreter', 'latex', 'location', 'southeast', 'FontSize', 12)
FigureTitle(3) = "1D_Lambda0_ANN";

% Optimize ANN hyperparameters
% Set hyperparameter objects
HyperParams_Yield = hyperparameters("fitrnet", lambda0, yD_lambda0);

% First object
HyperParams_Yield(1).Optimize = false;  % Number of layers
HyperParams_Yield(2).Optimize = false;   % Activation function
HyperParams_Yield(3).Optimize = false;  % Standardize
HyperParams_Yield(4).Optimize = true;   % Regularization strength
HyperParams_Yield(5).Optimize = false;  % Layer Weights Initializer
HyperParams_Yield(6).Optimize = false;  % Layer Biases Initializer
HyperParams_Yield(7).Optimize = true;   % Number of neurons layer 1
HyperParams_Yield(8).Optimize = true;   % Number of neurons layer 2
HyperParams_Yield(9).Optimize = false;  % Number of neurons layer 3
HyperParams_Yield(10).Optimize = false; % Number of neurons layer 4
HyperParams_Yield(11).Optimize = false; % Number of neurons layer 5

% Set hyperparameter ranges
HyperParams_Yield(7).Range = [10 30];   % Neurons layer 1
HyperParams_Yield(8).Range = [10 30];   % Neurons layer 2

% Set number of neurons in layers
NeuronsYield = [30 30];

% Optimize hyperparameters
ANN_Yield = fitrnet(lambda0, yD_lambda0,                   ...
                'LayerSizes', NeuronsYield,                     ...  
                'Activations', 'sigmoid',                       ...
                'Standardize', true,                            ...
                'OptimizeHyperparameters', HyperParams_Yield,   ...
                'HyperParameterOptimizationOptions',            ...
                struct("MaxObjectiveEvaluations", 50,           ...
                       'ShowPlots', false, 'KFold', 5));

% Fit final model
% ANN_Yield = fitrnet(lambda0, yD_lambda0, ...
%                     'LayerSizes', [15 13] , ...
%                     'Activations', 'sigmoid', ...
%                     'Standardize', true, ...
%                     'Lambda', 0.00018194);

% Predictions
ANN_Yield_pred = predict(ANN_Yield, lambda0_test);

figure; hold all
plot(lambda0, yD_lambda0*100, 'o')
plot(lambda0_test, yD_lambda0_test*100, '-', 'LineWidth', 2)
plot(lambda0_test, ANN_Yield_pred*100, '-', 'LineWidth', 2)
xlabel('Molar ratio of S and A, $\lambda_{0}$', 'interpreter', 'latex')
ylabel('Yield [$\%$]')
ylim([0 120])
legend("Data", "Ground truth", ...
    strcat("ANN ($n_{H_{1}} = ", num2str(ANN_Yield.LayerSizes(1)), ...
           ", n_{H_{2}} = ", num2str(ANN_Yield.LayerSizes(2)), ...
           ", \lambda = ", num2str(ANN_Yield.ModelParameters.Lambda), "$)"), ...
    'interpreter', 'latex', 'location', 'southeast', 'FontSize', 12)
FigureTitle(4) = "1D_Lambda0_ANN_best";

%% Fit MLR
MLR_Yield = fitlm(lambda0, yD_lambda0, 'quadratic');
MLR_Yield_pred = predict(MLR_Yield, lambda0_test);

%% Plot all best responses
figure; hold all
plot(lambda0, yD_lambda0*100, 'o')
plot(lambda0_test, yD_lambda0_test*100, '-', 'LineWidth', 2)
plot(lambda0_test, MLR_Yield_pred*100, '-', 'LineWidth', 2)
plot(lambda0_test, GPR_Yield_pred*100, '-', 'LineWidth', 2)
plot(lambda0_test, ANN_Yield_pred*100, '-', 'LineWidth', 2)
xlabel('Molar ratio of S and A, $\lambda_{0}$', 'interpreter', 'latex')
ylabel('Yield [$\%$]')
ylim([0 120])
legend("Data", "Ground truth",                              ...
    "MLR",                                                  ...
    strcat("ANN ($n_{H_{1}} = ", num2str(15),               ...
           ", n_{H_{2}} = ", num2str(13),                   ...
           ", \lambda = ", num2str(0.00018194), "$)"),      ...
    strcat("GPR ($\ell =", num2str(LengthScale), ", ",      ...
           "\sigma_{f} =", num2str(NoiseStandardDev),"$)"), ...
    'interpreter', 'latex', 'location', 'southeast')
FigureTitle(5) = "1D_Lambda0_MLR_GPR_ANN";

% Compute statistics
stats_MLR = rs_stats(yD_lambda0_test, MLR_Yield_pred);
stats_GPR = rs_stats(yD_lambda0_test, GPR_Yield_pred);
stats_ANN = rs_stats(yD_lambda0_test, ANN_Yield_pred);

%
disp('Statistics:')
disp([stats_MLR.R2 stats_GPR.R2 stats_ANN.R2; ...
      stats_MLR.SSE stats_GPR.SSE stats_ANN.SSE])

%% Save figures
if bool_SaveFigures
    for i = 1:length(FigureTitle)
        saveas(figure(i), FigureTitle(i), 'png')
    end
end

%% Using exact method directly
xnew    = lambda0_test;
xtrain  = lambda0;
beta    = GPR_Yield_1.Beta;
sigma_n = GPR_Yield_1.Sigma;
ell     = [0.01; 0.1; 10];
sigma_f = 1;
y       = yD_lambda0;

SqExpKernel = @(xi,xj,CL) sigma_f^2 * ...
                    exp(-1/2*((xi-xj)'*(xi-xj)) / CL^2);

K = zeros(length(xtrain), length(xtrain));
for i = 1:length(xtrain)
    for j = 1:length(xtrain)
        K(i,j) = SqExpKernel(xtrain(i), xtrain(j),ell(3));
    end
end
% 
% % alpha = (K + eye(length(xtrain))*sigma_n^2)^(-1) * (y-beta);
% alpha = (K + eye(length(xtrain))*sigma_n^2)^(-1) * (y-beta);
% 
% pred = zeros(length(xnew), length(xtrain));
% for i = 1:length(xnew)
%     for j = 1:length(xtrain)
%         pred(i,j) = alpha(j)*SqExpKernel(xnew(i),xtrain(j),ell(3));
%     end
% end
% 
% figure
% plot(lambda0_test, beta + sum(pred,2))

