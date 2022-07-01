%% Initialization
% Clear workspace, close all figures and clear command window
clear; close all; clc

% Add path to experimental plan function
addpath(strcat(erase(pwd, "Test scripts"), "\Functions"))
addpath(strcat(erase(pwd, "Test scripts"), "\easyGSA"))

%
figurer;

rng(4, 'twister')

%% User-defined decisions
bool_SaveFigures = false;
SimulationIndex = 1;
VariableNames = ["T", "pH", "$C_{A,0}$", "$\lambda_{0}$", "$t_{dose}$"];

% Set parameters for running experiment
setpoints = [20 11 20 2 30];
CV_input  = [inf inf inf inf; ...
             10 250 25 25];     % T pH Co CSdose
for i = 1:size(CV_input,1)
    sigma_input(i).T      = setpoints(1) / CV_input(i,1);
    sigma_input(i).pH     = setpoints(2) / CV_input(i,2);
    sigma_input(i).Co     = setpoints(3) / CV_input(i,3);
    sigma_input(i).CSdose = 333 / CV_input(i,4);
end
CV_output = [10 10 10 10 10];   % yA ymono ydi ytri yD

% Set parameters for GSA
N          = 1e5; 
pars       = {'T', 'pH', 'Co', 'lambda0', 'tdose'};
lbs        = [5 10 5 1 10];
ubs        = [40 12 50 3 60];
InputSpace = {'ParNames', pars, 'LowerBounds', lbs, 'UpperBounds', ubs};

%% One-factor-at-a-time method
% Set experimental parameters and generate experimental plan
OFAT.design    = 'OFAT';          % Experimental Design
OFAT.setpoints = [20 11 20 2 30]; % [T pH Co lambda0 tdose]
OFAT.lows      = [5 10 5 1 10];   % Low values
OFAT.highs     = [40 12 50 3 60]; % High values
OFAT.n         = 5;               % Number of factors
OFAT.seed      = 1;               % Set seed for random number generator
OFAT.plan      = exp_plan(OFAT);  % InstantLab lab plan

% Run experiments
OFAT.data = instantlab(OFAT.plan, sigma_input(SimulationIndex), OFAT.seed);

% Pre-allocate arrays
OFAT.T       = zeros(size(OFAT.data.out,2), 1);
OFAT.pH      = zeros(size(OFAT.data.out,2), 1);
OFAT.Co      = zeros(size(OFAT.data.out,2), 1);
OFAT.lambda0 = zeros(size(OFAT.data.out,2), 1);
OFAT.tdose   = zeros(size(OFAT.data.out,2), 1);
OFAT.yD      = zeros(size(OFAT.data.out,2), 1);
OFAT.yH      = zeros(size(OFAT.data.out,2), 1);
OFAT.D       = zeros(size(OFAT.data.out,2), 7);

% Populate arrays with InstantLab data
for i = 1:size(OFAT.data.out,2)
    OFAT.T(i, 1)       = OFAT.data.nom_input(i, 1);
    OFAT.pH(i, 1)      = OFAT.data.nom_input(i, 2);
    OFAT.Co(i, 1)      = OFAT.data.nom_input(i, 3);
    OFAT.lambda0(i, 1) = OFAT.data.nom_input(i, 4);
    OFAT.tdose(i, 1)   = OFAT.data.nom_input(i, 5);
    OFAT.yD(i, 1)      = OFAT.data.out{i}(6);
    OFAT.yH(i, 1)      = OFAT.data.out{i}(5);
    
    OFAT.D(i,:) = [OFAT.yD(i) OFAT.yH(i) OFAT.T(i) OFAT.pH(i) ...
                   OFAT.Co(i) OFAT.lambda0(i) OFAT.tdose(i)];
end

% Define design
OFAT.Yield = OFAT.D(:,1);   % Measured response
OFAT.Impurity = OFAT.D(:,2);
OFAT.X = OFAT.D(:,3:7); % Predictors

                          
% Perform stepwise multiple linear regression 
OFAT.ModelYield = stepwiselm(OFAT.X, OFAT.Yield, 'linear', ...
                                      'Criterion', 'AIC', ...
                                      'Upper', 'linear', ...
                'VarNames', {'T','pH','Co','lambda0','tdose','yield'});
OFAT.ModelYield

OFAT.ModelImpurity = stepwiselm(OFAT.X, OFAT.Impurity, 'linear', ...
                                      'Criterion', 'AIC', ...
                                      'Upper', 'linear', ...
                'VarNames', {'T','pH','Co','lambda0','tdose','impurity'});
OFAT.ModelImpurity

% Define prediction functions
f_OFAT_Yield    = @(x) predict(OFAT.ModelYield, x);
f_OFAT_Impurity = @(x) predict(OFAT.ModelImpurity, x);

% Compute Sobol indices
[Si_OFAT_Yield, STi_OFAT_Yield]       = easyGSA(f_OFAT_Yield, N, InputSpace{:});
[Si_OFAT_Impurity, STi_OFAT_Impurity] = easyGSA(f_OFAT_Impurity, N, InputSpace{:});

disp([Si_OFAT_Yield STi_OFAT_Yield])
disp([Si_OFAT_Impurity, STi_OFAT_Impurity])

% figure
% for i = 1:5
%     subplot(2,3,i); hold all
%     plot(OFAT.X(:,i), OFAT.Yield, 'o')
%     xlabel(VariableNames(i), 'interpreter', 'latex')
%     ylabel('Yield, $y_{D}$', 'interpreter', 'latex')
%     ylim([0 1])
% end
% sgtitle('One-factor-at-a-time')
% FigureTitle(1) = "DOE_Screening_OFAT";
% 
% figure
% for i = 1:5
%     subplot(2,3,i); hold all
%     plot(OFAT.X(:,i), OFAT.Impurity, 'o')
%     xlabel(VariableNames(i), 'interpreter', 'latex')
%     ylabel('Impurity, $y_{H}$', 'interpreter', 'latex')
%     ylim([0 1])
% end
% sgtitle('One-factor-at-a-time')
% FigureTitle(1) = "DOE_Screening_OFAT";

%% Two-level full factor design
% Set experimental parameters and generate experimental plan
ff2.design    = 'ff2';           % ff2 ff2frac ff3 dsd OAT13
ff2.setpoints = [20 11 20 2 30]; % [T pH lambda0 Co tdose]
ff2.lows      = [5 10 5 1 10];   % Low values
ff2.highs     = [40 12 50 3 60]; % High values

% Set seed for random number generator
ff2.seed = 2;

% 
ff2.plan = exp_plan(ff2);

% Run experiments
ff2.data = instantlab(ff2.plan, sigma_input(SimulationIndex), ff2.seed);

% Pre-allocate
ff2.T       = zeros(size(ff2.data.out,2), 1);
ff2.pH      = zeros(size(ff2.data.out,2), 1);
ff2.Co      = zeros(size(ff2.data.out,2), 1);
ff2.lambda0 = zeros(size(ff2.data.out,2), 1);
ff2.tdose   = zeros(size(ff2.data.out,2), 1);
ff2.D       = zeros(size(ff2.data.out,2), 7);
ff2.yD      = zeros(size(ff2.data.out,2), 1);
ff2.yH      = zeros(size(ff2.data.out,2), 1);

for i = 1:size(ff2.data.out,2)
    ff2.T(i, 1)       = ff2.data.nom_input(i, 1);
    ff2.pH(i, 1)      = ff2.data.nom_input(i, 2);
    ff2.Co(i, 1)      = ff2.data.nom_input(i, 3);
    ff2.lambda0(i, 1) = ff2.data.nom_input(i, 4);
    ff2.tdose(i, 1)   = ff2.data.nom_input(i, 5);
    ff2.yD(i, 1)      = ff2.data.out{i}(6);
    ff2.yH(i, 1)      = ff2.data.out{i}(5);
    
    ff2.D(i,:) = [ff2.yD(i) ff2.yH(i) ff2.T(i) ff2.pH(i) ff2.Co(i) ...
                  ff2.lambda0(i) ff2.tdose(i)];
end

% Compute response and design matrix
ff2.Yield    = ff2.D(:,1);      % Measured response
ff2.Impurity = ff2.D(:,2);      % Measured impurity
ff2.X        = ff2.D(:,3:7); % Predictors

                          
% Perform stepwise multiple linear regression 
ff2.ModelYield = stepwiselm(ff2.X, ff2.Yield, 'linear', 'Upper', 'interactions', ...
                'VarNames', {'T','pH','Co','lambda0','tdose','yield'}, ...
                'Criterion', 'AIC');
ff2.ModelYield

ff2.ModelImpurity = stepwiselm(ff2.X, ff2.Impurity, 'linear', 'Upper', 'interactions', ...
                'VarNames', {'T','pH','Co','lambda0','tdose','yield'}, ...
                'Criterion', 'AIC');
ff2.ModelImpurity

% Define prediction functions
f_ff2_Yield    = @(x) predict(ff2.ModelYield, x);
f_ff2_Impurity = @(x) predict(ff2.ModelImpurity, x);

% Compute Sobol indices
[Si_ff2_Yield, STi_ff2_Yield]       = easyGSA(f_ff2_Yield, N, InputSpace{:});
[Si_ff2_Impurity, STi_ff2_Impurity] = easyGSA(f_ff2_Impurity, N, InputSpace{:});

disp([Si_ff2_Yield STi_ff2_Yield])
disp([Si_ff2_Impurity STi_ff2_Impurity])

%% Latin Hypercube Sampling
% Set experimental parameters and generate experimental plan
lhs.design  = 'lhs';            % Latin hypercube sampling
lhs.lows    = [5 10 5 1 10];    % Low values
lhs.highs   = [40 12 50 3 60];  % High values
lhs.samples = 32;               % Number of samples to draw from distribution
lhs.seed    = 1;                % Seed for random number generator
lhs.plan    = exp_plan(lhs);    % Experimental plan

% Run experiments
lhs.data = instantlab(lhs.plan, sigma_input(SimulationIndex), lhs.seed);

% Pre-allocate
lhs.T       = zeros(size(lhs.data.out,2), 1);
lhs.pH      = zeros(size(lhs.data.out,2), 1);
lhs.Co      = zeros(size(lhs.data.out,2), 1);
lhs.lambda0 = zeros(size(lhs.data.out,2), 1);
lhs.tdose   = zeros(size(lhs.data.out,2), 1);
lhs.D       = zeros(size(lhs.data.out,2), 7);
lhs.yD      = zeros(size(lhs.data.out,2), 1);
lhs.yH      = zeros(size(lhs.data.out,2), 1);

for i = 1:size(lhs.data.out,2)
    lhs.T(i, 1)       = lhs.data.nom_input(i, 1);
    lhs.pH(i, 1)      = lhs.data.nom_input(i, 2);
    lhs.Co(i, 1)      = lhs.data.nom_input(i, 3);
    lhs.lambda0(i, 1) = lhs.data.nom_input(i, 4);
    lhs.tdose(i, 1)   = lhs.data.nom_input(i, 5);
    lhs.yD(i, 1)      = lhs.data.out{i}(6);
    lhs.yH(i, 1)      = lhs.data.out{i}(5);
    
    lhs.D(i,:) = [lhs.yD(i) lhs.yH(i) ...
                  lhs.T(i) lhs.pH(i) lhs.Co(i) lhs.lambda0(i) lhs.tdose(i)];
end

% Compute response and design matrix
lhs.Yield    = lhs.D(:,1);      % Measured response
lhs.Impurity = lhs.D(:,2);      % Measured impurity
lhs.X        = lhs.D(:,3:7);    % Predictors

% Perform stepwise multiple linear regression 
lhs.ModelYield = stepwiselm(lhs.X, lhs.Yield, 'linear', 'Upper', 'quadratic', ...
                                   'Criterion', 'AIC', ...
                'VarNames', {'T', 'pH','Co','lambda0','tdose','yield'});
lhs.ModelYield

lhs.ModelImpurity = stepwiselm(lhs.X, lhs.Impurity, 'linear', 'Upper', 'quadratic', ...
                                   'Criterion', 'AIC', ...
                'VarNames', {'T', 'pH','Co','lambda0','tdose','yield'});
lhs.ModelImpurity

% Define prediction functions
f_lhs_Yield    = @(x) predict(lhs.ModelYield, x);
f_lhs_Impurity = @(x) predict(lhs.ModelImpurity, x);

% Compute Sobol indices
[Si_lhs_Yield, STi_lhs_Yield]       = easyGSA(f_lhs_Yield, N, InputSpace{:});
[Si_lhs_Impurity, STi_lhs_Impurity] = easyGSA(f_lhs_Impurity, N, InputSpace{:});

disp([Si_lhs_Yield STi_lhs_Yield])
disp([Si_lhs_Impurity STi_lhs_Impurity])


%% Definitive screening dsign
% Set experimental parameters and generate experimental plan
dsd.design    = 'dsd';           % ff2 ff2frac ff3 dsd OAT13
dsd.setpoints = [20 11 20 2 30]; % [T pH lambda0 Co tdose]
dsd.lows      = [5 10 5 1 10];   % Low values
dsd.highs     = [40 12 50 3 60]; % High values

% Set seed for random number generator
dsd.seed = 2;

% 
dsd.plan = exp_plan(dsd);

% Run experiments
dsd.data = instantlab(dsd.plan, sigma_input(SimulationIndex), dsd.seed);

% Pre-allocate
dsd.T       = zeros(size(dsd.data.out,2), 1);
dsd.pH      = zeros(size(dsd.data.out,2), 1);
dsd.Co      = zeros(size(dsd.data.out,2), 1);
dsd.lambda0 = zeros(size(dsd.data.out,2), 1);
dsd.tdose   = zeros(size(dsd.data.out,2), 1);
dsd.D       = zeros(size(dsd.data.out,2), 7);
dsd.yD      = zeros(size(dsd.data.out,2), 1);
dsd.yH      = zeros(size(dsd.data.out,2), 1);

for i = 1:size(dsd.data.out,2)
    dsd.T(i, 1)       = dsd.data.nom_input(i, 1);
    dsd.pH(i, 1)      = dsd.data.nom_input(i, 2);
    dsd.Co(i, 1)      = dsd.data.nom_input(i, 3);
    dsd.lambda0(i, 1) = dsd.data.nom_input(i, 4);
    dsd.tdose(i, 1)   = dsd.data.nom_input(i, 5);
    dsd.yD(i, 1)      = dsd.data.out{i}(6);
    dsd.yH(i, 1)      = dsd.data.out{i}(5);
    
    dsd.D(i,:) = [dsd.yD(i) dsd.yH(i) dsd.T(i) dsd.pH(i) dsd.Co(i) ...
                  dsd.lambda0(i) dsd.tdose(i)];
end

% Compute response and design matrix
dsd.Yield    = dsd.D(:,1);      % Measured response
dsd.Impurity = dsd.D(:,2);      % Measured impurity
dsd.X        = dsd.D(:,3:7); % Predictors

                          
% Perform stepwise multiple linear regression 
dsd.ModelYield = stepwiselm(dsd.X, dsd.Yield, 'linear', 'Upper', 'quadratic', ...
                'VarNames', {'T','pH','Co','lambda0','tdose','yield'}, ...
                'Criterion', 'AIC');
dsd.ModelYield

dsd.ModelImpurity = stepwiselm(dsd.X, dsd.Impurity, 'linear', 'Upper', 'quadratic', ...
                'VarNames', {'T','pH','Co','lambda0','tdose','yield'}, ...
                'Criterion', 'AIC');
dsd.ModelImpurity

% Define prediction functions
f_dsd_Yield    = @(x) predict(dsd.ModelYield, x);
f_dsd_Impurity = @(x) predict(dsd.ModelImpurity, x);

% Compute Sobol indices
[Si_dsd_Yield, STi_dsd_Yield]       = easyGSA(f_dsd_Yield, N, InputSpace{:});
[Si_dsd_Impurity, STi_dsd_Impurity] = easyGSA(f_dsd_Impurity, N, InputSpace{:});

disp([Si_dsd_Yield STi_dsd_Yield])
disp([Si_dsd_Impurity STi_dsd_Impurity])


%% Save figures
if bool_SaveFigures
    for i = 1:length(FigureTitle) 
        saveas(figure(i), strcat(pwd, "\Images\DOE_Screening\", ...
            FigureTitle(i), "NoiseLevel_", num2str(SimulationIndex)), 'png')
    end
end


