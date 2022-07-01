%% Initialization
% Clear workspace, close all figures and clear command window
clear; close all; clc

% Add path to experimental plan function
addpath(strcat(erase(pwd, "Test scripts"), "\Functions"))

%
figurer;

%%
bool_SaveFigures = true;

% Set parameters
MWS         = 824.868;
MWA         = 3490.78;
lambda0     = 2;
CSdose_Mean = 333 / MWS;
tdose       = 30;
Co_Mean     = 20 / MWA;
N           = 1e4;       % Number of samples

% Compute dosing rate Q
Qdose = Co_Mean * lambda0 / (tdose * CSdose_Mean);

%
sigma_CSdose = CSdose_Mean / 25;
sigma_Co     = Co_Mean / 25;

rng(0, 'twister')
X = lhsdesign(N,2);

CSdose = icdf('Normal', X(:,1), CSdose_Mean, sigma_CSdose);
Co     = icdf('Normal', X(:,2), Co_Mean, sigma_Co);

lambda0 = Qdose * tdose * (CSdose ./ Co);

figure
histogram(lambda0, 'Normalization', 'probability')
xlabel('$\lambda_{0}$', 'interpreter', 'latex')
ylabel('Frequency')
FigureTitle(1) = "Lambda0Distribution";

if bool_SaveFigures
    for i = 1:length(FigureTitle) 
        saveas(gcf,strcat(pwd, '\Images\Lambda0distribution\', FigureTitle(i)), 'png')
    end
end