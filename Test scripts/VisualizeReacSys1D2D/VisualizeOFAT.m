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
bool_SaveFigures = false;

%%
reacstruc = reacstruccreate();

% Set point
reacstruc.process.T         = 20;
reacstruc.process.pH        = 11;
reacstruc.process.Co        = 20;
reacstruc.process.lambda0   = 2;
reacstruc.process.tdose     = 30;

%% Change lambda0
lambda0 = (1:0.1:3)';

yD_lambda0 = zeros(length(lambda0),1);
yH_lambda0 = zeros(length(lambda0),1);
for i = 1:length(lambda0)
    % Change molar ratio of sidechain S and component A
    reacstruc.process.lambda0 = lambda0(i);

    % Run simulation and save solutions in structure array
    reacstruc = reacsim(reacstruc);  

    % Record solution
    yD_lambda0(i) = reacstruc.out.y(end,7);        % Diacylated (product)
    yH_lambda0(i) = reacstruc.out.y(end,10);       % Triacylated (impurity)
end

figure
subplot(2,2,1)
plot(lambda0, yD_lambda0)
xlabel('lambda0')

%% Change temperature 
% Reset lambda0
reacstruc.process.lambda0   = 2;

% Set temperature
T = (5:1:40)';

yD_T = zeros(length(T),1);
yH_T = zeros(length(T),1);
for i = 1:length(T)
    % Change molar ratio of sidechain S and component A
    reacstruc.process.T = T(i);

    % Run simulation and save solutions in structure array
    reacstruc = reacsim(reacstruc);  

    % Record solution
    yD_T(i) = reacstruc.out.y(end,7);        % Diacylated (product)
    yH_T(i) = reacstruc.out.y(end,10);       % Triacylated (impurity)
end

subplot(2,2,2)
plot(T, yD_T, '-')
xlabel('Temperature')

%% Change Co
% Reset temperature
reacstruc.process.T   = 20;

% Set temperature
Co = (5:1:50)';

yD_Co = zeros(length(Co),1);
yH_Co = zeros(length(Co),1);
for i = 1:length(Co)
    % Change molar ratio of sidechain S and component A
    reacstruc.process.Co = Co(i);

    % Run simulation and save solutions in structure array
    reacstruc = reacsim(reacstruc);  

    % Record solution
    yD_Co(i) = reacstruc.out.y(end,7);        % Diacylated (product)
    yH_Co(i) = reacstruc.out.y(end,10);       % Triacylated (impurity)
end

subplot(2,2,3)
plot(Co, yD_Co, '-')
xlabel('Co')

%% Change pH
% Reset temperature
reacstruc.process.Co   = 20;

% Set temperature
pH = (10:0.5:12)';

yD_pH = zeros(length(pH),1);
yH_pH = zeros(length(pH),1);
for i = 1:length(pH)
    % Change molar ratio of sidechain S and component A
    reacstruc.process.pH = pH(i);

    % Run simulation and save solutions in structure array
    reacstruc = reacsim(reacstruc);  

    % Record solution
    yD_pH(i) = reacstruc.out.y(end,7);        % Diacylated (product)
    yH_pH(i) = reacstruc.out.y(end,10);       % Triacylated (impurity)
end

subplot(2,2,4)
plot(pH, yD_pH, '-')
xlabel('pH')

%%

figure
subplot(2,2,1)
plot(1 ./ (exp(-lambda0)+1), yD_lambda0)
subplot(2,2,2)
plot(T, yD_T)
subplot(2,2,3)
plot(1 ./ (exp(-0.1*Co)+1), yD_Co)
subplot(2,2,4)
plot(pH, yD_pH)



