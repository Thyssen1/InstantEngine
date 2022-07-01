%% Initialization
% Clear workspace, close all figures and clear command window
clear; close all; clc

% Add path to experimental plan function
addpath(strcat(erase(pwd, "Test scripts\VisualizeReacSys1D2D"), "\Functions"))

% Run custom figures script
figurer;
rng(0, 'twister');

%% User-defined decisions
bool_SaveFigures = true;

%% Model 1
% Initialize data
x = [-4; -2.5; -1; 0; 1.5];
y = [-2; 0; 1; 2; -1];

% Train model
gprMdl1 = fitrgp(x, y);

% Define points for which to predict
xtest = (-5:0.01:5)';

% Compute predictions and 95% confidence interval
[preds, ~, CI] = predict(gprMdl1, xtest);

% Visualize results
figure; hold all
plot(x,y,'o')
plot(xtest,preds,'-')
patch([xtest; flipud(xtest)], [CI(:,1); flipud(CI(:,2))], 'blue', 'FaceAlpha', .125) 
xlabel('x')
ylabel('y')
ylim([-4 4])
legend("Data", "Predictions", "$95 \%$ CI", 'location', 'northwest')
box on
FigureTitle(1) = "1D_Methods_GPR_Model_1";

%% Model 2
x = [-4; -2.5; -1; 0; 1.5];
y = [-2; 0; 1; 2; 2];

% Train model
gprMdl2 = fitrgp(x, y);

% Define points for which to predict
xtest = (-5:0.01:5)';

% Compute predictions and 95% confidence interval
[preds, ~, CI] = predict(gprMdl2, xtest);

% Visualize results
figure; hold all
plot(x,y,'o')
plot(xtest,preds,'-')
patch([xtest; flipud(xtest)], [CI(:,1); flipud(CI(:,2))], 'blue', 'FaceAlpha', .125) 
xlabel('x')
ylabel('y')
legend("Data", "Predictions", "$95 \%$ CI", 'location', 'northwest')
ylim([-4 4])
box on
FigureTitle(2) = "1D_Methods_GPR_Model_2";


%% Notes



%% Save figures
if bool_SaveFigures
    for i = 1:length(FigureTitle)
        saveas(figure(i), strcat(strcat(pwd, "\Images\"),FigureTitle(i)), 'png')
    end
end



