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

% 
lambda0 = (1:0.01:3)';
pH      = 10.0:0.5:12;

yA = zeros(length(lambda0), length(pH));
yB = zeros(length(lambda0), length(pH));
yC = zeros(length(lambda0), length(pH));
yD = zeros(length(lambda0), length(pH));
yE = zeros(length(lambda0), length(pH));
yF = zeros(length(lambda0), length(pH));
yG = zeros(length(lambda0), length(pH));
yH = zeros(length(lambda0), length(pH));
lambda_react = zeros(length(lambda0), length(pH));
for i = 1:length(pH)
    for j = 1:length(lambda0)
        % Change molar ratio of sidechain S and component A
        reacstruc.process.pH      = pH(i);
        reacstruc.process.lambda0 = lambda0(j);
    
        % Run simulation and save solutions in structure array
        reacstruc = reacsim(reacstruc);  
    
        % Record solution
        yS(i,j) = reacstruc.out.y(end,1);
        yA(i,j) = reacstruc.out.y(end,3);
        yB(i,j) = reacstruc.out.y(end,4);
        yC(i,j) = reacstruc.out.y(end,5);
        yE(i,j) = reacstruc.out.y(end,6);
        yD(i,j) = reacstruc.out.y(end,7);        % Diacylated (product)
        yF(i,j) = reacstruc.out.y(end,8);
        yG(i,j) = reacstruc.out.y(end,9);
        yH(i,j) = reacstruc.out.y(end,10);       % Triacylated (impurity)

        lambda_react(i,j) = yS(i,j) + (yB(i,j) + yC(i,j) + yE(i,j)) + ...
                            + 2*(yD(i,j) + yF(i,j) + yG(i,j)) + 3*yH(i,j);
    end

    lgn(i) = strcat("pH=",num2str(pH(i)));
end


figure; hold all
for i = 1:length(pH)
    plot(lambda0, yD(i,:), '-')
end
xlabel('$\lambda_{0}$')
ylabel('Yield ($\%$)')
legend(lgn, 'location', 'northwest')

saveas(gcf, "GSA_1D_lambda0_pH", 'png')

