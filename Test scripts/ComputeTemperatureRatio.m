%% Initialization
% Clear workspace, close all figures and clear command window
clear; close all; clc

% Add path to experimental plan function
addpath(strcat(erase(pwd, "Test scripts"), "\Functions"))
addpath(strcat(erase(pwd, "Test scripts"), "\easyGSA"))

%
figurer;

%%
bool_SaveFigures = true;

% Set up simulation scenario
reacstruc = reacstruccreate();
T = 5:1:40;
pH = 10:0.5:12;

k = zeros(length(T), length(pH), 3);
for i = 1:length(T)
    for j = 1:length(pH)
        %
        reacstruc.process.T = T(i);
        reacstruc.process.pH = pH(j);

        % Run simulation
        reacstruc = reacsim(reacstruc);

        %
        for p = 1:3
            k(i,j,p) = reacstruc.model.k(p);
        end
    end
end

figure; hold all
for i = 1:length(pH)
    s1(:,i) = k(:,i,1) ./ (k(:,i,1) + k(:,i,2) + k(:,i,3));
    plot(T, 2*k(:,i,1) ./ (k(:,i,1) + k(:,i,2) + k(:,i,3)), '-')
    lgn(i) = strcat("pH = ", num2str(pH(i)));
end
ylim(2*[0.475 0.5])
xlabel('Temperature [$^{o}$C]', 'interpreter', 'latex', 'FontSize', 14)
ylabel('Selectivity, $s_{1}+s_{2}$ (pH, T)', 'interpreter', 'latex', 'FontSize', 14)
legend(lgn, 'location', 'southeast', 'FontSize', 12)
FigureTitle(1) = "GSA_TempRatio_s1";

figure; hold all
for i = 1:length(pH)
    s3(:,i) = k(:,i,3) ./ (k(:,i,1) + k(:,i,2) + k(:,i,3));
    plot(T, k(:,i,3) ./ (k(:,i,1) + k(:,i,2) + k(:,i,3)), '-')
end
ylim([0 0.025])
xlabel('Temperature [$^{o}$C]', 'interpreter', 'latex', 'FontSize', 14)
ylabel('Selectivity, $s_{3}$ (pH, T)', 'interpreter', 'latex', 'FontSize', 14)
legend(lgn, 'location', 'northeast', 'FontSize', 12)
FigureTitle(2) = "GSA_TempRatio_s3";


if bool_SaveFigures
    for i = 1:length(FigureTitle) 
        saveas(figure(i), ...
            strcat(pwd, "\Images\TemperatureRatio\", FigureTitle(i)), 'png')
    end
end

