%% Initialization
% Clear workspace and close all figures
clear; close all; clc

%% Main code
% Initialize simulation scenario 
reacstruc = reacstruccreate();

% Read from excel sheet
[num, txt] = xlsread('reaction','input');

nsim = size(num,1);

%% Simulation
% Compute nsim number of simulations
for k = 1:nsim
    % k
    
    %%ID	Description	Temp	pH	startconc (g/L)	lamda	tdose
    reacstruc.process.T       = num(k,3); %C
    reacstruc.process.pH      = num(k,4); %
    reacstruc.process.Co      = num(k,5); % g/L
    reacstruc.process.lambda0 = num(k,6); % mol HEP/mol N9
    reacstruc.process.tdose   = num(k,7); % h
    
    T(k) = num(k,3);
    pH(k) = num(k,4);
    Co(k) = num(k,5);
    lambda0(k) = num(k,6);
    tdose(k) = num(k,7);
    
    % Run simulation and save solutions in structure array
    reacstruc = reacsim(reacstruc);
    
    % Unwrap solutions from structure array
    t = reacstruc.out.t;
    y = reacstruc.out.y;
    
    yD(k,1) = y(end,7);
end%k

min_range = min(min(yD));
max_range = min(1, max(max(yD))*1.2);

% Temperature
subplot(2,3,1); hold all
idx_T = [8 1 9];
plot(T(idx_T), yD(idx_T,1), 'ko-')
ylim([min_range max_range])
xlabel('Temperature [$^{o}$C]', 'interpreter', 'latex')
ylabel('$y_{D} [\%]$', 'interpreter', 'latex')
% legend("D")
% set(gca, 'YScale', 'log')

% pH
idx_pH = [10 1 11];
subplot(2,3,2); hold all
plot(sort(pH(idx_pH)), yD(idx_pH,1), 'ko-')
ylim([min_range max_range])
xlabel('pH', 'interpreter', 'latex')
ylabel('$y_{D} [\%]$', 'interpreter', 'latex')
% legend("D")
% set(gca, 'YScale', 'log')

% Co
idx_Co = [2 1 3];
subplot(2,3,3); hold all
plot(Co(idx_Co), yD(idx_Co), 'ko-')
ylim([min_range max_range])
xlabel('Co [g/L]', 'interpreter', 'latex')
ylabel('$y_{D} [\%]$', 'interpreter', 'latex')
% legend("D", 'location', 'southeast')
% set(gca, 'YScale', 'log')

% lambda0
idx_lambda0 = [4 1 5];
subplot(2,3,4); hold all
plot(lambda0(idx_lambda0), yD(idx_lambda0,1), 'ko-')
ylim([min_range max_range])
xlabel('$\lambda_{0}$', 'interpreter', 'latex')
ylabel('$y_{D} [\%]$', 'interpreter', 'latex')
% legend("D")
% set(gca, 'YScale', 'log')

% tdose
idx_tdose = [6 1 7];
subplot(2,3,5); hold all
plot(tdose(idx_tdose), yD(idx_tdose,1), 'ko-')
ylim([min_range max_range])
xlabel('t$_{dose}$ [min]', 'interpreter', 'latex')
ylabel('$y_{D} [\%]$', 'interpreter', 'latex')
% legend("D", 'location', 'southeast')
% set(gca, 'YScale', 'log')

saveas(gcf, "reac_sysphys", 'png')