function reacstruc=reacrun(varargin)

%% Initialization
% Clear workspace and close all figures
clear; close all; clc

%% Main code
% Initialize simulation scenario 
if (exist('reacstruc'))
    % Use reacstruct given by user as input
else
    reacstruc = reacstruccreate();
end

% Read from excel sheet
% [num, txt] = xlsread('reaction','input');
[num, txt] = xlsread(strcat(erase(pwd, "Functions"), "\Experimental plans\reaction"), 'input');
nsim = size(num,1);                             % Number of experiments

%% Simulation
% Compute nsim number of simulations
for k = 1:nsim
    %%ID	Description	Temp	pH	startconc (g/L)	lamda	tdose
    reacstruc.process.T       = num(k,3); %C
    reacstruc.process.pH      = num(k,4); %
    reacstruc.process.Co      = num(k,5); % g/L
    reacstruc.process.lambda0 = num(k,6); % mol HEP/mol N9
    reacstruc.process.tdose   = num(k,7); % h
    
    % Run simulation and save solutions in structure array
    reacstruc = reacsim(reacstruc);
    
    % Unwrap solutions from structure array
    t = reacstruc.out.t;
    y = reacstruc.out.y;
    
    % Unwrap molar ratio from structure array
    lambda0 = reacstruc.process.lambda0;
    
    % Compute amount (%) of mono-, di-, and triacylated components
    mono      = sum(y(:,4:6),2);    % Monoacylated 
    di        = sum(y(:,7:9),2);    % Diacylated
    tri       = sum(y(:,10),2);     % Triacylated
    xsim      = (1-y(:,3)-mono);    %
    yend(k,:) = y(end,:);           %
    
    % Compute selectivity of degradation
    reacstruc.out.sdeg = 1 - (mono(end) + 2*di(end) + 3*tri(end))/lambda0; 

    % Figures
    factor = 100;
    
    figure; hold all
    title('Time profiles')
    
    % Add plots to figure
    plot(t, factor * y(:,3),'-b')
    plot(t, factor * mono(:),'-g')
    plot(t, factor * di(:),'r-')
    plot(t, factor * tri(:),'-m')
    
    % Add limits and labels
    ylim([0 factor])
    xlabel('Time [min]')
    legend("Reagent", "Monoacylated", "Diacylated", "Triacylated")
    
    % Conversion plot
    figure; hold all
    title('Conversion plot')
    
    % Add plots to figure
    plot(factor * xsim, factor*y(:,3),'-b')
    plot(factor * xsim, factor*mono(:),'-g')
    plot(factor * xsim, factor*di(:),'r-')
    plot(factor * xsim, factor*tri(:),'-m')
    
    % Add limits
    ylim([0 factor])
    xlim([0 factor])
    xlabel('Time [min]')
    legend("Component A", "Mono", "Di", "Tri", 'location', 'north')
    
    % Unwrap parameters from structure array
    Price = reacstruc.optim.Price;
    cost = (1 + Price.SCrel*lambda0) / di(end);    % cost/cost of A
    sdeg = reacstruc.out.sdeg;

    out(k,:) = [k di(end) mono(end) tri(end) sdeg cost];
end
