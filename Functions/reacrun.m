function reacstruc=reacrun()

%% Initialization
% Clear workspace and close all figures
clear; close all; clc

%% Main code
% Initialize simulation scenario 
reacstruc = reacstruccreate();

% reacstruc.model.kAref = 1.6e5;

reacstruc.process.T       = 20;
reacstruc.process.pH      = 11;
reacstruc.process.lambda0 = 2;
reacstruc.process.Co      = 20;
reacstruc.process.tdose   = 30;

% Read from excel sheet
[num, txt] = xlsread(strcat(pwd, "\Experimental Plans\reaction.xls"),'input');
nsim = size(num,1);

%% Simulation
% Compute nsim number of simulations
% for k = 1:nsim
for k = 1
    % k
    
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
    
    % Unwrap parameter from structure array
    lambda0 = reacstruc.process.lambda0;
    
    mono      = sum(y(:,4:6),2);
    di        = sum(y(:,7:9),2);
    tri       = sum(y(:,10),2);
    xsim      = (1-y(:,3)-mono);
    yend(k,:) = y(end,:);
    sdeg      = 1 - (mono(end) + 2*di(end) + 3*tri(end))/lambda0;
    
    % Update structure array
    reacstruc.out.sdeg = sdeg;

    % Figures
    factor = 100;
    
    figure(1); hold all
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

    out(k,:) = [k di(end) mono(end) tri(end) sdeg cost];
end%k
% '[k yield sdeg cost]'
% out
