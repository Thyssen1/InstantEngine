%% Tutorial: Performing Sobol GSA on Ishigami function
% By Resul Al @DTU

% Model: Ishigami function [https://www.sfu.ca/~ssurjano/ishigami.html]
f = @(x) sin(x(:,1)) + 7.*sin(x(:,2)).^2 + 0.1.*x(:,3).^4.*sin(x(:,1));  
N = 1e4; % Number of MC samples. Minimum recommended: 1e3

% Uniform Input Space
pars = {'x1','x2','x3'}; % input parameter names
lbs  = -pi.*ones(1,3);   % lower bounds of input parameters
ubs  =  pi.*ones(1,3);   % upper bounds of input parameters
InputSpace = {'ParNames',pars,'LowerBounds',lbs,'UpperBounds',ubs};

% call easyGSA tool to perform Sobol sensitivity analysis with MC approach
[mcSi,mcSTi] = easyGSA(f,N,InputSpace{:})
         
                
% Analytical first order indices from doi:10.1016/j.ress.2008.07.008
Si_analytical = [0.3139 0.4424 0]';

% Plot comparative results
T = table(Si_analytical,mcSi,...
    'VariableNames', {'Analytical','MonteCarlo'}, ...
    'RowNames', pars);
fprintf("\n\nFirst Order Sensitivity Indices of Ishigami function\n\n")
disp(T)

H = [Si_analytical,mcSi]; c = categorical(pars);
bar(1:3,H(:,1)); legend({'Analytical','MonteCarlo'});
ylabel('First Order Sobol indices'); xlabel('Input Parameters');
print('Si_ishigami','-dpng','-r1200')

H = [mcSTi,gprSTi,annSTi]; c = categorical(pars);
bar(c,H); legend({'MonteCarlo'});
ylabel('Total Order Sobol indices'); xlabel('Input Parameters');
print('STi_ishigami','-dpng','-r1200')