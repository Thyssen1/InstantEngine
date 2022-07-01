function reacstruc=reacstruccreate()


%% Process parameters
% Define experiment name and name of components
reacstruc.process.name = 'NN9023 Acylation';
reacstruc.process.comps = {'SC','OH-','A','B','C','E','D','F','G','H','Vol'};

% Define temperature, pH and initial concentration
reacstruc.process.T  = 20;      % [C]   Temperature
reacstruc.process.pH = 11;    % [ ]   Concentration of H+ ions
reacstruc.process.Co = 20;      % [g/L] Initial concentration of material
                                % (backbone A)
reacstruc.model.sW1=0.5;

% Define purities of component A and sidechains, nominal mole ratio of side
% chains to component A, dosage time, Sialidase incubation time and dose
% concentration
% reacstruc.process.purityA  = 1; 
% reacstruc.process.SCpurity = 1;
reacstruc.process.lambda0  = 2.0;   % nominal mol SC/mol A
reacstruc.process.tdose    = 30;    % [min]     Dosing time
reacstruc.process.textra   = 10;    %           Sialidase incubation time
reacstruc.process.CSdose   = 333;   % dose conc in g/L, used for 
                                    % estimating dV/dt (Qdose)
                                    
%% Model parameters
% Define reference temperature and pH
reacstruc.model.Tref = 20;
reacstruc.model.pHref = 11;

% Define reference reaction rate for component A, reference selectivities
% for E and degradation of side chain, and exponent in OH 
reacstruc.model.kAref   = 2.0358e5;    % 1.17e5;%13.6e3; %9.9e4, kref=sum(kiR), reference T
reacstruc.model.sref(3) = 0.0029;   % 0.0013 for synth;% >0,  hvis triacyleret, s0(2) beregnes i sim
reacstruc.model.sdref   = 2e-4;     %  2.7832e-04kdi/k0 for sidec chain, pH 14, reference T
reacstruc.model.nd      = 1;        % OH exponent, n=1 default

% T dependency. Energy of activation for k1, k2, k3 and energy of
% activation for degradation
reacstruc.model.EA1 = 1e4*1.8; % Energy of activation for k1
% reacstruc.model.EA2 = 1e4*1.8; % Energy of activation for k1
reacstruc.model.EA2 = reacstruc.model.EA1;
reacstruc.model.EA3 = 1e4*1.4; % Energy of activation for k1
reacstruc.model.EA = [reacstruc.model.EA1 reacstruc.model.EA2 ...
                      reacstruc.model.EA3];
% reacstruc.model.EAd = 1e4*1.6;        % EA for kdeg
reacstruc.model.EAd = 1e4*1.6;          % EA for kdeg

%pH dep pH dependency
reacstruc.model.pkA1 = 10.75;
reacstruc.model.pkA2 = 10.75;
reacstruc.model.pkA3 = 9.2;
reacstruc.model.pkA = [reacstruc.model.pkA1 reacstruc.model.pkA2 ...
                       reacstruc.model.pkA3];    % Not fitted, Googled 
reacstruc.model.MWA = 3490.78;              % ELN 23642-075
reacstruc.model.MWS = 824.868;              % ELN 23642-075

% Define helper variables 
MWA = reacstruc.model.MWA;  % Molar weight of component A
MWS = reacstruc.model.MWS;  % Molar weight of component S (reagent)

% Compute molar weights
reacstruc.model.MW = [MWS 40 MWA (MWA+MWS)*[1 1 1] ...
                     (MWA + 2*MWS)*[1 1 1] (MWA + 3*MWS)]; 

%% solver options
reacstruc.solver.opts = odeset('initialstep',1e-9,'RelTol',1e-9,'AbsTol',1e-9);

%% Parameter fitting
% Define name of parameters for fitting 
reacstruc.parfit.par{1} = 'reacstruc.model.kAref'; 
reacstruc.parfit.par{2} = 'reacstruc.model.sref(3)'; 
reacstruc.parfit.par{3} = 'reacstruc.model.sdref'; 
reacstruc.parfit.par{4} = 'reacstruc.model.pkA1';
reacstruc.parfit.par{5} = 'reacstruc.model.pkA3'; 
reacstruc.parfit.par{6} = 'reacstruc.model.EA1'; 
reacstruc.parfit.par{7} = 'reacstruc.model.EA3'; 
reacstruc.parfit.par{8} = 'reacstruc.model.EAd';

% reacstruc.parfit.par{4} = 'reacstruc.model.pkA1'; 
% reacstruc.parfit.par{5} = 'reacstruc.model.pkA2'; 
% reacstruc.parfit.par{6} = 'reacstruc.model.pkA3'; 
% reacstruc.parfit.par{7} = 'reacstruc.model.EA1'; 
% reacstruc.parfit.par{8} = 'reacstruc.model.EA2'; 
% reacstruc.parfit.par{9} = 'reacstruc.model.EA3'; 
% reacstruc.parfit.par{10} = 'reacstruc.model.EAd';

% Define lower and upper bounds for parameter fitting
% reacstruc.parfit.LB = [0 0 0 0 0 0 0 0 0 0]';                  % Lower bounds
% reacstruc.parfit.UB = [1e7 .1 .01 14 14 14 1e5 1e5 1e5 1e5]';  % Upper bounds

reacstruc.parfit.LB = [0 0 0 0 0 0 0 0]';   % Lower bounds
reacstruc.parfit.UB = [1e7 .1 .01 14 14 1e5 1e5 1e5]';  % Upper bounds

%% Optimization parameters
% Compute price of sidechain (SC) per kg relative to precursor per kg
cost_relative = 50e3 / 22e3; 

reacstruc.optim.Price.SCrel = cost_relative * MWS / MWA; %SC cost/precursor cost , mol/mol
% reacstruc.optim.Price.SCrel = 2;
reacstruc.optim.var{1} = 'reacstruc.process.lambda0';   % Molar ratio
reacstruc.optim.var{2} = 'reacstruc.process.pH';
reacstruc.optim.var{3} = 'reacstruc.process.T';         % Temperature 
reacstruc.optim.var{4} = 'reacstruc.process.Co';        % Initial concentration
% reacstruc.optim.var{5} = 'reacstruc.process.tdose';     % Dosing time

% Define lower and upper bounds for optimization problem
% reacstruc.optim.LB = [1 10 5 5 10]';  % Lower bounds
% reacstruc.optim.UB = [3 12 40 50 60]';   % Upper bounds
reacstruc.optim.LB = [1 10 5 5]';  % Lower bounds
reacstruc.optim.UB = [3 12 40 50]';   % Upper bounds


%% InstantLab
reacstruc.instantlab.bool_Qdose = false;

