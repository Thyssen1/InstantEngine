function f = reacobj_yh(x,reacstruc)
%REACOBJ Summary of this function goes here
%   Detailed explanation goes here


%% Initialization

% Update reacstruc structure array with current parameter estimates
for i = 1:length(reacstruc.parfit.par)
    % Update parameter estimate with new x
    eval([reacstruc.parfit.par{i},'=',num2str(x(i),10),';']);
end

%% Main code
% Unwrap data from structure array
nom_input = reacstruc.data.input;

% Define helper variables
nsim = size(reacstruc.data.input,1);
Res(1:nsim, 1) = 0;

%ID	Description	Temp	pH	startconc (g/L)	lamda	tdose	Csdose	SC purity

% For nsim number each simulation do 
for k = 1:nsim
    % Update structure array for current simulation scenario
    reacstruc.process.T       = nom_input(k,1); % Celcius
    reacstruc.process.pH      = nom_input(k,2); % Logarithmic pH scale
    reacstruc.process.Co      = nom_input(k,3); % g/L
    reacstruc.process.lambda0 = nom_input(k,4); % mol /mol
    reacstruc.process.tdose   = nom_input(k,5); % g/L
    
    %Tid	A	mono	di	tri	eq.
    tdat = reacstruc.data.num{k}(1);
    yH   = reacstruc.data.num{k}(5);
    
    % Run simulation and save solutions in structure array
    reacstruc = reacsim(reacstruc);
    
    % Unwrap solutions
    t = reacstruc.out.t;
    y = reacstruc.out.y;
    
    yHout = interp1(t, y(:,10), tdat);
    
    Res(k,1) = yHout - yH;
end

% Compute value of objective function
f = Res;