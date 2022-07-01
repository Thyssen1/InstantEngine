function f = reacobj(x,reacstruc)
%REACOBJ Summary of this function goes here
%   Detailed explanation goes here


%% Initialization

% Update reacstruc structure array with current parameter estimates
for i = 1:length(reacstruc.parfit.par)
    LB = reacstruc.parfit.LB(i);
    UB = reacstruc.parfit.UB(i);
    
    % Scale back the initial guess
    z(i) = x(i) .* (UB - LB) + LB;
    
    % Update parameter estimate with new z
    eval([reacstruc.parfit.par{i},'=',num2str(z(i),10),';']);
end

% Reorder as column vector
z = z(:)';

%% Main code
% Unwrap data from structure array
nom_input = reacstruc.data.input;

% Define helper variables
nsim = size(reacstruc.data.input,1);
Res(1:nsim, 1) = 0;

%ID	Description	Temp	pH	startconc (g/L)	lamda	tdose	Csdose	SC purity

% For nsim number each simulation do 
for k = 1:nsim
%     if (isfield(reacstruc.process, 'Qdose'))
%         reacstruc.process = rmfield(reacstruc.process, 'Qdose');
%     end
    
    % Update structure array for current simulation scenario
    reacstruc.process.T       = nom_input(k,1); % Celcius
    reacstruc.process.pH      = nom_input(k,2); % Logarithmic pH scale
    reacstruc.process.Co      = nom_input(k,3); % g/L
    reacstruc.process.lambda0 = nom_input(k,4); % mol /mol
    reacstruc.process.tdose   = nom_input(k,5); % g/L
    
    %Tid	A	mono	di	tri	eq.
    num = reacstruc.data.num{k}(1:5);
    num = [[0 1 0 0 0];num]; %t and y for t=0
    tdat = num(:,1);
    
%     reacstruc.process.tdose = tdat(end); %min
    ydat = 0*num(:,2:5);                % Allocation
    ydat(:,1) = num(:,2);               % A
    ydat(:,2) = num(:,3);               % Mono
    ydat(:,3) = num(:,4);               % Di
    ydat(:,4) = num(:,5);               % Tri
    xdat = 1 - ydat(:,1) - ydat(:,2);
    
    % Define weight
    W = 1 + 0 * xdat;
    
    % Run simulation and save solutions in structure array
    reacstruc = reacsim(reacstruc);
    
    % Unwrap solutions
    t = reacstruc.out.t;
    y = reacstruc.out.y;
    
    % 
    ysum = 0*y(:,1:4);              % Allocation
    ysum(:,1) = y(:,3);             % Component A
    ysum(:,2) = sum(y(:,4:6),2);    % Monoacylated
    ysum(:,3) = sum(y(:,7:9),2);    % Diacylated
    ysum(:,4) = y(:,10);            % Triacylated
    
    %
    yAout = interp1(t, ysum(:,1), tdat);
    yBout = interp1(t, ysum(:,2), tdat);
    yDout = interp1(t, ysum(:,3), tdat);
    yHout = interp1(t, ysum(:,4), tdat);

    ResA = sum(W.*(yAout(:) - ydat(:,1)).^2);
    ResB = sum(W.*(yBout(:) - ydat(:,2)).^2);
    ResD = sum(W.*(yDout(:) - ydat(:,3)).^2);
    ResH = sum(W.*(yHout(:) - ydat(:,4)).^2);

    Res(k,1) = ResA + ResB + ResD + 1*ResH;
end

% Compute value of objective function
f = sum(Res);
% f = Res;


