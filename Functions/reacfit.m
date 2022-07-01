function reacstruc = reacfit(X, y)

% addpath(strcat(pwd, "\Functions"))

%% Main code
% Initialize simulation scenario 
reacstruc = reacstruccreate();

% Save data in data field of structure array 
for i = 1:size(X,1)
    reacstruc.data.input(i,:) = X(i,:);
    reacstruc.data.num{i}     = y(i,:);
end %i
con = [];

% For each parameter to be fitted do
for i = 1:length(reacstruc.parfit.par)
    % Define helper variable
    par   = reacstruc.parfit.par{i};    % Parameter name
    LB(i) = reacstruc.parfit.LB(i);     % Lower bound
    UB(i) = reacstruc.parfit.UB(i);     % Upper bound
    
    % Evaluate MATLAB expression
    z0(i) = eval(par);
    
    % Compute initial guess in deviation form
    x0(i) = (z0(i) - LB(i)) / (UB(i) - LB(i));
    
    % Set normalized constraints
    con(1,i) = 0;
    con(2,i) = 1;
end

% Convert vector of initial conditions to column vector
x0 = x0(:);

% Run start guess
f = reacobj(x0, reacstruc);

%%
% Define function handle and set optimizer tolerance
fun = @(x) reacobj(x, reacstruc);
options = optimset('TolX', 1e-9, 'TolFun', 1e-9);

% Compute solution to minimization problem with given constraints
% [sol, fval, exitflag, output] = fminsearchcon(fun, x0, ...
%                                               con(1,:), con(2,:), ...
%                                               [], [], [], options);
% ParameterMinimum = fminsearchcon(fun, x0,               ...
%                                  con(1,:), con(2,:),    ...
%                                  [], [], [], options);
ParameterMinimum = fmincon(fun, x0, [], [], [], [], con(1,:), con(2,:),    ...
                                 [], options);

% Define helper variable
par = reacstruc.parfit.par{i};  % Parameter name
LB  = reacstruc.parfit.LB;      % Lower bound
UB  = reacstruc.parfit.UB;      % Upper bound

x = ParameterMinimum;
z = zeros(length(reacstruc.parfit.par), 1);
for i = 1:length(reacstruc.parfit.par)
    z(i) = x(i) .* (UB(i) - LB(i)) + LB(i); %skalerer tilbage
    eval([reacstruc.parfit.par{i},'=',num2str(z(i),10),';']); %opdaterer model med ny z
end