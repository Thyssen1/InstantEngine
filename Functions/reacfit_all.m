function reacstruc = reacfit_all(X, y)

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
fun = @(x) reacobj_res(x, reacstruc);
% fun = @(x) reacobj_yd(x, reacstruc);
options = optimset('TolX', 1e-6, 'TolFun', 1e-6, 'maxfunevals', 1000);

% Compute solution to minimization problem with given constraints
% [sol, fval, exitflag, output] = fminsearchcon(fun, x0, ...
%                                               con(1,:), con(2,:), ...
%                                               [], [], [], options);
ParameterMinimum = lsqnonlin(fun, x0, con(1,:), con(2,:), options);

% Define helper variable
% par = reacstruc.parfit.par{i};  % Parameter name
LB  = reacstruc.parfit.LB;      % Lower bound
UB  = reacstruc.parfit.UB;      % Upper bound

x = ParameterMinimum;
z = zeros(length(reacstruc.parfit.par), 1);
for i = 1:length(reacstruc.parfit.par)
    z(i) = x(i) .* (UB(i) - LB(i)) + LB(i); %skalerer tilbage
    eval([reacstruc.parfit.par{i},'=',num2str(z(i),10),';']); %opdaterer model med ny z
end


%% Old

% % Read from excel sheet
% [num, txt] = xlsread(strcat(pwd, "\Experimental plans\", "reaction.xls"),'input');
% nsim = size(num,1);
% 
% %% Simulation
% figure; hold all
% for k = 1:nsim
%     reacstruc.process.T       = data.nom_input(k,1); % Celcius
%     reacstruc.process.pH      = data.nom_input(k,2); % Logarithmic pH scale
%     reacstruc.process.Co      = data.nom_input(k,3); % g/L
%     reacstruc.process.lambda0 = data.nom_input(k,4); % mol /mol
%     reacstruc.process.tdose   = data.nom_input(k,5); % g/L
%     
%     % Run sim with new parameters
%     reacstruc = reacsim(reacstruc);
%     
%     % Extract yA
%     yD(k,1) = reacstruc.out.y(end,7);
%     
%     yD_act(k,1) = data.act_out{k}(6);
% end
% 
% plot(yD_act, yD_act, '-o')
% plot(yD_act, yD, 'x')
% xlabel('Experimental')
% ylabel('Simulation')
% 
% % Save parameters tor data structure array
% data.params = z;
% data.yD = yD;
% data.yD_act = yD_act;
% save(strcat(pwd, "\Data (results)\", folder, "\data_", exp_plan), 'data')

end