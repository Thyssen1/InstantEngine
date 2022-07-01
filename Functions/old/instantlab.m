function reacstruc = instantlab(exp_plan, folder, seed, N)

%% Initialization
% Read experimental plans from excel sheet
tab_params = readtable(strcat(pwd,"\Experimental plans\",exp_plan,".xls"), ...
            'Sheet', 'Input', 'VariableNamingRule', 'preserve');
params = table2array(tab_params(1:end,3:end));

tab_uncert_input = readtable(strcat(pwd,"\Experimental plans\",exp_plan,".xls"), ...
            'Sheet', 'Input Uncertainty', 'VariableNamingRule', 'preserve');
uncert_input = table2array(tab_uncert_input);

tab_uncert_output = readtable(strcat(pwd,"\Experimental plans\",exp_plan,".xls"), ...
            'Sheet', 'Output Uncertainty', 'VariableNamingRule', 'preserve');
uncert_output = table2array(tab_uncert_output);

names = table2array(tab_params(:,2));

% Set helper variable
nexp = size(params,1);     % Number of experiments

% Set seed for random number generator
if nargin > 2
    % Set seed for reproducibility
    rng(seed);
else
    % Do not set a seed
end

%% Main code
% Compute simulation of each experiment with no noise
for k = 1:nexp    
    % Initialize simulation scenario 
    reacstruc = reacstruccreate();

    % ID Description Temp pH startconc (g/L) lamda tdose
    reacstruc.process.T       = params(k,1); % C
    reacstruc.process.pH      = params(k,2); %
    reacstruc.process.Co      = params(k,3); % g/L
    reacstruc.process.lambda0 = params(k,4); % mol HEP/mol N9
    reacstruc.process.tdose   = params(k,5); % h
    
    % Save nominal input (without noise)
    data.nom_input(k,:) = [reacstruc.process.T ...
                reacstruc.process.pH reacstruc.process.Co ...
                reacstruc.process.lambda0 reacstruc.process.tdose];
    
    % Compute dosing flow rate to get correct molar ratio
    CSdose   = reacstruc.process.CSdose;    % Nominal CSdose
    MWS      = reacstruc.model.MWS;         %  
    MWA      = reacstruc.model.MWA;         %
    lambda0  = reacstruc.process.lambda0;   % Nominal molar ratio
    Ao       = reacstruc.process.Co / MWA;  % Nominal CA at t=0
    tdose    = reacstruc.process.tdose;     % [min]     Dosing time
    
    % Compute dosing flow rate of sidechain S
    Qdose = (tdose>0)*Ao*lambda0 / (eps+tdose*(CSdose / MWS)); %[start vols]/min
    
                 
    % Unwrap noise input parameters
    sigma_T       = uncert_input(k,1);    % 
    sigma_pH      = uncert_input(k,2);    % 
    sigma_Co      = uncert_input(k,3);    % 
    sigma_CSdose  = uncert_input(k,4);    % 
    
    % Generate normally distributed noise
    x_T       = normrnd(0, params(k,1)*sigma_T);
    x_pH      = normrnd(0, params(k,2)*sigma_pH);
    x_Co      = normrnd(0, params(k,3)*sigma_Co);
    x_CSdose  = normrnd(0, CSdose*sigma_CSdose);
    
    % Initialize simulation scenario 
    reacstruc = reacstruccreate();
    reacstruc.process.Qdose = Qdose;

    % ID Description Temp pH startconc (g/L) lamda tdose
    reacstruc.process.T       = params(k,1) + x_T; % C
    reacstruc.process.pH      = params(k,2) + x_pH; %
    reacstruc.process.Co      = params(k,3) + x_Co; % g/L
    reacstruc.process.CSdose  = reacstruc.process.CSdose + x_CSdose;      
    reacstruc.process.tdose   = params(k,5); % h
    
    % Create helper variables
    CSdose = reacstruc.process.CSdose;  % Nominal CSdose
    tdose  = reacstruc.process.tdose;   
    Co     = reacstruc.process.Co / MWA;      
    
    % Compute molar ratio with uncertainties in CSdose and Co
    reacstruc.process.lambda0 = (Qdose * (CSdose/MWS) * tdose) / Co; % mol HEP/mol N9
    
    % Save actual input (with noise)
    data.act_input(k,:) = [reacstruc.process.T       ...
                           reacstruc.process.pH      ...
                           reacstruc.process.Co      ...
                           reacstruc.process.lambda0 ...
                           reacstruc.process.tdose];
    
    % Run instant lab simulation and save solutions in structure array
    reacstruc = reacsim(reacstruc);               
              
    % Unwrap from solution vectors
    tend    = reacstruc.out.t(end);
    yA      = reacstruc.out.y(end,3);           % Component A
    yD      = reacstruc.out.y(end,7);           % Component D, product
    ymono   = sum(reacstruc.out.y(end,4:6),2);  % Monoacylated
    ydi     = sum(reacstruc.out.y(end,7:9),2);  % Diacylated
    ytri    = reacstruc.out.y(end,10);          % Triacylated
    
    data.nom_out{k} = [tend yA ymono ydi ytri yD];  
    
    % Unwrap measurement noise parameters
    sigma_yA     = uncert_output(k,1);    % 
    sigma_yD     = uncert_output(k,2);
    sigma_ymono  = uncert_output(k,3);    %
    sigma_ydi    = uncert_output(k,4);    %
    sigma_ytri   = uncert_output(k,5);    %
    
    % Generate output (measurement) noise
    noise = [0 normrnd(0, yA*sigma_yA)          ...
               normrnd(0, ymono*sigma_ymono)    ...
               normrnd(0, ydi*sigma_ydi)        ...
               normrnd(0, ytri*sigma_ytri)      ...
               normrnd(0, yD*sigma_yD)];
    
    % Add to measurement
    data.act_out{k} = [tend yA ymono ydi ytri yD] + noise; 
end

% Save names of experiments in structure array
data.names = names;

% Save data 
save(strcat(pwd, "\Data (results)\", folder, "\data_", exp_plan), 'data')
end