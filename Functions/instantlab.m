function data = instantlab(exp_plan, sigma_input, seed)

%% Initialization
% Read experimental plans from excel sheet
% if ~exist('CV_input')
%     CV_input = [10 10 10 10];
% end
% if ~exist('CV_output')
%     CV_output = [10 10 10 10 10];
% end

% Set helper variable
nexp = size(exp_plan,1);     % Number of experiments

% Set seed for random number generator
if ~exist('seed')
    % Set seed for reproducibility
    rng(seed);
else
    % Do not set a seed
end

% Programmer debug mode
bool_debug = false;

%% Main code
% Compute simulation of each experiment with no noise
for k = 1:nexp    
    % Initialize simulation scenario 
    reacstruc = reacstruccreate();

    % ID Description Temp pH startconc (g/L) lamda tdose
    reacstruc.process.T       = exp_plan(k,1); % C
    reacstruc.process.pH      = exp_plan(k,2); %
    reacstruc.process.Co      = exp_plan(k,3); % g/L
    reacstruc.process.lambda0 = exp_plan(k,4); % mol HEP/mol N9
    reacstruc.process.tdose   = exp_plan(k,5); % h
    
    % Save nominal input (without noise)
    data.nom_input(k,:) = [reacstruc.process.T       ...
                           reacstruc.process.pH      ...
                           reacstruc.process.Co      ...
                           reacstruc.process.lambda0 ...
                           reacstruc.process.tdose];
    
    % Compute dosing flow rate to get correct molar ratio
    CSdose   = reacstruc.process.CSdose;    % Nominal CSdose
    MWS      = reacstruc.model.MWS;         %  
    MWA      = reacstruc.model.MWA;         %
    lambda0  = reacstruc.process.lambda0;   % Nominal molar ratio
    Ao       = reacstruc.process.Co / MWA;  % Nominal CA at t=0
    tdose    = reacstruc.process.tdose;     % [min]     Dosing time
    
    % Compute dosing flow rate of sidechain S with nominal values
    Qdose = (tdose>0)*Ao*lambda0 / (eps+tdose*(CSdose / MWS)); %[start vols]/min
    reacstruc.instantlab.bool_Qdose = true;
                 
    % Generate normally distributed noise
%     x_T       = normrnd(0, exp_plan(k,1) / (eps + CV_input(1)) );
%     x_pH      = normrnd(0, exp_plan(k,2) / (eps + CV_input(2)) );
%     x_Co      = normrnd(0, exp_plan(k,3) / (eps + CV_input(3)) );
%     x_CSdose  = normrnd(0, CSdose / (eps + CV_input(4)) );
    x_T      = normrnd(0, sigma_input.T);
    x_pH     = normrnd(0, sigma_input.pH);
    x_Co     = normrnd(0, sigma_input.Co);
    x_CSdose = normrnd(0, sigma_input.CSdose);
    
    % Initialize simulation scenario 
    reacstruc = reacstruccreate();
    reacstruc.process.Qdose = Qdose;    % Qdose with nominal values

    % ID Description Temp pH startconc (g/L) lamda tdose
    reacstruc.process.T       = exp_plan(k,1) + x_T; % C
    reacstruc.process.pH      = exp_plan(k,2) + x_pH; %
    reacstruc.process.Co      = exp_plan(k,3) + x_Co; % g/L
    reacstruc.process.CSdose  = reacstruc.process.CSdose + x_CSdose;      
    reacstruc.process.tdose   = exp_plan(k,5); % h
    
    % Create helper variables
    CSdose = reacstruc.process.CSdose;  % Nominal CSdose
    tdose  = reacstruc.process.tdose;   
    Co     = reacstruc.process.Co / MWA;      
    
    % Compute molar ratio with uncertainties in CSdose and Co
    reacstruc.process.lambda0 = (Qdose * (CSdose/MWS) * tdose) / Co; % mol HEP/mol N9
    
    % Save actual input (with noise)
    data.act_input(k,:) = [reacstruc.process.T          ...
                           reacstruc.process.pH         ...
                           reacstruc.process.Co         ...
                           reacstruc.process.lambda0    ...
                           reacstruc.process.tdose];
    
    % Run instant lab simulation and save solutions in structure array
    reacstruc = reacsim(reacstruc);       
    
    if (bool_debug)
        figure; hold all
        title('Time profiles')

        % Add plots to figure
        plot(reacstruc.out.t, 100 * reacstruc.out.y(:,3),'-b')
        plot(reacstruc.out.t, 100 * sum(reacstruc.out.y(:,4:6),2),'-g')
        plot(reacstruc.out.t, 100 * sum(reacstruc.out.y(:,7:9),2),'r-')
        plot(reacstruc.out.t, 100 * reacstruc.out.y(:,10),'-m')
    end
              
    % Unwrap from solution vectors
    tend = reacstruc.out.t(end);
    yA   = reacstruc.out.y(end,3);  % Component A
    yB   = reacstruc.out.y(end,4);  % Component B
    yC   = reacstruc.out.y(end,5);  % Component C
    yE   = reacstruc.out.y(end,6);  % Component E
    yD   = reacstruc.out.y(end,7);  % Component D, product
    yF   = reacstruc.out.y(end,8);  % Component F
    yG   = reacstruc.out.y(end,9);  % Component G
    yH   = reacstruc.out.y(end,10); % Component H, impurity   

    ymono   = sum(reacstruc.out.y(end,4:6),2);  % Monoacylated
    ydi     = sum(reacstruc.out.y(end,7:9),2);  % Diacylated
    ytri    = reacstruc.out.y(end,10);          % Triacylated
    
    data.out{k} = [tend yA ymono ydi ytri yD yB yC yD yE yF yG yH];  
end

% Save data 
% save(folder, 'data')
end