function reacstruc = reacsim(reacstruc)

%% Initialization
% Add relative path to functions folder
% addpath(strcat(pwd, '\Functions'))

%% Correlations?
% Unwrap reference selectivity from structure array
sref    = reacstruc.model.sref;

% Update structure array
reacstruc.model.sref(1) = (1 - sref(3)) / 2; %sW1 is weight between k1 and k2
reacstruc.model.sref(2) = reacstruc.model.sref(1);

% Update helper variable
sref = reacstruc.model.sref;

%%
% Update structure array
reacstruc.model.pkA(2) = reacstruc.model.pkA(1);    % Same pKa for 1&2 (Lys)
reacstruc.model.EA(2) = reacstruc.model.EA(1);      % Same Ea for 1&2  (Lys)

% Unwrap parameters from structure array
EA      = [reacstruc.model.EA1 reacstruc.model.EA2 reacstruc.model.EA3];
MWA     = reacstruc.model.MWA;
MWS     = reacstruc.model.MWS;
Co      = reacstruc.process.Co;
CSdose  = reacstruc.process.CSdose;
textra  = reacstruc.process.textra;
lambda0 = reacstruc.process.lambda0;
tdose   = reacstruc.process.tdose;
pH      = reacstruc.process.pH;

% process
Ao      = Co / MWA; %M
Qdose   = (tdose>0) * Ao * lambda0 / (eps + tdose * (CSdose / MWS));  %[start vols]/min
y0      = [lambda0*(tdose==0); 10^(pH-14); 1; 0; 0; 0; 0; 0; 0; 0; 1]; %y(2) is [OH-], unscaled, y(11) is scaled volume
tspan   = [0 tdose+textra];

%%
% Unwrap variables from structure array
T       = reacstruc.process.T;
kAref   = reacstruc.model.kAref;
sdref   = reacstruc.model.sdref;
EAd     = reacstruc.model.EAd;
Tref    = reacstruc.model.Tref;
pkA     = [reacstruc.model.pkA1; reacstruc.model.pkA2; reacstruc.model.pkA3];
pHref   = reacstruc.model.pHref;
nd      = reacstruc.model.nd;

% Temperature, kref-->k0(T)
k0  = kAref * sref .* exp(-EA * (1 / (T+273) - 1 / (Tref + 273))); 
kd0 = kAref * sdref * exp(-EAd * (1 / (T+273) - 1 / (Tref + 273))); 

% Afprotoniseringsgrad (pH), k0-->k(pH)
k     = k0 .* [(1 + 10.^(pkA - pHref)) ./ (1 + 10.^(pkA - pH))]';
OHref = 10^(pHref - 14);
OH    = 10^(pH - 14);
kd    = kd0 * (OH / OHref)^nd;

% Save parameters in structure array 
reacstruc.model.k       = k;
reacstruc.model.kd      = kd;%
reacstruc.process.y0    = y0;
reacstruc.process.tspan = tspan;
reacstruc.process.Qdose = Qdose;
reacstruc.process.Ao    = Ao;

%% Computations
% Set ODE solver options
opts = reacstruc.solver.opts;

% Compute solutions 
[t, y] = ode15s(@reacdt, tspan,y0,opts,reacstruc);

% Save solutions 
reacstruc.out.t = t;
reacstruc.out.y = y;
reacstruc.out.yreac = sum(y(end,4:6),2) + 2*sum(y(end,7:9),2) + ...
                      3*sum(y(end,10),2);
reacstruc.out.yield = sum(y(end,7:9));


