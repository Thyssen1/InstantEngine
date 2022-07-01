function reacstruc=reacsim(reacstruc)

%% correlations
reacstruc.model.sref(1) = (1-reacstruc.model.sref(3))/2; %sW1 is weight between k1 and k2
reacstruc.model.sref(2) = 1-reacstruc.model.sref(1)-reacstruc.model.sref(3);%sref sum to 1
% reacstruc.model.pkA(2)  = reacstruc.model.pkA(1); %same pKa for 1&2 (Lys)
% reacstruc.model.EA(2)   = reacstruc.model.EA(1); %same Ea for 1&2  (Lys)
reacstruc.model.pkA(1) = reacstruc.model.pkA1;
reacstruc.model.pkA(2) = reacstruc.model.pkA(1);
reacstruc.model.pkA(3) = reacstruc.model.pkA3;

reacstruc.model.EA(1) = reacstruc.model.EA1;
reacstruc.model.EA(2) = reacstruc.model.EA(1);
reacstruc.model.EA(3) = reacstruc.model.EA3;

%% process
Ao      = reacstruc.process.Co / reacstruc.model.MWA; %M
lambda0 = reacstruc.process.lambda0;
tdose   = reacstruc.process.tdose;

% if (~isfield(reacstruc.process, 'Qdose'))
if ~reacstruc.instantlab.bool_Qdose
    % Unwrap relevant variables
    CSdose = reacstruc.process.CSdose;
    MWS    = reacstruc.model.MWS;
    
    Qdose = (tdose>0)*Ao*lambda0 / (eps+tdose*(CSdose / MWS));  %[start vols]/min
    
    % Save in structure array
    reacstruc.process.Qdose = Qdose;
end

pH      = reacstruc.process.pH;
y0      = [lambda0*(tdose==0) 10^(pH-14) 1 0 0 0 0 0 0 0 1]; %y(2) is [OH-], unscaled, y(11) is scaled volume
y0      = y0(:);
tspan   = [0 tdose+reacstruc.process.textra];

%% Temperature, kref-->k0(T)
k0 = reacstruc.model.kAref*reacstruc.model.sref.*exp(-reacstruc.model.EA*(1/(reacstruc.process.T+273)-1/(reacstruc.model.Tref+273))); 
kd0 = reacstruc.model.kAref*reacstruc.model.sdref*exp(-reacstruc.model.EAd*(1/(reacstruc.process.T+273)-1/(reacstruc.model.Tref+273))); 

%% afprotoniseringsgrad (pH), k0-->k(pH)
pkA = reacstruc.model.pkA(:);
k = k0.*[(1+10.^(pkA-reacstruc.model.pHref))./(1+10.^(pkA-pH))]';
OHref = 10^(reacstruc.model.pHref-14);
OH = 10^(pH-14);
kd = kd0*(OH/OHref)^reacstruc.model.nd;

%[.0001*k(1) 0.5*k(3)/k(1)]
% kd
%% opdatering af model*********************************************
% reacstruc.model.s=s0.*xb/(s0(1)*xb(1)+s0(2)*xb(2)+s0(3)*xb(3)); %pH dep
% reacstruc.model.sd=sd0*y0(2)/(s0(1)*xb(1)+s0(2)*xb(2)+s0(3)*xb(3)); %pH dep
reacstruc.model.k       = k;
reacstruc.model.kd      = kd;%
reacstruc.process.y0    = y0;
reacstruc.process.tspan = tspan;
% reacstruc.process.Qdose = Qdose;

y0 = y0(:);

reacstruc.process.Ao=Ao;
%% Simulering
% Set ODE solver options
opts = reacstruc.solver.opts;

% Compute solutions
[t, y] = ode15s(@reacdt, tspan, y0, opts, reacstruc);

% Save output
reacstruc.out.t     = t;
reacstruc.out.y     = y;
reacstruc.out.yreac = sum(y(end,4:6),2)+2*sum(y(end,7:9),2)+3*sum(y(end,10),2);
reacstruc.out.yield = sum(y(end,7:9));


