function dydt = reacdt(t,y,reacstruc)

%% 
% CHANGED on 16-03-2022

%%
% REACDT    Model dy/dt = f(t,y,reacstruc) for acylation cube model
%
% This function implements the differential equation system for acylation 
% cube model. All parameters can be found in the reacstruc structure array.
%
% Syntax: ydot = CubeModel(t,y,reacstruc)
%
% Molar ratio of components:
% y(1):     SC
% y(2):     OH-
% y(3):     A 
% y(4):     B
% y(5):     C
% y(6):     E
% y(7):     D
% y(8):     F
% y(9):     G
% y(10):    H
% y(11):    Vol

% Unwrap parameters from structure array
k       = reacstruc.model.k;            % Reaction rate coefficient
kd      = reacstruc.model.kd;           % Degradation rate coefficient
Ao      = reacstruc.process.Ao;         % Initial concentration of component A
lambda0 = reacstruc.process.lambda0;    % Molar fraction
tdose   = reacstruc.process.tdose;      % Dosage time
Qdose   = reacstruc.process.Qdose;      % Dosage flow rate

%%
%scaled reaction rates
%r1=k(1)*xb(1)*y(3)*Ao*y(1)*Ao; dC/dt in M/min
%scaled w. Ao*S=Ao*y(1)*Ao
%r1=k(1)*xb(1)*y(3);

% Monoacylated
r1 = k(1) * y(1) * y(3) * (Ao / y(11))^2; % A -> B
r2 = k(2) * y(1) * y(3) * (Ao / y(11))^2; % A -> C
r5 = k(3) * y(1) * y(3) * (Ao / y(11))^2; % A -> E

% Diacylated
r4 = k(2) * y(1) * y(4) * (Ao / y(11))^2;   % B -> D
r8 = k(3) * y(1) * y(4) * (Ao / y(11))^2;   % B -> F
r3 = k(1) * y(1) * y(5) * (Ao / y(11))^2;   % C -> D
r7 = k(3) * y(1) * y(5) * (Ao / y(11))^2;   % C -> G
r6 = k(2) * y(1) * y(6) * (Ao / y(11))^2;   % E -> G
r9 = k(1) * y(1) * y(6) * (Ao / y(11))^2;   % E -> F

%triacyleret
r10 = k(3) * y(1) * y(7) * (Ao / y(11))^2;  % D -> H
r11 = k(2) * y(1) * y(8) * (Ao / y(11))^2;  % F -> H
r12 = k(1) * y(1) * y(9) * (Ao / y(11))^2;  % G -> H

%SC degredation (scaled)
rdeg = -kd * y(1) * (Ao / y(11)); %dCsdeg/dt in M/min, y(2) er [OH-], uskaleret, NB

%%
%Balances (scaled)
r(1,1)  = -(r1+r2+r3+r4+r5+r6+r7+r8+r9+r10+r11+r12); %Reagens
r(2,1)  = 0;           % OH-
r(3,1)  = -(r1+r2+r5); % Component A
r(4,1)  = r1-r4-r8;    % Component B
r(5,1)  = r2-r3-r7;    % Component C
r(6,1)  = r5-r9-r6;    % Component E
r(7,1)  = r4+r3-r10;   % Component D
r(8,1)  = r8+r9-r11;   % Component F
r(9,1)  = r6+r7-r12;   % Component G
r(10,1) = r10+r11+r12; % Component H

%%
%dn/dt=r*V, y=n*V/nA0, Z=y(11)=V/V0
% Scale back 
dydt = r / Ao* y(11);

%corection for dosing
dydt(1) = dydt(1) + (rdeg / Ao) * y(11) + (t<tdose) * lambda0 / (tdose+eps); 
%sidste led skal IKKE ganges med y(11),da lamda er i mol
dydt(11) = Qdose * (t<tdose);