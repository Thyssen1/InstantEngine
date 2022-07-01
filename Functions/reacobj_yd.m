function f = reacobj_yd(x,reacstruc)
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
% Res(1:nsim, 1) = 0;

%ID	Description	Temp	pH	startconc (g/L)	lamda	tdose	Csdose	SC purity

% For nsim number each simulation do 
Res = zeros(nsim*2,1);
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
    tdat = reacstruc.data.num{k}(1);
    yA   = reacstruc.data.num{k}(2);
    yB   = reacstruc.data.num{k}(7);
    yC   = reacstruc.data.num{k}(8);
    yD   = reacstruc.data.num{k}(6);
    yE   = reacstruc.data.num{k}(10);
    yF   = reacstruc.data.num{k}(11);
    yG   = reacstruc.data.num{k}(12);
    yH   = reacstruc.data.num{k}(13);

%      data.out{k} = [tend yA ymono ydi ytri yD yB yC yD yE yF yG yH]; 
    
    % Run simulation and save solutions in structure array
    reacstruc = reacsim(reacstruc);
    
    % Unwrap solutions
    t = reacstruc.out.t;
    y = reacstruc.out.y;
    
    yAout = interp1(t, y(:,3), tdat);
    yBout = interp1(t, y(:,4), tdat);
    yCout = interp1(t, y(:,5), tdat);
    yDout = interp1(t, y(:,7), tdat);
    yEout = interp1(t, y(:,6), tdat);
    yFout = interp1(t, y(:,8), tdat);
    yGout = interp1(t, y(:,9), tdat);
    yHout = interp1(t, y(:,10), tdat);
    
%     Res(k,1) = yDout - yD;
    Res(k,1)        = yAout - yA;
    Res(k+nsim,1)   = yBout - yB;
    Res(k+2*nsim,1) = yCout - yC;
    Res(k+3*nsim,1) = yDout - yD;
    Res(k+4*nsim,1) = yEout - yE;
    Res(k+5*nsim,1) = yFout - yF;
    Res(k+6*nsim,1) = yGout - yG;
    Res(k+7*nsim,1) = yHout - yH;

%     yA    = reacstruc.data.num{k}(2);
%     yMono = reacstruc.data.num{k}(3);
%     yDi   = reacstruc.data.num{k}(4);
%     yTri  = reacstruc.data.num{k}(5);
% 
%     yAout    = interp1(t, y(:,3), tdat);
%     yMonoOut = interp1(t, sum(y(:,4:6),2), tdat);
%     yDiOut   = interp1(t, sum(y(:,7:9),2), tdat);
%     yTriOut  = interp1(t, y(:,10), tdat);
% 
%     ResA    = (yAout(:) - yA)^2;
%     ResMono = (yMonoOut(:) - yMono)^2;
%     ResDi   = (yDiOut(:) - yDi)^2;
%     ResTri  = (yTriOut(:) - yTri)^2;
%     
%     Res(k,1) = sqrt(ResA + ResMono + ResDi + ResTri);
end

% Compute value of objective function
% f = sum(Res);
f = Res;