function [f,reacstruc] = reacoptim(x,reacstruc)

%% Initialization
Price=reacstruc.optim.Price;

%%
% Parameters for optimization
for i = 1:length(reacstruc.optim.var)
    % Declare helper variables
    LB(i) = reacstruc.optim.LB(i);  % Lower bound 
    UB(i) = reacstruc.optim.UB(i);  % Upper bound
    
    % Scale back and update optimization variable in structure array
    z(i) = x(i) .* (UB(i)-LB(i))+LB(i);                      % Scale back
    eval([reacstruc.optim.var{i},'=',num2str(z(i),10),';']); % Update
end
z=z(:)';

% Run simulation
reacstruc = reacsim(reacstruc);

% Unwrap  solution
t = reacstruc.out.t;
y = reacstruc.out.y;

%
ysum = 0*y(:,1:4); %allocation
ysum(:,1) = y(:,3);             % A
ysum(:,2) = sum(y(:,4:6),2);    % Monoacylated
ysum(:,3) = sum(y(:,7:9),2);    % Diacylated
ysum(:,4) = y(:,10);            % Triacylated

% Compute objective function
yield = ysum(end,3);
f = (1 + Price.SCrel*reacstruc.process.lambda0) / yield;    %cost/cost precursor

