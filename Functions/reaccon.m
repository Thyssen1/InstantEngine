function [c, ceq] = reaccon(x, params, reacstruc)

% reacstruc.process.T       = x(1); % C
% reacstruc.process.pH      = x(2); %
% reacstruc.process.Co      = x(3); % g/L
% reacstruc.process.lambda0 = x(4); % mol HEP/mol N9
% reacstruc.process.tdose   = x(5);

reacstruc.process.lambda0 = x(1); % mol HEP/mol N9
reacstruc.process.pH      = x(2); %
reacstruc.process.T       = x(3); % C
reacstruc.process.Co      = x(4); % g/L
% reacstruc.process.tdose   = x(5);
reacstruc = reacsim(reacstruc);

for i = 1:length(x)
    % Scale back from normalization
    z(i) = x(i) .* (params.UB(i) - params.LB(i)) +  params.LB(i);
end
z = z(:)';

[yield, tri] = reacsim_wrapper(z, reacstruc);

c = [0-yield; yield-1; 0-tri; tri-0.01];
ceq = [];

end