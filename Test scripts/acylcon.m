function [c, ceq] = acylcon(x, params)
% Constraints for statistical and machine learning models

for i = 1:length(x)
    % Scale back from normalization
    z(i) = x(i) .* (params.UB(i) - params.LB(i)) +  params.LB(i);
end
z = z(:)';

% Compute yield and impurity components
yield = predict(params.ModelYield, z);
tri   = predict(params.ModelImpurity, z);

% Compute inequality and equality constraints
% c   = [0-yield; yield-1; 0-tri; tri-0.01];
glb = params.glb;   % Upper constraint yield and impurity
gub = params.gub;   % Upper constraint yield and impurity

% glb = [0; 0];
% gub = [1; 0.01];

% c   = [0-yield; yield-1; 0-tri; tri-0.01];
c   = [glb(1)-yield; yield-gub(1); glb(2)-tri; tri-gub(2)];
ceq = [];

end