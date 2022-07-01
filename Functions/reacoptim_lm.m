function f = reacoptim_lm(x, params)

% Unwrap from structure array
for i = 1:length(x)
    % Scale back from normalization
    z(i) = x(i) .* (params.UB(i) - params.LB(i)) +  params.LB(i);
end

z = z(:)';
lambda0 = z(4);

% Compute objective function
yield = predict(params.ModelYield, z);
f = (1 + params.price * lambda0) / yield;    %cost/cost precursor

end
