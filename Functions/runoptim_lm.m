function sol = runoptim_lm(params, algorithm)

% Normalize initial condition for each optimization variable
for i = 1:length(params.z0)
    % Compute normalized initial condition
    x0(i) = (params.z0(i) - params.LB(i)) / (params.UB(i) - params.LB(i));

    % Set normalized bounds
    con(1,i) = 0;
    con(2,i) = 1;
end

% 
x0 = x0(:);

% Set function handle for objective function
fun = @(x) reacoptim_lm(x, params);

% Set optimizer tolerance
options = optimset('TolX',1e-9,'TolFun',1e-9);

% Set function handle for nonlinear constraints
nonlincon = @(x) acylcon(x, params);

% Compute optimum
switch algorithm
    case 'fminsearchcon'
        sol = fminsearchcon(fun, x0, con(1,:), con(2,:), [], [], nonlincon, options);
    case 'fmincon'
%         options = optimset('TolX',1e-9,'TolFun',1e-9, 'Algorithm', 'sqp');
        sol = fmincon(fun, x0, [], [], [], [], con(1,:), con(2,:), nonlincon, options);
end
% sol = fminsearchcon(fun, x0, con(1,:), con(2,:), [], [], nonlincon, options);
% sol = fmincon(fun, x0, [], [], [], [], con(1,:), con(2,:), nonlincon, options);

for i = 1:length(sol)
    % Scale back from normalization
    z_sol(i) = sol(i) .* (params.UB(i) - params.LB(i)) +  params.LB(i);
end
z_sol = z_sol(:)';

% Compute predicted yield
% yield_predicted = predict(params.mdl, z_sol);

end