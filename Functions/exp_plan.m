function plan = exp_plan(params)
%% EXP_PLAN
% Create experimental plan to be used by instantlab.m. 
% The user supplies a structure array of parameters to specify the type of 
% experimental plan, and the necessary input parameters.
% Designs:
%   'ff2': Two-level full factorial design
%       Required:
%           n:      Number of factors
%           reps:   Number of repetitions
%           lows:   Vector of low points
%           highs:  Vector of high points
%       Optional:
%           alpha:  Critical value (default: 0.05)
%           sigma:  Noise standard deviation (default: 0)
%           
%   'OATl3': Three-level one-variable-at-a-time (OAT)
%       Required:
%        a
%           
%       Optional:  
%       
%   'lhs': Latin hypercube sampling design
%       Required:
%        a
%           
%       Optional:    
%
%   'dsd': Definitive screening design based on publicly available MATLAB
%   script written by Jacob Albrecht, BMS based on the paper by 
%   Jones, B and Nachtsheim, C. (2011) 
%       Required:
%       
%       
% 

%% Main
design = params.design;

% Create experimental plan for InstantLab 
switch design
    case 'ff2' % Two-level full factorial design
        % Unwrap variables from structure array
        lows  = params.lows;
        highs = params.highs;
        
        % Two-level full factorial design
        n = length(lows);
        dFF2 = logical(ff2n(n));

        % Set up randomized experiments
        plan = [];
        temp = zeros(1, n);
        for i = 1:(2^n)
            temp(dFF2(i,:))  = highs(dFF2(i,:));
            temp(~dFF2(i,:)) = lows(~dFF2(i,:));

            plan = [plan; temp];
        end
        
    case 'ff2frac' % Two-level full factorial design
        % Unwrap variables from structure array
        n     = params.n;
        reps  = params.reps;
        lows  = params.lows;
        highs = params.highs;
        gen   = params.gen;

        % Two-level fractional factorial design
        plan = fracfact(gen);
        
        for i = 1:n
            plan(plan(:,i) == -1,i) = lows(i);
            plan(plan(:,i) == 1,i) = highs(i);
        end
        
        plan = repmat(plan, reps, 1);
        
    case 'ff3'
        % Unwrap variables from structure array
        lows  = params.lows;
        sets  = params.setpoints;
        highs = params.highs;
        
        plan = fullfact(3 * ones(1, length(sets)));
        
        for i = 1:length(sets)
            plan(plan(:,i) == 1,i)  = lows(i);
            plan(plan(:,i) == 2, i) = sets(i);
            plan(plan(:,i) == 3,i)  = highs(i);
        end
    case 'ff4'
        % Unwrap variables from structure array
        level_1 = params.level_1;
        level_2 = params.level_2;
        level_3 = params.level_3;
        level_4 = params.level_4;
        
        plan = fullfact(4 * ones(1,length(level_1)));
        
        for i = 1:length(level_1)
            plan(plan(:,i) == 1, i) = level_1(i);
            plan(plan(:,i) == 2, i) = level_2(i);
            plan(plan(:,i) == 3, i) = level_3(i);
            plan(plan(:,i) == 4, i) = level_4(i);
        end
    case 'ff5'
        % Unwrap variables from structure array
        level_1 = params.level_1;
        level_2 = params.level_2;
        level_3 = params.level_3;
        level_4 = params.level_4;
        level_5 = params.level_5;
        
        plan = fullfact(5 * ones(1, length(level_1)));
        
        for i = 1:length(level_1)
            plan(plan(:,i) == 1, i) = level_1(i);
            plan(plan(:,i) == 2, i) = level_2(i);
            plan(plan(:,i) == 3, i) = level_3(i);
            plan(plan(:,i) == 4, i) = level_4(i);
            plan(plan(:,i) == 5, i) = level_5(i);
        end  
    case 'ff6'
        % Unwrap variables from structure array
        level_1 = params.level_1;
        level_2 = params.level_2;
        level_3 = params.level_3;
        level_4 = params.level_4;
        level_5 = params.level_5;
        level_6 = params.level_6;
        
        plan = fullfact(6 * ones(1, length(level_1)));
        
        for i = 1:length(level_1)
            plan(plan(:,i) == 1, i) = level_1(i);
            plan(plan(:,i) == 2, i) = level_2(i);
            plan(plan(:,i) == 3, i) = level_3(i);
            plan(plan(:,i) == 4, i) = level_4(i);
            plan(plan(:,i) == 5, i) = level_5(i);
            plan(plan(:,i) == 6, i) = level_6(i);
        end
    case 'OFAT' % One-factor-at-a-time (OFAT) method
        % Unwrap parameters from structure array
        setpoints = params.setpoints;
        highs     = params.highs;
        lows      = params.lows;
        
        % Pre-allocate 
        plan = zeros(3 + 2*length(setpoints), length(setpoints));
        
        % Create plan
        plan(1:3,:) = repmat(setpoints, 3, 1);
        for i = 1:length(setpoints)
            plan(2*i+2, :) = setpoints;
            plan(2*i+2, i) = highs(i);

            plan(2*i+3, :) = setpoints;
            plan(2*i+3, i) = lows(i);
        end

    case 'lhs' % Latin hypercube sampling design
        % Unwrap parameters from structure array
        samples = params.samples;
        highs   = params.highs;
        lows    = params.lows;
        
        % Compute Latin Hypercube samples
        rng(params.seed) % Seed for reproducibility
        P = lhsdesign(samples, length(highs));
        plan = icdf('Uniform', P, repmat(lows, samples, 1), ...
                                  repmat(highs, samples, 1));
        
    case 'dsd'  % Definitive screening design
        % Unwrap parameters from structure array
        highs = params.highs;
        lows  = params.lows;
        
        %
        n = length(highs);
        setpoints = (highs + lows) ./ 2;
        
        plan = dsd(n, 0);
        idx_lows  = (plan == -1);
        idx_set   = (plan ==  0);
        idx_highs = (plan ==  1);
        
        for i = 1:size(plan,1)
            plan(i,idx_lows(i,:))  = lows(idx_lows(i,:));
            plan(i,idx_set(i,:))   = setpoints(idx_set(i,:));
            plan(i,idx_highs(i,:)) = highs(idx_highs(i,:));
        end
        
    case 'bbd' % Box-Behnken Design
        % Unwrap parameters from structure array
        highs     = params.highs;
        lows      = params.lows;
        centers   = (highs + lows) / 2;
        
        plan = bbdesign(length(highs));

        idx_lows   = (plan == -1);
        idx_center = (plan ==  0);
        idx_highs  = (plan ==  1);
        
        for i = 1:size(plan,1)
            plan(i,idx_lows(i,:))   = lows(idx_lows(i,:));
            plan(i,idx_center(i,:)) = centers(idx_center(i,:));
            plan(i,idx_highs(i,:))  = highs(idx_highs(i,:));
        end
  
    case 'ccd' % Central Co Design
        % Unwrap parameters from structure array
        setpoints = params.setpoints;
        highs = params.highs;
        lows = params.lows;
        
        plan = ccdesign(length(setpoints));

        idx_lows  = (plan == -1);
        idx_set   = (plan ==  0);
        idx_highs = (plan ==  1);
        
        for i = 1:size(plan,1)
            plan(i,idx_lows(i,:))  = lows(idx_lows(i,:));
            plan(i,idx_set(i,:))   = setpoints(idx_set(i,:));
            plan(i,idx_highs(i,:)) = highs(idx_highs(i,:));
        end
        
    otherwise
        plan = [];
end

end