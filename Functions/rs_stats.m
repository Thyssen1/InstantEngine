function stats = rs_stats(Yval, Yfit)
    % Calculate error stats for a given pair of Yval and Yfit.
    % By Resul Al @DTU.
    % Yval  - Measured output
    % Yfit  - Fitted output

    % Compare inputs
    if ~all(size(Yval)==size(Yfit))
        error 'Yval and Yfit must be the same size'
    end

    % Compute summary statistics
    err = Yval-Yfit;                    % Errors or residuals
    SSE = sum((err).^2);                % Sum of squared errors
    SST = sum((Yval-mean(Yval)).^2);    % TSS, calculated from data only.
%     R2 = 1 - SSE/SST; % another way
    R2 = corr(Yval,Yfit).^2;            % Coefficient of determination
    MSE = mean((err).^2);               % Mean squared error
    RMSE = sqrt(MSE);                   % Root mean squared error
    AAD = mean(abs(err));               % Average Absolute Deviation, also called MAE
    RE = abs(err)./abs(Yval);           % Relative Error
    ARE = mean(RE);                     % Average Relative Error

    % Save stats in structuer array
    stats.R2   = R2;
    stats.SSE  = SSE;
    stats.SST  = SST;
    stats.MSE  = MSE;
    stats.RMSE = RMSE;
    stats.AAD  = AAD;
    stats.ARE  = ARE;
    stats.RE   = RE;
end