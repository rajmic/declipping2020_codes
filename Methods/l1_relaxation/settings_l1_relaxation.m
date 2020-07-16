function paramsolver = settings_l1_relaxation(algorithm, reweighting)
% SETTINGS_l1_relaxation provides the parameter settings for algorithms based on l1 relaxation
%
% Pavel Záviška, Brno University of Technology, 2020


    switch algorithm
        case {'Douglas-Rachford', 'DR'}
            if reweighting % Rl1CC
                paramsolver.K = 6;             % maximum number of outer iterations
                paramsolver.epsilon = 0.0001;    % parameter of the reweighting step
                paramsolver.delta = 0.01;       % criterion for outer cycle

                paramsolver.maxit = 1000;    % maximum number of iterations
                paramsolver.lambda = 1;      % step size for DR algorithm (default 1)
                paramsolver.gamma = 1;       % DR parameter; here threshold for soft thresholding
                
            else
                paramsolver.maxit = 3000;    % maximum number of iterations
                paramsolver.lambda = 1;      % step size for DR algorithm (default 1)
                paramsolver.gamma = 1;       % DR parameter; here threshold for soft thresholding
                
            end

        case {'Chambolle-Pock', 'CP'}
            if reweighting % Rl1CC
                paramsolver.K = 6;             % maximum number of outer iterations
                paramsolver.epsilon = 0.0001;    % parameter of the reweighting step
                paramsolver.delta = 0.01;       % stopping criterion for outer cycle

                paramsolver.maxit = 1000;    % maximum number of iterations
                paramsolver.zeta = 1;      % CP parameter; here step for projection
                paramsolver.sigma = 1/paramsolver.zeta;     % CP parameter; here threshold for soft thresholding
                paramsolver.rho = 1;    % step size for CP algorithm rho = <0,1>
            else    
                paramsolver.maxit = 3000;    % maximum number of iterations
                paramsolver.zeta = 1;      % CP parameter; here step for projection
                paramsolver.sigma = 1/paramsolver.zeta;     % CP parameter; here threshold for soft thresholding
                paramsolver.rho = 1;    % step size for CP algorithm rho = <0,1>
                
            end

    end
    
end