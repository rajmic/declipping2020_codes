function paramsolver = settings_csl1(algorithm)
% SETTINGS_CSL1 provides the parameter settings of the Condat algorithm
%
% Pavel Záviška, Brno University of Technology, 2020

    switch algorithm
        case {'CSL1', 'csl1'}
            paramsolver.gamma = 0.01;  % regularization parameter; default 0.1
            paramsolver.sigma = 1;     % internal parameter of the optimization algorithm; default 0.5
            paramsolver.tau = (2/((1/paramsolver.gamma)+2*paramsolver.sigma)) - 0.001;  % internal parameter; here the threshold for soft thresholding  
            paramsolver.rho = 0.99;    % internal parameter; here the step size of the extrapolation step

            paramsolver.maxit = 500;   % maximum number of iterations

  
        case {'PCSL1', 'pcsl1'}
            paramsolver.gamma = 0.01;  % regularization parameter; default 0.1
            paramsolver.sigma = 1;     % internal parameter of the optimization algorithm; default 0.5
            paramsolver.tau = (2/((1/paramsolver.gamma)+2*paramsolver.sigma)) - 0.001;  % internal parameter; here the threshold for soft thresholding  
            paramsolver.rho = 0.99;    % internal parameter; here the step size of the extrapolation step

            paramsolver.maxit = 500;   % maximum number of iterations
            
        
         case {'PWCSL1', 'pwcsl1'}
            paramsolver.gamma = 0.01;  % regularization parameter; default 0.1
            paramsolver.sigma = 1;     % internal parameter of the optimization algorithm; default 0.5
            paramsolver.tau = (2/((1/paramsolver.gamma)+2*paramsolver.sigma)) - 0.001;  % internal parameter; here the threshold for soft thresholding  
            paramsolver.rho = 0.99;    % internal parameter; here the step size of the extrapolation step

            paramsolver.maxit = 5000;  % maximum number of iterations

    end
    
end