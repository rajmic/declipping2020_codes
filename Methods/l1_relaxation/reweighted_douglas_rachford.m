function [data_rec, dsdr_rec, obj_func] = reweighted_douglas_rachford(data_clipped, param, paramsolver, data)
% REWEIGHTED_DOUGLAS_RACHFORD is the implementation of Douglas-Rachford algorithm
% suited for the needs of signal declipping with reweighting introduced in
% "Recovering a Clipped Signal in Sparseland"
%
% Pavel Záviška, Brno University of Technology, 2020


% definition of soft thresholding
soft = @(z, T) sign(z).*max(abs(z)-T, 0);

% dsdr process
dsdr_rec = NaN(paramsolver.maxit, paramsolver.K);

% objective function process
obj_func = NaN(paramsolver.maxit, paramsolver.K);

% starting point
c = frana(param.F, data_clipped);

% initializations
c_tilde = zeros(length(c),1);
reweights = 1;

% outer cycle
for k = 1:paramsolver.K
    
    % inner iteration counter
    cnt = 1;
       
    c_tilde_old = c_tilde;
    
    % computation of weights
    T = paramsolver.gamma .* reweights;
    if param.weighting
        T = T.*param.weights;
    end
    
    
    % inner cycle
    while cnt <= paramsolver.maxit
        c_tilde = proj_parse_frame(c, param, data_clipped);
        c = c + paramsolver.lambda*(soft(2*c_tilde-c, T)-c_tilde);
        
        if paramsolver.comp_dsdr
            data_rec_tmp = postpad(frsyn(param.F, c_tilde), param.Ls); % reconstructed signal
            data_rec_tmp(param.masks.Mr) = data_clipped(param.masks.Mr); % replacing samples on reliable positions
            dsdr_rec(cnt, k) = sdr(data, data_rec_tmp) - sdr(data, data_clipped); % computing dSDR
        end
        
        if paramsolver.comp_obj
            obj_func(cnt, k) = norm(param.weights.*c_tilde, 1); % computing objective function (l1 norm of coefficients)
        end
        
        if paramsolver.verbose
            fprintf('Outer cycle no. %u/%u', k, paramsolver.K);
            fprintf('  Iteration number: %u', cnt);
            if paramsolver.comp_dsdr
                fprintf(' -- SDR improvement: %e', dsdr_rec(cnt));
            end
            if paramsolver.comp_obj
                fprintf(' -- Objective function value: %e', obj_func(cnt));
            end
            fprintf('\n')
        end
        
        cnt = cnt+1;
        
    end
    
    c_tilde = proj_parse_frame(c, param, data_clipped); % final projection into the constraints
    
    if paramsolver.verbose
        fprintf('Relative difference is: %d\n', norm(c_tilde_old - c_tilde));
    end
    if norm(c_tilde_old - c_tilde) < paramsolver.delta
        break
    end
    
    % computation of the re-weights
    reweights = 1./(abs(c_tilde) + paramsolver.epsilon);
    
    
end

data_rec = postpad(frsyn(param.F, c_tilde), param.Ls); % reconstructed signal

end

