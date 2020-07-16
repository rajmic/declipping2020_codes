function [data_rec, dsdr_rec, obj_func] = douglas_rachford(data_clipped, param, paramsolver, data)
% DOUGLAS_RACHFORD is the implementation of Douglas-Rachford algorithm suited
% for the needs of signal declipping
% 
% Pavel Záviška, Brno University of Technology, 2020


T = paramsolver.gamma;
% computation of weights
if param.weighting
    T = T.*param.weights;
end

% starting point
c = frana(param.F, data_clipped);
%c = zeros(length(c), 1);

% definition of soft thresholding
soft = @(z, T) sign(z).*max(abs(z)-T, 0);

% dsdr process
dsdr_rec = NaN(paramsolver.maxit,1);

% objective function process
obj_func = NaN(paramsolver.maxit, 1);

% iteration counter
cnt = 1;

while cnt <= paramsolver.maxit
    
    c_tilde = proj_parse_frame(c, param, data_clipped);
    c = c + paramsolver.lambda*(soft(2*c_tilde-c, T)-c_tilde);
    
    if paramsolver.comp_dsdr
        data_rec_tmp = postpad(frsyn(param.F, c_tilde), param.Ls); % reconstructed signal
        data_rec_tmp(param.masks.Mr) = data_clipped(param.masks.Mr); % replacing samples on reliable positions
        dsdr_rec(cnt) = sdr(data, data_rec_tmp) - sdr(data, data_clipped); % computing dSDR
    end
    
    if paramsolver.comp_obj
        obj_func(cnt) = norm(param.weights.*c_tilde, 1); % computing objective function (l1 norm of coefficients)
    end
    
    if paramsolver.verbose
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
data_rec = postpad(frsyn(param.F, c_tilde), param.Ls); % reconstructed signal

end

