function [data_rec, dsdr_rec, obj_func] = chambolle_pock(data_clipped, param, paramsolver, data)
% CHAMBOLLE_POCK is the implementation of Chamboll-Pock algorithm suited for
% the needs of signal declipping
%
% Pavel Záviška, Brno University of Technology, 2020


% computation of weights
T = 1;
if param.weighting
    T = T.*param.weights;
end

% starting point
sol.primal = data_clipped; %randi([-1000 1000], 80000,1)/1000;
sol.dual = frana(param.F, sol.primal);

data_rec = data_clipped;

% definition of clip function (result of the Fenchel-Rockafellar conjugate of soft thresholding)
clip = @(x,T) (sign(x).*min(abs(x), T));

% dsdr process
dsdr_rec = NaN(paramsolver.maxit, 1);

% objective function process
obj_func = NaN(paramsolver.maxit, 1);

% iteration counter
cnt = 1;

while cnt <= paramsolver.maxit

    sol.dual = clip(sol.dual + paramsolver.sigma.*frana(param.F, data_rec), T);
    
    sol.primal_old = sol.primal;
    sol.primal = proj_time(sol.primal - paramsolver.zeta*postpad(frsyn(param.F, sol.dual), param.Ls), param.masks, data_clipped);
    data_rec = sol.primal + paramsolver.rho*(sol.primal - sol.primal_old);
    %data_rec = f1.prox(data_rec);
    
    if paramsolver.comp_dsdr 
        dsdr_rec(cnt) = sdr(data, data_rec) - sdr(data, data_clipped); % computing dSDR
    end
    
    if paramsolver.comp_obj
        obj_func(cnt) = norm(param.weights .* frana(param.F, data_rec), 1); % computing objective function (l1 norm of coefficients)
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
    
    cnt = cnt + 1;
       
end

data_rec = proj_time(data_rec, param.masks, data_clipped); % final projection into the constraints

end