function [data_rec, dsdr_rec, obj_func] = reweighted_chambolle_pock(data_clipped, param, paramsolver, data)
% REWEIGHTED_CHAMBOLLE_POCK is the implementation of Chamboll-Pock algorithm suited for
% the needs of signal declipping with reweighting introduced in
% "Recovering a Clipped Signal in Sparseland"
%
% Pavel Záviška, Brno University of Technology, 2020


% starting point
sol.primal = data_clipped; %randi([-1000 1000], 80000,1)/1000;
sol.dual = frana(param.F, sol.primal);

% definition of clip function (result of the Fenchel-Rockafellar conjugate of soft thresholding)
clip = @(x,T) (sign(x).*min(abs(x), T));

% dsdr process
dsdr_rec = NaN(paramsolver.maxit, paramsolver.K);

% objective function process
obj_func = NaN(paramsolver.maxit, paramsolver.K);

% initialitions
data_rec = data_clipped;
reweights = 1;

% outer cycle
for k = 1:paramsolver.K
    
    % inner iteration counter
    cnt = 1;
      
    data_rec_old = data_rec;
    
    % computation of weights
    T = reweights;
    if param.weighting
        T = T.*param.weights;
    end
    
    
    % inner cycle
    while cnt <= paramsolver.maxit

        sol.dual = clip(sol.dual + paramsolver.sigma.*frana(param.F, data_rec), T);
        
        sol.primal_old = sol.primal;
        sol.primal = proj_time(sol.primal - paramsolver.zeta*postpad(frsyn(param.F, sol.dual), param.Ls), param.masks, data_clipped);
        data_rec = sol.primal + paramsolver.rho*(sol.primal - sol.primal_old);
        
        if paramsolver.comp_dsdr %|| paramsolver.dsdr_decrease_termination
            dsdr_rec(cnt, k) = sdr(data, data_rec) - sdr(data, data_clipped); % computing dSDR
        end
        
        if paramsolver.comp_obj
            obj_func(cnt, k) = norm(param.weights .* frana(param.F, data_rec), 1); % computing objective function (l1 norm of coefficients)
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
               
        cnt = cnt + 1;
        
    end
    
    if paramsolver.verbose
        fprintf('Relative difference is: %d\n', norm(data_rec_old - data_rec));
    end
    if norm(data_rec_old - data_rec) < paramsolver.delta
        break
    end    

    % computation of the re-weights
    reweights = 1./(abs(frana(param.F, data_rec)) + paramsolver.epsilon);
    
end

    data_rec = proj_time(data_rec, param.masks, data_clipped); % final projection into the constraints

end