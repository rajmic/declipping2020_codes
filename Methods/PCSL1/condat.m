function [data_rec, dsdr_rec, obj_func] = condat(data_clipped, param, paramsolver, data)
% CONDAT is the implementation of Condat algorithm suited for the needs of
% audio declipping based on the algorithm (P)CSL1 by Defraene et al.
%
% Pavel Záviška, Brno University of Technology, 2020


T = paramsolver.tau;
% computation of weights
if any(strcmp(param.algorithm, {'pcsl1', 'PCSL1'}))
    T = T .* param.gmt_weights;
elseif any(strcmp(param.algorithm, {'pwcsl1', 'PWCSL1'}))
    T = T .* param.parabola_weights;
end


% starting point
z = frana(param.F, data_clipped);
u_hl = zeros(length(data_clipped),1);
u_hl(param.masks.Mr) = data_clipped(param.masks.Mr);

% definition of soft thresholding
soft = @(z, T) sign(z).*max(abs(z)-T, 0);

% dsdr process
dsdr_rec = NaN(paramsolver.maxit,1);

% objective function process
obj_func = NaN(paramsolver.maxit, 1);

% iteration counter
cnt = 1;


while cnt <= paramsolver.maxit
    
    % soft thresholding step
    diff_r = (1/paramsolver.gamma) .* (param.masks.Mr .* (frsyn(param.F,z) - data_clipped));
    z_tilde_new = soft(z - paramsolver.tau*frana(param.F, diff_r + u_hl), T);
    
    % weighted average
    z_new = paramsolver.rho*z_tilde_new + (1-paramsolver.rho)*z;
    
    % projection on the clipped samples
    mu = u_hl + paramsolver.sigma * frsyn(param.F, (2*z_tilde_new) - z);
    u_hl_tilde_new = mu - paramsolver.sigma * proj_hl(mu / paramsolver.sigma, param.masks, data_clipped);
    
    % weighted average
    u_hl = paramsolver.rho * u_hl_tilde_new + (1 - paramsolver.rho) * u_hl;
    
    % compute the actual dSDR value
    if paramsolver.store_dsdr
        data_rec_tmp = frsyn(param.F, z_new); % reconstructed signal
        %data_rec_tmp(param.masks.Mr) = data_clipped(param.masks.Mr); % replacing samples on reliable positions
        dsdr_rec(cnt) = sdr(data, data_rec_tmp) - sdr(data, data_clipped); % computing dSDR
    end
    
    % compute and store the 
    if paramsolver.store_obj
        if any(strcmp(param.algorithm, {'csl1', 'CSL1'}))
            obj_func(cnt) = norm(z_new, 1); % computing objective function (l1 norm of coefficients)
        elseif any(strcmp(param.algorithm, {'pwcsl1', 'PWCSL1'}))
            obj_func(cnt) = norm(param.parabola_weights .* z_new, 1); % computing objective function (weighted-l1 norm of coefficients)
        elseif any(strcmp(param.algorithm, {'pcsl1', 'PCSL1'}))
            obj_func(cnt) = norm(param.gmt_weights .* z_new, 1); % computing objective function (weighted-l1 norm of coefficients)
        end
    end    
    
    % z update for the next iteration
    z = z_new;
    
    cnt = cnt + 1;
    
end

data_rec = frsyn(param.F, z); % reconstructed signal

end