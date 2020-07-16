function [data_rec, sdr_iter, obj_val] = aspade(data_clipped, param, paramsolver, data_orig)
% ASPADE computes the analysis version of declipping algorithm SPADE.
%
% Input parameters
%       data_clipped   vector of clipped signal
%       param          structure of parameters containing:
%                            Ls         length of the original signal
%                            theta      clipping threshold
%                            w          window length (in samples)
%                            a          window shift (in samples)
%                            wtype      string defining the type of the window, see http://ltfat.github.io/doc/sigproc/firwin.html
%                            F          frame (usually DFT frame)
%                            algorithm  string used to select A-SPADE or S-SPADE
%                            masks      structure containing clipping masks
%       paramsolver    structure of parameters for SPADE containing:
%                            s          relaxation stepsize
%                            r          relaxation steprate
%                            epsilon    stopping threshold of termination function
%                            maxit      maximal possible number of iterations with particular settings
%                            store_sdr  switch to enable computing SNR in each iteration
%                            store_obj  switch to enable storing the value of termination function in each iteration
%       data_orig      vector of original clean signal used to compute SNR during iterations
%
% Output parameters
%       data_rec       vector of restored (declipped) block of signal
%       sdr_iter       SDR values during iterations the blocks of signal
%       obj_val        Termination function values during iterations for the block of signal
%
%
%
%   x_hat (reconstructed signal - waveform)
%   z_bar (signal coefficients after hard thresholding)
%   u (dual variable, vector of signal coefficients)
%   k (required sparsity)
%   zEst = A.x_hat
%
% Pavel Záviška, Brno University of Technology, 2020



% initialization of variables
x_hat = data_clipped;
zEst = frana(param.F, x_hat);  % zEst = A.x_hat
u = zeros(length(zEst), 1);
k = paramsolver.s;
cnt = 1;
bestObj = Inf;

% initialization of vectors for SDR and termination function during iterations
sdr_iter = NaN(paramsolver.maxit, 1);
obj_val = NaN(paramsolver.maxit, 1);


while cnt <= paramsolver.maxit
    
    %set all but k largest coefficients to zero
    %(complex conjugate pairs are taken into consideration)
    z_bar = hard_thresholding(zEst + u, k);
    
    objVal = norm(zEst - z_bar); % update termination function
    
    if objVal <= bestObj
        data_rec = x_hat;
        bestObj = objVal;
    end
    
    % termination step
    if objVal <= paramsolver.epsilon
        break
    end
    
    % projection onto the set of feasible solutions
    b = z_bar - u;
    syn = frsyn(param.F, b);
    syn = syn(1:param.Ls);
    x_hat = proj_time(syn, param.masks, data_clipped);
    
    % computation and storing of SDR and termination function if requested
    if paramsolver.store_sdr
        sdr_iter(cnt) = sdr(data_orig,  x_hat);
    end
    if paramsolver.store_obj
        obj_val(cnt) = objVal;
    end
    
    % dual variable update
    zEst = frana(param.F, x_hat);
    u = u + zEst - z_bar;
    
    % iteration counter update
    cnt = cnt+1;
    
    % incrementation of variable k (require less sparse signal in next iteration)
    if mod(cnt,paramsolver.r) == 0
        k = k + paramsolver.s;
    end
    
end

end