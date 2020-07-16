function [data_rec_fin, sdr_iter, obj_iter] = csl1_segmentation(data_clipped, param, paramsolver, data_orig)
% PCSL1_SEGMENTATION serves as a middle layer between main_pcsl1.m and
% Condat algorithm used for the solution of the optimization problem.
% This fucntion performs input signal padding, computation of analysis and
% synthesis window, windowing the signal and after processing each block by
% selected variant of spade algorithm it foldes the blocks back together.
%
% Input parameters:
%       data_clipped   vector of clipped signal
%
%       param          structure of parameters containing:
%                            Ls         length of the original signal
%                            theta      clipping threshold
%                            w          window length (in samples)
%                            a          window shift (in samples)
%                            wtype      string defining the type of the window, see http://ltfat.github.io/doc/sigproc/firwin.html
%                            F          frame (usually DFT frame)
%                            algorithm  string used to select CSL1 or PCSL1
%                            masks      structure containing clipping masks
%                            fs         sampling frequency
%
%       paramsolver    structure of parameters for Condat algorithm containing:
%                            gamma      regularization parameter
%                            sigma      internal parameter of the optimiation algorithm
%                            tau        internal parameter; here the threshold for soft thresholding
%                            rho        internal parameter; here the step size of the extrapolation step
%
%       data_orig      vector of original clean signal used to compute SNR during iterations
%
%
% Output parameters:
%       data_rec_fin   vector of reconstructed (declipped) signal after OLA
%       sdr_iter       matrix of SDR values during iterations for all blocks of signal
%       obj_iter       matrix of termination function values during iterations for all blocks of signal
%
% Pavel Záviška, Brno University of Technology, 2020

% preparation for padding at the end of the signal
%L = dgtlength(param.Ls,param.a,param.w); % L is divisible by a and w
L = ceil(param.Ls/param.a)*param.a+(ceil(param.w/param.a)-1)*param.a; % L is divisible by a and minimum amount of zeros equals gl (window length). Zeros will be appended to avoid periodization of nonzero samples.
N = L/param.a; % number of signal blocks

% padding the signals and masks to length L
padding = zeros(L-param.Ls, 1);
data_clipped = [data_clipped; padding];
data_orig = [data_orig; padding];
param.masks.Mr = [param.masks.Mr; true(L-param.Ls,1)];
param.masks.Mh = [param.masks.Mh; false(L-param.Ls,1)];
param.masks.Ml = [param.masks.Ml; false(L-param.Ls,1)];

% construction of analysis and synthesis windows
g = gabwin(param.wtype, param.a, param.w, L);
gana = normalize(g,'peak');      % peak-normalization of the analysis window
gsyn = gabdual(gana, param.a, param.w)*param.w;  % computing the synthesis window

% this is substituting fftshift (computing indexes to swap left and right half of the windows)
idxrange = [0:ceil(param.w/2)-1,-floor(param.w/2):-1];
idxrange2 = idxrange+abs(min(idxrange))+1;

% computing the parabola weights in the case of PWCSL1
if any(strcmp(param.algorithm, {'pwcsl1', 'PWCSL1'}))
    param.parabola_weights = linspace(0, 1, floor(param.w*param.F.redundancy/2) + 1 )'.^2;
    if mod(param.w*param.F.redundancy, 2) == 0
        param.parabola_weights = [param.parabola_weights; param.parabola_weights(end-1:-1:2)];
    else
        param.parabola_weights = [param.parabola_weights; param.parabola_weights(end:-1:2)];
    end
end

% initialization of param_seg (parameters for one signal block)
param_seg = param;
param_seg.Ls = param.w;
param_seg.masks.Mr = true(param.w,1);
param_seg.masks.Mh = false(param.w,1);
param_seg.masks.Ml = false(param.w,1);

% initalization of the gmt_weights vector
param_seg.gmt_weights = 1;

% initialization of signal blocks
data_block = zeros(param.w,1);
data_orig_block = zeros(param.w,1);
data_rec_fin = zeros(L,1);

% initialization of matrices for recording SNR and for termination function during iterations
sdr_iter = NaN(paramsolver.maxit,N);
obj_iter = NaN(paramsolver.maxit,N);

for n=0:N-1
    % multiplying signal block with windows and choosing corresponding masks
    idx = mod(n*param.a + idxrange,L) + 1;
    data_block(idxrange2) = data_clipped(idx).*gana;
    data_orig_block(idxrange2) = data_orig(idx).*gana;
    param_seg.masks.Mr(idxrange2) = param.masks.Mr(idx);
    param_seg.masks.Mh(idxrange2) = param.masks.Mh(idx);
    param_seg.masks.Ml(idxrange2) = param.masks.Ml(idx);
     
    % compute the optimization problem using Condat algorithm
    switch param.algorithm
        case {'CSL1', 'csl1'}
            [data_rec_block, sdr_iter(:,n+1), obj_iter(:,n+1)] = condat(data_block, param_seg, paramsolver, data_orig_block);
        
        case {'PWCSL1', 'pwcsl1'}
            [data_rec_block, sdr_iter(:,n+1), obj_iter(:,n+1)] = condat(data_block, param_seg, paramsolver, data_orig_block);
        
        case {'PCSL1', 'pcsl1'}
            [data_rec_block, sdr_iter(:,n+1), obj_iter(:,n+1)] = condat(data_block, param_seg, paramsolver, data_orig_block);

            [gmt, ~] = masking2_ww([data_rec_block; zeros(param.w * (param.F.redundancy - 1), 1)], param.fs);
            gmt = gmt.';
            gmt(gmt>68) = 68;
            
            if mod(param.w*param.F.redundancy, 2) == 0
                gmt = [gmt; gmt(end-1:-1:2)];
            else
                gmt = [gmt; gmt(end:-1:2)];
            end
            
            param_seg.gmt_weights = 1./(gmt-min(gmt)+1);
            param_seg.gmt_weights = param_seg.gmt_weights./max(abs(param_seg.gmt_weights));
            
        otherwise
            error('Invalid algorithm!')
    end
    
    % Folding blocks together using Overlap-Add approach (OLA)
    data_rec_block = ifftshift(data_rec_block);
    data_rec_fin(idx) = data_rec_fin(idx) + data_rec_block.*gsyn;
    
    if paramsolver.verbose
        fprintf('  Processed signal blocks: %d/%d \n', n+1, N)
    end
end

% crop the padding of reconstructed signal
data_rec_fin = data_rec_fin(1:param.Ls);
