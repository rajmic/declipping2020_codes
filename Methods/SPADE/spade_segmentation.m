function [data_rec_fin, sdr_iter, obj_iter]=spade_segmentation(data_clipped, param, paramsolver, data_orig)
% SPADE_SEGMENTATION serves as a middle layer between declipping_main.m and
% SPADE algorithm (aspade.m or sspade_new.m)
% This fucntion performs input signal padding, computation of analysis and
% synthesis window, windowing the signal and after processing each block by
% selected variant of spade algorithm it foldes the blocks back together.
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
%
% Output parameters
%       data_rec_fin   vector of reconstructed (declipped) signal after OLA
%       sdr_iter       matrix of SNR values during iterations for all blocks of signal
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

% initialization of param_seg (parameters for one signal block)
param_seg = param;
param_seg.Ls = param.w;
param_seg.masks.Mr = true(param.w,1);
param_seg.masks.Mh = false(param.w,1);
param_seg.masks.Ml = false(param.w,1);

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
    
    % perform SPADE
    switch param.algorithm
        case {'ASPADE', 'aspade'}
           [data_rec_block, sdr_iter(:,n+1), obj_iter(:,n+1)] = aspade(data_block, param_seg, paramsolver, data_orig_block);
        case {'SSPADE_NEW', 'sspade_new'}
           [data_rec_block, sdr_iter(:,n+1), obj_iter(:,n+1)] = sspade_new(data_block, param_seg, paramsolver, data_orig_block); 
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
