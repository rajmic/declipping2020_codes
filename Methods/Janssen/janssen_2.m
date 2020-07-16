function [ data_rec_fin ] = janssen_2(data_gapped, param, paramsolver)
% JANSSEN_2 performs input signal padding, computation of analysis and
% synthesis window, windowing the signal and after processing each block by
% Janssens's regression it folds the blocks back together.
%
% Input parameters
%       data_gapped    vector of gapped signal
%       param          structure of parameters containing:
%                            Ls         length of the original signal
%                            w          window length (in samples)
%                            a          window shift (in samples)
%                            wtype      string defining the type of the window, see http://ltfat.github.io/doc/sigproc/firwin.html
%                            mask       logical vector indicating the reliable samples
%       paramsolver    structure of solver parameters containing:
%                            NIt        number of iterations
%                            p          order of the model
%                            verbose    display parameter
%
% The function is built in such a way that if NIt is a vector of possible
% numbers of iterations, it runs for max(NIt) iterations while saving the
% results for all the posibilities set by the vector NIt.
%
% Output parameters
%       data_rec_fin   vector of reconstructed signal after OLA
%
% Ondøej Mokrý and Pavel Záviška, Brno University of Technology, 2020


% preparation for padding at the end of the signal
L = ceil(param.Ls/param.a)*param.a+(ceil(param.w/param.a)-1)*param.a; % L is divisible by a and minimum amount of zeros equals gl (window length). Zeros will be appended to avoid periodization of nonzero samples.
N = L/param.a; % number of signal blocks

% padding the signals and mask to length L
padding = zeros(L-param.Ls, 1);
data_gapped = [data_gapped; padding];
param.mask = [param.mask; true(L-param.Ls,1)];

% construction of analysis and synthesis windows
g = gabwin(param.wtype, param.a, param.w, L);
gana = normalize(g,'peak'); % peak-normalization of the analysis window
gsyn = gabdual(gana, param.a, param.w)*param.w;  % computing the synthesis window
gana = fftshift(gana);
gsyn = fftshift(gsyn);

% initialization of param_seg (parameters for one signal block)
param_seg = param;
param_seg.Ls = param.w;
param_seg.mask = true(param.w,1);

% initialization of restored signal
data_rec_fin = zeros(L,length(paramsolver.NIt));

for n=0:N-1    
    % multiplying signal block with windows and choosing corresponding masks
    idx = mod(n*param.a + (0:param.w-1),L) + 1;
    data_block = data_gapped(idx).*gana;
    param_seg.mask = param.mask(idx);
    
    % running the implementation by Adler et.al.
    problemData.x = data_block;
    problemData.IMiss = ~param_seg.mask;
    data_rec_block = inpaintFrame_janssenInterpolation_2(problemData,paramsolver);
    
    for m = 1:length(paramsolver.NIt)
        
        % Folding blocks together using Overlap-add approach (OLA)
        data_rec_fin(idx,m) = data_rec_fin(idx,m) + data_rec_block(:,m).*gsyn;
    
    end
    
    if paramsolver.verbose
        fprintf('  Processed signal blocks: %d/%d \n', n+1, N)
    end
end 
    
% crop the padding of reconstructed signal
data_rec_fin = data_rec_fin(1:param.Ls,:);

end
