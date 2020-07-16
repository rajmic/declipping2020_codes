function [sig_rec, varargout] = thresholding(sig, varargin)

% THRESHOLDING serves as front end for the structured sparsity toolbox.
% It takes variable numbers of input and output and is incredibly user
% friendly!
%
% Input Options: 
% [sig_rec, varargout] = thresholding(sig) - sets up all settings by default
% [sig_rec, varargout] = thresholding(sig, settings) - settings must be 
% structure of structures as given out by thresholding as fourth output
% argument or the initialization function THRESH_DEFAULTS.M
% [sig_rec, varargout] = thresholding(sig, iter, trans, shrink) - iter,
% trans and shrink must be structures
%
% Output Options:
% [sig_rec] = thresholding(sig) - reconstructed signal only
% [sig_rec, Gs] = thresholding(sig) - Gs is shrunk coeff matrix
% [sig_rec, Gs, rel_err] = thresholding(sig) - relative iteration error
% [sig_rec, Gs, rel_err, settings] = thresholding(sig) - structure of
% settings


% Input:
if nargin == 1
    shrink.type = 'l'; trans.type = 's'; iter.type = 'fast';
elseif nargin == 2
    if ischar(varargin{1})
        defaultsonly = 1;
        iter = 'fast'; trans = 's'; shrink = 'fast';
    else
        iter = varargin{1}.iter; trans = varargin{1}.trans; shrink = varargin{1}.shrink;
    end
elseif nargin == 4
    iter = varargin{1}; trans = varargin{2}; shrink = varargin{3};
end
    
% DEFAULTS
% Iteration
if strcmp(trans.type, 'wmdct')
    iter = 0;
else
    if isfield(iter, {'type'}) == 0 
        iter.type = 'fast'; 
    end
    if isfield(iter, {'tol'}) == 0 
        iter.tol = 0.001; 
    end
    if isfield(iter, {'maxit'}) == 0 
        iter.maxit = 20; 
    end
    if isfield(iter, {'gamma'}) == 0
        iter.gamma = 1; % could sometime also be adjusted to the actual circumstances set by trans!
    end
    if isfield(iter, {'disp'}) == 0
        iter.disp = 0; % could sometime also be adjusted to the actual circumstances set by trans!
    end
end

% Transform
trans.Ls = length(sig);
if isfield(trans, {'type'}) == 0 
    trans.type = 's';
end
if strcmp(trans.type, 's')
    if isreal(sig)
        trans.type = 'sreal';
    end
end
if isfield(trans, {'g'}) == 0 && isfield(trans, {'M'}) == 0
    trans.M = 1024;
end
if isfield(trans, {'M'}) == 0 && isfield(trans, {'g'}) ~= 0
    trans.M = length(trans.g);
end
if isfield(trans, {'shift'}) == 0
    trans.shift = (trans.M)/4;
end
if isfield(trans, {'g'}) == 0
    trans.g = gabwin({'tight', 'hann'}, trans.shift, trans.M);
end
if isfield(trans, {'dg'}) == 0
    trans.dg = trans.g; % a little dangerous...
end

% Shrinkage
if isfield(shrink, {'type'}) == 0 
    shrink.type = 'l';
end
if isfield(shrink, {'lambda'}) == 0 
    shrink.lambda = 0.01;
end
if isfield(shrink, {'glabel'}) == 0 
    shrink.glabel = 'time';
end
if isfield(shrink, {'expo'}) == 0
    shrink.expo = 1;
end
if isfield(shrink, {'size'}) == 1 % the old format overwrites the new one...
    shrink.neigh = ones(shrink.size(1)+shrink.size(3)+1,...
        shrink.size(2)+shrink.size(4) +1);
        shrink.center = [shrink.size(1)+1, shrink.size(4)+1];
        if strcmp(shrink.type, 'l') % exept for lasso, where it should be:
            shrink.neigh = [1];
            shrink.center = [1,1];
        end
end
if isfield(shrink, {'neigh'}) == 0
    shrink.neigh = [1 1 1; 1 1 1; 1 1 1];
end
if isfield(shrink, {'center'}) == 0
    shrink.center = [2,2];
end
if strcmp(shrink.type, 'wgl')
    shrink.type = 'l';
elseif strcmp(shrink.type, 'pgl')
    shrink.type = 'gl';
elseif strcmp(shrink.type, 'pel')
    shrink.type = 'el';
end

% Display the Settings:
settings = struct('iter', iter, 'trans', trans, 'shrink', shrink);
%disp('Iteration Settings:'); disp(iter);
%disp('Transform Settings:'); disp(trans);
%disp('Shrinkage Settings:'); disp(shrink);


%if ischar(varargin{1})
%   if defaultsonly
%        sig_rec = 0; rel_err = 0; Gs = 0;
%        varargout{1} = Gs;
%        varargout{2} = rel_err;
%        varargout{3} = settings;
%   end
%   return
%end


% COMPUTATION:
G = trafo(sig, trans); % Compute the coefficients which are to be shrunk

if strcmp(trans.type, 'wmdct')
   Gs = gen_thresh(G, shrink);
   rel_err = 0; % Its orthogonal.
else
    if strcmp(iter.type, 'reg')
        [Gs, rel_err] = t_ista(G, settings);
    elseif strcmp(iter.type, 'fast')
        [Gs, rel_err] = t_fista(G, settings);
    end
end
sig_rec = itrafo(Gs, trans); % Reconstruct the signal
varargout{1} = Gs;
varargout{2} = rel_err;
varargout{3} = settings;
end

