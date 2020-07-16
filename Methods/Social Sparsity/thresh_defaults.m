% THRESH_DEFAULTS.M - quickly sets up parameter structures for
% thresholding.m
%
% [settings]= thresh_defaults(type, varargin)
% 
% INPUT
% type:     choices: 's' (stationary transform), 'wmdct' (windowed
%           modulated discrete cosine transform).
% varargin: - the signal, in case type = 'wmdct'
%           - time shift const. a and number of channels M, in case type =
%           's'
%           - in case type = wmdct: number of freq channels M, sig;
% 
% OUTPUT
% settings: including the below sub structures:
% iter:   iteration   settings as structure array
% trans:  transform   "-------------------------"
% shrink: shrinkage   "-------------------------"
%
% Edited by Kai Siedenburg 08.07.2011.

function [varargout]= thresh_defaults(varargin)

input = struct(varargin{:});

if isfield(input, {'type'}) == 0 
   disp('Please specify your favorite transform.')
end
if input.type == 's'
   input.type = 'gab';
elseif input.type == 'gabor'
    input.type = 'gab';
end
    
switch input.type
   case 'wmdct'
        if isfield(input, {'sig'}) == 0 
            disp('Please provide the signal to be processed.')
        end
        if isfield(input, {'M'}) == 0
           if isfield(input, {'sr'}) == 0
               M = 512;
           elseif input.sr == 8000
                M = 96;
           elseif input.sr == 16000
                M = 192;
           elseif input.sr == 44100
                M = 512;
           end
        else
           M = input.M;
        end
        Ls = length(input.sig);
        ls = ceil(Ls/(2*M))*2*M;
        g = wilorth(M, ls);
        trans = struct('type', 'wmdct', 'g', g, 'dg', g, 'M', M);
        iter.type = 'Welcome my friend to the wonderful world of orthogonality.';
        shrink = struct('type', 'wgl', 'lambda', 0.01, 'glabel', 'time', ...
                  'neigh', ones(1,9), 'center', [1,7], 'expo', 1);
              
    case 'gab'
        if isfield(input, {'M'}) == 0
            if isfield(input, {'sr'}) == 0
                M = 1024;
            elseif input.sr == 8000
                M = 192;
            elseif input.sr == 16000
                M = 384;
            end
        else
           M = input.M;
        end
        if isfield(input, {'shift'}) == 0 
           a = ceil(M/4);
        else
           a = input.shift;
        end
        g = gabwin({'tight', 'hann'}, a, M);
        iter = struct('type', 'fast', 'maxit', 20, 'tol', 1e-3, 'gamma', 1, 'disp', 1, 'method', 'reg', 'maxSubit', 10);
        trans = struct('type', 'sreal', 'g', g, 'dg', g, 'shift', a, 'M', M);
        shrink = struct('type', 'wgl', 'lambda', 0.01, 'glabel', 'time', ...
                        'neigh', ones(1,9), 'center', [1,7], 'expo', 1);
        if shrink.expo == 2
           iter.maxit = 1;
        end
end

if nargout == 1
    settings = struct('iter', iter, 'trans', trans, 'shrink', shrink);
    varargout{1} = settings;
elseif nargout == 3
    varargout{1} = iter;
    varargout{2} = trans;
    varargout{3} = shrink;
end
end


        
