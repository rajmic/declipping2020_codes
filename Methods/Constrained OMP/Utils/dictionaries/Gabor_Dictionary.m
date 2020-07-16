function D = Gabor_Dictionary(param)
% Windowed Gabor dictionary. In this implementation, the dictionary matrix
% is the concatenation of a DCT (left part of the matrix) and of a DST
% (right part).
% Note that one can use this dictionary 
%  - either by constraining the simulaneous selection of cosine and sine
%  atoms with same frequency in order to implement Gabor atoms;
%  - or, without any selection constraint, by considering that the
%  dictionary is not a Gabor dictionary but the concatenation of a DCT and 
%  of a DST.
%
% Usage:
%   D = Gabor_Dictionary(param)
%
% Inputs [and default values]:
%   - param.N: frame length [256]
%   - param.redundancyFactor: redundancy factor to adjust the number of
%   frequencies [1]. The number of atoms in the dictionary equals 
%   param.N*param.redundancyFactor
%   - param.wd: weigthing window function [@wSine]
%
% Output:
%   - Dictionary: D (cosine atoms followed by sine atoms)
%
%
% -------------------
%
% Audio Inpainting toolbox
% Date: June 28, 2011
% By Valentin Emiya, Amir Adler, Maria Jafari
% This code is distributed under the terms of the GNU Public License version 3 (http://www.gnu.org/licenses/gpl.txt).

% Check and load default parameters
defaultParam.N = 256;
defaultParam.redundancyFactor = 1;
defaultParam.wd = @wSine;

if nargin<1
    param = defaultParam;
else
    names = fieldnames(defaultParam);
    for k=1:length(names)
        if ~isfield(param,names{k}) || isempty(param.(names{k}))
            param.(names{k}) = defaultParam.(names{k});
        end
    end
end
K = param.N*param.redundancyFactor; % number of atoms
wd = param.wd(param.N); % weigthing window
u = 0:(param.N-1); % time
k=0:K/2-1; % frequency
D = diag(wd)*[cos(2*pi/K*(u.'+1/2)*(k+1/2)),sin(2*pi/K*(u.'+1/2)*(k+1/2))];

% normalisation
D = D*diag(1./sqrt(diag(D'*D)));

return
