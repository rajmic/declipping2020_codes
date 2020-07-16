function D = DCT_Dictionary(param)
% Windowed DCT dictionary
%
% Usage:
%   D = DCT_Dictionary(param)
%
%
% Inputs [and default values]:
%   - param.N: frame length [256]
%   - param.redundancyFactor: redundancy factor to adjust the number of
%   frequencies [1]. The number of atoms in the dictionary equals 
%   param.N*param.redundancyFactor
%   - param.wd: weigthing window function [@wSine]
%
% Output:
%   - Dictionary: D
%
%
% -------------------
%
% Audio Inpainting toolbox
% Date: June 28, 2011
% By Valentin Emiya, Amir Adler, Maria Jafari
% This code is distributed under the terms of the GNU Public License version 3 (http://www.gnu.org/licenses/gpl.txt).
% Windowed DCT dictionary
%

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
k=0:K-1; % frequency
D = diag(wd)*cos(pi/K*(u.'+1/2)*(k+1/2));

% normalisation
D = D*diag(1./sqrt(diag(D'*D)));
return
