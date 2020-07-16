function D = DCT_Dictionary(param)
% Windowed DCT dictionary
%
% Inputs:
%   - param.N: frame length
%   - param.M: number of atoms
%   - param.wa: analysis window function
%
% Output:
%   - Dictionary: D


wa = param.wa(param.N); % weigthing window
u = 0:(param.N-1); % time
k = 0:param.M-1; % frequency

D = diag(wa)*cos(pi/param.M*(u.'+1/2)*(k+1/2));

% normalisation
W = sqrt(sum(D.^2,1)); % 2-norm of each atom
D = bsxfun(@times,D,1./W);
return
