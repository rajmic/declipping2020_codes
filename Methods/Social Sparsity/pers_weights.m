function [W] = pers_weights(G, neigh, center)
% pers_weights smoothes the data for generating persistent shrinkage masks

% get data
%kk = 2; % neighborhood norm parameter - wel = 1, ow = 2;
[M, N] = size(G);
[MM, NN] = size(neigh);
c1 = center(1); c2 = center(2);
cc1 = ceil(MM/2); cc2 = ceil(NN/2); % take actual centers of neighborhood
neigh = neigh/norm(neigh(:),1); % normalize neighborhood weights

% check neighborhood - center point relations
if center(1) > MM 
    disp('This cannot be right... Check your center point.')
    return
elseif center(2) > NN
    disp('This cannot be right... Check your center point.')
    return
end

X2 = zeros(M+MM-1, N+NN-1); % make larger matrix with mirrowed borders:
X2(c1: M+c1-1, c2: N+c2-1) = abs(G).^2; % put x in the middle of x2
X2(:, 1:c2-1) =  fliplr(X2(:, c2 : 2*(c2-1))); % left border
X2(:, N+c2:end) = fliplr(X2(:, N - NN +2*c2: N+c2-1));% right border
X2(1:c1-1, :) = flipud(X2(c1:2*(c1-1), :)); % upper border
X2(M+c1:end, :) = flipud(X2(M-MM + 2*c1: M+c1-1, :)); % lower border

% compute the stuff via convolution!
X2 = (conv2(X2, neigh, 'same'));
W = X2(cc1 : M + cc1 - 1, cc2 : N + cc2 -1); % resize matrix 
W = sqrt(W); % should be included for WGL!
end


