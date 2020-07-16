function [Gs] = gen_thresh(G, shrink)


% computation of persistence weights
W = pers_weights(G, shrink.neigh, shrink.center);

% computation of threshold function values
switch shrink.type
    case 'l' % LASSO
        TF = l_threshfunc(W, shrink.lambda);
        
    case 'gl' % GROUP-LASSO
        shrink.lambda = shrink.lambda;
        TF = gl_threshfunc(W, shrink.lambda, shrink.glabel);
        
    case 'el' % ELITIST LASSO
        TF = el_threshfunc(W, shrink.lambda, shrink.glabel);
end

% compute shrinkage
TF(isinf(TF)) = 0;
Gs = 1 - TF.^(shrink.expo);
GGs = (Gs.*(Gs>0));

if shrink.orth==1
    [M, N] = size(GGs);
    [MM, NN] = size(shrink.neigh);
    c1 = shrink.center(1); c2 = shrink.center(2);
    cc1 = ceil(MM/2); cc2 = ceil(NN/2); % take actual centers of neighborhood
    neigh = shrink.neigh/norm(shrink.neigh(:),1); % normalize neighborhood weights
    
    X2 = zeros(M+MM-1, N+NN-1); % make larger matrix with mirrowed borders:
    X2(c1: M+c1-1, c2: N+c2-1) = GGs; % put x in the middle of x2
    X2(:, 1:c2-1) =  fliplr(X2(:, c2 : 2*(c2-1))); % left border
    X2(:, N+c2:end) = fliplr(X2(:, N - NN +2*c2: N+c2-1));% right border
    X2(1:c1-1, :) = flipud(X2(c1:2*(c1-1), :)); % upper border
    X2(M+c1:end, :) = flipud(X2(M-MM + 2*c1: M+c1-1, :)); % lower border
    
    % compute the stuff via convolution!
    X2 = (conv2(X2, neigh, 'same'));
    GGs = X2(cc1 : M + cc1 - 1, cc2 : N + cc2 -1); % resize matrix        
end

Gs = G.*GGs;


end
