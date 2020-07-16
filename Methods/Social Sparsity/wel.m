function Gs = wel(G, shrink)

lambda = shrink.lambda;
center = shrink.center;
neigh = shrink.neigh;
[MM, NN] = size(neigh);
c1 = center(1); c2 = center(2);
[M, N] = size(G);

% check neighborhood - center point relations
if center(1) > MM 
    disp('This cannot be right... Check your center point.')
    return
elseif center(2) > NN
    disp('This cannot be right... Check your center point.')
    return
end
X2 = zeros(M+MM-1, N+NN-1);
X2(c1: M+c1-1, c2: N+c2-1) = abs(G); % put x in the middle of x2
X2(:, 1:c2-1) =  fliplr(X2(:, c2 : 2*(c2-1))); % left border
X2(:, N+c2:end) = fliplr(X2(:, N - NN +2*c2: N+c2-1)); % right border
X2(1:c1-1, :) = flipud(X2(c1:2*(c1-1), :)); % upper border
X2(M+c1:end, :) = flipud(X2(M-MM + 2*c1: M+c1-1, :)); % lower border

neigh = shrink.neigh;
neighlong = neigh(:);
P = length(neighlong);
W = zeros(M,N,P);

i = c1 - 1;
j = NN - c2;
k = MM - c1;
l = c2 -1;

for m = 1:M-1
    for n = 1:N-1
        ind1 = MM-1+m-k : MM-1+m+i;
        ind2 = NN-1+n-l : NN-1+n+j;
        mat = X2(ind1, ind2).*neigh;
        W(m,n,:) = mat(:);
    end
end

dccSort = sort(W,3,'descend'); 
Z = zeros(1,1, P); Z(1,1,:) = (1:P)'; matind = repmat(Z,[M,N,1]); 
cumnorm1 = cumsum(dccSort,3); 
diff = W - lambda.*(cumnorm1-matind.*dccSort);
[ind] = find(diff<=0); 
matind(ind) = 0; 
indLok = max(matind,[],1);
ind = sub2ind([M,N,P],indLok,(1:N)); 
norm1 = cumnorm1(ind); 
LL = repmat(indLok,[M,1]);
norm1 = repmat(norm1,M,1);       

TF = (norm1.*lambda)./((1+lambda*LL).*W); 
TF(isinf(TF)) = 0;
Gs = 1 - TF.^(shrink.expo);
Gs = G.*(Gs.*(Gs>0));


%{
if isnumeric(x) == 1
    [M, N] = size(x);
    weights = abs(x);
    switch type % STATIONARY CASE
       case 'time' % groups correspond to the coeffs' time labels
            dccSort = sort(weights,1,'descend'); 
            matind = repmat((1:M)',[1,N]); 
            cumnorm1 = cumsum(dccSort,1); 
            diff = weights - lambda.*(cumnorm1-matind.*dccSort);
            [ind] = find(diff<=0); 
            matind(ind) = 0; 
            indLok = max(matind,[],1); 
            ind = sub2ind([M,N],indLok,(1:N)); 
            norm1 = cumnorm1(ind); 
            LL = repmat(indLok,[M,1]);
            norm1 = repmat(norm1,M,1);       

        case 'freq' % groups correspond to the coeffs' freq. labels 
            dccSort = sort(weights, 2, 'descend'); 
            matind = repmat((1:N),[M,1]); 
            cumnorm1 = cumsum(dccSort,2); 
            diff = weights - lambda.*(cumnorm1-matind.*dccSort);
            [ind] = find(diff<=0); 
            matind(ind) = 0; 
            indLok = max(matind,[],2); 
            ind = sub2ind([M,N],(1:M)', indLok); 
            norm1 = cumnorm1(ind); 
            LL = repmat(indLok,[1,N]);
            norm1 = repmat(norm1,1,N);     
    end
    shrink = 1 - (norm1.*lambda)./((1+lambda*LL).*weights); 
    shrink(isinf(shrink))=0;
    shrink = shrink .* (shrink>0);
    xs = x.*shrink;
    
elseif iscell(x) == 1 % NON-STATIONARY CASE
           for i = 1:length(x)
               M = length(x{i});
               weights = abs(x{i});
               dccSort = sort(weights,1,'descend'); 
               matind = (1:M)'; 
               cumnorm1 = cumsum(dccSort,1); 
               diff = weights - lambda.*(cumnorm1-matind.*dccSort);
               [ind] = find(diff<=0); 
               matind(ind) = 0; 
               indLok = max(matind,[],1); 
               ind = sub2ind([M,1],indLok,1); 
               norm1 = cumnorm1(ind); 
               LL = repmat(indLok,[M,1]);
               norm1 = repmat(norm1,M,1);
               shrink = 1 - (norm1.*lambda)./((1+lambda*LL).*weights); 
               shrink(isinf(shrink))=0;
               shrink = shrink .* (shrink>0);
               xs{i,1} = x{i}.*shrink;
           end
end
%}
end
