function threshfunc = el_threshfunc(W, lambda, type)

[M, N] = size(W);
    switch type % STATIONARY CASE
       case 'time' % groups correspond to the coeffs' time labels
            dccSort = sort(W,1,'descend'); 
            matind = repmat((1:M)',[1,N]); 
            cumnorm1 = cumsum(dccSort,1); 
            diff = W - lambda.*(cumnorm1-matind.*dccSort);
            [ind] = find(diff<=0); 
            matind(ind) = 0; 
            indLok = max(matind,[],1); 
            ind = sub2ind([M,N],indLok,(1:N)); 
            norm1 = cumnorm1(ind); 
            LL = repmat(indLok,[M,1]);
            norm1 = repmat(norm1,M,1);       

        case 'freq' % groups correspond to the coeffs' freq. labels 
            dccSort = sort(W, 2, 'descend'); 
            matind = repmat((1:N),[M,1]); 
            cumnorm1 = cumsum(dccSort,2); 
            diff = W - lambda.*(cumnorm1-matind.*dccSort);
            [ind] = find(diff<=0); 
            matind(ind) = 0; 
            indLok = max(matind,[],2); 
            ind = sub2ind([M,N],(1:M)', indLok); 
            norm1 = cumnorm1(ind); 
            LL = repmat(indLok,[1,N]);
            norm1 = repmat(norm1,1,N);     
    end
    threshfunc = (norm1.*lambda)./((1+lambda*LL).*W); 

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
