function X_ht = block_hard_thresholding(X,t,type)

[M,Nframes] = size(X);
X_ht = zeros(M,Nframes);

if strcmp(type,'column_constraint') % keeps t highest in each column
    X_sort = sort(abs(X),'descend');
    X_ht = X.*bsxfun(@ge,abs(X),X_sort(t,:));
elseif strcmp(type,'global_constraint') % keep t highest coefficients in total
    X_ht = hard_thresholding(X(:),t,'constraint');
    X_ht = reshape(X_ht,[M, Nframes]);
elseif strcmp(type,'regularizer') % thresholds components lower than t
    X_ht = X .* (abs(X_ht)>t);
else
    error('type unknown')
end

return
