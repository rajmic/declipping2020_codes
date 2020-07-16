function threshfunc = gl_threshfunc(W, lambda, type, arg1, arg2)

[M, N] = size(W);
switch type
    case 'time' % groups correspond to the coeffs' time labels 
        W = repmat(sqrt(sum(W.^2, 1)), M, 1);
    case 'freq' % groups correspond to the coeffs' freq. labels 
        W = repmat(sqrt(sum(W.^2, 2)), 1, N);
end
% compute threshold function
threshfunc = lambda./(W);

%{
if isnumeric(x) == 1 % STATIONARY CASE
    [M, N] = size(x);
    weights = zeros(M,N); % compute matrix of group energy weights
    switch type
        case 'time' % groups correspond to the coeffs' time labels 
            weights = repmat(sqrt(sum(abs(x).^2, 1)), M, 1);
        case 'exptime' % groups grow exponentially with frequency
            K = arg1; % number of groups per time index
            bin = round([M./((2*ones(K+1,1)).^([0:K]')); 0]);
            bin = M - bin;
            for i = 2:length(bin)
                weights(bin(i-1)+1:bin(i),:) = repmat(sqrt(sum(abs( ... 
                    x(bin(i-1)+1: bin(i),:)).^2, 1)), bin(i)-bin(i-1), 1);
            end
        case 'freq' % groups correspond to the coeffs' freq. labels 
            weights = repmat(sqrt(sum(abs(x).^2, 2)), 1, N);
        case 'onsetfreq' % group corr. to coeffs freq labels from onset to onset
            onsets = arg1;
            a = arg2;
            on = [floor(onsets/a)+1; N];
            for i = 2:length(on)
                weights(:, on(i-1):on(i)) = repmat(sqrt(sum(abs( ... 
                    x(:, on(i-1):on(i))).^2, 2)), 1, on(i)-on(i-1) +1);
            end
     end
    % compute shrinkage
    shrink = 1 - lambda./(weights);
    shrink(isinf(shrink))=0;
    shrink = shrink .* (shrink>0);
    xs = x.*shrink;
elseif iscell(x) == 1 % NON-STATIONARY CASE
    for i = 1:length(x)
        weight = norm(x{i},2); % compute matrix of group energy weights
        shrink = 1 - lambda./(weight); % compute shrinkage
        shrink(isinf(shrink))=0;
        shrink = shrink .* (shrink>0);
        xs{i,1} = x{i}.*shrink;
    end
end
%}
end
