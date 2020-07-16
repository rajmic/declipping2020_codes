function x_ht = hard_thresholding(x,t,type)
% keep the K highest elements out of x


if strcmp(type,'constraint') % keeps t highest components
    if t == length(x)
        x_ht = x;
    else
        a = sort(abs(x), 'descend');
        x_ht = x.*(abs(x)>a(t+1)); % keep only t highest components
        
    end
elseif strcmp(type,'regularizer') % thresholds components lower than t
        x_ht = x .* (abs(x)>t);
else
    error('type unknown')
end

return