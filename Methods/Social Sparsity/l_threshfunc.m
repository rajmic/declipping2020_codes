function threshfunc = l_threshfunc(W, lambda)

threshfunc = lambda./W;

%{
if isnumeric(x) == 1 % STATIONARY CASE
treshfunc = lambda./abs(x);
elseif iscell(x) == 1 % NON-STATIONARY CASE
    for i = 1:length(x)
        shrink = 1 - lambda./abs(x{i});
        shrink(isinf(shrink))=0;
        shrink = shrink .* (shrink>0);
        xs{i,1} = x{i}.*shrink;
    end        
end
%}
end
