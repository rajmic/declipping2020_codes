function weights = weights_ath(fmin, fmax, coefCount, F, Ls, weights_assingment)
% WEIGTHS_ATH function constructs the weighting vector for the needs of the
% l1-minimization based audio declipping using the Absolute Hearing
% Threshold (ATH) curve
%
% Pavel Záviška, Brno University of Technology, 2020


% settings
normalize = 1;         % vector normalization 
tau = 100;             % limitation of the curve's dynamic range

% construction of the ATH curve
f = linspace(fmin,fmax,coefCount);  
t = (3.64.*(f/1000).^(-0.8)-6.5.*exp(-0.6.*((f/1000)-3.3).^2)+10.^(-3).*(f/1000).^4); % ATH formula
t(t>tau) = tau; 

% weights assignment
switch weights_assingment
    case 1
        weights = 1./(t-min(t)+1);
    case 2
        weights = -t + tau;
    case 3
        weights = 2*10.^(-5)*10.^((-t+tau)/20);
    otherwise
        error('Invalid type of weights assingment! Select 1, 2 or 3.')
end

% normalization
if normalize == 1
    weights = weights/max(abs(weights));
end

% copying curves into one vector that fits DGT coefficients
Ncoef = frameclength(F, Ls);
Nrep = Ncoef/coefCount;
weights = repmat(weights', Nrep, 1);

end

