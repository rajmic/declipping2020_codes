function weights = weights_mpeg(w, M, a, F, data_clipped, fs, weights_assingment)
% WEIGTHS_MPEG function constructs the weighting vector for the needs of the
% l1-minimization based audio declipping MPEG Psychoacoustic model
%
% Pavel Záviška, Brno University of Technology, 2020


% settings
normalize = 1;         % vector normalization
tau = 100;             % limitation of the curve's dynamic range


% f = linspace(fmin,fmax,coefCount);
N = F.L / a;
data_clipped(F.L)=0;
weights = zeros(F.clength(F.L),1);

% slice a window of a signal and compute the masking2 threshold
% rectangular window or window used in STFT?

idxrange = [0:ceil(w/2)-1,-floor(w/2):-1];
dftreal_length = floor(M/2)+1;
% idxrange2 = idxrange+abs(min(idxrange))+1;

for n=0:N-1
    idx = mod(n*a + idxrange,F.L) + 1;
    
    % normalize each window
    [t, ~] = masking2(data_clipped(idx), fs);
    t(t>tau) = tau;      % cut dynamic range
    
    % weights assignment
    switch weights_assingment
        case 1
            w_temp = 1./(t-min(t)+1);
        case 2
            w_temp = -t + tau;
        case 3
            w_temp = 2*10.^(-5)*10.^((-t+tau)/20);
        otherwise
            error('Invalid type of weights assingment! Select 1, 2 or 3.')
    end
    
    % normalization
    if normalize
        w_temp = w_temp./max(w_temp);
    end
    weights(n*dftreal_length+1 : (n+1)*dftreal_length) = w_temp;
    
end

end

