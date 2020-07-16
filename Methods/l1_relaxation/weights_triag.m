function weights = weights_triag(fmin, fmax, coefCount, F, Ls, peakfreq, powercoef)
% WEIGHTS_TRIAG computes the 'triangular' normalized weights fitting the DGT
% coeffcients. When the peak frequency is the maximum frequency and
% powercoef is equal to 2, the function provides weights in the shape of
% parabola.
%
% Pavel Záviška, Brno University of Technology, 2020


if peakfreq < fmin || peakfreq > fmax
    error('Peak Frequency must be in the bounds of Fmin and Fmax!')
end

f = linspace(fmin, fmax, coefCount);
[~, peakindex] = min((abs(f - peakfreq)));

weights = [linspace(0, 1, peakindex) linspace(1, 0, coefCount - peakindex)].^powercoef;

% copying curves into one vector that fits the DGT coefficients
Ncoef = frameclength(F, Ls);
Nrep = Ncoef/coefCount;

weights = repmat(weights', Nrep, 1);

end