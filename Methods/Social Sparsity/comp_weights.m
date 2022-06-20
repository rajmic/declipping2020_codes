function weights = comp_weights(fmin, fmax, coefCount, peakfreq, powercoef)
% COMP_WEIGHTS computes the 'triangular' normalized weights  When the peak frequency is the maximum frequency and
% powercoef is equal to 2, the function provides weights in the shape of
% parabola.
%
% Pavel Záviška, Brno University of Technology, 2022


if peakfreq < fmin || peakfreq > fmax
    error('Peak Frequency must be in the bounds of Fmin and Fmax!')
end

f = linspace(fmin, fmax, coefCount);
[~, peakindex] = min((abs(f - peakfreq)));

weights = [linspace(0, 1, peakindex)'; linspace(1, 0, coefCount - peakindex)'].^powercoef;

end