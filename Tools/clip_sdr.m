function [clipped, masks, clipping_threshold, TrueSDR, percentage] = clip_sdr (signal, desiredSDR)
% CLIP_SDR function clip the input signal according to the desired SDR. 
%
% Pavel Záviška, Brno University of Technology, 2020

% Difference function
diffSDR = @(T) sdr(signal, hard_clip(signal, -T, T))-desiredSDR;

% Finding the optimal clipping threshold for given inputSDR
[clipping_threshold, diff_from_desiredSDR] = fzero(diffSDR,[eps 0.99*max(abs(signal))]);
TrueSDR = desiredSDR + diff_from_desiredSDR;

% Clipping the signal with the optimal clipping threshold
[clipped, masks] = hard_clip(signal, -clipping_threshold, clipping_threshold);

% Computing the the percentage of clipped samples
percentage = (sum(masks.Mh) + sum(masks.Ml)) / length(signal) * 100;

end
