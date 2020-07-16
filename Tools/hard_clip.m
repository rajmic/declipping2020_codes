function [clipped, masks] = hard_clip (signal, t_min, t_max)  
% HARD_CLIP performs hard clipping of input signal and returns clipped 
% signal and clipping masks, which are three logical vectors Mr, Mh and Ml
% stored in structure masks.
%
% Input parameters
%       signal      clean input signal
%       t_min       lower clipping threshold
%       t_max       upper clipping threshold
%
% Pavel Záviška, Brno University of Technology, 2020


%% test input parameters
if min(signal)>=t_min && max(signal)<=t_max
    warning('Clipping range too large. No clipping will occur!')
end
if t_min >= t_max
    error('Lower clipping level must be smaller than the upper one!')
end

%% hard clipping & computing clipping masks
clipped = signal;

masks.Mh = (signal>t_max);
masks.Ml = (signal<t_min);
masks.Mr = ~(masks.Mh|masks.Ml);

clipped(masks.Mh) = t_max;
clipped(masks.Ml) = t_min;

end

