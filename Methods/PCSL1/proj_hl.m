function proj = proj_hl(x, masks, data_clipped)
% PROJ_HL perfoms projection of vector x onto the set of feasible
% solutions (defined only for the clipped samples) for the declipping 
% problem in the time domain.
%
% Input parameters
%       x               vector of input signal
%       masks           structure of 3 logical masks Mr, Mh and Ml
%       data_clipped    clipped signal
%
% Pavel Záviška, Brno University of Technology, 2020


proj = x;
proj(masks.Mh) = max(x(masks.Mh), data_clipped(masks.Mh));
proj(masks.Ml) = min(x(masks.Ml), data_clipped(masks.Ml));


end