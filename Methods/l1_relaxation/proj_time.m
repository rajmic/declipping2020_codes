function proj = proj_time(x, masks, data_clipped)
% PROJ_TIME perfoms projection of vector x onto the set of feasible
% solutions for the declipping problem in time domain.
%
% Input parameters
%       x               vector of input signal
%       masks           structure of 3 logical masks Mr, Mh and Ml
%       data_clipped    clipped signal
%
% Pavel Záviška, Brno University of Technology, 2020

proj = x;
proj(masks.Mr) = data_clipped(masks.Mr);
proj(masks.Mh) = max(x(masks.Mh), data_clipped(masks.Mh));
proj(masks.Ml) = min(x(masks.Ml), data_clipped(masks.Ml));


end