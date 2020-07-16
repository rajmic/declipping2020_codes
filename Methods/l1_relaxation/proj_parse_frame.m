function [proj] = proj_parse_frame(c, param, data_clipped)
% PROJ_PARSE_FRAME perfoms projection of coefficients c
% onto a multidimensional interval for the use of audio declipping.
% 
% Input parameters
%       c               vector of input coefficients
%       param           structure containing clipping masks, clipping
%                       threshold theta and frame F
%       data_clipped    original clipped signal
%
% This projection
%       proj(z) = argmin_{u} ||z - u||_2 s.t. Au \in [b_1, b_2]
% can be evaluated as 
%       proj(z) = z-A^+(Az - proj(Az));  %here A^+ denotes pseudoinverse
% 
% The projection proj(Az) is computed by function proj_time.
% 
% Please note that this particular function works only for Parseval tight frame 
% (only in this case the pseudoinverse is identical to the analysis of the signal)
%
% Pavel Záviška, Brno University of Technology, 2020


% Synthesis of the signal 
syn = frsyn(param.F, c);

% Compute proj(Az)
proj_temp = proj_time(syn, param.masks, data_clipped);

% Final projection
proj = c - frana(param.F, syn-proj_temp);


end