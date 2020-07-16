function w = wSine(L)
% Symetric sine window with length L
%
% Usage:
%     w = wSine(L)
%
% Inputs:
%          - L - Window length
%
% Outputs:
%          - w - window
%
% -------------------
%
% Audio Inpainting toolbox
% Date: June 28, 2011
% By Valentin Emiya, Amir Adler, Maria Jafari
% This code is distributed under the terms of the GNU Public License version 3 (http://www.gnu.org/licenses/gpl.txt).

w = sin(((0:(L-1))+.5)/L*pi);
return
