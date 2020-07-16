function y = inpaintFrame_janssenInterpolation_2(problemData,param)
% Frame-level inpainting method based on the linear prediction by
% Janssen.
%
% Usage: xEst = inpaintFrame_janssenInterpolation(problemData,param)
%
% Inputs:
%          - problemData.x - observed signal to be inpainted
%          - problemData.Imiss - Indices of clean samples
%          - param.p - Order of the autoregressive model used for linear
%                      prediction
%          - param.NIt - number of iterations
%
% Outputs:
%          - y: estimated frame
%
% The function is modified in such a way that if NIt is a vector of possible
% numbers of iterations, it runs for max(NIt) iterations while saving the
% results for all the posibilities set by the vector NIt.
% 
% -------------------
%
% Audio Inpainting toolbox
% Date: June 28, 2011
% By Valentin Emiya, Amir Adler, Maria Jafari
% This code is distributed under the terms of the GNU Public License version 3 (http://www.gnu.org/licenses/gpl.txt).

s = problemData.x;
N = length(s);
Im = find(problemData.IMiss);
IObs = find(~problemData.IMiss);
M = length(Im);
Im = sort(Im); Im = Im(:); % Im: indexes of missing samples
s(Im) = 0;

if nargin<2 || ~isfield(param,'GR')
   param.GR = false;
end

NIt = param.NIt;
p = param.p;

y = NaN(length(s),length(NIt));

%IAA = abs(Im*ones(1,N)-ones(M,1)*(1:N));
IAA = abs(repmat(Im,1,N)-repmat(1:N,M,1));
IAA1 = IAA<=p;
AA = zeros(size(IAA));

if param.GR
   figure;
   hold on
end

variant = 1;

for k=1:max(NIt)
    
   % Re-estimation of LPC
   aEst = lpc(s,p).';
   
   % Re-estimation of the missing samples
   b = aEst.'*hankel(aEst.',[aEst(end),zeros(1,p)]);
   AA(IAA1) = b(IAA(IAA1)+1);
   [R, flagErr] = chol(AA(:,Im));
   if flagErr
      xEst = -inv(AA(:,Im))*AA(:,IObs)*s(IObs);
   else
      xEst = -R\(R'\(AA(:,IObs)*s(IObs)));
   end
   s(Im) = xEst;
   if param.GR
      e = filter(aEst,1,s);
      plot(k,10*log10(mean(e(p+1:end).^2)),'o')
   end
   
   if k == NIt(variant)
       y(:,variant) = s;
       variant = variant + 1;
   end
   
end

return
