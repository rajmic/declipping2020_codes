function y = inpaintFrame_consOMP(problemData,param)
% Inpainting method based on OMP with a constraint
% on the amplitude of the reconstructed samples an optional constraint
% on the maximum value of the clipped samples
%
% Usage: y = inpaintFrame_consOMP(problemData,param)
%
%
% Inputs:
%          - problemData.x: observed signal to be inpainted
%          - problemData.Imiss: Indices of clean samples
%          - param.D - the dictionary matrix (optional if param.D_fun is set)
%          - param.D_fun - a function handle that generates the dictionary 
%          matrix param.D if param.D is not given. See, e.g., DCT_Dictionary.m and Gabor_Dictionary.m
%          - param.wa - Analysis window
%          - param.Upper_Limit - if present and non-empty this fiels
%          indicates that an upper limit constraint is active and its
%          integer value is such that
%
% Outputs:
%          - y: estimated frame
%
% Note that the CVX library is needed.
%
% -------------------
%
% Audio Inpainting toolbox
% Date: June 28, 2011
% By Valentin Emiya, Amir Adler, Michael Elad, Maria Jafari
% This code is distributed under the terms of the GNU Public License version 3 (http://www.gnu.org/licenses/gpl.txt).
% ========================================================

%% Load data and parameters

x = problemData.x;
IObs = find(~problemData.IMiss);
p.N = length(x);
E2 = param.OMPerr^2;
E2M=E2*length(IObs);
wa = param.wa(param.N);

% build the dictionary matrix if only the dictionary generation function is given
if ~isfield(param,'D')
    param.D = param.D_fun(param);
end


% clipping level detection
clippingLevelEst = max(abs(x(:)./wa(:)));

IMiss = true(length(x),1);
IMiss(IObs) = false;
IMissPos = find(x>=0 & IMiss);
IMissNeg = find(x<0 & IMiss);

DictPos=param.D(IMissPos,:);
DictNeg=param.D(IMissNeg,:);

% Clipping level: take the analysis window into account
wa_pos = wa(IMissPos);
wa_neg = wa(IMissNeg);
b_ineq_pos = wa_pos(:)*clippingLevelEst;
b_ineq_neg = -wa_neg(:)*clippingLevelEst;
if isfield(param,'Upper_Limit') && ~isempty(param.Upper_Limit)
    b_ineq_pos_upper_limit = wa_pos(:)*param.Upper_Limit*clippingLevelEst;
    b_ineq_neg_upper_limit = -wa_neg(:)*param.Upper_Limit*clippingLevelEst;
else
    b_ineq_pos_upper_limit = Inf;
    b_ineq_neg_upper_limit = -Inf;
end

%%
Dict=param.D(IObs,:);
W=1./sqrt(diag(Dict'*Dict));
Dict=Dict*diag(W);
xObs=x(IObs);

residual=xObs;
maxNumCoef = param.sparsityDegree;
indx = [];
currResNorm2 = E2M*2; % set a value above the threshold in order to have/force at least one loop executed
j = 0;
a = [];
while currResNorm2>E2M && j < maxNumCoef,
    j = j+1;
    proj=Dict'*residual;
    [dum pos] = max(abs(proj));
    indx(j)=pos;
    a=pinv(Dict(:,indx(1:j)))*xObs;
    residual=xObs-Dict(:,indx(1:j))*a;
    currResNorm2=sum(residual.^2);
end;

if ~isempty(a)

    if isinf(b_ineq_pos_upper_limit)
        %% CVX code
        cvx_begin
        cvx_quiet(true)
        variable a(j)
        %minimize( sum(square(xObs-Dict*a)))
        minimize(norm(Dict(:,indx)*a-xObs))
        subject to
        DictPos(:,indx)*(W(indx).*a) >= b_ineq_pos
        DictNeg(:,indx)*(W(indx).*a) <= b_ineq_neg
        cvx_end
        fprintf('cvx_optval: %.3f\n', cvx_optval);
        if cvx_optval>1e3
            cvx_begin
            cvx_quiet(true)
            variable a(j)
            minimize(norm(Dict(:,indx)*a-xObs))
            cvx_end
        end
    else
        %% CVX code
        cvx_begin
        cvx_quiet(true)
        variable a(j)
        %minimize( sum(square(xObs-Dict*a)))
        minimize(norm(Dict(:,indx)*a-xObs))
        subject to
        DictPos(:,indx)*(W(indx).*a) >= b_ineq_pos
        DictNeg(:,indx)*(W(indx).*a) <= b_ineq_neg
        DictPos(:,indx)*(W(indx).*a) <= b_ineq_pos_upper_limit
        DictNeg(:,indx)*(W(indx).*a) >= b_ineq_neg_upper_limit
        cvx_end
        
        if cvx_optval>1e3
            cvx_begin
            cvx_quiet(true)
            variable a(j)
            minimize(norm(Dict(:,indx)*a-xObs))
            cvx_end
        end
    end

    %% Frame Reconstruction
    indx(length(a)+1:end) = [];

    Coeff = sparse(size(param.D,2),1);
    if (~isempty(indx))
        Coeff(indx) = a;
        Coeff = W.*Coeff;
    end
    y = param.D*Coeff;

else
    
    y = zeros(size(param.D,1),1);
    
end
    
    
return
