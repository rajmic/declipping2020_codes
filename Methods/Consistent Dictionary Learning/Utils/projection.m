function X_proj = projection(X, Y, distortion_param)

% X: current estimate
% Y: (nonlinear) measurements
% distortion_param: additional parameters

if strcmp(distortion_param.name, 'clean')

    X_proj = Y;

elseif strcmp(distortion_param.name, 'inpainting')

    X_proj = X;
    mask_mat = logical(distortion_param.reliable_samples_mat);
    
    X_proj(mask_mat) = Y(mask_mat);

elseif strcmp(distortion_param.name, 'clipping')

    pos_clipped = logical(distortion_param.positive_clipped);    
    neg_clipped = logical(distortion_param.negative_clipped);

    X_proj = Y;
    X_proj(pos_clipped) = max(Y(pos_clipped),X(pos_clipped));
    X_proj(neg_clipped) = min(Y(neg_clipped),X(neg_clipped));

elseif strcmp(distortion_param.name, 'quantization')
    
    Delta = 2/2^distortion_param.Nbits;
    
    pos = X>=Y+Delta/2; % samples that violate upper constraint
    neg = X<=Y-Delta/2; % samples that violate lower constraint

    X_proj = X;
    X_proj(pos) = Y(pos) + Delta/2; % project on upper constraint
    X_proj(neg) = Y(neg) - Delta/2; % project on lower constraint
    
        
elseif strcmp(distortion_param.name, 'zero_crossing_distortion')
    
    t = distortion_param.threshold;
    missing_samples = logical(distortion_param.missing_samples_mat);
    
    X_proj = X; % missing samples in [-t,t]
    X_proj(missing_samples & X>t) = t; % missing samples above t
    X_proj(missing_samples & X<-t) = -t; % below
%     
%     X_proj(~missing_samples & Y>0) = Y(~missing_samples & Y>0)+t;    
%     X_proj(~missing_samples & Y<0) = Y(~missing_samples & Y<0)-t;
    
    X_proj(~missing_samples) = Y(~missing_samples)+t*sign(Y(~missing_samples));
    
elseif strcmp(distortion_param.name, '1Bit_compression')
    
    % Projection on feasibility set:
    X_proj = X;
    X_proj((Y.*sign(X))<=0) = 0; % set to zero if Y and X have different signs
    
    % Projection on unit-norm set:
%     X_proj = bsxfun(@times,X_proj, distortion_param.frame_2norm./sqrt(sum(X_proj.^2,1)));

else
    
    error('Measurement function unknown')

end